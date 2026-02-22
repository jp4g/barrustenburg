//! AES-128 plookup tables.
//!
//! C++ source: plookup_tables/aes128.hpp
//!
//! Provides lookup tables for AES-128 operations:
//! - Sparse form (base-9) mapping for AES bytes
//! - Sparse normalization tables
//! - S-box substitution tables with MixColumns swizzle

use bbrs_ecc::curves::bn254::Fr;
use bbrs_numeric::bitop::pow64;

use super::sparse;
use super::types::{BasicTable, BasicTableId, MultiTable, MultiTableId};

// ---------------------------------------------------------------------------
// Helper: Fr exponentiation with u64 exponent
// ---------------------------------------------------------------------------

/// Raise an Fr element to a u64 power.
#[inline]
fn fr_pow(base: &Fr, exp: u64) -> Fr {
    base.pow(&[exp, 0, 0, 0])
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Base used for AES sparse representation.
const AES_BASE: u64 = 9;

/// AES normalization table (9 entries): only digit value 0 maps to 1; all
/// others map to 0. This is used to normalize sparse-form digits.
/// Retained for documentation purposes, matching C++ `aes_normalization_table`.
#[allow(dead_code)]
const AES_NORMALIZATION_TABLE: [u64; 9] = [1, 0, 0, 0, 0, 0, 0, 0, 0];

/// AES S-box lookup table (256 entries).
///
/// Standard AES SubBytes substitution table.
#[rustfmt::skip]
const AES_SBOX: [u8; 256] = [
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16,
];

// ---------------------------------------------------------------------------
// Sparse form helpers
// ---------------------------------------------------------------------------

/// Map an 8-bit binary integer into sparse base-9 form.
///
/// Each bit of `input` contributes `9^bit_position` to the output.
/// For 8-bit values, the result fits in u64 (9^7 = 4782969).
fn map_into_sparse_form_base9(input: u8) -> u64 {
    let mut out: u64 = 0;
    let mut base_power: u64 = 1;
    for i in 0..8 {
        if (input >> i) & 1 != 0 {
            out += base_power;
        }
        if i < 7 {
            base_power *= AES_BASE;
        }
    }
    out
}

/// Map a sparse base-9 value back to an 8-bit binary integer.
///
/// Decomposes the input into base-9 digits and checks if each is odd
/// (representing a set bit). Uses the same algorithm as the C++ version
/// which processes from the most significant digit down.
///
/// Mirrors C++ `numeric::map_from_sparse_form<AES_BASE>(input)`.
fn map_from_sparse_form_base9(input: u64) -> u8 {
    // Simple approach: decompose into base-9 digits, check odd/even.
    let mut target = input;
    let mut output: u8 = 0;
    for i in 0..8 {
        let digit = target % AES_BASE;
        if (digit & 1) == 1 {
            output |= 1u8 << i;
        }
        target /= AES_BASE;
    }
    output
}

// ---------------------------------------------------------------------------
// get_values_from_key functions
// ---------------------------------------------------------------------------

/// Compute sparse form of the key for AES sparse mapping.
///
/// Returns `[sparse_form(key[0]), 0]`.
///
/// Mirrors C++ `get_aes_sparse_values_from_key`.
pub fn get_aes_sparse_values_from_key(key: [u64; 2]) -> [Fr; 2] {
    let sparse = map_into_sparse_form_base9(key[0] as u8);
    [Fr::from(sparse), Fr::zero()]
}

/// Compute sparse normalization values for AES.
///
/// Converts sparse-form input back to binary, then re-encodes it in sparse form.
/// This normalizes the digit values (which may have been accumulated beyond {0,1}).
///
/// Returns `[sparse_form(map_from_sparse(key[0])), 0]`.
///
/// Mirrors C++ `get_aes_sparse_normalization_values_from_key`.
pub fn get_aes_sparse_normalization_values_from_key(key: [u64; 2]) -> [Fr; 2] {
    let byte = map_from_sparse_form_base9(key[0]);
    [Fr::from(map_into_sparse_form_base9(byte)), Fr::zero()]
}

/// Compute S-box output values for AES.
///
/// Given a sparse-form input, decodes the byte, looks up the S-box,
/// computes the MixColumns swizzle, and returns both in sparse form.
///
/// Returns `[sparse(sbox[byte]), sparse(sbox[byte] ^ swizzled)]`.
///
/// Mirrors C++ `get_aes_sbox_values_from_key`.
pub fn get_aes_sbox_values_from_key(key: [u64; 2]) -> [Fr; 2] {
    let byte = map_from_sparse_form_base9(key[0]);
    let sbox_value = AES_SBOX[byte as usize];
    // MixColumns: xtime(sbox_value) = (sbox_value << 1) ^ ((sbox_value >> 7) & 1) * 0x1b)
    let swizzled = (sbox_value << 1) ^ (((sbox_value >> 7) & 1) * 0x1b);
    [
        Fr::from(map_into_sparse_form_base9(sbox_value)),
        Fr::from(map_into_sparse_form_base9(sbox_value ^ swizzled)),
    ]
}

// ---------------------------------------------------------------------------
// BasicTable generators
// ---------------------------------------------------------------------------

/// Generate the AES sparse mapping table.
///
/// Maps each byte [0, 256) to its base-9 sparse form.
///
/// Mirrors C++ `generate_aes_sparse_table`.
pub fn generate_aes_sparse_table(id: BasicTableId, table_index: usize) -> BasicTable {
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = true;

    for i in 0u64..256 {
        let sparse = map_into_sparse_form_base9(i as u8);
        table.column_1.push(Fr::from(i));
        table.column_2.push(Fr::zero());
        table.column_3.push(Fr::from(sparse));
    }

    table.get_values_from_key = get_aes_sparse_values_from_key;

    table.column_1_step_size = Fr::from(256u64);
    table.column_2_step_size = Fr::zero();
    table.column_3_step_size = Fr::zero();

    table
}

/// Generate the AES sparse normalization table.
///
/// Enumerates all 9^4 = 6561 possible 4-digit base-9 values and maps each
/// to its normalized form (where each digit is replaced by digit & 1).
///
/// Mirrors C++ `generate_aes_sparse_normalization_table`.
pub fn generate_aes_sparse_normalization_table(
    id: BasicTableId,
    table_index: usize,
) -> BasicTable {
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    for i in 0..AES_BASE {
        let i_raw = i * AES_BASE * AES_BASE * AES_BASE;
        let i_normalized = (if (i & 1) == 1 { 1u64 } else { 0u64 }) * AES_BASE * AES_BASE * AES_BASE;
        for j in 0..AES_BASE {
            let j_raw = j * AES_BASE * AES_BASE;
            let j_normalized = (if (j & 1) == 1 { 1u64 } else { 0u64 }) * AES_BASE * AES_BASE;
            for k in 0..AES_BASE {
                let k_raw = k * AES_BASE;
                let k_normalized = (if (k & 1) == 1 { 1u64 } else { 0u64 }) * AES_BASE;
                for m in 0..AES_BASE {
                    let m_raw = m;
                    let m_normalized = if (m & 1) == 1 { 1u64 } else { 0u64 };
                    let left = i_raw + j_raw + k_raw + m_raw;
                    let right = i_normalized + j_normalized + k_normalized + m_normalized;
                    table.column_1.push(Fr::from(left));
                    table.column_2.push(Fr::from(right));
                    table.column_3.push(Fr::zero());
                }
            }
        }
    }

    table.get_values_from_key = get_aes_sparse_normalization_values_from_key;

    // 9^4 = 6561
    table.column_1_step_size = Fr::from(6561u64);
    table.column_2_step_size = Fr::from(6561u64);
    table.column_3_step_size = Fr::zero();

    table
}

/// Generate the AES S-box table.
///
/// For each byte [0, 256), stores the sparse-form input, S-box output,
/// and S-box-xor-swizzled output (for MixColumns).
///
/// Mirrors C++ `generate_aes_sbox_table`.
pub fn generate_aes_sbox_table(id: BasicTableId, table_index: usize) -> BasicTable {
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    for i in 0u64..256 {
        let byte = i as u8;
        let first = map_into_sparse_form_base9(byte);
        let sbox_value = AES_SBOX[byte as usize];
        let swizzled = (sbox_value << 1) ^ (((sbox_value >> 7) & 1) * 0x1b);
        let second = map_into_sparse_form_base9(sbox_value);
        let third = map_into_sparse_form_base9(sbox_value ^ swizzled);

        table.column_1.push(Fr::from(first));
        table.column_2.push(Fr::from(second));
        table.column_3.push(Fr::from(third));
    }

    table.get_values_from_key = get_aes_sbox_values_from_key;

    table.column_1_step_size = Fr::zero();
    table.column_2_step_size = Fr::zero();
    table.column_3_step_size = Fr::zero();

    table
}

// ---------------------------------------------------------------------------
// MultiTable constructors
// ---------------------------------------------------------------------------

/// Create the AES normalization multi-table.
///
/// Uses 2 lookups, each covering 4 base-9 digits (9^4 = 6561 entries).
///
/// Mirrors C++ `get_aes_normalization_table`.
pub fn get_aes_normalization_table(id: MultiTableId) -> MultiTable {
    let num_entries = 2;
    let base9_fr = Fr::from(AES_BASE);

    let mut column_1_coefficients = Vec::with_capacity(num_entries);
    let mut column_2_coefficients = Vec::with_capacity(num_entries);
    let mut column_3_coefficients = Vec::with_capacity(num_entries);

    for i in 0..num_entries {
        column_1_coefficients.push(fr_pow(&base9_fr, 4 * i as u64));
        column_2_coefficients.push(fr_pow(&base9_fr, 4 * i as u64));
        column_3_coefficients.push(Fr::zero());
    }

    let mut table = MultiTable::new_explicit(
        column_1_coefficients,
        column_2_coefficients,
        column_3_coefficients,
    );
    table.id = id;

    let slice_size = AES_BASE * AES_BASE * AES_BASE * AES_BASE; // 6561
    for _ in 0..num_entries {
        table.slice_sizes.push(slice_size);
        table.basic_table_ids.push(BasicTableId::AES_SPARSE_NORMALIZE);
        table
            .get_table_values
            .push(get_aes_sparse_normalization_values_from_key);
    }

    table
}

/// Create the AES input multi-table.
///
/// Uses 16 lookups to map 16 bytes into sparse base-9 form (one per byte).
///
/// Mirrors C++ `get_aes_input_table`.
pub fn get_aes_input_table(id: MultiTableId) -> MultiTable {
    let num_entries = 16;

    let mut table = MultiTable::new_repeated(
        Fr::from(256u64),
        Fr::zero(),
        Fr::zero(),
        num_entries,
    );
    table.id = id;

    for _ in 0..num_entries {
        table.slice_sizes.push(256);
        table.basic_table_ids.push(BasicTableId::AES_SPARSE_MAP);
        table
            .get_table_values
            .push(sparse::get_sparse_values_base9_rot0);
    }

    table
}

/// Create the AES S-box multi-table.
///
/// Single lookup covering the full S-box.
///
/// Mirrors C++ `get_aes_sbox_table`.
pub fn get_aes_sbox_table(id: MultiTableId) -> MultiTable {
    let mut table = MultiTable::new_repeated(
        Fr::zero(),
        Fr::zero(),
        Fr::zero(),
        1,
    );
    table.id = id;

    table.slice_sizes.push(pow64(AES_BASE, 8));
    table.basic_table_ids.push(BasicTableId::AES_SBOX_MAP);
    table.get_table_values.push(get_aes_sbox_values_from_key);

    table
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_map_into_sparse_form_base9() {
        // 0 -> 0
        assert_eq!(map_into_sparse_form_base9(0), 0);
        // 1 -> 9^0 = 1
        assert_eq!(map_into_sparse_form_base9(1), 1);
        // 2 -> 9^1 = 9
        assert_eq!(map_into_sparse_form_base9(2), 9);
        // 3 -> 9^0 + 9^1 = 10
        assert_eq!(map_into_sparse_form_base9(3), 10);
        // 0xFF -> sum of 9^i for i=0..7
        let mut expected = 0u64;
        for i in 0..8 {
            expected += pow64(AES_BASE, i);
        }
        assert_eq!(map_into_sparse_form_base9(0xFF), expected);
    }

    #[test]
    fn test_map_from_sparse_form_base9_roundtrip() {
        for byte in 0..=255u8 {
            let sparse = map_into_sparse_form_base9(byte);
            let recovered = map_from_sparse_form_base9(sparse);
            assert_eq!(recovered, byte, "roundtrip failed for byte {}", byte);
        }
    }

    #[test]
    fn test_aes_sbox_known_values() {
        // AES S-box(0x00) = 0x63
        assert_eq!(AES_SBOX[0x00], 0x63);
        // AES S-box(0x01) = 0x7c
        assert_eq!(AES_SBOX[0x01], 0x7c);
        // AES S-box(0x52) = 0x00 (inverse of first entry)
        assert_eq!(AES_SBOX[0x52], 0x00);
    }

    #[test]
    fn test_get_aes_sparse_values() {
        let result = get_aes_sparse_values_from_key([0, 0]);
        assert_eq!(result[0], Fr::from(0u64));
        assert_eq!(result[1], Fr::zero());

        let result = get_aes_sparse_values_from_key([1, 0]);
        assert_eq!(result[0], Fr::from(1u64));
    }

    #[test]
    fn test_get_aes_sbox_values() {
        // Input byte = 0 (sparse form of 0 is 0)
        let sparse_0 = map_into_sparse_form_base9(0);
        let result = get_aes_sbox_values_from_key([sparse_0, 0]);
        // S-box(0) = 0x63
        let sbox_val = 0x63u8;
        let swizzled = (sbox_val << 1) ^ (((sbox_val >> 7) & 1) * 0x1b);
        assert_eq!(result[0], Fr::from(map_into_sparse_form_base9(sbox_val)));
        assert_eq!(
            result[1],
            Fr::from(map_into_sparse_form_base9(sbox_val ^ swizzled))
        );
    }

    #[test]
    fn test_generate_aes_sparse_table() {
        let table = generate_aes_sparse_table(BasicTableId::AES_SPARSE_MAP, 0);
        assert_eq!(table.column_1.len(), 256);
        assert!(table.use_twin_keys);
        // First entry: 0 -> sparse(0) = 0
        assert_eq!(table.column_1[0], Fr::from(0u64));
        assert_eq!(table.column_3[0], Fr::from(0u64));
    }

    #[test]
    fn test_generate_aes_normalization_table() {
        let table = generate_aes_sparse_normalization_table(BasicTableId::AES_SPARSE_NORMALIZE, 0);
        // 9^4 = 6561 entries
        assert_eq!(table.column_1.len(), 6561);
        assert!(!table.use_twin_keys);
    }

    #[test]
    fn test_generate_aes_sbox_table() {
        let table = generate_aes_sbox_table(BasicTableId::AES_SBOX_MAP, 0);
        assert_eq!(table.column_1.len(), 256);
        assert!(!table.use_twin_keys);
    }

    #[test]
    fn test_aes_normalization_multitable() {
        let table = get_aes_normalization_table(MultiTableId::AES_NORMALIZE);
        assert_eq!(table.id, MultiTableId::AES_NORMALIZE);
        assert_eq!(table.slice_sizes.len(), 2);
        assert_eq!(table.slice_sizes[0], 6561);
        assert_eq!(table.slice_sizes[1], 6561);
    }

    #[test]
    fn test_aes_input_multitable() {
        let table = get_aes_input_table(MultiTableId::AES_INPUT);
        assert_eq!(table.id, MultiTableId::AES_INPUT);
        assert_eq!(table.slice_sizes.len(), 16);
        for &sz in &table.slice_sizes {
            assert_eq!(sz, 256);
        }
    }

    #[test]
    fn test_aes_sbox_multitable() {
        let table = get_aes_sbox_table(MultiTableId::AES_SBOX);
        assert_eq!(table.id, MultiTableId::AES_SBOX);
        assert_eq!(table.slice_sizes.len(), 1);
        assert_eq!(table.slice_sizes[0], pow64(AES_BASE, 8));
    }
}
