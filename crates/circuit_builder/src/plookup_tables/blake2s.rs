//! BLAKE2s plookup tables.
//!
//! C++ source: plookup_tables/blake2s.hpp
//!
//! Provides XOR-with-rotation lookup tables for the BLAKE2s hash algorithm.
//! BLAKE2s operates on 32-bit words and uses XOR combined with right-rotation
//! by 16, 12, 8, and 7 bits. These tables decompose 32-bit values into 6-bit
//! slices (with a 5-bit last slice) and compute XOR + rotation per slice.

use bbrs_ecc::curves::bn254::Fr;

use super::types::{BasicTable, BasicTableId, MultiTable, MultiTableId};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Number of bits in the last (residual) slice.
const BITS_IN_LAST_SLICE: u64 = 5;

/// Size of the last slice: 2^5 = 32.
const SIZE_OF_LAST_SLICE: u64 = 1u64 << BITS_IN_LAST_SLICE;

// ---------------------------------------------------------------------------
// Helper: 32-bit rotation
// ---------------------------------------------------------------------------

/// Rotate a 32-bit value right by `rotation` bits.
#[inline]
fn rotate32(value: u32, rotation: u32) -> u32 {
    if rotation == 0 {
        value
    } else {
        (value >> rotation) | (value << (32 - rotation))
    }
}

// ---------------------------------------------------------------------------
// get_values_from_key functions
// ---------------------------------------------------------------------------

/// Compute XOR-rotate values for 6-bit slices with no rotation.
pub fn get_xor_rotate_values_6_0(key: [u64; 2]) -> [Fr; 2] {
    let xored = (key[0] as u32) ^ (key[1] as u32);
    let rotated = rotate32(xored, 0);
    [Fr::from(rotated as u64), Fr::zero()]
}

/// Compute XOR-rotate values for 6-bit slices with 1-bit rotation.
pub fn get_xor_rotate_values_6_1(key: [u64; 2]) -> [Fr; 2] {
    let xored = (key[0] as u32) ^ (key[1] as u32);
    let rotated = rotate32(xored, 1);
    [Fr::from(rotated as u64), Fr::zero()]
}

/// Compute XOR-rotate values for 6-bit slices with 2-bit rotation.
pub fn get_xor_rotate_values_6_2(key: [u64; 2]) -> [Fr; 2] {
    let xored = (key[0] as u32) ^ (key[1] as u32);
    let rotated = rotate32(xored, 2);
    [Fr::from(rotated as u64), Fr::zero()]
}

/// Compute XOR-rotate values for 6-bit slices with 4-bit rotation.
pub fn get_xor_rotate_values_6_4(key: [u64; 2]) -> [Fr; 2] {
    let xored = (key[0] as u32) ^ (key[1] as u32);
    let rotated = rotate32(xored, 4);
    [Fr::from(rotated as u64), Fr::zero()]
}

/// Compute XOR-rotate values for the 5-bit last slice with mod-4 filter.
///
/// When filter=true, only the 2 least significant bits of each key participate
/// in the XOR, producing `ROTR^0((key[0] & 3) ^ (key[1] & 3))`.
pub fn get_xor_rotate_values_5_0_filtered(key: [u64; 2]) -> [Fr; 2] {
    let filtered0 = key[0] & 3;
    let filtered1 = key[1] & 3;
    let xored = (filtered0 as u32) ^ (filtered1 as u32);
    let rotated = rotate32(xored, 0);
    [Fr::from(rotated as u64), Fr::zero()]
}

// ---------------------------------------------------------------------------
// BasicTable generators
// ---------------------------------------------------------------------------

/// Generate a basic XOR-rotate table.
///
/// For each pair (i, j) in [0, 2^bits_per_slice), computes
/// `ROTR^num_rotated_bits(i ^ j)`. When `use_mod4_filter` is true,
/// the XOR operates only on the 2 least significant bits.
///
/// Mirrors C++ `generate_xor_rotate_table<bits_per_slice, num_rotated_bits, filter>`.
pub fn generate_xor_rotate_table(
    bits_per_slice: u64,
    num_rotated_bits: u32,
    use_mod4_filter: bool,
    id: BasicTableId,
    table_index: usize,
) -> BasicTable {
    let base = 1u64 << bits_per_slice;
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = true;

    for i in 0..base {
        for j in 0..base {
            table.column_1.push(Fr::from(i));
            table.column_2.push(Fr::from(j));

            let i_val = if use_mod4_filter { i & 3 } else { i };
            let j_val = if use_mod4_filter { j & 3 } else { j };
            let xored = (i_val as u32) ^ (j_val as u32);
            let result = rotate32(xored, num_rotated_bits);
            table.column_3.push(Fr::from(result as u64));
        }
    }

    // Select the appropriate get_values_from_key function
    table.get_values_from_key = select_xor_rotate_get_values(
        bits_per_slice,
        num_rotated_bits,
        use_mod4_filter,
    );

    table.column_1_step_size = Fr::from(base);
    table.column_2_step_size = Fr::from(base);
    table.column_3_step_size = Fr::from(base);

    table
}

/// Select the correct get_values_from_key function based on parameters.
fn select_xor_rotate_get_values(
    bits_per_slice: u64,
    num_rotated_bits: u32,
    use_mod4_filter: bool,
) -> fn([u64; 2]) -> [Fr; 2] {
    if use_mod4_filter {
        return get_xor_rotate_values_5_0_filtered;
    }
    match (bits_per_slice, num_rotated_bits) {
        (6, 0) => get_xor_rotate_values_6_0,
        (6, 1) => get_xor_rotate_values_6_1,
        (6, 2) => get_xor_rotate_values_6_2,
        (6, 4) => get_xor_rotate_values_6_4,
        _ => panic!(
            "No concrete get_values_from_key for blake2s XOR table with bits={}, rot={}",
            bits_per_slice, num_rotated_bits
        ),
    }
}

// ---------------------------------------------------------------------------
// MultiTable constructors
// ---------------------------------------------------------------------------

/// Create the BLAKE2s plain XOR multi-table (no rotation).
///
/// Uses 6 slices: 5 x 6-bit slices + 1 x 5-bit slice (with mod-4 filter).
/// Total: 5*6 + 5 = 35 >= 32 bits covered.
///
/// Mirrors C++ `get_blake2s_xor_table`.
pub fn get_blake2s_xor_table(id: MultiTableId) -> MultiTable {
    let num_entries = (32 + 2) / 6 + 1; // = 6
    let base = 1u64 << 6;

    let mut table = MultiTable::new_repeated(
        Fr::from(base),
        Fr::from(base),
        Fr::from(base),
        num_entries,
    );
    table.id = id;

    for _ in 0..(num_entries - 1) {
        table.slice_sizes.push(base);
        table.basic_table_ids.push(BasicTableId::BLAKE_XOR_ROTATE0);
        table.get_table_values.push(get_xor_rotate_values_6_0);
    }

    // Last slice: 5-bit with mod-4 filter
    table.slice_sizes.push(SIZE_OF_LAST_SLICE);
    table.basic_table_ids.push(BasicTableId::BLAKE_XOR_ROTATE0_SLICE5_MOD4);
    table.get_table_values.push(get_xor_rotate_values_5_0_filtered);

    table
}

/// Create the BLAKE2s XOR-rotate-16 multi-table.
///
/// Computes ROTR^16(a ^ b) using 6 slices with appropriate coefficients.
///
/// Mirrors C++ `get_blake2s_xor_rotate_16_table`.
pub fn get_blake2s_xor_rotate_16_table(id: MultiTableId) -> MultiTable {
    let base = 1u64 << 6;
    let coefficient_16 = Fr::from(1u64 << 16).invert();

    let column_1_coefficients = vec![
        Fr::from(1u64),
        Fr::from(1u64 << 6),
        Fr::from(1u64 << 12),
        Fr::from(1u64 << 18),
        Fr::from(1u64 << 24),
        Fr::from(1u64 << 30),
    ];

    let column_3_coefficients = vec![
        Fr::one(),
        Fr::from(1u64 << 6),
        coefficient_16,
        coefficient_16 * Fr::from(1u64 << 2),
        coefficient_16 * Fr::from(1u64 << 8),
        coefficient_16 * Fr::from(1u64 << 14),
    ];

    let mut table = MultiTable::new_explicit(
        column_1_coefficients.clone(),
        column_1_coefficients,
        column_3_coefficients,
    );
    table.id = id;
    table.slice_sizes = vec![base, base, base, base, base, SIZE_OF_LAST_SLICE];
    table.basic_table_ids = vec![
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE4,
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE0_SLICE5_MOD4,
    ];
    table.get_table_values = vec![
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_6_4,
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_5_0_filtered,
    ];

    table
}

/// Create the BLAKE2s XOR-rotate-8 multi-table.
///
/// Computes ROTR^8(a ^ b) using 6 slices with appropriate coefficients.
///
/// Mirrors C++ `get_blake2s_xor_rotate_8_table`.
pub fn get_blake2s_xor_rotate_8_table(id: MultiTableId) -> MultiTable {
    let base = 1u64 << 6;
    let coefficient_24 = Fr::from(1u64 << 24).invert();

    let column_1_coefficients = vec![
        Fr::from(1u64),
        Fr::from(1u64 << 6),
        Fr::from(1u64 << 12),
        Fr::from(1u64 << 18),
        Fr::from(1u64 << 24),
        Fr::from(1u64 << 30),
    ];

    let column_3_coefficients = vec![
        Fr::one(),
        coefficient_24,
        coefficient_24 * Fr::from(1u64 << 4),
        coefficient_24 * Fr::from(1u64 << (4 + 6)),
        coefficient_24 * Fr::from(1u64 << (4 + 12)),
        coefficient_24 * Fr::from(1u64 << (4 + 18)),
    ];

    let mut table = MultiTable::new_explicit(
        column_1_coefficients.clone(),
        column_1_coefficients,
        column_3_coefficients,
    );
    table.id = id;
    table.slice_sizes = vec![base, base, base, base, base, SIZE_OF_LAST_SLICE];
    table.basic_table_ids = vec![
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE2,
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE0_SLICE5_MOD4,
    ];
    table.get_table_values = vec![
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_6_2,
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_5_0_filtered,
    ];

    table
}

/// Create the BLAKE2s XOR-rotate-7 multi-table.
///
/// Computes ROTR^7(a ^ b) using 6 slices with appropriate coefficients.
///
/// Mirrors C++ `get_blake2s_xor_rotate_7_table`.
pub fn get_blake2s_xor_rotate_7_table(id: MultiTableId) -> MultiTable {
    let base = 1u64 << 6;
    let coefficient_25 = Fr::from(1u64 << 25).invert();

    let column_1_coefficients = vec![
        Fr::from(1u64),
        Fr::from(1u64 << 6),
        Fr::from(1u64 << 12),
        Fr::from(1u64 << 18),
        Fr::from(1u64 << 24),
        Fr::from(1u64 << 30),
    ];

    let column_3_coefficients = vec![
        Fr::one(),
        coefficient_25,
        coefficient_25 * Fr::from(1u64 << 5),
        coefficient_25 * Fr::from(1u64 << (5 + 6)),
        coefficient_25 * Fr::from(1u64 << (5 + 12)),
        coefficient_25 * Fr::from(1u64 << (5 + 18)),
    ];

    let mut table = MultiTable::new_explicit(
        column_1_coefficients.clone(),
        column_1_coefficients,
        column_3_coefficients,
    );
    table.id = id;
    table.slice_sizes = vec![base, base, base, base, base, SIZE_OF_LAST_SLICE];
    table.basic_table_ids = vec![
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE1,
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE0,
        BasicTableId::BLAKE_XOR_ROTATE0_SLICE5_MOD4,
    ];
    table.get_table_values = vec![
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_6_1,
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_6_0,
        get_xor_rotate_values_5_0_filtered,
    ];

    table
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rotate32() {
        assert_eq!(rotate32(0x12345678, 0), 0x12345678);
        assert_eq!(rotate32(1, 1), 0x80000000);
        assert_eq!(rotate32(0x80000001, 1), 0xC0000000);
    }

    #[test]
    fn test_xor_rotate_values_no_rotation() {
        let result = get_xor_rotate_values_6_0([0b101010, 0b010101]);
        // 0b101010 ^ 0b010101 = 0b111111 = 63
        assert_eq!(result[0], Fr::from(63u64));
        assert_eq!(result[1], Fr::zero());
    }

    #[test]
    fn test_xor_rotate_values_filtered() {
        // Only the bottom 2 bits of each key matter
        let result = get_xor_rotate_values_5_0_filtered([0b11111, 0b11101]);
        // (0b11111 & 3) ^ (0b11101 & 3) = 0b11 ^ 0b01 = 0b10 = 2
        assert_eq!(result[0], Fr::from(2u64));
    }

    #[test]
    fn test_generate_xor_rotate_table_basic() {
        let table = generate_xor_rotate_table(
            6, 0, false,
            BasicTableId::BLAKE_XOR_ROTATE0,
            0,
        );
        // 2^6 * 2^6 = 4096 entries
        assert_eq!(table.column_1.len(), 4096);
        assert!(table.use_twin_keys);
    }

    #[test]
    fn test_generate_xor_rotate_table_filtered() {
        let table = generate_xor_rotate_table(
            BITS_IN_LAST_SLICE, 0, true,
            BasicTableId::BLAKE_XOR_ROTATE0_SLICE5_MOD4,
            0,
        );
        // 2^5 * 2^5 = 1024 entries
        assert_eq!(table.column_1.len(), 1024);

        // For i=7, j=5 with filter: (7&3)^(5&3) = 3^1 = 2
        let idx = 7 * 32 + 5;
        assert_eq!(table.column_3[idx], Fr::from(2u64));
    }

    #[test]
    fn test_blake2s_xor_table() {
        let table = get_blake2s_xor_table(MultiTableId::BLAKE_XOR);
        assert_eq!(table.id, MultiTableId::BLAKE_XOR);
        assert_eq!(table.slice_sizes.len(), 6);
        // First 5 slices are 64, last is 32
        for i in 0..5 {
            assert_eq!(table.slice_sizes[i], 64);
        }
        assert_eq!(table.slice_sizes[5], SIZE_OF_LAST_SLICE);
    }

    #[test]
    fn test_blake2s_xor_rotate_16_table() {
        let table = get_blake2s_xor_rotate_16_table(MultiTableId::BLAKE_XOR_ROTATE_16);
        assert_eq!(table.id, MultiTableId::BLAKE_XOR_ROTATE_16);
        assert_eq!(table.slice_sizes.len(), 6);
        assert_eq!(table.basic_table_ids.len(), 6);
        assert_eq!(table.get_table_values.len(), 6);
    }

    #[test]
    fn test_blake2s_xor_rotate_8_table() {
        let table = get_blake2s_xor_rotate_8_table(MultiTableId::BLAKE_XOR_ROTATE_8);
        assert_eq!(table.id, MultiTableId::BLAKE_XOR_ROTATE_8);
        assert_eq!(table.slice_sizes.len(), 6);
    }

    #[test]
    fn test_blake2s_xor_rotate_7_table() {
        let table = get_blake2s_xor_rotate_7_table(MultiTableId::BLAKE_XOR_ROTATE_7);
        assert_eq!(table.id, MultiTableId::BLAKE_XOR_ROTATE_7);
        assert_eq!(table.slice_sizes.len(), 6);
    }
}
