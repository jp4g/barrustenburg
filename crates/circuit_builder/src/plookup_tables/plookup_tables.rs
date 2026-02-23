//! Native plookup table computation (builder-agnostic).
//!
//! C++ source: barretenberg/cpp/src/barretenberg/stdlib_circuit_builders/plookup_tables/plookup_tables.{hpp,cpp}
//!
//! Key functions:
//!   - `get_lookup_accumulators()`: Slices inputs, computes table values, builds accumulator arrays.
//!   - `get_multitable()`: Returns a MultiTable by ID (lazy initialization of all tables).
//!   - `create_basic_table()`: Factory function for BasicTable creation by ID.

use std::sync::LazyLock;

use bbrs_ecc::curves::bn254::Fr;
use bbrs_numeric::uint256::U256Ext;
use bbrs_numeric::U256;

use super::types::*;
use super::{aes128, blake2s, dummy, fixed_base, keccak, non_native_group_generator, sha256, sparse, uint};

// ---------------------------------------------------------------------------
// Input slicing
// ---------------------------------------------------------------------------

/// Decompose a field element into variable-base slices.
///
/// For each `slice_size` in the list, extracts `accumulator % slice_size` as one slice,
/// then divides `accumulator /= slice_size` for the next iteration.
///
/// Port of C++ `numeric::slice_input_using_variable_bases`.
pub fn slice_input_using_variable_bases(input: Fr, slice_sizes: &[u64]) -> Vec<u64> {
    let standard = input.from_montgomery_form();
    let mut accumulator = U256::from_limbs(standard.data);
    let mut result = Vec::with_capacity(slice_sizes.len());

    for &size in slice_sizes {
        let divisor = U256::from_limbs([size, 0, 0, 0]);
        let (quotient, remainder) = accumulator.div_rem(&divisor.to_nz().unwrap());
        result.push(remainder.as_words()[0]);
        accumulator = quotient;
    }

    result
}

// ---------------------------------------------------------------------------
// Multi-table registry (lazy static)
// ---------------------------------------------------------------------------

static MULTI_TABLES: LazyLock<Vec<MultiTable>> = LazyLock::new(|| {
    let mut tables = vec![MultiTable::empty(); MultiTableId::NUM_MULTI_TABLES];

    tables[MultiTableId::SHA256_CH_INPUT.0] =
        sha256::get_choose_input_table(MultiTableId::SHA256_CH_INPUT);
    tables[MultiTableId::SHA256_CH_OUTPUT.0] =
        sha256::get_choose_output_table(MultiTableId::SHA256_CH_OUTPUT);
    tables[MultiTableId::SHA256_MAJ_INPUT.0] =
        sha256::get_majority_input_table(MultiTableId::SHA256_MAJ_INPUT);
    tables[MultiTableId::SHA256_MAJ_OUTPUT.0] =
        sha256::get_majority_output_table(MultiTableId::SHA256_MAJ_OUTPUT);
    tables[MultiTableId::SHA256_WITNESS_INPUT.0] =
        sha256::get_witness_extension_input_table(MultiTableId::SHA256_WITNESS_INPUT);
    tables[MultiTableId::SHA256_WITNESS_OUTPUT.0] =
        sha256::get_witness_extension_output_table(MultiTableId::SHA256_WITNESS_OUTPUT);

    tables[MultiTableId::AES_NORMALIZE.0] =
        aes128::get_aes_normalization_table(MultiTableId::AES_NORMALIZE);
    tables[MultiTableId::AES_INPUT.0] = aes128::get_aes_input_table(MultiTableId::AES_INPUT);
    tables[MultiTableId::AES_SBOX.0] = aes128::get_aes_sbox_table(MultiTableId::AES_SBOX);

    tables[MultiTableId::UINT8_XOR.0] = uint::get_uint_xor_table(8, MultiTableId::UINT8_XOR);
    tables[MultiTableId::UINT16_XOR.0] = uint::get_uint_xor_table(16, MultiTableId::UINT16_XOR);
    tables[MultiTableId::UINT32_XOR.0] = uint::get_uint_xor_table(32, MultiTableId::UINT32_XOR);
    tables[MultiTableId::UINT64_XOR.0] = uint::get_uint_xor_table(64, MultiTableId::UINT64_XOR);
    tables[MultiTableId::UINT8_AND.0] = uint::get_uint_and_table(8, MultiTableId::UINT8_AND);
    tables[MultiTableId::UINT16_AND.0] = uint::get_uint_and_table(16, MultiTableId::UINT16_AND);
    tables[MultiTableId::UINT32_AND.0] = uint::get_uint_and_table(32, MultiTableId::UINT32_AND);
    tables[MultiTableId::UINT64_AND.0] = uint::get_uint_and_table(64, MultiTableId::UINT64_AND);

    tables[MultiTableId::SECP256K1_XLO.0] = non_native_group_generator::get_xlo_multi_table(
        MultiTableId::SECP256K1_XLO,
        BasicTableId::SECP256K1_XLO_BASIC,
    );
    tables[MultiTableId::SECP256K1_XHI.0] = non_native_group_generator::get_xhi_multi_table(
        MultiTableId::SECP256K1_XHI,
        BasicTableId::SECP256K1_XHI_BASIC,
    );
    tables[MultiTableId::SECP256K1_YLO.0] = non_native_group_generator::get_ylo_multi_table(
        MultiTableId::SECP256K1_YLO,
        BasicTableId::SECP256K1_YLO_BASIC,
    );
    tables[MultiTableId::SECP256K1_YHI.0] = non_native_group_generator::get_yhi_multi_table(
        MultiTableId::SECP256K1_YHI,
        BasicTableId::SECP256K1_YHI_BASIC,
    );
    tables[MultiTableId::SECP256K1_XYPRIME.0] =
        non_native_group_generator::get_xyprime_multi_table(
            MultiTableId::SECP256K1_XYPRIME,
            BasicTableId::SECP256K1_XYPRIME_BASIC,
        );
    tables[MultiTableId::SECP256K1_XLO_ENDO.0] =
        non_native_group_generator::get_xlo_endo_multi_table(
            MultiTableId::SECP256K1_XLO_ENDO,
            BasicTableId::SECP256K1_XLO_ENDO_BASIC,
        );
    tables[MultiTableId::SECP256K1_XHI_ENDO.0] =
        non_native_group_generator::get_xhi_endo_multi_table(
            MultiTableId::SECP256K1_XHI_ENDO,
            BasicTableId::SECP256K1_XHI_ENDO_BASIC,
        );
    tables[MultiTableId::SECP256K1_XYPRIME_ENDO.0] =
        non_native_group_generator::get_xyprime_endo_multi_table(
            MultiTableId::SECP256K1_XYPRIME_ENDO,
            BasicTableId::SECP256K1_XYPRIME_ENDO_BASIC,
        );

    tables[MultiTableId::BLAKE_XOR.0] =
        blake2s::get_blake2s_xor_table(MultiTableId::BLAKE_XOR);
    tables[MultiTableId::BLAKE_XOR_ROTATE_16.0] =
        blake2s::get_blake2s_xor_rotate_16_table(MultiTableId::BLAKE_XOR_ROTATE_16);
    tables[MultiTableId::BLAKE_XOR_ROTATE_8.0] =
        blake2s::get_blake2s_xor_rotate_8_table(MultiTableId::BLAKE_XOR_ROTATE_8);
    tables[MultiTableId::BLAKE_XOR_ROTATE_7.0] =
        blake2s::get_blake2s_xor_rotate_7_table(MultiTableId::BLAKE_XOR_ROTATE_7);

    tables[MultiTableId::KECCAK_FORMAT_INPUT.0] =
        keccak::input::get_keccak_input_table(MultiTableId::KECCAK_FORMAT_INPUT);
    tables[MultiTableId::KECCAK_THETA_OUTPUT.0] =
        keccak::theta::get_theta_output_table(MultiTableId::KECCAK_THETA_OUTPUT);
    tables[MultiTableId::KECCAK_CHI_OUTPUT.0] =
        keccak::chi::get_chi_output_table(MultiTableId::KECCAK_CHI_OUTPUT);
    tables[MultiTableId::KECCAK_FORMAT_OUTPUT.0] =
        keccak::output::get_keccak_output_table(MultiTableId::KECCAK_FORMAT_OUTPUT);

    tables[MultiTableId::FIXED_BASE_LEFT_LO.0] =
        fixed_base::get_fixed_base_table(0, 128, MultiTableId::FIXED_BASE_LEFT_LO);
    tables[MultiTableId::FIXED_BASE_LEFT_HI.0] =
        fixed_base::get_fixed_base_table(1, 126, MultiTableId::FIXED_BASE_LEFT_HI);
    tables[MultiTableId::FIXED_BASE_RIGHT_LO.0] =
        fixed_base::get_fixed_base_table(2, 128, MultiTableId::FIXED_BASE_RIGHT_LO);
    tables[MultiTableId::FIXED_BASE_RIGHT_HI.0] =
        fixed_base::get_fixed_base_table(3, 126, MultiTableId::FIXED_BASE_RIGHT_HI);

    // Keccak normalize-and-rotate: 25 lanes
    for i in 0..25 {
        tables[MultiTableId::KECCAK_NORMALIZE_AND_ROTATE.0 + i] =
            keccak::rho::get_rho_output_table(MultiTableId::KECCAK_NORMALIZE_AND_ROTATE, i);
    }

    tables[MultiTableId::HONK_DUMMY_MULTI.0] = dummy::get_honk_dummy_multitable();

    tables
});

/// Return the multitable with the provided ID.
///
/// Constructs all MultiTables on first access (they're lightweight metadata objects).
///
/// Port of C++ `plookup::get_multitable`.
pub fn get_multitable(id: MultiTableId) -> &'static MultiTable {
    &MULTI_TABLES[id.0]
}

// ---------------------------------------------------------------------------
// Lookup accumulator computation
// ---------------------------------------------------------------------------

/// Given a table ID and key(s) for a key-value lookup, return the lookup accumulators.
///
/// Slices inputs into variable-base limbs, looks up per-limb values, then
/// reconstructs running accumulators so that the first entry contains the
/// full accumulated value. This is the core native lookup logic.
///
/// Port of C++ `plookup::get_lookup_accumulators`.
pub fn get_lookup_accumulators(
    id: MultiTableId,
    key_a: Fr,
    key_b: Fr,
    is_2_to_1_lookup: bool,
) -> ReadData<Fr> {
    let multi_table = get_multitable(id);
    let num_lookups = multi_table.basic_table_ids.len();

    let key_a_slices = slice_input_using_variable_bases(key_a, &multi_table.slice_sizes);
    let key_b_slices = slice_input_using_variable_bases(key_b, &multi_table.slice_sizes);

    let mut column_1_raw_values = Vec::with_capacity(num_lookups);
    let mut column_2_raw_values = Vec::with_capacity(num_lookups);
    let mut column_3_raw_values = Vec::with_capacity(num_lookups);

    let mut lookup = ReadData::<Fr>::new();

    for i in 0..num_lookups {
        let values = (multi_table.get_table_values[i])([key_a_slices[i], key_b_slices[i]]);

        column_1_raw_values.push(Fr::from(key_a_slices[i]));
        column_2_raw_values.push(if is_2_to_1_lookup {
            Fr::from(key_b_slices[i])
        } else {
            values[0]
        });
        column_3_raw_values.push(if is_2_to_1_lookup {
            values[0]
        } else {
            values[1]
        });

        let lookup_entry = LookupEntry {
            key: [
                U256::from_limbs([key_a_slices[i], 0, 0, 0]),
                U256::from_limbs([key_b_slices[i], 0, 0, 0]),
            ],
            value: values,
        };
        lookup.lookup_entries.push(lookup_entry);
    }

    // Build accumulators: row[i-1] = raw[i-1] + row[i] * step_size[i]
    let c1 = lookup.column_mut(ColumnIdx::C1);
    c1.resize(num_lookups, Fr::zero());
    let c2 = lookup.column_mut(ColumnIdx::C2);
    c2.resize(num_lookups, Fr::zero());
    let c3 = lookup.column_mut(ColumnIdx::C3);
    c3.resize(num_lookups, Fr::zero());

    lookup.column_mut(ColumnIdx::C1)[num_lookups - 1] = column_1_raw_values[num_lookups - 1];
    lookup.column_mut(ColumnIdx::C2)[num_lookups - 1] = column_2_raw_values[num_lookups - 1];
    lookup.column_mut(ColumnIdx::C3)[num_lookups - 1] = column_3_raw_values[num_lookups - 1];

    for i in (1..num_lookups).rev() {
        let prev_c1 = lookup.column(ColumnIdx::C1)[i];
        lookup.column_mut(ColumnIdx::C1)[i - 1] =
            column_1_raw_values[i - 1] + prev_c1 * multi_table.column_1_step_sizes[i];

        let prev_c2 = lookup.column(ColumnIdx::C2)[i];
        lookup.column_mut(ColumnIdx::C2)[i - 1] =
            column_2_raw_values[i - 1] + prev_c2 * multi_table.column_2_step_sizes[i];

        let prev_c3 = lookup.column(ColumnIdx::C3)[i];
        lookup.column_mut(ColumnIdx::C3)[i - 1] =
            column_3_raw_values[i - 1] + prev_c3 * multi_table.column_3_step_sizes[i];
    }

    lookup
}

// ---------------------------------------------------------------------------
// Basic table factory
// ---------------------------------------------------------------------------

/// Create a basic table from its ID and table index.
///
/// Port of C++ `plookup::create_basic_table`.
pub fn create_basic_table(id: BasicTableId, index: usize) -> BasicTable {
    let id_val = id.0;

    // Fixed-base table ranges
    if id_val >= BasicTableId::FIXED_BASE_0_0.0 && id_val < BasicTableId::FIXED_BASE_1_0.0 {
        return fixed_base::generate_basic_fixed_base_table(
            0,
            id,
            index,
            id_val - BasicTableId::FIXED_BASE_0_0.0,
        );
    }
    if id_val >= BasicTableId::FIXED_BASE_1_0.0 && id_val < BasicTableId::FIXED_BASE_2_0.0 {
        return fixed_base::generate_basic_fixed_base_table(
            1,
            id,
            index,
            id_val - BasicTableId::FIXED_BASE_1_0.0,
        );
    }
    if id_val >= BasicTableId::FIXED_BASE_2_0.0 && id_val < BasicTableId::FIXED_BASE_3_0.0 {
        return fixed_base::generate_basic_fixed_base_table(
            2,
            id,
            index,
            id_val - BasicTableId::FIXED_BASE_2_0.0,
        );
    }
    if id_val >= BasicTableId::FIXED_BASE_3_0.0 && id_val < BasicTableId::HONK_DUMMY_BASIC1.0 {
        return fixed_base::generate_basic_fixed_base_table(
            3,
            id,
            index,
            id_val - BasicTableId::FIXED_BASE_3_0.0,
        );
    }

    match id {
        BasicTableId::AES_SPARSE_MAP => {
            sparse::generate_sparse_table_with_rotation(9, 8, 0, id, index)
        }
        BasicTableId::AES_SBOX_MAP => aes128::generate_aes_sbox_table(id, index),
        BasicTableId::AES_SPARSE_NORMALIZE => {
            aes128::generate_aes_sparse_normalization_table(id, index)
        }
        BasicTableId::SHA256_WITNESS_NORMALIZE => {
            sha256::generate_witness_extension_normalization_table(id, index)
        }
        BasicTableId::SHA256_WITNESS_SLICE_3 => {
            sparse::generate_sparse_table_with_rotation(16, 3, 0, id, index)
        }
        BasicTableId::SHA256_WITNESS_SLICE_7_ROTATE_4 => {
            sparse::generate_sparse_table_with_rotation(16, 7, 4, id, index)
        }
        BasicTableId::SHA256_WITNESS_SLICE_8_ROTATE_7 => {
            sparse::generate_sparse_table_with_rotation(16, 8, 7, id, index)
        }
        BasicTableId::SHA256_WITNESS_SLICE_14_ROTATE_1 => {
            sparse::generate_sparse_table_with_rotation(16, 14, 1, id, index)
        }
        BasicTableId::SHA256_CH_NORMALIZE => {
            sha256::generate_choose_normalization_table(id, index)
        }
        BasicTableId::SHA256_MAJ_NORMALIZE => {
            sha256::generate_majority_normalization_table(id, index)
        }
        BasicTableId::SHA256_BASE28 => {
            sparse::generate_sparse_table_with_rotation(28, 11, 0, id, index)
        }
        BasicTableId::SHA256_BASE28_ROTATE6 => {
            sparse::generate_sparse_table_with_rotation(28, 11, 6, id, index)
        }
        BasicTableId::SHA256_BASE28_ROTATE3 => {
            sparse::generate_sparse_table_with_rotation(28, 11, 3, id, index)
        }
        BasicTableId::SHA256_BASE16 => {
            sparse::generate_sparse_table_with_rotation(16, 11, 0, id, index)
        }
        BasicTableId::SHA256_BASE16_ROTATE2 => {
            sparse::generate_sparse_table_with_rotation(16, 11, 2, id, index)
        }
        BasicTableId::UINT_XOR_SLICE_6_ROTATE_0 => {
            uint::generate_xor_rotate_table(6, id, index)
        }
        BasicTableId::UINT_XOR_SLICE_4_ROTATE_0 => {
            uint::generate_xor_rotate_table(4, id, index)
        }
        BasicTableId::UINT_XOR_SLICE_2_ROTATE_0 => {
            uint::generate_xor_rotate_table(2, id, index)
        }
        BasicTableId::UINT_AND_SLICE_6_ROTATE_0 => {
            uint::generate_and_rotate_table(6, id, index)
        }
        BasicTableId::UINT_AND_SLICE_4_ROTATE_0 => {
            uint::generate_and_rotate_table(4, id, index)
        }
        BasicTableId::UINT_AND_SLICE_2_ROTATE_0 => {
            uint::generate_and_rotate_table(2, id, index)
        }
        BasicTableId::SECP256K1_XLO_BASIC => {
            non_native_group_generator::generate_xlo_table(id, index)
        }
        BasicTableId::SECP256K1_XHI_BASIC => {
            non_native_group_generator::generate_xhi_table(id, index)
        }
        BasicTableId::SECP256K1_YLO_BASIC => {
            non_native_group_generator::generate_ylo_table(id, index)
        }
        BasicTableId::SECP256K1_YHI_BASIC => {
            non_native_group_generator::generate_yhi_table(id, index)
        }
        BasicTableId::SECP256K1_XYPRIME_BASIC => {
            non_native_group_generator::generate_xyprime_table(id, index)
        }
        BasicTableId::SECP256K1_XLO_ENDO_BASIC => {
            non_native_group_generator::generate_xlo_endo_table(id, index)
        }
        BasicTableId::SECP256K1_XHI_ENDO_BASIC => {
            non_native_group_generator::generate_xhi_endo_table(id, index)
        }
        BasicTableId::SECP256K1_XYPRIME_ENDO_BASIC => {
            non_native_group_generator::generate_xyprime_endo_table(id, index)
        }
        BasicTableId::BLAKE_XOR_ROTATE0 => {
            blake2s::generate_xor_rotate_table(6, 0, false, id, index)
        }
        BasicTableId::BLAKE_XOR_ROTATE0_SLICE5_MOD4 => {
            blake2s::generate_xor_rotate_table(5, 0, true, id, index)
        }
        BasicTableId::BLAKE_XOR_ROTATE1 => {
            blake2s::generate_xor_rotate_table(6, 1, false, id, index)
        }
        BasicTableId::BLAKE_XOR_ROTATE2 => {
            blake2s::generate_xor_rotate_table(6, 2, false, id, index)
        }
        BasicTableId::BLAKE_XOR_ROTATE4 => {
            blake2s::generate_xor_rotate_table(6, 4, false, id, index)
        }
        BasicTableId::HONK_DUMMY_BASIC1 => {
            dummy::generate_honk_dummy_table_basic1(id, index)
        }
        BasicTableId::HONK_DUMMY_BASIC2 => {
            dummy::generate_honk_dummy_table_basic2(id, index)
        }
        BasicTableId::KECCAK_INPUT => keccak::input::generate_keccak_input_table(id, index),
        BasicTableId::KECCAK_THETA => {
            keccak::theta::generate_theta_renormalization_table(id, index)
        }
        BasicTableId::KECCAK_CHI => keccak::chi::generate_chi_renormalization_table(id, index),
        BasicTableId::KECCAK_OUTPUT => {
            keccak::output::generate_keccak_output_table(id, index)
        }
        BasicTableId::KECCAK_RHO_1 => {
            keccak::rho::generate_rho_renormalization_table(id, index, 1)
        }
        BasicTableId::KECCAK_RHO_2 => {
            keccak::rho::generate_rho_renormalization_table(id, index, 2)
        }
        BasicTableId::KECCAK_RHO_3 => {
            keccak::rho::generate_rho_renormalization_table(id, index, 3)
        }
        BasicTableId::KECCAK_RHO_4 => {
            keccak::rho::generate_rho_renormalization_table(id, index, 4)
        }
        BasicTableId::KECCAK_RHO_5 => {
            keccak::rho::generate_rho_renormalization_table(id, index, 5)
        }
        BasicTableId::KECCAK_RHO_6 => {
            keccak::rho::generate_rho_renormalization_table(id, index, 6)
        }
        BasicTableId::KECCAK_RHO_7 => {
            keccak::rho::generate_rho_renormalization_table(id, index, 7)
        }
        BasicTableId::KECCAK_RHO_8 => {
            keccak::rho::generate_rho_renormalization_table(id, index, 8)
        }
        _ => panic!("table id {} does not exist", id.0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_slice_input_using_variable_bases_uniform() {
        // Decompose 239 with base [4, 16, 4]: 239 % 4 = 3, (239/4)%16 = 11, (239/64)%4 = 3
        let input = Fr::from(239u64);
        let slices = slice_input_using_variable_bases(input, &[4, 16, 4]);
        assert_eq!(slices, vec![3, 11, 3]);
    }

    #[test]
    fn test_slice_input_using_variable_bases_powers_of_2() {
        // Decompose with uniform base of 64 (6-bit slices)
        let input = Fr::from(0xABCDu64);
        let slices = slice_input_using_variable_bases(input, &[64, 64, 64]);
        // 0xABCD = 43981 = 13 + 47*64 + 10*4096
        assert_eq!(slices[0], 43981 % 64);
        assert_eq!(slices[1], (43981 / 64) % 64);
        assert_eq!(slices[2], (43981 / 4096) % 64);
    }

    #[test]
    fn test_get_multitable_uint32_xor() {
        let table = get_multitable(MultiTableId::UINT32_XOR);
        // UINT32_XOR should have ceil(32/6) = 6 basic tables
        assert_eq!(table.basic_table_ids.len(), 6);
    }

    #[test]
    fn test_get_lookup_accumulators_uint32_xor() {
        let left = Fr::from(0x12345678u64);
        let right = Fr::from(0xDEADBEEFu64);
        let lookup = get_lookup_accumulators(
            MultiTableId::UINT32_XOR,
            left,
            right,
            true,
        );

        // First accumulator entry should be the full accumulated value
        let expected_xor = 0x12345678u64 ^ 0xDEADBEEFu64;
        assert_eq!(
            lookup.column(ColumnIdx::C3)[0],
            Fr::from(expected_xor)
        );
        assert_eq!(lookup.column(ColumnIdx::C1)[0], left);
        assert_eq!(lookup.column(ColumnIdx::C2)[0], right);
    }

    #[test]
    fn test_get_lookup_accumulators_uint32_and() {
        let left = Fr::from(0x12345678u64);
        let right = Fr::from(0xDEADBEEFu64);
        let lookup = get_lookup_accumulators(
            MultiTableId::UINT32_AND,
            left,
            right,
            true,
        );

        let expected_and = 0x12345678u64 & 0xDEADBEEFu64;
        assert_eq!(
            lookup.column(ColumnIdx::C3)[0],
            Fr::from(expected_and)
        );
    }

    #[test]
    fn test_create_basic_table_xor() {
        let table = create_basic_table(BasicTableId::UINT_XOR_SLICE_6_ROTATE_0, 0);
        assert_eq!(table.id, BasicTableId::UINT_XOR_SLICE_6_ROTATE_0);
        // 6-bit XOR table: 64 * 64 = 4096 entries
        assert_eq!(table.size(), 4096);
    }

    #[test]
    fn test_create_basic_table_honk_dummy() {
        let table = create_basic_table(BasicTableId::HONK_DUMMY_BASIC1, 0);
        assert_eq!(table.id, BasicTableId::HONK_DUMMY_BASIC1);
        assert!(table.size() > 0);
    }
}
