//! SHA-256 plookup tables.
//!
//! C++ source: plookup_tables/sha256.hpp
//!
//! Provides normalization tables for SHA-256 choose, majority, and witness
//! extension operations, as well as multi-tables that assemble sparse-form
//! representations with rotation coefficients.

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
// Normalization tables
// ---------------------------------------------------------------------------

/// Choose normalization table (28 entries).
///
/// Maps the weighted sum `e + 2f + 3g` (combined with XOR result digit)
/// to the normalized choose output bit.
const CHOOSE_NORMALIZATION_TABLE: [u64; 28] = [
    // xor result = 0
    0, // e + 2f + 3g = 0
    0, // e + 2f + 3g = 1
    0, // e + 2f + 3g = 2
    1, // e + 2f + 3g = 3
    0, // e + 2f + 3g = 4
    1, // e + 2f + 3g = 5
    1, // e + 2f + 3g = 6
    // xor result = 1
    1, // e + 2f + 3g = 0
    1, // e + 2f + 3g = 1
    1, // e + 2f + 3g = 2
    2, // e + 2f + 3g = 3
    1, // e + 2f + 3g = 4
    2, // e + 2f + 3g = 5
    2, // e + 2f + 3g = 6
    // xor result = 2
    0, // e + 2f + 3g = 0
    0, // e + 2f + 3g = 1
    0, // e + 2f + 3g = 2
    1, // e + 2f + 3g = 3
    0, // e + 2f + 3g = 4
    1, // e + 2f + 3g = 5
    1, // e + 2f + 3g = 6
    1, // e + 2f + 3g = 0 (xor result = 3 start)
    // xor result = 3
    1, // e + 2f + 3g = 1
    1, // e + 2f + 3g = 2
    2, // e + 2f + 3g = 3
    1, // e + 2f + 3g = 4
    2, // e + 2f + 3g = 5
    2, // e + 2f + 3g = 6
];

/// Majority normalization table (16 entries).
///
/// Maps `a + b + c` (combined with XOR result digit) to the majority output bit.
const MAJORITY_NORMALIZATION_TABLE: [u64; 16] = [
    // xor result = 0
    0, 0, 1, 1,
    // xor result = 1
    1, 1, 2, 2,
    // xor result = 2
    0, 0, 1, 1,
    // xor result = 3
    1, 1, 2, 2,
];

/// Witness extension normalization table (16 entries).
const WITNESS_EXTENSION_NORMALIZATION_TABLE: [u64; 16] = [
    // xor result = 0
    0, 1, 0, 1,
    // xor result = 1
    1, 2, 1, 2,
    // xor result = 2
    0, 1, 0, 1,
    // xor result = 3
    1, 2, 1, 2,
];

// ---------------------------------------------------------------------------
// Concrete get_values_from_key functions for normalization
// ---------------------------------------------------------------------------

/// Get normalization values for witness extension (base=16).
pub fn get_witness_normalization_values(key: [u64; 2]) -> [Fr; 2] {
    sparse::get_sparse_normalization_values(16, &WITNESS_EXTENSION_NORMALIZATION_TABLE, key)
}

/// Get normalization values for choose (base=28).
pub fn get_choose_normalization_values(key: [u64; 2]) -> [Fr; 2] {
    sparse::get_sparse_normalization_values(28, &CHOOSE_NORMALIZATION_TABLE, key)
}

/// Get normalization values for majority (base=16).
pub fn get_majority_normalization_values(key: [u64; 2]) -> [Fr; 2] {
    sparse::get_sparse_normalization_values(16, &MAJORITY_NORMALIZATION_TABLE, key)
}

// ---------------------------------------------------------------------------
// Generate normalization BasicTables
// ---------------------------------------------------------------------------

/// Generate the witness extension normalization table.
///
/// Mirrors C++ `generate_witness_extension_normalization_table`.
pub fn generate_witness_extension_normalization_table(
    id: BasicTableId,
    table_index: usize,
) -> BasicTable {
    sparse::generate_sparse_normalization_table(
        16,
        3,
        &WITNESS_EXTENSION_NORMALIZATION_TABLE,
        id,
        table_index,
        get_witness_normalization_values,
    )
}

/// Generate the choose normalization table.
///
/// Mirrors C++ `generate_choose_normalization_table`.
pub fn generate_choose_normalization_table(
    id: BasicTableId,
    table_index: usize,
) -> BasicTable {
    sparse::generate_sparse_normalization_table(
        28,
        2,
        &CHOOSE_NORMALIZATION_TABLE,
        id,
        table_index,
        get_choose_normalization_values,
    )
}

/// Generate the majority normalization table.
///
/// Mirrors C++ `generate_majority_normalization_table`.
pub fn generate_majority_normalization_table(
    id: BasicTableId,
    table_index: usize,
) -> BasicTable {
    sparse::generate_sparse_normalization_table(
        16,
        3,
        &MAJORITY_NORMALIZATION_TABLE,
        id,
        table_index,
        get_majority_normalization_values,
    )
}

// ---------------------------------------------------------------------------
// Rotation multiplier helpers
// ---------------------------------------------------------------------------

/// Compute the majority rotation multipliers for SHA-256.
///
/// Rotations: 2, 13, 22 using base 16 sparse representation.
/// Returns the multipliers used to combine sparse limbs with rotation
/// coefficients in the multi-table accumulator.
///
/// Mirrors C++ `get_majority_rotation_multipliers()`.
pub fn get_majority_rotation_multipliers() -> [Fr; 3] {
    let base = Fr::from(16u64);

    // Scaling factors applied to a's sparse limbs, excluding the rotated limb.
    let rot2_coefficients = [Fr::zero(), fr_pow(&base, 11 - 2), fr_pow(&base, 22 - 2)];
    let rot13_coefficients = [fr_pow(&base, 32 - 13), Fr::zero(), fr_pow(&base, 22 - 13)];
    let rot22_coefficients = [fr_pow(&base, 32 - 22), fr_pow(&base, 32 - 22 + 11), Fr::zero()];

    let target_rotation_coefficients = [
        rot2_coefficients[0] + rot13_coefficients[0] + rot22_coefficients[0],
        rot2_coefficients[1] + rot13_coefficients[1] + rot22_coefficients[1],
        rot2_coefficients[2] + rot13_coefficients[2] + rot22_coefficients[2],
    ];

    let column_2_row_1_multiplier = target_rotation_coefficients[0];
    let column_2_row_2_multiplier =
        target_rotation_coefficients[0] * (-fr_pow(&base, 11)) + target_rotation_coefficients[1];

    [column_2_row_1_multiplier, column_2_row_2_multiplier, Fr::zero()]
}

/// Compute the choose rotation multipliers for SHA-256.
///
/// Rotations: 6, 11, 25 using base 28 sparse representation.
///
/// Mirrors C++ `get_choose_rotation_multipliers()`.
pub fn get_choose_rotation_multipliers() -> [Fr; 3] {
    let base28 = Fr::from(28u64);

    let _column_2_row_3_coefficients = [
        Fr::one(),
        fr_pow(&base28, 11),
        fr_pow(&base28, 22),
    ];

    // Scaling factors applied to a's sparse limbs, excluding the rotated limb.
    let rot6_coefficients = [Fr::zero(), fr_pow(&base28, 11 - 6), fr_pow(&base28, 22 - 6)];
    let rot11_coefficients = [fr_pow(&base28, 32 - 11), Fr::zero(), fr_pow(&base28, 22 - 11)];
    let rot25_coefficients = [fr_pow(&base28, 32 - 25), fr_pow(&base28, 32 - 25 + 11), Fr::zero()];

    let target_rotation_coefficients = [
        rot6_coefficients[0] + rot11_coefficients[0] + rot25_coefficients[0],
        rot6_coefficients[1] + rot11_coefficients[1] + rot25_coefficients[1],
        rot6_coefficients[2] + rot11_coefficients[2] + rot25_coefficients[2],
    ];

    let column_2_row_1_multiplier = target_rotation_coefficients[0];

    // Compute the correct scaling factor for a0's 1st limb
    let current_coefficients = [
        Fr::one() * column_2_row_1_multiplier,
        fr_pow(&base28, 11) * column_2_row_1_multiplier,
        fr_pow(&base28, 22) * column_2_row_1_multiplier,
    ];

    let column_2_row_3_multiplier =
        -(current_coefficients[2]) + target_rotation_coefficients[2];

    [column_2_row_1_multiplier, Fr::zero(), column_2_row_3_multiplier]
}

// ---------------------------------------------------------------------------
// MultiTable constructors: output (normalization) tables
// ---------------------------------------------------------------------------

/// Create the witness extension output multi-table.
///
/// Mirrors C++ `get_witness_extension_output_table`.
pub fn get_witness_extension_output_table(id: MultiTableId) -> MultiTable {
    let num_entries = 11;

    let mut table = MultiTable::new_repeated(
        Fr::from(pow64(16, 3)),
        Fr::from(1u64 << 3),
        Fr::zero(),
        num_entries,
    );
    table.id = id;

    for _ in 0..num_entries {
        table.slice_sizes.push(pow64(16, 3));
        table.basic_table_ids.push(BasicTableId::SHA256_WITNESS_NORMALIZE);
        table.get_table_values.push(get_witness_normalization_values);
    }

    table
}

/// Create the choose output multi-table.
///
/// Mirrors C++ `get_choose_output_table`.
pub fn get_choose_output_table(id: MultiTableId) -> MultiTable {
    let num_entries = 16;

    let mut table = MultiTable::new_repeated(
        Fr::from(pow64(28, 2)),
        Fr::from(1u64 << 2),
        Fr::zero(),
        num_entries,
    );
    table.id = id;

    for _ in 0..num_entries {
        table.slice_sizes.push(pow64(28, 2));
        table.basic_table_ids.push(BasicTableId::SHA256_CH_NORMALIZE);
        table.get_table_values.push(get_choose_normalization_values);
    }

    table
}

/// Create the majority output multi-table.
///
/// Mirrors C++ `get_majority_output_table`.
pub fn get_majority_output_table(id: MultiTableId) -> MultiTable {
    let num_entries = 11;

    let mut table = MultiTable::new_repeated(
        Fr::from(pow64(16, 3)),
        Fr::from(1u64 << 3),
        Fr::zero(),
        num_entries,
    );
    table.id = id;

    for _ in 0..num_entries {
        table.slice_sizes.push(pow64(16, 3));
        table.basic_table_ids.push(BasicTableId::SHA256_MAJ_NORMALIZE);
        table.get_table_values.push(get_majority_normalization_values);
    }

    table
}

// ---------------------------------------------------------------------------
// MultiTable constructors: input (sparse conversion) tables
// ---------------------------------------------------------------------------

/// Create the witness extension input multi-table.
///
/// Decomposes a 32-bit scalar into slices (3, 7, 8, 14 bits) and converts
/// each to sparse base-16 form with optional rotation.
///
/// Mirrors C++ `get_witness_extension_input_table`.
pub fn get_witness_extension_input_table(id: MultiTableId) -> MultiTable {
    let column_1_coefficients = vec![
        Fr::from(1u64),
        Fr::from(1u64 << 3),
        Fr::from(1u64 << 10),
        Fr::from(1u64 << 18),
    ];
    let column_2_coefficients = vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()];
    let column_3_coefficients = vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()];

    let mut table = MultiTable::new_explicit(
        column_1_coefficients,
        column_2_coefficients,
        column_3_coefficients,
    );
    table.id = id;
    table.slice_sizes = vec![1 << 3, 1 << 7, 1 << 8, 1 << 18];
    table.basic_table_ids = vec![
        BasicTableId::SHA256_WITNESS_SLICE_3,
        BasicTableId::SHA256_WITNESS_SLICE_7_ROTATE_4,
        BasicTableId::SHA256_WITNESS_SLICE_8_ROTATE_7,
        BasicTableId::SHA256_WITNESS_SLICE_14_ROTATE_1,
    ];
    table.get_table_values = vec![
        sparse::get_sparse_values_base16_rot0,
        sparse::get_sparse_values_base16_rot4,
        sparse::get_sparse_values_base16_rot7,
        sparse::get_sparse_values_base16_rot1,
    ];

    table
}

/// Create the choose input multi-table.
///
/// Decomposes a 32-bit scalar into 3 slices (11, 11, 10 bits) and converts
/// each to sparse base-28 form. Column 3 contains rotation combination terms.
///
/// Mirrors C++ `get_choose_input_table`.
pub fn get_choose_input_table(id: MultiTableId) -> MultiTable {
    let base28 = Fr::from(28u64);

    // Scaling factors for rotations 6, 11, 25
    let rot6_coefficients = [Fr::zero(), fr_pow(&base28, 11 - 6), fr_pow(&base28, 22 - 6)];
    let rot11_coefficients = [fr_pow(&base28, 32 - 11), Fr::zero(), fr_pow(&base28, 22 - 11)];
    let rot25_coefficients = [fr_pow(&base28, 32 - 25), fr_pow(&base28, 32 - 25 + 11), Fr::zero()];

    let target_rotation_coefficients = [
        rot6_coefficients[0] + rot11_coefficients[0] + rot25_coefficients[0],
        rot6_coefficients[1] + rot11_coefficients[1] + rot25_coefficients[1],
        rot6_coefficients[2] + rot11_coefficients[2] + rot25_coefficients[2],
    ];

    let column_2_row_1_multiplier = target_rotation_coefficients[0];

    let current_coefficients = [
        column_2_row_1_multiplier,
        fr_pow(&base28, 11) * column_2_row_1_multiplier,
        fr_pow(&base28, 22) * column_2_row_1_multiplier,
    ];

    let column_3_row_2_multiplier =
        -(current_coefficients[1]) + target_rotation_coefficients[1];

    let column_1_coefficients = vec![
        Fr::from(1u64),
        Fr::from(1u64 << 11),
        Fr::from(1u64 << 22),
    ];
    let column_2_coefficients = vec![Fr::one(), fr_pow(&base28, 11), fr_pow(&base28, 22)];
    let column_3_coefficients = vec![
        Fr::one(),
        column_3_row_2_multiplier + Fr::one(),
        Fr::one(),
    ];

    let mut table = MultiTable::new_explicit(
        column_1_coefficients,
        column_2_coefficients,
        column_3_coefficients,
    );
    table.id = id;
    table.slice_sizes = vec![1 << 11, 1 << 11, 1 << 10];
    table.basic_table_ids = vec![
        BasicTableId::SHA256_BASE28_ROTATE6,
        BasicTableId::SHA256_BASE28,
        BasicTableId::SHA256_BASE28_ROTATE3,
    ];
    table.get_table_values = vec![
        sparse::get_sparse_values_base28_rot6,
        sparse::get_sparse_values_base28_rot0,
        sparse::get_sparse_values_base28_rot3,
    ];

    table
}

/// Create the majority input multi-table.
///
/// Decomposes a 32-bit scalar into 3 slices (11, 11, 10 bits) and converts
/// each to sparse base-16 form. Column 3 contains rotation combination terms
/// for rotations 2, 13, 22.
///
/// Mirrors C++ `get_majority_input_table`.
pub fn get_majority_input_table(id: MultiTableId) -> MultiTable {
    let base = Fr::from(16u64);

    // Scaling factors for rotations 2, 13, 22
    let rot2_coefficients = [Fr::zero(), fr_pow(&base, 11 - 2), fr_pow(&base, 22 - 2)];
    let rot13_coefficients = [fr_pow(&base, 32 - 13), Fr::zero(), fr_pow(&base, 22 - 13)];
    let rot22_coefficients = [fr_pow(&base, 32 - 22), fr_pow(&base, 32 - 22 + 11), Fr::zero()];

    let target_rotation_coefficients = [
        rot2_coefficients[0] + rot13_coefficients[0] + rot22_coefficients[0],
        rot2_coefficients[1] + rot13_coefficients[1] + rot22_coefficients[1],
        rot2_coefficients[2] + rot13_coefficients[2] + rot22_coefficients[2],
    ];

    let column_2_row_3_multiplier =
        target_rotation_coefficients[1] * (-fr_pow(&base, 11)) + target_rotation_coefficients[2];

    let column_1_coefficients = vec![
        Fr::from(1u64),
        Fr::from(1u64 << 11),
        Fr::from(1u64 << 22),
    ];
    let column_2_coefficients = vec![Fr::one(), fr_pow(&base, 11), fr_pow(&base, 22)];
    let column_3_coefficients = vec![
        Fr::one(),
        Fr::one(),
        Fr::one() + column_2_row_3_multiplier,
    ];

    let mut table = MultiTable::new_explicit(
        column_1_coefficients,
        column_2_coefficients,
        column_3_coefficients,
    );
    table.id = id;
    table.slice_sizes = vec![1 << 11, 1 << 11, 1 << 10];
    table.basic_table_ids = vec![
        BasicTableId::SHA256_BASE16_ROTATE2,
        BasicTableId::SHA256_BASE16_ROTATE2,
        BasicTableId::SHA256_BASE16,
    ];
    table.get_table_values = vec![
        sparse::get_sparse_values_base16_rot2,
        sparse::get_sparse_values_base16_rot2,
        sparse::get_sparse_values_base16_rot0,
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
    fn test_normalization_table_sizes() {
        assert_eq!(CHOOSE_NORMALIZATION_TABLE.len(), 28);
        assert_eq!(MAJORITY_NORMALIZATION_TABLE.len(), 16);
        assert_eq!(WITNESS_EXTENSION_NORMALIZATION_TABLE.len(), 16);
    }

    #[test]
    fn test_generate_choose_normalization() {
        let table = generate_choose_normalization_table(BasicTableId::SHA256_CH_NORMALIZE, 0);
        // 28^2 = 784 entries
        assert_eq!(table.column_1.len(), pow64(28, 2) as usize);
        assert!(!table.use_twin_keys);
    }

    #[test]
    fn test_generate_majority_normalization() {
        let table = generate_majority_normalization_table(BasicTableId::SHA256_MAJ_NORMALIZE, 0);
        // 16^3 = 4096 entries
        assert_eq!(table.column_1.len(), pow64(16, 3) as usize);
    }

    #[test]
    fn test_generate_witness_extension_normalization() {
        let table = generate_witness_extension_normalization_table(
            BasicTableId::SHA256_WITNESS_NORMALIZE,
            0,
        );
        // 16^3 = 4096 entries
        assert_eq!(table.column_1.len(), pow64(16, 3) as usize);
    }

    #[test]
    fn test_choose_output_table() {
        let table = get_choose_output_table(MultiTableId::SHA256_CH_OUTPUT);
        assert_eq!(table.id, MultiTableId::SHA256_CH_OUTPUT);
        assert_eq!(table.slice_sizes.len(), 16);
        assert_eq!(table.basic_table_ids.len(), 16);
        assert_eq!(table.get_table_values.len(), 16);
    }

    #[test]
    fn test_majority_output_table() {
        let table = get_majority_output_table(MultiTableId::SHA256_MAJ_OUTPUT);
        assert_eq!(table.id, MultiTableId::SHA256_MAJ_OUTPUT);
        assert_eq!(table.slice_sizes.len(), 11);
    }

    #[test]
    fn test_witness_extension_output_table() {
        let table = get_witness_extension_output_table(MultiTableId::SHA256_WITNESS_OUTPUT);
        assert_eq!(table.id, MultiTableId::SHA256_WITNESS_OUTPUT);
        assert_eq!(table.slice_sizes.len(), 11);
    }

    #[test]
    fn test_choose_input_table() {
        let table = get_choose_input_table(MultiTableId::SHA256_CH_INPUT);
        assert_eq!(table.id, MultiTableId::SHA256_CH_INPUT);
        assert_eq!(table.slice_sizes.len(), 3);
        assert_eq!(table.slice_sizes, vec![1 << 11, 1 << 11, 1 << 10]);
        assert_eq!(table.basic_table_ids.len(), 3);
        assert_eq!(table.get_table_values.len(), 3);
    }

    #[test]
    fn test_majority_input_table() {
        let table = get_majority_input_table(MultiTableId::SHA256_MAJ_INPUT);
        assert_eq!(table.id, MultiTableId::SHA256_MAJ_INPUT);
        assert_eq!(table.slice_sizes.len(), 3);
        assert_eq!(table.slice_sizes, vec![1 << 11, 1 << 11, 1 << 10]);
    }

    #[test]
    fn test_witness_extension_input_table() {
        let table = get_witness_extension_input_table(MultiTableId::SHA256_WITNESS_INPUT);
        assert_eq!(table.id, MultiTableId::SHA256_WITNESS_INPUT);
        assert_eq!(table.slice_sizes.len(), 4);
        assert_eq!(table.slice_sizes, vec![1 << 3, 1 << 7, 1 << 8, 1 << 18]);
    }

    #[test]
    fn test_rotation_multipliers_length() {
        let maj_rot = get_majority_rotation_multipliers();
        assert_eq!(maj_rot.len(), 3);
        // Third element should be zero
        assert_eq!(maj_rot[2], Fr::zero());

        let ch_rot = get_choose_rotation_multipliers();
        assert_eq!(ch_rot.len(), 3);
        // Second element should be zero
        assert_eq!(ch_rot[1], Fr::zero());
    }
}
