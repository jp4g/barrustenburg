//! Keccak RHO round plookup tables.
//!
//! C++ source: plookup_tables/keccak/keccak_rho.hpp
//!
//! The RHO round performs left-rotation on each of the 25 hash lanes by a
//! fixed rotation value. At this point in the algorithm each quasi-bit is
//! in the range [0, 1, 2].
//!
//! The RHO lookup tables:
//!   1. Normalize each quasi-bit so that P_out = sum (b_i mod 2) * 11^i
//!   2. Perform left-rotation by lane-specific offsets
//!   3. Extract the MSB of the non-rotated normalized output
//!
//! In C++ this is templated on TABLE_BITS and LANE_INDEX. In Rust we use
//! runtime parameters and concrete function pointers for TABLE_BITS 1..8.

use bbrs_ecc::curves::bn254::Fr;
use bbrs_numeric::bitop::pow64;

use crate::plookup_tables::types::{BasicTable, BasicTableId, MultiTable, MultiTableId};

/// Base for sparse representation.
const BASE: u64 = 11;

/// Effective base: maximum quasi-bit value at this stage.
const EFFECTIVE_BASE: u64 = 3;

/// Maximum bits per lookup sub-table.
const MAXIMUM_MULTITABLE_BITS: usize = 8;

/// Rotation offsets per lane, indexed by y * 5 + x.
const ROTATIONS: [usize; 25] = [
    0, 1, 62, 28, 27, 36, 44, 6, 55, 20, 3, 10, 43, 25, 39, 41, 45, 15, 21, 8, 18, 2, 61, 56,
    14,
];

/// Normalization table: 0->0, 1->1, 2->0.
const RHO_NORMALIZATION_TABLE: [u64; 3] = [0, 1, 0];

/// Precompute base multipliers for a given table_bits.
fn get_scaled_bases(table_bits: usize) -> Vec<u64> {
    let mut result = Vec::with_capacity(table_bits);
    let mut acc = 1u64;
    for _ in 0..table_bits {
        result.push(acc);
        acc *= BASE;
    }
    result
}

/// Increment `counts` (quasi-bits in [0, EFFECTIVE_BASE-1]) for a given
/// table_bits and return (value, normalized_value, msb_extraction).
fn get_column_values_for_next_row(
    counts: &mut Vec<usize>,
    table_bits: usize,
) -> [u64; 3] {
    let scaled_bases = get_scaled_bases(table_bits);
    let divisor = pow64(BASE, (table_bits as u64) - 1);

    for i in 0..table_bits {
        if counts[i] == (EFFECTIVE_BASE as usize) - 1 {
            counts[i] = 0;
        } else {
            counts[i] += 1;
            break;
        }
    }

    let mut value = 0u64;
    let mut normalized_value = 0u64;
    for i in 0..table_bits {
        value += (counts[i] as u64) * scaled_bases[i];
        normalized_value += RHO_NORMALIZATION_TABLE[counts[i]] * scaled_bases[i];
    }
    [value, normalized_value, normalized_value / divisor]
}

/// Core renormalization logic parameterized by table_bits.
fn get_rho_renormalization_values_impl(key: [u64; 2], table_bits: u64) -> [Fr; 2] {
    let mut accumulator = 0u64;
    let mut input = key[0];
    let mut base_shift = 1u64;
    let divisor_exponent = table_bits - 1;
    let divisor = pow64(BASE, divisor_exponent);

    while input > 0 {
        let slice = input % BASE;
        let bit = RHO_NORMALIZATION_TABLE[slice as usize];
        accumulator += bit * base_shift;
        input /= BASE;
        base_shift *= BASE;
    }

    [Fr::from(accumulator), Fr::from(accumulator / divisor)]
}

// Concrete function pointers for TABLE_BITS = 1..8.
// Required because BasicTable.get_values_from_key is a fn pointer, not a closure.

/// Rho renormalization for TABLE_BITS = 1.
pub fn get_rho_values_1(key: [u64; 2]) -> [Fr; 2] {
    get_rho_renormalization_values_impl(key, 1)
}

/// Rho renormalization for TABLE_BITS = 2.
pub fn get_rho_values_2(key: [u64; 2]) -> [Fr; 2] {
    get_rho_renormalization_values_impl(key, 2)
}

/// Rho renormalization for TABLE_BITS = 3.
pub fn get_rho_values_3(key: [u64; 2]) -> [Fr; 2] {
    get_rho_renormalization_values_impl(key, 3)
}

/// Rho renormalization for TABLE_BITS = 4.
pub fn get_rho_values_4(key: [u64; 2]) -> [Fr; 2] {
    get_rho_renormalization_values_impl(key, 4)
}

/// Rho renormalization for TABLE_BITS = 5.
pub fn get_rho_values_5(key: [u64; 2]) -> [Fr; 2] {
    get_rho_renormalization_values_impl(key, 5)
}

/// Rho renormalization for TABLE_BITS = 6.
pub fn get_rho_values_6(key: [u64; 2]) -> [Fr; 2] {
    get_rho_renormalization_values_impl(key, 6)
}

/// Rho renormalization for TABLE_BITS = 7.
pub fn get_rho_values_7(key: [u64; 2]) -> [Fr; 2] {
    get_rho_renormalization_values_impl(key, 7)
}

/// Rho renormalization for TABLE_BITS = 8.
pub fn get_rho_values_8(key: [u64; 2]) -> [Fr; 2] {
    get_rho_renormalization_values_impl(key, 8)
}

/// Look up the concrete function pointer for a given table_bits (1..=8).
fn get_rho_fn(table_bits: usize) -> fn([u64; 2]) -> [Fr; 2] {
    match table_bits {
        1 => get_rho_values_1,
        2 => get_rho_values_2,
        3 => get_rho_values_3,
        4 => get_rho_values_4,
        5 => get_rho_values_5,
        6 => get_rho_values_6,
        7 => get_rho_values_7,
        8 => get_rho_values_8,
        _ => panic!("Rho TABLE_BITS must be in range 1..=8, got {}", table_bits),
    }
}

/// Generate a Rho renormalization BasicTable for a given `table_bits`.
///
/// Table size = EFFECTIVE_BASE^table_bits.
pub fn generate_rho_renormalization_table(
    id: BasicTableId,
    table_index: usize,
    table_bits: usize,
) -> BasicTable {
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    let table_size = pow64(EFFECTIVE_BASE, table_bits as u64);

    let mut counts = vec![0usize; table_bits];
    let mut column_values = [0u64; 3];

    for _ in 0..table_size {
        table.column_1.push(Fr::from(column_values[0]));
        table.column_2.push(Fr::from(column_values[1]));
        table.column_3.push(Fr::from(column_values[2]));
        column_values = get_column_values_for_next_row(&mut counts, table_bits);
    }

    table.get_values_from_key = get_rho_fn(table_bits);

    let step_size = pow64(BASE, table_bits as u64);
    table.column_1_step_size = Fr::from(step_size);
    table.column_2_step_size = Fr::from(step_size);
    table.column_3_step_size = Fr::zero();

    table
}

/// Create the Rho MultiTable for a given lane index.
///
/// Rotations are performed by splitting the input into left and right bit
/// slices. Both slices are fed into sub-tables of varying sizes to correctly
/// range-constrain the input.
///
/// `lane_index` selects the rotation value from the ROTATIONS array (0..24).
pub fn get_rho_output_table(id: MultiTableId, lane_index: usize) -> MultiTable {
    let left_bits = ROTATIONS[lane_index];
    let right_bits = 64 - ROTATIONS[lane_index];

    let num_right_tables =
        right_bits / MAXIMUM_MULTITABLE_BITS + if right_bits % MAXIMUM_MULTITABLE_BITS > 0 { 1 } else { 0 };
    let num_left_tables =
        left_bits / MAXIMUM_MULTITABLE_BITS + if left_bits % MAXIMUM_MULTITABLE_BITS > 0 { 1 } else { 0 };

    let mut table = MultiTable::empty();
    table.id = id;

    table.column_1_step_sizes.push(Fr::one());
    table.column_2_step_sizes.push(Fr::one());
    table.column_3_step_sizes.push(Fr::one());

    // Generate table selector values for the 'right' slice.
    for i in 0..num_right_tables {
        let num_bits_processed = i * MAXIMUM_MULTITABLE_BITS;
        let bit_slice = if num_bits_processed + MAXIMUM_MULTITABLE_BITS > right_bits {
            right_bits % MAXIMUM_MULTITABLE_BITS
        } else {
            MAXIMUM_MULTITABLE_BITS
        };

        let scaled_base = pow64(BASE, bit_slice as u64);

        if i == num_right_tables - 1 {
            table.column_1_step_sizes.push(Fr::from(scaled_base));
            table.column_2_step_sizes.push(Fr::zero());
            table.column_3_step_sizes.push(Fr::zero());
        } else {
            table.column_1_step_sizes.push(Fr::from(scaled_base));
            table.column_2_step_sizes.push(Fr::from(scaled_base));
            table.column_3_step_sizes.push(Fr::zero());
        }

        table.slice_sizes.push(scaled_base);
        table.get_table_values.push(get_rho_fn(bit_slice));
        // BasicTableId for KECCAK_RHO_1 + (bit_slice - 1)
        table
            .basic_table_ids
            .push(BasicTableId(BasicTableId::KECCAK_RHO_1.0 + (bit_slice - 1)));
    }

    // Generate table selector values for the 'left' slice.
    for i in 0..num_left_tables {
        let num_bits_processed = i * MAXIMUM_MULTITABLE_BITS;
        let bit_slice = if num_bits_processed + MAXIMUM_MULTITABLE_BITS > left_bits {
            left_bits % MAXIMUM_MULTITABLE_BITS
        } else {
            MAXIMUM_MULTITABLE_BITS
        };

        let scaled_base = pow64(BASE, bit_slice as u64);

        if i != num_left_tables - 1 {
            table.column_1_step_sizes.push(Fr::from(scaled_base));
            table.column_2_step_sizes.push(Fr::from(scaled_base));
            table.column_3_step_sizes.push(Fr::zero());
        }

        table.slice_sizes.push(scaled_base);
        table.get_table_values.push(get_rho_fn(bit_slice));
        table
            .basic_table_ids
            .push(BasicTableId(BasicTableId::KECCAK_RHO_1.0 + (bit_slice - 1)));
    }

    // Compute coefficients from step sizes.
    let num_entries = table.column_1_step_sizes.len();
    table.column_1_coefficients = Vec::with_capacity(num_entries);
    table.column_2_coefficients = Vec::with_capacity(num_entries);
    table.column_3_coefficients = Vec::with_capacity(num_entries);

    let mut c1_acc = Fr::one();
    let mut c2_acc = Fr::one();
    let mut c3_acc = Fr::one();
    for i in 0..num_entries {
        table.column_1_coefficients.push(c1_acc);
        table.column_2_coefficients.push(c2_acc);
        table.column_3_coefficients.push(c3_acc);
        c1_acc = c1_acc * table.column_1_step_sizes[i];
        c2_acc = c2_acc * table.column_2_step_sizes[i];
        c3_acc = c3_acc * table.column_3_step_sizes[i];
    }

    table
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rho_normalization() {
        // 0 -> 0
        let r = get_rho_values_8([0, 0]);
        assert_eq!(r[0], Fr::zero());

        // 1 -> 1
        let r = get_rho_values_8([1, 0]);
        assert_eq!(r[0], Fr::from(1u64));

        // 2 -> 0
        let r = get_rho_values_8([2, 0]);
        assert_eq!(r[0], Fr::zero());
    }

    #[test]
    fn test_generate_rho_table() {
        // TABLE_BITS = 8 -> table_size = 3^8 = 6561
        let table = generate_rho_renormalization_table(
            BasicTableId::KECCAK_RHO_8,
            0,
            8,
        );
        assert_eq!(table.column_1.len(), 6561);
    }

    #[test]
    fn test_rho_output_table_lane0() {
        // Lane 0 has rotation = 0 -> right_bits = 64, left_bits = 0
        let multi = get_rho_output_table(MultiTableId::KECCAK_NORMALIZE_AND_ROTATE, 0);
        // right_bits=64, 64/8=8 tables, left_bits=0, 0 tables -> 8 total
        assert_eq!(multi.slice_sizes.len(), 8);
    }

    #[test]
    fn test_rho_output_table_lane1() {
        // Lane 1 has rotation = 1 -> right_bits = 63, left_bits = 1
        let multi = get_rho_output_table(
            MultiTableId(MultiTableId::KECCAK_NORMALIZE_AND_ROTATE.0 + 1),
            1,
        );
        // right_bits=63, ceil(63/8)=8, left_bits=1, ceil(1/8)=1 -> 9 total
        assert_eq!(multi.slice_sizes.len(), 9);
    }

    #[test]
    fn test_all_rho_fns() {
        // Ensure all 8 concrete fns work
        for bits in 1..=8 {
            let f = get_rho_fn(bits);
            let r = f([0, 0]);
            assert_eq!(r[0], Fr::zero());
        }
    }
}
