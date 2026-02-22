//! Keccak THETA round plookup tables.
//!
//! C++ source: plookup_tables/keccak/keccak_theta.hpp
//!
//! The THETA round performs:
//!   C0 = A0 ^ A1 ^ A2 ^ A3 ^ A4
//!   C1 = B0 ^ B1 ^ B2 ^ B3 ^ B4
//!   D  = C0 ^ ROTATE_LEFT(C1, 1)
//!
//! In the base-11 sparse representation, XOR becomes addition.
//! Each quasi-bit of D can reach a maximum value of 10. The THETA
//! normalization table maps each quasi-bit back to binary (even=0, odd=1).

use bbrs_ecc::curves::bn254::Fr;
use bbrs_numeric::bitop::pow64;

use crate::plookup_tables::types::{BasicTable, BasicTableId, MultiTable, MultiTableId};

/// Base for sparse representation.
const BASE: u64 = 11;

/// Number of quasi-bits per lookup table slice.
const TABLE_BITS: usize = 4;

/// Normalization table: even values -> 0, odd values -> 1.
const THETA_NORMALIZATION_TABLE: [u64; 11] = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0];

/// Precompute an array of base multipliers (11^i for i = 0..TABLE_BITS-1).
fn get_scaled_bases() -> [u64; TABLE_BITS] {
    let mut result = [0u64; TABLE_BITS];
    let mut acc = 1u64;
    for i in 0..TABLE_BITS {
        result[i] = acc;
        acc *= BASE;
    }
    result
}

/// Increment `counts` (array of quasi-bits in [0, BASE-1]) and return
/// the (value, normalized_value) column pair for the new row.
fn get_column_values_for_next_row(counts: &mut [usize; TABLE_BITS]) -> [u64; 2] {
    let scaled_bases = get_scaled_bases();

    for i in 0..TABLE_BITS {
        if counts[i] == (BASE as usize) - 1 {
            counts[i] = 0;
        } else {
            counts[i] += 1;
            break;
        }
    }

    let mut value = 0u64;
    let mut normalized_value = 0u64;
    for i in 0..TABLE_BITS {
        value += (counts[i] as u64) * scaled_bases[i];
        normalized_value += THETA_NORMALIZATION_TABLE[counts[i]] * scaled_bases[i];
    }
    [value, normalized_value]
}

/// Given a table input value, return the normalized output.
///
/// Decomposes the input in base-11 and applies the normalization table
/// to each quasi-bit. Returns [normalized_accumulator, 0].
pub fn get_theta_values(key: [u64; 2]) -> [Fr; 2] {
    let mut accumulator = 0u64;
    let mut input = key[0];
    let mut base_shift = 1u64;
    while input > 0 {
        let slice = input % BASE;
        let bit = THETA_NORMALIZATION_TABLE[slice as usize];
        accumulator += bit * base_shift;
        input /= BASE;
        base_shift *= BASE;
    }
    [Fr::from(accumulator), Fr::zero()]
}

/// Generate the THETA normalization BasicTable.
///
/// Table size = BASE^TABLE_BITS. Each row maps a base-11 quasi-bit pattern
/// to its normalized form.
pub fn generate_theta_renormalization_table(id: BasicTableId, table_index: usize) -> BasicTable {
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    let table_size = pow64(BASE, TABLE_BITS as u64);

    let mut counts = [0usize; TABLE_BITS];
    let mut column_values = [0u64; 2];

    for _ in 0..table_size {
        table.column_1.push(Fr::from(column_values[0]));
        table.column_2.push(Fr::from(column_values[1]));
        table.column_3.push(Fr::zero());
        column_values = get_column_values_for_next_row(&mut counts);
    }

    table.get_values_from_key = get_theta_values;

    let step_size = pow64(BASE, TABLE_BITS as u64);
    table.column_1_step_size = Fr::from(step_size);
    table.column_2_step_size = Fr::from(step_size);
    table.column_3_step_size = Fr::zero();

    table
}

/// Create the THETA MultiTable used by plookup for a sequence of lookups.
///
/// 64-bit integers require ceil(64 / TABLE_BITS) lookups.
pub fn get_theta_output_table(id: MultiTableId) -> MultiTable {
    let num_tables: usize = (64 / TABLE_BITS) + if 64 % TABLE_BITS == 0 { 0 } else { 1 };

    let column_multiplier = pow64(BASE, TABLE_BITS as u64);
    let mut table = MultiTable::new_repeated(
        Fr::from(column_multiplier),
        Fr::from(column_multiplier),
        Fr::zero(),
        num_tables,
    );

    table.id = id;
    for _ in 0..num_tables {
        table.slice_sizes.push(pow64(BASE, TABLE_BITS as u64));
        table.basic_table_ids.push(BasicTableId::KECCAK_THETA);
        table.get_table_values.push(get_theta_values);
    }

    table
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_theta_normalization() {
        // 0 -> 0
        let r = get_theta_values([0, 0]);
        assert_eq!(r[0], Fr::zero());

        // 1 -> 1 (odd)
        let r = get_theta_values([1, 0]);
        assert_eq!(r[0], Fr::from(1u64));

        // 2 -> 0 (even)
        let r = get_theta_values([2, 0]);
        assert_eq!(r[0], Fr::zero());

        // 3 -> 1 (odd)
        let r = get_theta_values([3, 0]);
        assert_eq!(r[0], Fr::from(1u64));
    }

    #[test]
    fn test_generate_theta_table() {
        let table =
            generate_theta_renormalization_table(BasicTableId::KECCAK_THETA, 0);
        // table_size = 11^4 = 14641
        assert_eq!(table.column_1.len(), 14641);
        assert_eq!(table.column_1[0], Fr::zero());
        assert_eq!(table.column_2[0], Fr::zero());
    }

    #[test]
    fn test_get_theta_output_table() {
        let multi = get_theta_output_table(MultiTableId::KECCAK_THETA_OUTPUT);
        assert_eq!(multi.id, MultiTableId::KECCAK_THETA_OUTPUT);
        // 64 / 4 = 16
        assert_eq!(multi.slice_sizes.len(), 16);
    }
}
