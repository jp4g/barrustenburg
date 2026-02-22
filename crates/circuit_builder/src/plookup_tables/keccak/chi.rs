//! Keccak CHI round plookup tables.
//!
//! C++ source: plookup_tables/keccak/keccak_chi.hpp
//!
//! The CHI round performs A ^ (~B & C) on 3 hash lanes.
//! In the base-11 sparse representation this becomes the linear expression
//! 2*A - B + C + Q where Q = sum_{i=0}^{63} 11^i.
//!
//! The normalization table maps the algebraic output back to binary:
//!   0 -> 0, 1 -> 0, 2 -> 1, 3 -> 1, 4 -> 0
//!
//! The table also extracts the MSB of the normalized output for the most
//! significant slice.

use bbrs_ecc::curves::bn254::Fr;
use bbrs_numeric::bitop::pow64;

use crate::plookup_tables::types::{BasicTable, BasicTableId, MultiTable, MultiTableId};

/// Base for sparse representation.
const BASE: u64 = 11;

/// Effective base: maximum value each quasi-bit can reach at this stage.
const EFFECTIVE_BASE: u64 = 5;

/// Number of quasi-bits per lookup table slice.
const TABLE_BITS: usize = 6;

/// CHI normalization: maps algebraic output to binary.
const CHI_NORMALIZATION_TABLE: [u64; 5] = [0, 0, 1, 1, 0];

/// Precompute array of base multipliers (11^i for i = 0..TABLE_BITS-1).
fn get_scaled_bases() -> [u64; TABLE_BITS] {
    let mut result = [0u64; TABLE_BITS];
    let mut acc = 1u64;
    for i in 0..TABLE_BITS {
        result[i] = acc;
        acc *= BASE;
    }
    result
}

/// Increment `counts` (array of quasi-bits in [0, EFFECTIVE_BASE-1]) and
/// return the (value, normalized_value, msb) columns for the new row.
fn get_column_values_for_next_row(counts: &mut [usize; TABLE_BITS]) -> [u64; 3] {
    let scaled_bases = get_scaled_bases();

    let divisor_exponent: u64 = if 64 % TABLE_BITS == 0 {
        (TABLE_BITS as u64) - 1
    } else {
        (64 % TABLE_BITS) as u64 - 1
    };
    let divisor = pow64(BASE, divisor_exponent);

    for i in 0..TABLE_BITS {
        if counts[i] == (EFFECTIVE_BASE as usize) - 1 {
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
        normalized_value += CHI_NORMALIZATION_TABLE[counts[i]] * scaled_bases[i];
    }
    [value, normalized_value, normalized_value / divisor]
}

/// Given a table input value, return the normalized output and MSB extraction.
///
/// Decomposes the input in base-11 and applies the CHI normalization table.
/// Returns [normalized_accumulator, normalized_accumulator / divisor].
pub fn get_chi_values(key: [u64; 2]) -> [Fr; 2] {
    let mut accumulator = 0u64;
    let mut input = key[0];
    let mut base_shift = 1u64;
    let divisor_exponent: u64 = if 64 % TABLE_BITS == 0 {
        (TABLE_BITS as u64) - 1
    } else {
        (64 % TABLE_BITS) as u64 - 1
    };
    let divisor = pow64(BASE, divisor_exponent);

    while input > 0 {
        let slice = input % BASE;
        let bit = CHI_NORMALIZATION_TABLE[slice as usize];
        accumulator += bit * base_shift;
        input /= BASE;
        base_shift *= BASE;
    }

    [Fr::from(accumulator), Fr::from(accumulator / divisor)]
}

/// Generate the CHI normalization BasicTable.
///
/// Table size = EFFECTIVE_BASE^TABLE_BITS. Each row maps a base-11 quasi-bit
/// pattern to its normalized form and extracts the MSB.
pub fn generate_chi_renormalization_table(id: BasicTableId, table_index: usize) -> BasicTable {
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    let table_size = pow64(EFFECTIVE_BASE, TABLE_BITS as u64);

    let mut counts = [0usize; TABLE_BITS];
    let mut column_values = [0u64; 3];

    for _ in 0..table_size {
        table.column_1.push(Fr::from(column_values[0]));
        table.column_2.push(Fr::from(column_values[1]));
        table.column_3.push(Fr::from(column_values[2]));
        column_values = get_column_values_for_next_row(&mut counts);
    }

    table.get_values_from_key = get_chi_values;

    let step_size = pow64(BASE, TABLE_BITS as u64);
    table.column_1_step_size = Fr::from(step_size);
    table.column_2_step_size = Fr::from(step_size);
    table.column_3_step_size = Fr::zero();

    table
}

/// Create the CHI MultiTable used by plookup for a sequence of lookups.
///
/// 64-bit integers require ceil(64 / TABLE_BITS) lookups.
pub fn get_chi_output_table(id: MultiTableId) -> MultiTable {
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
        table.basic_table_ids.push(BasicTableId::KECCAK_CHI);
        table.get_table_values.push(get_chi_values);
    }

    table
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chi_normalization() {
        // 0 -> 0
        let r = get_chi_values([0, 0]);
        assert_eq!(r[0], Fr::zero());

        // 1 -> 0
        let r = get_chi_values([1, 0]);
        assert_eq!(r[0], Fr::zero());

        // 2 -> 1 (in base-11 position 0, so value = 1)
        let r = get_chi_values([2, 0]);
        assert_eq!(r[0], Fr::from(1u64));

        // 3 -> 1
        let r = get_chi_values([3, 0]);
        assert_eq!(r[0], Fr::from(1u64));

        // 4 -> 0
        let r = get_chi_values([4, 0]);
        assert_eq!(r[0], Fr::zero());
    }

    #[test]
    fn test_generate_chi_table() {
        let table = generate_chi_renormalization_table(BasicTableId::KECCAK_CHI, 0);
        // table_size = 5^6 = 15625
        assert_eq!(table.column_1.len(), 15625);
        assert_eq!(table.column_1[0], Fr::zero());
        assert_eq!(table.column_2[0], Fr::zero());
        assert_eq!(table.column_3[0], Fr::zero());
    }

    #[test]
    fn test_get_chi_output_table() {
        let multi = get_chi_output_table(MultiTableId::KECCAK_CHI_OUTPUT);
        assert_eq!(multi.id, MultiTableId::KECCAK_CHI_OUTPUT);
        // 64 / 6 = 10, 64 % 6 = 4 != 0, so 11
        assert_eq!(multi.slice_sizes.len(), 11);
    }
}
