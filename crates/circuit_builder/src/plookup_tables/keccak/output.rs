//! Keccak output plookup tables.
//!
//! C++ source: plookup_tables/keccak/keccak_output.hpp
//!
//! Converts a base-11 sparse integer representation back into a regular
//! base-2 binary integer. Used by the Keccak hash algorithm to convert
//! the output of the algorithm into a regular integer.
//!
//! At this point in the algorithm each quasi-bit can only take values [0, 1],
//! so EFFECTIVE_BASE = 2 even though we use base-11 representation.
//! The output column is a proper binary integer (using bit shifts, not base-11).

use bbrs_ecc::curves::bn254::Fr;
use bbrs_numeric::bitop::pow64;

use crate::plookup_tables::types::{BasicTable, BasicTableId, MultiTable, MultiTableId};

/// Base for sparse representation.
const BASE: u64 = 11;

/// Effective base at output stage: quasi-bits can only be 0 or 1.
const EFFECTIVE_BASE: u64 = 2;

/// Number of quasi-bits per lookup table slice.
const TABLE_BITS: usize = 8;

/// Output normalization table (identity for binary values).
const OUTPUT_NORMALIZATION_TABLE: [u64; 2] = [0, 1];

/// Precompute base multipliers (11^i for i = 0..TABLE_BITS-1).
fn get_scaled_bases() -> [u64; TABLE_BITS] {
    let mut result = [0u64; TABLE_BITS];
    let mut acc = 1u64;
    for i in 0..TABLE_BITS {
        result[i] = acc;
        acc *= BASE;
    }
    result
}

/// Increment `counts` (quasi-bits in [0, EFFECTIVE_BASE-1]) and return
/// (sparse_value, binary_value) for the new row.
///
/// Note: the output column uses bit-shifting (base-2) rather than base-11
/// for the normalized value, since we're converting back to binary.
fn get_column_values_for_next_row(counts: &mut [usize; TABLE_BITS]) -> [u64; 2] {
    let scaled_bases = get_scaled_bases();

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
        normalized_value += OUTPUT_NORMALIZATION_TABLE[counts[i]] << i;
    }
    [value, normalized_value]
}

/// Given a table input (base-11 sparse integer), return the binary output.
///
/// Mirrors C++ `sparse_tables::get_sparse_normalization_values<BASE, OUTPUT_NORMALIZATION_TABLE>`.
/// Decomposes input in base-11 and maps each quasi-bit through the
/// normalization table, accumulating with bit-shifts (base-2 output).
///
/// Returns [binary_accumulator, 0].
pub fn get_keccak_output_normalization_values(key: [u64; 2]) -> [Fr; 2] {
    let mut accumulator = 0u64;
    let mut input = key[0];
    let mut count = 0u64;
    while input > 0 {
        let slice = input % BASE;
        let bit = OUTPUT_NORMALIZATION_TABLE[slice as usize];
        accumulator += bit << count;
        input -= slice;
        input /= BASE;
        count += 1;
    }
    [Fr::from(accumulator), Fr::zero()]
}

/// Generate the Keccak output BasicTable.
///
/// Table size = EFFECTIVE_BASE^TABLE_BITS = 2^8 = 256.
/// Column 1 is the base-11 sparse representation; column 2 is the binary
/// output; column 3 is always 0.
pub fn generate_keccak_output_table(id: BasicTableId, table_index: usize) -> BasicTable {
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    let table_size = pow64(EFFECTIVE_BASE, TABLE_BITS as u64);

    let mut counts = [0usize; TABLE_BITS];
    let mut column_values = [0u64; 2];

    for _ in 0..table_size {
        table.column_1.push(Fr::from(column_values[0]));
        table.column_2.push(Fr::from(column_values[1]));
        table.column_3.push(Fr::zero());
        column_values = get_column_values_for_next_row(&mut counts);
    }

    table.get_values_from_key = get_keccak_output_normalization_values;

    table.column_1_step_size = Fr::from(pow64(BASE, TABLE_BITS as u64));
    table.column_2_step_size = Fr::from(1u64 << TABLE_BITS);
    table.column_3_step_size = Fr::zero();

    table
}

/// Create the KeccakOutput MultiTable for converting sparse integers to binary.
///
/// 64-bit integers require ceil(64 / TABLE_BITS) lookups.
/// Column 1 uses base-11 step sizes; column 2 uses base-2 step sizes.
pub fn get_keccak_output_table(id: MultiTableId) -> MultiTable {
    let num_tables: usize = 64 / TABLE_BITS + if 64 % TABLE_BITS == 0 { 0 } else { 1 };

    let column_1_multiplier = pow64(BASE, TABLE_BITS as u64);
    let column_2_multiplier = 1u64 << TABLE_BITS;

    let mut table = MultiTable::new_repeated(
        Fr::from(column_1_multiplier),
        Fr::from(column_2_multiplier),
        Fr::zero(),
        num_tables,
    );

    table.id = id;
    for _ in 0..num_tables {
        table.slice_sizes.push(pow64(BASE, TABLE_BITS as u64));
        table.basic_table_ids.push(BasicTableId::KECCAK_OUTPUT);
        table
            .get_table_values
            .push(get_keccak_output_normalization_values);
    }

    table
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_output_normalization() {
        // 0 -> 0
        let r = get_keccak_output_normalization_values([0, 0]);
        assert_eq!(r[0], Fr::zero());
        assert_eq!(r[1], Fr::zero());

        // 1 (base-11 position 0 = 1) -> binary 1
        let r = get_keccak_output_normalization_values([1, 0]);
        assert_eq!(r[0], Fr::from(1u64));

        // 11 (base-11 position 1 = 1) -> binary 2
        let r = get_keccak_output_normalization_values([11, 0]);
        assert_eq!(r[0], Fr::from(2u64));

        // 12 (base-11: positions 0 and 1 both 1) -> binary 3
        let r = get_keccak_output_normalization_values([12, 0]);
        assert_eq!(r[0], Fr::from(3u64));
    }

    #[test]
    fn test_generate_keccak_output_table() {
        let table = generate_keccak_output_table(BasicTableId::KECCAK_OUTPUT, 0);
        // 2^8 = 256
        assert_eq!(table.column_1.len(), 256);
        assert_eq!(table.column_1[0], Fr::zero());
        assert_eq!(table.column_2[0], Fr::zero());
        assert_eq!(table.column_3[0], Fr::zero());
    }

    #[test]
    fn test_get_keccak_output_table() {
        let multi = get_keccak_output_table(MultiTableId::KECCAK_FORMAT_OUTPUT);
        assert_eq!(multi.id, MultiTableId::KECCAK_FORMAT_OUTPUT);
        // 64 / 8 = 8
        assert_eq!(multi.slice_sizes.len(), 8);
    }
}
