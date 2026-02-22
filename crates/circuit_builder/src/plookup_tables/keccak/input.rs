//! Keccak input plookup tables.
//!
//! C++ source: plookup_tables/keccak/keccak_input.hpp
//!
//! Generates plookup tables that convert 64-bit integers into a base-11 sparse
//! representation used by the Keccak hash algorithm circuit.
//!
//! Each 64-bit hash lane is represented as P = sum_{j=0}^{63} b_j * 11^j.
//! This table maps binary integer slices into base-11 integer slices.
//! It also extracts the most significant bit of the input slice (used by
//! stdlib::keccak to efficiently left-rotate by 1 bit).

use bbrs_ecc::curves::bn254::Fr;
use bbrs_numeric::bitop::pow64;

use crate::plookup_tables::types::{BasicTable, BasicTableId, MultiTable, MultiTableId};

/// Base used for sparse representation.
const BASE: u64 = 11;

/// Number of bits per lookup table slice.
const TABLE_BITS: u64 = 8;

/// Map a binary integer into base-11 sparse form.
///
/// For each bit i of `input`, if bit i is set, add BASE^i to the result.
/// This mirrors C++ `numeric::map_into_sparse_form<BASE>`.
fn map_into_sparse_form(input: u64) -> u64 {
    let mut out: u64 = 0;
    let mut base_power: u64 = 1;
    for i in 0..32 {
        let sparse_bit = (input >> i) & 1;
        if sparse_bit != 0 {
            out += base_power;
        }
        if i < 31 {
            base_power = base_power.wrapping_mul(BASE);
        }
    }
    out
}

/// Given a table input value, return the table output value.
///
/// key[0] = table input (binary integer slice).
/// Returns [sparse_form, msb_extraction].
///
/// The MSB extraction shifts by `(64 % TABLE_BITS == 0) ? TABLE_BITS - 1 : (64 % TABLE_BITS) - 1`.
pub fn get_keccak_input_values(key: [u64; 2]) -> [Fr; 2] {
    let t0 = map_into_sparse_form(key[0]);
    let msb_shift: u64 = if 64 % TABLE_BITS == 0 {
        TABLE_BITS - 1
    } else {
        (64 % TABLE_BITS) - 1
    };
    let t1 = key[0] >> msb_shift;
    [Fr::from(t0), Fr::from(t1)]
}

/// Generate the basic plookup table that maps a TABLE_BITS-slice of a base-2
/// integer into a base-11 representation.
pub fn generate_keccak_input_table(id: BasicTableId, table_index: usize) -> BasicTable {
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    let table_size = 1u64 << TABLE_BITS;
    let msb_shift: u64 = if 64 % TABLE_BITS == 0 {
        TABLE_BITS - 1
    } else {
        (64 % TABLE_BITS) - 1
    };

    for i in 0..table_size {
        let source = i;
        let target = map_into_sparse_form(source);
        table.column_1.push(Fr::from(source));
        table.column_2.push(Fr::from(target));
        table.column_3.push(Fr::from(source >> msb_shift));
    }

    table.get_values_from_key = get_keccak_input_values;

    let sparse_step_size = pow64(BASE, TABLE_BITS);
    table.column_1_step_size = Fr::from(1u64 << TABLE_BITS);
    table.column_2_step_size = Fr::from(sparse_step_size);
    table.column_3_step_size = Fr::from(sparse_step_size);

    table
}

/// Create the KeccakInput MultiTable used by plookup to generate a sequence
/// of lookups.
///
/// 64-bit integers require 8 lookups of 8 bits each. The MultiTable defines
/// coefficient values that allow plookup to derive column entries from
/// relative differences between wire values.
pub fn get_keccak_input_table(id: MultiTableId) -> MultiTable {
    let num_entries = 8;

    let col1_coeff = Fr::from(1u64 << 8);
    let col2_coeff = Fr::from(pow64(11, 8));
    let col3_coeff = Fr::zero();

    let mut table = MultiTable::new_repeated(col1_coeff, col2_coeff, col3_coeff, num_entries);
    table.id = id;

    for _ in 0..num_entries {
        table.slice_sizes.push(1u64 << 8);
        table.basic_table_ids.push(BasicTableId::KECCAK_INPUT);
        table.get_table_values.push(get_keccak_input_values);
    }

    table
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_map_into_sparse_form() {
        // 0 -> 0
        assert_eq!(map_into_sparse_form(0), 0);
        // 1 -> 11^0 = 1
        assert_eq!(map_into_sparse_form(1), 1);
        // 2 (bit 1 set) -> 11^1 = 11
        assert_eq!(map_into_sparse_form(2), 11);
        // 3 (bits 0,1 set) -> 11^0 + 11^1 = 1 + 11 = 12
        assert_eq!(map_into_sparse_form(3), 12);
        // 4 (bit 2 set) -> 11^2 = 121
        assert_eq!(map_into_sparse_form(4), 121);
    }

    #[test]
    fn test_get_keccak_input_values() {
        let result = get_keccak_input_values([0, 0]);
        assert_eq!(result[0], Fr::zero());

        let result = get_keccak_input_values([1, 0]);
        assert_eq!(result[0], Fr::from(1u64));
    }

    #[test]
    fn test_generate_keccak_input_table() {
        let table = generate_keccak_input_table(BasicTableId::KECCAK_INPUT, 0);
        assert_eq!(table.column_1.len(), 256); // 2^8
        assert_eq!(table.column_2.len(), 256);
        assert_eq!(table.column_3.len(), 256);
        // First row should be all zeros
        assert_eq!(table.column_1[0], Fr::zero());
        assert_eq!(table.column_2[0], Fr::zero());
        assert_eq!(table.column_3[0], Fr::zero());
    }

    #[test]
    fn test_get_keccak_input_table() {
        let multi = get_keccak_input_table(MultiTableId::KECCAK_FORMAT_INPUT);
        assert_eq!(multi.id, MultiTableId::KECCAK_FORMAT_INPUT);
        assert_eq!(multi.slice_sizes.len(), 8);
        assert_eq!(multi.basic_table_ids.len(), 8);
        assert_eq!(multi.get_table_values.len(), 8);
    }
}
