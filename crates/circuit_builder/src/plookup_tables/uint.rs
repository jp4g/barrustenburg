//! Unsigned integer XOR and AND lookup tables.
//!
//! C++ source: plookup_tables/uint.hpp
//!
//! Provides lookup tables for 8/16/32/64-bit XOR and AND operations,
//! decomposed into 6-bit slices (with 2-bit or 4-bit remainder).

use bbrs_ecc::curves::bn254::Fr;

use super::types::{BasicTable, BasicTableId, MultiTable, MultiTableId};

// rotate64 is trivially: val.rotate_left(num_bits as u32) for 64-bit rotation.
// But for sub-64-bit slices, rotation doesn't apply (num_rotated_output_bits=0 in all uses).
// The C++ rotate64 with 0 bits just returns the value unchanged.

/// XOR value function for 6-bit slices, no rotation.
pub fn get_xor_rotate_values_6_0(key: [u64; 2]) -> [Fr; 2] {
    [Fr::from(key[0] ^ key[1]), Fr::zero()]
}

/// XOR value function for 2-bit slices, no rotation.
pub fn get_xor_rotate_values_2_0(key: [u64; 2]) -> [Fr; 2] {
    [Fr::from(key[0] ^ key[1]), Fr::zero()]
}

/// XOR value function for 4-bit slices, no rotation.
pub fn get_xor_rotate_values_4_0(key: [u64; 2]) -> [Fr; 2] {
    [Fr::from(key[0] ^ key[1]), Fr::zero()]
}

/// AND value function for 6-bit slices, no rotation.
pub fn get_and_rotate_values_6_0(key: [u64; 2]) -> [Fr; 2] {
    [Fr::from(key[0] & key[1]), Fr::zero()]
}

/// AND value function for 2-bit slices, no rotation.
pub fn get_and_rotate_values_2_0(key: [u64; 2]) -> [Fr; 2] {
    [Fr::from(key[0] & key[1]), Fr::zero()]
}

/// AND value function for 4-bit slices, no rotation.
pub fn get_and_rotate_values_4_0(key: [u64; 2]) -> [Fr; 2] {
    [Fr::from(key[0] & key[1]), Fr::zero()]
}

/// Generate XOR rotate table with given slice size.
pub fn generate_xor_rotate_table(
    bits_per_slice: u64,
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
            table.column_3.push(Fr::from(i ^ j));
        }
    }

    table.get_values_from_key = match bits_per_slice {
        2 => get_xor_rotate_values_2_0,
        4 => get_xor_rotate_values_4_0,
        6 => get_xor_rotate_values_6_0,
        _ => get_xor_rotate_values_6_0,
    };

    table.column_1_step_size = Fr::from(base);
    table.column_2_step_size = Fr::from(base);
    table.column_3_step_size = Fr::from(base);

    table
}

/// Generate AND rotate table with given slice size.
pub fn generate_and_rotate_table(
    bits_per_slice: u64,
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
            table.column_3.push(Fr::from(i & j));
        }
    }

    table.get_values_from_key = match bits_per_slice {
        2 => get_and_rotate_values_2_0,
        4 => get_and_rotate_values_4_0,
        6 => get_and_rotate_values_6_0,
        _ => get_and_rotate_values_6_0,
    };

    table.column_1_step_size = Fr::from(base);
    table.column_2_step_size = Fr::from(base);
    table.column_3_step_size = Fr::from(base);

    table
}

/// Create XOR MultiTable for the given uint size (8, 16, 32, or 64).
pub fn get_uint_xor_table(uint_size: usize, id: MultiTableId) -> MultiTable {
    assert!(
        uint_size == 8 || uint_size == 16 || uint_size == 32 || uint_size == 64,
        "unsupported uint size for XOR table"
    );

    let table_bit_size = 6usize;
    let num_entries = uint_size / table_bit_size;
    let base = 1u64 << table_bit_size;

    let mut table = MultiTable::new_repeated(Fr::from(base), Fr::from(base), Fr::from(base), num_entries);
    table.id = id;

    for _ in 0..num_entries {
        table.slice_sizes.push(base);
        table
            .basic_table_ids
            .push(BasicTableId::UINT_XOR_SLICE_6_ROTATE_0);
        table.get_table_values.push(get_xor_rotate_values_6_0);
    }

    // Remaining bits
    let last_table_bit_size = uint_size - table_bit_size * num_entries;
    let last_slice_size = 1u64 << last_table_bit_size;
    table.slice_sizes.push(last_slice_size);

    if uint_size == 8 || uint_size == 32 {
        // 2 bits remaining (8 % 6 = 2, 32 % 6 = 2)
        table
            .basic_table_ids
            .push(BasicTableId::UINT_XOR_SLICE_2_ROTATE_0);
        table.get_table_values.push(get_xor_rotate_values_2_0);
    } else {
        // 4 bits remaining (16 % 6 = 4, 64 % 6 = 4)
        table
            .basic_table_ids
            .push(BasicTableId::UINT_XOR_SLICE_4_ROTATE_0);
        table.get_table_values.push(get_xor_rotate_values_4_0);
    }

    table
}

/// Create AND MultiTable for the given uint size (8, 16, 32, or 64).
pub fn get_uint_and_table(uint_size: usize, id: MultiTableId) -> MultiTable {
    assert!(
        uint_size == 8 || uint_size == 16 || uint_size == 32 || uint_size == 64,
        "unsupported uint size for AND table"
    );

    let table_bit_size = 6usize;
    let num_entries = uint_size / table_bit_size;
    let base = 1u64 << table_bit_size;

    let mut table = MultiTable::new_repeated(Fr::from(base), Fr::from(base), Fr::from(base), num_entries);
    table.id = id;

    for _ in 0..num_entries {
        table.slice_sizes.push(base);
        table
            .basic_table_ids
            .push(BasicTableId::UINT_AND_SLICE_6_ROTATE_0);
        table.get_table_values.push(get_and_rotate_values_6_0);
    }

    // Remaining bits
    let last_table_bit_size = uint_size - table_bit_size * num_entries;
    let last_slice_size = 1u64 << last_table_bit_size;
    table.slice_sizes.push(last_slice_size);

    if uint_size == 8 || uint_size == 32 {
        table
            .basic_table_ids
            .push(BasicTableId::UINT_AND_SLICE_2_ROTATE_0);
        table.get_table_values.push(get_and_rotate_values_2_0);
    } else {
        table
            .basic_table_ids
            .push(BasicTableId::UINT_AND_SLICE_4_ROTATE_0);
        table.get_table_values.push(get_and_rotate_values_4_0);
    }

    table
}
