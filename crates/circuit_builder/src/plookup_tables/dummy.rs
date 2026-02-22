//! Dummy plookup tables for UltraHonk.
//!
//! C++ source: plookup_tables/dummy.hpp
//!
//! These tables ensure that the table, sorted, and lookup selector polynomials
//! are non-zero in UltraHonk proofs.

use bbrs_ecc::curves::bn254::Fr;

use super::types::{BasicTable, BasicTableId, MultiTable, MultiTableId};

/// Get value from key for dummy table 1 (HONK_DUMMY_BASIC1).
pub fn get_value_from_key_basic1(key: [u64; 2]) -> [Fr; 2] {
    let id_val = BasicTableId::HONK_DUMMY_BASIC1.0 as u64;
    [Fr::from(key[0] * 3 + key[1] * 4 + id_val * 0x1337), Fr::zero()]
}

/// Get value from key for dummy table 2 (HONK_DUMMY_BASIC2).
pub fn get_value_from_key_basic2(key: [u64; 2]) -> [Fr; 2] {
    let id_val = BasicTableId::HONK_DUMMY_BASIC2.0 as u64;
    [Fr::from(key[0] * 3 + key[1] * 4 + id_val * 0x1337), Fr::zero()]
}

/// Generate a dummy basic table.
fn generate_honk_dummy_table(
    id: BasicTableId,
    table_index: usize,
    get_values: fn([u64; 2]) -> [Fr; 2],
) -> BasicTable {
    let base: u64 = 1 << 1; // 2
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = true;

    for i in 0..base {
        for j in 0..base {
            table.column_1.push(Fr::from(i));
            table.column_2.push(Fr::from(j));
            table
                .column_3
                .push(Fr::from(i * 3 + j * 4 + id.0 as u64 * 0x1337));
        }
    }

    table.get_values_from_key = get_values;
    table.column_1_step_size = Fr::from(base);
    table.column_2_step_size = Fr::from(base);
    table.column_3_step_size = Fr::from(base);

    table
}

/// Generate HONK_DUMMY_BASIC1 table.
pub fn generate_honk_dummy_table_basic1(id: BasicTableId, table_index: usize) -> BasicTable {
    generate_honk_dummy_table(id, table_index, get_value_from_key_basic1)
}

/// Generate HONK_DUMMY_BASIC2 table.
pub fn generate_honk_dummy_table_basic2(id: BasicTableId, table_index: usize) -> BasicTable {
    generate_honk_dummy_table(id, table_index, get_value_from_key_basic2)
}

/// Create the Honk dummy MultiTable (2 lookups on 2 basic tables).
pub fn get_honk_dummy_multitable() -> MultiTable {
    let num_elements = 1u64 << 1; // 2
    let num_lookups = 2;

    let mut table = MultiTable::new_repeated(
        Fr::from(num_elements),
        Fr::from(num_elements),
        Fr::from(num_elements),
        num_lookups,
    );
    table.id = MultiTableId::HONK_DUMMY_MULTI;

    table.slice_sizes.push(num_elements);
    table.basic_table_ids.push(BasicTableId::HONK_DUMMY_BASIC1);
    table
        .get_table_values
        .push(get_value_from_key_basic1);

    table.slice_sizes.push(num_elements);
    table.basic_table_ids.push(BasicTableId::HONK_DUMMY_BASIC2);
    table
        .get_table_values
        .push(get_value_from_key_basic2);

    table
}
