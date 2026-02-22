//! ECC generator lookup tables for non-native field operations.
//!
//! C++ source: plookup_tables/non_native_group_generator.{hpp,cpp}
//!
//! Precomputes 256-entry tables of secp256k1 generator point multiples,
//! storing x/y coordinates as 68-bit limbs for non-native field simulation.

use std::sync::OnceLock;

use bbrs_ecc::curves::bn254::Fr;
use bbrs_ecc::curves::secp256k1::{G1Element, Secp256k1FqParams};
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_ecc::groups::element::Element;

use super::types::{BasicTable, BasicTableId, MultiTable, MultiTableId};

/// Number of limb bits in field simulation (68 bits per limb, 4 limbs for 272-bit capacity).
const NUM_LIMB_BITS: u32 = 68;

/// Precomputed generator tables for secp256k1.
struct GeneratorTables {
    xlo: [(Fr, Fr); 256],
    xhi: [(Fr, Fr); 256],
    ylo: [(Fr, Fr); 256],
    yhi: [(Fr, Fr); 256],
    xyprime: [(Fr, Fr); 256],
    endo_xlo: [(Fr, Fr); 256],
    endo_xhi: [(Fr, Fr); 256],
    endo_xyprime: [(Fr, Fr); 256],
}

static GENERATOR_TABLES: OnceLock<GeneratorTables> = OnceLock::new();

fn get_tables() -> &'static GeneratorTables {
    GENERATOR_TABLES.get_or_init(init_generator_tables)
}

fn init_generator_tables() -> GeneratorTables {
    let base_point = G1Element::one();
    let d2 = base_point.dbl();

    let mut point_table = vec![G1Element::one(); 256];
    point_table[128] = base_point;
    for i in 1..128 {
        point_table[i + 128] = point_table[i + 127] + d2;
    }
    for i in 0..128 {
        point_table[127 - i] = -point_table[128 + i];
    }
    Element::batch_normalize(&mut point_table);

    // cube root of unity for endomorphism: beta
    let beta = Field::<Secp256k1FqParams>::from_raw(Secp256k1FqParams::CUBE_ROOT);

    let mut tables = GeneratorTables {
        xlo: [(Fr::zero(), Fr::zero()); 256],
        xhi: [(Fr::zero(), Fr::zero()); 256],
        ylo: [(Fr::zero(), Fr::zero()); 256],
        yhi: [(Fr::zero(), Fr::zero()); 256],
        xyprime: [(Fr::zero(), Fr::zero()); 256],
        endo_xlo: [(Fr::zero(), Fr::zero()); 256],
        endo_xhi: [(Fr::zero(), Fr::zero()); 256],
        endo_xyprime: [(Fr::zero(), Fr::zero()); 256],
    };

    for i in 0..256 {
        let affine = point_table[i].to_affine();
        let x_field = affine.x;
        let y_field = affine.y;
        let endo_x_field = x_field * beta;

        // Get standard-form (non-Montgomery) limbs
        let x_std = x_field.from_montgomery_form();
        let y_std = y_field.from_montgomery_form();
        let endo_x_std = endo_x_field.from_montgomery_form();

        // Extract 68-bit limbs from x
        let x0 = extract_limb_68(&x_std.data, 0);
        let x1 = extract_limb_68(&x_std.data, 1);
        let x2 = extract_limb_68(&x_std.data, 2);
        let x3 = extract_limb_68(&x_std.data, 3);

        // Extract 68-bit limbs from endo_x
        let ex0 = extract_limb_68(&endo_x_std.data, 0);
        let ex1 = extract_limb_68(&endo_x_std.data, 1);
        let ex2 = extract_limb_68(&endo_x_std.data, 2);
        let ex3 = extract_limb_68(&endo_x_std.data, 3);

        // Extract 68-bit limbs from y
        let y0 = extract_limb_68(&y_std.data, 0);
        let y1 = extract_limb_68(&y_std.data, 1);
        let y2 = extract_limb_68(&y_std.data, 2);
        let y3 = extract_limb_68(&y_std.data, 3);

        tables.xlo[i] = (Fr::from(x0), Fr::from(x1));
        tables.xhi[i] = (Fr::from(x2), Fr::from(x3));
        tables.ylo[i] = (Fr::from(y0), Fr::from(y1));
        tables.yhi[i] = (Fr::from(y2), Fr::from(y3));

        // Prime-basis: store as Fr values derived from the full coordinate
        // fr(x) and fr(y) — these are the coordinates reduced mod BN254 Fr
        tables.xyprime[i] = (
            Fr::from_limbs(x_std.data),
            Fr::from_limbs(y_std.data),
        );
        tables.endo_xyprime[i] = (
            Fr::from_limbs(endo_x_std.data),
            Fr::from_limbs(y_std.data),
        );

        tables.endo_xlo[i] = (Fr::from(ex0), Fr::from(ex1));
        tables.endo_xhi[i] = (Fr::from(ex2), Fr::from(ex3));
    }

    tables
}

/// Extract a 68-bit limb from a [u64; 4] array at the given limb index (0-3).
fn extract_limb_68(data: &[u64; 4], limb_idx: u32) -> u64 {
    let bit_offset = limb_idx * NUM_LIMB_BITS;
    let mask = (1u128 << NUM_LIMB_BITS) - 1;

    // Reconstruct relevant portion as u128
    let word_idx = (bit_offset / 64) as usize;
    let bit_in_word = bit_offset % 64;

    if word_idx >= 4 {
        return 0;
    }

    let mut val = data[word_idx] as u128;
    if word_idx + 1 < 4 {
        val |= (data[word_idx + 1] as u128) << 64;
    }
    ((val >> bit_in_word) & mask) as u64
}

// ---------------------------------------------------------------------------
// Value lookup functions (fn pointers for BasicTable)
// ---------------------------------------------------------------------------

pub fn get_xlo_values(key: [u64; 2]) -> [Fr; 2] {
    let t = get_tables();
    let idx = key[0] as usize;
    [t.xlo[idx].0, t.xlo[idx].1]
}

pub fn get_xhi_values(key: [u64; 2]) -> [Fr; 2] {
    let t = get_tables();
    let idx = key[0] as usize;
    [t.xhi[idx].0, t.xhi[idx].1]
}

pub fn get_ylo_values(key: [u64; 2]) -> [Fr; 2] {
    let t = get_tables();
    let idx = key[0] as usize;
    [t.ylo[idx].0, t.ylo[idx].1]
}

pub fn get_yhi_values(key: [u64; 2]) -> [Fr; 2] {
    let t = get_tables();
    let idx = key[0] as usize;
    [t.yhi[idx].0, t.yhi[idx].1]
}

pub fn get_xyprime_values(key: [u64; 2]) -> [Fr; 2] {
    let t = get_tables();
    let idx = key[0] as usize;
    [t.xyprime[idx].0, t.xyprime[idx].1]
}

pub fn get_xlo_endo_values(key: [u64; 2]) -> [Fr; 2] {
    let t = get_tables();
    let idx = key[0] as usize;
    [t.endo_xlo[idx].0, t.endo_xlo[idx].1]
}

pub fn get_xhi_endo_values(key: [u64; 2]) -> [Fr; 2] {
    let t = get_tables();
    let idx = key[0] as usize;
    [t.endo_xhi[idx].0, t.endo_xhi[idx].1]
}

pub fn get_xyprime_endo_values(key: [u64; 2]) -> [Fr; 2] {
    let t = get_tables();
    let idx = key[0] as usize;
    [t.endo_xyprime[idx].0, t.endo_xyprime[idx].1]
}

// ---------------------------------------------------------------------------
// BasicTable generators
// ---------------------------------------------------------------------------

fn generate_table_from_precomputed(
    id: BasicTableId,
    table_index: usize,
    data: &[(Fr, Fr); 256],
    get_values: fn([u64; 2]) -> [Fr; 2],
) -> BasicTable {
    let _ = get_tables(); // ensure initialized
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    for i in 0..256 {
        table.column_1.push(Fr::from(i as u64));
        table.column_2.push(data[i].0);
        table.column_3.push(data[i].1);
    }

    table.get_values_from_key = get_values;
    table.column_1_step_size = Fr::zero();
    table.column_2_step_size = Fr::zero();
    table.column_3_step_size = Fr::zero();

    table
}

pub fn generate_xlo_table(id: BasicTableId, table_index: usize) -> BasicTable {
    generate_table_from_precomputed(id, table_index, &get_tables().xlo, get_xlo_values)
}

pub fn generate_xhi_table(id: BasicTableId, table_index: usize) -> BasicTable {
    generate_table_from_precomputed(id, table_index, &get_tables().xhi, get_xhi_values)
}

pub fn generate_ylo_table(id: BasicTableId, table_index: usize) -> BasicTable {
    generate_table_from_precomputed(id, table_index, &get_tables().ylo, get_ylo_values)
}

pub fn generate_yhi_table(id: BasicTableId, table_index: usize) -> BasicTable {
    generate_table_from_precomputed(id, table_index, &get_tables().yhi, get_yhi_values)
}

pub fn generate_xyprime_table(id: BasicTableId, table_index: usize) -> BasicTable {
    generate_table_from_precomputed(id, table_index, &get_tables().xyprime, get_xyprime_values)
}

pub fn generate_xlo_endo_table(id: BasicTableId, table_index: usize) -> BasicTable {
    generate_table_from_precomputed(id, table_index, &get_tables().endo_xlo, get_xlo_endo_values)
}

pub fn generate_xhi_endo_table(id: BasicTableId, table_index: usize) -> BasicTable {
    generate_table_from_precomputed(id, table_index, &get_tables().endo_xhi, get_xhi_endo_values)
}

pub fn generate_xyprime_endo_table(id: BasicTableId, table_index: usize) -> BasicTable {
    generate_table_from_precomputed(
        id,
        table_index,
        &get_tables().endo_xyprime,
        get_xyprime_endo_values,
    )
}

// ---------------------------------------------------------------------------
// MultiTable generators
// ---------------------------------------------------------------------------

fn get_ecc_multi_table(
    id: MultiTableId,
    basic_id: BasicTableId,
    get_values: fn([u64; 2]) -> [Fr; 2],
) -> MultiTable {
    let mut table = MultiTable::new_repeated(Fr::from(256u64), Fr::zero(), Fr::zero(), 1);
    table.id = id;
    table.slice_sizes.push(512);
    table.basic_table_ids.push(basic_id);
    table.get_table_values.push(get_values);
    table
}

pub fn get_xlo_multi_table(id: MultiTableId, basic_id: BasicTableId) -> MultiTable {
    get_ecc_multi_table(id, basic_id, get_xlo_values)
}

pub fn get_xhi_multi_table(id: MultiTableId, basic_id: BasicTableId) -> MultiTable {
    get_ecc_multi_table(id, basic_id, get_xhi_values)
}

pub fn get_ylo_multi_table(id: MultiTableId, basic_id: BasicTableId) -> MultiTable {
    get_ecc_multi_table(id, basic_id, get_ylo_values)
}

pub fn get_yhi_multi_table(id: MultiTableId, basic_id: BasicTableId) -> MultiTable {
    get_ecc_multi_table(id, basic_id, get_yhi_values)
}

pub fn get_xyprime_multi_table(id: MultiTableId, basic_id: BasicTableId) -> MultiTable {
    get_ecc_multi_table(id, basic_id, get_xyprime_values)
}

pub fn get_xlo_endo_multi_table(id: MultiTableId, basic_id: BasicTableId) -> MultiTable {
    get_ecc_multi_table(id, basic_id, get_xlo_endo_values)
}

pub fn get_xhi_endo_multi_table(id: MultiTableId, basic_id: BasicTableId) -> MultiTable {
    get_ecc_multi_table(id, basic_id, get_xhi_endo_values)
}

pub fn get_xyprime_endo_multi_table(id: MultiTableId, basic_id: BasicTableId) -> MultiTable {
    get_ecc_multi_table(id, basic_id, get_xyprime_endo_values)
}
