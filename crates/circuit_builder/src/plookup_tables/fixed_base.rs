//! Fixed-base scalar multiplication plookup tables.
//!
//! C++ source: plookup_tables/fixed_base/fixed_base.{hpp,cpp}
//!
//! Precomputes lookup tables for fixed-base scalar multiplication over Grumpkin.
//! Each table stores affine points of the form [G_offset] + j * [P] for j in [0, 512).
//! Offset generators prevent point-at-infinity edge cases.

use std::sync::OnceLock;

use bbrs_crypto::generators::{derive_generators, precomputed::DEFAULT_GENERATORS};
use bbrs_ecc::curves::bn254::Fr;
use bbrs_ecc::curves::grumpkin::{self, GrumpkinFr};
use bbrs_ecc::groups::element::Element;

use super::fixed_base_params::FixedBaseParams;
use super::types::{BasicTable, BasicTableId, GetValuesFromKey, MultiTable, MultiTableId};

type GrumpkinAffine = grumpkin::G1Affine;
type GrumpkinElement = grumpkin::G1Element;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const MAX_TABLE_SIZE: usize = FixedBaseParams::MAX_TABLE_SIZE as usize;
const BITS_PER_TABLE: u64 = FixedBaseParams::BITS_PER_TABLE;
const NUM_FIXED_BASE_MULTI_TABLES: usize = FixedBaseParams::NUM_FIXED_BASE_MULTI_TABLES as usize;
#[cfg(test)]
const NUM_TABLES_PER_LO: usize = FixedBaseParams::NUM_TABLES_PER_LO_MULTITABLE as usize;
#[cfg(test)]
const NUM_TABLES_PER_HI: usize = FixedBaseParams::NUM_TABLES_PER_HI_MULTITABLE as usize;
const MAX_NUM_TABLES: usize = FixedBaseParams::MAX_NUM_TABLES_IN_MULTITABLE as usize;

/// MAX_LO_SCALAR = 2^128, used to compute hi base points.
fn max_lo_scalar() -> GrumpkinFr {
    GrumpkinFr::from_limbs([0, 0, 1, 0])
}

// ---------------------------------------------------------------------------
// Generator points
// ---------------------------------------------------------------------------

/// Left generator point (precomputed_generators[0]).
pub fn lhs_generator_point() -> GrumpkinAffine {
    DEFAULT_GENERATORS[0]
}

/// Right generator point (precomputed_generators[1]).
pub fn rhs_generator_point() -> GrumpkinAffine {
    DEFAULT_GENERATORS[1]
}

// ---------------------------------------------------------------------------
// Precomputed tables (lazy init)
// ---------------------------------------------------------------------------

/// All fixed-base scalar multiplication tables.
/// `tables[multitable_index][table_index]` = Vec of affine points.
///
/// Layout:
///   [0] = LHS LO (15 tables)
///   [1] = LHS HI (14 tables)
///   [2] = RHS LO (15 tables)
///   [3] = RHS HI (14 tables)
struct FixedBaseTables {
    tables: [Vec<Vec<GrumpkinAffine>>; NUM_FIXED_BASE_MULTI_TABLES],
    offset_generators: [GrumpkinAffine; NUM_FIXED_BASE_MULTI_TABLES],
}

static FIXED_BASE_TABLES: OnceLock<FixedBaseTables> = OnceLock::new();

fn get_all_tables() -> &'static FixedBaseTables {
    FIXED_BASE_TABLES.get_or_init(init_fixed_base_tables)
}

/// Serialize a Grumpkin affine point to bytes (32 bytes x big-endian + 32 bytes y big-endian).
/// Matches C++ `write(buf, affine_element)`.
fn serialize_affine(point: &GrumpkinAffine) -> Vec<u8> {
    let mut buf = Vec::with_capacity(64);
    let x_std = point.x.from_montgomery_form();
    let y_std = point.y.from_montgomery_form();
    for field_data in &[x_std.data, y_std.data] {
        for i in (0..4).rev() {
            for j in (0..8).rev() {
                buf.push((field_data[i] >> (j * 8)) as u8);
            }
        }
    }
    buf
}

/// Generate a single lookup table: { [G] + j*[P] : j = 0..MAX_TABLE_SIZE-1 }.
fn generate_single_lookup_table(
    base_point: &GrumpkinAffine,
    offset_generator: &GrumpkinAffine,
) -> Vec<GrumpkinAffine> {
    let mut table_raw = Vec::with_capacity(MAX_TABLE_SIZE);
    let mut accumulator = GrumpkinElement::from_affine(offset_generator);
    let base_proj = GrumpkinElement::from_affine(base_point);

    for _ in 0..MAX_TABLE_SIZE {
        table_raw.push(accumulator);
        accumulator = accumulator + base_proj;
    }

    Element::batch_normalize(&mut table_raw);

    table_raw.iter().map(|e| e.to_affine()).collect()
}

/// Generate all basic tables for a scalar mul of `num_bits` bits.
fn generate_tables(input: &GrumpkinAffine, num_bits: u64) -> Vec<Vec<GrumpkinAffine>> {
    let num_tables = ((num_bits + BITS_PER_TABLE - 1) / BITS_PER_TABLE) as usize;

    let input_buf = serialize_affine(input);
    let offset_generators = derive_generators(&input_buf, num_tables, 0);

    let mut accumulator = GrumpkinElement::from_affine(input);
    let mut tables = Vec::with_capacity(num_tables);

    for i in 0..num_tables {
        tables.push(generate_single_lookup_table(
            &accumulator.to_affine(),
            &offset_generators[i],
        ));
        for _ in 0..BITS_PER_TABLE {
            accumulator = accumulator.dbl();
        }
    }

    tables
}

/// Compute the sum of offset generators for a scalar mul of `num_bits` bits.
fn compute_generator_offset(input: &GrumpkinAffine, num_bits: u64) -> GrumpkinAffine {
    let num_tables = ((num_bits + BITS_PER_TABLE - 1) / BITS_PER_TABLE) as usize;

    let input_buf = serialize_affine(input);
    let offset_generators = derive_generators(&input_buf, num_tables, 0);

    let mut total = GrumpkinElement::infinity();
    for g in &offset_generators {
        total = total + GrumpkinElement::from_affine(g);
    }

    total.to_affine()
}

fn init_fixed_base_tables() -> FixedBaseTables {
    let lhs_lo = lhs_generator_point();
    let rhs_lo = rhs_generator_point();
    let lhs_hi_proj = GrumpkinElement::from_affine(&lhs_lo).mul(&max_lo_scalar());
    let rhs_hi_proj = GrumpkinElement::from_affine(&rhs_lo).mul(&max_lo_scalar());
    let lhs_hi = lhs_hi_proj.to_affine();
    let rhs_hi = rhs_hi_proj.to_affine();

    let lo_bits = FixedBaseParams::BITS_PER_LO_SCALAR;
    let hi_bits = FixedBaseParams::BITS_PER_HI_SCALAR;

    FixedBaseTables {
        tables: [
            generate_tables(&lhs_lo, lo_bits),
            generate_tables(&lhs_hi, hi_bits),
            generate_tables(&rhs_lo, lo_bits),
            generate_tables(&rhs_hi, hi_bits),
        ],
        offset_generators: [
            compute_generator_offset(&lhs_lo, lo_bits),
            compute_generator_offset(&lhs_hi, hi_bits),
            compute_generator_offset(&rhs_lo, lo_bits),
            compute_generator_offset(&rhs_hi, hi_bits),
        ],
    }
}

// ---------------------------------------------------------------------------
// Value lookup functions (const-generic dispatch)
// ---------------------------------------------------------------------------

/// Look up (x, y) coordinates from precomputed table.
/// MI = multitable index, TI = table index within multitable.
fn get_fixed_base_values<const MI: usize, const TI: usize>(key: [u64; 2]) -> [Fr; 2] {
    let tables = get_all_tables();
    let idx = key[0] as usize;
    let p = tables.tables[MI][TI][idx];
    // GrumpkinFq = BN254 Fr, so point.x and point.y are already Fr
    [p.x, p.y]
}

/// Select the correct monomorphized function pointer for (mi, ti).
pub fn select_fixed_base_get_values(mi: usize, ti: usize) -> GetValuesFromKey {
    macro_rules! check {
        ($mi:literal, $ti:literal) => {
            if mi == $mi && ti == $ti {
                return get_fixed_base_values::<$mi, $ti>;
            }
        };
    }

    // Multitable 0 (LHS LO): 15 tables
    check!(0, 0); check!(0, 1); check!(0, 2); check!(0, 3); check!(0, 4);
    check!(0, 5); check!(0, 6); check!(0, 7); check!(0, 8); check!(0, 9);
    check!(0, 10); check!(0, 11); check!(0, 12); check!(0, 13); check!(0, 14);

    // Multitable 1 (LHS HI): 14 tables
    check!(1, 0); check!(1, 1); check!(1, 2); check!(1, 3); check!(1, 4);
    check!(1, 5); check!(1, 6); check!(1, 7); check!(1, 8); check!(1, 9);
    check!(1, 10); check!(1, 11); check!(1, 12); check!(1, 13);

    // Multitable 2 (RHS LO): 15 tables
    check!(2, 0); check!(2, 1); check!(2, 2); check!(2, 3); check!(2, 4);
    check!(2, 5); check!(2, 6); check!(2, 7); check!(2, 8); check!(2, 9);
    check!(2, 10); check!(2, 11); check!(2, 12); check!(2, 13); check!(2, 14);

    // Multitable 3 (RHS HI): 14 tables
    check!(3, 0); check!(3, 1); check!(3, 2); check!(3, 3); check!(3, 4);
    check!(3, 5); check!(3, 6); check!(3, 7); check!(3, 8); check!(3, 9);
    check!(3, 10); check!(3, 11); check!(3, 12); check!(3, 13);

    panic!(
        "No fixed base get_values function for multitable_index={}, table_index={}",
        mi, ti
    );
}

// ---------------------------------------------------------------------------
// BasicTable generator
// ---------------------------------------------------------------------------

/// Generate a BasicTable for a specific bit-slice of fixed-base scalar multiplication.
///
/// Mirrors C++ `table::generate_basic_fixed_base_table<multitable_index>`.
pub fn generate_basic_fixed_base_table(
    multitable_index: usize,
    id: BasicTableId,
    basic_table_index: usize,
    table_index: usize,
) -> BasicTable {
    assert!(multitable_index < NUM_FIXED_BASE_MULTI_TABLES);
    assert!(table_index < MAX_NUM_TABLES);

    let multitable_bits = if multitable_index % 2 == 0 {
        FixedBaseParams::BITS_PER_LO_SCALAR
    } else {
        FixedBaseParams::BITS_PER_HI_SCALAR
    };
    let bits_covered = BITS_PER_TABLE * table_index as u64;
    let is_small = (multitable_bits - bits_covered) < BITS_PER_TABLE;
    let table_bits = if is_small {
        multitable_bits - bits_covered
    } else {
        BITS_PER_TABLE
    };
    let table_size = 1usize << table_bits;

    let all = get_all_tables();
    let basic_table = &all.tables[multitable_index][table_index];

    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = basic_table_index;
    table.use_twin_keys = false;

    for i in 0..table_size {
        table.column_1.push(Fr::from(i as u64));
        // GrumpkinFq = BN254 Fr, direct field reuse
        table.column_2.push(basic_table[i].x);
        table.column_3.push(basic_table[i].y);
    }

    table.get_values_from_key = select_fixed_base_get_values(multitable_index, table_index);

    table.column_1_step_size = Fr::from(table_size as u64);
    table.column_2_step_size = Fr::zero(); // No accumulation for coordinates
    table.column_3_step_size = Fr::zero();

    table
}

// ---------------------------------------------------------------------------
// MultiTable generator
// ---------------------------------------------------------------------------

/// Create a fixed-base scalar multiplication MultiTable.
///
/// Mirrors C++ `table::get_fixed_base_table<multitable_index, num_bits>`.
pub fn get_fixed_base_table(
    multitable_index: usize,
    num_bits: u64,
    id: MultiTableId,
) -> MultiTable {
    assert!(
        num_bits == FixedBaseParams::BITS_PER_LO_SCALAR
            || num_bits == FixedBaseParams::BITS_PER_HI_SCALAR
    );
    let num_tables = ((num_bits + BITS_PER_TABLE - 1) / BITS_PER_TABLE) as usize;

    let basic_table_id_bases = [
        BasicTableId::FIXED_BASE_0_0,
        BasicTableId::FIXED_BASE_1_0,
        BasicTableId::FIXED_BASE_2_0,
        BasicTableId::FIXED_BASE_3_0,
    ];

    // Column 1 accumulates (scalar slices), columns 2&3 do NOT accumulate (individual points)
    let mut table = MultiTable::new_repeated(
        Fr::from(MAX_TABLE_SIZE as u64),
        Fr::zero(),
        Fr::zero(),
        num_tables,
    );
    table.id = id;
    table.get_table_values.resize(num_tables, super::types::default_get_values);
    table.basic_table_ids.resize(num_tables, BasicTableId::XOR);

    for i in 0..num_tables {
        table.slice_sizes.push(MAX_TABLE_SIZE as u64);
        table.get_table_values[i] = select_fixed_base_get_values(multitable_index, i);
        let base_id = basic_table_id_bases[multitable_index].0;
        table.basic_table_ids[i] = BasicTableId(base_id + i);
    }

    table
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Returns true if a precomputed table exists for the given point.
pub fn lookup_table_exists_for_point(input: &GrumpkinAffine) -> bool {
    *input == lhs_generator_point() || *input == rhs_generator_point()
}

/// Get the (LO, HI) MultiTable IDs for a given generator point.
pub fn get_lookup_table_ids_for_point(input: &GrumpkinAffine) -> [MultiTableId; 2] {
    assert!(
        lookup_table_exists_for_point(input),
        "No fixed-base table exists for input point"
    );
    if *input == lhs_generator_point() {
        [MultiTableId::FIXED_BASE_LEFT_LO, MultiTableId::FIXED_BASE_LEFT_HI]
    } else {
        [MultiTableId::FIXED_BASE_RIGHT_LO, MultiTableId::FIXED_BASE_RIGHT_HI]
    }
}

/// Get the offset generator for a given table ID.
pub fn get_generator_offset_for_table_id(table_id: MultiTableId) -> GrumpkinAffine {
    let all = get_all_tables();
    if table_id == MultiTableId::FIXED_BASE_LEFT_LO {
        all.offset_generators[0]
    } else if table_id == MultiTableId::FIXED_BASE_LEFT_HI {
        all.offset_generators[1]
    } else if table_id == MultiTableId::FIXED_BASE_RIGHT_LO {
        all.offset_generators[2]
    } else if table_id == MultiTableId::FIXED_BASE_RIGHT_HI {
        all.offset_generators[3]
    } else {
        panic!("Invalid table_id for fixed base offset generator")
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generator_points_on_curve() {
        let lhs = lhs_generator_point();
        let rhs = rhs_generator_point();
        assert!(lhs.on_curve(), "LHS generator not on curve");
        assert!(rhs.on_curve(), "RHS generator not on curve");
        assert_ne!(lhs, rhs, "LHS and RHS generators should differ");
    }

    #[test]
    fn test_lookup_table_exists() {
        let lhs = lhs_generator_point();
        let rhs = rhs_generator_point();
        assert!(lookup_table_exists_for_point(&lhs));
        assert!(lookup_table_exists_for_point(&rhs));
    }

    #[test]
    fn test_tables_initialized() {
        let tables = get_all_tables();
        // LHS LO: 15 tables, each with up to 512 entries
        assert_eq!(tables.tables[0].len(), NUM_TABLES_PER_LO);
        assert_eq!(tables.tables[0][0].len(), MAX_TABLE_SIZE);
        // LHS HI: 14 tables
        assert_eq!(tables.tables[1].len(), NUM_TABLES_PER_HI);
        // RHS LO: 15 tables
        assert_eq!(tables.tables[2].len(), NUM_TABLES_PER_LO);
        // RHS HI: 14 tables
        assert_eq!(tables.tables[3].len(), NUM_TABLES_PER_HI);
    }

    #[test]
    fn test_table_points_on_curve() {
        let tables = get_all_tables();
        // Spot check: first table of LHS LO, first 10 entries
        for i in 0..10 {
            assert!(
                tables.tables[0][0][i].on_curve(),
                "Point {} in table[0][0] not on curve",
                i
            );
        }
    }

    #[test]
    fn test_offset_generators_on_curve() {
        let tables = get_all_tables();
        for (i, g) in tables.offset_generators.iter().enumerate() {
            assert!(g.on_curve(), "Offset generator {} not on curve", i);
        }
    }

    #[test]
    fn test_get_values_returns_point_coords() {
        let tables = get_all_tables();
        let fn_ptr = select_fixed_base_get_values(0, 0);
        let result = fn_ptr([0, 0]);
        let expected_point = tables.tables[0][0][0];
        assert_eq!(result[0], expected_point.x);
        assert_eq!(result[1], expected_point.y);
    }

    #[test]
    fn test_generate_basic_table() {
        let table = generate_basic_fixed_base_table(
            0,
            BasicTableId::FIXED_BASE_0_0,
            0,
            0,
        );
        assert_eq!(table.column_1.len(), MAX_TABLE_SIZE);
        assert_eq!(table.column_2.len(), MAX_TABLE_SIZE);
        assert_eq!(table.column_3.len(), MAX_TABLE_SIZE);
        assert!(!table.use_twin_keys);
    }

    #[test]
    fn test_get_fixed_base_multi_table() {
        let table = get_fixed_base_table(
            0,
            FixedBaseParams::BITS_PER_LO_SCALAR,
            MultiTableId::FIXED_BASE_LEFT_LO,
        );
        assert_eq!(table.id, MultiTableId::FIXED_BASE_LEFT_LO);
        assert_eq!(table.slice_sizes.len(), NUM_TABLES_PER_LO);
        assert_eq!(table.basic_table_ids.len(), NUM_TABLES_PER_LO);
        assert_eq!(table.get_table_values.len(), NUM_TABLES_PER_LO);
    }
}
