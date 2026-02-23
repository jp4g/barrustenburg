//! Circuit plookup interface.
//!
//! Port of C++ `barretenberg/stdlib/primitives/plookup/plookup.{hpp,cpp}`.
//!
//! Provides the stdlib (circuit) interface for plookup operations. The key function is
//! `get_lookup_accumulators()`, which has two code paths:
//!
//!   1. **Constant path** (builder-agnostic): If all inputs are constants, wraps native
//!      values as `FieldT` constants. No witnesses or gates are created.
//!
//!   2. **Variable path** (builder-modifying): If any input is a witness, calls
//!      `create_gates_from_plookup_accumulators()` to create actual lookup gates.

use bbrs_circuit_builder::gate_data;
use bbrs_circuit_builder::plookup_tables::plookup_tables as native;
use bbrs_circuit_builder::plookup_tables::types::{ColumnIdx, MultiTableId};
use bbrs_ecc::curves::bn254::Bn254FrParams;
use bbrs_ecc::fields::field::Field;

use super::field::FieldT;
use super::witness::{BuilderRef, IS_CONSTANT};

type P = Bn254FrParams;
type Fr = Field<P>;

/// Read a pair of values (C2, C3) from a 1-to-2 lookup table.
///
/// Port of C++ `plookup_read::read_pair_from_table`.
pub fn read_pair_from_table(
    id: MultiTableId,
    key: &FieldT<P>,
) -> (FieldT<P>, FieldT<P>) {
    let lookup = get_lookup_accumulators(id, key, &FieldT::default(), false);
    (lookup.columns[1][0].clone(), lookup.columns[2][0].clone())
}

/// Read a single value from a 2-to-1 lookup table (two keys → one output in C3).
///
/// Port of C++ `plookup_read::read_from_2_to_1_table`.
pub fn read_from_2_to_1_table(
    id: MultiTableId,
    key_a: &FieldT<P>,
    key_b: &FieldT<P>,
) -> FieldT<P> {
    let lookup = get_lookup_accumulators(id, key_a, key_b, true);
    lookup.columns[2][0].clone()
}

/// Read a single value from a 1-to-2 lookup table (returns just C2).
///
/// Port of C++ `plookup_read::read_from_1_to_2_table`.
pub fn read_from_1_to_2_table(
    id: MultiTableId,
    key_a: &FieldT<P>,
) -> FieldT<P> {
    let lookup = get_lookup_accumulators(id, key_a, &FieldT::default(), false);
    lookup.columns[1][0].clone()
}

/// Container for stdlib-level lookup read data.
///
/// Stores `FieldT` accumulator values for each of the three columns.
#[derive(Clone)]
pub struct ReadData {
    pub columns: [Vec<FieldT<P>>; 3],
}

impl ReadData {
    fn new() -> Self {
        Self {
            columns: [Vec::new(), Vec::new(), Vec::new()],
        }
    }
}

/// Get lookup accumulators from a plookup table.
///
/// This is the main entry point for stdlib code that needs table lookups. It bridges
/// the native plookup computation to circuit gate creation.
///
/// Port of C++ `plookup_read::get_lookup_accumulators`.
pub fn get_lookup_accumulators(
    id: MultiTableId,
    key_a: &FieldT<P>,
    key_b: &FieldT<P>,
    is_2_to_1_lookup: bool,
) -> ReadData {
    let ctx: BuilderRef<P> = key_a
        .get_context()
        .clone()
        .or_else(|| key_b.get_context().clone())
        .expect("at least one key must have a builder context");

    // Compute native lookup accumulators
    let native_data = native::get_lookup_accumulators(
        id,
        key_a.get_value(),
        key_b.get_value(),
        is_2_to_1_lookup,
    );

    let num_lookups = native_data.column(ColumnIdx::C1).len();
    let mut lookup = ReadData::new();

    // Constant path: both keys are constants (no gates needed)
    if key_a.is_constant() && (key_b.is_constant() || !is_2_to_1_lookup) {
        for i in 0..num_lookups {
            lookup.columns[0].push(FieldT::constant_with_context(
                ctx.clone(),
                native_data.column(ColumnIdx::C1)[i],
            ));
            lookup.columns[1].push(FieldT::constant_with_context(
                ctx.clone(),
                native_data.column(ColumnIdx::C2)[i],
            ));
            lookup.columns[2].push(FieldT::constant_with_context(
                ctx.clone(),
                native_data.column(ColumnIdx::C3)[i],
            ));
        }
        return lookup;
    }

    // Variable path: at least one key is a witness, create lookup gates
    let mut lhs_index = key_a.get_witness_index();
    let mut rhs_index = key_b.get_witness_index();

    if key_a.is_constant() {
        lhs_index = ctx.borrow_mut().put_constant_variable(key_a.get_value());
    }
    if key_b.is_constant() && is_2_to_1_lookup {
        rhs_index = ctx.borrow_mut().put_constant_variable(key_b.get_value());
    }

    let key_b_witness = if rhs_index == IS_CONSTANT {
        None
    } else {
        Some(rhs_index)
    };

    // Convert native ReadData to gate_data ReadData for the builder
    let native_multi_table = native::get_multitable(id);
    let gate_multi_table = to_gate_multi_table(native_multi_table);
    let gate_read_data = to_gate_read_data(&native_data, is_2_to_1_lookup);

    let accumulator_witnesses = ctx.borrow_mut().create_gates_from_plookup_accumulators(
        &gate_multi_table,
        &gate_read_data,
        lhs_index,
        key_b_witness,
    );

    for i in 0..num_lookups {
        lookup.columns[0].push(FieldT::from_witness_index(
            ctx.clone(),
            accumulator_witnesses[gate_data::ColumnIdx::C1][i],
        ));
        lookup.columns[1].push(FieldT::from_witness_index(
            ctx.clone(),
            accumulator_witnesses[gate_data::ColumnIdx::C2][i],
        ));
        lookup.columns[2].push(FieldT::from_witness_index(
            ctx.clone(),
            accumulator_witnesses[gate_data::ColumnIdx::C3][i],
        ));
    }

    lookup
}

/// Convert a plookup_tables MultiTable to a gate_data MultiTable.
fn to_gate_multi_table(
    src: &bbrs_circuit_builder::plookup_tables::types::MultiTable,
) -> gate_data::MultiTable {
    gate_data::MultiTable {
        id: gate_data::MultiTableId(src.id.0 as u64),
        basic_table_ids: src
            .basic_table_ids
            .iter()
            .map(|id| gate_data::BasicTableId(id.0 as u64))
            .collect(),
        column_1_step_sizes: src
            .column_1_step_sizes
            .iter()
            .map(|f| fr_to_u64(f))
            .collect(),
        column_2_step_sizes: src
            .column_2_step_sizes
            .iter()
            .map(|f| fr_to_u64(f))
            .collect(),
        column_3_step_sizes: src
            .column_3_step_sizes
            .iter()
            .map(|f| fr_to_u64(f))
            .collect(),
    }
}

/// Convert a plookup_tables ReadData<Fr> to a gate_data ReadData<Fr>.
fn to_gate_read_data(
    src: &bbrs_circuit_builder::plookup_tables::types::ReadData<Fr>,
    is_2_to_1_lookup: bool,
) -> gate_data::ReadData<Fr> {
    let mut dst = gate_data::ReadData::<Fr>::default();
    dst[gate_data::ColumnIdx::C1] = src.column(ColumnIdx::C1).clone();
    dst[gate_data::ColumnIdx::C2] = src.column(ColumnIdx::C2).clone();
    dst[gate_data::ColumnIdx::C3] = src.column(ColumnIdx::C3).clone();

    // Convert LookupEntry to 3-column Vec<Fr> format used by gate_data
    for entry in &src.lookup_entries {
        let components = entry.to_table_components(is_2_to_1_lookup);
        dst.lookup_entries
            .push(vec![components[0], components[1], components[2]]);
    }
    dst
}

/// Convert an Fr field element to u64 (extracts low limb from standard form).
fn fr_to_u64(f: &Fr) -> u64 {
    let standard = f.from_montgomery_form();
    standard.data[0]
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use std::cell::RefCell;
    use std::rc::Rc;

    fn make_builder() -> BuilderRef<P> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    #[test]
    fn test_uint32_xor() {
        let builder = make_builder();

        let left_value: u64 = 0x12345678;
        let right_value: u64 = 0xDEADBEEF;

        let left = FieldT::from_witness(builder.clone(), Fr::from(left_value));
        let right = FieldT::from_witness(builder.clone(), Fr::from(right_value));

        let lookup = get_lookup_accumulators(
            MultiTableId::UINT32_XOR,
            &left,
            &right,
            true,
        );

        let expected_xor = left_value ^ right_value;
        assert_eq!(lookup.columns[2][0].get_value(), Fr::from(expected_xor));
        assert_eq!(lookup.columns[0][0].get_value(), Fr::from(left_value));
        assert_eq!(lookup.columns[1][0].get_value(), Fr::from(right_value));

        UltraCircuitChecker::check(&mut builder.borrow_mut())
            .expect("circuit check failed");
    }

    #[test]
    fn test_uint32_and() {
        let builder = make_builder();

        let left_value: u64 = 0x12345678;
        let right_value: u64 = 0xDEADBEEF;

        let left = FieldT::from_witness(builder.clone(), Fr::from(left_value));
        let right = FieldT::from_witness(builder.clone(), Fr::from(right_value));

        let lookup = get_lookup_accumulators(
            MultiTableId::UINT32_AND,
            &left,
            &right,
            true,
        );

        let expected_and = left_value & right_value;
        assert_eq!(lookup.columns[2][0].get_value(), Fr::from(expected_and));
        assert_eq!(lookup.columns[0][0].get_value(), Fr::from(left_value));
        assert_eq!(lookup.columns[1][0].get_value(), Fr::from(right_value));

        UltraCircuitChecker::check(&mut builder.borrow_mut())
            .expect("circuit check failed");
    }

    #[test]
    fn test_read_pair_from_table() {
        let builder = make_builder();
        // Use a 1-to-2 table (SECP256K1_XLO) for read_pair_from_table.
        // Key must be in range [0, 512) for this single-lookup table.
        let key = FieldT::from_witness(builder.clone(), Fr::from(42u64));

        let (c2, c3) = read_pair_from_table(MultiTableId::SECP256K1_XLO, &key);

        // Verify we get valid field elements and the circuit is satisfied
        assert!(!c2.is_constant());
        assert!(!c3.is_constant());
        UltraCircuitChecker::check(&mut builder.borrow_mut())
            .expect("circuit check failed");
    }

    #[test]
    fn test_constant_inputs_constant_outputs() {
        let builder = make_builder();

        let left = FieldT::constant_with_context(builder.clone(), Fr::from(0x12345678u64));
        let right = FieldT::constant_with_context(builder.clone(), Fr::from(0xDEADBEEFu64));

        assert!(left.is_constant());
        assert!(right.is_constant());

        let lookup = get_lookup_accumulators(
            MultiTableId::UINT32_XOR,
            &left,
            &right,
            true,
        );

        // Result should be constant
        assert!(lookup.columns[2][0].is_constant());

        let expected = 0x12345678u64 ^ 0xDEADBEEFu64;
        assert_eq!(lookup.columns[2][0].get_value(), Fr::from(expected));
    }

    #[test]
    fn test_variable_inputs_variable_outputs() {
        let builder = make_builder();

        let left = FieldT::from_witness(builder.clone(), Fr::from(0x12345678u64));
        let right = FieldT::from_witness(builder.clone(), Fr::from(0xDEADBEEFu64));

        assert!(!left.is_constant());
        assert!(!right.is_constant());

        let lookup = get_lookup_accumulators(
            MultiTableId::UINT32_XOR,
            &left,
            &right,
            true,
        );

        // Result should NOT be constant
        assert!(!lookup.columns[2][0].is_constant());

        let expected = 0x12345678u64 ^ 0xDEADBEEFu64;
        assert_eq!(lookup.columns[2][0].get_value(), Fr::from(expected));

        UltraCircuitChecker::check(&mut builder.borrow_mut())
            .expect("circuit check failed");
    }

    #[test]
    fn test_mixed_constant_variable_inputs() {
        let builder = make_builder();

        let left = FieldT::constant_with_context(builder.clone(), Fr::from(0x12345678u64));
        let right = FieldT::from_witness(builder.clone(), Fr::from(0xDEADBEEFu64));

        assert!(left.is_constant());
        assert!(!right.is_constant());

        let lookup = get_lookup_accumulators(
            MultiTableId::UINT32_XOR,
            &left,
            &right,
            true,
        );

        // Result should NOT be constant (one input is variable)
        assert!(!lookup.columns[2][0].is_constant());

        let expected = 0x12345678u64 ^ 0xDEADBEEFu64;
        assert_eq!(lookup.columns[2][0].get_value(), Fr::from(expected));

        UltraCircuitChecker::check(&mut builder.borrow_mut())
            .expect("circuit check failed");
    }

    #[test]
    fn test_read_from_2_to_1_table() {
        let builder = make_builder();

        let left = FieldT::from_witness(builder.clone(), Fr::from(0xAABBu64));
        let right = FieldT::from_witness(builder.clone(), Fr::from(0xCCDDu64));

        let result = read_from_2_to_1_table(MultiTableId::UINT32_XOR, &left, &right);

        let expected = 0xAABBu64 ^ 0xCCDDu64;
        assert_eq!(result.get_value(), Fr::from(expected));
        UltraCircuitChecker::check(&mut builder.borrow_mut())
            .expect("circuit check failed");
    }
}
