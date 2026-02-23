//! In-circuit Poseidon2 hash function.
//!
//! Port of `barretenberg/stdlib/hash/poseidon2/poseidon2.hpp` and
//! `poseidon2_permutation.hpp`.
//!
//! Provides circuit-level Poseidon2 hashing using dedicated Poseidon2 gates.
//! Each round of the permutation creates a single gate in the corresponding
//! execution trace block (external or internal), matching the C++ implementation.

use bbrs_circuit_builder::gate_data::{Poseidon2ExternalGate, Poseidon2InternalGate};
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use crate::primitives::field::FieldT;
use crate::primitives::witness::BuilderRef;

const T: usize = 4;
const ROUNDS_F: usize = 8;
const ROUNDS_P: usize = 56;
const RATE: usize = 3;

// ════════════════════════════════════════════════════════════════════════
//  Native helpers (compute expected output values for witness creation)
// ════════════════════════════════════════════════════════════════════════

/// Get round constants for a given round, converted to `Field<P>`.
///
/// Uses `from_raw` since ROUND_CONSTANTS are already in Montgomery form.
/// NOTE: Only correct when `P = Bn254FrParams`.
fn get_round_constants<P: FieldParams>(round_idx: usize) -> [Field<P>; T] {
    use bbrs_crypto::poseidon2::params::ROUND_CONSTANTS;
    let rc = &(*ROUND_CONSTANTS)[round_idx];
    [
        Field::from_raw(rc[0].data),
        Field::from_raw(rc[1].data),
        Field::from_raw(rc[2].data),
        Field::from_raw(rc[3].data),
    ]
}

/// Get internal matrix diagonal values, converted to `Field<P>`.
///
/// Uses `from_raw` since INTERNAL_DIAGONAL values are already in Montgomery form.
fn get_internal_diagonal<P: FieldParams>() -> [Field<P>; T] {
    use bbrs_crypto::poseidon2::params::INTERNAL_DIAGONAL;
    let diag = &*INTERNAL_DIAGONAL;
    [
        Field::from_raw(diag[0].data),
        Field::from_raw(diag[1].data),
        Field::from_raw(diag[2].data),
        Field::from_raw(diag[3].data),
    ]
}

/// Native S-box: x -> x^5.
#[inline]
fn native_sbox<P: FieldParams>(x: Field<P>) -> Field<P> {
    let xx = x * x;
    let xxxx = xx * xx;
    xxxx * x
}

/// Native external MDS matrix multiplication.
///
/// Matrix:
/// | 5 7 1 3 |
/// | 4 6 1 1 |
/// | 1 3 5 7 |
/// | 1 1 4 6 |
fn native_matrix_mul_external<P: FieldParams>(state: &mut [Field<P>; T]) {
    let t0 = state[0] + state[1];
    let t1 = state[2] + state[3];
    let mut t2 = state[1] + state[1];
    t2 = t2 + t1;
    let mut t3 = state[3] + state[3];
    t3 = t3 + t0;
    let mut t4 = t1 + t1;
    t4 = t4 + t4;
    t4 = t4 + t3;
    let mut t5 = t0 + t0;
    t5 = t5 + t5;
    t5 = t5 + t2;
    let t6 = t3 + t5;
    let t7 = t2 + t4;
    state[0] = t6;
    state[1] = t5;
    state[2] = t7;
    state[3] = t4;
}

/// Native internal MDS matrix multiplication.
fn native_matrix_mul_internal<P: FieldParams>(state: &mut [Field<P>; T]) {
    let diag = get_internal_diagonal::<P>();
    let mut sum = state[0];
    for i in 1..T {
        sum = sum + state[i];
    }
    for i in 0..T {
        state[i] = state[i] * diag[i] + sum;
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Circuit-level MDS (for initial linear layer, uses FieldT arithmetic)
// ════════════════════════════════════════════════════════════════════════

/// Apply external MDS matrix to circuit field elements using FieldT arithmetic.
fn circuit_matrix_mul_external<P: FieldParams>(state: &mut [FieldT<P>; T]) {
    let t0 = &state[0] + &state[1];
    let t1 = &state[2] + &state[3];
    let t2 = {
        let doubled = &state[1] + &state[1];
        doubled + &t1
    };
    let t3 = {
        let doubled = &state[3] + &state[3];
        doubled + &t0
    };
    let t4 = {
        let doubled = &t1 + &t1;
        let quadrupled = &doubled + &doubled;
        quadrupled + &t3
    };
    let t5 = {
        let doubled = &t0 + &t0;
        let quadrupled = &doubled + &doubled;
        quadrupled + &t2
    };
    let t6 = &t3 + &t5;
    let t7 = &t2 + &t4;
    state[0] = t6;
    state[1] = t5;
    state[2] = t7;
    state[3] = t4;
}

// ════════════════════════════════════════════════════════════════════════
//  Gate creation helpers
// ════════════════════════════════════════════════════════════════════════

/// Create a Poseidon2 external round gate and update state to output values.
fn poseidon2_external_gate<P: FieldParams>(
    ctx: &BuilderRef<P>,
    state: &mut [FieldT<P>; T],
    round_idx: usize,
) {
    // Capture native values before normalization (which may create gates)
    let mut native_state = [
        state[0].get_value(),
        state[1].get_value(),
        state[2].get_value(),
        state[3].get_value(),
    ];

    // Normalize to get clean witness indices
    for s in state.iter_mut() {
        *s = s.normalize();
    }

    // Create the gate
    ctx.borrow_mut()
        .create_poseidon2_external_gate(&Poseidon2ExternalGate {
            a: state[0].witness_index,
            b: state[1].witness_index,
            c: state[2].witness_index,
            d: state[3].witness_index,
            round_idx,
        });

    // Compute native output: add_RC -> sbox -> external MDS
    let rc = get_round_constants::<P>(round_idx);
    for i in 0..T {
        native_state[i] = native_state[i] + rc[i];
    }
    for i in 0..T {
        native_state[i] = native_sbox(native_state[i]);
    }
    native_matrix_mul_external(&mut native_state);

    // Create new witness variables for output
    let indices = {
        let mut builder = ctx.borrow_mut();
        [
            builder.base.add_variable(native_state[0]),
            builder.base.add_variable(native_state[1]),
            builder.base.add_variable(native_state[2]),
            builder.base.add_variable(native_state[3]),
        ]
    };
    for i in 0..T {
        state[i] = FieldT::from_witness_index(ctx.clone(), indices[i]);
    }
}

/// Create a Poseidon2 internal round gate and update state to output values.
fn poseidon2_internal_gate<P: FieldParams>(
    ctx: &BuilderRef<P>,
    state: &mut [FieldT<P>; T],
    round_idx: usize,
) {
    let mut native_state = [
        state[0].get_value(),
        state[1].get_value(),
        state[2].get_value(),
        state[3].get_value(),
    ];

    for s in state.iter_mut() {
        *s = s.normalize();
    }

    ctx.borrow_mut()
        .create_poseidon2_internal_gate(&Poseidon2InternalGate {
            a: state[0].witness_index,
            b: state[1].witness_index,
            c: state[2].witness_index,
            d: state[3].witness_index,
            round_idx,
        });

    // Compute native output: add_RC[0] to first element -> sbox first element -> internal MDS
    let rc = get_round_constants::<P>(round_idx);
    native_state[0] = native_state[0] + rc[0];
    native_state[0] = native_sbox(native_state[0]);
    native_matrix_mul_internal(&mut native_state);

    let indices = {
        let mut builder = ctx.borrow_mut();
        [
            builder.base.add_variable(native_state[0]),
            builder.base.add_variable(native_state[1]),
            builder.base.add_variable(native_state[2]),
            builder.base.add_variable(native_state[3]),
        ]
    };
    for i in 0..T {
        state[i] = FieldT::from_witness_index(ctx.clone(), indices[i]);
    }
}

/// Add an end row to the external block (selector = 0, holds output state).
fn add_external_end_row<P: FieldParams>(ctx: &BuilderRef<P>, state: &[FieldT<P>; T]) {
    ctx.borrow_mut().create_poseidon2_external_end_row(
        state[0].witness_index,
        state[1].witness_index,
        state[2].witness_index,
        state[3].witness_index,
    );
}

/// Add an end row to the internal block (selector = 0, holds output state).
fn add_internal_end_row<P: FieldParams>(ctx: &BuilderRef<P>, state: &[FieldT<P>; T]) {
    ctx.borrow_mut().create_poseidon2_internal_end_row(
        state[0].witness_index,
        state[1].witness_index,
        state[2].witness_index,
        state[3].witness_index,
    );
}

// ════════════════════════════════════════════════════════════════════════
//  Public API
// ════════════════════════════════════════════════════════════════════════

/// In-circuit Poseidon2 permutation over 4 field elements.
///
/// Applies the full Poseidon2 permutation using dedicated circuit gates:
/// - Initial external MDS linear layer
/// - 4 external (full) rounds
/// - 56 internal (partial) rounds
/// - 4 external (full) rounds
///
/// Modifies `state` in place with the permuted values.
pub fn permutation<P: FieldParams>(ctx: &BuilderRef<P>, state: &mut [FieldT<P>; T]) {
    // Initial linear layer (uses FieldT arithmetic → regular arithmetic gates)
    circuit_matrix_mul_external(state);

    // Ensure all elements are witnesses (constants can't be used as gate wire indices)
    for s in state.iter_mut() {
        if s.is_constant() {
            s.convert_constant_to_fixed_witness(ctx.clone());
        }
    }

    // First 4 external rounds (rounds 0-3)
    let rounds_f_beginning = ROUNDS_F / 2;
    for i in 0..rounds_f_beginning {
        poseidon2_external_gate(ctx, state, i);
    }
    add_external_end_row(ctx, state);

    // 56 internal rounds (rounds 4-59)
    let p_end = rounds_f_beginning + ROUNDS_P;
    for i in rounds_f_beginning..p_end {
        poseidon2_internal_gate(ctx, state, i);
    }
    add_internal_end_row(ctx, state);

    // Last 4 external rounds (rounds 60-63)
    let num_rounds = ROUNDS_F + ROUNDS_P;
    for i in p_end..num_rounds {
        poseidon2_external_gate(ctx, state, i);
    }
    add_external_end_row(ctx, state);
}

/// In-circuit Poseidon2 sponge hash.
///
/// Hashes a slice of circuit field elements into a single field element using
/// the Poseidon2 sponge construction with rate=3, capacity=1.
///
/// The capacity element is initialized with `IV = input_length << 64`.
pub fn hash<P: FieldParams>(ctx: BuilderRef<P>, input: &[FieldT<P>]) -> FieldT<P> {
    // IV = input_length << 64
    let iv = Field::<P>::from_limbs([0, input.len() as u64, 0, 0]);

    let mut state: [FieldT<P>; T] = [
        FieldT::with_context(ctx.clone()),
        FieldT::with_context(ctx.clone()),
        FieldT::with_context(ctx.clone()),
        FieldT::constant_with_context(ctx.clone(), iv),
    ];

    let mut cache: [FieldT<P>; RATE] = [
        FieldT::with_context(ctx.clone()),
        FieldT::with_context(ctx.clone()),
        FieldT::with_context(ctx.clone()),
    ];
    let mut cache_size: usize = 0;

    // Absorb all input elements
    for elem in input {
        if cache_size == RATE {
            // Cache full: duplex (absorb cache into state, then permute)
            for i in 0..RATE {
                state[i] = state[i].clone() + cache[i].clone();
            }
            permutation(&ctx, &mut state);
            cache = [
                FieldT::with_context(ctx.clone()),
                FieldT::with_context(ctx.clone()),
                FieldT::with_context(ctx.clone()),
            ];
            cache[0] = elem.clone();
            cache_size = 1;
        } else {
            cache[cache_size] = elem.clone();
            cache_size += 1;
        }
    }

    // Squeeze: final duplex
    for i in 0..RATE {
        state[i] = state[i].clone() + cache[i].clone();
    }
    permutation(&ctx, &mut state);

    state[0].normalize()
}

// ════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use std::rc::Rc;

    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type FrField = FieldT<Bn254FrParams>;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    fn check_circuit(builder: &BuilderRef<Bn254FrParams>) -> Result<(), String> {
        UltraCircuitChecker::check(&mut builder.borrow_mut())
    }

    /// Test that in-circuit permutation produces the same output as native.
    #[test]
    fn test_permutation_matches_native() {
        let builder = make_builder();

        let input = [
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
        ];

        // In-circuit
        let mut circuit_state: [FrField; T] = [
            FrField::from_witness(builder.clone(), input[0]),
            FrField::from_witness(builder.clone(), input[1]),
            FrField::from_witness(builder.clone(), input[2]),
            FrField::from_witness(builder.clone(), input[3]),
        ];
        permutation(&builder, &mut circuit_state);

        // Native
        let native_output = bbrs_crypto::poseidon2::permutation::permutation(&input);

        for i in 0..T {
            assert_eq!(
                circuit_state[i].get_value(),
                native_output[i],
                "permutation output mismatch at index {}",
                i
            );
        }

        assert!(check_circuit(&builder).is_ok());
    }

    /// Test permutation against known test vector [0, 1, 2, 3].
    #[test]
    fn test_permutation_test_vector() {
        let builder = make_builder();

        let input = *bbrs_crypto::poseidon2::params::TEST_VECTOR_INPUT;
        let expected = *bbrs_crypto::poseidon2::params::TEST_VECTOR_OUTPUT;

        let mut circuit_state: [FrField; T] = [
            FrField::from_witness(builder.clone(), input[0]),
            FrField::from_witness(builder.clone(), input[1]),
            FrField::from_witness(builder.clone(), input[2]),
            FrField::from_witness(builder.clone(), input[3]),
        ];
        permutation(&builder, &mut circuit_state);

        for i in 0..T {
            assert_eq!(
                circuit_state[i].get_value(),
                expected[i],
                "test vector mismatch at index {}",
                i
            );
        }

        assert!(check_circuit(&builder).is_ok());
    }

    /// Test in-circuit hash matches native hash for a single element.
    #[test]
    fn test_hash_single_element() {
        let builder = make_builder();
        let val = Fr::from(42u64);
        let input = [FrField::from_witness(builder.clone(), val)];

        let result = hash(builder.clone(), &input);
        let native_result = bbrs_crypto::poseidon2::sponge::hash(&[val]);

        assert_eq!(
            result.get_value(),
            native_result,
        );
        assert!(check_circuit(&builder).is_ok());
    }

    /// Test hash of empty input.
    #[test]
    fn test_hash_empty_input() {
        let builder = make_builder();
        let input: &[FrField] = &[];

        let result = hash(builder.clone(), input);
        let native_result = bbrs_crypto::poseidon2::sponge::hash(&[]);

        assert_eq!(
            result.get_value(),
            native_result,
        );
        assert!(check_circuit(&builder).is_ok());
    }

    /// Test hash of two elements.
    #[test]
    fn test_hash_two_elements() {
        let builder = make_builder();
        let vals = [Fr::from(7u64), Fr::from(13u64)];
        let input: Vec<FrField> = vals
            .iter()
            .map(|v| FrField::from_witness(builder.clone(), *v))
            .collect();

        let result = hash(builder.clone(), &input);
        let native_result = bbrs_crypto::poseidon2::sponge::hash(&vals);

        assert_eq!(
            result.get_value(),
            native_result,
        );
        assert!(check_circuit(&builder).is_ok());
    }

    /// Test hash with exactly RATE (3) elements.
    #[test]
    fn test_hash_rate_elements() {
        let builder = make_builder();
        let vals = [Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let input: Vec<FrField> = vals
            .iter()
            .map(|v| FrField::from_witness(builder.clone(), *v))
            .collect();

        let result = hash(builder.clone(), &input);
        let native_result = bbrs_crypto::poseidon2::sponge::hash(&vals);

        assert_eq!(
            result.get_value(),
            native_result,
        );
        assert!(check_circuit(&builder).is_ok());
    }

    /// Test hash with more than RATE elements (triggers multi-block absorption).
    #[test]
    fn test_hash_many_elements() {
        let builder = make_builder();
        let vals: Vec<Fr> = (1..=7).map(|i| Fr::from(i as u64)).collect();
        let input: Vec<FrField> = vals
            .iter()
            .map(|v| FrField::from_witness(builder.clone(), *v))
            .collect();

        let result = hash(builder.clone(), &input);
        let native_result = bbrs_crypto::poseidon2::sponge::hash(&vals);

        assert_eq!(
            result.get_value(),
            native_result,
        );
        assert!(check_circuit(&builder).is_ok());
    }

    /// Test hash with known consistency vector from C++ tests.
    #[test]
    fn test_hash_consistency_vector() {
        let builder = make_builder();
        let val = Fr::from_limbs([
            0xabcdef0123456789,
            0xfedcba9876543210,
            0xa0b1c2d3e4f56789,
            0x9a807b615c4d3e2f,
        ]);
        let vals = [val, val, val, val];
        let input: Vec<FrField> = vals
            .iter()
            .map(|v| FrField::from_witness(builder.clone(), *v))
            .collect();

        let result = hash(builder.clone(), &input);
        let native_result = bbrs_crypto::poseidon2::sponge::hash(&vals);

        assert_eq!(
            result.get_value(),
            native_result,
        );
        assert!(check_circuit(&builder).is_ok());
    }

    /// Test that different inputs produce different hashes.
    #[test]
    fn test_hash_different_inputs() {
        let builder = make_builder();
        let vals1 = [Fr::from(1u64), Fr::from(2u64)];
        let vals2 = [Fr::from(2u64), Fr::from(1u64)];

        let input1: Vec<FrField> = vals1
            .iter()
            .map(|v| FrField::from_witness(builder.clone(), *v))
            .collect();
        let input2: Vec<FrField> = vals2
            .iter()
            .map(|v| FrField::from_witness(builder.clone(), *v))
            .collect();

        let r1 = hash(builder.clone(), &input1);
        let r2 = hash(builder.clone(), &input2);

        assert_ne!(
            r1.get_value(),
            r2.get_value(),
            "different input order must produce different hash"
        );
        assert!(check_circuit(&builder).is_ok());
    }

    /// Test circuit failure: assert wrong output value.
    #[test]
    fn test_circuit_failure_wrong_output() {
        let builder = make_builder();
        let vals = [Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let input: Vec<FrField> = vals
            .iter()
            .map(|v| FrField::from_witness(builder.clone(), *v))
            .collect();

        let result = hash(builder.clone(), &input);

        // Assert the result equals a deliberately wrong value
        let wrong_value = FrField::from_witness(builder.clone(), Fr::from(999u64));
        result.assert_equal(&wrong_value, "deliberately wrong");

        // Builder should detect the value mismatch
        assert!(builder.borrow().base.failed());
    }

    /// Test circuit failure: wrong permutation output.
    #[test]
    fn test_circuit_failure_wrong_permutation() {
        let builder = make_builder();

        let input = [
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
        ];
        let mut circuit_state: [FrField; T] = [
            FrField::from_witness(builder.clone(), input[0]),
            FrField::from_witness(builder.clone(), input[1]),
            FrField::from_witness(builder.clone(), input[2]),
            FrField::from_witness(builder.clone(), input[3]),
        ];
        permutation(&builder, &mut circuit_state);

        // Assert output[0] equals a wrong constant
        let wrong = FrField::from_witness(builder.clone(), Fr::from(0u64));
        circuit_state[0].assert_equal(&wrong, "deliberately wrong");

        // Builder should detect the value mismatch
        assert!(builder.borrow().base.failed());
    }
}
