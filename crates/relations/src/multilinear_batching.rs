//! Port of `multilinear_batching_relation.hpp` — Multilinear Batching Relations for HyperNova.
//!
//! Two relation implementations for batching polynomial evaluation claims:
//! - `accumulator`: accumulator contribution (batched_unshifted_accumulator * eq_accumulator, etc.)
//! - `instance`: instance contribution (batched_unshifted_instance * eq_instance, etc.)

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use crate::relation_parameters::RelationParameters;

/// Input elements for multilinear batching relations.
///
/// Contains 6 field elements representing the batched polynomials and eq selectors
/// for both accumulator and instance contributions.
pub struct InputElements<P: FieldParams> {
    pub batched_unshifted_accumulator: Field<P>,
    pub batched_unshifted_instance: Field<P>,
    pub eq_accumulator: Field<P>,
    pub eq_instance: Field<P>,
    pub batched_shifted_accumulator: Field<P>,
    pub batched_shifted_instance: Field<P>,
}

impl<P: FieldParams> InputElements<P> {
    /// Deterministic test values: 1, 2, 3, 4, 5, 6.
    pub fn get_special() -> Self {
        Self {
            batched_unshifted_accumulator: Field::from(1u64),
            batched_unshifted_instance: Field::from(2u64),
            eq_accumulator: Field::from(3u64),
            eq_instance: Field::from(4u64),
            batched_shifted_accumulator: Field::from(5u64),
            batched_shifted_instance: Field::from(6u64),
        }
    }

    /// Random test values.
    pub fn get_random() -> Self {
        Self {
            batched_unshifted_accumulator: Field::random_element(),
            batched_unshifted_instance: Field::random_element(),
            eq_accumulator: Field::random_element(),
            eq_instance: Field::random_element(),
            batched_shifted_accumulator: Field::random_element(),
            batched_shifted_instance: Field::random_element(),
        }
    }
}

/// Accumulator contribution to the multilinear batching sumcheck.
///
/// Subrelation 0: batched_unshifted_accumulator * eq_accumulator (non-shifted)
/// Subrelation 1: batched_shifted_accumulator * eq_accumulator (shifted)
///
/// Both subrelations are linearly dependent (not scaled by separator challenges).
pub mod accumulator {
    use super::*;

    pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 2] = [3, 3];
    pub const SUBRELATION_LINEARLY_INDEPENDENT: [bool; 2] = [false, false];
    pub const NUM_SUBRELATIONS: usize = 2;

    pub fn accumulate<P: FieldParams>(
        evals: &mut [Field<P>; NUM_SUBRELATIONS],
        input: &InputElements<P>,
        _params: &RelationParameters<Field<P>>,
        _scaling_factor: &Field<P>,
    ) {
        evals[0] = evals[0] + input.batched_unshifted_accumulator * input.eq_accumulator;
        evals[1] = evals[1] + input.batched_shifted_accumulator * input.eq_accumulator;
    }

    pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
        (input.batched_unshifted_accumulator.is_zero()
            && input.batched_shifted_accumulator.is_zero())
            || input.eq_accumulator.is_zero()
    }
}

/// Instance contribution to the multilinear batching sumcheck.
///
/// Subrelation 0: batched_unshifted_instance * eq_instance (non-shifted)
/// Subrelation 1: batched_shifted_instance * eq_instance (shifted)
///
/// Both subrelations are linearly dependent (not scaled by separator challenges).
pub mod instance {
    use super::*;

    pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 2] = [3, 3];
    pub const SUBRELATION_LINEARLY_INDEPENDENT: [bool; 2] = [false, false];
    pub const NUM_SUBRELATIONS: usize = 2;

    pub fn accumulate<P: FieldParams>(
        evals: &mut [Field<P>; NUM_SUBRELATIONS],
        input: &InputElements<P>,
        _params: &RelationParameters<Field<P>>,
        _scaling_factor: &Field<P>,
    ) {
        evals[0] = evals[0] + input.batched_unshifted_instance * input.eq_instance;
        evals[1] = evals[1] + input.batched_shifted_instance * input.eq_instance;
    }

    pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
        (input.batched_unshifted_accumulator.is_zero()
            && input.batched_unshifted_instance.is_zero()
            && input.batched_shifted_accumulator.is_zero()
            && input.batched_shifted_instance.is_zero())
            || (input.eq_accumulator.is_zero() && input.eq_instance.is_zero())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_multilinear_batching_accumulator_consistency() {
        let run_case = |input: &InputElements<P>, seed: [Fr; 2]| {
            let mut accum = seed;
            let mut expected = seed;

            expected[0] =
                expected[0] + input.batched_unshifted_accumulator * input.eq_accumulator;
            expected[1] =
                expected[1] + input.batched_shifted_accumulator * input.eq_accumulator;

            let params = RelationParameters::<Fr>::new();
            let scaling = Fr::one();
            accumulator::accumulate(&mut accum, input, &params, &scaling);

            assert_eq!(accum, expected);
        };

        // Zero seed with special inputs
        run_case(
            &InputElements::<P>::get_special(),
            [Fr::zero(), Fr::zero()],
        );

        // Random seed with random inputs
        run_case(
            &InputElements::<P>::get_random(),
            [Fr::random_element(), Fr::random_element()],
        );
    }

    #[test]
    fn test_multilinear_batching_instance_consistency() {
        let run_case = |input: &InputElements<P>, seed: [Fr; 2]| {
            let mut accum = seed;
            let mut expected = seed;

            expected[0] = expected[0] + input.batched_unshifted_instance * input.eq_instance;
            expected[1] = expected[1] + input.batched_shifted_instance * input.eq_instance;

            let params = RelationParameters::<Fr>::new();
            let scaling = Fr::one();
            instance::accumulate(&mut accum, input, &params, &scaling);

            assert_eq!(accum, expected);
        };

        // Zero seed with special inputs
        run_case(
            &InputElements::<P>::get_special(),
            [Fr::zero(), Fr::zero()],
        );

        // Random seed with random inputs
        run_case(
            &InputElements::<P>::get_random(),
            [Fr::random_element(), Fr::random_element()],
        );
    }

    #[test]
    fn test_multilinear_batching_accumulator_skip() {
        // Case 1: eq_accumulator is zero → should skip
        let mut input = InputElements::<P>::get_random();
        input.eq_accumulator = Fr::zero();
        assert!(accumulator::skip(&input));

        // Case 2: both batched accumulators are zero → should skip
        let mut input = InputElements::<P>::get_random();
        input.batched_unshifted_accumulator = Fr::zero();
        input.batched_shifted_accumulator = Fr::zero();
        assert!(accumulator::skip(&input));

        // Case 3: non-zero batched and eq → should not skip
        let mut input = InputElements::<P>::get_random();
        input.batched_unshifted_accumulator = Fr::from(1u64);
        input.eq_accumulator = Fr::from(1u64);
        assert!(!accumulator::skip(&input));
    }

    #[test]
    fn test_multilinear_batching_instance_skip() {
        // Case 1: both eq values are zero → should skip
        let mut input = InputElements::<P>::get_random();
        input.eq_accumulator = Fr::zero();
        input.eq_instance = Fr::zero();
        assert!(instance::skip(&input));

        // Case 2: all batched values are zero → should skip
        let mut input = InputElements::<P>::get_random();
        input.batched_unshifted_accumulator = Fr::zero();
        input.batched_unshifted_instance = Fr::zero();
        input.batched_shifted_accumulator = Fr::zero();
        input.batched_shifted_instance = Fr::zero();
        assert!(instance::skip(&input));

        // Case 3: eq_accumulator zero but eq_instance non-zero → should NOT skip
        let mut input = InputElements::<P>::get_random();
        input.eq_accumulator = Fr::zero();
        input.eq_instance = Fr::from(1u64);
        assert!(!instance::skip(&input));

        // Case 4: eq_instance zero but eq_accumulator non-zero → should NOT skip
        let mut input = InputElements::<P>::get_random();
        input.eq_accumulator = Fr::from(1u64);
        input.eq_instance = Fr::zero();
        assert!(!instance::skip(&input));

        // Case 5: all non-zero → should NOT skip
        let input = InputElements {
            batched_unshifted_accumulator: Fr::from(1u64),
            batched_unshifted_instance: Fr::from(1u64),
            eq_accumulator: Fr::from(1u64),
            eq_instance: Fr::from(1u64),
            batched_shifted_accumulator: Fr::from(1u64),
            batched_shifted_instance: Fr::from(1u64),
        };
        assert!(!instance::skip(&input));
    }
}
