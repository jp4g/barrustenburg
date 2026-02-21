//! Port of `relation_types.hpp` — Relation wrapper trait and concepts.
//!
//! The C++ code uses a `Relation<RelationImpl>` wrapper class that computes
//! `RELATION_LENGTH` from the max of `SUBRELATION_PARTIAL_LENGTHS` and defines
//! associated type aliases for sumcheck containers.
//!
//! In Rust we use a trait hierarchy: `RelationImpl` defines the core relation
//! logic, and `Relation` is a supertrait that adds computed constants.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

/// Core trait that each relation must implement.
///
/// Port of the C++ `RelationImpl` pattern where each relation defines:
/// - `SUBRELATION_PARTIAL_LENGTHS`: array of lengths for each subrelation
/// - `SUBRELATION_LINEARLY_INDEPENDENT`: optional array indicating independence
/// - `accumulate()`: the relation arithmetic
/// - `skip()`: optional fast-skip check
pub trait RelationImpl {
    /// The field type used by this relation.
    type P: FieldParams;

    /// Number of subrelations in this relation.
    const NUM_SUBRELATIONS: usize;

    /// Partial lengths for each subrelation (1 + degree of each subrelation polynomial).
    /// Length = NUM_SUBRELATIONS.
    fn subrelation_partial_lengths() -> &'static [usize];

    /// Whether each subrelation is linearly independent from the others.
    /// Default: all true (each gets scaled by its own power of the separator challenge).
    /// When false for a subrelation, it "merges" with another and doesn't need separate scaling.
    fn subrelation_linearly_independent() -> &'static [bool] {
        // Default: all subrelations are linearly independent
        // We return a static slice of `true` values. Each concrete impl can override.
        &[true; 32][..Self::NUM_SUBRELATIONS]
    }

    /// Accumulate the contribution of this relation into `evals`.
    ///
    /// - `evals`: container of accumulators (one per subrelation), either Univariates or field values
    /// - `in_values`: the wire/selector values at the current row (AllEntities pattern)
    /// - `params`: relation parameters (eta, beta, gamma, etc.)
    /// - `scaling_factor`: multiplier for this row's contribution
    fn accumulate<Acc: AccumulatorTypes<Self::P>>(
        evals: &mut Acc::Accumulators,
        in_values: &Acc::AllValues,
        params: &RelationParameters<Field<Self::P>>,
        scaling_factor: &Field<Self::P>,
    );

    /// Returns true if the contribution from all subrelations for the provided inputs
    /// is identically zero, allowing the caller to skip accumulation.
    ///
    /// Default: never skip.
    fn skip<T>(_in_values: &T) -> bool {
        false
    }
}

/// Marker trait for types that can serve as accumulator containers.
///
/// In the C++ code, this is the split between prover (Univariate accumulators)
/// and verifier (FF value accumulators). The prover accumulates `Univariate<FF, N>`
/// for each subrelation; the verifier accumulates plain `FF` values.
pub trait AccumulatorTypes<P: FieldParams> {
    /// The container of accumulators — one per subrelation.
    /// For prover: tuple of Univariates with different lengths.
    /// For verifier: array of FF values.
    type Accumulators;

    /// The type providing wire/selector values for a single row.
    type AllValues;
}

/// Computed trait that mirrors `Relation<RelationImpl>` in C++.
///
/// Adds `RELATION_LENGTH` (max of subrelation partial lengths) and
/// the sumcheck container type aliases.
pub trait Relation: RelationImpl {
    /// Maximum partial length across all subrelations.
    /// In C++: `*std::max_element(SUBRELATION_PARTIAL_LENGTHS.begin(), ...end())`
    fn relation_length() -> usize {
        *Self::subrelation_partial_lengths()
            .iter()
            .max()
            .unwrap_or(&0)
    }
}

// Blanket impl: every RelationImpl is also a Relation.
impl<R: RelationImpl> Relation for R {}

/// Check whether a given subrelation is linearly independent.
///
/// Port of C++ `subrelation_is_linearly_independent<Relation, subrelation_index>()`.
#[inline]
pub fn subrelation_is_linearly_independent<R: RelationImpl>(subrelation_index: usize) -> bool {
    let independent = R::subrelation_linearly_independent();
    if subrelation_index < independent.len() {
        independent[subrelation_index]
    } else {
        true // default: independent
    }
}

use crate::relation_parameters::RelationParameters;

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};
    use bbrs_polynomials::univariate::Univariate;

    // ---- Test: CreateSumcheckTupleOfTuplesOfUnivariates ----
    // Port of C++ `RelationTypes::CreateSumcheckTupleOfTuplesOfUnivariates`
    //
    // The C++ test creates a tuple-of-tuples where:
    //   Relation1 has SUBRELATION_PARTIAL_LENGTHS = [3]  → one Univariate<FF,3>
    //   Relation2 has SUBRELATION_PARTIAL_LENGTHS = [2,5] → (Univariate<FF,2>, Univariate<FF,5>)
    // And checks they all initialize to zero.
    //
    // In Rust we represent this with concrete types since we don't have variadic tuples.
    // The nested_containers module provides the type-construction pattern.

    #[test]
    fn test_create_sumcheck_tuple_of_tuples_of_univariates() {
        // Relation1: single subrelation with partial length 3
        type Univ3 = Univariate<Bn254FrParams, 3>;
        // Relation2: two subrelations with partial lengths 2, 5
        type Univ2 = Univariate<Bn254FrParams, 2>;
        type Univ5 = Univariate<Bn254FrParams, 5>;

        // Create the "tuple of tuples" as nested Vecs/tuples
        // Relation1's accumulators: (Univariate<FF,3>,)
        let rel1_accum: (Univ3,) = (Univ3::zero(),);
        // Relation2's accumulators: (Univariate<FF,2>, Univariate<FF,5>)
        let rel2_accum: (Univ2, Univ5) = (Univ2::zero(), Univ5::zero());

        // Check all are zero-initialized
        assert_eq!(rel1_accum.0, Univ3::zero());
        assert_eq!(rel2_accum.0, Univ2::zero());
        assert_eq!(rel2_accum.1, Univ5::zero());

        // Test creating from type (value-initialization)
        let rel1_from_type: (Univ3,) = Default::default();
        let rel2_from_type: (Univ2, Univ5) = Default::default();

        assert_eq!(rel1_from_type.0, Univ3::zero());
        assert_eq!(rel2_from_type.0, Univ2::zero());
        assert_eq!(rel2_from_type.1, Univ5::zero());
    }

    // ---- Test: IsSkippableConcept ----
    // Port of C++ `RelationTypes::IsSkippableConcept`
    //
    // Tests the `isSkippable` concept which checks if a relation has a `skip` method.
    // In Rust we test this via the trait's default vs overridden `skip()`.

    struct SkippableInputs {
        _input: i32,
    }

    // A relation that overrides skip (always returns false)
    struct Relation1Skippable;
    impl RelationImpl for Relation1Skippable {
        type P = Bn254FrParams;
        const NUM_SUBRELATIONS: usize = 1;

        fn subrelation_partial_lengths() -> &'static [usize] {
            &[3]
        }

        fn accumulate<Acc: AccumulatorTypes<Self::P>>(
            _evals: &mut Acc::Accumulators,
            _in_values: &Acc::AllValues,
            _params: &RelationParameters<Fr>,
            _scaling_factor: &Fr,
        ) {
        }

        fn skip<T>(_in_values: &T) -> bool {
            false
        }
    }

    // A relation that does NOT override skip (uses default, which returns false)
    struct Relation2NoSkip;
    impl RelationImpl for Relation2NoSkip {
        type P = Bn254FrParams;
        const NUM_SUBRELATIONS: usize = 1;

        fn subrelation_partial_lengths() -> &'static [usize] {
            &[3]
        }

        fn accumulate<Acc: AccumulatorTypes<Self::P>>(
            _evals: &mut Acc::Accumulators,
            _in_values: &Acc::AllValues,
            _params: &RelationParameters<Fr>,
            _scaling_factor: &Fr,
        ) {
        }
        // skip is NOT overridden — uses default
    }

    #[test]
    fn test_is_skippable_concept() {
        let inputs = SkippableInputs { _input: 42 };

        // Relation1 has skip → returns false
        assert!(!Relation1Skippable::skip(&inputs));

        // Relation2 uses default skip → also returns false
        assert!(!Relation2NoSkip::skip(&inputs));

        // In C++, the concept checks compile-time whether skip *exists*.
        // In Rust, every RelationImpl has skip (via default impl), so the
        // distinction is about whether the relation has meaningful skip logic.
        // The test verifies the trait mechanism works correctly.

        // Test that subrelation_is_linearly_independent defaults to true
        assert!(subrelation_is_linearly_independent::<Relation1Skippable>(0));
        assert!(subrelation_is_linearly_independent::<Relation2NoSkip>(0));
    }
}
