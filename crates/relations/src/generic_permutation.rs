//! Port of `generic_permutation_relation.hpp` — Generic Log-Derivative Set Permutation Relation.
//!
//! Provides a trait-based generic permutation mechanism using the log-derivative technique.
//! Concrete permutation relations specialize `GenericPermutationSettings` to define which
//! columns participate in the permutation check.
//!
//! The relation enforces two subrelations:
//! 1. **Inverse correctness**: `I * read_term * write_term - inverse_exists = 0`
//! 2. **Log-derivative sum**: `Σ [read_pred / read_term - write_pred / write_term] = 0`
//!
//! Permutations always have exactly 1 read term and 1 write term (unlike lookups which
//! can have multiple).

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use crate::relation_parameters::RelationParameters;

/// Settings trait for specializing the generic permutation relation.
///
/// Each concrete permutation check implements this trait to define:
/// - How many columns per set (COLUMNS_PER_SET)
/// - How to extract the inverse polynomial and predicates from the row
/// - How to compute read/write terms (γ + Σ column_j · β^j)
///
/// Port of C++ `GenericPermutationRelationImpl<Settings, FF>`.
pub trait GenericPermutationSettings {
    /// The field parameters for this permutation.
    type P: FieldParams;

    /// The type representing all entity values at a single row.
    type AllValues;

    /// Number of columns bundled into each permutation set tuple.
    const COLUMNS_PER_SET: usize;

    /// Subrelation partial lengths. Default: [5, 5] (READ_TERMS + WRITE_TERMS + 3 = 1 + 1 + 3).
    fn subrelation_partial_lengths() -> [usize; 2] {
        [5, 5]
    }

    /// Whether each subrelation is linearly independent.
    /// Default: first is independent, second is dependent (summed globally).
    fn subrelation_linearly_independent() -> [bool; 2] {
        [true, false]
    }

    /// Get the value of the inverse polynomial at this row.
    fn get_inverse_polynomial(input: &Self::AllValues) -> Field<Self::P>;

    /// Compute the "inverse exists" predicate.
    /// Logical OR of the two set predicates: `a + b - a*b`.
    fn compute_inverse_exists(input: &Self::AllValues) -> Field<Self::P>;

    /// Compute the read term (first set): γ + Σ(column_j · β^j) for columns in set 1.
    fn compute_read_term(
        input: &Self::AllValues,
        params: &RelationParameters<Field<Self::P>>,
    ) -> Field<Self::P>;

    /// Compute the write term (second set): γ + Σ(column_j · β^j) for columns in set 2.
    fn compute_write_term(
        input: &Self::AllValues,
        params: &RelationParameters<Field<Self::P>>,
    ) -> Field<Self::P>;

    /// Compute the predicate enabling the first (read) set.
    fn compute_read_term_predicate(input: &Self::AllValues) -> Field<Self::P>;

    /// Compute the predicate enabling the second (write) set.
    fn compute_write_term_predicate(input: &Self::AllValues) -> Field<Self::P>;
}

/// Accumulate the generic log-derivative permutation subrelation contributions.
///
/// Port of C++ `accumulate_logderivative_permutation_subrelation_contributions`.
///
/// Given the inverse polynomial I (precomputed as I = 1 / (read_term · write_term)),
/// this function derives individual inverse terms and accumulates:
///
/// - Subrelation 0 (inverse correctness, scaled by `scaling_factor`):
///   `(read_term · write_term · I - inverse_exists) · scaling_factor`
///
/// - Subrelation 1 (log-derivative, NOT scaled — linearly dependent):
///   `read_pred · I · write_term - write_pred · I · read_term`
///   = `read_pred / read_term - write_pred / write_term`
pub fn accumulate<S: GenericPermutationSettings>(
    evals: &mut [Field<S::P>; 2],
    input: &S::AllValues,
    params: &RelationParameters<Field<S::P>>,
    scaling_factor: &Field<S::P>,
) {
    let inverse = S::get_inverse_polynomial(input);

    let read_term = S::compute_read_term(input, params);
    let write_term = S::compute_write_term(input, params);

    // Product of both terms
    let product = read_term * write_term;

    let inverse_exists = S::compute_inverse_exists(input);

    // Subrelation 0: inverse correctness
    // (read_term * write_term * inverse - inverse_exists) * scaling_factor
    evals[0] = evals[0] + (product * inverse - inverse_exists) * *scaling_factor;

    // Derive individual inverses from the product inverse:
    // 1/read_term  = inverse * write_term
    // 1/write_term = inverse * read_term
    let inv_read = inverse * write_term;
    let inv_write = inverse * read_term;

    // Subrelation 1: log-derivative permutation (not scaled by scaling_factor)
    // read_pred / read_term - write_pred / write_term
    evals[1] = evals[1] + S::compute_read_term_predicate(input) * inv_read;
    evals[1] = evals[1] - S::compute_write_term_predicate(input) * inv_write;
}
