//! Port of `generic_lookup_relation.hpp` — Generic Log-Derivative Lookup Relation.
//!
//! Provides a trait-based generic lookup mechanism using the log-derivative technique.
//! Concrete lookup relations specialize `GenericLookupSettings` to define which
//! columns participate in the lookup and how read/write terms are computed.
//!
//! The relation enforces two subrelations:
//! 1. **Inverse correctness**: `I * ∏(read_terms) * ∏(write_terms) - inverse_exists = 0`
//! 2. **Log-derivative sum**: `Σ [read_pred * 1/read_term - write_pred * read_count / write_term] = 0`

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use crate::relation_parameters::RelationParameters;

/// Settings trait for specializing the generic lookup relation.
///
/// Each concrete lookup (e.g., a specific table lookup in Ultra circuits) implements
/// this trait to define:
/// - How many read/write terms exist per row
/// - How to extract the inverse polynomial and predicates from the row
/// - How to compute the read/write terms (γ + Σ column_j · β^j)
///
/// Port of C++ `GenericLookupRelationImpl<Settings, FF>`.
pub trait GenericLookupSettings {
    /// The field parameters for this lookup.
    type P: FieldParams;

    /// The type representing all entity values at a single row.
    type AllValues;

    /// Number of read terms (lookups performed) per row.
    const READ_TERMS: usize;

    /// Number of write terms (table insertions) per row.
    const WRITE_TERMS: usize;

    /// Partial lengths for the two subrelations: [inverse_correctness, log_derivative].
    fn subrelation_partial_lengths() -> [usize; 2];

    /// Whether each subrelation is linearly independent.
    /// Default: first is independent (scaled by separator), second is dependent (summed globally).
    fn subrelation_linearly_independent() -> [bool; 2] {
        [true, false]
    }

    /// Get the value of the inverse polynomial at this row.
    fn get_inverse_polynomial(input: &Self::AllValues) -> Field<Self::P>;

    /// Compute the "inverse exists" predicate.
    /// Non-zero when at least one read or write operation is active at this row.
    fn compute_inverse_exists(input: &Self::AllValues) -> Field<Self::P>;

    /// Compute the `read_index`-th read term: γ + Σ(column_j · β^j).
    fn compute_read_term(
        input: &Self::AllValues,
        params: &RelationParameters<Field<Self::P>>,
        read_index: usize,
    ) -> Field<Self::P>;

    /// Compute the `write_index`-th write term: γ + Σ(column_j · β^j).
    fn compute_write_term(
        input: &Self::AllValues,
        params: &RelationParameters<Field<Self::P>>,
        write_index: usize,
    ) -> Field<Self::P>;

    /// Compute the predicate enabling the `read_index`-th read term.
    fn compute_read_term_predicate(
        input: &Self::AllValues,
        read_index: usize,
    ) -> Field<Self::P>;

    /// Compute the predicate enabling the `write_index`-th write term.
    fn compute_write_term_predicate(
        input: &Self::AllValues,
        write_index: usize,
    ) -> Field<Self::P>;

    /// Get the read count for the `write_index`-th write term.
    /// This is the number of times the corresponding table entry has been looked up.
    fn lookup_read_counts(
        input: &Self::AllValues,
        write_index: usize,
    ) -> Field<Self::P>;
}

/// Accumulate the generic log-derivative lookup subrelation contributions.
///
/// Port of C++ `accumulate_logderivative_lookup_subrelation_contributions`.
///
/// Given the inverse polynomial I (precomputed as I = 1 / ∏(read_terms) · ∏(write_terms)),
/// this function derives individual inverse terms and accumulates:
///
/// - Subrelation 0 (inverse correctness, scaled by `scaling_factor`):
///   `(∏(all_terms) · I - inverse_exists) · scaling_factor`
///
/// - Subrelation 1 (log-derivative, NOT scaled — linearly dependent):
///   `Σ_i read_pred_i / read_term_i - Σ_j write_pred_j · read_count_j / write_term_j`
pub fn accumulate<S: GenericLookupSettings>(
    evals: &mut [Field<S::P>; 2],
    input: &S::AllValues,
    params: &RelationParameters<Field<S::P>>,
    scaling_factor: &Field<S::P>,
) {
    let num_total_terms = S::READ_TERMS + S::WRITE_TERMS;

    let inverse = S::get_inverse_polynomial(input);

    // Compute all read and write terms
    let mut terms: Vec<Field<S::P>> = Vec::with_capacity(num_total_terms);
    for i in 0..S::READ_TERMS {
        terms.push(S::compute_read_term(input, params, i));
    }
    for i in 0..S::WRITE_TERMS {
        terms.push(S::compute_write_term(input, params, i));
    }

    // Build prefix product: denom_acc[i] = ∏_{j=0}^{i} terms[j]
    let mut denom_acc = terms.clone();
    for i in 0..num_total_terms - 1 {
        let prev = denom_acc[i];
        denom_acc[i + 1] = denom_acc[i + 1] * prev;
    }

    let inverse_exists = S::compute_inverse_exists(input);

    // Subrelation 0: inverse correctness
    // product_of_all_terms * inverse - inverse_exists (when I is correct, product * I = 1)
    evals[0] = evals[0]
        + (denom_acc[num_total_terms - 1] * inverse - inverse_exists) * *scaling_factor;

    // Derive individual inverse terms: denom_acc[i] = I · ∏_{j≠i} terms[j] = 1/terms[i]
    let mut inv_acc = inverse;
    for i in 0..num_total_terms - 1 {
        let idx = num_total_terms - 1 - i;
        denom_acc[idx] = denom_acc[idx - 1] * inv_acc;
        inv_acc = inv_acc * terms[idx];
    }
    denom_acc[0] = inv_acc;

    // Subrelation 1: log-derivative accumulation (not scaled by scaling_factor)
    // Read contributions: +read_predicate[i] * (1/read_term[i])
    for i in 0..S::READ_TERMS {
        evals[1] = evals[1] + S::compute_read_term_predicate(input, i) * denom_acc[i];
    }

    // Write contributions: -write_predicate[i] * read_count[i] * (1/write_term[i])
    for i in 0..S::WRITE_TERMS {
        let predicate = S::compute_write_term_predicate(input, i);
        let read_count = S::lookup_read_counts(input, i);
        evals[1] = evals[1] - predicate * (denom_acc[i + S::READ_TERMS] * read_count);
    }
}
