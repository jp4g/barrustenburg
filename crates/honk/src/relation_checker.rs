//! Port of `relation_checker.hpp` — debugging utility for checking relation satisfaction.
//!
//! Provides a generic function to check whether a set of polynomials satisfies a given
//! relation at every row of the execution trace.

use std::collections::BTreeMap;

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use bbrs_relations::relation_parameters::RelationParameters;

/// Map from subrelation index to the first row index where it failed.
pub type FirstSubrelationFailures = BTreeMap<usize, u32>;

/// Map from relation name to its first subrelation failures.
pub type AllSubrelationFailures = BTreeMap<String, FirstSubrelationFailures>;

/// Check that a single relation is satisfied for a set of polynomial values.
///
/// Port of C++ `RelationChecker<void>::check<Relation>()`.
///
/// Parameters:
/// - `accumulate_fn`: function that accumulates the relation into `evals` given row values and params
/// - `get_row_fn`: function that extracts row `i` from the polynomials
/// - `num_subrelations`: number of subrelations
/// - `linearly_independent`: slice indicating which subrelations are linearly independent
/// - `polynomial_size`: number of rows in the execution trace
/// - `params`: relation parameters
/// - `has_linearly_dependent`: whether to handle linearly dependent subrelations specially
///
/// Returns a map of first failures per subrelation (empty if all pass).
pub fn check_relation<P, Row, F, G>(
    accumulate_fn: F,
    get_row_fn: G,
    num_subrelations: usize,
    linearly_independent: &[bool],
    polynomial_size: usize,
    params: &RelationParameters<Field<P>>,
    has_linearly_dependent: bool,
) -> FirstSubrelationFailures
where
    P: FieldParams,
    F: Fn(&mut Vec<Field<P>>, &Row, &RelationParameters<Field<P>>, &Field<P>),
    G: Fn(usize) -> Row,
{
    let mut first_failures = FirstSubrelationFailures::new();
    let mut result = vec![Field::<P>::zero(); num_subrelations];
    let one = Field::<P>::one();

    for i in 0..polynomial_size {
        let row = get_row_fn(i);
        accumulate_fn(&mut result, &row, params, &one);

        // Check linearly independent subrelations at each row
        for (sub_idx, element) in result.iter().enumerate() {
            if has_linearly_dependent {
                if !element.is_zero()
                    && sub_idx < linearly_independent.len()
                    && linearly_independent[sub_idx]
                {
                    first_failures.entry(sub_idx).or_insert(i as u32);
                }
            } else if !element.is_zero() {
                first_failures.entry(sub_idx).or_insert(i as u32);
            }
        }
    }

    // For linearly dependent subrelations, check the accumulated sum over all rows
    if has_linearly_dependent {
        for (sub_idx, element) in result.iter().enumerate() {
            if !element.is_zero()
                && sub_idx < linearly_independent.len()
                && !linearly_independent[sub_idx]
            {
                first_failures.entry(sub_idx).or_insert(0);
            }
        }
    }

    first_failures
}
