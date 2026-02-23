//! Port of `logderiv_lookup_relation.hpp` — Log-Derivative Lookup Relation.
//!
//! Three subrelations with partial lengths [5, 5, 3].
//! Linearly independent: [true, false, true].

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::input_elements::InputElements;
use crate::relation_parameters::RelationParameters;

pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 3] = [5, 5, 3];
pub const NUM_SUBRELATIONS: usize = 3;

/// Compute the write term for the lookup relation.
/// write_term = 1 / (table_1 + gamma + table_2*eta + table_3*eta^2 + table_4*eta^3)
fn compute_write_term<P: FieldParams>(
    input: &InputElements<P>,
    params: &RelationParameters<Field<P>>,
) -> Field<P> {
    let table_1 = input.table_1();
    let table_2 = input.table_2();
    let table_3 = input.table_3();
    let table_4 = input.table_4();

    let mut write_term = table_1 + params.gamma;
    write_term = write_term + table_2 * params.eta;
    write_term = write_term + table_3 * params.eta_two;
    write_term = write_term + table_4 * params.eta_three;

    write_term.invert()
}

/// Compute the read term for the lookup relation.
///
/// The wire values for lookup gates are accumulators structured so that the differences
/// w_i + step_size*w_i_shift yield entries present in column i of the corresponding table.
///
/// read_term = 1 / (derived_entry_1 + derived_entry_2*eta + derived_entry_3*eta^2 + table_index*eta^3)
///
/// where:
///   derived_entry_1 = w_1 + q_r*w_1_shift + gamma
///   derived_entry_2 = w_2 + q_m*w_2_shift
///   derived_entry_3 = w_3 + q_c*w_3_shift
///   table_index     = q_o
fn compute_read_term<P: FieldParams>(
    input: &InputElements<P>,
    params: &RelationParameters<Field<P>>,
) -> Field<P> {
    let w_1 = input.w_l();
    let w_2 = input.w_r();
    let w_3 = input.w_o();
    let w_1_shift = input.w_l_shift();
    let w_2_shift = input.w_r_shift();
    let w_3_shift = input.w_o_shift();

    let negative_column_1_step_size = input.q_r();
    let negative_column_2_step_size = input.q_m();
    let negative_column_3_step_size = input.q_c();
    let table_index = input.q_o();

    let derived_entry_1 = w_1 + negative_column_1_step_size * w_1_shift + params.gamma;
    let derived_entry_2 = w_2 + negative_column_2_step_size * w_2_shift;
    let derived_entry_3 = w_3 + negative_column_3_step_size * w_3_shift;

    let mut read_term = derived_entry_1;
    read_term = read_term + derived_entry_2 * params.eta;
    read_term = read_term + derived_entry_3 * params.eta_two;
    read_term = read_term + table_index * params.eta_three;

    read_term.invert()
}

/// Accumulate the LogDerivLookup relation into the output array.
///
/// Port of C++ `LogDerivLookupRelationImpl::accumulate` for the verifier (FF) path.
pub fn accumulate<P: FieldParams>(
    evals: &mut [Field<P>; NUM_SUBRELATIONS],
    input: &InputElements<P>,
    params: &RelationParameters<Field<P>>,
    scaling_factor: &Field<P>,
) {
    let q_lookup = input.q_lookup();
    let lookup_read_counts = input.lookup_read_counts();
    let lookup_read_tags = input.lookup_read_tags();
    let lookup_inverses = input.lookup_inverses();

    let write_term = compute_write_term(input, params);
    let read_term = compute_read_term(input, params);

    // Determine whether the inverse is well-defined
    let inverse_exists = lookup_inverses * write_term * read_term;

    // Subrelation 0 (length 5, linearly independent):
    // (read_term * write_term * lookup_inverses - inverse_exists) * scaling_factor
    evals[0] = evals[0]
        + (read_term * write_term * lookup_inverses - inverse_exists) * *scaling_factor;

    // Subrelation 1 (length 5, NOT linearly independent):
    // (q_lookup * write_term - lookup_read_counts * read_term) * lookup_inverses
    // Note: no scaling_factor for non-linearly-independent subrelations
    evals[1] = evals[1]
        + (q_lookup * write_term - lookup_read_counts * read_term) * lookup_inverses;

    // Subrelation 2 (length 3, linearly independent):
    // (lookup_read_tags^2 - lookup_read_tags) * scaling_factor
    evals[2] = evals[2]
        + (lookup_read_tags.sqr() - lookup_read_tags) * *scaling_factor;
}

/// Check if this relation can be skipped for the given inputs.
pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
    input.q_lookup().is_zero() && input.lookup_read_counts().is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    fn validate_relation(input: &InputElements<P>) {
        let params = RelationParameters::<Fr>::get_random();
        let mut accum = [Fr::zero(); NUM_SUBRELATIONS];
        let one = Fr::one();
        accumulate(&mut accum, input, &params, &one);

        // Compute expected independently
        let q_lookup = input.q_lookup();
        let lookup_read_counts = input.lookup_read_counts();
        let lookup_read_tags = input.lookup_read_tags();
        let lookup_inverses = input.lookup_inverses();

        let write_term = compute_write_term(input, &params);
        let read_term = compute_read_term(input, &params);
        let inverse_exists = lookup_inverses * write_term * read_term;

        let mut exp = [Fr::zero(); NUM_SUBRELATIONS];
        exp[0] = read_term * write_term * lookup_inverses - inverse_exists;
        exp[1] = (q_lookup * write_term - lookup_read_counts * read_term) * lookup_inverses;
        exp[2] = lookup_read_tags.sqr() - lookup_read_tags;

        assert_eq!(accum, exp, "accumulate should match hand-computed expected values");
    }

    /// Port of C++ `UltraRelationConsistency::LogDerivLookupRelation`.
    #[test]
    fn test_logderiv_lookup_relation_consistency() {
        // Test with random inputs
        let input = InputElements::<P>::get_random();
        validate_relation(&input);

        // Test with special (deterministic) inputs
        let input_special = InputElements::<P>::get_special();
        validate_relation(&input_special);
    }

    #[test]
    fn test_logderiv_lookup_skip() {
        // When q_lookup = 0 and lookup_read_counts = 0, skip should return true
        let mut input = InputElements::<P>::get_random();
        input.data[11] = Fr::zero(); // q_lookup
        input.data[42] = Fr::zero(); // lookup_read_counts
        assert!(skip(&input), "skip should be true when q_lookup=0 and read_counts=0");

        // When q_lookup != 0, skip should return false
        let mut input2 = InputElements::<P>::get_random();
        input2.data[11] = Fr::one();
        assert!(!skip(&input2), "skip should be false when q_lookup != 0");

        // When read_counts != 0, skip should return false
        let mut input3 = InputElements::<P>::get_random();
        input3.data[11] = Fr::zero();
        input3.data[42] = Fr::one();
        assert!(!skip(&input3), "skip should be false when read_counts != 0");
    }
}
