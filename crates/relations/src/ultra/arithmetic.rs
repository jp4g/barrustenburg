//! Port of `ultra_arithmetic_relation.hpp` — Ultra Arithmetic Relation.
//!
//! Two subrelations with partial lengths [6, 5].
//! Toggled by q_arith ∈ {0, 1, 2, 3}.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::input_elements::InputElements;
use crate::relation_parameters::RelationParameters;

pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 2] = [6, 5];
pub const NUM_SUBRELATIONS: usize = 2;

/// Accumulate the Ultra Arithmetic relation into the output array.
///
/// Port of C++ `ArithmeticRelationImpl::accumulate` for the verifier (FF) path.
pub fn accumulate<P: FieldParams>(
    evals: &mut [Field<P>; NUM_SUBRELATIONS],
    input: &InputElements<P>,
    _params: &RelationParameters<Field<P>>,
    scaling_factor: &Field<P>,
) {
    let w_l = input.w_l();
    let w_r = input.w_r();
    let w_o = input.w_o();
    let w_4 = input.w_4();
    let w_4_shift = input.w_4_shift();
    let w_l_shift = input.w_l_shift();

    let q_m = input.q_m();
    let q_l = input.q_l();
    let q_r = input.q_r();
    let q_o = input.q_o();
    let q_4 = input.q_4();
    let q_c = input.q_c();
    let q_arith = input.q_arith();

    let neg_half = (Field::<P>::zero() - Field::from(2u64)).invert();

    let q_arith_sub_1 = q_arith - Field::one();
    let scaled_q_arith = q_arith * *scaling_factor;

    // Subrelation 1
    {
        let tmp0 = w_r * w_l * neg_half * (q_arith - Field::from(3u64)) * q_m;
        let mut tmp1 = q_l * w_l + q_r * w_r + q_o * w_o + q_4 * w_4 + q_c;
        tmp1 = tmp1 + q_arith_sub_1 * w_4_shift;

        evals[0] = evals[0] + (tmp0 + tmp1) * scaled_q_arith;
    }

    // Subrelation 2
    {
        let tmp_0 = w_l + w_4 - w_l_shift + q_m;
        let tmp_1 = tmp_0 * (q_arith - Field::from(2u64));
        let tmp_2 = q_arith_sub_1 * scaled_q_arith;
        evals[1] = evals[1] + tmp_1 * tmp_2;
    }
}

/// Check if this relation can be skipped for the given inputs.
pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
    input.q_arith().is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    fn validate_relation(expected: &[Fr; NUM_SUBRELATIONS], input: &InputElements<P>) {
        let params = RelationParameters::<Fr>::get_random();
        let mut accum = [Fr::zero(); NUM_SUBRELATIONS];
        let one = Fr::one();
        accumulate(&mut accum, input, &params, &one);
        assert_eq!(accum, *expected);
    }

    /// Port of C++ `UltraRelationConsistency::ArithmeticRelation`.
    #[test]
    fn test_arithmetic_relation_consistency() {
        let run_test = |random_inputs: bool, q_arith_value: Fr| {
            let mut input = if random_inputs {
                InputElements::<P>::get_random()
            } else {
                InputElements::<P>::get_special()
            };

            let w_1 = input.w_l();
            let w_2 = input.w_r();
            let w_3 = input.w_o();
            let w_4 = input.w_4();
            let w_4_shift = input.w_4_shift();
            let w_1_shift = input.w_l_shift();
            let q_m = input.q_m();
            let q_l = input.q_l();
            let q_r = input.q_r();
            let q_o = input.q_o();
            let q_4 = input.q_4();
            let q_c = input.q_c();

            input.set_q_arith(q_arith_value);
            let q_arith = input.q_arith();

            let neg_half = (Fr::zero() - Fr::from(2u64)).invert();

            let mut expected = [Fr::zero(); NUM_SUBRELATIONS];

            let one = Fr::one();
            let two = Fr::from(2u64);
            let three = Fr::from(3u64);
            let six = Fr::from(6u64);

            if q_arith == one {
                expected[0] = q_m * w_2 * w_1 + q_l * w_1 + q_r * w_2 + q_o * w_3 + q_4 * w_4 + q_c;
                // contribution_2: None (stays zero)
            } else if q_arith == two {
                expected[0] = q_m * w_2 * w_1;
                expected[0] = expected[0]
                    + (q_l * w_1 + q_r * w_2 + q_o * w_3 + q_4 * w_4 + w_4_shift + q_c) * two;
                // contribution_2: None
            } else if q_arith == three {
                expected[0] = q_l * w_1 + q_r * w_2 + q_o * w_3 + q_4 * w_4 + q_c;
                expected[0] = expected[0] + w_4_shift * two;
                expected[0] = expected[0] * three;

                expected[1] = (w_1 + w_4 - w_1_shift + q_m) * six;
            } else {
                expected[0] = (q_arith - three) * (q_m * w_2 * w_1) * neg_half;
                expected[0] = expected[0] + q_l * w_1 + q_r * w_2 + q_o * w_3 + q_4 * w_4 + q_c;
                expected[0] = expected[0] + (q_arith - one) * w_4_shift;
                expected[0] = expected[0] * q_arith;

                expected[1] = w_1 + w_4 - w_1_shift + q_m;
                expected[1] = expected[1] * (q_arith - two) * (q_arith - one) * q_arith;
            }

            validate_relation(&expected, &input);
        };

        run_test(false, Fr::random_element());
        run_test(true, Fr::random_element());
        run_test(true, Fr::one());
        run_test(true, Fr::from(2u64));
        run_test(true, Fr::from(3u64));
    }
}
