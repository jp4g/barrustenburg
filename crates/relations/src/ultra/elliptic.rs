//! Port of `elliptic_relation.hpp` — Elliptic Curve Point Addition/Doubling Relation.
//!
//! Two subrelations with partial lengths [6, 6].
//! Handles both point addition and point doubling via q_is_double selector.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_ecc::curves::bn254::Bn254FrParams;

use super::input_elements::InputElements;
use crate::relation_parameters::RelationParameters;

pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 2] = [6, 6];
pub const NUM_SUBRELATIONS: usize = 2;

/// Get the curve b parameter for short Weierstrass: y^2 = x^3 + b.
///
/// For BN254 G1: b = 3, for Grumpkin: b = -17.
/// Selects based on modulus comparison (matching C++ behavior).
fn get_curve_b<P: FieldParams>() -> Field<P> {
    if P::MODULUS == Bn254FrParams::MODULUS {
        // Grumpkin curve: b = -17
        Field::zero() - Field::from(17u64)
    } else {
        // BN254 G1 curve: b = 3
        Field::from(3u64)
    }
}

/// Accumulate the Elliptic relation.
pub fn accumulate<P: FieldParams>(
    evals: &mut [Field<P>; NUM_SUBRELATIONS],
    input: &InputElements<P>,
    _params: &RelationParameters<Field<P>>,
    scaling_factor: &Field<P>,
) {
    // Wire assignments for elliptic gate:
    // x1 = w_r, y1 = w_o, x2 = w_l_shift, y2 = w_4_shift,
    // x3 = w_r_shift, y3 = w_o_shift
    // q_sign = q_l, q_is_double = q_m, q_elliptic = q_elliptic
    let x_1 = input.w_r();
    let y_1 = input.w_o();
    let x_2 = input.w_l_shift();
    let y_2 = input.w_4_shift();
    let x_3 = input.w_r_shift();
    let y_3 = input.w_o_shift();
    let q_sign = input.q_l();
    let q_is_double = input.q_m();
    let q_elliptic = input.q_elliptic();

    let x2_sub_x1 = x_2 - x_1;
    let x1_mul_3 = x_1 + x_1 + x_1;
    let x3_sub_x1 = x_3 - x_1;
    let x3_plus_two_x1 = x3_sub_x1 + x1_mul_3; // x3 + 2*x1
    let x3_plus_x2_plus_x1 = x3_plus_two_x1 + x2_sub_x1; // x3 + x2 + x1

    // Point addition, x-coordinate:
    // (x3 + x2 + x1)(x2 - x1)^2 - y2^2 - y1^2 + 2*y2*q_sign*y1 = 0
    let y2_sqr = y_2 * y_2;
    let y1_sqr = y_1 * y_1;
    let y2_mul_q_sign = y_2 * q_sign;
    let x_add_identity = x3_plus_x2_plus_x1 * (x2_sub_x1 * x2_sub_x1)
        - (y2_sqr + y1_sqr)
        + (y2_mul_q_sign + y2_mul_q_sign) * y_1;

    let q_elliptic_by_scaling = q_elliptic * *scaling_factor;
    let q_elliptic_q_double_scaling = q_elliptic_by_scaling * q_is_double;
    // a * q_elliptic * (q_is_double - 1) = (q_elliptic_q_double_scaling - q_elliptic_by_scaling)
    let neg_q_elliptic_not_double_scaling = q_elliptic_q_double_scaling - q_elliptic_by_scaling;

    evals[0] = evals[0] - x_add_identity * neg_q_elliptic_not_double_scaling;

    // Point addition, y-coordinate:
    // (y1 + y3)(x2 - x1) + (x3 - x1)(q_sign*y2 - y1) = 0
    let y1_plus_y3 = y_1 + y_3;
    let y_diff = y2_mul_q_sign - y_1;
    let y_add_identity = y1_plus_y3 * x2_sub_x1 + x3_sub_x1 * y_diff;
    evals[1] = evals[1] - y_add_identity * neg_q_elliptic_not_double_scaling;

    // Point doubling, x-coordinate:
    // (x3 + 2*x1)*4*y1^2 - 9*x1*(y1^2 - b) = 0
    let curve_b = get_curve_b::<P>();
    let x_pow_4_mul_3 = (y1_sqr - curve_b) * x1_mul_3;
    let mut y1_sqr_mul_4 = y1_sqr + y1_sqr;
    y1_sqr_mul_4 = y1_sqr_mul_4 + y1_sqr_mul_4;
    let x1_pow_4_mul_9 = x_pow_4_mul_3 + x_pow_4_mul_3 + x_pow_4_mul_3;
    let x_double_identity = x3_plus_two_x1 * y1_sqr_mul_4 - x1_pow_4_mul_9;
    evals[0] = evals[0] + x_double_identity * q_elliptic_q_double_scaling;

    // Point doubling, y-coordinate:
    // (y1 + y3)(2*y1) - (3*x1^2)(x1 - x3) = 0
    // Negated: (3*x1^2)(x3 - x1) + (2*y1)(y1 + y3)
    let x1_sqr_mul_3 = x1_mul_3 * x_1;
    let neg_y_double_identity = x1_sqr_mul_3 * x3_sub_x1 + (y_1 + y_1) * y1_plus_y3;
    evals[1] = evals[1] - neg_y_double_identity * q_elliptic_q_double_scaling;
}

pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
    input.q_elliptic().is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_elliptic_relation_consistency() {
        let run_test = |random_inputs: bool| {
            let mut input = if random_inputs {
                InputElements::<P>::get_random()
            } else {
                InputElements::<P>::get_special()
            };

            // Set q_sign = -1 (the test explicitly sets this)
            input.set_q_l(Fr::zero() - Fr::one());

            let x_1 = input.w_r();
            let y_1 = input.w_o();
            let x_2 = input.w_l_shift();
            let y_2 = input.w_4_shift();
            let x_3 = input.w_r_shift();
            let y_3 = input.w_o_shift();
            let q_sign = input.q_l();
            let q_elliptic = input.q_elliptic();
            let q_is_double = input.q_m();

            let y_diff = q_sign * y_2 - y_1;
            let x_diff = x_2 - x_1;
            let x_diff_sqr = x_diff * x_diff;
            let lambda = y_diff * x_diff.invert();
            let lambda_sqr = lambda * lambda;

            // Point addition x: (x3 - lambda^2 + x1 + x2) * (x2 - x1)^2 = 0
            let x_add_identity = (x_3 - lambda_sqr + x_1 + x_2) * x_diff_sqr;
            // Point addition y: (y3 - lambda*(x1 - x3) + y1) * (x2 - x1)
            let y_add_identity = (y_3 - lambda * (x_1 - x_3) + y_1) * x_diff;

            // Point doubling
            let curve_b = get_curve_b::<P>();
            let y1_sqr = y_1 * y_1;
            let x_pow_4 = (y1_sqr - curve_b) * x_1;
            let lambda_sqr_d = x_pow_4 * Fr::from(9u64) * (y1_sqr * Fr::from(4u64)).invert();
            let lambda_d = (x_1 * x_1 * Fr::from(3u64)) * (y_1 * Fr::from(2u64)).invert();

            let x_double_identity = (x_3 - lambda_sqr_d + x_1 * Fr::from(2u64)) * (y1_sqr * Fr::from(4u64));
            let y_double_identity = (y_3 - lambda_d * (x_1 - x_3) + y_1) * (y_1 * Fr::from(2u64)) * (Fr::zero() - Fr::one());

            let mut expected = [Fr::zero(); NUM_SUBRELATIONS];
            expected[0] = (x_add_identity * (Fr::zero() - q_is_double + Fr::one()) + x_double_identity * q_is_double) * q_elliptic;
            expected[1] = (y_add_identity * (Fr::zero() - q_is_double + Fr::one()) + y_double_identity * q_is_double) * q_elliptic;

            let params = RelationParameters::<Fr>::get_random();
            let mut accum = [Fr::zero(); NUM_SUBRELATIONS];
            let one = Fr::one();
            accumulate(&mut accum, &input, &params, &one);
            assert_eq!(accum, expected);
        };

        run_test(false);
        run_test(true);
    }
}
