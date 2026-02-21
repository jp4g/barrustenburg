//! Port of `delta_range_constraint_relation.hpp` — Delta Range Constraint Relation.
//!
//! Four subrelations with partial lengths [6, 6, 6, 6].
//! Enforces: D(D-1)(D-2)(D-3) = 0 for D = wire differences.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::input_elements::InputElements;
use crate::relation_parameters::RelationParameters;

pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 4] = [6, 6, 6, 6];
pub const NUM_SUBRELATIONS: usize = 4;

/// Accumulate the Delta Range Constraint relation.
pub fn accumulate<P: FieldParams>(
    evals: &mut [Field<P>; NUM_SUBRELATIONS],
    input: &InputElements<P>,
    _params: &RelationParameters<Field<P>>,
    scaling_factor: &Field<P>,
) {
    let w_1 = input.w_l();
    let w_2 = input.w_r();
    let w_3 = input.w_o();
    let w_4 = input.w_4();
    let w_1_shift = input.w_l_shift();
    let q_delta_range = input.q_delta_range();

    let q_delta_range_scaled = q_delta_range * *scaling_factor;

    let one = Field::<P>::one();
    let two = Field::from(2u64);
    let three = Field::from(3u64);

    // D(D-1)(D-2)(D-3) = ((D-3)*D) * ((D-3)*D + 2) = (D^2 - 3D)(D^2 - 3D + 2)
    let compute_delta_constraint = |delta: Field<P>| -> Field<P> {
        delta * (delta - one) * (delta - two) * (delta - three)
    };

    let delta_1 = w_2 - w_1;
    let delta_2 = w_3 - w_2;
    let delta_3 = w_4 - w_3;
    let delta_4 = w_1_shift - w_4;

    evals[0] = evals[0] + compute_delta_constraint(delta_1) * q_delta_range_scaled;
    evals[1] = evals[1] + compute_delta_constraint(delta_2) * q_delta_range_scaled;
    evals[2] = evals[2] + compute_delta_constraint(delta_3) * q_delta_range_scaled;
    evals[3] = evals[3] + compute_delta_constraint(delta_4) * q_delta_range_scaled;
}

pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
    input.q_delta_range().is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_delta_range_constraint_relation_consistency() {
        let run_test = |random_inputs: bool| {
            let input = if random_inputs {
                InputElements::<P>::get_random()
            } else {
                InputElements::<P>::get_special()
            };

            let w_1 = input.w_l();
            let w_2 = input.w_r();
            let w_3 = input.w_o();
            let w_4 = input.w_4();
            let w_1_shift = input.w_l_shift();
            let q_delta_range = input.q_delta_range();

            let one = Fr::one();
            let two = Fr::from(2u64);
            let three = Fr::from(3u64);

            let delta_1 = w_2 - w_1;
            let delta_2 = w_3 - w_2;
            let delta_3 = w_4 - w_3;
            let delta_4 = w_1_shift - w_4;

            let dc = |d: Fr| d * (d - one) * (d - two) * (d - three);

            let mut expected = [Fr::zero(); NUM_SUBRELATIONS];
            expected[0] = dc(delta_1) * q_delta_range;
            expected[1] = dc(delta_2) * q_delta_range;
            expected[2] = dc(delta_3) * q_delta_range;
            expected[3] = dc(delta_4) * q_delta_range;

            let params = RelationParameters::<Fr>::get_random();
            let mut accum = [Fr::zero(); NUM_SUBRELATIONS];
            accumulate(&mut accum, &input, &params, &one);
            assert_eq!(accum, expected);
        };

        run_test(false);
        run_test(true);
    }
}
