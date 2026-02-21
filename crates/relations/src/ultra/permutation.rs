//! Port of `permutation_relation.hpp` — Ultra Permutation Relation.
//!
//! Two subrelations with partial lengths [6, 3].

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::input_elements::InputElements;
use crate::relation_parameters::RelationParameters;

pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 2] = [6, 3];
pub const NUM_SUBRELATIONS: usize = 2;

/// Accumulate the Ultra Permutation relation into the output array.
pub fn accumulate<P: FieldParams>(
    evals: &mut [Field<P>; NUM_SUBRELATIONS],
    input: &InputElements<P>,
    params: &RelationParameters<Field<P>>,
    scaling_factor: &Field<P>,
) {
    let w_1 = input.w_l();
    let w_2 = input.w_r();
    let w_3 = input.w_o();
    let w_4 = input.w_4();

    let sigma_1 = input.sigma_1();
    let sigma_2 = input.sigma_2();
    let sigma_3 = input.sigma_3();
    let sigma_4 = input.sigma_4();

    let id_1 = input.id_1();
    let id_2 = input.id_2();
    let id_3 = input.id_3();
    let id_4 = input.id_4();

    let z_perm = input.z_perm();
    let z_perm_shift = input.z_perm_shift();
    let lagrange_first = input.lagrange_first();
    let lagrange_last = input.lagrange_last();

    let beta = params.beta;
    let gamma = params.gamma;
    let public_input_delta = params.public_input_delta;

    // Contribution (1): grand product construction
    let numerator = (w_1 + id_1 * beta + gamma)
        * (w_2 + id_2 * beta + gamma)
        * (w_3 + id_3 * beta + gamma)
        * (w_4 + id_4 * beta + gamma);

    let denominator = (w_1 + sigma_1 * beta + gamma)
        * (w_2 + sigma_2 * beta + gamma)
        * (w_3 + sigma_3 * beta + gamma)
        * (w_4 + sigma_4 * beta + gamma);

    let contribution_1 = (z_perm + lagrange_first) * numerator
        - (z_perm_shift + lagrange_last * public_input_delta) * denominator;
    evals[0] = evals[0] + contribution_1 * *scaling_factor;

    // Contribution (2): left-shiftable polynomial
    let contribution_2 = z_perm_shift * lagrange_last;
    evals[1] = evals[1] + contribution_2 * *scaling_factor;
}

pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
    (input.z_perm() - input.z_perm_shift()).is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_permutation_relation_consistency() {
        let run_test = |random_inputs: bool| {
            let input = if random_inputs {
                InputElements::<P>::get_random()
            } else {
                InputElements::<P>::get_special()
            };

            let params = RelationParameters::<Fr>::get_random();
            let beta = params.beta;
            let gamma = params.gamma;
            let public_input_delta = params.public_input_delta;

            let w_1 = input.w_l();
            let w_2 = input.w_r();
            let w_3 = input.w_o();
            let w_4 = input.w_4();
            let sigma_1 = input.sigma_1();
            let sigma_2 = input.sigma_2();
            let sigma_3 = input.sigma_3();
            let sigma_4 = input.sigma_4();
            let id_1 = input.id_1();
            let id_2 = input.id_2();
            let id_3 = input.id_3();
            let id_4 = input.id_4();
            let z_perm = input.z_perm();
            let z_perm_shift = input.z_perm_shift();
            let lagrange_first = input.lagrange_first();
            let lagrange_last = input.lagrange_last();

            let mut expected = [Fr::zero(); NUM_SUBRELATIONS];

            expected[0] = (z_perm + lagrange_first)
                * (w_1 + id_1 * beta + gamma)
                * (w_2 + id_2 * beta + gamma)
                * (w_3 + id_3 * beta + gamma)
                * (w_4 + id_4 * beta + gamma)
                - (z_perm_shift + lagrange_last * public_input_delta)
                    * (w_1 + sigma_1 * beta + gamma)
                    * (w_2 + sigma_2 * beta + gamma)
                    * (w_3 + sigma_3 * beta + gamma)
                    * (w_4 + sigma_4 * beta + gamma);

            expected[1] = z_perm_shift * lagrange_last;

            let mut accum = [Fr::zero(); NUM_SUBRELATIONS];
            let one = Fr::one();
            accumulate(&mut accum, &input, &params, &one);
            assert_eq!(accum, expected);
        };

        run_test(false);
        run_test(true);
    }
}
