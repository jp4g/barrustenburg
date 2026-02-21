//! Port of `poseidon2_external_relation.hpp` — Poseidon2 External Round Relation.
//!
//! Four subrelations with partial lengths [7, 7, 7, 7].

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::input_elements::InputElements;
use crate::relation_parameters::RelationParameters;

pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 4] = [7, 7, 7, 7];
pub const NUM_SUBRELATIONS: usize = 4;

/// Accumulate the Poseidon2 External Round relation.
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
    let w_2_shift = input.w_r_shift();
    let w_3_shift = input.w_o_shift();
    let w_4_shift = input.w_4_shift();
    let q_1 = input.q_l();    // c_1 (round constant)
    let q_2 = input.q_r();    // c_2
    let q_3 = input.q_o();    // c_3
    let q_4 = input.q_4();    // c_4
    let q_poseidon2_external = input.q_poseidon2_external();

    // Add round constants
    let s1 = w_1 + q_1;
    let s2 = w_2 + q_2;
    let s3 = w_3 + q_3;
    let s4 = w_4 + q_4;

    // Apply S-box: x^5
    let sbox = |x: Field<P>| -> Field<P> {
        let t2 = x * x;
        let t4 = t2 * t2;
        t4 * x
    };

    let u1 = sbox(s1);
    let u2 = sbox(s2);
    let u3 = sbox(s3);
    let u4 = sbox(s4);

    // Matrix mul v = M_E * u with 14 additions
    let t0 = u1 + u2;           // u_1 + u_2
    let t1 = u3 + u4;           // u_3 + u_4
    let mut t2 = u2 + u2;       // 2u_2
    t2 = t2 + t1;               // 2u_2 + u_3 + u_4
    let mut t3 = u4 + u4;       // 2u_4
    t3 = t3 + t0;               // u_1 + u_2 + 2u_4

    let mut v4 = t1 + t1;
    v4 = v4 + v4;
    v4 = v4 + t3;               // u_1 + u_2 + 4u_3 + 6u_4

    let mut v2 = t0 + t0;
    v2 = v2 + v2;
    v2 = v2 + t2;               // 4u_1 + 6u_2 + u_3 + u_4

    let v1 = t3 + v2;           // 5u_1 + 7u_2 + u_3 + 3u_4
    let v3 = t2 + v4;           // u_1 + 3u_2 + 5u_3 + 7u_4

    let q_pos_by_scaling = q_poseidon2_external * *scaling_factor;

    evals[0] = evals[0] + q_pos_by_scaling * (v1 - w_1_shift);
    evals[1] = evals[1] + q_pos_by_scaling * (v2 - w_2_shift);
    evals[2] = evals[2] + q_pos_by_scaling * (v3 - w_3_shift);
    evals[3] = evals[3] + q_pos_by_scaling * (v4 - w_4_shift);
}

pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
    input.q_poseidon2_external().is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_poseidon2_external_relation_consistency() {
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
            let w_2_shift = input.w_r_shift();
            let w_3_shift = input.w_o_shift();
            let w_4_shift = input.w_4_shift();
            let q_1 = input.q_l();
            let q_2 = input.q_r();
            let q_3 = input.q_o();
            let q_4 = input.q_4();
            let q_poseidon2_external = input.q_poseidon2_external();

            let s1 = w_1 + q_1;
            let s2 = w_2 + q_2;
            let s3 = w_3 + q_3;
            let s4 = w_4 + q_4;

            let sbox = |x: Fr| -> Fr {
                let t2 = x * x;
                let t4 = t2 * t2;
                t4 * x
            };

            let u1 = sbox(s1);
            let u2 = sbox(s2);
            let u3 = sbox(s3);
            let u4 = sbox(s4);

            let t0 = u1 + u2;
            let t1 = u3 + u4;
            let mut t2 = u2 + u2;
            t2 = t2 + t1;
            let mut t3 = u4 + u4;
            t3 = t3 + t0;
            let mut v4 = t1 + t1;
            v4 = v4 + v4;
            v4 = v4 + t3;
            let mut v2 = t0 + t0;
            v2 = v2 + v2;
            v2 = v2 + t2;
            let v1 = t3 + v2;
            let v3 = t2 + v4;

            let mut expected = [Fr::zero(); NUM_SUBRELATIONS];
            expected[0] = q_poseidon2_external * (v1 - w_1_shift);
            expected[1] = q_poseidon2_external * (v2 - w_2_shift);
            expected[2] = q_poseidon2_external * (v3 - w_3_shift);
            expected[3] = q_poseidon2_external * (v4 - w_4_shift);

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
