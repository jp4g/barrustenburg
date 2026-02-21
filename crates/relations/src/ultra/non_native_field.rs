//! Port of `non_native_field_relation.hpp` — Non-Native Field Arithmetic Relation.
//!
//! One subrelation with partial length [6].

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::input_elements::InputElements;
use crate::relation_parameters::RelationParameters;

pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 1] = [6];
pub const NUM_SUBRELATIONS: usize = 1;

/// Accumulate the Non-Native Field relation.
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

    let q_2 = input.q_r();
    let q_3 = input.q_o();
    let q_4 = input.q_4();
    let q_m = input.q_m();
    let q_nnf = input.q_nnf();

    // 2^68 = [0, 16, 0, 0] in limb form (68 = 64 + 4, so second limb = 2^4 = 16)
    let limb_size = Field::<P>::from_limbs([0, 16, 0, 0]);
    // 2^14 = 16384
    let sublimb_shift = Field::<P>::from(16384u64);
    let sublimb_shift_2 = sublimb_shift * sublimb_shift;
    let sublimb_shift_3 = sublimb_shift_2 * sublimb_shift;
    let sublimb_shift_4 = sublimb_shift_3 * sublimb_shift;

    // Gate 1: [(w_1 * w_2_shift) + (w_1_shift * w_2)] * LIMB_SIZE + (w_1_shift * w_2_shift) - (w_3 + w_4)
    let mut nnf_gate_1 = (w_1 * w_2_shift + w_1_shift * w_2) * limb_size;
    nnf_gate_1 = nnf_gate_1 + w_1_shift * w_2_shift;
    nnf_gate_1 = nnf_gate_1 - (w_3 + w_4);

    // Gate 2: [(w_1 * w_4) + (w_2 * w_3) - w_3_shift] * LIMB_SIZE - w_4_shift + (w_1 * w_2_shift) + (w_1_shift * w_2)
    let mut nnf_gate_2 = (w_1 * w_4 + w_2 * w_3 - w_3_shift) * limb_size;
    nnf_gate_2 = nnf_gate_2 - w_4_shift;
    nnf_gate_2 = nnf_gate_2 + w_1 * w_2_shift + w_1_shift * w_2;

    // Gate 3: [(w_1 * w_2_shift) + (w_1_shift * w_2)] * LIMB_SIZE + (w_1_shift * w_2_shift) + w_4 - (w_3_shift + w_4_shift)
    let mut nnf_gate_3 = (w_1 * w_2_shift + w_1_shift * w_2) * limb_size;
    nnf_gate_3 = nnf_gate_3 + w_1_shift * w_2_shift;
    nnf_gate_3 = nnf_gate_3 + w_4;
    nnf_gate_3 = nnf_gate_3 - (w_3_shift + w_4_shift);

    // Limb accumulator 1
    let limb_accumulator_1 = w_1 + w_2 * sublimb_shift + w_3 * sublimb_shift_2
        + w_1_shift * sublimb_shift_3 + w_2_shift * sublimb_shift_4 - w_4;

    // Limb accumulator 2
    let limb_accumulator_2 = w_3 + w_4 * sublimb_shift + w_1_shift * sublimb_shift_2
        + w_2_shift * sublimb_shift_3 + w_3_shift * sublimb_shift_4 - w_4_shift;

    // Multiply by selector products
    let nnf_gate_1 = nnf_gate_1 * (q_2 * q_3);
    let nnf_gate_2 = nnf_gate_2 * (q_2 * q_4);
    let nnf_gate_3 = nnf_gate_3 * (q_2 * q_m);
    let limb_accumulator_1 = limb_accumulator_1 * (q_3 * q_4);
    let limb_accumulator_2 = limb_accumulator_2 * (q_3 * q_m);

    let non_native_field_identity = nnf_gate_1 + nnf_gate_2 + nnf_gate_3;
    let limb_accumulator_identity = limb_accumulator_1 + limb_accumulator_2;

    evals[0] = evals[0] + (non_native_field_identity + limb_accumulator_identity) * q_nnf * *scaling_factor;
}

pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
    input.q_nnf().is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_non_native_field_relation_consistency() {
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
            let q_2 = input.q_r();
            let q_3 = input.q_o();
            let q_4 = input.q_4();
            let q_m = input.q_m();
            let q_nnf = input.q_nnf();

            let limb_size = Fr::from_limbs([0, 16, 0, 0]); // 2^68
            let sublimb_shift = Fr::from(16384u64); // 2^14
            let sublimb_shift_2 = sublimb_shift * sublimb_shift;
            let sublimb_shift_3 = sublimb_shift_2 * sublimb_shift;
            let sublimb_shift_4 = sublimb_shift_3 * sublimb_shift;

            let mut nnf_gate_1 = (w_1 * w_2_shift + w_1_shift * w_2) * limb_size;
            nnf_gate_1 = nnf_gate_1 + w_1_shift * w_2_shift;
            nnf_gate_1 = nnf_gate_1 - (w_3 + w_4);

            let mut nnf_gate_2 = (w_1 * w_4 + w_2 * w_3 - w_3_shift) * limb_size;
            nnf_gate_2 = nnf_gate_2 - w_4_shift;
            nnf_gate_2 = nnf_gate_2 + w_1 * w_2_shift + w_1_shift * w_2;

            let mut nnf_gate_3 = (w_1 * w_2_shift + w_1_shift * w_2) * limb_size;
            nnf_gate_3 = nnf_gate_3 + w_1_shift * w_2_shift;
            nnf_gate_3 = nnf_gate_3 + w_4;
            nnf_gate_3 = nnf_gate_3 - (w_3_shift + w_4_shift);

            let limb_accumulator_1 = w_1 + w_2 * sublimb_shift + w_3 * sublimb_shift_2
                + w_1_shift * sublimb_shift_3 + w_2_shift * sublimb_shift_4 - w_4;

            let limb_accumulator_2 = w_3 + w_4 * sublimb_shift + w_1_shift * sublimb_shift_2
                + w_2_shift * sublimb_shift_3 + w_3_shift * sublimb_shift_4 - w_4_shift;

            let nnf_gate_1 = nnf_gate_1 * (q_2 * q_3);
            let nnf_gate_2 = nnf_gate_2 * (q_2 * q_4);
            let nnf_gate_3 = nnf_gate_3 * (q_2 * q_m);
            let limb_accumulator_1 = limb_accumulator_1 * (q_3 * q_4);
            let limb_accumulator_2 = limb_accumulator_2 * (q_3 * q_m);

            let non_native_field_identity = nnf_gate_1 + nnf_gate_2 + nnf_gate_3;
            let limb_accumulator_identity = limb_accumulator_1 + limb_accumulator_2;

            let mut expected = [Fr::zero(); NUM_SUBRELATIONS];
            expected[0] = (non_native_field_identity + limb_accumulator_identity) * q_nnf;

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
