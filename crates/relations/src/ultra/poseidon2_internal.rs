//! Port of `poseidon2_internal_relation.hpp` — Poseidon2 Internal Round Relation.
//!
//! Four subrelations with partial lengths [7, 7, 7, 7].

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::input_elements::InputElements;
use crate::relation_parameters::RelationParameters;

pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 4] = [7, 7, 7, 7];
pub const NUM_SUBRELATIONS: usize = 4;

/// Get the Poseidon2 internal matrix diagonal values.
///
/// These are the `internal_matrix_diagonal` constants from the C++ code.
fn get_internal_matrix_diagonal<P: FieldParams>() -> [Field<P>; 4] {
    // Hardcoded from poseidon2_params.hpp (these are d_i - 1, i.e. the diagonal minus one)
    // Hex:  10dc6e9c006ea38b 04b1e03b4bd9490c 0d03f98929ca1d7f b56821fd19d3b6e7
    //       0c28145b6a44df3e 0149b3d0a30b3bb5 99df9756d4dd9b84 a86b38cfb45a740b
    //       00544b8338791518 b2c7645a50392798 b21f75bb60e35961 70067d00141cac15
    //       222c01175718386f 2e2e82eb122789e3 52e105a3b8fa8526 13bc534433ee428b
    [
        Field::from_limbs([0xb56821fd19d3b6e7, 0x0d03f98929ca1d7f, 0x04b1e03b4bd9490c, 0x10dc6e9c006ea38b]),
        Field::from_limbs([0xa86b38cfb45a740b, 0x99df9756d4dd9b84, 0x0149b3d0a30b3bb5, 0x0c28145b6a44df3e]),
        Field::from_limbs([0x70067d00141cac15, 0xb21f75bb60e35961, 0xb2c7645a50392798, 0x00544b8338791518]),
        Field::from_limbs([0x13bc534433ee428b, 0x52e105a3b8fa8526, 0x2e2e82eb122789e3, 0x222c01175718386f]),
    ]
}

/// Accumulate the Poseidon2 Internal Round relation.
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
    let q_1 = input.q_l();  // round constant c_0
    let q_poseidon2_internal = input.q_poseidon2_internal();

    let diag = get_internal_matrix_diagonal::<P>();

    // Add round constant and apply S-box to first element only
    let v1 = w_1 + q_1;
    let mut u1 = v1 * v1; // v1^2
    u1 = u1 * u1;          // v1^4
    u1 = u1 * v1;          // v1^5

    // Internal matrix multiplication: M_I * [u1, w_2, w_3, w_4]
    // M_I has diagonal = diag, off-diagonal = 1
    // v_i = diag[i] * u_i + sum(all)
    let sum = u1 + w_2 + w_3 + w_4;

    let t0 = u1 * diag[0] + sum;
    let t1 = w_2 * diag[1] + sum;
    let t2 = w_3 * diag[2] + sum;
    let t3 = w_4 * diag[3] + sum;

    let q_pos_by_scaling = q_poseidon2_internal * *scaling_factor;

    evals[0] = evals[0] + q_pos_by_scaling * (t0 - w_1_shift);
    evals[1] = evals[1] + q_pos_by_scaling * (t1 - w_2_shift);
    evals[2] = evals[2] + q_pos_by_scaling * (t2 - w_3_shift);
    evals[3] = evals[3] + q_pos_by_scaling * (t3 - w_4_shift);
}

pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
    input.q_poseidon2_internal().is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_poseidon2_internal_relation_consistency() {
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
            let q_poseidon2_internal = input.q_poseidon2_internal();

            let diag = get_internal_matrix_diagonal::<P>();

            let v1 = w_1 + q_1;
            let mut u1 = v1 * v1;
            u1 = u1 * u1;
            u1 = u1 * v1;

            let sum = u1 + w_2 + w_3 + w_4;
            let t0 = u1 * diag[0] + sum;
            let t1 = w_2 * diag[1] + sum;
            let t2 = w_3 * diag[2] + sum;
            let t3 = w_4 * diag[3] + sum;

            let mut expected = [Fr::zero(); NUM_SUBRELATIONS];
            expected[0] = q_poseidon2_internal * (t0 - w_1_shift);
            expected[1] = q_poseidon2_internal * (t1 - w_2_shift);
            expected[2] = q_poseidon2_internal * (t2 - w_3_shift);
            expected[3] = q_poseidon2_internal * (t3 - w_4_shift);

            let params = RelationParameters::<Fr>::get_random();
            let mut accum = [Fr::zero(); NUM_SUBRELATIONS];
            let one = Fr::one();
            accumulate(&mut accum, &input, &params, &one);
            assert_eq!(accum, expected);
        };

        run_test(false);
        run_test(true);
    }

    /// Port of C++ `RelationManual::Poseidon2InternalRelationZeros`
    ///
    /// Deterministic vector test: all wires zero, selector enabled → output should be zero.
    #[test]
    fn test_poseidon2_internal_relation_zeros() {
        let mut input = InputElements::<P>::get_special();
        for e in input.data.iter_mut() {
            *e = Fr::zero();
        }
        // Set q_poseidon2_internal = 1
        input.data[13] = Fr::one();

        let params = RelationParameters::<Fr>::get_random();
        let mut accum = [Fr::zero(); NUM_SUBRELATIONS];
        let one = Fr::one();
        accumulate(&mut accum, &input, &params, &one);

        for i in 0..NUM_SUBRELATIONS {
            assert!(accum[i].is_zero(), "Subrelation {i} should be zero when all wires are zero");
        }
    }

    /// Port of C++ `RelationManual::Poseidon2InternalRelationRandom`
    ///
    /// Deterministic vector test with known C++ values:
    ///   w = (1, 2, 3, 4), q_l = 5
    ///   u1 = (w_1 + q_1)^5 = (1+5)^5 = 7776
    ///   sum = u1 + w_2 + w_3 + w_4 = 7776 + 2 + 3 + 4 = 7785
    ///   M_I * [u1, w_2, w_3, w_4]:
    ///     w_l_shift = 0x122d9ce41e83c533318954d77a4ebc40eb729f6543ebd5f2e4ecb175ced3bc74
    ///     w_r_shift = 0x185028b6d489be7c029367a14616776b33bf2eada9bb370950d6719f68b5067f
    ///     w_o_shift = 0x00fce289a96b3f4a18562d0ef0ab76ca165e613222aa0c24501377003c5622a8
    ///     w_4_shift = 0x27e7677799fda1694819803f459b76d2fb1c45fdf0773375c72d61e8efb92893
    #[test]
    fn test_poseidon2_internal_relation_random() {
        let mut input = InputElements::<P>::get_special();
        for e in input.data.iter_mut() {
            *e = Fr::zero();
        }
        input.data[13] = Fr::one();                    // q_poseidon2_internal
        input.data[28] = Fr::from(1u64);               // w_l
        input.data[29] = Fr::from(2u64);               // w_r
        input.data[30] = Fr::from(3u64);               // w_o
        input.data[31] = Fr::from(4u64);               // w_4
        input.data[1]  = Fr::from(5u64);               // q_l

        // Expected M_I * [u1, w_2, w_3, w_4] values from C++ (big-endian hex → little-endian limbs)
        // 0x122d9ce41e83c533_318954d77a4ebc40_eb729f6543ebd5f2_e4ecb175ced3bc74
        input.data[35] = Fr::from_limbs([0xe4ecb175ced3bc74, 0xeb729f6543ebd5f2, 0x318954d77a4ebc40, 0x122d9ce41e83c533]);
        // 0x185028b6d489be7c_029367a14616776b_33bf2eada9bb3709_50d6719f68b5067f
        input.data[36] = Fr::from_limbs([0x50d6719f68b5067f, 0x33bf2eada9bb3709, 0x029367a14616776b, 0x185028b6d489be7c]);
        // 0x00fce289a96b3f4a_18562d0ef0ab76ca_165e613222aa0c24_501377003c5622a8
        input.data[37] = Fr::from_limbs([0x501377003c5622a8, 0x165e613222aa0c24, 0x18562d0ef0ab76ca, 0x00fce289a96b3f4a]);
        // 0x27e7677799fda169_4819803f459b76d2_fb1c45fdf0773375_c72d61e8efb92893
        input.data[38] = Fr::from_limbs([0xc72d61e8efb92893, 0xfb1c45fdf0773375, 0x4819803f459b76d2, 0x27e7677799fda169]);

        let params = RelationParameters::<Fr>::get_random();
        let mut accum = [Fr::zero(); NUM_SUBRELATIONS];
        let one = Fr::one();
        accumulate(&mut accum, &input, &params, &one);

        for i in 0..NUM_SUBRELATIONS {
            assert!(accum[i].is_zero(), "Subrelation {i} should be zero for known Poseidon2 internal test vector");
        }
    }
}
