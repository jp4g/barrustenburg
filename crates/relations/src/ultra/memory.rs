//! Port of `memory_relation.hpp` — RAM/ROM Memory Relation.
//!
//! Six subrelations with partial lengths [6, 6, 6, 6, 6, 6].

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::input_elements::InputElements;
use crate::relation_parameters::RelationParameters;

pub const SUBRELATION_PARTIAL_LENGTHS: [usize; 6] = [6, 6, 6, 6, 6, 6];
pub const NUM_SUBRELATIONS: usize = 6;

/// Accumulate the Memory relation.
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
    let w_1_shift = input.w_l_shift();
    let w_2_shift = input.w_r_shift();
    let w_3_shift = input.w_o_shift();
    let w_4_shift = input.w_4_shift();

    let q_1 = input.q_l();
    let q_2 = input.q_r();
    let q_3 = input.q_o();
    let q_4 = input.q_4();
    let q_m = input.q_m();
    let q_c = input.q_c();
    let q_memory = input.q_memory();

    let eta = params.eta;
    let eta_two = params.eta_two;
    let eta_three = params.eta_three;

    let one = Field::<P>::one();

    // Memory Record Check: w3*eta3 + w2*eta2 + w1*eta + qc - w4
    let mut memory_record_check = w_3 * eta_three;
    memory_record_check = memory_record_check + w_2 * eta_two;
    memory_record_check = memory_record_check + w_1 * eta;
    memory_record_check = memory_record_check + q_c;
    let partial_record_check = memory_record_check;
    memory_record_check = memory_record_check - w_4;

    // ROM Consistency Check
    let index_delta = w_1_shift - w_1;
    let record_delta = w_4_shift - w_4;

    let index_is_monotonically_increasing = index_delta * index_delta - index_delta;

    let adjacent_values_match_if_adjacent_indices_match =
        (Field::zero() - index_delta + one) * record_delta;

    let q_memory_by_scaling = q_memory * *scaling_factor;

    evals[1] = evals[1] + adjacent_values_match_if_adjacent_indices_match * (q_1 * q_2) * q_memory_by_scaling;
    evals[2] = evals[2] + index_is_monotonically_increasing * (q_1 * q_2) * q_memory_by_scaling;

    let rom_consistency_check_identity = memory_record_check * (q_1 * q_2);

    // RAM Consistency Check
    let access_type = w_4 - partial_record_check;
    let access_check = access_type * access_type - access_type;

    let mut next_gate_access_type = w_3_shift * eta_three;
    next_gate_access_type = next_gate_access_type + w_2_shift * eta_two;
    next_gate_access_type = next_gate_access_type + w_1_shift * eta;
    next_gate_access_type = w_4_shift - next_gate_access_type;

    let value_delta = w_3_shift - w_3;
    let adjacent_values_match_if_adjacent_indices_match_and_next_access_is_a_read_operation =
        (Field::zero() - index_delta + one) * value_delta * (Field::zero() - next_gate_access_type + one);

    let next_gate_access_type_is_boolean =
        next_gate_access_type * next_gate_access_type - next_gate_access_type;

    evals[3] = evals[3] + adjacent_values_match_if_adjacent_indices_match_and_next_access_is_a_read_operation * q_3 * q_memory_by_scaling;
    evals[4] = evals[4] + index_is_monotonically_increasing * q_3 * q_memory_by_scaling;
    evals[5] = evals[5] + next_gate_access_type_is_boolean * q_3 * q_memory_by_scaling;
    let ram_consistency_check_identity = access_check * q_3;

    // RAM/ROM access check gate
    let memory_record_check_scaled = memory_record_check * (q_1 * q_m);

    // RAM Timestamp Consistency Check
    let timestamp_delta = w_2_shift - w_2;
    let ram_timestamp_check_identity =
        ((Field::zero() - index_delta + one) * timestamp_delta - w_3) * (q_1 * q_4);

    // Complete memory identity
    let mut memory_identity = rom_consistency_check_identity;
    memory_identity = memory_identity + ram_timestamp_check_identity;
    memory_identity = memory_identity + memory_record_check_scaled;
    memory_identity = memory_identity + ram_consistency_check_identity;

    evals[0] = evals[0] + memory_identity * q_memory_by_scaling;
}

pub fn skip<P: FieldParams>(input: &InputElements<P>) -> bool {
    input.q_memory().is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_memory_relation_consistency() {
        let run_test = |random_inputs: bool| {
            let input = if random_inputs {
                InputElements::<P>::get_random()
            } else {
                InputElements::<P>::get_special()
            };

            let params = RelationParameters::<Fr>::get_random();
            let eta = params.eta;
            let eta_two = params.eta_two;
            let eta_three = params.eta_three;

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
            let q_m = input.q_m();
            let q_c = input.q_c();
            let q_memory = input.q_memory();

            let one = Fr::one();
            let neg_one = Fr::zero() - one;

            // Memory Record Check
            let mut memory_record_check = w_3 * eta_three;
            memory_record_check = memory_record_check + w_2 * eta_two;
            memory_record_check = memory_record_check + w_1 * eta;
            memory_record_check = memory_record_check + q_c;
            let partial_record_check = memory_record_check;
            memory_record_check = memory_record_check - w_4;

            // ROM
            let index_delta = w_1_shift - w_1;
            let record_delta = w_4_shift - w_4;
            let index_is_monotonically_increasing = index_delta * index_delta - index_delta;
            let adjacent_values_match_if_adjacent_indices_match =
                (index_delta * neg_one + one) * record_delta;

            let mut expected = [Fr::zero(); NUM_SUBRELATIONS];
            expected[1] = adjacent_values_match_if_adjacent_indices_match * (q_1 * q_2);
            expected[2] = index_is_monotonically_increasing * (q_1 * q_2);
            let rom_consistency_check_identity = memory_record_check * (q_1 * q_2);

            // RAM
            let access_type = w_4 - partial_record_check;
            let access_check = access_type * access_type - access_type;

            let mut next_gate_access_type = w_3_shift * eta_three;
            next_gate_access_type = next_gate_access_type + w_2_shift * eta_two;
            next_gate_access_type = next_gate_access_type + w_1_shift * eta;
            next_gate_access_type = w_4_shift - next_gate_access_type;

            let value_delta = w_3_shift - w_3;
            let adjacent_ram = (index_delta * neg_one + one) * value_delta
                * (next_gate_access_type * neg_one + one);
            let next_gate_access_type_is_boolean =
                next_gate_access_type * next_gate_access_type - next_gate_access_type;

            expected[3] = adjacent_ram * q_3;
            expected[4] = index_is_monotonically_increasing * q_3;
            expected[5] = next_gate_access_type_is_boolean * q_3;
            let ram_consistency_check_identity = access_check * q_3;

            memory_record_check = memory_record_check * (q_1 * q_m);

            let timestamp_delta = w_2_shift - w_2;
            let ram_timestamp_check_identity =
                ((index_delta * neg_one + one) * timestamp_delta - w_3) * (q_1 * q_4);

            let mut memory_identity = rom_consistency_check_identity;
            memory_identity = memory_identity + ram_timestamp_check_identity;
            memory_identity = memory_identity + memory_record_check;
            memory_identity = memory_identity + ram_consistency_check_identity;

            expected[0] = memory_identity * q_memory;
            expected[1] = expected[1] * q_memory;
            expected[2] = expected[2] * q_memory;
            expected[3] = expected[3] * q_memory;
            expected[4] = expected[4] * q_memory;
            expected[5] = expected[5] * q_memory;

            let mut accum = [Fr::zero(); NUM_SUBRELATIONS];
            let scaling = Fr::one();
            accumulate(&mut accum, &input, &params, &scaling);
            assert_eq!(accum, expected);
        };

        run_test(false);
        run_test(true);
    }
}
