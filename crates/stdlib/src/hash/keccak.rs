//! In-circuit Keccak hash function.
//!
//! Port of C++ `barretenberg/stdlib/hash/keccak/keccak.{hpp,cpp}`.
//!
//! Creates constraints that evaluate the Keccak-f[1600] permutation using
//! plookup tables for efficient base-11 sparse representation operations.
//! Ultra only due to heavy lookup table use.

use bbrs_circuit_builder::plookup_tables::types::MultiTableId;
use bbrs_ecc::curves::bn254::Bn254FrParams;
use bbrs_ecc::fields::field::Field;
use bbrs_numeric::uint256::U256Ext;
use bbrs_numeric::U256;
use std::sync::LazyLock;

use crate::primitives::field::FieldT;
use crate::primitives::plookup;
use crate::primitives::witness::BuilderRef;

type P = Bn254FrParams;
type Fr = Field<P>;

// ════════════════════════════════════════════════════════════════════════
//  Constants
// ════════════════════════════════════════════════════════════════════════

const BASE: u64 = 11;
const NUM_KECCAK_ROUNDS: usize = 24;
const NUM_KECCAK_LANES: usize = 25;
const KECCAK_LANE_SIZE: usize = 64;
const MAXIMUM_MULTITABLE_BITS: usize = 8;

/// Keccak round constants (IOTA round).
const RC: [u64; NUM_KECCAK_ROUNDS] = [
    0x0000000000000001, 0x0000000000008082, 0x800000000000808a, 0x8000000080008000,
    0x000000000000808b, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
    0x000000000000008a, 0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
    0x000000008000808b, 0x800000000000008b, 0x8000000000008089, 0x8000000000008003,
    0x8000000000008002, 0x8000000000000080, 0x000000000000800a, 0x800000008000000a,
    0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008,
];

/// Rotation offsets, y vertically, x horizontally: r[y * 5 + x].
const ROTATIONS: [usize; NUM_KECCAK_LANES] = [
    0, 1, 62, 28, 27, 36, 44, 6, 55, 20, 3, 10, 43, 25, 39, 41, 45, 15, 21, 8, 18, 2, 61, 56, 14,
];

/// Sparse round constants: RC values converted to base-11 representation.
static SPARSE_RC: LazyLock<[U256; NUM_KECCAK_ROUNDS]> = LazyLock::new(|| {
    let mut out = [U256::ZERO; NUM_KECCAK_ROUNDS];
    for i in 0..NUM_KECCAK_ROUNDS {
        out[i] = convert_to_sparse(U256::from_words([RC[i], 0, 0, 0]));
    }
    out
});

/// Chi offset: sum_{i=0}^{63} 11^i. Added in the CHI round to compute 1 + 2A - B + C.
static CHI_OFFSET: LazyLock<U256> = LazyLock::new(|| {
    let mut result = U256::ZERO;
    let base = U256::from_words([BASE, 0, 0, 0]);
    for _ in 0..KECCAK_LANE_SIZE {
        result = result.wrapping_mul(&base);
        result = result.wrapping_add(&U256::ONE);
    }
    result
});

// ════════════════════════════════════════════════════════════════════════
//  Sparse representation helpers
// ════════════════════════════════════════════════════════════════════════

/// Convert a binary integer into a base-11 integer.
///
/// Input  = sum_{i=0}^{63} b_i * 2^i
/// Output = sum_{i=0}^{63} b_i * 11^i
pub fn convert_to_sparse(input: U256) -> U256 {
    let mut bits = [0u64; 64];
    let mut count = 0usize;
    let mut val = input;
    while val != U256::ZERO {
        bits[count] = val.as_words()[0] & 1;
        count += 1;
        val = val.wrapping_shr_vartime(1);
    }
    let mut output = U256::ZERO;
    let base = U256::from_words([BASE, 0, 0, 0]);
    for i in (0..count).rev() {
        output = output.wrapping_mul(&base);
        output = output.wrapping_add(&U256::from_words([bits[i], 0, 0, 0]));
    }
    output
}

/// Normalize a base-11 integer where each base value can be > 1.
///
/// Input  = sum_{i=0}^{63} b_i * 11^i
/// Output = sum_{i=0}^{63} (b_i & 1) * 11^i
fn normalize_sparse(input: U256) -> U256 {
    let mut slices = [0u64; 64];
    let mut count = 0usize;
    let mut val = input;
    let base = U256::from_words([BASE, 0, 0, 0]);
    let base_nz = base.to_nz().unwrap();
    while val != U256::ZERO {
        let (quotient, remainder) = val.div_rem(&base_nz);
        slices[count] = remainder.as_words()[0] & 1;
        count += 1;
        val = quotient;
    }
    let mut out = U256::ZERO;
    for i in (0..count).rev() {
        out = out.wrapping_mul(&base);
        out = out.wrapping_add(&U256::from_words([slices[i], 0, 0, 0]));
    }
    out
}

// ════════════════════════════════════════════════════════════════════════
//  Internal state
// ════════════════════════════════════════════════════════════════════════

struct KeccakState {
    state: Vec<FieldT<P>>,
    state_msb: Vec<FieldT<P>>,
    twisted_state: Vec<FieldT<P>>,
    context: BuilderRef<P>,
}

// ════════════════════════════════════════════════════════════════════════
//  normalize_and_rotate
// ════════════════════════════════════════════════════════════════════════

/// Normalize a base-11 limb and left-rotate by ROTATIONS[lane_index] bits.
/// Also extracts the most significant bit of the normalized limb.
///
/// Uses the KECCAK_NORMALIZE_AND_ROTATE plookup table.
///
/// Returns (rotated_normalized_limb, msb).
fn normalize_and_rotate(lane_index: usize, limb: &FieldT<P>) -> (FieldT<P>, FieldT<P>) {
    let left_bits = ROTATIONS[lane_index];
    let right_bits = KECCAK_LANE_SIZE - left_bits;

    let num_right_tables = (right_bits + MAXIMUM_MULTITABLE_BITS - 1) / MAXIMUM_MULTITABLE_BITS;
    let num_left_tables = if left_bits == 0 {
        0
    } else {
        (left_bits + MAXIMUM_MULTITABLE_BITS - 1) / MAXIMUM_MULTITABLE_BITS
    };
    let total_tables = num_right_tables + num_left_tables;

    let table_id = MultiTableId(MultiTableId::KECCAK_NORMALIZE_AND_ROTATE.0 + lane_index);

    let accumulators = plookup::get_lookup_accumulators(
        table_id,
        limb,
        &FieldT::default(),
        false,
    );

    // MSB is in the last C3 entry
    let msb = accumulators.columns[2][total_tables - 1].clone();

    // Right output is the full accumulator in C2[0]
    let right_output = accumulators.columns[1][0].clone();

    if num_left_tables == 0 {
        (right_output, msb)
    } else {
        // Left output starts at C2[num_right_tables]
        let left_output = accumulators.columns[1][num_right_tables].clone();

        // Stitch: rotated = left + right * 11^rotation
        let shift = u256_pow_base11(ROTATIONS[lane_index]);
        let shift_fr = Fr::from_limbs(shift.limbs());
        let result = left_output + right_output * FieldT::from_field(shift_fr);
        (result, msb)
    }
}

/// Compute 11^exp as U256.
fn u256_pow_base11(exp: usize) -> U256 {
    let base = U256::from_words([BASE, 0, 0, 0]);
    let mut result = U256::ONE;
    for _ in 0..exp {
        result = result.wrapping_mul(&base);
    }
    result
}

// ════════════════════════════════════════════════════════════════════════
//  Round functions
// ════════════════════════════════════════════════════════════════════════

/// Compute twisted representation of hash lanes.
///
/// If the bit slices for a regular variable are [b63, ..., b0],
/// the twisted representation is [b63, ..., b0, b63] (65-bit).
///
/// twisted_limb = state[i] * 11 + state_msb[i]
fn compute_twisted_state(internal: &mut KeccakState) {
    let base_ft = FieldT::from_field(Fr::from(BASE));
    for i in 0..NUM_KECCAK_LANES {
        internal.twisted_state[i] =
            (internal.state[i].clone() * base_ft.clone() + internal.state_msb[i].clone())
                .normalize();
    }
}

/// THETA round.
///
/// Evaluates XOR operations via additions in base-11 representation.
/// Uses twisted representation for cheap left-rotation by 1 bit.
fn theta(internal: &mut KeccakState) {
    let base_ft = FieldT::from_field(Fr::from(BASE));

    // Compute column parity C[i] = XOR of all 5 lanes in column i
    let mut c = Vec::with_capacity(5);
    for i in 0..5 {
        c.push(FieldT::accumulate(&[
            internal.twisted_state[i].clone(),
            internal.twisted_state[5 + i].clone(),
            internal.twisted_state[10 + i].clone(),
            internal.twisted_state[15 + i].clone(),
            internal.twisted_state[20 + i].clone(),
        ]));
    }

    // Compute D[i] using twisted representation for rotation by 1
    let mut d = Vec::with_capacity(5);
    for i in 0..5 {
        let non_shifted = c[(i + 4) % 5].clone();
        let shifted = c[(i + 1) % 5].clone() * base_ft.clone();
        d.push(non_shifted + shifted);
    }

    // D contains 66 base-11 slices. Remove top and bottom slices, normalize.
    let base_u256 = U256::from_words([BASE, 0, 0, 0]);
    let base_nz = base_u256.to_nz().unwrap();
    let divisor = u256_pow_base11(KECCAK_LANE_SIZE);
    let divisor_nz = divisor.to_nz().unwrap();
    let multiplicand = u256_pow_base11(KECCAK_LANE_SIZE + 1);
    let multiplicand_fr = Fr::from_limbs(multiplicand.limbs());
    let base_fr = Fr::from(BASE);

    for i in 0..5 {
        let d_native = d[i].get_value().from_montgomery_form();
        let d_val = U256::from_words(d_native.data);

        let (d_quotient, lo_native) = d_val.div_rem(&base_nz);
        let (hi_native, mid_native) = d_quotient.div_rem(&divisor_nz);

        let hi = FieldT::from_witness(
            internal.context.clone(),
            Fr::from_limbs(hi_native.limbs()),
        );
        let mid = FieldT::from_witness(
            internal.context.clone(),
            Fr::from_limbs(mid_native.limbs()),
        );
        let lo = FieldT::from_witness(
            internal.context.clone(),
            Fr::from_limbs(lo_native.limbs()),
        );

        // D[i] == hi * 11^65 + mid * 11 + lo
        let reconstructed = (hi.clone() * FieldT::from_field(multiplicand_fr))
            .add_two(&(mid.clone() * FieldT::from_field(base_fr)), &lo);
        d[i].assert_equal(&reconstructed, "keccak theta: D decomposition mismatch");

        // Range-constrain hi and lo to [0, 11]
        {
            let mut ctx = internal.context.borrow_mut();
            ctx.create_new_range_constraint(
                hi.get_witness_index(),
                BASE,
                "keccak theta: hi out of range",
            );
            ctx.create_new_range_constraint(
                lo.get_witness_index(),
                BASE,
                "keccak theta: lo out of range",
            );
        }

        // Normalize mid via KECCAK_THETA_OUTPUT table
        d[i] = plookup::read_from_1_to_2_table(MultiTableId::KECCAK_THETA_OUTPUT, &mid);
    }

    // XOR D into state
    for i in 0..5 {
        for j in 0..5 {
            internal.state[j * 5 + i] =
                internal.state[j * 5 + i].clone() + d[i].clone();
        }
    }
}

/// RHO round.
///
/// Normalizes and left-rotates each lane by its rotation offset.
fn rho(internal: &mut KeccakState) {
    for i in 0..NUM_KECCAK_LANES {
        let (rotated, msb) = normalize_and_rotate(i, &internal.state[i]);
        internal.state[i] = rotated;
        internal.state_msb[i] = msb;
    }
}

/// PI round.
///
/// Permutes the keccak lanes. Zero constraints (just reordering).
fn pi(internal: &mut KeccakState) {
    let mut b = Vec::with_capacity(NUM_KECCAK_LANES);
    for i in 0..NUM_KECCAK_LANES {
        b.push(internal.state[i].clone());
    }

    for y in 0..5 {
        for x in 0..5 {
            let u = (0 * x + 1 * y) % 5;
            let v = (2 * x + 3 * y) % 5;
            internal.state[v * 5 + u] = b[5 * y + x].clone();
        }
    }
}

/// CHI round.
///
/// Applies: A XOR (~B AND C) via linear operation 1 + 2A - B + C in base-11.
fn chi(internal: &mut KeccakState) {
    let chi_offset_fr = Fr::from_limbs(CHI_OFFSET.limbs());

    for y in 0..5 {
        let mut lane_outputs = Vec::with_capacity(5);
        for x in 0..5 {
            let a = internal.state[y * 5 + x].clone();
            let b = internal.state[y * 5 + ((x + 1) % 5)].clone();
            let c = internal.state[y * 5 + ((x + 2) % 5)].clone();

            // 1 + 2A - B + C = (A + A + CHI_OFFSET) + (-B) + C
            let two_a_plus_offset = a.clone() + a + FieldT::from_field(chi_offset_fr);
            lane_outputs.push(two_a_plus_offset.add_two(&(-b), &c));
        }

        for x in 0..5 {
            let accumulators = plookup::get_lookup_accumulators(
                MultiTableId::KECCAK_CHI_OUTPUT,
                &lane_outputs[x],
                &FieldT::default(),
                false,
            );
            internal.state[y * 5 + x] = accumulators.columns[1][0].clone();
            let last_c3 = accumulators.columns[2].len() - 1;
            internal.state_msb[y * 5 + x] = accumulators.columns[2][last_c3].clone();
        }
    }
}

/// IOTA round.
///
/// XORs first hash limb with a precomputed round constant, then normalizes.
fn iota(internal: &mut KeccakState, round: usize) {
    let sparse_rc_fr = Fr::from_limbs(SPARSE_RC[round].limbs());
    let xor_result = internal.state[0].clone() + FieldT::from_field(sparse_rc_fr);

    let (normalized, msb) = normalize_and_rotate(0, &xor_result);
    internal.state[0] = normalized;
    internal.state_msb[0] = msb;

    // Compute twisted state for next round (skip for last round)
    if round != NUM_KECCAK_ROUNDS - 1 {
        compute_twisted_state(internal);
    }
}

/// Run the full Keccak-f[1600] permutation (24 rounds).
fn keccakf1600(internal: &mut KeccakState) {
    for i in 0..NUM_KECCAK_ROUNDS {
        theta(internal);
        rho(internal);
        pi(internal);
        chi(internal);
        iota(internal, i);
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Extended ↔ Normal conversion
// ════════════════════════════════════════════════════════════════════════

/// Convert extended base-11 representation back to normal 64-bit lanes.
fn extended_to_normal(internal: &KeccakState) -> Vec<FieldT<P>> {
    let mut conversion = Vec::with_capacity(NUM_KECCAK_LANES);
    for i in 0..NUM_KECCAK_LANES {
        let output_limb = plookup::read_from_1_to_2_table(
            MultiTableId::KECCAK_FORMAT_OUTPUT,
            &internal.state[i],
        );
        conversion.push(output_limb);
    }
    conversion
}

// ════════════════════════════════════════════════════════════════════════
//  Public API
// ════════════════════════════════════════════════════════════════════════

/// Compute the Keccak-f[1600] permutation on a state of 25 64-bit lanes.
///
/// Converts the state into extended base-11 representation, runs keccakf1600,
/// then converts back to normal 64-bit lanes.
///
/// Port of C++ `keccak::permutation_opcode`.
pub fn permutation_opcode(state: &[FieldT<P>; NUM_KECCAK_LANES], ctx: BuilderRef<P>) -> Vec<FieldT<P>> {
    let mut internal = KeccakState {
        state: Vec::with_capacity(NUM_KECCAK_LANES),
        state_msb: Vec::with_capacity(NUM_KECCAK_LANES),
        twisted_state: vec![FieldT::default(); NUM_KECCAK_LANES],
        context: ctx,
    };

    // Convert 64-bit lanes into extended base-11 representation
    for i in 0..NUM_KECCAK_LANES {
        let accumulators = plookup::get_lookup_accumulators(
            MultiTableId::KECCAK_FORMAT_INPUT,
            &state[i],
            &FieldT::default(),
            false,
        );
        internal.state.push(accumulators.columns[1][0].clone());
        let last_c3 = accumulators.columns[2].len() - 1;
        internal.state_msb.push(accumulators.columns[2][last_c3].clone());
    }

    compute_twisted_state(&mut internal);
    keccakf1600(&mut internal);
    extended_to_normal(&internal)
}

// ════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use std::cell::RefCell;
    use std::rc::Rc;

    fn make_builder() -> BuilderRef<P> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    /// Native Keccak-f[1600] permutation for reference testing.
    fn native_keccakf1600(state: &mut [u64; 25]) {
        const KECCAK_RC: [u64; 24] = [
            0x0000000000000001, 0x0000000000008082, 0x800000000000808a,
            0x8000000080008000, 0x000000000000808b, 0x0000000080000001,
            0x8000000080008081, 0x8000000000008009, 0x000000000000008a,
            0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
            0x000000008000808b, 0x800000000000008b, 0x8000000000008089,
            0x8000000000008003, 0x8000000000008002, 0x8000000000000080,
            0x000000000000800a, 0x800000008000000a, 0x8000000080008081,
            0x8000000000008080, 0x0000000080000001, 0x8000000080008008,
        ];
        const KECCAK_ROT: [u32; 25] = [
            0, 1, 62, 28, 27, 36, 44, 6, 55, 20, 3, 10, 43, 25, 39,
            41, 45, 15, 21, 8, 18, 2, 61, 56, 14,
        ];

        for round in 0..24 {
            // Theta
            let mut c = [0u64; 5];
            for x in 0..5 {
                c[x] = state[x] ^ state[x + 5] ^ state[x + 10] ^ state[x + 15] ^ state[x + 20];
            }
            let mut d = [0u64; 5];
            for x in 0..5 {
                d[x] = c[(x + 4) % 5] ^ c[(x + 1) % 5].rotate_left(1);
            }
            for i in 0..25 {
                state[i] ^= d[i % 5];
            }

            // Rho
            let mut rotated = [0u64; 25];
            for i in 0..25 {
                rotated[i] = state[i].rotate_left(KECCAK_ROT[i]);
            }

            // Pi
            let mut permuted = [0u64; 25];
            for y in 0..5 {
                for x in 0..5 {
                    let u = y % 5;
                    let v = (2 * x + 3 * y) % 5;
                    permuted[v * 5 + u] = rotated[5 * y + x];
                }
            }

            // Chi
            let mut chi_out = [0u64; 25];
            for y in 0..5 {
                for x in 0..5 {
                    chi_out[y * 5 + x] = permuted[y * 5 + x]
                        ^ (!permuted[y * 5 + (x + 1) % 5] & permuted[y * 5 + (x + 2) % 5]);
                }
            }

            // Iota
            chi_out[0] ^= KECCAK_RC[round];
            *state = chi_out;
        }
    }

    // ── Table tests ────────────────────────────────────────────────────
    // Note: Individual table tests verify value correctness. Full circuit
    // correctness is validated by the permutation_opcode test which exercises
    // all keccak lookup tables together.

    /// Test KECCAK_FORMAT_INPUT table.
    ///
    /// Port of C++ `stdlib_keccak::keccak_format_input_table`.
    #[test]
    fn test_keccak_format_input_table() {
        let builder = make_builder();

        let test_values: [u64; 5] = [
            0x123456789ABCDEF0,
            0x0000000000000001,
            0xFFFFFFFFFFFFFFFF,
            0x8000000000000000,
            0xDEADBEEFCAFEBABE,
        ];

        for &limb_native in &test_values {
            let limb = FieldT::from_witness(builder.clone(), Fr::from(limb_native));
            let accumulators = plookup::get_lookup_accumulators(
                MultiTableId::KECCAK_FORMAT_INPUT,
                &limb,
                &FieldT::default(),
                false,
            );

            let sparse_limb = &accumulators.columns[1][0];
            let last_c3 = accumulators.columns[2].len() - 1;
            let msb = &accumulators.columns[2][last_c3];

            let expected_sparse = convert_to_sparse(U256::from_words([limb_native, 0, 0, 0]));
            let expected_msb = (limb_native >> 63) & 1;

            let sparse_val = sparse_limb.get_value().from_montgomery_form();
            assert_eq!(
                U256::from_words(sparse_val.data),
                expected_sparse,
                "sparse mismatch for {limb_native:#x}"
            );
            assert_eq!(
                msb.get_value(),
                Fr::from(expected_msb),
                "msb mismatch for {limb_native:#x}"
            );
        }
    }

    /// Test KECCAK_FORMAT_OUTPUT table.
    ///
    /// Port of C++ `stdlib_keccak::keccak_format_output_table`.
    #[test]
    fn test_keccak_format_output_table() {
        let builder = make_builder();

        let test_values: [u64; 5] = [
            0x123456789ABCDEF0,
            0x0000000000000001,
            0xFFFFFFFFFFFFFFFF,
            0x8000000000000000,
            0xDEADBEEFCAFEBABE,
        ];

        for &limb_native in &test_values {
            let extended_native = convert_to_sparse(U256::from_words([limb_native, 0, 0, 0]));
            let extended_fr = Fr::from_limbs(extended_native.limbs());
            let limb = FieldT::from_witness(builder.clone(), extended_fr);

            let normalized_limb = plookup::read_from_1_to_2_table(
                MultiTableId::KECCAK_FORMAT_OUTPUT,
                &limb,
            );

            assert_eq!(
                normalized_limb.get_value(),
                Fr::from(limb_native),
                "output mismatch for {limb_native:#x}"
            );
        }
    }

    /// Test KECCAK_THETA_OUTPUT table.
    ///
    /// Port of C++ `stdlib_keccak::keccak_theta_output_table`.
    #[test]
    fn test_keccak_theta_output_table() {
        let builder = make_builder();

        // Test values: base-11 numbers where each "digit" is 0..10
        let base_values: [[u64; 8]; 3] = [
            [0, 1, 2, 3, 4, 5, 6, 7],
            [10, 9, 8, 7, 6, 5, 4, 3],
            [1, 0, 1, 0, 1, 0, 1, 0],
        ];

        for test in &base_values {
            let mut extended_native = 0u64;
            let mut expected_normalized = 0u64;
            for j in 0..8 {
                extended_native = extended_native * 11 + test[j];
                expected_normalized = expected_normalized * 11 + (test[j] & 1);
            }

            let limb = FieldT::from_witness(builder.clone(), Fr::from(extended_native));
            let normalized = plookup::read_from_1_to_2_table(
                MultiTableId::KECCAK_THETA_OUTPUT,
                &limb,
            );

            assert_eq!(
                normalized.get_value(),
                Fr::from(expected_normalized),
                "theta normalization mismatch"
            );
        }
    }

    /// Test normalize_and_rotate (RHO round).
    ///
    /// Port of C++ `stdlib_keccak::keccak_rho_output_table`.
    #[test]
    fn test_keccak_rho_output_table() {
        let builder = make_builder();

        // Test a subset of lanes to keep test time reasonable
        for lane in [0, 1, 5, 12, 24] {
            // Build a quasi-bit extended value (each digit in [0, 1, 2])
            let mut extended_native = U256::ZERO;
            let mut binary_native = 0u64;
            let base = U256::from_words([11, 0, 0, 0]);

            // Use deterministic values based on lane index
            for j in 0..64 {
                extended_native = extended_native.wrapping_mul(&base);
                binary_native = binary_native << 1;
                let base_value = ((lane * 7 + j * 13 + 3) % 3) as u64;
                extended_native = extended_native.wrapping_add(
                    &U256::from_words([base_value, 0, 0, 0]),
                );
                binary_native += base_value & 1;
            }

            let left_bits = ROTATIONS[lane];
            let right_bits = 64 - left_bits;

            let binary_rotated = if left_bits == 0 {
                // No rotation
                binary_native
            } else {
                let left = binary_native >> right_bits;
                let right = binary_native & ((1u64 << right_bits) - 1);
                left | (right << left_bits)
            };

            let expected_limb = convert_to_sparse(U256::from_words([binary_rotated, 0, 0, 0]));
            let expected_msb = binary_native >> 63;

            let extended_fr = Fr::from_limbs(extended_native.limbs());
            let limb = FieldT::from_witness(builder.clone(), extended_fr);
            let (result_limb, result_msb) = normalize_and_rotate(lane, &limb);

            let result_val = result_limb.get_value().from_montgomery_form();
            assert_eq!(
                U256::from_words(result_val.data),
                expected_limb,
                "rho lane {lane}: limb mismatch"
            );
            assert_eq!(
                result_msb.get_value(),
                Fr::from(expected_msb),
                "rho lane {lane}: msb mismatch"
            );
        }
    }

    /// Test CHI output table.
    ///
    /// Port of C++ `stdlib_keccak::keccak_chi_output_table`.
    #[test]
    fn test_keccak_chi_output_table() {
        let builder = make_builder();
        let chi_normalization_table: [u64; 5] = [0, 0, 1, 1, 0];

        let test_inputs: [[u64; 8]; 3] = [
            [0, 1, 2, 3, 4, 0, 1, 2],
            [4, 3, 2, 1, 0, 4, 3, 2],
            [1, 0, 3, 2, 4, 1, 0, 3],
        ];

        for test in &test_inputs {
            let mut normalized_native = 0u64;
            let mut extended_native = 0u64;
            let mut binary_native = 0u64;

            for j in 0..8 {
                extended_native = extended_native * 11 + test[j];
                normalized_native = normalized_native * 11 + chi_normalization_table[test[j] as usize];
                binary_native = (binary_native << 1) + chi_normalization_table[test[j] as usize];
            }

            let limb = FieldT::from_witness(builder.clone(), Fr::from(extended_native));
            let accumulators = plookup::get_lookup_accumulators(
                MultiTableId::KECCAK_CHI_OUTPUT,
                &limb,
                &FieldT::default(),
                false,
            );

            let normalized = &accumulators.columns[1][0];
            let last_c3 = accumulators.columns[2].len() - 1;
            let msb = &accumulators.columns[2][last_c3];

            assert_eq!(
                normalized.get_value(),
                Fr::from(normalized_native),
                "chi normalization mismatch"
            );
            assert_eq!(
                msb.get_value(),
                Fr::from(binary_native >> 63),
                "chi msb mismatch"
            );
        }
    }

    /// Test convert_to_sparse and normalize_sparse helpers.
    #[test]
    fn test_sparse_conversion() {
        // 0 -> 0
        assert_eq!(convert_to_sparse(U256::ZERO), U256::ZERO);

        // 1 -> 1
        assert_eq!(
            convert_to_sparse(U256::ONE),
            U256::ONE
        );

        // 3 (0b11) -> 1 + 1*11 = 12
        let three = U256::from_words([3, 0, 0, 0]);
        let expected = U256::from_words([12, 0, 0, 0]);
        assert_eq!(convert_to_sparse(three), expected);

        // Round-trip: normalize_sparse(convert_to_sparse(x)) == convert_to_sparse(x)
        let test = U256::from_words([0xDEADBEEF, 0, 0, 0]);
        let sparse = convert_to_sparse(test);
        assert_eq!(normalize_sparse(sparse), sparse);
    }

    /// Test full permutation with circuit correctness verification.
    ///
    /// Port of C++ `stdlib_keccak::permutation_opcode`.
    #[test]
    fn test_permutation_opcode() {
        let builder = make_builder();

        // Use deterministic test state
        let mut native_state = [0u64; 25];
        for i in 0..25 {
            native_state[i] = (i as u64 + 1).wrapping_mul(0x0123456789ABCDEFu64);
        }

        // Run native permutation
        let mut expected_state = native_state;
        native_keccakf1600(&mut expected_state);

        // Create circuit state
        let circuit_state: [FieldT<P>; NUM_KECCAK_LANES] =
            std::array::from_fn(|i| {
                FieldT::from_witness(builder.clone(), Fr::from(native_state[i]))
            });

        // Run circuit permutation
        let circuit_output = permutation_opcode(&circuit_state, builder.clone());

        // Compare outputs
        for i in 0..25 {
            let circuit_value = circuit_output[i].get_value();
            assert_eq!(
                circuit_value,
                Fr::from(expected_state[i]),
                "permutation output mismatch at lane {i}"
            );
        }

        // Verify circuit correctness
        UltraCircuitChecker::check(&mut builder.borrow_mut())
            .expect("circuit check failed");
    }
}
