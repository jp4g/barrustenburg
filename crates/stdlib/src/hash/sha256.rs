//! In-circuit SHA-256 compression function.
//!
//! Port of `barretenberg/stdlib/hash/sha256/sha256.{hpp,cpp}`.
//!
//! Provides circuit-level SHA-256 compression using plookup tables for the
//! choose (Ch), majority (Maj), and witness extension (sigma) operations.
//! Only the compression function is exposed; this is sufficient for DSL use.

use bbrs_circuit_builder::plookup_tables::sha256 as sha256_tables;
use bbrs_circuit_builder::plookup_tables::sparse;
use bbrs_circuit_builder::plookup_tables::types::MultiTableId;
use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};
use crate::primitives::field::FieldT;
use crate::primitives::plookup;

type P = Bn254FrParams;

// ════════════════════════════════════════════════════════════════════════
//  Constants
// ════════════════════════════════════════════════════════════════════════

#[allow(dead_code)]
const INIT_CONSTANTS: [u64; 8] = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
];

const ROUND_CONSTANTS: [u64; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

// ════════════════════════════════════════════════════════════════════════
//  Helpers
// ════════════════════════════════════════════════════════════════════════

/// Raise Fr to a u64 power.
#[inline]
fn fr_pow(base: &Fr, exp: u64) -> Fr {
    base.pow(&[exp, 0, 0, 0])
}

/// Convert a u64 value to sparse base-16 form as an Fr.
#[inline]
fn sparse_base16(value: u64) -> Fr {
    sparse::get_sparse_table_with_rotation_values(16, 0, [value, 0])[0]
}

/// Compute the left multipliers for witness extension (sigma0).
///
/// These correspond to the C++ constexpr `left_multipliers` computed from base=16:
///   [0]: 16^25 + 16^14
///   [1]: 16^17 + 1
///   [2]: 16^24 + 16^3 + 16^7
///   [3]: 16^11 + 16^15 + 1
fn compute_left_multipliers() -> [Fr; 4] {
    let base = Fr::from(16u64);
    [
        fr_pow(&base, 25) + fr_pow(&base, 14),
        fr_pow(&base, 17) + Fr::one(),
        fr_pow(&base, 24) + fr_pow(&base, 3) + fr_pow(&base, 7),
        fr_pow(&base, 11) + fr_pow(&base, 15) + Fr::one(),
    ]
}

/// Compute the right multipliers for witness extension (sigma1).
///
/// These correspond to the C++ constexpr `right_multipliers` computed from base=16:
///   [0]: 16^15 + 16^13
///   [1]: 16^18 + 16^16
///   [2]: 16^23 + 1
///   [3]: 16^1 + 16^8
fn compute_right_multipliers() -> [Fr; 4] {
    let base = Fr::from(16u64);
    [
        fr_pow(&base, 15) + fr_pow(&base, 13),
        fr_pow(&base, 18) + fr_pow(&base, 16),
        fr_pow(&base, 23) + Fr::one(),
        fr_pow(&base, 1) + fr_pow(&base, 8),
    ]
}

// ════════════════════════════════════════════════════════════════════════
//  Data structures
// ════════════════════════════════════════════════════════════════════════

/// Sparse representation of a 32-bit value.
///
/// Tracks both the normal (binary) and sparse (base-16) representations.
/// For constants, the sparse form is computed eagerly. For witnesses,
/// the sparse form is computed lazily via plookup tables.
#[derive(Clone)]
struct SparseValue {
    normal: FieldT<P>,
    sparse: FieldT<P>,
}

impl SparseValue {
    fn new(input: &FieldT<P>) -> Self {
        let normal = input.clone();
        let sparse = if normal.is_constant() {
            let val = normal.get_value().from_montgomery_form().data[0];
            let ctx_opt = normal.get_context().clone();
            if let Some(ctx) = ctx_opt {
                FieldT::constant_with_context(ctx, sparse_base16(val))
            } else {
                FieldT::from_field(sparse_base16(val))
            }
        } else {
            FieldT::default()
        };
        Self { normal, sparse }
    }
}

/// Sparse witness limbs for message schedule extension.
///
/// When a witness is decomposed via the SHA256_WITNESS_INPUT plookup table,
/// it produces 4 sparse limbs and 4 rotated limbs used for sigma computations.
#[derive(Clone)]
struct SparseWitnessLimbs {
    normal: FieldT<P>,
    sparse_limbs: [FieldT<P>; 4],
    rotated_limbs: [FieldT<P>; 4],
    has_sparse_limbs: bool,
}

impl SparseWitnessLimbs {
    fn new(normal: FieldT<P>) -> Self {
        Self {
            normal,
            sparse_limbs: [
                FieldT::default(),
                FieldT::default(),
                FieldT::default(),
                FieldT::default(),
            ],
            rotated_limbs: [
                FieldT::default(),
                FieldT::default(),
                FieldT::default(),
                FieldT::default(),
            ],
            has_sparse_limbs: false,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Internal functions
// ════════════════════════════════════════════════════════════════════════

/// Convert a witness to sparse form using the SHA256_WITNESS_INPUT plookup table.
fn convert_witness(w: &FieldT<P>) -> SparseWitnessLimbs {
    let mut result = SparseWitnessLimbs::new(w.clone());

    let lookup = plookup::get_lookup_accumulators(
        MultiTableId::SHA256_WITNESS_INPUT,
        w,
        &FieldT::default(),
        false,
    );

    result.sparse_limbs = [
        lookup.columns[1][0].clone(), // C2[0]
        lookup.columns[1][1].clone(), // C2[1]
        lookup.columns[1][2].clone(), // C2[2]
        lookup.columns[1][3].clone(), // C2[3]
    ];
    result.rotated_limbs = [
        lookup.columns[2][0].clone(), // C3[0]
        lookup.columns[2][1].clone(), // C3[1]
        lookup.columns[2][2].clone(), // C3[2]
        lookup.columns[2][3].clone(), // C3[3]
    ];
    result.has_sparse_limbs = true;

    result
}

/// Map a field element into choose sparse form using SHA256_CH_INPUT.
fn map_into_choose_sparse_form(e: &FieldT<P>) -> SparseValue {
    let sparse = plookup::read_from_1_to_2_table(MultiTableId::SHA256_CH_INPUT, e);
    SparseValue {
        normal: e.clone(),
        sparse,
    }
}

/// Map a field element into majority sparse form using SHA256_MAJ_INPUT.
fn map_into_maj_sparse_form(e: &FieldT<P>) -> SparseValue {
    let sparse = plookup::read_from_1_to_2_table(MultiTableId::SHA256_MAJ_INPUT, e);
    SparseValue {
        normal: e.clone(),
        sparse,
    }
}

/// SHA-256 choose function: Ch(e, f, g) = (e AND f) XOR (NOT e AND g).
///
/// Uses plookup tables to compute the sparse-form choose result.
/// Also decomposes `e` via SHA256_CH_INPUT to obtain rotation components.
fn choose(e: &mut SparseValue, f: &SparseValue, g: &SparseValue) -> FieldT<P> {
    let lookup = plookup::get_lookup_accumulators(
        MultiTableId::SHA256_CH_INPUT,
        &e.normal,
        &FieldT::default(),
        false,
    );
    let rotation_coefficients = sha256_tables::get_choose_rotation_multipliers();

    let rotation_result = lookup.columns[2][0].clone(); // C3[0]
    e.sparse = lookup.columns[1][0].clone(); // C2[0]
    let sparse_limb_3 = lookup.columns[1][2].clone(); // C2[2]

    // Compute XOR result using rotation coefficients.
    // xor_result = (rotation_result * 7) + e.sparse * (rotation_coefficients[0] * 7 + 1)
    //              + sparse_limb_3 * (rotation_coefficients[2] * 7)
    let fr_7 = FieldT::from_field(Fr::from(7u64));
    let rot_times_7 = &rotation_result * &fr_7;
    let coeff_0 = FieldT::from_field(rotation_coefficients[0] * Fr::from(7u64) + Fr::one());
    let coeff_2 = FieldT::from_field(rotation_coefficients[2] * Fr::from(7u64));
    let e_term = &e.sparse * &coeff_0;
    let limb3_term = &sparse_limb_3 * &coeff_2;
    let xor_result = rot_times_7.add_two(&e_term, &limb3_term);

    // choose_result_sparse = xor_result + 2*f.sparse + 3*g.sparse
    let f_doubled = &f.sparse + &f.sparse;
    let g_tripled = &(&g.sparse + &g.sparse) + &g.sparse;
    let choose_result_sparse = xor_result.add_two(&f_doubled, &g_tripled);

    // Normalize through SHA256_CH_OUTPUT table
    plookup::read_from_1_to_2_table(MultiTableId::SHA256_CH_OUTPUT, &choose_result_sparse)
}

/// SHA-256 majority function: Maj(a, b, c) = (a AND b) XOR (a AND c) XOR (b AND c).
///
/// Uses plookup tables to compute the sparse-form majority result.
/// Also decomposes `a` via SHA256_MAJ_INPUT to obtain rotation components.
fn majority(a: &mut SparseValue, b: &SparseValue, c: &SparseValue) -> FieldT<P> {
    let lookup = plookup::get_lookup_accumulators(
        MultiTableId::SHA256_MAJ_INPUT,
        &a.normal,
        &FieldT::default(),
        false,
    );
    let rotation_coefficients = sha256_tables::get_majority_rotation_multipliers();

    let rotation_result = lookup.columns[2][0].clone(); // C3[0]
    a.sparse = lookup.columns[1][0].clone(); // C2[0]
    let sparse_accumulator_2 = lookup.columns[1][1].clone(); // C2[1]

    // xor_result = (rotation_result * 4) + a.sparse * (rotation_coefficients[0] * 4 + 1)
    //              + sparse_accumulator_2 * (rotation_coefficients[1] * 4)
    let fr_4 = FieldT::from_field(Fr::from(4u64));
    let rot_times_4 = &rotation_result * &fr_4;
    let coeff_0 = FieldT::from_field(rotation_coefficients[0] * Fr::from(4u64) + Fr::one());
    let coeff_1 = FieldT::from_field(rotation_coefficients[1] * Fr::from(4u64));
    let a_term = &a.sparse * &coeff_0;
    let acc2_term = &sparse_accumulator_2 * &coeff_1;
    let xor_result = rot_times_4.add_two(&a_term, &acc2_term);

    // majority_result_sparse = xor_result + b.sparse + c.sparse
    let majority_result_sparse = xor_result.add_two(&b.sparse, &c.sparse);

    // Normalize through SHA256_MAJ_OUTPUT table
    plookup::read_from_1_to_2_table(MultiTableId::SHA256_MAJ_OUTPUT, &majority_result_sparse)
}

/// Modular 32-bit addition: (a + b) mod 2^32, with an overflow range constraint.
fn add_normalize(a: &FieldT<P>, b: &FieldT<P>) -> FieldT<P> {
    let ctx_opt = a
        .get_context()
        .clone()
        .or_else(|| b.get_context().clone());

    let sum_val = a.get_value() + b.get_value();
    let sum_u256 = sum_val.from_montgomery_form();
    let normalized_sum = sum_u256.data[0] & 0xFFFF_FFFF;

    if a.is_constant() && b.is_constant() {
        if let Some(ctx) = ctx_opt {
            return FieldT::constant_with_context(ctx, Fr::from(normalized_sum));
        } else {
            return FieldT::from_field(Fr::from(normalized_sum));
        }
    }

    let ctx = ctx_opt.expect("add_normalize: at least one input must have context");
    let overflow_val = (sum_u256.data[0] - normalized_sum) >> 32;
    let overflow = FieldT::from_witness(ctx.clone(), Fr::from(overflow_val));

    // result = a + b + overflow * (-2^32) = a + b - overflow * 2^32
    let neg_pow_32 = FieldT::constant_with_context(ctx, -Fr::from(1u64 << 32));
    let overflow_correction = &overflow * &neg_pow_32;
    let result = a.add_two(b, &overflow_correction);

    // Overflow must fit in 3 bits (max sum < 3 * 2^32)
    overflow.create_range_constraint(3, "sha256 add_normalize overflow");

    result
}

/// Extend the 16-word message schedule to 64 words.
///
/// For each word w[i] (16 <= i < 64):
///   sigma0 = ROTR7(w[i-15]) XOR ROTR18(w[i-15]) XOR SHR3(w[i-15])
///   sigma1 = ROTR17(w[i-2]) XOR ROTR19(w[i-2]) XOR SHR10(w[i-2])
///   w[i] = (sigma0 + sigma1 + w[i-16] + w[i-7]) mod 2^32
fn extend_witness(w_in: &[FieldT<P>; 16]) -> [FieldT<P>; 64] {
    let left_muls = compute_left_multipliers();
    let right_muls = compute_right_multipliers();

    let mut w_sparse: Vec<SparseWitnessLimbs> = Vec::with_capacity(64);
    for i in 0..16 {
        w_sparse.push(SparseWitnessLimbs::new(w_in[i].clone()));
    }
    for _ in 16..64 {
        w_sparse.push(SparseWitnessLimbs::new(FieldT::default()));
    }

    for i in 16..64 {
        // Ensure w[i-15] has sparse limbs
        if !w_sparse[i - 15].has_sparse_limbs {
            w_sparse[i - 15] = convert_witness(&w_sparse[i - 15].normal);
        }
        // Ensure w[i-2] has sparse limbs
        if !w_sparse[i - 2].has_sparse_limbs {
            w_sparse[i - 2] = convert_witness(&w_sparse[i - 2].normal);
        }

        // Compute sigma0 sparse terms using left multipliers
        let w_left = &w_sparse[i - 15];
        let left = [
            &w_left.sparse_limbs[0] * &FieldT::from_field(left_muls[0]),
            &w_left.sparse_limbs[1] * &FieldT::from_field(left_muls[1]),
            &w_left.sparse_limbs[2] * &FieldT::from_field(left_muls[2]),
            &w_left.sparse_limbs[3] * &FieldT::from_field(left_muls[3]),
        ];

        // Compute sigma1 sparse terms using right multipliers
        let w_right = &w_sparse[i - 2];
        let right = [
            &w_right.sparse_limbs[0] * &FieldT::from_field(right_muls[0]),
            &w_right.sparse_limbs[1] * &FieldT::from_field(right_muls[1]),
            &w_right.sparse_limbs[2] * &FieldT::from_field(right_muls[2]),
            &w_right.sparse_limbs[3] * &FieldT::from_field(right_muls[3]),
        ];

        // Combine sigma0: left[0] + left[1] + left[2] + left[3] + rotated_limbs[1]
        // then multiply by 4
        let left_xor_sparse = {
            let sum_01_2 = left[0].add_two(&left[1], &left[2]);
            let sum_3_rot = sum_01_2.add_two(&left[3], &w_sparse[i - 15].rotated_limbs[1]);
            &sum_3_rot * &FieldT::from_field(Fr::from(4u64))
        };

        // Combine sigma1 + sigma0:
        // right[0] + right[1] + right[2] + right[3] + rotated_limbs[2] + rotated_limbs[3] + left_xor_sparse
        let xor_result_sparse = {
            let sum_01_2 = right[0].add_two(&right[1], &right[2]);
            let sum_3_rot2 = sum_01_2.add_two(&right[3], &w_sparse[i - 2].rotated_limbs[2]);
            sum_3_rot2.add_two(&w_sparse[i - 2].rotated_limbs[3], &left_xor_sparse)
        };

        // Normalize through witness extension output table
        let xor_result = plookup::read_from_1_to_2_table(
            MultiTableId::SHA256_WITNESS_OUTPUT,
            &xor_result_sparse,
        );

        // w_out_raw = xor_result + w[i-16] + w[i-7]
        let w_out_raw = xor_result.add_two(&w_sparse[i - 16].normal, &w_sparse[i - 7].normal);

        // Truncate to 32 bits with range constraint
        let w_out = if w_out_raw.is_constant() {
            let val = w_out_raw.get_value().from_montgomery_form().data[0] & 0xFFFF_FFFF;
            let ctx = w_out_raw
                .get_context()
                .clone()
                .unwrap_or_else(|| w_in[0].get_context().clone().unwrap());
            FieldT::constant_with_context(ctx, Fr::from(val))
        } else {
            let ctx = w_out_raw.get_context().clone().unwrap();
            let val = w_out_raw.get_value().from_montgomery_form().data[0] & 0xFFFF_FFFF;
            let w_out = FieldT::from_witness(ctx, Fr::from(val));

            // Range-check the overflow: (w_out_raw - w_out) / 2^32 must fit in 3 bits
            let inv_pow_two = Fr::from(2u64).pow(&[32, 0, 0, 0]).invert();
            let w_out_raw_scaled = &w_out_raw * &FieldT::from_field(inv_pow_two);
            let w_out_scaled = &w_out * &FieldT::from_field(inv_pow_two);
            let divisor = &w_out_raw_scaled - &w_out_scaled;
            divisor.create_range_constraint(3, "sha256 extend_witness overflow");
            w_out
        };

        w_sparse[i] = SparseWitnessLimbs::new(w_out);
    }

    // Extract the normal values
    let mut w_extended: [FieldT<P>; 64] = std::array::from_fn(|_| FieldT::default());
    for i in 0..64 {
        w_extended[i] = w_sparse[i].normal.clone();
    }
    w_extended
}

// ════════════════════════════════════════════════════════════════════════
//  Public API
// ════════════════════════════════════════════════════════════════════════

/// Apply the SHA-256 compression function to a single 512-bit message block.
///
/// This is the only public entry point for the stdlib SHA-256 implementation.
/// We implement only the compression function (rather than a full hash) because
/// this is all that is required in DSL.
///
/// # Arguments
/// - `h_init`: The 8-word (256-bit) initial hash state. For the first block,
///   this should be the standard SHA-256 IV. For subsequent blocks, this is
///   the output of the previous compression.
/// - `input`: The 16-word (512-bit) message block to compress.
///
/// # Returns
/// The updated 8-word hash state after compression.
pub fn sha256_block(
    h_init: &[FieldT<P>; 8],
    input: &[FieldT<P>; 16],
) -> [FieldT<P>; 8] {
    // Initialize round variables with previous block output.
    // a and e use plain SparseValue (their sparse form will be computed
    // on first use in majority/choose respectively).
    let mut a = SparseValue::new(&h_init[0]);
    let mut b = map_into_maj_sparse_form(&h_init[1]);
    let mut c = map_into_maj_sparse_form(&h_init[2]);
    let mut d = SparseValue::new(&h_init[3]);
    let mut e = SparseValue::new(&h_init[4]);
    let mut f = map_into_choose_sparse_form(&h_init[5]);
    let mut g = map_into_choose_sparse_form(&h_init[6]);
    let mut h = SparseValue::new(&h_init[7]);

    // Extend the 16-word message to 64 words
    let w = extend_witness(input);

    // Apply SHA-256 compression rounds
    for i in 0..64 {
        let ch = choose(&mut e, &f, &g);
        let maj = majority(&mut a, &b, &c);

        // temp1 = ch + h + w[i] + round_constant[i]
        let w_plus_k = &w[i] + &FieldT::from_field(Fr::from(ROUND_CONSTANTS[i]));
        let temp1 = ch.add_two(&h.normal, &w_plus_k);

        h = g;
        g = f;
        f = SparseValue {
            normal: e.normal.clone(),
            sparse: e.sparse.clone(),
        };
        e.normal = add_normalize(&d.normal, &temp1);
        e.sparse = FieldT::default();
        d = SparseValue {
            normal: c.normal.clone(),
            sparse: c.sparse.clone(),
        };
        c = SparseValue {
            normal: b.normal.clone(),
            sparse: b.sparse.clone(),
        };
        b = SparseValue {
            normal: a.normal.clone(),
            sparse: a.sparse.clone(),
        };
        a.normal = add_normalize(&temp1, &maj);
        a.sparse = FieldT::default();
    }

    // Add into previous block output
    let mut output: [FieldT<P>; 8] = std::array::from_fn(|_| FieldT::default());
    output[0] = add_normalize(&a.normal, &h_init[0]);
    output[1] = add_normalize(&b.normal, &h_init[1]);
    output[2] = add_normalize(&c.normal, &h_init[2]);
    output[3] = add_normalize(&d.normal, &h_init[3]);
    output[4] = add_normalize(&e.normal, &h_init[4]);
    output[5] = add_normalize(&f.normal, &h_init[5]);
    output[6] = add_normalize(&g.normal, &h_init[6]);
    output[7] = add_normalize(&h.normal, &h_init[7]);

    // Range-check outputs to 32 bits (prevent overflow attacks on add_normalize)
    for out in &output {
        out.create_range_constraint(32, "sha256_block output range check");
    }

    output
}

// ════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use std::rc::Rc;

    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use crate::primitives::witness::{BuilderRef, IS_CONSTANT};

    fn make_builder() -> BuilderRef<P> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    fn check_circuit(builder: &BuilderRef<P>) -> Result<(), String> {
        UltraCircuitChecker::check(&mut builder.borrow_mut())
    }

    /// Minimal test: does a single SHA256 plookup table read work?
    #[test]
    fn test_sha256_basic_plookup() {
        let builder = make_builder();
        let val = FieldT::from_witness(builder.clone(), Fr::from(0x12345678u64));

        // Try reading from SHA256_CH_INPUT (1-to-2 table)
        let _sparse = plookup::read_from_1_to_2_table(MultiTableId::SHA256_CH_INPUT, &val);

        check_circuit(&builder).expect("basic SHA256 plookup failed");
    }

    /// Test sha256_block against NIST vector one ("abc").
    ///
    /// Padded block: "abc" + 0x80 + zeros + 64-bit length (24 bits).
    /// Single block since message fits in 55 bytes.
    #[test]
    fn test_sha256_block_nist_vector_one() {
        let builder = make_builder();

        let h_init_vals: [u32; 8] = [
            0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
            0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
        ];
        let padded_block: [u32; 16] = [
            0x61626380, // "abc" + padding bit
            0x00000000, 0x00000000, 0x00000000,
            0x00000000, 0x00000000, 0x00000000, 0x00000000,
            0x00000000, 0x00000000, 0x00000000, 0x00000000,
            0x00000000, 0x00000000, 0x00000000,
            0x00000018, // length in bits (24)
        ];
        let expected: [u32; 8] = [
            0xba7816bf, 0x8f01cfea, 0x414140de, 0x5dae2223,
            0xb00361a3, 0x96177a9c, 0xb410ff61, 0xf20015ad,
        ];

        // Verify native implementation first
        let native_output = bbrs_crypto::sha256::sha256_block(&h_init_vals, &padded_block);
        for i in 0..8 {
            assert_eq!(native_output[i], expected[i], "Native mismatch at index {}", i);
        }

        // Create circuit witnesses
        let h_init: [FieldT<P>; 8] = std::array::from_fn(|i| {
            FieldT::from_witness(builder.clone(), Fr::from(h_init_vals[i] as u64))
        });
        let block: [FieldT<P>; 16] = std::array::from_fn(|i| {
            FieldT::from_witness(builder.clone(), Fr::from(padded_block[i] as u64))
        });

        let circuit_output = sha256_block(&h_init, &block);

        // Compare outputs first (before circuit check, to isolate logic vs constraint issues)
        for i in 0..8 {
            let circuit_val = circuit_output[i].get_value().from_montgomery_form().data[0] as u32;
            assert_eq!(circuit_val, expected[i], "Circuit output mismatch at index {}", i);
        }

        // Verify circuit correctness
        check_circuit(&builder).expect("circuit check failed");
    }

    /// Test sha256_block against NIST vector two (56-byte message, two blocks).
    ///
    /// Message: "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"
    /// Tests chained compression across two blocks.
    #[test]
    fn test_sha256_block_nist_vector_two() {
        let builder = make_builder();

        let h_init_vals: [u32; 8] = [
            0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
            0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
        ];
        let block_1: [u32; 16] = [
            0x61626364, 0x62636465, 0x63646566, 0x64656667,
            0x65666768, 0x66676869, 0x6768696a, 0x68696a6b,
            0x696a6b6c, 0x6a6b6c6d, 0x6b6c6d6e, 0x6c6d6e6f,
            0x6d6e6f70, 0x6e6f7071, 0x80000000, 0x00000000,
        ];
        let block_2: [u32; 16] = [
            0x00000000, 0x00000000, 0x00000000, 0x00000000,
            0x00000000, 0x00000000, 0x00000000, 0x00000000,
            0x00000000, 0x00000000, 0x00000000, 0x00000000,
            0x00000000, 0x00000000, 0x00000000, 0x000001c0,
        ];
        let expected: [u32; 8] = [
            0x248d6a61, 0xd20638b8, 0xe5c02693, 0x0c3e6039,
            0xa33ce459, 0x64ff2167, 0xf6ecedd4, 0x19db06c1,
        ];

        // Verify native implementation
        let h_after_block1 = bbrs_crypto::sha256::sha256_block(&h_init_vals, &block_1);
        let native_output = bbrs_crypto::sha256::sha256_block(&h_after_block1, &block_2);
        for i in 0..8 {
            assert_eq!(native_output[i], expected[i], "Native mismatch at index {}", i);
        }

        // Circuit: first block
        let h_init: [FieldT<P>; 8] = std::array::from_fn(|i| {
            FieldT::from_witness(builder.clone(), Fr::from(h_init_vals[i] as u64))
        });
        let b1: [FieldT<P>; 16] = std::array::from_fn(|i| {
            FieldT::from_witness(builder.clone(), Fr::from(block_1[i] as u64))
        });

        let h_mid = sha256_block(&h_init, &b1);

        // Circuit: second block
        let b2: [FieldT<P>; 16] = std::array::from_fn(|i| {
            FieldT::from_witness(builder.clone(), Fr::from(block_2[i] as u64))
        });

        let circuit_output = sha256_block(&h_mid, &b2);

        // Verify circuit correctness
        check_circuit(&builder).expect("circuit check failed");

        // Compare outputs
        for i in 0..8 {
            let circuit_val = circuit_output[i].get_value().from_montgomery_form().data[0] as u32;
            assert_eq!(circuit_val, expected[i], "Circuit mismatch at index {}", i);
        }
    }

    /// Test extend_witness constraints (boomerang attack regression).
    ///
    /// Verifies that modifying any extended witness word causes circuit failure.
    #[test]
    fn test_extend_witness_constraints() {
        let builder = make_builder();

        // Create random-ish deterministic input witnesses
        let mut input: [FieldT<P>; 16] = std::array::from_fn(|_| FieldT::default());
        let seeds: [u32; 16] = [
            0x12345678, 0x9ABCDEF0, 0x13579BDF, 0x2468ACE0,
            0xDEADBEEF, 0xCAFEBABE, 0x01234567, 0x89ABCDEF,
            0xFEDCBA98, 0x76543210, 0xAAAA5555, 0x5555AAAA,
            0xFFFF0000, 0x0000FFFF, 0xA5A5A5A5, 0x5A5A5A5A,
        ];
        for i in 0..16 {
            let mut elt = FieldT::from_witness(builder.clone(), Fr::from(seeds[i] as u64));
            elt.fix_witness();
            input[i] = elt;
        }

        // Extend the witness
        let w_ext = extend_witness(&input);

        // Verify circuit is initially valid
        check_circuit(&builder).expect("circuit should be valid before modification");

        // Try modifying each extended witness and verify circuit fails
        let mut any_modification_passed = false;
        for ext_word in &w_ext {
            let variable_index = ext_word.witness_index;
            if variable_index == IS_CONSTANT {
                continue;
            }

            let backup = builder.borrow().base.get_variable(variable_index);

            // Set to a different value
            let mut modified_val = Fr::from(0x42424242u64);
            if modified_val == backup {
                modified_val = Fr::from(0x13131313u64);
            }

            builder.borrow_mut().base.set_variable_unchecked(variable_index, modified_val);

            if check_circuit(&builder).is_ok() {
                any_modification_passed = true;
            }

            builder.borrow_mut().base.set_variable_unchecked(variable_index, backup);
        }

        assert!(!any_modification_passed, "modifying an extended witness should cause circuit failure");
    }
}
