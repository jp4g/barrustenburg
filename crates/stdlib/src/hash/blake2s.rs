//! Circuit Blake2s hash function.
//!
//! Port of `barretenberg/stdlib/hash/blake2s/blake2s_plookup.{hpp,cpp}`.
//!
//! Implements the BLAKE2s-256 hash function in-circuit using plookup tables
//! for XOR-rotation operations and field arithmetic for modular addition.

use bbrs_ecc::curves::bn254::Bn254FrParams;
use bbrs_ecc::fields::field::Field;

use bbrs_circuit_builder::plookup_tables::types::MultiTableId;

use crate::primitives::byte_array::ByteArrayT;
use crate::primitives::field::FieldT;
use crate::primitives::plookup;
use crate::primitives::witness::{BuilderRef, WitnessT};

type P = Bn254FrParams;
type Fr = Field<P>;

// ════════════════════════════════════════════════════════════════════════
//  Constants
// ════════════════════════════════════════════════════════════════════════

/// Blake2s initialization vector.
const IV: [u32; 8] = [
    0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
    0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
];

/// Blake2s message permutation schedule.
const SIGMA: [[usize; 16]; 10] = [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
    [11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4],
    [7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8],
    [9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13],
    [2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9],
    [12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11],
    [13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10],
    [6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5],
    [10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0],
];

// ════════════════════════════════════════════════════════════════════════
//  Helpers: modular 32-bit addition
// ════════════════════════════════════════════════════════════════════════

/// Add two 32-bit circuit values modulo 2^32.
fn add_mod_32(a: &FieldT<P>, b: &FieldT<P>) -> FieldT<P> {
    let ctx = a
        .get_context()
        .clone()
        .or_else(|| b.get_context().clone());

    if a.is_constant() && b.is_constant() {
        let a_val = a.get_value().from_montgomery_form().data[0] as u32;
        let b_val = b.get_value().from_montgomery_form().data[0] as u32;
        let result = Fr::from(a_val.wrapping_add(b_val) as u64);
        return match ctx {
            Some(c) => FieldT::constant_with_context(c, result),
            None => FieldT::from_field(result),
        };
    }

    let ctx = ctx.expect("at least one must have context");

    let a_val = a.get_value().from_montgomery_form().data[0];
    let b_val = b.get_value().from_montgomery_form().data[0];
    let sum = a_val + b_val;
    let carry = sum >> 32;
    let result = sum & 0xFFFF_FFFF;

    let result_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(result)));
    result_wit.create_range_constraint(32, "blake2s add_mod_32: result");

    let carry_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(carry)));
    carry_wit.create_range_constraint(1, "blake2s add_mod_32: carry");

    // Constrain: a + b = result + carry * 2^32
    let shift = FieldT::from_field(Fr::from(1u64 << 32));
    let reconstructed = result_wit.clone() + carry_wit * shift;
    let sum_field = a.clone() + b.clone();
    sum_field.assert_equal(&reconstructed, "blake2s add_mod_32: constraint");

    result_wit
}

/// Add three 32-bit circuit values modulo 2^32.
fn add_mod_32_three(a: &FieldT<P>, b: &FieldT<P>, c: &FieldT<P>) -> FieldT<P> {
    let ctx = a
        .get_context()
        .clone()
        .or_else(|| b.get_context().clone())
        .or_else(|| c.get_context().clone());

    if a.is_constant() && b.is_constant() && c.is_constant() {
        let a_val = a.get_value().from_montgomery_form().data[0] as u32;
        let b_val = b.get_value().from_montgomery_form().data[0] as u32;
        let c_val = c.get_value().from_montgomery_form().data[0] as u32;
        let result = Fr::from(a_val.wrapping_add(b_val).wrapping_add(c_val) as u64);
        return match ctx {
            Some(c) => FieldT::constant_with_context(c, result),
            None => FieldT::from_field(result),
        };
    }

    let ctx = ctx.expect("at least one must have context");

    let a_val = a.get_value().from_montgomery_form().data[0];
    let b_val = b.get_value().from_montgomery_form().data[0];
    let c_val = c.get_value().from_montgomery_form().data[0];
    let sum = a_val + b_val + c_val;
    let carry = sum >> 32;
    let result = sum & 0xFFFF_FFFF;

    let result_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(result)));
    result_wit.create_range_constraint(32, "blake2s add_mod_32_three: result");

    let carry_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(carry)));
    carry_wit.create_range_constraint(2, "blake2s add_mod_32_three: carry");

    // Constrain: a + b + c = result + carry * 2^32
    let shift = FieldT::from_field(Fr::from(1u64 << 32));
    let reconstructed = result_wit.clone() + carry_wit * shift;
    let sum_field = a.add_two(b, c);
    sum_field.assert_equal(&reconstructed, "blake2s add_mod_32_three: constraint");

    result_wit
}

// ════════════════════════════════════════════════════════════════════════
//  Helpers: XOR with rotation
// ════════════════════════════════════════════════════════════════════════

/// XOR two 32-bit values and rotate right by the specified amount.
///
/// Uses `create_logic_constraint` (UINT32_XOR plookup) for XOR, then manual
/// bit decomposition for the rotation.
fn xor_rotate(a: &FieldT<P>, b: &FieldT<P>, rotation: u32) -> FieldT<P> {
    let xor_val = plookup::read_from_2_to_1_table(MultiTableId::UINT32_XOR, a, b);

    if rotation == 0 {
        return xor_val;
    }

    if xor_val.is_constant() {
        let native = xor_val.get_value().from_montgomery_form().data[0] as u32;
        let rotated = native.rotate_right(rotation);
        let result = Fr::from(rotated as u64);
        return match xor_val.get_context().clone() {
            Some(c) => FieldT::constant_with_context(c, result),
            None => FieldT::from_field(result),
        };
    }

    let ctx = xor_val.get_context().as_ref().unwrap().clone();
    let native_xor = xor_val.get_value().from_montgomery_form().data[0] as u32;

    // Split into lo (rotation bits) and hi (32 - rotation bits)
    let lo_bits = rotation;
    let hi_bits = 32 - rotation;
    let lo_mask = (1u32 << lo_bits) - 1;

    let lo = native_xor & lo_mask;
    let hi = native_xor >> lo_bits;

    let lo_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(lo as u64)));
    lo_wit.create_range_constraint(lo_bits as usize, "blake2s xor_rotate: lo");

    let hi_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(hi as u64)));
    hi_wit.create_range_constraint(hi_bits as usize, "blake2s xor_rotate: hi");

    // Constrain: xor_val == lo + hi * 2^lo_bits
    let shift_lo = FieldT::from_field(Fr::from(1u64 << lo_bits));
    let reconstructed = lo_wit.clone() + hi_wit.clone() * shift_lo;
    xor_val.assert_equal(&reconstructed, "blake2s xor_rotate: decomposition");

    // ROTR(xor_val, rotation) = hi + lo * 2^hi_bits
    let shift_hi = FieldT::from_field(Fr::from(1u64 << hi_bits));
    hi_wit + lo_wit * shift_hi
}

// ════════════════════════════════════════════════════════════════════════
//  G mixing function
// ════════════════════════════════════════════════════════════════════════

/// BLAKE2s G mixing function.
fn blake2s_g(
    v: &mut [FieldT<P>; 16],
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    x: &FieldT<P>,
    y: &FieldT<P>,
) {
    let new_a = add_mod_32_three(&v[a], &v[b], x);
    v[a] = new_a;

    let new_d = xor_rotate(&v[d], &v[a], 16);
    v[d] = new_d;

    let new_c = add_mod_32(&v[c], &v[d]);
    v[c] = new_c;

    let new_b = xor_rotate(&v[b], &v[c], 12);
    v[b] = new_b;

    let new_a = add_mod_32_three(&v[a], &v[b], y);
    v[a] = new_a;

    let new_d = xor_rotate(&v[d], &v[a], 8);
    v[d] = new_d;

    let new_c = add_mod_32(&v[c], &v[d]);
    v[c] = new_c;

    let new_b = xor_rotate(&v[b], &v[c], 7);
    v[b] = new_b;
}

// ════════════════════════════════════════════════════════════════════════
//  Compression function
// ════════════════════════════════════════════════════════════════════════

/// BLAKE2s compression function.
fn blake2s_compress(
    h: &mut [FieldT<P>; 8],
    m: &[FieldT<P>; 16],
    t_lo: u32,
    t_hi: u32,
    f: u32,
) {
    let ctx = h[0]
        .get_context()
        .clone()
        .expect("state must have context");

    // Initialize working variables
    let mut v: [FieldT<P>; 16] = [
        h[0].clone(),
        h[1].clone(),
        h[2].clone(),
        h[3].clone(),
        h[4].clone(),
        h[5].clone(),
        h[6].clone(),
        h[7].clone(),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[0] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[1] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[2] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[3] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from((IV[4] ^ t_lo) as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from((IV[5] ^ t_hi) as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from((IV[6] ^ f) as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[7] as u64)),
    ];

    // 10 rounds of mixing
    for round in 0..10 {
        let s = &SIGMA[round];

        // Column step
        blake2s_g(&mut v, 0, 4, 8, 12, &m[s[0]], &m[s[1]]);
        blake2s_g(&mut v, 1, 5, 9, 13, &m[s[2]], &m[s[3]]);
        blake2s_g(&mut v, 2, 6, 10, 14, &m[s[4]], &m[s[5]]);
        blake2s_g(&mut v, 3, 7, 11, 15, &m[s[6]], &m[s[7]]);

        // Diagonal step
        blake2s_g(&mut v, 0, 5, 10, 15, &m[s[8]], &m[s[9]]);
        blake2s_g(&mut v, 1, 6, 11, 12, &m[s[10]], &m[s[11]]);
        blake2s_g(&mut v, 2, 7, 8, 13, &m[s[12]], &m[s[13]]);
        blake2s_g(&mut v, 3, 4, 9, 14, &m[s[14]], &m[s[15]]);
    }

    // Finalize: h[i] = h[i] ^ v[i] ^ v[i+8]
    for i in 0..8 {
        let tmp = plookup::read_from_2_to_1_table(MultiTableId::UINT32_XOR, &h[i], &v[i]);
        h[i] = plookup::read_from_2_to_1_table(MultiTableId::UINT32_XOR, &tmp, &v[i + 8]);
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Byte <-> word conversion helpers
// ════════════════════════════════════════════════════════════════════════

/// Convert 4 consecutive bytes (at `offset`) to a 32-bit little-endian word.
fn bytes_to_word(bytes: &[FieldT<P>], offset: usize) -> FieldT<P> {
    // word = b[0] + b[1]*256 + b[2]*65536 + b[3]*16777216
    let s8 = FieldT::<P>::from_field(Fr::from(256u64));
    let s16 = FieldT::<P>::from_field(Fr::from(65536u64));
    let s24 = FieldT::<P>::from_field(Fr::from(16777216u64));

    let scaled = [
        bytes[offset].clone(),
        bytes[offset + 1].clone() * s8,
        bytes[offset + 2].clone() * s16,
        bytes[offset + 3].clone() * s24,
    ];
    FieldT::accumulate(&scaled)
}

/// Decompose a 32-bit word into 4 bytes (little-endian).
fn word_to_bytes(ctx: &BuilderRef<P>, word: &FieldT<P>) -> [FieldT<P>; 4] {
    let native = word.get_value().from_montgomery_form().data[0] as u32;

    let b0 = (native & 0xFF) as u64;
    let b1 = ((native >> 8) & 0xFF) as u64;
    let b2 = ((native >> 16) & 0xFF) as u64;
    let b3 = ((native >> 24) & 0xFF) as u64;

    if word.is_constant() {
        return [
            FieldT::from_field(Fr::from(b0)),
            FieldT::from_field(Fr::from(b1)),
            FieldT::from_field(Fr::from(b2)),
            FieldT::from_field(Fr::from(b3)),
        ];
    }

    let b0_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(b0)));
    b0_wit.create_range_constraint(8, "blake2s word_to_bytes: b0");
    let b1_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(b1)));
    b1_wit.create_range_constraint(8, "blake2s word_to_bytes: b1");
    let b2_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(b2)));
    b2_wit.create_range_constraint(8, "blake2s word_to_bytes: b2");
    let b3_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(b3)));
    b3_wit.create_range_constraint(8, "blake2s word_to_bytes: b3");

    // Constrain: word = b0 + b1*256 + b2*65536 + b3*16777216
    let s8 = FieldT::<P>::from_field(Fr::from(256u64));
    let s16 = FieldT::<P>::from_field(Fr::from(65536u64));
    let s24 = FieldT::<P>::from_field(Fr::from(16777216u64));
    let reconstructed = FieldT::accumulate(&[
        b0_wit.clone(),
        b1_wit.clone() * s8,
        b2_wit.clone() * s16,
        b3_wit.clone() * s24,
    ]);
    word.assert_equal(&reconstructed, "blake2s word_to_bytes: reconstruction");

    [b0_wit, b1_wit, b2_wit, b3_wit]
}

// ════════════════════════════════════════════════════════════════════════
//  Public API
// ════════════════════════════════════════════════════════════════════════

/// Compute the BLAKE2s-256 hash of a byte array in-circuit.
///
/// Takes a byte array of arbitrary length and returns a 32-byte hash result.
/// Each output byte is range-constrained to 8 bits.
pub fn blake2s(input: ByteArrayT<P>) -> ByteArrayT<P> {
    let ctx = input.get_context().as_ref().unwrap().clone();
    let input_bytes = input.bytes();
    let input_len = input_bytes.len();

    // Initialize hash state with parameter block:
    // h[0] = IV[0] ^ 0x01010020 (fanout=1, depth=1, digest_length=32)
    let mut h: [FieldT<P>; 8] = [
        FieldT::constant_with_context(ctx.clone(), Fr::from((IV[0] ^ 0x0101_0020) as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[1] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[2] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[3] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[4] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[5] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[6] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[7] as u64)),
    ];

    // Pad input to a multiple of 64 bytes
    let num_blocks = if input_len == 0 { 1 } else { (input_len + 63) / 64 };
    let padded_len = num_blocks * 64;

    let mut padded_bytes: Vec<FieldT<P>> = Vec::with_capacity(padded_len);
    padded_bytes.extend_from_slice(input_bytes);
    for _ in input_len..padded_len {
        padded_bytes.push(FieldT::constant_with_context(ctx.clone(), Fr::zero()));
    }

    // Process each block
    let mut bytes_compressed: u32 = 0;

    for block_idx in 0..num_blocks {
        let block_start = block_idx * 64;
        let is_last = block_idx == num_blocks - 1;

        if is_last {
            bytes_compressed = input_len as u32;
        } else {
            bytes_compressed += 64;
        }

        // Convert 64 bytes to 16 32-bit words (little-endian)
        let m: [FieldT<P>; 16] = std::array::from_fn(|i| {
            bytes_to_word(&padded_bytes, block_start + i * 4)
        });

        let f = if is_last { 0xFFFF_FFFFu32 } else { 0u32 };
        blake2s_compress(&mut h, &m, bytes_compressed, 0, f);
    }

    // Convert state to 32 output bytes (little-endian word order)
    let mut output_bytes: Vec<FieldT<P>> = Vec::with_capacity(32);
    for word in &h {
        let wb = word_to_bytes(&ctx, word);
        output_bytes.extend(wb);
    }

    ByteArrayT::from_values(Some(ctx), output_bytes)
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

    fn check_circuit(builder: &BuilderRef<P>) -> Result<(), String> {
        UltraCircuitChecker::check(&mut builder.borrow_mut())
    }

    /// Hash the given bytes in-circuit and verify against the expected output.
    fn run_blake2s_test(input: &[u8], expected: &[u8; 32]) {
        let builder = make_builder();
        let input_arr = ByteArrayT::from_bytes(builder.clone(), input);
        let output = blake2s(input_arr);

        assert_eq!(output.size(), 32);
        assert_eq!(output.get_value(), expected.to_vec(), "hash mismatch for input len={}", input.len());
        check_circuit(&builder).expect("circuit check failed");
    }

    #[test]
    fn test_blake2s_empty() {
        run_blake2s_test(
            b"",
            &[
                0x69, 0x21, 0x7A, 0x30, 0x79, 0x90, 0x80, 0x94, 0xE1, 0x11, 0x21, 0xD0, 0x42,
                0x35, 0x4A, 0x7C, 0x1F, 0x55, 0xB6, 0x48, 0x2C, 0xA1, 0xA5, 0x1E, 0x1B, 0x25,
                0x0D, 0xFD, 0x1E, 0xD0, 0xEE, 0xF9,
            ],
        );
    }

    #[test]
    fn test_blake2s_single_byte() {
        run_blake2s_test(
            b"a",
            &[
                0x4A, 0x0D, 0x12, 0x98, 0x73, 0x40, 0x30, 0x37, 0xC2, 0xCD, 0x9B, 0x90, 0x48,
                0x20, 0x36, 0x87, 0xF6, 0x23, 0x3F, 0xB6, 0x73, 0x89, 0x56, 0xE0, 0x34, 0x9B,
                0xD4, 0x32, 0x0F, 0xEC, 0x3E, 0x90,
            ],
        );
    }

    #[test]
    fn test_blake2s_short_string() {
        run_blake2s_test(
            b"abc",
            &[
                0x50, 0x8C, 0x5E, 0x8C, 0x32, 0x7C, 0x14, 0xE2, 0xE1, 0xA7, 0x2B, 0xA3, 0x4E,
                0xEB, 0x45, 0x2F, 0x37, 0x45, 0x8B, 0x20, 0x9E, 0xD6, 0x3A, 0x29, 0x4D, 0x99,
                0x9B, 0x4C, 0x86, 0x67, 0x59, 0x82,
            ],
        );
    }

    #[test]
    fn test_blake2s_one_block() {
        run_blake2s_test(
            b"abcdefghijklmnop",
            &[
                0xB6, 0x77, 0x5F, 0xD6, 0x8A, 0x7B, 0x03, 0xF1, 0x77, 0x42, 0x6A, 0x0E, 0xF1,
                0xAC, 0xEF, 0x97, 0xAC, 0x07, 0x0A, 0xE0, 0xD3, 0x30, 0xBB, 0x46, 0x2E, 0xBB,
                0x52, 0x93, 0x16, 0xA6, 0x1C, 0xF7,
            ],
        );
    }

    #[test]
    fn test_blake2s_multi_block() {
        run_blake2s_test(
            b"abcdefghijklmnopqrstuvwxyz0123456789abcdefghijklmnopqrstuvwxyz0123456789",
            &[
                0x44, 0xDD, 0xDB, 0x39, 0xBD, 0xB2, 0xAF, 0x80, 0xC1, 0x47, 0x89, 0x4C, 0x1D,
                0x75, 0x6A, 0xDA, 0x3D, 0x1C, 0x2A, 0xC2, 0xB1, 0x00, 0x54, 0x1E, 0x04, 0xFE,
                0x87, 0xB4, 0xA5, 0x9E, 0x12, 0x43,
            ],
        );
    }

    #[test]
    fn test_blake2s_native_consistency() {
        let builder = make_builder();
        let input = b"abcdefghijklmnopqrstuvwxyz";
        let input_arr = ByteArrayT::from_bytes(builder.clone(), input);

        let circuit_output = blake2s(input_arr);
        let native_output = bbrs_crypto::blake2s::blake2s(input);

        assert_eq!(circuit_output.size(), 32);
        assert_eq!(circuit_output.get_value(), native_output.to_vec());
        check_circuit(&builder).expect("circuit check failed");
    }
}
