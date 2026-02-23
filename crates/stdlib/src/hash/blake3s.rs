//! Circuit Blake3s hash function.
//!
//! Port of `barretenberg/stdlib/hash/blake3s/blake3s_plookup.{hpp,cpp}`.
//!
//! Implements the BLAKE3 hash function in-circuit using plookup tables
//! for XOR-rotation operations and field arithmetic for modular addition.
//!
//! Supports hashing of up to 1024 bytes (one chunk) and produces a 32-byte digest.
//! Does not implement tree-hashing mode, keyed hashing, or key derivation.

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

/// BLAKE3 initialization vector (same as SHA-256 H0 / BLAKE2s IV).
const IV: [u32; 8] = [
    0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
    0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
];

/// BLAKE3 domain separation flags.
const CHUNK_START: u8 = 1;
const CHUNK_END: u8 = 2;
const ROOT: u8 = 8;

/// BLAKE3 message schedule (7 rounds).
///
/// Each row gives the indices into the 16-word message block for that round.
/// Derived by repeatedly applying the BLAKE3 message permutation.
const MSG_SCHEDULE: [[usize; 16]; 7] = [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8],
    [3, 4, 10, 12, 13, 2, 7, 14, 6, 5, 9, 0, 11, 15, 8, 1],
    [10, 7, 12, 9, 14, 3, 13, 15, 4, 0, 11, 2, 5, 8, 1, 6],
    [12, 13, 9, 11, 15, 10, 14, 8, 7, 2, 5, 3, 0, 1, 6, 4],
    [9, 14, 11, 5, 8, 12, 15, 1, 13, 3, 0, 10, 2, 6, 4, 7],
    [11, 15, 5, 0, 1, 9, 8, 6, 14, 10, 2, 12, 3, 4, 7, 13],
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
    result_wit.create_range_constraint(32, "blake3s add_mod_32: result");

    let carry_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(carry)));
    carry_wit.create_range_constraint(1, "blake3s add_mod_32: carry");

    let shift = FieldT::from_field(Fr::from(1u64 << 32));
    let reconstructed = result_wit.clone() + carry_wit * shift;
    let sum_field = a.clone() + b.clone();
    sum_field.assert_equal(&reconstructed, "blake3s add_mod_32: constraint");

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
    result_wit.create_range_constraint(32, "blake3s add_mod_32_three: result");

    let carry_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(carry)));
    carry_wit.create_range_constraint(2, "blake3s add_mod_32_three: carry");

    let shift = FieldT::from_field(Fr::from(1u64 << 32));
    let reconstructed = result_wit.clone() + carry_wit * shift;
    let sum_field = a.add_two(b, c);
    sum_field.assert_equal(&reconstructed, "blake3s add_mod_32_three: constraint");

    result_wit
}

// ════════════════════════════════════════════════════════════════════════
//  Helpers: XOR with rotation
// ════════════════════════════════════════════════════════════════════════

/// XOR two 32-bit values and rotate right by the specified amount.
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

    let lo_bits = rotation;
    let hi_bits = 32 - rotation;
    let lo_mask = (1u32 << lo_bits) - 1;

    let lo = native_xor & lo_mask;
    let hi = native_xor >> lo_bits;

    let lo_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(lo as u64)));
    lo_wit.create_range_constraint(lo_bits as usize, "blake3s xor_rotate: lo");

    let hi_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(hi as u64)));
    hi_wit.create_range_constraint(hi_bits as usize, "blake3s xor_rotate: hi");

    let shift_lo = FieldT::from_field(Fr::from(1u64 << lo_bits));
    let reconstructed = lo_wit.clone() + hi_wit.clone() * shift_lo;
    xor_val.assert_equal(&reconstructed, "blake3s xor_rotate: decomposition");

    let shift_hi = FieldT::from_field(Fr::from(1u64 << hi_bits));
    hi_wit + lo_wit * shift_hi
}

// ════════════════════════════════════════════════════════════════════════
//  G mixing function
// ════════════════════════════════════════════════════════════════════════

/// BLAKE3 G mixing function (same rotation constants as BLAKE2s: 16, 12, 8, 7).
fn blake3s_g(
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
//  Compression functions
// ════════════════════════════════════════════════════════════════════════

/// BLAKE3 core compression (pre-finalization).
///
/// Initializes the 16-word state from the chaining value, IV, counter,
/// block length, and flags, then runs 7 rounds of mixing.
fn compress_pre(
    cv: &[FieldT<P>; 8],
    msg: &[FieldT<P>; 16],
    counter: u64,
    block_len: u32,
    flags: u8,
    ctx: &BuilderRef<P>,
) -> [FieldT<P>; 16] {
    let counter_lo = (counter & 0xFFFF_FFFF) as u64;
    let counter_hi = (counter >> 32) as u64;

    let mut state: [FieldT<P>; 16] = [
        cv[0].clone(),
        cv[1].clone(),
        cv[2].clone(),
        cv[3].clone(),
        cv[4].clone(),
        cv[5].clone(),
        cv[6].clone(),
        cv[7].clone(),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[0] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[1] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[2] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[3] as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(counter_lo)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(counter_hi)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(block_len as u64)),
        FieldT::constant_with_context(ctx.clone(), Fr::from(flags as u64)),
    ];

    // 7 rounds of mixing with BLAKE3 message schedule
    for round in 0..7 {
        let s = &MSG_SCHEDULE[round];

        // Column step
        blake3s_g(&mut state, 0, 4, 8, 12, &msg[s[0]], &msg[s[1]]);
        blake3s_g(&mut state, 1, 5, 9, 13, &msg[s[2]], &msg[s[3]]);
        blake3s_g(&mut state, 2, 6, 10, 14, &msg[s[4]], &msg[s[5]]);
        blake3s_g(&mut state, 3, 7, 11, 15, &msg[s[6]], &msg[s[7]]);

        // Diagonal step
        blake3s_g(&mut state, 0, 5, 10, 15, &msg[s[8]], &msg[s[9]]);
        blake3s_g(&mut state, 1, 6, 11, 12, &msg[s[10]], &msg[s[11]]);
        blake3s_g(&mut state, 2, 7, 8, 13, &msg[s[12]], &msg[s[13]]);
        blake3s_g(&mut state, 3, 4, 9, 14, &msg[s[14]], &msg[s[15]]);
    }

    state
}

/// BLAKE3 compress-in-place: updates the chaining value.
///
/// cv[i] = state[i] ^ state[i+8] for i in 0..8
fn compress_in_place(
    cv: &mut [FieldT<P>; 8],
    msg: &[FieldT<P>; 16],
    counter: u64,
    block_len: u32,
    flags: u8,
    ctx: &BuilderRef<P>,
) {
    let state = compress_pre(cv, msg, counter, block_len, flags, ctx);

    for i in 0..8 {
        cv[i] = plookup::read_from_2_to_1_table(
            MultiTableId::UINT32_XOR,
            &state[i],
            &state[i + 8],
        );
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Byte <-> word conversion helpers
// ════════════════════════════════════════════════════════════════════════

/// Convert 4 consecutive bytes (at `offset`) to a 32-bit little-endian word.
fn bytes_to_word(bytes: &[FieldT<P>], offset: usize) -> FieldT<P> {
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
    b0_wit.create_range_constraint(8, "blake3s word_to_bytes: b0");
    let b1_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(b1)));
    b1_wit.create_range_constraint(8, "blake3s word_to_bytes: b1");
    let b2_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(b2)));
    b2_wit.create_range_constraint(8, "blake3s word_to_bytes: b2");
    let b3_wit = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Fr::from(b3)));
    b3_wit.create_range_constraint(8, "blake3s word_to_bytes: b3");

    let s8 = FieldT::<P>::from_field(Fr::from(256u64));
    let s16 = FieldT::<P>::from_field(Fr::from(65536u64));
    let s24 = FieldT::<P>::from_field(Fr::from(16777216u64));
    let reconstructed = FieldT::accumulate(&[
        b0_wit.clone(),
        b1_wit.clone() * s8,
        b2_wit.clone() * s16,
        b3_wit.clone() * s24,
    ]);
    word.assert_equal(&reconstructed, "blake3s word_to_bytes: reconstruction");

    [b0_wit, b1_wit, b2_wit, b3_wit]
}

// ════════════════════════════════════════════════════════════════════════
//  Public API
// ════════════════════════════════════════════════════════════════════════

/// Compute the BLAKE3 hash of a byte array in-circuit.
///
/// Takes a byte array of up to 1024 bytes and returns a 32-byte hash result.
/// Each output byte is range-constrained to 8 bits.
///
/// This implements single-chunk BLAKE3 hashing (no tree mode).
pub fn blake3s(input: ByteArrayT<P>) -> ByteArrayT<P> {
    let ctx = input.get_context().as_ref().unwrap().clone();
    let input_bytes = input.bytes();
    let input_len = input_bytes.len();

    assert!(
        input_len <= 1024,
        "blake3s: input must be at most 1024 bytes (one chunk)"
    );

    // Initialize chaining value with IV (no parameter block XOR, unlike BLAKE2s)
    let mut cv: [FieldT<P>; 8] = std::array::from_fn(|i| {
        FieldT::constant_with_context(ctx.clone(), Fr::from(IV[i] as u64))
    });

    // Determine number of 64-byte blocks.
    // Empty input still requires one block (zero-padded, block_len=0).
    let num_blocks = if input_len == 0 {
        1
    } else {
        (input_len + 63) / 64
    };

    // Pad input to fill all blocks with zeros
    let padded_len = num_blocks * 64;
    let mut padded: Vec<FieldT<P>> = Vec::with_capacity(padded_len);
    padded.extend_from_slice(input_bytes);
    for _ in input_len..padded_len {
        padded.push(FieldT::constant_with_context(ctx.clone(), Fr::zero()));
    }

    // Process blocks.
    // For single-chunk BLAKE3, the counter is always 0 (chunk counter).
    // blocks_compressed tracks how many blocks have been processed, used only
    // for the CHUNK_START flag (set on the first block only).
    let chunk_counter: u64 = 0;
    let mut blocks_compressed: u64 = 0;

    for block_idx in 0..num_blocks {
        let is_last = block_idx == num_blocks - 1;
        let block_start = block_idx * 64;

        // Convert 64 bytes to 16 32-bit words (little-endian)
        let msg: [FieldT<P>; 16] =
            std::array::from_fn(|i| bytes_to_word(&padded, block_start + i * 4));

        if is_last {
            // Final block: compute flags and block length
            let mut flags: u8 = CHUNK_END | ROOT;
            if blocks_compressed == 0 {
                flags |= CHUNK_START;
            }

            let block_len = if input_len == 0 {
                0u32
            } else {
                (input_len - block_start) as u32
            };

            // Compress and produce output
            let state = compress_pre(&cv, &msg, chunk_counter, block_len, flags, &ctx);

            // Output: state[i] ^ state[i+8] for i in 0..8, converted to bytes
            let mut output_bytes: Vec<FieldT<P>> = Vec::with_capacity(32);
            for i in 0..8 {
                let word = plookup::read_from_2_to_1_table(
                    MultiTableId::UINT32_XOR,
                    &state[i],
                    &state[i + 8],
                );
                let wb = word_to_bytes(&ctx, &word);
                output_bytes.extend(wb);
            }

            return ByteArrayT::from_values(Some(ctx), output_bytes);
        } else {
            // Intermediate block: compress in place to update chaining value
            let mut flags: u8 = 0;
            if blocks_compressed == 0 {
                flags |= CHUNK_START;
            }

            compress_in_place(&mut cv, &msg, chunk_counter, 64, flags, &ctx);
            blocks_compressed += 1;
        }
    }

    unreachable!()
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
    fn run_blake3s_test(input: &[u8], expected: &[u8; 32]) {
        let builder = make_builder();
        let input_arr = ByteArrayT::from_bytes(builder.clone(), input);
        let output = blake3s(input_arr);

        assert_eq!(output.size(), 32);
        assert_eq!(
            output.get_value(),
            expected.to_vec(),
            "hash mismatch for input len={}",
            input.len()
        );
        check_circuit(&builder).expect("circuit check failed");
    }

    #[test]
    fn test_blake3s_empty() {
        let expected = bbrs_crypto::blake3s::blake3s(b"");
        run_blake3s_test(b"", &expected);
    }

    #[test]
    fn test_blake3s_single_byte() {
        let input = b"\x00";
        let expected = bbrs_crypto::blake3s::blake3s(input);
        run_blake3s_test(input, &expected);
    }

    #[test]
    fn test_blake3s_short_string() {
        let input = b"abc";
        let expected = bbrs_crypto::blake3s::blake3s(input);
        run_blake3s_test(input, &expected);
    }

    #[test]
    fn test_blake3s_one_block() {
        let input = b"abcdefghijklmnop";
        let expected = bbrs_crypto::blake3s::blake3s(input);
        run_blake3s_test(input, &expected);
    }

    #[test]
    fn test_blake3s_multi_block() {
        // 65 bytes: spans two 64-byte blocks
        let input: Vec<u8> = (0..65).map(|i| (i % 251) as u8).collect();
        let expected = bbrs_crypto::blake3s::blake3s(&input);
        run_blake3s_test(&input, &expected);
    }

    #[test]
    fn test_blake3s_native_consistency() {
        let builder = make_builder();
        let input = b"abcdefghijklmnopqrstuvwxyz";
        let input_arr = ByteArrayT::from_bytes(builder.clone(), input);

        let circuit_output = blake3s(input_arr);
        let native_output = bbrs_crypto::blake3s::blake3s(input);

        assert_eq!(circuit_output.size(), 32);
        assert_eq!(circuit_output.get_value(), native_output.to_vec());
        check_circuit(&builder).expect("circuit check failed");
    }

    #[test]
    fn test_blake3s_edge_case_addition_overflow() {
        // Edge case from C++ tests: input bytes that trigger addition overflow
        // in the G function. Tests that normalization handles overflow correctly.
        let input: Vec<u8> = vec![
            0xC3, 0x2B, 0xC3, 0x91, 0x23, 0xFF, 0xFF, 0xFF, 0xFF, 0xC3, 0xFF, 0xFF, 0xFF, 0xFF,
            0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
            0xC3, 0x03, 0x83, 0x83, 0x83, 0x40,
        ];
        let expected = bbrs_crypto::blake3s::blake3s(&input);
        run_blake3s_test(&input, &expected);
    }
}
