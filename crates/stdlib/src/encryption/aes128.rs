//! In-circuit AES-128-CBC encryption.
//!
//! Port of C++ `barretenberg/stdlib/encryption/aes128/aes128.{hpp,cpp}`.
//!
//! Implements AES-128-CBC encryption using plookup tables for in-circuit
//! evaluation. The core operations (SubBytes, ShiftRows, MixColumns,
//! AddRoundKey) are computed in sparse base-9 form where XOR becomes
//! addition followed by normalization.

use bbrs_circuit_builder::plookup_tables::types::{ColumnIdx, MultiTableId};
use bbrs_ecc::curves::bn254::Bn254FrParams;
use bbrs_ecc::fields::field::Field;

use crate::primitives::field::FieldT;
use crate::primitives::plookup;
use crate::primitives::witness::BuilderRef;

type P = Bn254FrParams;
type Fr = Field<P>;

// ════════════════════════════════════════════════════════════════════════
//  Constants
// ════════════════════════════════════════════════════════════════════════

const AES128_BASE: u64 = 9;

/// Convert a u128 into an Fr field element.
#[inline]
fn fr_from_u128(val: u128) -> Fr {
    Fr::from_limbs([val as u64, (val >> 64) as u64, 0, 0])
}
const BLOCK_SIZE: usize = 16;
const EXTENDED_KEY_LENGTH: usize = 176; // 11 round keys * 16 bytes
const NUM_ROUNDS: usize = 10;
const COLUMN_SIZE: usize = 4;

// ════════════════════════════════════════════════════════════════════════
//  Sparse form helpers
// ════════════════════════════════════════════════════════════════════════

/// Map an 8-bit value into sparse base-9 form.
fn map_into_sparse_form_base9(input: u8) -> u64 {
    let mut out: u64 = 0;
    let mut base_power: u64 = 1;
    for i in 0..8 {
        if (input >> i) & 1 != 0 {
            out += base_power;
        }
        if i < 7 {
            base_power *= AES128_BASE;
        }
    }
    out
}

/// Map a sparse base-9 value back to an 8-bit value.
fn map_from_sparse_form_base9(input: u64) -> u8 {
    let mut target = input;
    let mut output: u8 = 0;
    for i in 0..8 {
        let digit = target % AES128_BASE;
        if (digit & 1) == 1 {
            output |= 1u8 << i;
        }
        target /= AES128_BASE;
    }
    output
}

// ════════════════════════════════════════════════════════════════════════
//  Internal type aliases
// ════════════════════════════════════════════════════════════════════════

/// A byte in the AES state: `.0` is the S-box output (×1), `.1` is the
/// S-box output XOR xtime (×3), used by MixColumns.
type BytePair = (FieldT<P>, FieldT<P>);

// ════════════════════════════════════════════════════════════════════════
//  Plookup wrappers
// ════════════════════════════════════════════════════════════════════════

/// Normalize a sparse-form byte via the AES_NORMALIZE plookup table.
fn normalize_sparse_form(byte: &FieldT<P>) -> FieldT<P> {
    plookup::read_from_1_to_2_table(MultiTableId::AES_NORMALIZE, byte)
}

/// Apply the AES S-box via the AES_SBOX plookup table.
/// Returns `(S(x), S(x) ^ xtime(S(x)))`.
fn apply_aes_sbox_map(input: &FieldT<P>) -> BytePair {
    plookup::read_pair_from_table(MultiTableId::AES_SBOX, input)
}

// ════════════════════════════════════════════════════════════════════════
//  Sparse form conversions
// ════════════════════════════════════════════════════════════════════════

/// Convert a 128-bit block into 16 sparse-form bytes via the AES_INPUT plookup table.
pub fn convert_into_sparse_bytes(
    ctx: BuilderRef<P>,
    block_data: &FieldT<P>,
) -> [FieldT<P>; BLOCK_SIZE] {
    let mut block_data_copy = block_data.clone();
    if block_data_copy.is_constant() {
        block_data_copy.convert_constant_to_fixed_witness(ctx);
    }
    let lookup = plookup::get_lookup_accumulators(
        MultiTableId::AES_INPUT,
        &block_data_copy,
        &FieldT::default(),
        false,
    );
    let mut sparse_bytes: [FieldT<P>; BLOCK_SIZE] = std::array::from_fn(|_| FieldT::default());
    for i in 0..BLOCK_SIZE {
        sparse_bytes[BLOCK_SIZE - 1 - i] = lookup.columns[ColumnIdx::C2 as usize][i].clone();
    }
    sparse_bytes
}

/// Convert 16 sparse-form bytes back to a 128-bit field element.
pub fn convert_from_sparse_bytes(
    ctx: BuilderRef<P>,
    sparse_bytes: &[FieldT<P>; BLOCK_SIZE],
) -> FieldT<P> {
    // Compute the native value by decoding each sparse byte
    let mut accumulator: u128 = 0;
    for i in 0..BLOCK_SIZE {
        let sparse_byte = {
            let standard = sparse_bytes[i].get_value().from_montgomery_form();
            standard.data[0]
        };
        let byte = map_from_sparse_form_base9(sparse_byte) as u128;
        accumulator = (accumulator << 8) + byte;
    }

    // Create witness for the result and constrain via AES_INPUT lookup
    let result = FieldT::from_witness(ctx, fr_from_u128(accumulator));
    let lookup = plookup::get_lookup_accumulators(
        MultiTableId::AES_INPUT,
        &result,
        &FieldT::default(),
        false,
    );
    for i in 0..BLOCK_SIZE {
        sparse_bytes[BLOCK_SIZE - 1 - i].assert_equal(
            &lookup.columns[ColumnIdx::C2 as usize][i],
            "convert_from_sparse_bytes: sparse byte mismatch",
        );
    }
    result
}

// ════════════════════════════════════════════════════════════════════════
//  Key expansion
// ════════════════════════════════════════════════════════════════════════

/// AES-128 key expansion (FIPS 197 Section 5.2).
///
/// Expands a 128-bit key into 176 sparse-form bytes (11 round keys).
fn expand_key(ctx: BuilderRef<P>, key: &FieldT<P>) -> [FieldT<P>; EXTENDED_KEY_LENGTH] {
    const ROUND_CONSTANTS: [u8; 11] = [0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36];

    let sparse_round_constants: [FieldT<P>; 11] = std::array::from_fn(|i| {
        FieldT::constant_with_context(
            ctx.clone(),
            Fr::from(map_into_sparse_form_base9(ROUND_CONSTANTS[i])),
        )
    });

    let mut round_key: [FieldT<P>; EXTENDED_KEY_LENGTH] =
        std::array::from_fn(|_| FieldT::default());
    let sparse_key = convert_into_sparse_bytes(ctx.clone(), key);

    let mut add_counts = [1u64; EXTENDED_KEY_LENGTH];

    // First 16 bytes = original key
    for i in 0..16 {
        round_key[i] = sparse_key[i].clone();
    }

    // Iterate over the 40 words (4 words per round for 10 rounds)
    for i in 4..44 {
        let k = (i - 1) * 4;
        let mut temp: [FieldT<P>; 4] = [
            round_key[k].clone(),
            round_key[k + 1].clone(),
            round_key[k + 2].clone(),
            round_key[k + 3].clone(),
        ];
        let mut temp_add_counts: [u64; 4] = [
            add_counts[k],
            add_counts[k + 1],
            add_counts[k + 2],
            add_counts[k + 3],
        ];

        if (i & 0x03) == 0 {
            // RotWord
            temp.rotate_left(1);
            temp_add_counts.rotate_left(1);

            // SubWord via S-box lookup
            for t in temp.iter_mut() {
                *t = apply_aes_sbox_map(t).0;
            }

            // Add round constant to first byte
            temp[0] = &temp[0] + &sparse_round_constants[i >> 2];
            temp_add_counts[0] += 1;
        }

        let j = i * 4;
        let prev_k = (i - 4) * 4;
        for m in 0..4 {
            round_key[j + m] = &round_key[prev_k + m] + &temp[m];
            add_counts[j + m] = add_counts[prev_k + m] + temp_add_counts[m];
        }

        // Normalize when accumulated additions exceed threshold
        const TARGET: u64 = 3;
        for m in 0..4 {
            let byte_index = j + m;
            if add_counts[byte_index] > TARGET
                || (add_counts[byte_index] > 1 && (byte_index & 12) == 12)
            {
                round_key[byte_index] = normalize_sparse_form(&round_key[byte_index]);
                add_counts[byte_index] = 1;
            }
        }
    }

    round_key
}

// ════════════════════════════════════════════════════════════════════════
//  AES round operations
// ════════════════════════════════════════════════════════════════════════

/// ShiftRows (FIPS 197 Section 5.1.2).
fn shift_rows(state: &mut [BytePair; BLOCK_SIZE]) {
    // Row 1: shift left by 1
    let temp = state[1].clone();
    state[1] = state[5].clone();
    state[5] = state[9].clone();
    state[9] = state[13].clone();
    state[13] = temp;

    // Row 2: shift left by 2
    let temp = state[2].clone();
    state[2] = state[10].clone();
    state[10] = temp;
    let temp = state[6].clone();
    state[6] = state[14].clone();
    state[14] = temp;

    // Row 3: shift left by 3 (= right by 1)
    let temp = state[3].clone();
    state[3] = state[15].clone();
    state[15] = state[11].clone();
    state[11] = state[7].clone();
    state[7] = temp;
}

/// MixColumns on a single column + AddRoundKey (FIPS 197 Sections 5.1.3 & 5.1.4).
///
/// Uses the byte_pair structure:
///   - `.0` = S(x)                    (×1)
///   - `.1` = S(x) ^ xtime(S(x))     (×3, precomputed)
///
/// 2·x = 3·x ^ x = `.1` + `.0` (in sparse form, addition = XOR before normalization)
fn mix_column_and_add_round_key(
    column: &mut [BytePair; COLUMN_SIZE],
    round_key: &[FieldT<P>; EXTENDED_KEY_LENGTH],
    round: usize,
    col_idx: usize,
) {
    // t0 = s0 + s3 + 3·s1
    let t0 = column[0].0.add_two(&column[3].0, &column[1].1);
    // t1 = s1 + s2 + 3·s3
    let t1 = column[1].0.add_two(&column[2].0, &column[3].1);

    // r0 = 2·s0 + 3·s1 + s2 + s3 = t0 + s2 + 3·s0
    let r0 = t0.add_two(&column[2].0, &column[0].1);
    // r1 = s0 + 2·s1 + 3·s2 + s3 = t0 + s1 + 3·s2
    let r1 = t0.add_two(&column[1].0, &column[2].1);
    // r2 = s0 + s1 + 2·s2 + 3·s3 = t1 + s0 + 3·s2
    let r2 = t1.add_two(&column[0].0, &column[2].1);
    // r3 = 3·s0 + s1 + s2 + 2·s3 = t1 + 3·s0 + s3
    let r3 = t1.add_two(&column[0].1, &column[3].0);

    let key_offset = round * BLOCK_SIZE + col_idx * COLUMN_SIZE;

    column[0].0 = &r0 + &round_key[key_offset];
    column[1].0 = &r1 + &round_key[key_offset + 1];
    column[2].0 = &r2 + &round_key[key_offset + 2];
    column[3].0 = &r3 + &round_key[key_offset + 3];
}

/// MixColumns + AddRoundKey for all 4 columns.
fn mix_columns_and_add_round_key(
    state: &mut [BytePair; BLOCK_SIZE],
    round_key: &[FieldT<P>; EXTENDED_KEY_LENGTH],
    round: usize,
) {
    for col in 0..4 {
        let base = col * COLUMN_SIZE;
        let mut column: [BytePair; COLUMN_SIZE] = [
            state[base].clone(),
            state[base + 1].clone(),
            state[base + 2].clone(),
            state[base + 3].clone(),
        ];
        mix_column_and_add_round_key(&mut column, round_key, round, col);
        for j in 0..COLUMN_SIZE {
            state[base + j] = column[j].clone();
        }
    }
}

/// SubBytes: apply S-box to each byte of the state.
fn sub_bytes(state: &mut [BytePair; BLOCK_SIZE]) {
    for i in 0..BLOCK_SIZE {
        state[i] = apply_aes_sbox_map(&state[i].0);
    }
}

/// AddRoundKey: XOR state with round key (via sparse-form addition).
fn add_round_key(
    state: &mut [BytePair; BLOCK_SIZE],
    round_key: &[FieldT<P>; EXTENDED_KEY_LENGTH],
    round: usize,
) {
    let key_offset = round * BLOCK_SIZE;
    for i in 0..BLOCK_SIZE {
        state[i].0 = &state[i].0 + &round_key[key_offset + i];
    }
}

/// XOR state with IV (via sparse-form addition).
fn xor_with_iv(state: &mut [BytePair; BLOCK_SIZE], iv: &[FieldT<P>; BLOCK_SIZE]) {
    for i in 0..BLOCK_SIZE {
        state[i].0 = &state[i].0 + &iv[i];
    }
}

/// Single-block AES-128 cipher.
fn aes128_cipher(
    state: &mut [BytePair; BLOCK_SIZE],
    round_key: &[FieldT<P>; EXTENDED_KEY_LENGTH],
) {
    // Initial round key addition
    add_round_key(state, round_key, 0);
    for i in 0..BLOCK_SIZE {
        state[i].0 = normalize_sparse_form(&state[i].0);
    }

    // Rounds 1-9: SubBytes, ShiftRows, MixColumns, AddRoundKey
    for round in 1..NUM_ROUNDS {
        sub_bytes(state);
        shift_rows(state);
        mix_columns_and_add_round_key(state, round_key, round);
        for i in 0..BLOCK_SIZE {
            state[i].0 = normalize_sparse_form(&state[i].0);
        }
    }

    // Final round (no MixColumns)
    sub_bytes(state);
    shift_rows(state);
    add_round_key(state, round_key, NUM_ROUNDS);
}

// ════════════════════════════════════════════════════════════════════════
//  Public interface
// ════════════════════════════════════════════════════════════════════════

/// AES-128-CBC encryption in-circuit.
///
/// Each element in `input` is a 128-bit block packed into a field element (big-endian).
/// `iv` and `key` are also 128-bit field elements.
///
/// Returns the ciphertext as a vector of 128-bit field elements.
pub fn encrypt_buffer_cbc(
    input: &[FieldT<P>],
    iv: &FieldT<P>,
    key: &FieldT<P>,
) -> Vec<FieldT<P>> {
    // Check if all inputs are constants
    let all_constants = key.is_constant()
        && iv.is_constant()
        && input.iter().all(|b| b.is_constant());

    if all_constants {
        return encrypt_buffer_cbc_constant(input, iv, key);
    }

    // Find a valid builder context
    let ctx = key
        .get_context()
        .clone()
        .or_else(|| iv.get_context().clone())
        .or_else(|| {
            input
                .iter()
                .find_map(|b| b.get_context().clone())
        })
        .expect("at least one input must have a builder context");

    let round_key = expand_key(ctx.clone(), key);
    let num_blocks = input.len();

    // Convert all input blocks to sparse form
    let mut sparse_state: Vec<BytePair> = Vec::with_capacity(num_blocks * BLOCK_SIZE);
    for block in input {
        let bytes = convert_into_sparse_bytes(ctx.clone(), block);
        for byte in bytes {
            sparse_state.push((byte, FieldT::constant_with_context(ctx.clone(), Fr::zero())));
        }
    }

    let mut sparse_iv = convert_into_sparse_bytes(ctx.clone(), iv);

    // Encrypt each block in CBC mode
    for i in 0..num_blocks {
        let base = i * BLOCK_SIZE;
        let mut round_state: [BytePair; BLOCK_SIZE] = std::array::from_fn(|j| {
            sparse_state[base + j].clone()
        });

        xor_with_iv(&mut round_state, &sparse_iv);
        aes128_cipher(&mut round_state, &round_key);

        // Update IV for next block (ciphertext feeds forward)
        for j in 0..BLOCK_SIZE {
            sparse_iv[j] = round_state[j].0.clone();
            sparse_state[base + j] = round_state[j].clone();
        }
    }

    // Normalize and convert back to field elements
    let mut sparse_output: Vec<FieldT<P>> = Vec::with_capacity(num_blocks * BLOCK_SIZE);
    for pair in &sparse_state {
        sparse_output.push(normalize_sparse_form(&pair.0));
    }

    let mut output = Vec::with_capacity(num_blocks);
    for i in 0..num_blocks {
        let block: [FieldT<P>; BLOCK_SIZE] = std::array::from_fn(|j| {
            sparse_output[i * BLOCK_SIZE + j].clone()
        });
        output.push(convert_from_sparse_bytes(ctx.clone(), &block));
    }
    output
}

/// Constant-only path: use native crypto to compute AES-128-CBC.
fn encrypt_buffer_cbc_constant(
    input: &[FieldT<P>],
    iv: &FieldT<P>,
    key: &FieldT<P>,
) -> Vec<FieldT<P>> {
    // Extract key bytes
    let key_value = key.get_value().from_montgomery_form();
    let mut key_bytes = [0u8; 16];
    for i in 0..16 {
        key_bytes[15 - i] = ((key_value.data[i / 8] >> ((i % 8) * 8)) & 0xFF) as u8;
    }

    // Extract IV bytes
    let iv_value = iv.get_value().from_montgomery_form();
    let mut iv_bytes = [0u8; 16];
    for i in 0..16 {
        iv_bytes[15 - i] = ((iv_value.data[i / 8] >> ((i % 8) * 8)) & 0xFF) as u8;
    }

    // Extract input bytes
    let mut input_bytes = vec![0u8; input.len() * 16];
    for (block_idx, block) in input.iter().enumerate() {
        let block_value = block.get_value().from_montgomery_form();
        for i in 0..16 {
            input_bytes[block_idx * 16 + 15 - i] =
                ((block_value.data[i / 8] >> ((i % 8) * 8)) & 0xFF) as u8;
        }
    }

    // Run native AES encryption
    let ciphertext = bbrs_crypto::aes128::encrypt_buffer_cbc(&input_bytes, &iv_bytes, &key_bytes);

    // Convert result back to field elements
    let mut result = Vec::with_capacity(input.len());
    for block_idx in 0..input.len() {
        let mut value: u128 = 0;
        for i in 0..16 {
            value = (value << 8) | (ciphertext[block_idx * 16 + i] as u128);
        }
        result.push(FieldT::from_field(fr_from_u128(value)));
    }
    result
}

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

    /// Convert 16 bytes to a u128 value (big-endian).
    fn convert_bytes_to_u128(data: &[u8]) -> u128 {
        let mut result: u128 = 0;
        for &b in &data[..16] {
            result = (result << 8) | (b as u128);
        }
        result
    }

    /// Create a FieldT as either constant or witness.
    fn create_field_element(
        builder: BuilderRef<P>,
        value: u128,
        as_witness: bool,
    ) -> FieldT<P> {
        if as_witness {
            FieldT::from_witness(builder, fr_from_u128(value))
        } else {
            FieldT::from_field(fr_from_u128(value))
        }
    }

    // NIST test vectors for AES-128-CBC
    const KEY: [u8; 16] = [
        0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6,
        0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c,
    ];
    const IV: [u8; 16] = [
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
        0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
    ];
    const PLAINTEXT: [u8; 64] = [
        0x6b, 0xc1, 0xbe, 0xe2, 0x2e, 0x40, 0x9f, 0x96,
        0xe9, 0x3d, 0x7e, 0x11, 0x73, 0x93, 0x17, 0x2a,
        0xae, 0x2d, 0x8a, 0x57, 0x1e, 0x03, 0xac, 0x9c,
        0x9e, 0xb7, 0x6f, 0xac, 0x45, 0xaf, 0x8e, 0x51,
        0x30, 0xc8, 0x1c, 0x46, 0xa3, 0x5c, 0xe4, 0x11,
        0xe5, 0xfb, 0xc1, 0x19, 0x1a, 0x0a, 0x52, 0xef,
        0xf6, 0x9f, 0x24, 0x45, 0xdf, 0x4f, 0x9b, 0x17,
        0xad, 0x2b, 0x41, 0x7b, 0xe6, 0x6c, 0x37, 0x10,
    ];
    const CIPHERTEXT: [u8; 64] = [
        0x76, 0x49, 0xab, 0xac, 0x81, 0x19, 0xb2, 0x46,
        0xce, 0xe9, 0x8e, 0x9b, 0x12, 0xe9, 0x19, 0x7d,
        0x50, 0x86, 0xcb, 0x9b, 0x50, 0x72, 0x19, 0xee,
        0x95, 0xdb, 0x11, 0x3a, 0x91, 0x76, 0x78, 0xb2,
        0x73, 0xbe, 0xd6, 0xb8, 0xe3, 0xc1, 0x74, 0x3b,
        0x71, 0x16, 0xe6, 0x9e, 0x22, 0x22, 0x95, 0x16,
        0x3f, 0xf1, 0xca, 0xa1, 0x68, 0x1f, 0xac, 0x09,
        0x12, 0x0e, 0xca, 0x30, 0x75, 0x86, 0xe1, 0xa7,
    ];

    /// Get expected ciphertext blocks as Fr values.
    fn expected_blocks() -> [Fr; 4] {
        [
            fr_from_u128(convert_bytes_to_u128(&CIPHERTEXT[0..])),
            fr_from_u128(convert_bytes_to_u128(&CIPHERTEXT[16..])),
            fr_from_u128(convert_bytes_to_u128(&CIPHERTEXT[32..])),
            fr_from_u128(convert_bytes_to_u128(&CIPHERTEXT[48..])),
        ]
    }

    /// Run AES-128-CBC test with specified witness/constant configuration.
    fn test_aes128_combination(
        key_as_witness: bool,
        iv_as_witness: bool,
        input_as_witness: bool,
    ) {
        let builder = make_builder();

        let mut in_field = Vec::with_capacity(4);
        for i in 0..4 {
            in_field.push(create_field_element(
                builder.clone(),
                convert_bytes_to_u128(&PLAINTEXT[i * 16..]),
                input_as_witness,
            ));
        }

        let key_field = create_field_element(
            builder.clone(),
            convert_bytes_to_u128(&KEY),
            key_as_witness,
        );
        let iv_field = create_field_element(
            builder.clone(),
            convert_bytes_to_u128(&IV),
            iv_as_witness,
        );

        let expected = expected_blocks();
        let result = encrypt_buffer_cbc(&in_field, &iv_field, &key_field);

        for i in 0..4 {
            assert_eq!(
                result[i].get_value(),
                expected[i],
                "Block {} mismatch (key={}, iv={}, input={})",
                i,
                if key_as_witness { "witness" } else { "constant" },
                if iv_as_witness { "witness" } else { "constant" },
                if input_as_witness { "witness" } else { "constant" },
            );
        }

        if key_as_witness || iv_as_witness || input_as_witness {
            UltraCircuitChecker::check(&mut builder.borrow_mut())
                .expect("circuit check failed");
        }
    }

    /// Run AES-128-CBC test with mixed per-block witness/constant configuration.
    fn test_aes128_mixed_input(
        key_as_witness: bool,
        iv_as_witness: bool,
        block_config: &[bool],
    ) {
        let builder = make_builder();

        let mut in_field = Vec::with_capacity(block_config.len());
        for (i, &as_witness) in block_config.iter().enumerate() {
            in_field.push(create_field_element(
                builder.clone(),
                convert_bytes_to_u128(&PLAINTEXT[i * 16..]),
                as_witness,
            ));
        }

        let key_field = create_field_element(
            builder.clone(),
            convert_bytes_to_u128(&KEY),
            key_as_witness,
        );
        let iv_field = create_field_element(
            builder.clone(),
            convert_bytes_to_u128(&IV),
            iv_as_witness,
        );

        let expected = expected_blocks();
        let result = encrypt_buffer_cbc(&in_field, &iv_field, &key_field);

        for i in 0..result.len() {
            assert_eq!(result[i].get_value(), expected[i], "Block {} mismatch", i);
        }

        let has_witness =
            key_as_witness || iv_as_witness || block_config.iter().any(|&w| w);

        if has_witness {
            UltraCircuitChecker::check(&mut builder.borrow_mut())
                .expect("circuit check failed");
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Test 1: All witness
    // ════════════════════════════════════════════════════════════════════

    #[test]
    fn encrypt_64_bytes_all_witness() {
        test_aes128_combination(true, true, true);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Test 2: All constant
    // ════════════════════════════════════════════════════════════════════

    #[test]
    fn encrypt_64_bytes_all_constant() {
        test_aes128_combination(false, false, false);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Test 3: Key witness, IV/input constant
    // ════════════════════════════════════════════════════════════════════

    #[test]
    fn encrypt_64_bytes_key_witness_iv_constant_input_constant() {
        test_aes128_combination(true, false, false);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Test 4: IV witness, key/input constant
    // ════════════════════════════════════════════════════════════════════

    #[test]
    fn encrypt_64_bytes_key_constant_iv_witness_input_constant() {
        test_aes128_combination(false, true, false);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Test 5: Input witness, key/IV constant
    // ════════════════════════════════════════════════════════════════════

    #[test]
    fn encrypt_64_bytes_key_constant_iv_constant_input_witness() {
        test_aes128_combination(false, false, true);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Test 6: Key+IV witness, input constant
    // ════════════════════════════════════════════════════════════════════

    #[test]
    fn encrypt_64_bytes_key_witness_iv_witness_input_constant() {
        test_aes128_combination(true, true, false);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Test 7: Key+input witness, IV constant
    // ════════════════════════════════════════════════════════════════════

    #[test]
    fn encrypt_64_bytes_key_witness_iv_constant_input_witness() {
        test_aes128_combination(true, false, true);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Test 8: IV+input witness, key constant
    // ════════════════════════════════════════════════════════════════════

    #[test]
    fn encrypt_64_bytes_key_constant_iv_witness_input_witness() {
        test_aes128_combination(false, true, true);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Test 9: Original test (backward compatibility)
    // ════════════════════════════════════════════════════════════════════

    #[test]
    fn encrypt_64_bytes_original() {
        let builder = make_builder();

        let in_field: Vec<FieldT<P>> = (0..4)
            .map(|i| {
                FieldT::from_witness(
                    builder.clone(),
                    fr_from_u128(convert_bytes_to_u128(&PLAINTEXT[i * 16..])),
                )
            })
            .collect();

        let key_field = FieldT::from_witness(
            builder.clone(),
            fr_from_u128(convert_bytes_to_u128(&KEY)),
        );
        let iv_field = FieldT::from_witness(
            builder.clone(),
            fr_from_u128(convert_bytes_to_u128(&IV)),
        );

        let expected = expected_blocks();
        let result = encrypt_buffer_cbc(&in_field, &iv_field, &key_field);

        for i in 0..4 {
            assert_eq!(result[i].get_value(), expected[i]);
        }

        UltraCircuitChecker::check(&mut builder.borrow_mut())
            .expect("circuit check failed");
    }

    // ════════════════════════════════════════════════════════════════════
    //  Tests 10-17: Mixed input block configurations
    // ════════════════════════════════════════════════════════════════════

    #[test]
    fn encrypt_64_bytes_mixed_input_first_witness_rest_constant() {
        test_aes128_mixed_input(false, false, &[true, false, false, false]);
    }

    #[test]
    fn encrypt_64_bytes_mixed_input_alternating_witness_constant() {
        test_aes128_mixed_input(false, false, &[true, false, true, false]);
    }

    #[test]
    fn encrypt_64_bytes_mixed_input_first_constant_rest_witness() {
        test_aes128_mixed_input(false, false, &[false, true, true, true]);
    }

    #[test]
    fn encrypt_64_bytes_mixed_input_key_witness_mixed_blocks() {
        test_aes128_mixed_input(true, false, &[true, false, true, false]);
    }

    #[test]
    fn encrypt_64_bytes_mixed_input_iv_witness_mixed_blocks() {
        test_aes128_mixed_input(false, true, &[false, true, false, true]);
    }

    #[test]
    fn encrypt_64_bytes_mixed_input_key_iv_witness_mixed_blocks() {
        test_aes128_mixed_input(true, true, &[true, false, true, false]);
    }

    #[test]
    fn encrypt_64_bytes_mixed_input_all_witness_blocks() {
        test_aes128_mixed_input(false, false, &[true, true, true, true]);
    }

    #[test]
    fn encrypt_64_bytes_mixed_input_all_constant_blocks() {
        test_aes128_mixed_input(false, false, &[false, false, false, false]);
    }
}
