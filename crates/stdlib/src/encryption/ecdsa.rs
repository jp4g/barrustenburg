//! In-circuit ECDSA signature verification.
//!
//! Port of C++ `barretenberg/stdlib/encryption/ecdsa/ecdsa.{hpp,cpp}`.
//!
//! Verifies ECDSA signatures over non-native curves (e.g. secp256k1, secp256r1)
//! inside a BN254 arithmetic circuit using BigFieldT for non-native field
//! arithmetic and BigGroupT for non-native elliptic curve operations.

use bbrs_ecc::curves::bn254::Bn254FrParams;
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_ecc::groups::curve_params::CurveParams;
use bbrs_numeric::{U256, U256Ext};

use crate::hash::sha256::sha256_block;
use crate::primitives::bigfield::BigFieldT;
use crate::primitives::biggroup::BigGroupT;
use crate::primitives::byte_array::ByteArrayT;
use crate::primitives::field::FieldT;

type P = Bn254FrParams;
type Fr = Field<P>;

/// Non-native scalar field element (signature scalars live here).
type ScalarFieldT<C> = BigFieldT<P, <C as CurveParams>::ScalarFieldParams>;

/// Non-native base field element (point coordinates live here).
type BaseFieldT<C> = BigFieldT<P, <C as CurveParams>::BaseFieldParams>;

// ════════════════════════════════════════════════════════════════════════
//  ECDSA signature (circuit type)
// ════════════════════════════════════════════════════════════════════════

/// An ECDSA signature in-circuit, with r, s as 32-byte arrays and v as 1 byte.
pub struct EcdsaSignatureT {
    pub r: ByteArrayT<P>,
    pub s: ByteArrayT<P>,
    pub v: ByteArrayT<P>,
}

// ════════════════════════════════════════════════════════════════════════
//  Full SHA-256 hash (padding + compression)
// ════════════════════════════════════════════════════════════════════════

/// SHA-256 initial hash values.
const SHA256_INIT: [u32; 8] = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
];

/// Compute the full SHA-256 hash of a byte array in-circuit.
///
/// Applies standard PKCS padding (append 0x80, zeros, 64-bit length) and
/// iterates the SHA-256 compression function over all 512-bit blocks.
/// Returns a 32-byte ByteArrayT containing the hash output.
fn sha256_full(input: &ByteArrayT<P>) -> ByteArrayT<P> {
    let ctx = input.get_context().clone().expect("sha256_full: need circuit context");
    let msg_bytes = input.bytes();
    let msg_len = msg_bytes.len();
    let msg_len_bits = (msg_len as u64) * 8;

    // Padding: append 0x80, then zeros, then 8-byte big-endian length
    // Total padded length must be a multiple of 64 bytes
    let padded_len = {
        let len_with_padding = msg_len + 1 + 8; // 1 for 0x80, 8 for length
        ((len_with_padding + 63) / 64) * 64
    };
    let num_blocks = padded_len / 64;

    // Build padded message as FieldT values
    let mut padded: Vec<FieldT<P>> = Vec::with_capacity(padded_len);
    // Copy message bytes
    for byte in msg_bytes {
        padded.push(byte.clone());
    }
    // Append 0x80
    padded.push(FieldT::constant_with_context(ctx.clone(), Fr::from(0x80u64)));
    // Append zeros until 8 bytes before end
    let zeros_needed = padded_len - msg_len - 1 - 8;
    for _ in 0..zeros_needed {
        padded.push(FieldT::constant_with_context(ctx.clone(), Fr::zero()));
    }
    // Append 8-byte big-endian message length in bits
    for i in (0..8).rev() {
        let byte_val = ((msg_len_bits >> (i * 8)) & 0xFF) as u64;
        padded.push(FieldT::constant_with_context(ctx.clone(), Fr::from(byte_val)));
    }

    assert_eq!(padded.len(), padded_len);

    // Process blocks
    let mut h: [FieldT<P>; 8] = std::array::from_fn(|i| {
        FieldT::constant_with_context(ctx.clone(), Fr::from(SHA256_INIT[i] as u64))
    });

    for block_idx in 0..num_blocks {
        let block_start = block_idx * 64;

        // Convert 64 bytes into 16 × 32-bit words (big-endian)
        let block: [FieldT<P>; 16] = std::array::from_fn(|word_idx| {
            let byte_offset = block_start + word_idx * 4;
            // word = byte[0]*2^24 + byte[1]*2^16 + byte[2]*2^8 + byte[3]
            let b0 = &padded[byte_offset];
            let b1 = &padded[byte_offset + 1];
            let b2 = &padded[byte_offset + 2];
            let b3 = &padded[byte_offset + 3];

            let s24 = FieldT::from_field(Fr::from(1u64 << 24));
            let s16 = FieldT::from_field(Fr::from(1u64 << 16));
            let s8 = FieldT::from_field(Fr::from(1u64 << 8));

            let word = &(b0 * &s24) + &(&(b1 * &s16) + &(&(b2 * &s8) + b3));
            word
        });

        h = sha256_block(&h, &block);
    }

    // Convert 8 × 32-bit words back to 32 bytes (big-endian)
    let mut result_bytes: Vec<FieldT<P>> = Vec::with_capacity(32);
    for word in &h {
        let word_bytes = ByteArrayT::from_field(word, 4, None);
        for byte in word_bytes.bytes() {
            result_bytes.push(byte.clone());
        }
    }

    ByteArrayT::from_values(Some(ctx), result_bytes)
}

// ════════════════════════════════════════════════════════════════════════
//  Byte array → BigFieldT conversion
// ════════════════════════════════════════════════════════════════════════

/// Number of bits per BigFieldT limb.
const NUM_LIMB_BITS: u64 = 68;

/// Create a BigFieldT from a 32-byte ByteArrayT, constraining them to
/// represent the same 256-bit integer.
///
/// The approach:
/// 1. Extract native byte values and compute the U256
/// 2. Create a BigFieldT witness from that U256 (with limb range constraints)
/// 3. Constrain equality: split into two halves (lo 136 bits, hi 120 bits)
///    and assert each matches between the byte representation and the
///    BigFieldT limb representation
fn bigfield_from_bytes<T: FieldParams>(byte_array: &ByteArrayT<P>) -> BigFieldT<P, T> {
    assert_eq!(byte_array.size(), 32, "bigfield_from_bytes: expected 32 bytes");

    let ctx = byte_array.get_context().clone().expect("bigfield_from_bytes: need context");
    let bytes = byte_array.bytes();

    // Compute native U256 value from big-endian bytes
    let native_bytes = byte_array.get_value();
    let mut u256_val = U256::ZERO;
    for &b in &native_bytes {
        u256_val = u256_val.wrapping_shl_vartime(8);
        u256_val = u256_val.wrapping_add(&U256::from(b as u64));
    }

    // Create BigFieldT witness from U256 (this adds limb range constraints)
    let bigfield = BigFieldT::<P, T>::from_witness(ctx.clone(), u256_val);

    // Constrain: bytes represent the same integer as the bigfield limbs.
    //
    // Split at bit 136 (= 2 * NUM_LIMB_BITS = 17 bytes from LSB):
    //   lo = bytes[15..32] interpreted as bits 0-135  (136 bits, fits in native field)
    //   hi = bytes[0..15]  interpreted as bits 136-255 (120 bits, fits in native field)
    //
    // From bigfield:
    //   lo_limbs = limb[0] + limb[1] * 2^68
    //   hi_limbs = limb[2] + limb[3] * 2^68

    // Accumulate lo from bytes: byte[15]*2^128 + byte[16]*2^120 + ... + byte[31]*2^0
    let mut lo_terms: Vec<FieldT<P>> = Vec::with_capacity(17);
    for i in 15..32 {
        let bit_pos = (31 - i) * 8; // byte[31] → bit 0, byte[15] → bit 128
        let scale = Fr::from_limbs(*U256::from(1u64).wrapping_shl_vartime(bit_pos as u32).as_words());
        lo_terms.push(&bytes[i] * &FieldT::from_field(scale));
    }
    let lo_from_bytes = FieldT::accumulate(&lo_terms);

    // Accumulate hi from bytes: byte[0]*2^112 + byte[1]*2^104 + ... + byte[14]*2^0
    let mut hi_terms: Vec<FieldT<P>> = Vec::with_capacity(15);
    for i in 0..15 {
        let bit_pos = (14 - i) * 8; // byte[14] → bit 0, byte[0] → bit 112
        let scale = Fr::from_limbs(*U256::from(1u64).wrapping_shl_vartime(bit_pos as u32).as_words());
        hi_terms.push(&bytes[i] * &FieldT::from_field(scale));
    }
    let hi_from_bytes = FieldT::accumulate(&hi_terms);

    // Reconstruct lo and hi from BigFieldT limbs
    let shift_68 = Fr::from_limbs(*U256::from(1u64).wrapping_shl_vartime(NUM_LIMB_BITS as u32).as_words());
    let lo_from_limbs = &bigfield.binary_basis_limbs[0].element
        + &(&bigfield.binary_basis_limbs[1].element * &FieldT::from_field(shift_68));
    let hi_from_limbs = &bigfield.binary_basis_limbs[2].element
        + &(&bigfield.binary_basis_limbs[3].element * &FieldT::from_field(shift_68));

    // Assert equality (both sides < 2^136 and native field > 2^253, so mod-field comparison is exact)
    lo_from_bytes.assert_equal(&lo_from_limbs, "ecdsa: byte→bigfield lo mismatch");
    hi_from_bytes.assert_equal(&hi_from_limbs, "ecdsa: byte→bigfield hi mismatch");

    bigfield
}

// ════════════════════════════════════════════════════════════════════════
//  ECDSA signature verification
// ════════════════════════════════════════════════════════════════════════

/// Verify an ECDSA signature in-circuit.
///
/// Port of C++ `ecdsa::verify_signature<Builder, Curve, Fq, Fr, G1>`.
///
/// The verification algorithm:
/// 1. Hash the message with SHA-256 to get `z` (32 bytes)
/// 2. Convert `z`, `r`, `s` to non-native scalar field elements
/// 3. Compute `u1 = z / s` and `u2 = r / s`
/// 4. Compute `R = u1 * G + u2 * PK` using batch_mul
/// 5. Assert `R.x == r` (comparing base field coordinate against scalar field value)
///
/// For secp256k1 and secp256r1, the scalar field modulus is smaller than the
/// base field modulus, so `r` (which is `R.x mod n`) can be directly compared
/// as a base field element.
pub fn verify_signature<C: CurveParams>(
    message: &ByteArrayT<P>,
    public_key: &BigGroupT<P, C>,
    sig: &EcdsaSignatureT,
) {
    let ctx = message.get_context().clone().expect("verify_signature: need context");

    // Step 1: Hash message with SHA-256
    let hashed_message = sha256_full(message);

    // Step 2: Convert hash and signature components to non-native scalar field elements
    let z: ScalarFieldT<C> = bigfield_from_bytes(&hashed_message);
    let r: ScalarFieldT<C> = bigfield_from_bytes(&sig.r);
    let s: ScalarFieldT<C> = bigfield_from_bytes(&sig.s);

    // Step 3: Compute u1 = z / s, u2 = r / s (in the non-native scalar field)
    let u1 = z.div(&s);
    let u2 = r.div(&s);

    // Step 4: Compute R = u1 * G + u2 * PK using batch_mul
    // For ECDSA-specific curves (secp256k1, secp256r1), we use with_edgecases=false
    // because the scalars u1, u2 are essentially random, making edge cases negligible.
    let generator = BigGroupT::<P, C>::one(Some(ctx.clone()));
    let result = BigGroupT::<P, C>::batch_mul(
        &[public_key.clone(), generator],
        &[u2, u1],
        0,     // max_num_bits = 0 (use full scalar field width)
        false, // with_edgecases = false for ECDSA
    );

    // Step 5: Assert R.x == r
    // r was parsed as a scalar field element, but R.x is a base field element.
    // Since r < Fr_mod < Fq_mod for secp256k1/r1, we can compare in the base field.
    // Create an Fq element from the r bytes and assert it equals result.x.
    let r_as_base_field: BaseFieldT<C> = bigfield_from_bytes(&sig.r);
    result.x().assert_equal(&r_as_base_field, "ecdsa: R.x != r");
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
    use bbrs_crypto::ecdsa::{
        ecdsa_construct_signature, ecdsa_verify_signature, EcdsaKeyPair, EcdsaSignature,
    };
    use bbrs_crypto::hashers::Sha256Hasher;
    use bbrs_ecc::curves::secp256k1::Secp256k1G1Params;
    use bbrs_ecc::curves::secp256r1::Secp256r1G1Params;
    use bbrs_ecc::groups::curve_params::ScalarField;
    use bbrs_ecc::groups::element::Element;
    use crate::primitives::witness::BuilderRef;

    fn make_builder() -> BuilderRef<P> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    fn check_circuit(builder: &BuilderRef<P>) -> Result<(), String> {
        UltraCircuitChecker::check(&mut builder.borrow_mut())
    }

    /// Helper: create a random ECDSA key pair.
    fn make_keypair<C: CurveParams>() -> EcdsaKeyPair<C> {
        let private_key = ScalarField::<C>::random_element();
        let generator = Element::<C>::one();
        let public_key = generator
            .mul_without_endomorphism(&private_key)
            .to_affine();
        EcdsaKeyPair {
            private_key,
            public_key,
        }
    }

    /// Helper: create circuit witnesses for an ECDSA verification.
    fn setup_ecdsa_circuit<C: CurveParams>(
        message: &str,
        account: &EcdsaKeyPair<C>,
        sig: &EcdsaSignature,
    ) -> (BuilderRef<P>, ByteArrayT<P>, BigGroupT<P, C>, EcdsaSignatureT) {
        let builder = make_builder();

        // Create message witness
        let msg_bytes = ByteArrayT::from_bytes(builder.clone(), message.as_bytes());

        // Create public key witness
        let pk = BigGroupT::<P, C>::from_witness(builder.clone(), &account.public_key);

        // Create signature witnesses
        let sig_t = EcdsaSignatureT {
            r: ByteArrayT::from_bytes(builder.clone(), &sig.r),
            s: ByteArrayT::from_bytes(builder.clone(), &sig.s),
            v: ByteArrayT::from_bytes(builder.clone(), &[sig.v]),
        };

        (builder, msg_bytes, pk, sig_t)
    }

    // ── Test: SHA-256 full hash matches native ─────────────────────────

    #[test]
    fn test_sha256_full_matches_native() {
        let builder = make_builder();
        let message = b"test ECDSA message";

        let msg_bytes = ByteArrayT::from_bytes(builder.clone(), message);
        let circuit_hash = sha256_full(&msg_bytes);

        let native_hash = bbrs_crypto::sha256::sha256(message);
        let circuit_hash_bytes = circuit_hash.get_value();

        // Compare each byte
        for i in 0..32 {
            assert_eq!(
                circuit_hash_bytes[i], native_hash[i],
                "SHA-256 mismatch at byte {i}"
            );
        }

        check_circuit(&builder).expect("SHA-256 full hash circuit check failed");
    }

    // ── Test: bigfield_from_bytes roundtrip ─────────────────────────────

    #[test]
    fn test_bigfield_from_bytes() {
        use bbrs_ecc::curves::secp256k1::Secp256k1FrParams;

        let builder = make_builder();

        // Create a random 32-byte value
        let mut test_bytes = [0u8; 32];
        let random_val = Field::<Secp256k1FrParams>::random_element().from_montgomery_form();
        let val_bytes = {
            let mut bytes = [0u8; 32];
            for i in 0..4 {
                let limb = random_val.data[3 - i];
                bytes[i * 8..(i + 1) * 8].copy_from_slice(&limb.to_be_bytes());
            }
            bytes
        };
        test_bytes = val_bytes;

        let byte_array = ByteArrayT::from_bytes(builder.clone(), &test_bytes);
        let _bigfield: BigFieldT<P, Secp256k1FrParams> = bigfield_from_bytes(&byte_array);

        check_circuit(&builder).expect("bigfield_from_bytes circuit check failed");
    }

    // ── Test: ECDSA verify signature secp256k1 ─────────────────────────

    #[test]
    fn test_ecdsa_verify_signature_secp256k1() {
        let account = make_keypair::<Secp256k1G1Params>();
        let message = "test ECDSA message for secp256k1";

        // Sign natively
        let sig = ecdsa_construct_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account);

        // Verify natively first
        assert!(
            ecdsa_verify_signature::<Sha256Hasher, Secp256k1G1Params>(
                message,
                &account.public_key,
                &sig
            ),
            "Native ECDSA verification failed"
        );

        // Verify in circuit
        let (builder, msg, pk, sig_t) =
            setup_ecdsa_circuit::<Secp256k1G1Params>(message, &account, &sig);
        verify_signature::<Secp256k1G1Params>(&msg, &pk, &sig_t);

        check_circuit(&builder).expect("ECDSA secp256k1 circuit check failed");
    }

    // ── Test: ECDSA verify signature secp256r1 ─────────────────────────

    #[test]
    fn test_ecdsa_verify_signature_secp256r1() {
        let account = make_keypair::<Secp256r1G1Params>();
        let message = "test ECDSA message for secp256r1";

        let sig = ecdsa_construct_signature::<Sha256Hasher, Secp256r1G1Params>(message, &account);

        assert!(
            ecdsa_verify_signature::<Sha256Hasher, Secp256r1G1Params>(
                message,
                &account.public_key,
                &sig
            ),
            "Native ECDSA verification failed"
        );

        let (builder, msg, pk, sig_t) =
            setup_ecdsa_circuit::<Secp256r1G1Params>(message, &account, &sig);
        verify_signature::<Secp256r1G1Params>(&msg, &pk, &sig_t);

        check_circuit(&builder).expect("ECDSA secp256r1 circuit check failed");
    }

    // ── Test: ECDSA fails with wrong message ───────────────────────────

    #[test]
    fn test_ecdsa_fails_with_wrong_message() {
        let account = make_keypair::<Secp256k1G1Params>();
        let message = "original message";
        let wrong_message = "tampered message";

        let sig =
            ecdsa_construct_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account);

        // Verify with wrong message in circuit — should fail circuit check
        let (builder, wrong_msg, pk, sig_t) =
            setup_ecdsa_circuit::<Secp256k1G1Params>(wrong_message, &account, &sig);
        verify_signature::<Secp256k1G1Params>(&wrong_msg, &pk, &sig_t);

        assert!(
            check_circuit(&builder).is_err(),
            "ECDSA should fail with tampered message"
        );
    }

    // ── Test: ECDSA fails with wrong public key ────────────────────────

    #[test]
    fn test_ecdsa_fails_with_wrong_pubkey() {
        let account = make_keypair::<Secp256k1G1Params>();
        let wrong_account = make_keypair::<Secp256k1G1Params>();
        let message = "test message";

        let sig =
            ecdsa_construct_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account);

        // Verify with wrong public key
        let builder = make_builder();
        let msg = ByteArrayT::from_bytes(builder.clone(), message.as_bytes());
        let wrong_pk = BigGroupT::<P, Secp256k1G1Params>::from_witness(
            builder.clone(),
            &wrong_account.public_key,
        );
        let sig_t = EcdsaSignatureT {
            r: ByteArrayT::from_bytes(builder.clone(), &sig.r),
            s: ByteArrayT::from_bytes(builder.clone(), &sig.s),
            v: ByteArrayT::from_bytes(builder.clone(), &[sig.v]),
        };

        verify_signature::<Secp256k1G1Params>(&msg, &wrong_pk, &sig_t);

        assert!(
            check_circuit(&builder).is_err(),
            "ECDSA should fail with wrong public key"
        );
    }

    // ── Test: ECDSA fails with tampered signature r ────────────────────

    #[test]
    fn test_ecdsa_fails_with_tampered_r() {
        let account = make_keypair::<Secp256k1G1Params>();
        let message = "test message for tampered r";

        let mut sig =
            ecdsa_construct_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account);

        // Tamper with r (flip a bit)
        sig.r[0] ^= 0x01;

        let (builder, msg, pk, sig_t) =
            setup_ecdsa_circuit::<Secp256k1G1Params>(message, &account, &sig);
        verify_signature::<Secp256k1G1Params>(&msg, &pk, &sig_t);

        assert!(
            check_circuit(&builder).is_err(),
            "ECDSA should fail with tampered r"
        );
    }

    // ── Test: ECDSA fails with tampered signature s ────────────────────

    #[test]
    fn test_ecdsa_fails_with_tampered_s() {
        let account = make_keypair::<Secp256k1G1Params>();
        let message = "test message for tampered s";

        let mut sig =
            ecdsa_construct_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account);

        // Tamper with s (flip a bit)
        sig.s[0] ^= 0x01;

        let (builder, msg, pk, sig_t) =
            setup_ecdsa_circuit::<Secp256k1G1Params>(message, &account, &sig);
        verify_signature::<Secp256k1G1Params>(&msg, &pk, &sig_t);

        assert!(
            check_circuit(&builder).is_err(),
            "ECDSA should fail with tampered s"
        );
    }

    // ── Test: SHA-256 empty message ────────────────────────────────────

    #[test]
    fn test_sha256_full_empty_message() {
        let builder = make_builder();

        let msg_bytes = ByteArrayT::from_bytes(builder.clone(), b"");
        let circuit_hash = sha256_full(&msg_bytes);

        let native_hash = bbrs_crypto::sha256::sha256(b"");
        let circuit_hash_bytes = circuit_hash.get_value();

        for i in 0..32 {
            assert_eq!(
                circuit_hash_bytes[i], native_hash[i],
                "SHA-256 empty message mismatch at byte {i}"
            );
        }

        check_circuit(&builder).expect("SHA-256 empty message circuit check failed");
    }

    // ── Test: SHA-256 multi-block message ──────────────────────────────

    #[test]
    fn test_sha256_full_multi_block() {
        let builder = make_builder();

        // 56+ bytes requires 2 blocks
        let message = b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
        let msg_bytes = ByteArrayT::from_bytes(builder.clone(), message);
        let circuit_hash = sha256_full(&msg_bytes);

        let native_hash = bbrs_crypto::sha256::sha256(message);
        let circuit_hash_bytes = circuit_hash.get_value();

        for i in 0..32 {
            assert_eq!(
                circuit_hash_bytes[i], native_hash[i],
                "SHA-256 multi-block mismatch at byte {i}"
            );
        }

        check_circuit(&builder).expect("SHA-256 multi-block circuit check failed");
    }

    // ── Test: SHA-256 exact block boundary ─────────────────────────────

    #[test]
    fn test_sha256_full_55_bytes() {
        let builder = make_builder();

        // 55 bytes: padding fits in one block (55 + 1 + 8 = 64)
        let message = &[0x61u8; 55]; // 55 'a's
        let msg_bytes = ByteArrayT::from_bytes(builder.clone(), message);
        let circuit_hash = sha256_full(&msg_bytes);

        let native_hash = bbrs_crypto::sha256::sha256(message);
        let circuit_hash_bytes = circuit_hash.get_value();

        for i in 0..32 {
            assert_eq!(
                circuit_hash_bytes[i], native_hash[i],
                "SHA-256 55-byte mismatch at byte {i}"
            );
        }

        check_circuit(&builder).expect("SHA-256 55-byte circuit check failed");
    }

    #[test]
    fn test_sha256_full_56_bytes() {
        let builder = make_builder();

        // 56 bytes: padding doesn't fit, needs 2 blocks (56 + 1 + 8 = 65 > 64)
        let message = &[0x61u8; 56]; // 56 'a's
        let msg_bytes = ByteArrayT::from_bytes(builder.clone(), message);
        let circuit_hash = sha256_full(&msg_bytes);

        let native_hash = bbrs_crypto::sha256::sha256(message);
        let circuit_hash_bytes = circuit_hash.get_value();

        for i in 0..32 {
            assert_eq!(
                circuit_hash_bytes[i], native_hash[i],
                "SHA-256 56-byte mismatch at byte {i}"
            );
        }

        check_circuit(&builder).expect("SHA-256 56-byte circuit check failed");
    }

    // ── Test: ECDSA secp256k1 with short message ───────────────────────

    #[test]
    fn test_ecdsa_verify_short_message_secp256k1() {
        let account = make_keypair::<Secp256k1G1Params>();
        let message = "hi";

        let sig = ecdsa_construct_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account);

        let (builder, msg, pk, sig_t) =
            setup_ecdsa_circuit::<Secp256k1G1Params>(message, &account, &sig);
        verify_signature::<Secp256k1G1Params>(&msg, &pk, &sig_t);

        check_circuit(&builder).expect("ECDSA short message circuit check failed");
    }

    // ── Test: ECDSA secp256k1 with long message ────────────────────────

    #[test]
    fn test_ecdsa_verify_long_message_secp256k1() {
        let account = make_keypair::<Secp256k1G1Params>();
        // Message longer than one SHA-256 block
        let message = "This is a longer test message that exceeds a single SHA-256 block boundary for ECDSA verification testing purposes.";

        let sig = ecdsa_construct_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account);

        let (builder, msg, pk, sig_t) =
            setup_ecdsa_circuit::<Secp256k1G1Params>(message, &account, &sig);
        verify_signature::<Secp256k1G1Params>(&msg, &pk, &sig_t);

        check_circuit(&builder).expect("ECDSA long message circuit check failed");
    }
}
