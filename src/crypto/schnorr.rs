use crate::crypto::hashers::Hasher;
use crate::crypto::generators::GeneratorContext;
use crate::crypto::pedersen_hash;
use crate::ecc::curves::grumpkin::GrumpkinG1Params;
use crate::ecc::groups::affine_element::AffineElement;
use crate::ecc::groups::curve_params::{CurveParams, ScalarField};
use crate::ecc::groups::element::Element;

/// A Schnorr key pair.
pub struct SchnorrKeyPair<C: CurveParams> {
    pub private_key: ScalarField<C>,
    pub public_key: AffineElement<C>,
}

/// A Schnorr signature (BB-custom format).
pub struct SchnorrSignature {
    pub s: [u8; 32],
    pub e: [u8; 32],
}

/// BB's custom Schnorr challenge generation.
///
/// e = H(pedersen_hash(R.x, pk.x, pk.y) || message)
///
/// This is NOT standard Schnorr — it uses Pedersen hash for key/nonce
/// compression to keep the preimage small.
fn schnorr_generate_challenge<H: Hasher>(
    message: &str,
    pubkey: &AffineElement<GrumpkinG1Params>,
    r: &AffineElement<GrumpkinG1Params>,
) -> [u8; 32] {
    // Pedersen hash of [R.x, pubkey.x, pubkey.y] as Grumpkin Fq elements
    let ctx = GeneratorContext::default();
    let compressed = pedersen_hash::hash(&[r.x, pubkey.x, pubkey.y], &ctx);

    // buffer = compressed.to_be_bytes() || message.as_bytes()
    let compressed_bytes = compressed.to_be_bytes();
    let mut buffer = Vec::with_capacity(32 + message.len());
    buffer.extend_from_slice(&compressed_bytes);
    buffer.extend_from_slice(message.as_bytes());

    // e = H::hash(buffer)
    let hash_result = H::hash(&buffer);
    let mut e = [0u8; 32];
    e.copy_from_slice(&hash_result[..32]);
    e
}

/// Construct a Schnorr signature on Grumpkin.
///
/// Matches C++ `schnorr_construct_signature<Hash, Fq, Fr, G1>`.
pub fn schnorr_construct_signature<H: Hasher>(
    message: &str,
    account: &SchnorrKeyPair<GrumpkinG1Params>,
) -> SchnorrSignature {
    // k = random nonce
    let k = ScalarField::<GrumpkinG1Params>::random_element();

    // R = G * k
    let generator = Element::<GrumpkinG1Params>::one();
    let r_point = generator.mul_without_endomorphism(&k).to_affine();

    // Generate challenge e (raw hash bytes)
    let e_raw = schnorr_generate_challenge::<H>(message, &account.public_key, &r_point);

    // Interpret e as field element (with modular reduction)
    let e_fr = ScalarField::<GrumpkinG1Params>::from_be_bytes(&e_raw);

    // s = k - private_key * e
    let s_fr = k - (account.private_key * e_fr);

    SchnorrSignature {
        s: s_fr.to_be_bytes(),
        e: e_raw, // Store raw hash bytes, NOT reduced field element
    }
}

/// Verify a Schnorr signature on Grumpkin.
///
/// Matches C++ `schnorr_verify_signature<Hash, Fq, Fr, G1>`.
pub fn schnorr_verify_signature<H: Hasher>(
    message: &str,
    public_key: &AffineElement<GrumpkinG1Params>,
    sig: &SchnorrSignature,
) -> bool {
    // Validate public key
    if public_key.is_point_at_infinity() || !public_key.on_curve() {
        return false;
    }

    // Deserialize s and e with modular reduction
    let e_fr = ScalarField::<GrumpkinG1Params>::from_be_bytes(&sig.e);
    let s_fr = ScalarField::<GrumpkinG1Params>::from_be_bytes(&sig.s);

    // Reject if either is zero
    if e_fr.is_zero() || s_fr.is_zero() {
        return false;
    }

    // R = pubkey * e + G * s
    let generator = Element::<GrumpkinG1Params>::one();
    let pubkey_proj = Element::<GrumpkinG1Params>::from_affine(public_key);
    let r_point = pubkey_proj.double_scalar_mul(&e_fr, &generator, &s_fr);

    // Reject if R is infinity
    if r_point.is_point_at_infinity() {
        return false;
    }

    let r_affine = r_point.to_affine();

    // Regenerate challenge and compare raw bytes
    let target_e = schnorr_generate_challenge::<H>(message, public_key, &r_affine);

    // Byte comparison (not field comparison!)
    sig.e == target_e
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crypto::hashers::Blake2sHasher;

    fn make_schnorr_keypair() -> SchnorrKeyPair<GrumpkinG1Params> {
        let private_key = ScalarField::<GrumpkinG1Params>::random_element();
        let generator = Element::<GrumpkinG1Params>::one();
        let public_key = generator
            .mul_without_endomorphism(&private_key)
            .to_affine();
        SchnorrKeyPair {
            private_key,
            public_key,
        }
    }

    #[test]
    fn schnorr_sign_verify_grumpkin() {
        let account = make_schnorr_keypair();
        let message = "test Schnorr message on Grumpkin";
        let sig = schnorr_construct_signature::<Blake2sHasher>(message, &account);
        let valid =
            schnorr_verify_signature::<Blake2sHasher>(message, &account.public_key, &sig);
        assert!(valid, "Schnorr signature should verify on Grumpkin");
    }

    #[test]
    fn schnorr_verify_rejects_tampered_message() {
        let account = make_schnorr_keypair();
        let message = "original Schnorr message";
        let sig = schnorr_construct_signature::<Blake2sHasher>(message, &account);
        let valid = schnorr_verify_signature::<Blake2sHasher>(
            "tampered Schnorr message",
            &account.public_key,
            &sig,
        );
        assert!(!valid, "Schnorr should reject tampered message");
    }

    #[test]
    fn schnorr_verify_rejects_wrong_pubkey() {
        let account = make_schnorr_keypair();
        let wrong_account = make_schnorr_keypair();
        let message = "Schnorr wrong key test";
        let sig = schnorr_construct_signature::<Blake2sHasher>(message, &account);
        let valid = schnorr_verify_signature::<Blake2sHasher>(
            message,
            &wrong_account.public_key,
            &sig,
        );
        assert!(!valid, "Schnorr should reject wrong public key");
    }
}
