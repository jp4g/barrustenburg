use crate::crypto::hashers::Hasher;
use crate::crypto::hmac::get_unbiased_field_from_hmac;
use crate::ecc::fields::field_params::FieldParams;
use crate::ecc::groups::affine_element::AffineElement;
use crate::ecc::groups::curve_params::{BaseField, CurveParams, ScalarField};
use crate::ecc::groups::element::Element;

/// An ECDSA key pair.
pub struct EcdsaKeyPair<C: CurveParams> {
    pub private_key: ScalarField<C>,
    pub public_key: AffineElement<C>,
}

/// An ECDSA signature with recovery ID.
pub struct EcdsaSignature {
    pub r: [u8; 32],
    pub s: [u8; 32],
    pub v: u8,
}

/// Construct an ECDSA signature.
///
/// Matches C++ `ecdsa_construct_signature<Hash, Fq, Fr, G1>`.
pub fn ecdsa_construct_signature<H: Hasher, C: CurveParams>(
    message: &str,
    account: &EcdsaKeyPair<C>,
) -> EcdsaSignature {
    // Serialize private key to bytes
    let pk_bytes = account.private_key.to_be_bytes();

    // Derive nonce k via HMAC
    let k: ScalarField<C> =
        get_unbiased_field_from_hmac::<H, C::ScalarFieldParams>(message.as_bytes(), &pk_bytes);

    // R = G * k
    let generator = Element::<C>::one();
    let r_point = generator.mul_without_endomorphism(&k).to_affine();

    // sig.r = R.x serialized as base field bytes
    let r_fq: BaseField<C> = r_point.x;
    let sig_r = r_fq.to_be_bytes();

    // Hash message -> z
    let z_bytes = H::hash(message.as_bytes());
    let z_fr = ScalarField::<C>::from_be_bytes(z_bytes[..32].try_into().unwrap());

    // Convert r to scalar field: r_fr = Fr::from_be_bytes(sig.r)
    let r_fr = ScalarField::<C>::from_be_bytes(&sig_r);

    // s = (z + r * priv) / k
    let s_fr = (z_fr + r_fr * account.private_key) * k.invert();

    // Check if s is "low" (s < modulus/2)
    let s_reduced = s_fr.from_montgomery_form();
    let modulus = C::ScalarFieldParams::MODULUS;
    // modulus / 2 (shift right by 1)
    let half_mod = [
        (modulus[0] >> 1) | (modulus[1] << 63),
        (modulus[1] >> 1) | (modulus[2] << 63),
        (modulus[2] >> 1) | (modulus[3] << 63),
        modulus[3] >> 1,
    ];
    let is_s_low = cmp_u256(&s_reduced.data, &half_mod) == std::cmp::Ordering::Less;

    let final_s = if is_s_low { s_fr } else { s_fr.negate() };
    let sig_s = final_s.to_be_bytes();

    // Recovery ID computation
    // is_r_finite: check if R.x (as Fq integer) < Fr::modulus
    // In other words, the Fq representation and Fr representation of R.x are the same integer
    let r_fq_reduced = r_fq.from_montgomery_form();
    let r_fr_reduced = r_fr.from_montgomery_form();
    let is_r_finite = r_fq_reduced.data == r_fr_reduced.data;

    // y_parity: LSB of R.y
    let r_y_reduced = r_point.y.from_montgomery_form();
    let y_parity = r_y_reduced.data[0] & 1 == 1;

    let recovery_bit = y_parity ^ is_s_low;

    let v = 27u8 + recovery_bit as u8 + 2 * (!is_r_finite) as u8;

    EcdsaSignature {
        r: sig_r,
        s: sig_s,
        v,
    }
}

/// Verify an ECDSA signature.
///
/// Matches C++ `ecdsa_verify_signature<Hash, Fq, Fr, G1>`.
pub fn ecdsa_verify_signature<H: Hasher, C: CurveParams>(
    message: &str,
    public_key: &AffineElement<C>,
    signature: &EcdsaSignature,
) -> bool {
    // Validate public key
    if public_key.is_point_at_infinity() || !public_key.on_curve() {
        return false;
    }

    // Parse r and s from bytes
    let r_fr = ScalarField::<C>::from_be_bytes(&signature.r);
    let s_fr = ScalarField::<C>::from_be_bytes(&signature.s);

    // Validate r, s non-zero
    if r_fr.is_zero() || s_fr.is_zero() {
        return false;
    }

    // Validate s < modulus/2
    let s_reduced = s_fr.from_montgomery_form();
    let modulus = C::ScalarFieldParams::MODULUS;
    let half_mod = [
        (modulus[0] >> 1) | (modulus[1] << 63),
        (modulus[1] >> 1) | (modulus[2] << 63),
        (modulus[2] >> 1) | (modulus[3] << 63),
        modulus[3] >> 1,
    ];
    if cmp_u256(&s_reduced.data, &half_mod) != std::cmp::Ordering::Less {
        return false;
    }

    // z = hash(message) as Fr
    let z_bytes = H::hash(message.as_bytes());
    let z = ScalarField::<C>::from_be_bytes(z_bytes[..32].try_into().unwrap());

    // s_inv, u1, u2
    let s_inv = s_fr.invert();
    let u1 = z * s_inv;
    let u2 = r_fr * s_inv;

    // R = u1 * G + u2 * pubkey
    let generator = Element::<C>::one();
    let pubkey_proj = Element::<C>::from_affine(public_key);
    let r_point = generator.double_scalar_mul(&u1, &pubkey_proj, &u2);

    if r_point.is_point_at_infinity() {
        return false;
    }

    // Check R.x mod Fr == r
    let r_affine = r_point.to_affine();
    let r_x_bytes = r_affine.x.to_be_bytes();
    let r_x_fr = ScalarField::<C>::from_be_bytes(&r_x_bytes);

    r_x_fr == r_fr
}

/// Recover the public key from an ECDSA signature.
///
/// Matches C++ `ecdsa_recover_public_key<Hash, Fq, Fr, G1>`.
pub fn ecdsa_recover_public_key<H: Hasher, C: CurveParams>(
    message: &str,
    sig: &EcdsaSignature,
) -> Option<AffineElement<C>> {
    // Parse r, s
    let r_fr = ScalarField::<C>::from_be_bytes(&sig.r);
    let s_fr = ScalarField::<C>::from_be_bytes(&sig.s);

    if r_fr.is_zero() || s_fr.is_zero() {
        return None;
    }

    // Validate s < modulus/2
    let s_reduced = s_fr.from_montgomery_form();
    let modulus = C::ScalarFieldParams::MODULUS;
    let half_mod = [
        (modulus[0] >> 1) | (modulus[1] << 63),
        (modulus[1] >> 1) | (modulus[2] << 63),
        (modulus[2] >> 1) | (modulus[3] << 63),
        modulus[3] >> 1,
    ];
    if cmp_u256(&s_reduced.data, &half_mod) != std::cmp::Ordering::Less {
        return None;
    }

    // Determine y-parity from v
    let recovery_id = sig.v - 27;
    let is_r_finite = (recovery_id & 2) == 0;

    // Recover R from r (x-coordinate)
    // If !is_r_finite, we need to add Fr::modulus to get the actual x in Fq
    let r_x: BaseField<C> = if is_r_finite {
        BaseField::<C>::from_be_bytes(&sig.r)
    } else {
        // R.x = r + Fr::modulus (as Fq element)
        let r_int = r_fr.from_montgomery_form();
        let fr_mod = C::ScalarFieldParams::MODULUS;
        // Add Fr::modulus to r and interpret as Fq
        let mut carry = 0u64;
        let (d0, c) = add_with_carry(r_int.data[0], fr_mod[0], carry);
        carry = c;
        let (d1, c) = add_with_carry(r_int.data[1], fr_mod[1], carry);
        carry = c;
        let (d2, c) = add_with_carry(r_int.data[2], fr_mod[2], carry);
        carry = c;
        let (d3, _) = add_with_carry(r_int.data[3], fr_mod[3], carry);
        BaseField::<C>::from_limbs([d0, d1, d2, d3])
    };

    // Determine which y to use based on recovery_bit
    // recovery_bit = y_parity XOR is_s_low
    // From verify perspective: recovery_id & 1 gives the recovery_bit
    let recovery_bit = (recovery_id & 1) != 0;

    // We need to figure out the sign bit for y.
    // recovery_bit = y_parity XOR is_s_low
    // is_s_low is true (we validated s < modulus/2)
    // so y_parity = recovery_bit XOR true = !recovery_bit
    let y_parity = !recovery_bit;

    let r_point = AffineElement::<C>::from_x_coordinate(r_x, y_parity)?;

    // Compute public key: pk = r_inv * (s * R - z * G)
    let z_bytes = H::hash(message.as_bytes());
    let z = ScalarField::<C>::from_be_bytes(z_bytes[..32].try_into().unwrap());

    let r_inv = r_fr.invert();
    let u1 = (z * r_inv).negate(); // -z * r_inv
    let u2 = s_fr * r_inv;

    let generator = Element::<C>::one();
    let r_proj = Element::<C>::from_affine(&r_point);
    let pk = r_proj.double_scalar_mul(&u2, &generator, &u1);

    if pk.is_point_at_infinity() {
        return None;
    }

    Some(pk.to_affine())
}

/// Compare two u256 values represented as [u64; 4] limbs (little-endian).
fn cmp_u256(a: &[u64; 4], b: &[u64; 4]) -> std::cmp::Ordering {
    for i in (0..4).rev() {
        match a[i].cmp(&b[i]) {
            std::cmp::Ordering::Equal => continue,
            ord => return ord,
        }
    }
    std::cmp::Ordering::Equal
}

/// Add with carry helper.
fn add_with_carry(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let result = a as u128 + b as u128 + carry as u128;
    (result as u64, (result >> 64) as u64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crypto::hashers::Sha256Hasher;
    use crate::ecc::curves::secp256k1::Secp256k1G1Params;
    use crate::ecc::curves::secp256r1::Secp256r1G1Params;

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

    #[test]
    fn ecdsa_sign_verify_secp256k1() {
        let account = make_keypair::<Secp256k1G1Params>();
        let message = "test ECDSA message for secp256k1";
        let sig = ecdsa_construct_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account);
        let valid =
            ecdsa_verify_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account.public_key, &sig);
        assert!(valid, "ECDSA signature should verify on secp256k1");
    }

    #[test]
    fn ecdsa_sign_verify_secp256r1() {
        let account = make_keypair::<Secp256r1G1Params>();
        let message = "test ECDSA message for secp256r1";
        let sig = ecdsa_construct_signature::<Sha256Hasher, Secp256r1G1Params>(message, &account);
        let valid =
            ecdsa_verify_signature::<Sha256Hasher, Secp256r1G1Params>(message, &account.public_key, &sig);
        assert!(valid, "ECDSA signature should verify on secp256r1");
    }

    #[test]
    fn ecdsa_recover_public_key_secp256k1() {
        let account = make_keypair::<Secp256k1G1Params>();
        let message = "test ECDSA recovery on secp256k1";
        let sig = ecdsa_construct_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account);
        let recovered =
            ecdsa_recover_public_key::<Sha256Hasher, Secp256k1G1Params>(message, &sig);
        assert!(recovered.is_some(), "should recover public key");
        assert_eq!(
            recovered.unwrap(),
            account.public_key,
            "recovered key should match original"
        );
    }

    #[test]
    fn ecdsa_verify_rejects_tampered_message() {
        let account = make_keypair::<Secp256k1G1Params>();
        let message = "original message";
        let sig = ecdsa_construct_signature::<Sha256Hasher, Secp256k1G1Params>(message, &account);
        let valid = ecdsa_verify_signature::<Sha256Hasher, Secp256k1G1Params>(
            "tampered message",
            &account.public_key,
            &sig,
        );
        assert!(!valid, "ECDSA should reject tampered message");
    }
}
