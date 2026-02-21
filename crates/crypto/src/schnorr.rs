use crate::hashers::Hasher;
use crate::generators::GeneratorContext;
use crate::pedersen_hash;
use bbrs_ecc::curves::grumpkin::GrumpkinG1Params;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::{CurveParams, ScalarField};
use bbrs_ecc::groups::element::Element;

type GrumpkinFr = ScalarField<GrumpkinG1Params>;
type GrumpkinAffine = AffineElement<GrumpkinG1Params>;

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
    pubkey: &GrumpkinAffine,
    r: &GrumpkinAffine,
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
    let k = GrumpkinFr::random_element();

    // R = G * k
    let generator = Element::<GrumpkinG1Params>::one();
    let r_point = generator.mul_without_endomorphism(&k).to_affine();

    // Generate challenge e (raw hash bytes)
    let e_raw = schnorr_generate_challenge::<H>(message, &account.public_key, &r_point);

    // Interpret e as field element (with modular reduction)
    let e_fr = GrumpkinFr::from_be_bytes(&e_raw);

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
    public_key: &GrumpkinAffine,
    sig: &SchnorrSignature,
) -> bool {
    // Validate public key
    if public_key.is_point_at_infinity() || !public_key.on_curve() {
        return false;
    }

    // Deserialize s and e with modular reduction
    let e_fr = GrumpkinFr::from_be_bytes(&sig.e);
    let s_fr = GrumpkinFr::from_be_bytes(&sig.s);

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

// ---------------------------------------------------------------------------
// Proof of Possession
// ---------------------------------------------------------------------------

/// A Schnorr proof of knowledge of a secret key corresponding to a given public key.
///
/// Follows the SpeedyMuSig specification (https://eprint.iacr.org/2021/1375.pdf).
/// Prevents rogue-key attacks in multisig protocols.
pub struct SchnorrProofOfPossession {
    /// challenge = e = H_reg(G, pk, pk, R)
    pub challenge: [u8; 32],
    /// response = z = k - e * sk
    pub response: GrumpkinFr,
}

impl SchnorrProofOfPossession {
    /// Create a new proof of possession for a given account.
    pub fn new<H: Hasher>(account: &SchnorrKeyPair<GrumpkinG1Params>) -> Self {
        let k = GrumpkinFr::random_element();
        let r_point = (Element::<GrumpkinG1Params>::one()
            .mul_without_endomorphism(&k))
        .to_affine();

        let challenge = pop_generate_challenge::<H>(&account.public_key, &r_point);
        let challenge_fr = GrumpkinFr::from_be_bytes(&challenge);
        let response = k - challenge_fr * account.private_key;

        Self {
            challenge,
            response,
        }
    }

    /// Verify that this proof is valid for the given public key.
    pub fn verify<H: Hasher>(&self, public_key: &GrumpkinAffine) -> bool {
        if self.response.is_zero() {
            return false;
        }
        if !public_key.on_curve() || public_key.is_point_at_infinity() {
            return false;
        }

        let challenge_fr = GrumpkinFr::from_be_bytes(&self.challenge);

        // R = e*pk + z*G
        let pk_proj = Element::<GrumpkinG1Params>::from_affine(public_key);
        let generator = Element::<GrumpkinG1Params>::one();
        let r_point = pk_proj.double_scalar_mul(&challenge_fr, &generator, &self.response);

        if r_point.is_point_at_infinity() {
            return false;
        }

        let r_affine = r_point.to_affine();
        let challenge_computed = pop_generate_challenge::<H>(public_key, &r_affine);
        self.challenge == challenge_computed
    }
}

/// Generate the Fiat-Shamir challenge for proof of possession.
///
/// e = H_reg("h_reg" || G || pk || pk || R)
///
/// Uses x-first (legacy) byte order for affine point serialization,
/// matching the C++ `write(buf, affine_element)` function.
fn pop_generate_challenge<H: Hasher>(
    public_key: &GrumpkinAffine,
    r: &GrumpkinAffine,
) -> [u8; 32] {
    let domain_sep = b"h_reg";
    let g = GrumpkinAffine::one();

    let mut buf = Vec::with_capacity(5 + 64 * 4);
    buf.extend_from_slice(domain_sep);
    buf.extend_from_slice(&g.to_buffer());
    buf.extend_from_slice(&public_key.to_buffer());
    buf.extend_from_slice(&public_key.to_buffer());
    buf.extend_from_slice(&r.to_buffer());

    let hash_result = H::hash(&buf);
    let mut challenge = [0u8; 32];
    challenge.copy_from_slice(&hash_result[..32]);
    challenge
}

// ---------------------------------------------------------------------------
// SpeedyMuSig Multisig Protocol
// ---------------------------------------------------------------------------

/// A signer's public key bundled with a proof of possession.
pub struct MultiSigPublicKey<HRegNon: Hasher> {
    pub public_key: GrumpkinAffine,
    pub proof_of_possession: SchnorrProofOfPossession,
    _phantom: std::marker::PhantomData<HRegNon>,
}

impl<HRegNon: Hasher> MultiSigPublicKey<HRegNon> {
    pub fn new(account: &SchnorrKeyPair<GrumpkinG1Params>) -> Self {
        Self {
            public_key: account.public_key,
            proof_of_possession: SchnorrProofOfPossession::new::<HRegNon>(account),
            _phantom: std::marker::PhantomData,
        }
    }
}

/// Round 1 private output: nonce scalars.
pub struct RoundOnePrivateOutput {
    pub r: GrumpkinFr,
    pub s: GrumpkinFr,
}

/// Round 1 public output: nonce group elements.
#[derive(Clone, Copy)]
pub struct RoundOnePublicOutput {
    pub r: GrumpkinAffine,
    pub s: GrumpkinAffine,
}

impl PartialEq for RoundOnePublicOutput {
    fn eq(&self, other: &Self) -> bool {
        self.r == other.r && self.s == other.s
    }
}

/// Round 2 public output: a single field element (the signer's signature share).
pub type RoundTwoPublicOutput = GrumpkinFr;

/// SpeedyMuSig 2-round interactive multisignature protocol.
///
/// - `HRegNon`: Hash function for H_reg (PoP) and H_non (nonce challenge)
/// - `HSig`: Hash function for the final Schnorr signature challenge
///
/// HRegNon and HSig must be different types for proper domain separation.
pub struct SchnorrMultisig<HRegNon: Hasher, HSig: Hasher> {
    _reg: std::marker::PhantomData<HRegNon>,
    _sig: std::marker::PhantomData<HSig>,
}

impl<HRegNon: Hasher, HSig: Hasher> SchnorrMultisig<HRegNon, HSig> {
    /// Validate signer public keys and combine into an aggregate key.
    ///
    /// Checks: all on_curve, no infinity, no duplicates, all PoPs valid.
    /// Returns `X_agg = pk_1 + pk_2 + ... + pk_n`.
    pub fn validate_and_combine_signer_pubkeys(
        signer_pubkeys: &[MultiSigPublicKey<HRegNon>],
    ) -> Option<GrumpkinAffine> {
        // Check for duplicate public keys
        for i in 0..signer_pubkeys.len() {
            for j in (i + 1)..signer_pubkeys.len() {
                if signer_pubkeys[i].public_key == signer_pubkeys[j].public_key {
                    return None;
                }
            }
        }

        let mut aggregate = Element::<GrumpkinG1Params>::from_affine(
            &AffineElement::infinity(),
        );

        for spk in signer_pubkeys {
            if !spk.public_key.on_curve() || spk.public_key.is_point_at_infinity() {
                return None;
            }
            if !spk.proof_of_possession.verify::<HRegNon>(&spk.public_key) {
                return None;
            }
            aggregate = aggregate + Element::<GrumpkinG1Params>::from_affine(&spk.public_key);
        }

        let aggregate_affine = aggregate.to_affine();
        if aggregate_affine.is_point_at_infinity() {
            return None;
        }
        Some(aggregate_affine)
    }

    /// Round 1: Generate random nonce keypairs.
    pub fn construct_signature_round_1() -> (RoundOnePublicOutput, RoundOnePrivateOutput) {
        let r_scalar = GrumpkinFr::random_element();
        let r_point = Element::<GrumpkinG1Params>::one()
            .mul_without_endomorphism(&r_scalar)
            .to_affine();

        let s_scalar = GrumpkinFr::random_element();
        let s_point = Element::<GrumpkinG1Params>::one()
            .mul_without_endomorphism(&s_scalar)
            .to_affine();

        (
            RoundOnePublicOutput {
                r: r_point,
                s: s_point,
            },
            RoundOnePrivateOutput {
                r: r_scalar,
                s: s_scalar,
            },
        )
    }

    /// Round 2: Compute the signer's signature share.
    pub fn construct_signature_round_2(
        message: &str,
        signer: &SchnorrKeyPair<GrumpkinG1Params>,
        signer_round_1_private: &RoundOnePrivateOutput,
        signer_pubkeys: &[MultiSigPublicKey<HRegNon>],
        round_1_nonces: &[RoundOnePublicOutput],
    ) -> Option<RoundTwoPublicOutput> {
        if round_1_nonces.len() != signer_pubkeys.len() {
            return None;
        }
        if !valid_round1_nonces(round_1_nonces) {
            return None;
        }

        let aggregate_pubkey = Self::validate_and_combine_signer_pubkeys(signer_pubkeys)?;

        let a = multisig_generate_nonce_challenge::<HRegNon>(
            message,
            &aggregate_pubkey,
            round_1_nonces,
        );
        let _r = construct_multisig_nonce(&a, round_1_nonces);

        // Schnorr challenge using HSig
        let e_buf = schnorr_generate_challenge::<HSig>(message, &aggregate_pubkey, &_r);
        let e = GrumpkinFr::from_be_bytes(&e_buf);

        // z = r + s*a - sk*e
        let z = signer_round_1_private.r + signer_round_1_private.s * a - signer.private_key * e;
        Some(z)
    }

    /// Combine signature shares into a final Schnorr signature.
    ///
    /// Verifies the combined signature before returning.
    pub fn combine_signatures(
        message: &str,
        signer_pubkeys: &[MultiSigPublicKey<HRegNon>],
        round_1_nonces: &[RoundOnePublicOutput],
        round_2_shares: &[RoundTwoPublicOutput],
    ) -> Option<SchnorrSignature> {
        let num_signers = signer_pubkeys.len();
        if round_1_nonces.len() != num_signers || round_2_shares.len() != num_signers {
            return None;
        }
        if !valid_round1_nonces(round_1_nonces) {
            return None;
        }

        let aggregate_pubkey = Self::validate_and_combine_signer_pubkeys(signer_pubkeys)?;

        let a = multisig_generate_nonce_challenge::<HRegNon>(
            message,
            &aggregate_pubkey,
            round_1_nonces,
        );
        let r = construct_multisig_nonce(&a, round_1_nonces);

        let e_buf = schnorr_generate_challenge::<HSig>(message, &aggregate_pubkey, &r);

        // s = sum(z_i)
        let mut s = GrumpkinFr::zero();
        for z in round_2_shares {
            s = s + *z;
        }

        let sig = SchnorrSignature {
            s: s.to_be_bytes(),
            e: e_buf,
        };

        // Verify before returning
        if !schnorr_verify_signature::<HSig>(message, &aggregate_pubkey, &sig) {
            return None;
        }

        Some(sig)
    }
}

/// Generate the Fiat-Shamir nonce challenge for multisig.
///
/// a = H_non("h_nonce" || G || X_agg || "m_start" || len(m) || m || "m_end" || R1 || S1 || ... || Rn || Sn)
fn multisig_generate_nonce_challenge<H: Hasher>(
    message: &str,
    aggregate_pubkey: &GrumpkinAffine,
    round_1_nonces: &[RoundOnePublicOutput],
) -> GrumpkinFr {
    let mut buf = Vec::new();

    // Domain separator
    buf.extend_from_slice(b"h_nonce");

    // Generator point
    let g = GrumpkinAffine::one();
    buf.extend_from_slice(&g.to_buffer());

    // Aggregate public key
    buf.extend_from_slice(&aggregate_pubkey.to_buffer());

    // Message with framing
    buf.extend_from_slice(b"m_start");
    buf.extend_from_slice(&(message.len() as u32).to_be_bytes());
    buf.extend_from_slice(message.as_bytes());
    buf.extend_from_slice(b"m_end");

    // All nonces
    for nonce in round_1_nonces {
        buf.extend_from_slice(&nonce.r.to_buffer());
        buf.extend_from_slice(&nonce.s.to_buffer());
    }

    let hash_result = H::hash(&buf);
    let mut challenge_bytes = [0u8; 32];
    challenge_bytes.copy_from_slice(&hash_result[..32]);
    GrumpkinFr::from_be_bytes(&challenge_bytes)
}

/// Compute the combined multisig nonce: R = (R_1 + ... + R_n) + (S_1 + ... + S_n) * a
fn construct_multisig_nonce(
    a: &GrumpkinFr,
    round_1_nonces: &[RoundOnePublicOutput],
) -> GrumpkinAffine {
    let mut r_sum = Element::<GrumpkinG1Params>::from_affine(&round_1_nonces[0].r);
    let mut s_sum = Element::<GrumpkinG1Params>::from_affine(&round_1_nonces[0].s);

    for nonce in &round_1_nonces[1..] {
        r_sum = r_sum + Element::<GrumpkinG1Params>::from_affine(&nonce.r);
        s_sum = s_sum + Element::<GrumpkinG1Params>::from_affine(&nonce.s);
    }

    (r_sum + s_sum.mul_without_endomorphism(a)).to_affine()
}

/// Validate that all round-1 nonces are on-curve, non-infinity, and no duplicates.
fn valid_round1_nonces(nonces: &[RoundOnePublicOutput]) -> bool {
    for nonce in nonces {
        if !nonce.r.on_curve() || nonce.r.is_point_at_infinity() {
            return false;
        }
        if !nonce.s.on_curve() || nonce.s.is_point_at_infinity() {
            return false;
        }
    }
    // Check for duplicates
    for i in 0..nonces.len() {
        for j in (i + 1)..nonces.len() {
            if nonces[i] == nonces[j] {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hashers::{Blake2sHasher, KeccakHasher, Sha256Hasher};
    use bbrs_ecc::fields::field::Field;

    fn make_schnorr_keypair() -> SchnorrKeyPair<GrumpkinG1Params> {
        let private_key = GrumpkinFr::random_element();
        let generator = Element::<GrumpkinG1Params>::one();
        let public_key = generator
            .mul_without_endomorphism(&private_key)
            .to_affine();
        SchnorrKeyPair {
            private_key,
            public_key,
        }
    }

    // -----------------------------------------------------------------------
    // Basic Schnorr sign/verify tests
    // -----------------------------------------------------------------------

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

    #[test]
    fn schnorr_sign_verify_keccak() {
        let account = make_schnorr_keypair();
        let message = "The quick brown fox jumped over the lazy dog.";
        let sig = schnorr_construct_signature::<KeccakHasher>(message, &account);
        let valid =
            schnorr_verify_signature::<KeccakHasher>(message, &account.public_key, &sig);
        assert!(valid);
    }

    #[test]
    fn schnorr_sign_verify_sha256() {
        let account = make_schnorr_keypair();
        let message = "The quick brown dog jumped over the lazy fox.";
        let sig = schnorr_construct_signature::<Sha256Hasher>(message, &account);
        let valid =
            schnorr_verify_signature::<Sha256Hasher>(message, &account.public_key, &sig);
        assert!(valid);
    }

    #[test]
    fn schnorr_signature_consistency() {
        let message_a = "The quick brown fox jumped over the lazy dog.";
        let message_b = "The quick brown dog jumped over the lazy fox.";

        let account_a = make_schnorr_keypair();
        let account_b = make_schnorr_keypair();

        assert!(account_a.private_key != account_b.private_key);
        assert!(account_a.public_key != account_b.public_key);

        // Same message, same account — random k means different sigs
        let sig_a =
            schnorr_construct_signature::<KeccakHasher>(message_a, &account_a);
        let sig_b =
            schnorr_construct_signature::<KeccakHasher>(message_a, &account_a);
        assert!(sig_a.e != sig_b.e);
        assert!(sig_a.s != sig_b.s);

        // Same message, different accounts
        let sig_c =
            schnorr_construct_signature::<KeccakHasher>(message_a, &account_a);
        let sig_d =
            schnorr_construct_signature::<KeccakHasher>(message_a, &account_b);
        assert!(sig_c.e != sig_d.e);
        assert!(sig_c.s != sig_d.s);

        // Different message, same account
        let sig_e =
            schnorr_construct_signature::<KeccakHasher>(message_a, &account_a);
        let sig_f =
            schnorr_construct_signature::<KeccakHasher>(message_b, &account_a);
        assert!(sig_e.e != sig_f.e);
        assert!(sig_e.s != sig_f.s);

        // Different message, different accounts
        let sig_g =
            schnorr_construct_signature::<KeccakHasher>(message_a, &account_a);
        let sig_h =
            schnorr_construct_signature::<KeccakHasher>(message_b, &account_b);
        assert!(sig_g.e != sig_h.e);
        assert!(sig_g.s != sig_h.s);

        // All signatures verify
        assert!(schnorr_verify_signature::<KeccakHasher>(
            message_a, &account_a.public_key, &sig_a
        ));
        assert!(schnorr_verify_signature::<KeccakHasher>(
            message_a, &account_a.public_key, &sig_b
        ));
        assert!(schnorr_verify_signature::<KeccakHasher>(
            message_a, &account_a.public_key, &sig_c
        ));
        assert!(schnorr_verify_signature::<KeccakHasher>(
            message_a, &account_b.public_key, &sig_d
        ));
        assert!(schnorr_verify_signature::<KeccakHasher>(
            message_a, &account_a.public_key, &sig_e
        ));
        assert!(schnorr_verify_signature::<KeccakHasher>(
            message_b, &account_a.public_key, &sig_f
        ));
        assert!(schnorr_verify_signature::<KeccakHasher>(
            message_a, &account_a.public_key, &sig_g
        ));
        assert!(schnorr_verify_signature::<KeccakHasher>(
            message_b, &account_b.public_key, &sig_h
        ));
    }

    // -----------------------------------------------------------------------
    // Proof of Possession tests
    // -----------------------------------------------------------------------

    macro_rules! pop_tests {
        ($hasher:ty, $suffix:ident) => {
            paste::paste! {
                #[test]
                fn [<pop_valid_proof_ $suffix>]() {
                    let account = make_schnorr_keypair();
                    let pop = SchnorrProofOfPossession::new::<$hasher>(&account);
                    assert!(pop.verify::<$hasher>(&account.public_key));
                }

                #[test]
                fn [<pop_invalid_empty_ $suffix>]() {
                    let account = make_schnorr_keypair();
                    let pop = SchnorrProofOfPossession {
                        challenge: [0u8; 32],
                        response: GrumpkinFr::zero(),
                    };
                    assert!(!pop.verify::<$hasher>(&account.public_key));
                }

                #[test]
                fn [<pop_wrong_account_ $suffix>]() {
                    let account1 = make_schnorr_keypair();
                    let account2 = make_schnorr_keypair();
                    let pop = SchnorrProofOfPossession::new::<$hasher>(&account1);
                    assert!(!pop.verify::<$hasher>(&account2.public_key));
                }

                #[test]
                fn [<pop_zero_challenge_ $suffix>]() {
                    let account = make_schnorr_keypair();
                    let mut pop = SchnorrProofOfPossession::new::<$hasher>(&account);
                    pop.challenge = [0u8; 32];
                    assert!(!pop.verify::<$hasher>(&account.public_key));
                }

                #[test]
                fn [<pop_zero_response_ $suffix>]() {
                    let account = make_schnorr_keypair();
                    let mut pop = SchnorrProofOfPossession::new::<$hasher>(&account);
                    pop.response = GrumpkinFr::zero();
                    assert!(!pop.verify::<$hasher>(&account.public_key));
                }
            }
        };
    }

    pop_tests!(KeccakHasher, keccak);
    pop_tests!(Sha256Hasher, sha256);
    pop_tests!(Blake2sHasher, blake2s);

    // -----------------------------------------------------------------------
    // Multisig tests
    // -----------------------------------------------------------------------

    fn create_multisig<HRegNon: Hasher>(
        message: &str,
        accounts: &[SchnorrKeyPair<GrumpkinG1Params>],
        tamper_pop: bool,
    ) -> Option<SchnorrSignature> {
        let num_signers = accounts.len();

        let mut signer_pubkeys: Vec<MultiSigPublicKey<HRegNon>> = accounts
            .iter()
            .map(|a| MultiSigPublicKey::<HRegNon>::new(a))
            .collect();

        if tamper_pop {
            signer_pubkeys[0].proof_of_possession.response =
                signer_pubkeys[0].proof_of_possession.response + Field::one();
        }

        let mut round1_pub = Vec::with_capacity(num_signers);
        let mut round1_priv = Vec::with_capacity(num_signers);

        for _ in 0..num_signers {
            let (pub_out, priv_out) =
                SchnorrMultisig::<HRegNon, Blake2sHasher>::construct_signature_round_1();
            round1_pub.push(pub_out);
            round1_priv.push(priv_out);
        }

        let mut round2 = Vec::with_capacity(num_signers);
        for i in 0..num_signers {
            if let Some(z) =
                SchnorrMultisig::<HRegNon, Blake2sHasher>::construct_signature_round_2(
                    message,
                    &accounts[i],
                    &round1_priv[i],
                    &signer_pubkeys,
                    &round1_pub,
                )
            {
                round2.push(z);
            }
        }

        SchnorrMultisig::<HRegNon, Blake2sHasher>::combine_signatures(
            message,
            &signer_pubkeys,
            &round1_pub,
            &round2,
        )
    }

    macro_rules! multisig_tests {
        ($hasher:ty, $suffix:ident) => {
            paste::paste! {
                #[test]
                fn [<multisig_5_signers_ $suffix>]() {
                    let message = "The quick brown dog jumped over the lazy fox.";
                    let num_signers = 5;
                    let accounts: Vec<_> =
                        (0..num_signers).map(|_| make_schnorr_keypair()).collect();

                    let sig = create_multisig::<$hasher>(message, &accounts, false);
                    assert!(sig.is_some(), "Multisig should produce a valid signature");

                    // Verify with aggregate key
                    let signer_pubkeys: Vec<_> = accounts
                        .iter()
                        .map(|a| MultiSigPublicKey::<$hasher>::new(a))
                        .collect();
                    let agg_key = SchnorrMultisig::<$hasher, Blake2sHasher>::validate_and_combine_signer_pubkeys(
                        &signer_pubkeys,
                    );
                    assert!(agg_key.is_some());
                    let valid = schnorr_verify_signature::<Blake2sHasher>(
                        message,
                        &agg_key.unwrap(),
                        &sig.unwrap(),
                    );
                    assert!(valid, "Combined multisig should verify as regular Schnorr");
                }

                #[test]
                fn [<multisig_invalid_pop_ $suffix>]() {
                    let message = "The quick brown dog jumped over the lazy fox.";
                    let num_signers = 5;
                    let accounts: Vec<_> =
                        (0..num_signers).map(|_| make_schnorr_keypair()).collect();

                    let sig = create_multisig::<$hasher>(message, &accounts, true);
                    assert!(sig.is_none(), "Multisig should fail with invalid PoP");
                }

                #[test]
                fn [<multisig_duplicate_keys_ $suffix>]() {
                    let message = "The quick brown dog jumped over the lazy fox.";
                    let num_signers = 5;
                    let mut accounts: Vec<_> =
                        (0..num_signers).map(|_| make_schnorr_keypair()).collect();

                    // Duplicate account[4] into account[2]
                    accounts[2] = SchnorrKeyPair {
                        private_key: accounts[4].private_key,
                        public_key: accounts[4].public_key,
                    };

                    let sig = create_multisig::<$hasher>(message, &accounts, false);
                    assert!(sig.is_none(), "Multisig should fail with duplicate keys");
                }
            }
        };
    }

    multisig_tests!(KeccakHasher, keccak);
    multisig_tests!(Sha256Hasher, sha256);
}
