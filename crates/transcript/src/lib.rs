//! Barretenberg transcript crate.
//!
//! Port of C++ `transcript/transcript.hpp` and `ecc/fields/field_conversion.hpp`.
//! Provides Fiat-Shamir transcript functionality for prover-verifier communication.

pub mod codec;
pub mod manifest;
pub mod transcript;

use bbrs_ecc::curves::bn254::Fr;

/// Poseidon2-based hasher for the native transcript.
pub struct Poseidon2Hasher;

impl transcript::TranscriptHasher for Poseidon2Hasher {
    fn hash(input: &[Fr]) -> Fr {
        bbrs_crypto::poseidon2::Poseidon2::hash(input)
    }
}

/// NativeTranscript: BaseTranscript using FrCodec + Poseidon2.
///
/// This is the primary transcript type for native (non-circuit) proving/verification.
/// Equivalent to C++ `NativeTranscript = BaseTranscript<FrCodec, Poseidon2<Bn254ScalarFieldParams>>`.
pub type NativeTranscript = transcript::BaseTranscript<Poseidon2Hasher>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codec::FieldSerializable;

    /// Test: prover sends data, exports proof, verifier loads and receives same data,
    /// both generate identical challenges.
    #[test]
    fn test_prover_verifier_deterministic_challenge() {
        // Prover side
        let mut prover = NativeTranscript::new();

        let val_a = Fr::from(100u64);
        let val_b = Fr::from(200u64);
        let val_c: u32 = 42;

        prover.send_to_verifier("a", &val_a);
        prover.send_to_verifier("b", &val_b);
        prover.send_to_verifier("c", &val_c);

        let prover_challenge = prover.get_challenge("alpha");

        // Export proof
        let proof = prover.export_proof();

        // Verifier side
        let mut verifier = NativeTranscript::from_proof(&proof);

        let recv_a: Fr = verifier.receive_from_prover("a");
        let recv_b: Fr = verifier.receive_from_prover("b");
        let recv_c: u32 = verifier.receive_from_prover("c");

        assert_eq!(recv_a, val_a);
        assert_eq!(recv_b, val_b);
        assert_eq!(recv_c, val_c);

        let verifier_challenge = verifier.get_challenge("alpha");

        assert_eq!(
            prover_challenge, verifier_challenge,
            "Prover and verifier must derive identical challenges"
        );

        // Challenge must not be zero (extremely unlikely for Poseidon2)
        assert_ne!(prover_challenge, Fr::zero());
    }

    /// Test: multiple rounds of send + challenge generation.
    #[test]
    fn test_multi_round_challenges() {
        let mut prover = NativeTranscript::new();
        let mut challenges_prover = Vec::new();

        // Round 0
        prover.send_to_verifier("poly_0", &Fr::from(1u64));
        challenges_prover.push(prover.get_challenge("r0"));

        // Round 1
        prover.send_to_verifier("poly_1", &Fr::from(2u64));
        challenges_prover.push(prover.get_challenge("r1"));

        // Round 2
        prover.send_to_verifier("poly_2", &Fr::from(3u64));
        challenges_prover.push(prover.get_challenge("r2"));

        let proof = prover.export_proof();

        // Verifier
        let mut verifier = NativeTranscript::from_proof(&proof);
        let mut challenges_verifier = Vec::new();

        let _: Fr = verifier.receive_from_prover("poly_0");
        challenges_verifier.push(verifier.get_challenge("r0"));

        let _: Fr = verifier.receive_from_prover("poly_1");
        challenges_verifier.push(verifier.get_challenge("r1"));

        let _: Fr = verifier.receive_from_prover("poly_2");
        challenges_verifier.push(verifier.get_challenge("r2"));

        for i in 0..3 {
            assert_eq!(
                challenges_prover[i], challenges_verifier[i],
                "Challenge mismatch at round {i}"
            );
        }

        // All challenges should be distinct
        assert_ne!(challenges_prover[0], challenges_prover[1]);
        assert_ne!(challenges_prover[1], challenges_prover[2]);
    }

    /// Test: multiple challenges per round (get_challenges with multiple labels).
    #[test]
    fn test_multiple_challenges_per_round() {
        let mut prover = NativeTranscript::new();
        prover.send_to_verifier("data", &Fr::from(999u64));
        let p_challenges = prover.get_challenges(&["alpha", "beta", "gamma"]);

        let proof = prover.export_proof();
        let mut verifier = NativeTranscript::from_proof(&proof);
        let _: Fr = verifier.receive_from_prover("data");
        let v_challenges = verifier.get_challenges(&["alpha", "beta", "gamma"]);

        assert_eq!(p_challenges.len(), 3);
        assert_eq!(v_challenges.len(), 3);
        for i in 0..3 {
            assert_eq!(p_challenges[i], v_challenges[i]);
        }
        // All three should be distinct
        assert_ne!(p_challenges[0], p_challenges[1]);
        assert_ne!(p_challenges[1], p_challenges[2]);
        assert_ne!(p_challenges[0], p_challenges[2]);
    }

    /// Test: prover_init_empty + verifier_init_empty convenience functions.
    #[test]
    fn test_init_empty_roundtrip() {
        let mut prover = NativeTranscript::prover_init_empty();
        let p_challenge = prover.get_challenge("test");

        let mut verifier = NativeTranscript::verifier_init_empty(&mut prover);
        let v_challenge = verifier.get_challenge("test");

        assert_eq!(p_challenge, v_challenge);
        assert_ne!(p_challenge, Fr::zero());
    }

    /// Test: export_proof returns correct slices for shared transcript.
    #[test]
    fn test_export_proof_slicing() {
        let mut transcript = NativeTranscript::new();

        // First prover sends 2 elements
        transcript.send_to_verifier("a", &Fr::from(1u64));
        transcript.send_to_verifier("b", &Fr::from(2u64));
        let slice1 = transcript.export_proof();
        assert_eq!(slice1.len(), 2);

        // Second prover sends 3 elements
        transcript.send_to_verifier("c", &Fr::from(3u64));
        transcript.send_to_verifier("d", &Fr::from(4u64));
        transcript.send_to_verifier("e", &Fr::from(5u64));
        let slice2 = transcript.export_proof();
        assert_eq!(slice2.len(), 3);
    }

    /// Test: Grumpkin Fr (BN254 Fq) roundtrip through transcript.
    #[test]
    fn test_grumpkin_fr_through_transcript() {
        use bbrs_ecc::curves::grumpkin::GrumpkinFr;

        let val = GrumpkinFr::from(0xDEADBEEF_CAFEBABEu64);
        let mut prover = NativeTranscript::new();
        prover.send_to_verifier("grumpkin_val", &val);
        let challenge = prover.get_challenge("ch");

        let proof = prover.export_proof();
        let mut verifier = NativeTranscript::from_proof(&proof);
        let recv: GrumpkinFr = verifier.receive_from_prover("grumpkin_val");
        let v_challenge = verifier.get_challenge("ch");

        assert_eq!(val, recv);
        assert_eq!(challenge, v_challenge);
    }

    /// Test: challenge determinism — same inputs always produce same challenges.
    #[test]
    fn test_challenge_determinism() {
        let make_transcript = || {
            let mut t = NativeTranscript::new();
            t.send_to_verifier("x", &Fr::from(42u64));
            t.send_to_verifier("y", &Fr::from(99u64));
            t.get_challenge("alpha")
        };

        let c1 = make_transcript();
        let c2 = make_transcript();
        assert_eq!(c1, c2, "Same transcript operations must produce same challenge");
    }

    /// Test: different inputs produce different challenges.
    #[test]
    fn test_challenge_sensitivity() {
        let mut t1 = NativeTranscript::new();
        t1.send_to_verifier("x", &Fr::from(1u64));
        let c1 = t1.get_challenge("alpha");

        let mut t2 = NativeTranscript::new();
        t2.send_to_verifier("x", &Fr::from(2u64));
        let c2 = t2.get_challenge("alpha");

        assert_ne!(c1, c2, "Different inputs must produce different challenges");
    }

    /// Test: AffineElement serialization through transcript.
    #[test]
    fn test_grumpkin_affine_through_transcript() {
        use bbrs_ecc::curves::grumpkin;

        // Use the Grumpkin generator point (base field = BN254 Fr, so 1 Fr per coord = 2 Fr total)
        let point = grumpkin::G1Affine::one();
        assert_eq!(<grumpkin::G1Affine as FieldSerializable>::NUM_FR, 2);

        let mut prover = NativeTranscript::new();
        prover.send_to_verifier("commit", &point);
        let challenge = prover.get_challenge("ch");

        let proof = prover.export_proof();
        let mut verifier = NativeTranscript::from_proof(&proof);
        let recv: grumpkin::G1Affine = verifier.receive_from_prover("commit");
        let v_challenge = verifier.get_challenge("ch");

        assert_eq!(point, recv);
        assert_eq!(challenge, v_challenge);
    }

    /// Test: BN254 affine point through transcript.
    #[test]
    fn test_bn254_affine_through_transcript() {
        use bbrs_ecc::curves::bn254;

        let point = bn254::G1Affine::one();
        // BN254 G1 base field = Fq (needs 2 Fr per coord) → 4 Fr total
        assert_eq!(<bn254::G1Affine as FieldSerializable>::NUM_FR, 4);

        let mut prover = NativeTranscript::new();
        prover.send_to_verifier("commit", &point);
        let challenge = prover.get_challenge("ch");

        let proof = prover.export_proof();
        let mut verifier = NativeTranscript::from_proof(&proof);
        let recv: bn254::G1Affine = verifier.receive_from_prover("commit");
        let v_challenge = verifier.get_challenge("ch");

        assert_eq!(point, recv);
        assert_eq!(challenge, v_challenge);
    }
}
