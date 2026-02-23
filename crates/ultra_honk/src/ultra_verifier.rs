//! Port of `ultra_verifier.hpp`/`.cpp` — Ultra Honk Verifier.
//!
//! The Ultra verifier orchestrates the full proof verification:
//! 1. Oink phase: receive commitments, generate challenges
//! 2. Sumcheck phase: run sumcheck verifier
//! 3. PCS phase: run Shplemini verification + KZG pairing check

use bbrs_commitment_schemes::claim_batcher::{Batch, ClaimBatcher};
use bbrs_commitment_schemes::kzg::KZG;
use bbrs_commitment_schemes::shplonk::shplemini::ShpleminiVerifier;
use bbrs_ecc::curves::bn254::{Bn254G1Params, Fr, G1Affine};
use bbrs_transcript::NativeTranscript;

use crate::oink_verifier::OinkVerifier;
use crate::ultra_sumcheck::UltraSumcheckVerifier;
use crate::verification_key::VerificationKey;

/// Output of the Ultra Honk verifier.
pub struct UltraVerifierOutput {
    /// Whether the proof verified successfully.
    pub result: bool,
}

/// Ultra Honk verifier.
///
/// Port of C++ `UltraVerifier_<UltraFlavor>`.
pub struct UltraVerifier {
    pub verification_key: VerificationKey,
}

impl UltraVerifier {
    pub fn new(verification_key: VerificationKey) -> Self {
        Self { verification_key }
    }

    /// Verify a proof.
    pub fn verify(self, proof: &[Fr]) -> UltraVerifierOutput {
        let mut transcript = NativeTranscript::from_proof(proof);

        // Phase 1: Oink — receive commitments, generate challenges
        let oink_output = {
            let oink = OinkVerifier::new(&self.verification_key);
            oink.verify(&mut transcript)
        };

        let log_circuit_size = self.verification_key.log_circuit_size;

        // Get gate challenges for pow-polynomial
        let gate_challenges: Vec<Fr> = (0..log_circuit_size)
            .map(|i| transcript.get_challenge(&format!("Sumcheck:gate_challenge_{}", i)))
            .collect();

        // Phase 2: Sumcheck
        let sumcheck_output = {
            let sumcheck_verifier = UltraSumcheckVerifier::new(
                oink_output.alpha,
                log_circuit_size,
            );
            sumcheck_verifier.verify(
                &mut transcript,
                &oink_output.relation_parameters,
                &gate_challenges,
            )
        };

        if !sumcheck_output.verified {
            #[cfg(test)]
            eprintln!("[UltraVerifier] SUMCHECK FAILED");
            return UltraVerifierOutput { result: false };
        }
        #[cfg(test)]
        eprintln!("[UltraVerifier] Sumcheck passed");

        // Phase 3: PCS verification (Shplemini + KZG)
        // Build claim batcher with all commitments and evaluations
        let evals = &sumcheck_output.claimed_evaluations;

        // Collect unshifted commitments: precomputed (from VK) + witness (from oink)
        let vk = &self.verification_key;
        let precomputed_commitments = vk.get_all_commitments();
        let witness_commitments = [
            oink_output.witness_commitments.w_l,
            oink_output.witness_commitments.w_r,
            oink_output.witness_commitments.w_o,
            oink_output.witness_commitments.w_4,
            oink_output.witness_commitments.z_perm,
            oink_output.witness_commitments.lookup_inverses,
            oink_output.witness_commitments.lookup_read_counts,
            oink_output.witness_commitments.lookup_read_tags,
        ];

        let all_unshifted_commitments: Vec<G1Affine> = precomputed_commitments
            .iter()
            .chain(witness_commitments.iter())
            .copied()
            .collect();

        // Unshifted evaluations (precomputed + witness)
        let unshifted_evals: Vec<Fr> = evals
            .get_precomputed()
            .iter()
            .chain(evals.get_witness().iter())
            .map(|v| **v)
            .collect();

        // Shifted commitments (same polynomials as to_be_shifted, reuse their commitments)
        let shifted_commitments: Vec<G1Affine> = vec![
            oink_output.witness_commitments.w_l,
            oink_output.witness_commitments.w_r,
            oink_output.witness_commitments.w_o,
            oink_output.witness_commitments.w_4,
            oink_output.witness_commitments.z_perm,
        ];

        // Shifted evaluations
        let shifted_evals: Vec<Fr> = evals
            .get_shifted()
            .iter()
            .map(|v| **v)
            .collect();

        // Build claim batcher
        let mut claim_batcher = ClaimBatcher::<Bn254G1Params>::new();
        claim_batcher.unshifted = Some(Batch {
            commitments: all_unshifted_commitments,
            evaluations: unshifted_evals,
            scalar: Fr::one(),
        });
        claim_batcher.shifted = Some(Batch {
            commitments: shifted_commitments,
            evaluations: shifted_evals,
            scalar: Fr::one(),
        });

        // Padding indicator array (all ones for non-ZK)
        let padding_indicator_array: Vec<Fr> = vec![Fr::one(); log_circuit_size];

        // Get G1 identity point
        let mut vk_pcs = self.verification_key.pcs_verification_key;
        vk_pcs.initialize();
        let g1_identity = vk_pcs.get_g1_identity();

        // Run ShpleminiVerifier
        let batch_opening_claim = ShpleminiVerifier::compute_batch_opening_claim_without_zk(
            &padding_indicator_array,
            &mut claim_batcher,
            &sumcheck_output.challenge,
            &g1_identity,
            &mut transcript,
        );

        // Run KZG batch verification
        let pairing_points = KZG::reduce_verify_batch_opening_claim(
            batch_opening_claim,
            &mut transcript,
        );

        // Final pairing check
        let result = pairing_points.check(&mut vk_pcs);
        #[cfg(test)]
        eprintln!("[UltraVerifier] Pairing check result: {}", result);

        UltraVerifierOutput { result }
    }
}
