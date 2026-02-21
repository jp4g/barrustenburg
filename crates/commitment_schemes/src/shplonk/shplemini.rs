//! Shplemini: combined Gemini + Shplonk prover and verifier.
//!
//! C++ source: barretenberg/commitment_schemes/shplonk/shplemini.hpp
//!
//! Combines verifiers from four protocols:
//! 1. Batch opening protocol: reduces multilinear claims to a single batched polynomial
//! 2. Gemini protocol: reduces the batched polynomial to Gemini univariate openings
//! 3. Shplonk protocol: reduces Gemini univariate openings to a single batched univariate
//! 4. KZG or IPA protocol: verifies the final univariate opening

use bbrs_ecc::curves::bn254::{Bn254FrParams, Bn254G1Params, Fr, G1Affine as Bn254G1Affine};
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;
use bbrs_transcript::NativeTranscript;

use crate::claim::{BatchOpeningClaim, OpeningPair};
use crate::claim_batcher::ClaimBatcher;
use crate::commitment_key::CommitmentKey;
use crate::gemini::{self, GeminiProver, GeminiVerifier, PolynomialBatcher};

use super::{
    compute_shplonk_batching_challenge_powers, ProverOpeningClaim, ShplonkProver, ShplonkVerifier,
};

// ── ShpleminiProver ──────────────────────────────────────────────────────────

/// Combined Gemini + Shplonk prover.
///
/// C++ source: `ShpleminiProver_<Curve>`
pub struct ShpleminiProver;

impl ShpleminiProver {
    /// Full Shplemini prove: Gemini reduction + Shplonk batching.
    ///
    /// C++ source: `ShpleminiProver_::prove`
    ///
    /// Specialized for BN254 because NativeTranscript methods are BN254-specific.
    pub fn prove(
        circuit_size: usize,
        batcher: &mut PolynomialBatcher<Bn254FrParams>,
        multilinear_challenge: &[Fr],
        commitment_key: &CommitmentKey<Bn254G1Params>,
        transcript: &mut NativeTranscript,
    ) -> ProverOpeningClaim<Bn254FrParams> {
        let virtual_log_n = multilinear_challenge.len();
        let log_n = (circuit_size as f64).log2() as usize;

        // Get rho challenge for batching multilinear polynomials
        let rho: Fr = transcript.get_challenge("rho");

        // Compute A_0 = F + G/X (batched polynomial)
        let a_0 = batcher.compute_batched(rho);

        // Compute Gemini fold polynomials
        let fold_polynomials =
            GeminiProver::compute_fold_polynomials(log_n, multilinear_challenge, &a_0);

        // Commit and send fold polynomial commitments
        for (l, fold_poly) in fold_polynomials.iter().enumerate().take(log_n - 1) {
            let commitment = commitment_key.commit(fold_poly);
            let label = format!("Gemini:FOLD_{}", l + 1);
            transcript.send_to_verifier(&label, &commitment);
        }

        let r: Fr = transcript.get_challenge("Gemini:r");

        // Run Gemini prover to get opening claims
        let gemini_claims =
            GeminiProver::prove(log_n, batcher, multilinear_challenge, rho, r);

        // Convert gemini claims to shplonk claims
        let mut opening_claims: Vec<ProverOpeningClaim<Bn254FrParams>> = gemini_claims
            .into_iter()
            .map(|c| ProverOpeningClaim {
                polynomial: c.polynomial,
                opening_pair: OpeningPair {
                    challenge: c.opening_pair.challenge,
                    evaluation: c.opening_pair.evaluation,
                },
                gemini_fold: c.gemini_fold,
            })
            .collect();

        // Send Gemini evaluations to transcript
        for claim in opening_claims.iter() {
            if !claim.gemini_fold {
                // A_0+(r) and A_0-(-r) evaluations - sent as part of the protocol
            }
        }

        let mut libra_opening_claims = Vec::new();
        let mut sumcheck_round_claims = Vec::new();

        ShplonkProver::prove(
            commitment_key,
            &mut opening_claims,
            transcript,
            &mut libra_opening_claims,
            &mut sumcheck_round_claims,
            virtual_log_n,
        )
    }
}

// ── ShpleminiVerifier ────────────────────────────────────────────────────────

/// Combined Gemini + Shplonk verifier for efficient batch verification.
///
/// C++ source: `ShpleminiVerifier_<Curve>`
pub struct ShpleminiVerifier;

impl ShpleminiVerifier {
    /// Compute the batch opening claim combining Gemini, Shplonk, and optional ZK data.
    ///
    /// C++ source: `ShpleminiVerifier_::compute_batch_opening_claim`
    ///
    /// Specialized for BN254 because NativeTranscript methods are BN254-specific.
    pub fn compute_batch_opening_claim(
        padding_indicator_array: &[Fr],
        claim_batcher: &mut ClaimBatcher<Bn254G1Params>,
        multivariate_challenge: &[Fr],
        g1_identity: &Bn254G1Affine,
        transcript: &mut NativeTranscript,
    ) -> BatchOpeningClaim<Bn254G1Params> {
        let virtual_log_n = multivariate_challenge.len();
        let mut batched_evaluation = Fr::zero();

        // Get rho challenge for batching multilinear polynomials
        let gemini_batching_challenge: Fr = transcript.get_challenge("rho");

        // Get Gemini fold commitments
        let mut fold_commitments: Vec<Bn254G1Affine> = Vec::with_capacity(virtual_log_n - 1);
        for l in 1..virtual_log_n {
            let label = format!("Gemini:FOLD_{}", l);
            let comm: Bn254G1Affine = transcript.receive_from_prover(&label);
            fold_commitments.push(comm);
        }

        // Get Gemini evaluation challenge
        let gemini_evaluation_challenge: Fr = transcript.get_challenge("Gemini:r");

        // Get Gemini fold negative evaluations
        let mut gemini_fold_neg_evaluations: Vec<Fr> = Vec::with_capacity(virtual_log_n);
        for l in 0..virtual_log_n {
            let label = format!("Gemini:a_{}", l);
            let eval: Fr = transcript.receive_from_prover(&label);
            gemini_fold_neg_evaluations.push(eval);
        }

        // Compute powers of r: (r, r^2, r^4, ...)
        let gemini_eval_challenge_powers =
            gemini::powers_of_evaluation_challenge(gemini_evaluation_challenge, virtual_log_n);

        // Process Shplonk transcript data
        let shplonk_batching_challenge: Fr = transcript.get_challenge("Shplonk:nu");

        let shplonk_batching_challenge_powers = compute_shplonk_batching_challenge_powers(
            &shplonk_batching_challenge,
            virtual_log_n,
            false,
            false,
        );

        let q_commitment: Bn254G1Affine = transcript.receive_from_prover("Shplonk:Q");

        // Start populating commitments with Q
        let mut commitments: Vec<Bn254G1Affine> = vec![q_commitment];

        // Get Shplonk opening point z
        let shplonk_evaluation_challenge: Fr = transcript.get_challenge("Shplonk:z");

        let mut constant_term_accumulator = Fr::zero();
        let mut scalars: Vec<Fr> = vec![Fr::one()];

        // Compute inverted Gemini denominators
        let inverse_vanishing_evals =
            ShplonkVerifier::<Bn254G1Params>::compute_inverted_gemini_denominators(
                &shplonk_evaluation_challenge,
                &gemini_eval_challenge_powers,
            );

        // Compute scalars for each batch
        claim_batcher.compute_scalars_for_each_batch(
            &inverse_vanishing_evals,
            &shplonk_batching_challenge,
            &gemini_evaluation_challenge,
        );

        // No interleaving in the basic case
        let shplonk_batching_pos = Fr::zero();
        let shplonk_batching_neg = Fr::zero();

        // Update commitments and scalars from claim batcher
        claim_batcher.update_batch_mul_inputs_and_batched_evaluation(
            &mut commitments,
            &mut scalars,
            &mut batched_evaluation,
            &gemini_batching_challenge,
            shplonk_batching_pos,
            shplonk_batching_neg,
        );

        // Reconstruct fold positive evaluations
        let p_neg = Fr::zero();
        let gemini_fold_pos_evaluations = GeminiVerifier::compute_fold_pos_evaluations(
            padding_indicator_array,
            batched_evaluation,
            multivariate_challenge,
            &gemini_eval_challenge_powers,
            &gemini_fold_neg_evaluations,
            p_neg,
        );

        // Batch Gemini claims from prover
        Self::batch_gemini_claims_received_from_prover::<Bn254G1Params>(
            padding_indicator_array,
            &fold_commitments,
            &gemini_fold_neg_evaluations,
            &gemini_fold_pos_evaluations,
            &inverse_vanishing_evals,
            &shplonk_batching_challenge_powers,
            &mut commitments,
            &mut scalars,
            &mut constant_term_accumulator,
        );

        // A_0 contributions
        let full_a_0_pos = gemini_fold_pos_evaluations[0];
        let a_0_pos = full_a_0_pos; // No interleaving: p_pos = 0

        // Add A_0+(r)/(z-r) to constant term
        constant_term_accumulator =
            constant_term_accumulator + a_0_pos * inverse_vanishing_evals[0];
        // Add A_0-(-r)/(z+r) * nu to constant term
        constant_term_accumulator = constant_term_accumulator
            + gemini_fold_neg_evaluations[0]
                * shplonk_batching_challenge
                * inverse_vanishing_evals[1];

        // Finalize
        commitments.push(*g1_identity);
        scalars.push(constant_term_accumulator);

        BatchOpeningClaim {
            commitments,
            scalars,
            evaluation_point: shplonk_evaluation_challenge,
        }
    }

    /// Place fold polynomial commitments and compute corresponding scalar multipliers.
    ///
    /// C++ source: `ShpleminiVerifier_::batch_gemini_claims_received_from_prover`
    pub fn batch_gemini_claims_received_from_prover<C: CurveParams>(
        padding_indicator_array: &[Field<C::ScalarFieldParams>],
        fold_commitments: &[AffineElement<C>],
        gemini_neg_evaluations: &[Field<C::ScalarFieldParams>],
        gemini_pos_evaluations: &[Field<C::ScalarFieldParams>],
        inverse_vanishing_evals: &[Field<C::ScalarFieldParams>],
        shplonk_batching_challenge_powers: &[Field<C::ScalarFieldParams>],
        commitments: &mut Vec<AffineElement<C>>,
        scalars: &mut Vec<Field<C::ScalarFieldParams>>,
        constant_term_accumulator: &mut Field<C::ScalarFieldParams>,
    ) {
        let virtual_log_n = gemini_neg_evaluations.len();

        // Start from j=1 (A_0 is handled separately)
        for j in 1..virtual_log_n {
            let pos_index = 2 * j;
            let neg_index = 2 * j + 1;

            // (nu^{2j}) / (z - r^{2^j})
            let scaling_factor_pos = shplonk_batching_challenge_powers[pos_index]
                * inverse_vanishing_evals[pos_index];
            // (nu^{2j+1}) / (z + r^{2^j})
            let scaling_factor_neg = shplonk_batching_challenge_powers[neg_index]
                * inverse_vanishing_evals[neg_index];

            // Accumulate constant term contribution
            *constant_term_accumulator = *constant_term_accumulator
                + scaling_factor_neg * gemini_neg_evaluations[j]
                + scaling_factor_pos * gemini_pos_evaluations[j];

            // Place scaling factor and commitment
            scalars.push(
                -padding_indicator_array[j] * (scaling_factor_neg + scaling_factor_pos),
            );
            commitments.push(fold_commitments[j - 1]);
        }
    }
}
