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
use bbrs_polynomials::polynomial::Polynomial;
use bbrs_transcript::NativeTranscript;

use crate::claim::{BatchOpeningClaim, OpeningPair};
use crate::claim_batcher::ClaimBatcher;
use crate::commitment_key::CommitmentKey;
use crate::gemini::{self, GeminiProver, GeminiVerifier, PolynomialBatcher};
use crate::small_subgroup_ipa::{SmallSubgroupCurveConfig, SmallSubgroupIPAVerifier};

use super::{
    compute_shplonk_batching_challenge_powers, ProverOpeningClaim, ShplonkProver, ShplonkVerifier,
    NUM_INTERLEAVING_CLAIMS, NUM_LIBRA_COMMITMENTS, NUM_SMALL_IPA_EVALUATIONS,
};

// ── RepeatedCommitmentsData ─────────────────────────────────────────────────

/// Data for combining scalars of repeating commitments.
///
/// Port of C++ `RepeatedCommitmentsData` from `flavor/repeated_commitments_data.hpp`.
#[derive(Default)]
pub struct RepeatedCommitmentsData {
    pub first_range_to_be_shifted_start: usize,
    pub first_range_shifted_start: usize,
    pub first_range_size: usize,
    pub second_range_to_be_shifted_start: usize,
    pub second_range_shifted_start: usize,
    pub second_range_size: usize,
}

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
        libra_polynomials: &[Polynomial<Bn254FrParams>; NUM_SMALL_IPA_EVALUATIONS],
        sumcheck_round_univariates: &[Polynomial<Bn254FrParams>],
        sumcheck_round_evaluations: &[[Fr; 3]],
    ) -> ProverOpeningClaim<Bn254FrParams> {
        let has_zk = libra_polynomials[0].size() > 0;
        let virtual_log_n = multilinear_challenge.len();
        let log_n = (circuit_size as f64).log2() as usize;

        // Get rho challenge for batching multilinear polynomials
        let rho: Fr = transcript.get_challenge("rho");

        // Compute A_0 = F + G/X (batched polynomial)
        let a_0 = batcher.compute_batched(rho);

        // Compute Gemini fold polynomials
        let fold_polynomials =
            GeminiProver::compute_fold_polynomials(log_n, multilinear_challenge, &a_0);

        // Commit and send fold polynomial commitments (virtual_log_n - 1, matching C++)
        for (l, fold_poly) in fold_polynomials.iter().enumerate().take(virtual_log_n - 1) {
            let commitment = commitment_key.commit(fold_poly);
            let label = format!("Gemini:FOLD_{}", l + 1);
            transcript.send_to_verifier(&label, &commitment);
        }

        let r: Fr = transcript.get_challenge("Gemini:r");

        // Compute partially evaluated batch polynomials A_0+(X) and A_0-(X)
        let (a_0_pos, a_0_neg) = batcher.compute_partially_evaluated_batch_polynomials(r);

        // Construct univariate opening claims
        let gemini_claims = GeminiProver::construct_univariate_opening_claims(
            virtual_log_n,
            a_0_pos,
            a_0_neg,
            fold_polynomials,
            r,
        );

        // Send Gemini fold negative evaluations to verifier
        // claims[0] = A_0+(r), claims[1] = A_0-(-r), claims[2..] = fold evaluations
        // C++ sends claims[1..=virtual_log_n]; verifier reads virtual_log_n evaluations
        for l in 0..virtual_log_n {
            let label = format!("Gemini:a_{}", l);
            transcript.send_to_verifier(
                &label,
                &gemini_claims[l + 1].opening_pair.evaluation,
            );
        }

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

        // Create opening claims for Libra masking univariates
        let mut libra_opening_claims = if has_zk {
            let gemini_r = opening_claims[0].opening_pair.challenge;
            Self::compute_libra_opening_claims(gemini_r, libra_polynomials, transcript)
        } else {
            Vec::new()
        };

        // Create opening claims for Sumcheck Round Univariates (used in ECCVM)
        let mut sumcheck_round_claims = if !sumcheck_round_univariates.is_empty() {
            Self::compute_sumcheck_round_claims(
                circuit_size,
                multilinear_challenge,
                sumcheck_round_univariates,
                sumcheck_round_evaluations,
            )
        } else {
            Vec::new()
        };

        ShplonkProver::prove(
            commitment_key,
            &mut opening_claims,
            transcript,
            &mut libra_opening_claims,
            &mut sumcheck_round_claims,
            virtual_log_n,
        )
    }

    /// Convenience wrapper that calls `prove` without ZK or committed sumcheck data.
    pub fn prove_without_zk(
        circuit_size: usize,
        batcher: &mut PolynomialBatcher<Bn254FrParams>,
        multilinear_challenge: &[Fr],
        commitment_key: &CommitmentKey<Bn254G1Params>,
        transcript: &mut NativeTranscript,
    ) -> ProverOpeningClaim<Bn254FrParams> {
        let empty_libra: [Polynomial<Bn254FrParams>; NUM_SMALL_IPA_EVALUATIONS] =
            core::array::from_fn(|_| Polynomial::new(0, 0, 0));
        Self::prove(
            circuit_size,
            batcher,
            multilinear_challenge,
            commitment_key,
            transcript,
            &empty_libra,
            &[],
            &[],
        )
    }

    /// For ZK Flavors: Evaluate the polynomials used in SmallSubgroupIPA argument,
    /// send the evaluations to the verifier, and populate a vector of opening claims.
    ///
    /// C++ source: `ShpleminiProver_::compute_libra_opening_claims`
    pub fn compute_libra_opening_claims(
        gemini_r: Fr,
        libra_polynomials: &[Polynomial<Bn254FrParams>; NUM_SMALL_IPA_EVALUATIONS],
        transcript: &mut NativeTranscript,
    ) -> Vec<ProverOpeningClaim<Bn254FrParams>> {
        let subgroup_generator = <Bn254G1Params as SmallSubgroupCurveConfig>::subgroup_generator();

        let libra_eval_labels = [
            "Libra:concatenation_eval",
            "Libra:shifted_grand_sum_eval",
            "Libra:grand_sum_eval",
            "Libra:quotient_eval",
        ];
        let evaluation_points = [
            gemini_r,
            gemini_r * subgroup_generator,
            gemini_r,
            gemini_r,
        ];

        let mut libra_opening_claims = Vec::with_capacity(NUM_SMALL_IPA_EVALUATIONS);
        for idx in 0..NUM_SMALL_IPA_EVALUATIONS {
            let evaluation = libra_polynomials[idx].evaluate(&evaluation_points[idx]);
            transcript.send_to_verifier(libra_eval_labels[idx], &evaluation);
            libra_opening_claims.push(ProverOpeningClaim {
                polynomial: libra_polynomials[idx].clone(),
                opening_pair: OpeningPair {
                    challenge: evaluation_points[idx],
                    evaluation,
                },
                gemini_fold: false,
            });
        }

        libra_opening_claims
    }

    /// Create a vector of 3*log_n opening claims for the evaluations of Sumcheck
    /// Round Univariates at 0, 1, and a round challenge.
    ///
    /// C++ source: `ShpleminiProver_::compute_sumcheck_round_claims`
    pub fn compute_sumcheck_round_claims(
        circuit_size: usize,
        multilinear_challenge: &[Fr],
        sumcheck_round_univariates: &[Polynomial<Bn254FrParams>],
        sumcheck_round_evaluations: &[[Fr; 3]],
    ) -> Vec<ProverOpeningClaim<Bn254FrParams>> {
        let log_n = (circuit_size as f64).log2() as usize;
        let mut sumcheck_round_claims = Vec::new();

        for idx in 0..log_n {
            let evaluation_points = [Fr::zero(), Fr::one(), multilinear_challenge[idx]];
            let round_univariate = sumcheck_round_univariates[idx].clone();

            for (eval_idx, eval_point) in evaluation_points.iter().enumerate() {
                sumcheck_round_claims.push(ProverOpeningClaim {
                    polynomial: round_univariate.clone(),
                    opening_pair: OpeningPair {
                        challenge: *eval_point,
                        evaluation: sumcheck_round_evaluations[idx][eval_idx],
                    },
                    gemini_fold: false,
                });
            }
        }

        sumcheck_round_claims
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
    #[allow(clippy::too_many_arguments)]
    pub fn compute_batch_opening_claim(
        padding_indicator_array: &[Fr],
        claim_batcher: &mut ClaimBatcher<Bn254G1Params>,
        multivariate_challenge: &[Fr],
        g1_identity: &Bn254G1Affine,
        transcript: &mut NativeTranscript,
        repeated_commitments: &RepeatedCommitmentsData,
        has_zk: bool,
        libra_commitments: &[Bn254G1Affine; NUM_LIBRA_COMMITMENTS],
        libra_univariate_evaluation: Fr,
        sumcheck_round_commitments: &[Bn254G1Affine],
        sumcheck_round_evaluations: &[[Fr; 3]],
    ) -> (BatchOpeningClaim<Bn254G1Params>, bool) {
        let virtual_log_n = multivariate_challenge.len();
        let committed_sumcheck = !sumcheck_round_evaluations.is_empty();
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

        // Get interleaving evaluations (if any)
        let mut p_pos = Fr::zero();
        let mut p_neg = Fr::zero();
        if claim_batcher.interleaved.is_some() {
            p_pos = transcript.receive_from_prover("Gemini:P_pos");
            p_neg = transcript.receive_from_prover("Gemini:P_neg");
        }

        // Compute powers of r: (r, r^2, r^4, ...)
        let gemini_eval_challenge_powers =
            gemini::powers_of_evaluation_challenge(gemini_evaluation_challenge, virtual_log_n);

        // Read Libra evaluations for ZK flavors
        let mut libra_evaluations = [Fr::zero(); NUM_SMALL_IPA_EVALUATIONS];
        if has_zk {
            libra_evaluations[0] = transcript.receive_from_prover("Libra:concatenation_eval");
            libra_evaluations[1] = transcript.receive_from_prover("Libra:shifted_grand_sum_eval");
            libra_evaluations[2] = transcript.receive_from_prover("Libra:grand_sum_eval");
            libra_evaluations[3] = transcript.receive_from_prover("Libra:quotient_eval");
        }

        // Process Shplonk transcript data
        let shplonk_batching_challenge: Fr = transcript.get_challenge("Shplonk:nu");

        let shplonk_batching_challenge_powers = compute_shplonk_batching_challenge_powers(
            &shplonk_batching_challenge,
            virtual_log_n,
            has_zk,
            committed_sumcheck,
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

        // Handle interleaving
        let mut shplonk_batching_pos = Fr::zero();
        let mut shplonk_batching_neg = Fr::zero();
        if claim_batcher.interleaved.is_some() {
            let interleaved_pos_index = 2 * virtual_log_n;
            let interleaved_neg_index = interleaved_pos_index + 1;
            shplonk_batching_pos = shplonk_batching_challenge_powers[interleaved_pos_index];
            shplonk_batching_neg = shplonk_batching_challenge_powers[interleaved_neg_index];
            if let Some(ref il) = claim_batcher.interleaved {
                constant_term_accumulator = constant_term_accumulator
                    + il.shplonk_denominator
                        * (p_pos * shplonk_batching_pos + p_neg * shplonk_batching_neg);
            }
        }

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
        let a_0_pos = full_a_0_pos - p_pos;

        // Add A_0+(r)/(z-r) to constant term
        constant_term_accumulator =
            constant_term_accumulator + a_0_pos * inverse_vanishing_evals[0];
        // Add A_0-(-r)/(z+r) * nu to constant term
        constant_term_accumulator = constant_term_accumulator
            + gemini_fold_neg_evaluations[0]
                * shplonk_batching_challenge
                * inverse_vanishing_evals[1];

        // Remove repeated commitments
        Self::remove_repeated_commitments(
            &mut commitments,
            &mut scalars,
            repeated_commitments,
            has_zk,
        );

        // Add ZK data (Libra commitments and evaluations)
        let mut consistency_checked = true;
        if has_zk {
            Self::add_zk_data(
                virtual_log_n,
                &mut commitments,
                &mut scalars,
                &mut constant_term_accumulator,
                libra_commitments,
                &libra_evaluations,
                &gemini_evaluation_challenge,
                &shplonk_batching_challenge_powers,
                &shplonk_evaluation_challenge,
            );

            consistency_checked =
                SmallSubgroupIPAVerifier::<Bn254G1Params>::check_libra_evaluations_consistency(
                    &libra_evaluations,
                    &gemini_evaluation_challenge,
                    multivariate_challenge,
                    &libra_univariate_evaluation,
                );
        }

        // Add committed sumcheck data
        if committed_sumcheck {
            Self::batch_sumcheck_round_claims(
                &mut commitments,
                &mut scalars,
                &mut constant_term_accumulator,
                multivariate_challenge,
                &shplonk_batching_challenge_powers,
                &shplonk_evaluation_challenge,
                sumcheck_round_commitments,
                sumcheck_round_evaluations,
            );
        }

        // Finalize
        commitments.push(*g1_identity);
        scalars.push(constant_term_accumulator);

        (
            BatchOpeningClaim {
                commitments,
                scalars,
                evaluation_point: shplonk_evaluation_challenge,
            },
            consistency_checked,
        )
    }

    /// Convenience wrapper without ZK or committed sumcheck parameters.
    pub fn compute_batch_opening_claim_without_zk(
        padding_indicator_array: &[Fr],
        claim_batcher: &mut ClaimBatcher<Bn254G1Params>,
        multivariate_challenge: &[Fr],
        g1_identity: &Bn254G1Affine,
        transcript: &mut NativeTranscript,
    ) -> BatchOpeningClaim<Bn254G1Params> {
        let default_repeated = RepeatedCommitmentsData::default();
        let default_libra = [Bn254G1Affine::infinity(); NUM_LIBRA_COMMITMENTS];
        let (claim, _) = Self::compute_batch_opening_claim(
            padding_indicator_array,
            claim_batcher,
            multivariate_challenge,
            g1_identity,
            transcript,
            &default_repeated,
            false,
            &default_libra,
            Fr::zero(),
            &[],
            &[],
        );
        claim
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

    /// Combines scalars of repeating commitments to reduce the number of
    /// scalar multiplications performed by the verifier.
    ///
    /// C++ source: `ShpleminiVerifier_::remove_repeated_commitments`
    pub fn remove_repeated_commitments(
        commitments: &mut Vec<Bn254G1Affine>,
        scalars: &mut Vec<Fr>,
        repeated_commitments: &RepeatedCommitmentsData,
        has_zk: bool,
    ) {
        // Offset: Q commitment + optional masking poly commitment
        let offset = if has_zk { 2 } else { 1 };

        let first_range_to_be_shifted_start =
            repeated_commitments.first_range_to_be_shifted_start + offset;
        let first_range_shifted_start = repeated_commitments.first_range_shifted_start + offset;
        let first_range_size = repeated_commitments.first_range_size;

        let second_range_to_be_shifted_start =
            repeated_commitments.second_range_to_be_shifted_start + offset;
        let second_range_shifted_start = repeated_commitments.second_range_shifted_start + offset;
        let second_range_size = repeated_commitments.second_range_size;

        // Combine first range scalars
        for i in 0..first_range_size {
            let idx_to_be_shifted = i + first_range_to_be_shifted_start;
            let idx_shifted = i + first_range_shifted_start;
            scalars[idx_to_be_shifted] = scalars[idx_to_be_shifted] + scalars[idx_shifted];
        }

        // Combine second range scalars
        for i in 0..second_range_size {
            let idx_to_be_shifted = i + second_range_to_be_shifted_start;
            let idx_shifted = i + second_range_shifted_start;
            scalars[idx_to_be_shifted] = scalars[idx_to_be_shifted] + scalars[idx_shifted];
        }

        // Erase shifted entries (from higher index first to preserve indices)
        if second_range_shifted_start > first_range_shifted_start {
            for _ in 0..second_range_size {
                scalars.remove(second_range_shifted_start);
                commitments.remove(second_range_shifted_start);
            }
            for _ in 0..first_range_size {
                scalars.remove(first_range_shifted_start);
                commitments.remove(first_range_shifted_start);
            }
        } else {
            for _ in 0..first_range_size {
                scalars.remove(first_range_shifted_start);
                commitments.remove(first_range_shifted_start);
            }
            for _ in 0..second_range_size {
                scalars.remove(second_range_shifted_start);
                commitments.remove(second_range_shifted_start);
            }
        }
    }

    /// Add the opening data corresponding to Libra masking univariates to the
    /// batched opening claim.
    ///
    /// C++ source: `ShpleminiVerifier_::add_zk_data`
    #[allow(clippy::too_many_arguments)]
    pub fn add_zk_data(
        virtual_log_n: usize,
        commitments: &mut Vec<Bn254G1Affine>,
        scalars: &mut Vec<Fr>,
        constant_term_accumulator: &mut Fr,
        libra_commitments: &[Bn254G1Affine; NUM_LIBRA_COMMITMENTS],
        libra_evaluations: &[Fr; NUM_SMALL_IPA_EVALUATIONS],
        gemini_evaluation_challenge: &Fr,
        shplonk_batching_challenge_powers: &[Fr],
        shplonk_evaluation_challenge: &Fr,
    ) {
        // Add Libra commitments
        for comm in libra_commitments.iter() {
            commitments.push(*comm);
        }

        // Compute Shplonk denominators
        let subgroup_generator = <Bn254G1Params as SmallSubgroupCurveConfig>::subgroup_generator();
        let mut denominators = [Fr::zero(); NUM_SMALL_IPA_EVALUATIONS];
        denominators[0] =
            (*shplonk_evaluation_challenge - *gemini_evaluation_challenge).invert();
        denominators[1] = (*shplonk_evaluation_challenge
            - subgroup_generator * *gemini_evaluation_challenge)
            .invert();
        denominators[2] = denominators[0];
        denominators[3] = denominators[0];

        // Compute batching scalars
        let mut batching_scalars = [Fr::zero(); NUM_SMALL_IPA_EVALUATIONS];
        for idx in 0..NUM_SMALL_IPA_EVALUATIONS {
            let scaling_factor = denominators[idx]
                * shplonk_batching_challenge_powers
                    [2 * virtual_log_n + NUM_INTERLEAVING_CLAIMS + idx];
            batching_scalars[idx] = -scaling_factor;
            *constant_term_accumulator =
                *constant_term_accumulator + scaling_factor * libra_evaluations[idx];
        }

        // To save a scalar mul, add the sum of batching scalars for the big sum evaluations
        scalars.push(batching_scalars[0]);
        scalars.push(batching_scalars[1] + batching_scalars[2]);
        scalars.push(batching_scalars[3]);
    }

    /// Adds the Sumcheck data into the Shplemini BatchOpeningClaim.
    ///
    /// C++ source: `ShpleminiVerifier_::batch_sumcheck_round_claims`
    #[allow(clippy::too_many_arguments)]
    pub fn batch_sumcheck_round_claims(
        commitments: &mut Vec<Bn254G1Affine>,
        scalars: &mut Vec<Fr>,
        constant_term_accumulator: &mut Fr,
        multilinear_challenge: &[Fr],
        shplonk_batching_challenge_powers: &[Fr],
        shplonk_evaluation_challenge: &Fr,
        sumcheck_round_commitments: &[Bn254G1Affine],
        sumcheck_round_evaluations: &[[Fr; 3]],
    ) {
        let num_gemini_claims = 2 * multilinear_challenge.len();

        // Denominators for evaluation claims at 0 and 1
        let const_denominators = [
            shplonk_evaluation_challenge.invert(),
            (*shplonk_evaluation_challenge - Fr::one()).invert(),
        ];

        // Compute denominators for round challenges and add commitments
        let mut denominators: Vec<Fr> = Vec::with_capacity(multilinear_challenge.len());
        for (challenge, comm) in multilinear_challenge
            .iter()
            .zip(sumcheck_round_commitments.iter())
        {
            denominators.push(*shplonk_evaluation_challenge - *challenge);
            commitments.push(*comm);
        }

        // Batch invert denominators
        Field::batch_invert(&mut denominators);

        // Compute scalars for each round
        let mut power = num_gemini_claims + NUM_INTERLEAVING_CLAIMS + NUM_SMALL_IPA_EVALUATIONS;
        for (eval_array, denominator) in sumcheck_round_evaluations.iter().zip(denominators.iter())
        {
            let mut batched_scalar = Fr::zero();
            let mut const_term_contribution = Fr::zero();

            // Contributions from evaluations at 0 and 1
            for idx in 0..2 {
                let current_scaling_factor =
                    const_denominators[idx] * shplonk_batching_challenge_powers[power];
                power += 1;
                batched_scalar = batched_scalar - current_scaling_factor;
                const_term_contribution =
                    const_term_contribution + current_scaling_factor * eval_array[idx];
            }

            // Contribution from evaluation at challenge u_i
            let current_scaling_factor =
                *denominator * shplonk_batching_challenge_powers[power];
            power += 1;
            batched_scalar = batched_scalar - current_scaling_factor;
            const_term_contribution =
                const_term_contribution + current_scaling_factor * eval_array[2];

            *constant_term_accumulator = *constant_term_accumulator + const_term_contribution;
            scalars.push(batched_scalar);
        }
    }
}
