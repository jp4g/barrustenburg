//! Port of `sumcheck.hpp` — SumcheckProver and SumcheckVerifier.
//!
//! Implements the full sumcheck protocol orchestration for SumcheckTestFlavor (non-ZK).

use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};
use bbrs_polynomials::gate_separator::GateSeparatorPolynomial;
use bbrs_polynomials::univariate::Univariate;
use bbrs_relations::relation_parameters::RelationParameters;
use bbrs_transcript::codec::FieldSerializable;
use bbrs_transcript::transcript::BaseTranscript;

use bbrs_flavor::sumcheck_test_flavor::{
    AllEntities, AllValues, PartiallyEvaluatedMultivariates, ProverPolynomials,
    BATCHED_RELATION_PARTIAL_LENGTH, NUM_ALL_ENTITIES,
};

use crate::sumcheck_output::SumcheckOutput;
use crate::sumcheck_round::{
    initialize_relation_separator, SumcheckProverRound, SumcheckVerifierRound,
    SubrelationSeparators,
};

// ============================================================================
// FieldSerializable for Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH>
// ============================================================================

/// Wrapper to send/receive Univariates via the transcript.
///
/// The C++ transcript serializes Univariate evaluations as N consecutive Fr elements.
struct UnivariateFrWrapper(Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH>);

impl FieldSerializable for UnivariateFrWrapper {
    const NUM_FR: usize = BATCHED_RELATION_PARTIAL_LENGTH;

    fn serialize_to_frs(&self) -> Vec<Fr> {
        self.0.evaluations.to_vec()
    }

    fn deserialize_from_frs(frs: &[Fr]) -> Self {
        assert_eq!(frs.len(), BATCHED_RELATION_PARTIAL_LENGTH);
        let mut evals = [Fr::zero(); BATCHED_RELATION_PARTIAL_LENGTH];
        evals.copy_from_slice(frs);
        Self(Univariate::new(evals))
    }
}

// ============================================================================
// SumcheckProver
// ============================================================================

/// Sumcheck prover for SumcheckTestFlavor (non-ZK).
///
/// Port of C++ `SumcheckProver<Flavor>` for the non-ZK path.
pub struct SumcheckProver<'a, H: bbrs_transcript::transcript::TranscriptHasher> {
    /// The size of the hypercube, 2^d.
    pub multivariate_n: usize,
    /// The number of variables d.
    pub multivariate_d: usize,
    /// Reference to the prover polynomials.
    pub full_polynomials: &'a ProverPolynomials<Bn254FrParams>,
    /// Transcript for Fiat-Shamir.
    pub transcript: &'a mut BaseTranscript<H>,
    /// Round computation state.
    pub round: SumcheckProverRound<Bn254FrParams>,
    /// Subrelation separator challenges.
    pub alphas: SubrelationSeparators<Bn254FrParams>,
    /// Gate challenges for pow-polynomial.
    pub gate_challenges: Vec<Fr>,
    /// Relation parameters.
    pub relation_parameters: RelationParameters<Fr>,
}

impl<'a, H: bbrs_transcript::transcript::TranscriptHasher> SumcheckProver<'a, H> {
    /// Create a new SumcheckProver.
    pub fn new(
        multivariate_n: usize,
        full_polynomials: &'a ProverPolynomials<Bn254FrParams>,
        transcript: &'a mut BaseTranscript<H>,
        alpha: Fr,
        gate_challenges: Vec<Fr>,
        relation_parameters: RelationParameters<Fr>,
    ) -> Self {
        let multivariate_d = log2(multivariate_n);
        Self {
            multivariate_n,
            multivariate_d,
            full_polynomials,
            transcript,
            round: SumcheckProverRound::new(multivariate_n),
            alphas: initialize_relation_separator(alpha),
            gate_challenges,
            relation_parameters,
        }
    }

    /// Run the sumcheck protocol (non-ZK).
    ///
    /// Returns the challenge vector and claimed evaluations.
    ///
    /// Port of C++ `SumcheckProver::prove()` (non-ZK version).
    pub fn prove(mut self) -> SumcheckOutput<Bn254FrParams, AllValues<Bn254FrParams>> {
        let gate_separators =
            GateSeparatorPolynomial::new(self.gate_challenges.clone(), self.multivariate_d);

        let mut multivariate_challenge = Vec::with_capacity(self.multivariate_d);

        // Round 0: compute univariate from full polynomials
        let round_univariate = self.round.compute_univariate(
            self.full_polynomials,
            &self.relation_parameters,
            &gate_separators,
            &self.alphas,
        );

        // Initialize partially evaluated polynomials
        let mut partially_evaluated_polynomials =
            PartiallyEvaluatedMultivariates::from_prover_polynomials(
                self.full_polynomials,
                self.multivariate_n,
            );

        // Send round 0 univariate, get challenge, partially evaluate
        self.transcript.send_to_verifier(
            "Sumcheck:univariate_0",
            &UnivariateFrWrapper(round_univariate),
        );
        let round_challenge = self.transcript.get_challenge("Sumcheck:u_0");
        multivariate_challenge.push(round_challenge);

        Self::partially_evaluate_from_full(
            &mut partially_evaluated_polynomials,
            self.full_polynomials,
            round_challenge,
        );
        let mut gate_separators = gate_separators;
        gate_separators.partially_evaluate(round_challenge);
        self.round.round_size >>= 1;

        // Rounds 1 through d-1
        for round_idx in 1..self.multivariate_d {
            let round_univariate = self.round.compute_univariate(
                &partially_evaluated_polynomials,
                &self.relation_parameters,
                &gate_separators,
                &self.alphas,
            );

            self.transcript.send_to_verifier(
                &format!("Sumcheck:univariate_{}", round_idx),
                &UnivariateFrWrapper(round_univariate),
            );
            let round_challenge = self
                .transcript
                .get_challenge(&format!("Sumcheck:u_{}", round_idx));
            multivariate_challenge.push(round_challenge);

            Self::partially_evaluate_in_place(
                &mut partially_evaluated_polynomials,
                round_challenge,
            );
            gate_separators.partially_evaluate(round_challenge);
            self.round.round_size >>= 1;
        }

        // Extract claimed evaluations from the top row
        let multivariate_evaluations =
            Self::extract_claimed_evaluations(&partially_evaluated_polynomials);

        // Send evaluations to verifier
        let eval_array: [Fr; NUM_ALL_ENTITIES] = {
            let all = multivariate_evaluations.get_all();
            let mut arr = [Fr::zero(); NUM_ALL_ENTITIES];
            for (i, val) in all.iter().enumerate() {
                arr[i] = **val;
            }
            arr
        };
        self.transcript
            .send_to_verifier("Sumcheck:evaluations", &eval_array);

        SumcheckOutput {
            challenge: multivariate_challenge,
            claimed_evaluations: multivariate_evaluations,
            verified: false,
        }
    }

    /// Partially evaluate: fold full_polynomials into partially_evaluated_polynomials.
    ///
    /// pep[i/2] = poly[2i] + challenge * (poly[2i+1] - poly[2i])
    fn partially_evaluate_from_full(
        pep: &mut PartiallyEvaluatedMultivariates<Bn254FrParams>,
        full_polynomials: &ProverPolynomials<Bn254FrParams>,
        round_challenge: Fr,
    ) {
        let pep_all = pep.get_all_mut();
        let poly_all = full_polynomials.get_all();

        for (pep_poly, src_poly) in pep_all.into_iter().zip(poly_all.iter()) {
            let limit = src_poly.end_index();
            let mut i = 0;
            while i < limit {
                let val_0 = src_poly.get(i);
                let val_1 = src_poly.get(i + 1);
                *pep_poly.at_mut(i >> 1) = val_0 + round_challenge * (val_1 - val_0);
                i += 2;
            }
            let new_end = (limit / 2) + (limit % 2);
            pep_poly.shrink_end_index(new_end);
        }
    }

    /// Partially evaluate: fold partially_evaluated_polynomials in place.
    fn partially_evaluate_in_place(
        pep: &mut PartiallyEvaluatedMultivariates<Bn254FrParams>,
        round_challenge: Fr,
    ) {
        let pep_all = pep.get_all_mut();
        for pep_poly in pep_all.into_iter() {
            let limit = pep_poly.end_index();
            let mut i = 0;
            while i < limit {
                let val_0 = pep_poly.get(i);
                let val_1 = pep_poly.get(i + 1);
                *pep_poly.at_mut(i >> 1) = val_0 + round_challenge * (val_1 - val_0);
                i += 2;
            }
            let new_end = (limit / 2) + (limit % 2);
            pep_poly.shrink_end_index(new_end);
        }
    }

    /// Extract claimed evaluations from the top row of partially evaluated polynomials.
    fn extract_claimed_evaluations(
        pep: &PartiallyEvaluatedMultivariates<Bn254FrParams>,
    ) -> AllValues<Bn254FrParams> {
        let all = pep.get_all();
        AllEntities {
            q_m: all[0].get(0),
            q_l: all[1].get(0),
            q_r: all[2].get(0),
            q_o: all[3].get(0),
            q_4: all[4].get(0),
            q_c: all[5].get(0),
            q_arith: all[6].get(0),
            q_test: all[7].get(0),
            w_l: all[8].get(0),
            w_r: all[9].get(0),
            w_o: all[10].get(0),
            w_4: all[11].get(0),
            w_test_1: all[12].get(0),
            w_test_2: all[13].get(0),
            w_l_shift: all[14].get(0),
            w_4_shift: all[15].get(0),
        }
    }
}

// ============================================================================
// SumcheckVerifier
// ============================================================================

/// Sumcheck verifier for SumcheckTestFlavor (non-ZK).
///
/// Port of C++ `SumcheckVerifier<Flavor>`.
pub struct SumcheckVerifier<'a, H: bbrs_transcript::transcript::TranscriptHasher> {
    /// Transcript for Fiat-Shamir.
    pub transcript: &'a mut BaseTranscript<H>,
    /// Verifier round state.
    pub round: SumcheckVerifierRound<Bn254FrParams>,
    /// Subrelation separator challenges.
    pub alphas: SubrelationSeparators<Bn254FrParams>,
    /// Number of sumcheck rounds.
    pub multivariate_d: usize,
}

impl<'a, H: bbrs_transcript::transcript::TranscriptHasher> SumcheckVerifier<'a, H> {
    /// Create a new SumcheckVerifier.
    pub fn new(
        transcript: &'a mut BaseTranscript<H>,
        alpha: Fr,
        multivariate_d: usize,
        target_sum: Fr,
    ) -> Self {
        Self {
            transcript,
            round: SumcheckVerifierRound::new(target_sum),
            alphas: initialize_relation_separator(alpha),
            multivariate_d,
        }
    }

    /// Run the sumcheck verification protocol.
    ///
    /// Port of C++ `SumcheckVerifier::verify()` (non-ZK, non-Grumpkin).
    pub fn verify(
        mut self,
        relation_parameters: &RelationParameters<Fr>,
        gate_challenges: &[Fr],
        padding_indicator_array: &[Fr],
    ) -> SumcheckOutput<Bn254FrParams, AllValues<Bn254FrParams>> {
        let mut gate_separators = GateSeparatorPolynomial::new_verifier(gate_challenges.to_vec());

        let mut multivariate_challenge = Vec::with_capacity(self.multivariate_d);
        let mut verified = true;

        // Process each round
        for round_idx in 0..self.multivariate_d {
            // Receive round univariate from transcript
            let round_univariate: UnivariateFrWrapper = self
                .transcript
                .receive_from_prover(&format!("Sumcheck:univariate_{}", round_idx));
            let round_univariate = round_univariate.0;

            let indicator = padding_indicator_array[round_idx];

            // Check sum
            self.round.check_sum(&round_univariate, indicator);
            verified = verified && !self.round.round_failed;

            // Get challenge
            let round_challenge = self
                .transcript
                .get_challenge(&format!("Sumcheck:u_{}", round_idx));
            multivariate_challenge.push(round_challenge);

            // Compute next target sum
            self.round
                .compute_next_target_sum(&round_univariate, round_challenge, indicator);

            gate_separators.partially_evaluate(round_challenge);
        }

        // Receive claimed evaluations
        let transcript_evaluations: [Fr; NUM_ALL_ENTITIES] = self
            .transcript
            .receive_from_prover("Sumcheck:evaluations");

        let purported_evaluations = AllEntities {
            q_m: transcript_evaluations[0],
            q_l: transcript_evaluations[1],
            q_r: transcript_evaluations[2],
            q_o: transcript_evaluations[3],
            q_4: transcript_evaluations[4],
            q_c: transcript_evaluations[5],
            q_arith: transcript_evaluations[6],
            q_test: transcript_evaluations[7],
            w_l: transcript_evaluations[8],
            w_r: transcript_evaluations[9],
            w_o: transcript_evaluations[10],
            w_4: transcript_evaluations[11],
            w_test_1: transcript_evaluations[12],
            w_test_2: transcript_evaluations[13],
            w_l_shift: transcript_evaluations[14],
            w_4_shift: transcript_evaluations[15],
        };

        // Compute full relation purported value
        let full_honk_purported_value = self.round.compute_full_relation_purported_value(
            &purported_evaluations,
            relation_parameters,
            &gate_separators,
            &self.alphas,
        );

        // Final verification
        verified = self
            .round
            .perform_final_verification(full_honk_purported_value)
            && verified;

        SumcheckOutput {
            challenge: multivariate_challenge,
            claimed_evaluations: purported_evaluations,
            verified,
        }
    }
}

// ============================================================================
// Helpers
// ============================================================================

/// Compute log2 of a power of 2.
fn log2(n: usize) -> usize {
    assert!(n.is_power_of_two(), "n must be a power of 2, got {}", n);
    n.trailing_zeros() as usize
}
