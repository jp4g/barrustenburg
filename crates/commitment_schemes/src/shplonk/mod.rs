//! SHPLONK: batched polynomial opening via quotient polynomial construction.
//!
//! C++ source: barretenberg/commitment_schemes/shplonk/shplonk.hpp
//!
//! Reduces multiple claims about commitments, each opened at a single point,
//! into a single claim for a single polynomial opened at a single point.
//!
//! Terminology:
//! - B_k(X) is a random linear combination of all polynomials opened at Omega_k
//! - T_k(X) interpolates B_k(X) over Omega_k
//! - z_k(X) = product of (X - x) for x in Omega_k
//! - nu is the batching challenge, z is the random evaluation challenge

pub mod shplemini;

use bbrs_ecc::curves::bn254::{Bn254FrParams, Bn254G1Params, Fr, G1Affine as Bn254G1Affine};
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;
use bbrs_ecc::groups::element::Element;
use bbrs_numeric::bitop::round_up_power_2;
use bbrs_polynomials::polynomial::Polynomial;
use bbrs_transcript::NativeTranscript;

use crate::claim::{BatchOpeningClaim, OpeningClaim, OpeningPair};
use crate::commitment_key::CommitmentKey;

// ── Constants ────────────────────────────────────────────────────────────────

/// Number of extra interleaving claims (used by Translator).
///
/// C++ source: `NUM_INTERLEAVING_CLAIMS` in constants.hpp.
pub const NUM_INTERLEAVING_CLAIMS: usize = 2;

/// Number of SmallSubgroupIPA evaluations sent by the prover (ZK flavors).
///
/// C++ source: `NUM_SMALL_IPA_EVALUATIONS` in constants.hpp.
pub const NUM_SMALL_IPA_EVALUATIONS: usize = 4;

/// Number of Libra commitments (ZK flavors).
///
/// C++ source: `NUM_LIBRA_COMMITMENTS` in constants.hpp.
pub const NUM_LIBRA_COMMITMENTS: usize = 3;

// ── ProverOpeningClaim (field-level, used internally) ────────────────────────

/// Polynomial p and an opening pair (r, v) such that p(r) = v.
///
/// This is the field-level version used internally by Shplonk.
/// Structurally identical to `crate::claim::ProverOpeningClaim` but parameterized
/// on FieldParams rather than CurveParams.
pub struct ProverOpeningClaim<P: FieldParams> {
    pub polynomial: Polynomial<P>,
    pub opening_pair: OpeningPair<P>,
    pub gemini_fold: bool,
}

// ── ShplonkProver ────────────────────────────────────────────────────────────

/// Shplonk prover: batches multiple opening claims into a single claim.
///
/// C++ source: `ShplonkProver_<Curve>`
pub struct ShplonkProver;

impl ShplonkProver {
    /// Compute batched quotient polynomial Q(X) = sum_j nu^j * (f_j(X) - v_j) / (X - x_j).
    ///
    /// C++ source: `ShplonkProver_::compute_batched_quotient`
    pub fn compute_batched_quotient<P: FieldParams>(
        virtual_log_n: usize,
        opening_claims: &[ProverOpeningClaim<P>],
        nu: &Field<P>,
        gemini_fold_pos_evaluations: &mut [Field<P>],
        libra_opening_claims: &[ProverOpeningClaim<P>],
        sumcheck_round_claims: &[ProverOpeningClaim<P>],
    ) -> Polynomial<P> {
        // Find maximum polynomial size among all claims
        let mut max_poly_size: usize = 0;
        for claim in opening_claims
            .iter()
            .chain(libra_opening_claims.iter())
            .chain(sumcheck_round_claims.iter())
        {
            max_poly_size = max_poly_size.max(claim.polynomial.size());
        }
        // Round up to next power of 2
        max_poly_size = round_up_power_2(max_poly_size as u64) as usize;

        // Q(X) = sum_j nu^j * (f_j(X) - v_j) / (X - x_j)
        let mut q = Polynomial::<P>::new(max_poly_size, max_poly_size, 0);

        let mut current_nu = Field::<P>::one();
        let mut fold_idx = 0;

        for claim in opening_claims {
            // Gemini Fold Polynomials opened at -r^{2^j} and r^{2^j}
            if claim.gemini_fold {
                let mut tmp = claim.polynomial.clone();
                *tmp.at_mut(0) = tmp.get(0) - gemini_fold_pos_evaluations[fold_idx];
                fold_idx += 1;
                tmp.factor_roots(&(-claim.opening_pair.challenge));
                q.add_scaled(&tmp.as_span(), current_nu);
                current_nu = current_nu * *nu;
            }

            // Compute individual claim quotient: (f_j(X) - v_j) / (X - x_j)
            let mut tmp = claim.polynomial.clone();
            *tmp.at_mut(0) = tmp.get(0) - claim.opening_pair.evaluation;
            tmp.factor_roots(&claim.opening_pair.challenge);
            q.add_scaled(&tmp.as_span(), current_nu);
            current_nu = current_nu * *nu;
        }

        // Libra opening claims use a fixed batching power
        if !libra_opening_claims.is_empty() {
            current_nu = nu.pow(&[(2 * virtual_log_n + NUM_INTERLEAVING_CLAIMS) as u64, 0, 0, 0]);
        }

        for claim in libra_opening_claims {
            let mut tmp = claim.polynomial.clone();
            *tmp.at_mut(0) = tmp.get(0) - claim.opening_pair.evaluation;
            tmp.factor_roots(&claim.opening_pair.challenge);
            q.add_scaled(&tmp.as_span(), current_nu);
            current_nu = current_nu * *nu;
        }

        for claim in sumcheck_round_claims {
            let mut tmp = claim.polynomial.clone();
            *tmp.at_mut(0) = tmp.get(0) - claim.opening_pair.evaluation;
            tmp.factor_roots(&claim.opening_pair.challenge);
            q.add_scaled(&tmp.as_span(), current_nu);
            current_nu = current_nu * *nu;
        }

        q
    }

    /// Compute partially evaluated batched quotient polynomial difference
    /// G(X) = Q(X) - Q_z(X), such that G(z) = 0.
    ///
    /// C++ source: `ShplonkProver_::compute_partially_evaluated_batched_quotient`
    pub fn compute_partially_evaluated_batched_quotient<P: FieldParams>(
        virtual_log_n: usize,
        opening_claims: &mut [ProverOpeningClaim<P>],
        batched_quotient_q: Polynomial<P>,
        nu_challenge: &Field<P>,
        z_challenge: &Field<P>,
        gemini_fold_pos_evaluations: &[Field<P>],
        libra_opening_claims: &mut [ProverOpeningClaim<P>],
        sumcheck_opening_claims: &mut [ProverOpeningClaim<P>],
    ) -> ProverOpeningClaim<P> {
        // Count total opening claims (Gemini folds produce 2 claims each)
        let num_gemini_opening_claims = 2 * opening_claims.len();
        let num_opening_claims =
            num_gemini_opening_claims + libra_opening_claims.len() + sumcheck_opening_claims.len();

        // Compute inverse vanishing evaluations: 1 / (z - x_j)
        let mut inverse_vanishing_evals: Vec<Field<P>> = Vec::with_capacity(num_opening_claims);

        for claim in opening_claims.iter() {
            if claim.gemini_fold {
                // 1 / (z + r^{2^j}) for the positive evaluation point
                inverse_vanishing_evals.push(*z_challenge + claim.opening_pair.challenge);
            }
            // 1 / (z - x_j) for the negative evaluation point
            inverse_vanishing_evals.push(*z_challenge - claim.opening_pair.challenge);
        }

        for claim in libra_opening_claims.iter() {
            inverse_vanishing_evals.push(*z_challenge - claim.opening_pair.challenge);
        }

        for claim in sumcheck_opening_claims.iter() {
            inverse_vanishing_evals.push(*z_challenge - claim.opening_pair.challenge);
        }

        Field::batch_invert(&mut inverse_vanishing_evals);

        // G(X) = Q(X) - sum_j nu^j * (f_j(X) - v_j) / (z - x_j)
        let mut g = batched_quotient_q;
        let mut current_nu = Field::<P>::one();
        let mut idx = 0;
        let mut fold_idx = 0;

        for claim in opening_claims.iter_mut() {
            if claim.gemini_fold {
                let mut tmp = claim.polynomial.clone();
                *tmp.at_mut(0) = tmp.get(0) - gemini_fold_pos_evaluations[fold_idx];
                fold_idx += 1;
                let scaling_factor = current_nu * inverse_vanishing_evals[idx];
                idx += 1;
                g.add_scaled(&tmp.as_span(), -scaling_factor);
                current_nu = current_nu * *nu_challenge;
            }

            *claim.polynomial.at_mut(0) =
                claim.polynomial.get(0) - claim.opening_pair.evaluation;
            let scaling_factor = current_nu * inverse_vanishing_evals[idx];
            idx += 1;
            g.add_scaled(&claim.polynomial.as_span(), -scaling_factor);
            current_nu = current_nu * *nu_challenge;
        }

        // Libra claims
        if !libra_opening_claims.is_empty() {
            current_nu = nu_challenge.pow(&[(2 * virtual_log_n + NUM_INTERLEAVING_CLAIMS) as u64, 0, 0, 0]);
        }

        for claim in libra_opening_claims.iter_mut() {
            *claim.polynomial.at_mut(0) =
                claim.polynomial.get(0) - claim.opening_pair.evaluation;
            let scaling_factor = current_nu * inverse_vanishing_evals[idx];
            idx += 1;
            g.add_scaled(&claim.polynomial.as_span(), -scaling_factor);
            current_nu = current_nu * *nu_challenge;
        }

        for claim in sumcheck_opening_claims.iter_mut() {
            *claim.polynomial.at_mut(0) =
                claim.polynomial.get(0) - claim.opening_pair.evaluation;
            let scaling_factor = current_nu * inverse_vanishing_evals[idx];
            idx += 1;
            g.add_scaled(&claim.polynomial.as_span(), -scaling_factor);
            current_nu = current_nu * *nu_challenge;
        }

        ProverOpeningClaim {
            polynomial: g,
            opening_pair: OpeningPair {
                challenge: *z_challenge,
                evaluation: Field::zero(),
            },
            gemini_fold: false,
        }
    }

    /// Compute evaluations of fold polynomials Fold_i at r^{2^i} for i > 0.
    ///
    /// C++ source: `ShplonkProver_::compute_gemini_fold_pos_evaluations`
    pub fn compute_gemini_fold_pos_evaluations<P: FieldParams>(
        opening_claims: &[ProverOpeningClaim<P>],
    ) -> Vec<Field<P>> {
        let mut evals = Vec::new();
        for claim in opening_claims {
            if claim.gemini_fold {
                let eval_point = -claim.opening_pair.challenge;
                let evaluation = claim.polynomial.evaluate(&eval_point);
                evals.push(evaluation);
            }
        }
        evals
    }

    /// Returns a batched opening claim equivalent to a set of opening claims
    /// consisting of polynomials, each opened at a single point.
    ///
    /// C++ source: `ShplonkProver_::prove`
    ///
    /// Specialized for BN254 because NativeTranscript methods are BN254-specific.
    pub fn prove(
        commitment_key: &CommitmentKey<Bn254G1Params>,
        opening_claims: &mut [ProverOpeningClaim<Bn254FrParams>],
        transcript: &mut NativeTranscript,
        libra_opening_claims: &mut [ProverOpeningClaim<Bn254FrParams>],
        sumcheck_round_claims: &mut [ProverOpeningClaim<Bn254FrParams>],
        virtual_log_n: usize,
    ) -> ProverOpeningClaim<Bn254FrParams> {
        let nu: Fr = transcript.get_challenge("Shplonk:nu");

        let mut gemini_fold_pos_evaluations =
            Self::compute_gemini_fold_pos_evaluations(opening_claims);

        let batched_quotient = Self::compute_batched_quotient(
            virtual_log_n,
            opening_claims,
            &nu,
            &mut gemini_fold_pos_evaluations,
            libra_opening_claims,
            sumcheck_round_claims,
        );

        let batched_quotient_commitment = commitment_key.commit(&batched_quotient);
        transcript.send_to_verifier("Shplonk:Q", &batched_quotient_commitment);

        let z: Fr = transcript.get_challenge("Shplonk:z");

        Self::compute_partially_evaluated_batched_quotient(
            virtual_log_n,
            opening_claims,
            batched_quotient,
            &nu,
            &z,
            &gemini_fold_pos_evaluations,
            libra_opening_claims,
            sumcheck_round_claims,
        )
    }
}

// ── ShplonkVerifier ──────────────────────────────────────────────────────────

/// A claim constructed as a linear combination of commitments.
///
/// Used to update the internal state of the Shplonk verifier.
/// The state is updated to add the check: sum_l a_{j_l} * f_{j_l}(x) = v.
///
/// C++ source: `ShplonkVerifier_::LinearCombinationOfClaims`
pub struct LinearCombinationOfClaims<P: FieldParams> {
    pub indices: Vec<usize>,
    pub scalars: Vec<Field<P>>,
    pub opening_pair: OpeningPair<P>,
}

/// Shplonk verifier: reduces multiple opening claims to a single KZG/IPA claim.
///
/// C++ source: `ShplonkVerifier_<Curve>`
pub struct ShplonkVerifier<C: CurveParams> {
    pows_of_nu: Vec<Field<C::ScalarFieldParams>>,
    pow_idx: usize,
    quotient: AffineElement<C>,
    z_challenge: Field<C::ScalarFieldParams>,
    commitments: Vec<AffineElement<C>>,
    scalars: Vec<Field<C::ScalarFieldParams>>,
    identity_scalar_coefficient: Field<C::ScalarFieldParams>,
    evaluation: Field<C::ScalarFieldParams>,
}

impl<C: CurveParams> ShplonkVerifier<C> {
    /// Update the internal state with a new claim.
    ///
    /// C++ source: `ShplonkVerifier_::update`
    pub fn update(
        &mut self,
        update_data: &LinearCombinationOfClaims<C::ScalarFieldParams>,
        inverse_vanishing_eval: &Field<C::ScalarFieldParams>,
    ) {
        let scalar_factor = self.pows_of_nu[self.pow_idx] * *inverse_vanishing_eval;

        for (index, coefficient) in update_data.indices.iter().zip(update_data.scalars.iter()) {
            let scaling_factor = scalar_factor * *coefficient;
            self.scalars[index + 1] = self.scalars[index + 1] - scaling_factor;
        }

        self.identity_scalar_coefficient =
            self.identity_scalar_coefficient + scalar_factor * update_data.opening_pair.evaluation;

        self.pow_idx += 1;
    }

    /// Finalize the Shplonk verification and return the opening claim.
    ///
    /// Computes: [Q] - sum_i s_i * [f_i] + theta * [1]
    ///
    /// C++ source: `ShplonkVerifier_::finalize`
    pub fn finalize(
        &mut self,
        g1_identity: &AffineElement<C>,
    ) -> OpeningClaim<C> {
        self.commitments.push(*g1_identity);
        self.scalars.push(self.identity_scalar_coefficient);

        let mut result = Element::<C>::infinity();
        for (commitment, scalar) in self.commitments.iter().zip(self.scalars.iter()) {
            result = result + Element::<C>::from_affine(commitment).mul(scalar);
        }

        OpeningClaim {
            opening_pair: OpeningPair {
                challenge: self.z_challenge,
                evaluation: self.evaluation,
            },
            commitment: result.to_affine(),
        }
    }

    /// Export a BatchOpeningClaim instead of performing final batch_mul.
    ///
    /// C++ source: `ShplonkVerifier_::export_batch_opening_claim`
    pub fn export_batch_opening_claim(
        &mut self,
        g1_identity: &AffineElement<C>,
    ) -> BatchOpeningClaim<C> {
        self.commitments.push(*g1_identity);
        self.scalars.push(self.identity_scalar_coefficient);

        BatchOpeningClaim {
            commitments: self.commitments.clone(),
            scalars: self.scalars.clone(),
            evaluation_point: self.z_challenge,
        }
    }

    /// Update state with vector claims (no finalize).
    ///
    /// C++ source: `ShplonkVerifier_::reduce_verification_vector_claims_no_finalize`
    pub fn reduce_verification_vector_claims_no_finalize(
        &mut self,
        claims: &[LinearCombinationOfClaims<C::ScalarFieldParams>],
    ) {
        let mut inverse_vanishing_evals: Vec<Field<C::ScalarFieldParams>> = claims
            .iter()
            .map(|c| self.z_challenge - c.opening_pair.challenge)
            .collect();
        Field::batch_invert(&mut inverse_vanishing_evals);

        for (claim, inv) in claims.iter().zip(inverse_vanishing_evals.iter()) {
            self.update(claim, inv);
        }
    }

    /// Reduce verification with vector claims and finalize.
    ///
    /// C++ source: `ShplonkVerifier_::reduce_verification_vector_claims`
    pub fn reduce_verification_vector_claims(
        &mut self,
        g1_identity: &AffineElement<C>,
        claims: &[LinearCombinationOfClaims<C::ScalarFieldParams>],
    ) -> OpeningClaim<C> {
        self.reduce_verification_vector_claims_no_finalize(claims);
        self.finalize(g1_identity)
    }

    /// Compute inverted Gemini denominators:
    /// 1/(z-r), 1/(z+r), ..., 1/(z-r^{2^{d-1}}), 1/(z+r^{2^{d-1}})
    ///
    /// C++ source: `ShplonkVerifier_::compute_inverted_gemini_denominators`
    pub fn compute_inverted_gemini_denominators(
        shplonk_eval_challenge: &Field<C::ScalarFieldParams>,
        gemini_eval_challenge_powers: &[Field<C::ScalarFieldParams>],
    ) -> Vec<Field<C::ScalarFieldParams>> {
        let virtual_log_n = gemini_eval_challenge_powers.len();
        let num_gemini_claims = 2 * virtual_log_n;
        let mut denominators = Vec::with_capacity(num_gemini_claims);

        for power in gemini_eval_challenge_powers {
            // 1/(z - r^{2^j})
            denominators.push(*shplonk_eval_challenge - *power);
            // 1/(z + r^{2^j})
            denominators.push(*shplonk_eval_challenge + *power);
        }

        Field::batch_invert(&mut denominators);
        denominators
    }
}

// ── ShplonkVerifier BN254-specific methods (transcript-dependent) ────────────

impl ShplonkVerifier<Bn254G1Params> {
    /// Construct a new Shplonk verifier from polynomial commitments and transcript.
    ///
    /// C++ source: `ShplonkVerifier_::ShplonkVerifier_()`
    ///
    /// Specialized for BN254 because NativeTranscript methods are BN254-specific.
    pub fn new(
        polynomial_commitments: &mut Vec<Bn254G1Affine>,
        transcript: &mut NativeTranscript,
        num_claims: usize,
    ) -> Self {
        assert!(
            num_claims > 1,
            "Using Shplonk with just one claim. Should use batch reduction."
        );

        let nu: Fr = transcript.get_challenge("Shplonk:nu");
        let quotient: Bn254G1Affine = transcript.receive_from_prover("Shplonk:Q");
        let z_challenge: Fr = transcript.get_challenge("Shplonk:z");

        // Start with [Q] as first commitment, scalar = 1
        let mut commitments = vec![quotient];
        let mut scalars = vec![Fr::one()];

        // Precompute powers of nu
        let mut pows_of_nu = vec![Fr::one(), nu];
        for _ in 2..num_claims {
            pows_of_nu.push(*pows_of_nu.last().unwrap() * nu);
        }

        // Add polynomial commitments
        commitments.append(&mut polynomial_commitments.clone());
        // Initialize corresponding scalars to zero
        let num_poly_commitments = commitments.len() - 1;
        scalars.resize(1 + num_poly_commitments, Fr::zero());

        Self {
            pows_of_nu,
            pow_idx: 0,
            quotient,
            z_challenge,
            commitments,
            scalars,
            identity_scalar_coefficient: Fr::zero(),
            evaluation: Fr::zero(),
        }
    }

    /// Instantiate a Shplonk verifier and update its state with the provided claims.
    ///
    /// C++ source: `ShplonkVerifier_::reduce_verification_no_finalize`
    pub fn reduce_verification_no_finalize(
        claims: &[OpeningClaim<Bn254G1Params>],
        transcript: &mut NativeTranscript,
    ) -> Self {
        let num_claims = claims.len();
        let mut polynomial_commitments: Vec<Bn254G1Affine> =
            claims.iter().map(|c| c.commitment).collect();

        let mut verifier = Self::new(&mut polynomial_commitments, transcript, num_claims);

        // Compute { 1 / (z - x_i) }
        let mut inverse_vanishing_evals: Vec<Fr> = claims
            .iter()
            .map(|c| verifier.z_challenge - c.opening_pair.challenge)
            .collect();
        Fr::batch_invert(&mut inverse_vanishing_evals);

        for (idx, claim) in claims.iter().enumerate() {
            let lcc = LinearCombinationOfClaims {
                indices: vec![idx],
                scalars: vec![Fr::one()],
                opening_pair: OpeningPair {
                    challenge: claim.opening_pair.challenge,
                    evaluation: claim.opening_pair.evaluation,
                },
            };
            verifier.update(&lcc, &inverse_vanishing_evals[idx]);
        }

        verifier
    }

    /// Instantiate, update, and finalize in one call.
    ///
    /// C++ source: `ShplonkVerifier_::reduce_verification`
    pub fn reduce_verification(
        g1_identity: &Bn254G1Affine,
        claims: &[OpeningClaim<Bn254G1Params>],
        transcript: &mut NativeTranscript,
    ) -> OpeningClaim<Bn254G1Params> {
        let mut verifier = Self::reduce_verification_no_finalize(claims, transcript);
        verifier.finalize(g1_identity)
    }
}

// ── Free functions ───────────────────────────────────────────────────────────

/// Precompute powers of nu needed to batch all univariate claims.
///
/// C++ source: `compute_shplonk_batching_challenge_powers`
pub fn compute_shplonk_batching_challenge_powers<P: FieldParams>(
    shplonk_batching_challenge: &Field<P>,
    virtual_log_n: usize,
    has_zk: bool,
    committed_sumcheck: bool,
) -> Vec<Field<P>> {
    let mut num_powers = 2 * virtual_log_n + NUM_INTERLEAVING_CLAIMS;
    const NUM_COMMITTED_SUMCHECK_CLAIMS_PER_ROUND: usize = 3;

    if has_zk {
        num_powers += NUM_SMALL_IPA_EVALUATIONS;
    }

    if committed_sumcheck {
        num_powers += NUM_COMMITTED_SUMCHECK_CLAIMS_PER_ROUND * virtual_log_n;
    }

    let mut result = Vec::with_capacity(num_powers);
    result.push(Field::one());
    for idx in 1..num_powers {
        result.push(result[idx - 1] * *shplonk_batching_challenge);
    }
    result
}

#[cfg(test)]
mod tests;
