//! Small Subgroup IPA: a ZK protocol to prove inner products of small vectors.
//!
//! Port of C++ `commitment_schemes/small_subgroup_ipa/`.
//!
//! The protocol proves that `<F, G> = s` for a witness polynomial G and a public
//! challenge polynomial F, where both are evaluated in the Lagrange basis over a
//! small multiplicative subgroup H of the scalar field.
//!
//! The verifier checks the "grand sum identity":
//!   L_1(r) * A(r) + (r - g^{-1}) * (A(g*r) - A(r) - F(r)*G(r)) + L_{|H|}(r) * (A(r) - s)
//!     = Z_H(r) * Q(r)
//!
//! where A is the grand sum polynomial, Q is the quotient, and L_1, L_{|H|} are
//! Lagrange polynomials for the first and last elements of H.

#[cfg(test)]
mod tests;

use std::marker::PhantomData;

use bbrs_ecc::curves::bn254::{Bn254FqParams, Bn254FrParams, Bn254G1Params};
use bbrs_ecc::curves::grumpkin::GrumpkinG1Params;
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;

// ── Global constants from C++ constants.hpp ──────────────────────────────────

/// Number of evaluations sent by the SmallSubgroupIPA prover.
/// Evaluations: G(r), A(g*r), A(r), Q(r).
pub const NUM_SMALL_IPA_EVALUATIONS: usize = 4;

/// Number of translation evaluations (ECCVM-specific).
pub const NUM_TRANSLATION_EVALUATIONS: usize = 5;

/// Number of masked rows in sumcheck.
pub const NUM_MASKED_ROWS: usize = 3;

/// Number of disabled rows in sumcheck = NUM_MASKED_ROWS + 1.
pub const NUM_DISABLED_ROWS_IN_SUMCHECK: usize = NUM_MASKED_ROWS + 1;

/// Max log circuit size (fixed proof size).
pub const CONST_PROOF_SIZE_LOG_N: usize = 28;

/// Number of barycentric evaluations computed by the verifier (F(r), L_1(r), L_{|H|}(r)).
const NUM_BARYCENTRIC_EVALUATIONS: usize = 3;

/// Convenience type alias: scalar field of a curve.
type FF<C> = Field<<C as CurveParams>::ScalarFieldParams>;

// ── SmallSubgroupCurveConfig trait ───────────────────────────────────────────

/// Trait providing curve-specific constants for the Small Subgroup IPA protocol.
///
/// Each curve has a multiplicative subgroup H of its scalar field with a fixed
/// generator `g`. The subgroup size and generator values are protocol constants.
pub trait SmallSubgroupCurveConfig: CurveParams {
    /// Size of the multiplicative subgroup H in the scalar field.
    const SUBGROUP_SIZE: usize;

    /// Length of random polynomials masking Sumcheck univariates.
    const LIBRA_UNIVARIATES_LENGTH: usize;

    /// Fixed generator of the subgroup H.
    fn subgroup_generator() -> FF<Self>;

    /// Inverse of the subgroup generator.
    fn subgroup_generator_inverse() -> FF<Self>;
}

// ── BN254 implementation ─────────────────────────────────────────────────────

impl SmallSubgroupCurveConfig for Bn254G1Params {
    const SUBGROUP_SIZE: usize = 256;
    const LIBRA_UNIVARIATES_LENGTH: usize = 9;

    fn subgroup_generator() -> Field<Bn254FrParams> {
        // C++: 0x07b0c561a6148404f086204a9f36ffb0617942546750f230c893619174a57a76
        Field::from_limbs([
            0xc893619174a57a76,
            0x617942546750f230,
            0xf086204a9f36ffb0,
            0x07b0c561a6148404,
        ])
    }

    fn subgroup_generator_inverse() -> Field<Bn254FrParams> {
        // C++: 0x204bd3277422fad364751ad938e2b5e6a54cf8c68712848a692c553d0329f5d6
        Field::from_limbs([
            0x692c553d0329f5d6,
            0xa54cf8c68712848a,
            0x64751ad938e2b5e6,
            0x204bd3277422fad3,
        ])
    }
}

// ── Grumpkin implementation ──────────────────────────────────────────────────

impl SmallSubgroupCurveConfig for GrumpkinG1Params {
    const SUBGROUP_SIZE: usize = 87;
    const LIBRA_UNIVARIATES_LENGTH: usize = 3;

    fn subgroup_generator() -> Field<Bn254FqParams> {
        // C++: 0x147c647c09fb639514909e9f0513f31ec1a523bf8a0880bc7c24fbc962a9586b
        Field::from_limbs([
            0x7c24fbc962a9586b,
            0xc1a523bf8a0880bc,
            0x14909e9f0513f31e,
            0x147c647c09fb6395,
        ])
    }

    fn subgroup_generator_inverse() -> Field<Bn254FqParams> {
        // C++: 0x0c68e27477b5e78cfab790bd3b59806fa871771f71ec7452cde5384f6e3a1988
        Field::from_limbs([
            0xcde5384f6e3a1988,
            0xa871771f71ec7452,
            0xfab790bd3b59806f,
            0x0c68e27477b5e78c,
        ])
    }
}

// ── Utility functions (from small_subgroup_ipa_utils.hpp) ────────────────────

/// Get the evaluation labels for SmallSubgroupIPA transcript entries.
///
/// Port of C++ `get_evaluation_labels`.
pub fn get_evaluation_labels(label_prefix: &str) -> [String; NUM_SMALL_IPA_EVALUATIONS] {
    [
        format!("{}concatenation_eval", label_prefix),
        format!("{}grand_sum_shift_eval", label_prefix),
        format!("{}grand_sum_eval", label_prefix),
        format!("{}quotient_eval", label_prefix),
    ]
}

/// Compute the evaluation points for SmallSubgroupIPA.
///
/// Returns [r, r*g, r, r] where g is the subgroup generator.
/// Port of C++ `compute_evaluation_points`.
pub fn compute_evaluation_points<P: FieldParams>(
    small_ipa_eval_challenge: &Field<P>,
    subgroup_gen: &Field<P>,
) -> [Field<P>; NUM_SMALL_IPA_EVALUATIONS] {
    [
        *small_ipa_eval_challenge,
        *small_ipa_eval_challenge * *subgroup_gen,
        *small_ipa_eval_challenge,
        *small_ipa_eval_challenge,
    ]
}

/// Commitments to the SmallSubgroupIPA witness polynomials [G], [A], [Q].
///
/// Port of C++ `SmallSubgroupIPACommitments`.
pub struct SmallSubgroupIPACommitments<C: CurveParams> {
    pub concatenated: AffineElement<C>,
    pub grand_sum: AffineElement<C>,
    pub quotient: AffineElement<C>,
}

impl<C: CurveParams> SmallSubgroupIPACommitments<C> {
    /// Returns [G], [A], [A], [Q] — grand_sum appears twice (opened at 2 points).
    pub fn get_all(&self) -> [&AffineElement<C>; NUM_SMALL_IPA_EVALUATIONS] {
        [
            &self.concatenated,
            &self.grand_sum,
            &self.grand_sum,
            &self.quotient,
        ]
    }
}

// ── Challenge polynomial construction ────────────────────────────────────────

/// Compute Lagrange-basis coefficients of the challenge polynomial F for ZK-Sumcheck.
///
/// Given multivariate challenge (u_0, ..., u_{D-1}), constructs the vector:
///   (1, u_0, u_0^2, ..., u_0^{L-1}, 1, u_1, ..., u_{D-1}^{L-1}, 0, ..., 0)
/// where L = LIBRA_UNIVARIATES_LENGTH.
///
/// Port of C++ `compute_challenge_polynomial_coeffs`.
pub fn compute_challenge_polynomial_coeffs<C: SmallSubgroupCurveConfig>(
    multivariate_challenge: &[FF<C>],
) -> Vec<FF<C>> {
    let mut coeffs = vec![FF::<C>::zero(); C::SUBGROUP_SIZE];
    let libra_len = C::LIBRA_UNIVARIATES_LENGTH;

    // First coefficient is 1
    coeffs[0] = FF::<C>::one();

    // Populate powers of each challenge
    for (round_idx, challenge) in multivariate_challenge.iter().enumerate() {
        let current_idx = 1 + libra_len * round_idx;
        coeffs[current_idx] = FF::<C>::one();
        for idx in (current_idx + 1)..(current_idx + libra_len) {
            coeffs[idx] = coeffs[idx - 1] * *challenge;
        }
    }

    coeffs
}

/// Compute Lagrange-basis coefficients of the challenge polynomial F for ECCVM translation.
///
/// Constructs coefficients: (1, x, x^2, ..., x^{M-1}, v, v*x, ..., v^{N-1}*x^{M-1}, 0, ..., 0)
/// where M = NUM_DISABLED_ROWS_IN_SUMCHECK and N = NUM_TRANSLATION_EVALUATIONS.
///
/// Port of C++ `compute_eccvm_challenge_coeffs`.
pub fn compute_eccvm_challenge_coeffs<C: SmallSubgroupCurveConfig>(
    evaluation_challenge_x: &FF<C>,
    batching_challenge_v: &FF<C>,
) -> Vec<FF<C>> {
    let mut coeffs = vec![FF::<C>::zero(); C::SUBGROUP_SIZE];
    let mut v_power = FF::<C>::one();

    for poly_idx in 0..NUM_TRANSLATION_EVALUATIONS {
        let start = NUM_DISABLED_ROWS_IN_SUMCHECK * poly_idx;
        coeffs[start] = v_power;

        for idx in (start + 1)..(start + NUM_DISABLED_ROWS_IN_SUMCHECK) {
            coeffs[idx] = coeffs[idx - 1] * *evaluation_challenge_x;
        }

        v_power = v_power * *batching_challenge_v;
    }

    coeffs
}

// ── SmallSubgroupIPAVerifier ─────────────────────────────────────────────────

/// Verifier for the Small Subgroup IPA protocol.
///
/// Checks the grand sum identity:
///   L_1(r)*A(r) + (r - g^{-1})*(A(g*r) - A(r) - F(r)*G(r)) + L_{|H|}(r)*(A(r) - s)
///     = Z_H(r) * Q(r)
///
/// Port of C++ `SmallSubgroupIPAVerifier<Curve>`.
pub struct SmallSubgroupIPAVerifier<C: SmallSubgroupCurveConfig> {
    _phantom: PhantomData<C>,
}

impl<C: SmallSubgroupCurveConfig> SmallSubgroupIPAVerifier<C> {
    /// Generic consistency check agnostic to the challenge polynomial F.
    ///
    /// Verifies: L_1(r)*A(r) + (r - g^{-1})*(A(g*r) - A(r) - F(r)*G(r))
    ///           + L_{|H|}(r)*(A(r) - s) - Z_H(r)*Q(r) == 0
    ///
    /// Port of C++ `SmallSubgroupIPAVerifier::check_consistency`.
    pub fn check_consistency(
        small_ipa_evaluations: &[FF<C>; NUM_SMALL_IPA_EVALUATIONS],
        small_ipa_eval_challenge: &FF<C>,
        challenge_polynomial: &[FF<C>],
        inner_product_eval_claim: &FF<C>,
        vanishing_poly_eval: &FF<C>,
    ) -> bool {
        // Check edge case: evaluation challenge in subgroup
        Self::handle_edge_cases(vanishing_poly_eval);

        // Compute evaluations at r of F, L_1, L_{|H|}
        let [challenge_poly, lagrange_first, lagrange_last] =
            Self::compute_batched_barycentric_evaluations(
                challenge_polynomial,
                small_ipa_eval_challenge,
                vanishing_poly_eval,
            );

        let concatenated_at_r = &small_ipa_evaluations[0];
        let grand_sum_shifted_eval = &small_ipa_evaluations[1];
        let grand_sum_eval = &small_ipa_evaluations[2];
        let quotient_eval = &small_ipa_evaluations[3];

        // Compute the grand sum identity check:
        // L_1(r)*A(r) + (r - 1/g)*(A(g*r) - A(r) - F(r)*G(r)) + L_{|H|}(r)*(A(r) - s) - Z_H(r)*Q(r)
        let mut diff = lagrange_first * *grand_sum_eval;
        diff = diff
            + (*small_ipa_eval_challenge - C::subgroup_generator_inverse())
                * (*grand_sum_shifted_eval - *grand_sum_eval - *concatenated_at_r * challenge_poly);
        diff = diff + lagrange_last * (*grand_sum_eval - *inner_product_eval_claim)
            - *vanishing_poly_eval * *quotient_eval;

        diff == FF::<C>::zero()
    }

    /// Consistency check for ZK-Sumcheck (Libra) evaluations.
    ///
    /// The challenge polynomial is constructed from concatenated powers of
    /// the sumcheck multivariate challenges.
    ///
    /// Port of C++ `check_libra_evaluations_consistency`.
    pub fn check_libra_evaluations_consistency(
        libra_evaluations: &[FF<C>; NUM_SMALL_IPA_EVALUATIONS],
        gemini_evaluation_challenge: &FF<C>,
        multilinear_challenge: &[FF<C>],
        inner_product_eval_claim: &FF<C>,
    ) -> bool {
        let vanishing_poly_eval =
            gemini_evaluation_challenge.pow(&[C::SUBGROUP_SIZE as u64, 0, 0, 0])
                - FF::<C>::one();

        Self::check_consistency(
            libra_evaluations,
            gemini_evaluation_challenge,
            &compute_challenge_polynomial_coeffs::<C>(multilinear_challenge),
            inner_product_eval_claim,
            &vanishing_poly_eval,
        )
    }

    /// Consistency check for ECCVM translation evaluations.
    ///
    /// The challenge polynomial is constructed from products of x^i * v^j.
    ///
    /// Port of C++ `check_eccvm_evaluations_consistency`.
    pub fn check_eccvm_evaluations_consistency(
        small_ipa_evaluations: &[FF<C>; NUM_SMALL_IPA_EVALUATIONS],
        evaluation_challenge: &FF<C>,
        evaluation_challenge_x: &FF<C>,
        batching_challenge_v: &FF<C>,
        inner_product_eval_claim: &FF<C>,
    ) -> bool {
        let vanishing_poly_eval =
            evaluation_challenge.pow(&[C::SUBGROUP_SIZE as u64, 0, 0, 0]) - FF::<C>::one();

        Self::check_consistency(
            small_ipa_evaluations,
            evaluation_challenge,
            &compute_eccvm_challenge_coeffs::<C>(evaluation_challenge_x, batching_challenge_v),
            inner_product_eval_claim,
            &vanishing_poly_eval,
        )
    }

    /// Panic if the evaluation challenge lands in the subgroup (negligible probability).
    ///
    /// Port of C++ `handle_edge_cases`.
    fn handle_edge_cases(vanishing_poly_eval: &FF<C>) {
        assert!(
            *vanishing_poly_eval != FF::<C>::zero(),
            "Evaluation challenge is in the SmallSubgroup."
        );
    }

    /// Efficient batch evaluation of the challenge polynomial, Lagrange first, and Lagrange last.
    ///
    /// Uses the barycentric formula over the interpolation domain (1, g, g^2, ..., g^{|H|-1}).
    /// Returns [F(r), L_1(r), L_{|H|}(r)].
    ///
    /// Port of C++ `compute_batched_barycentric_evaluations`.
    pub fn compute_batched_barycentric_evaluations(
        coeffs: &[FF<C>],
        r: &FF<C>,
        vanishing_poly_eval: &FF<C>,
    ) -> [FF<C>; NUM_BARYCENTRIC_EVALUATIONS] {
        let n = C::SUBGROUP_SIZE;

        // Construct denominators: d_i = r * g^{-i} - 1
        let mut denominators = vec![FF::<C>::zero(); n];
        let g_inv = C::subgroup_generator_inverse();
        let mut running_power = FF::<C>::one();

        for i in 0..n {
            denominators[i] = running_power * *r - FF::<C>::one();
            running_power = running_power * g_inv;
        }

        // Batch invert denominators
        FF::<C>::batch_invert(&mut denominators);

        // Compute inner product of coeffs and inverted denominators
        let challenge_poly_eval: FF<C> = coeffs
            .iter()
            .zip(denominators.iter())
            .fold(FF::<C>::zero(), |acc, (c, d)| acc + *c * *d);

        // numerator = Z_H(r) / |H| = (r^n - 1) / n
        let n_inv = FF::<C>::from(n as u64).invert();
        let numerator = *vanishing_poly_eval * n_inv;

        [
            challenge_poly_eval * numerator,
            denominators[0] * numerator,
            denominators[n - 1] * numerator,
        ]
    }

    /// Compute the evaluation points for this curve's SmallSubgroupIPA.
    pub fn evaluation_points(
        small_ipa_eval_challenge: &FF<C>,
    ) -> [FF<C>; NUM_SMALL_IPA_EVALUATIONS] {
        compute_evaluation_points(small_ipa_eval_challenge, &C::subgroup_generator())
    }

    /// Get the evaluation labels for this curve's SmallSubgroupIPA.
    pub fn evaluation_labels(label_prefix: &str) -> [String; NUM_SMALL_IPA_EVALUATIONS] {
        get_evaluation_labels(label_prefix)
    }
}
