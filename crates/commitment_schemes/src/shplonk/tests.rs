//! Tests for Shplonk and Shplemini commitment schemes.
//!
//! Port of C++ shplonk.test.cpp and shplemini.test.cpp.

use bbrs_ecc::curves::bn254::{
    Bn254FrParams, Bn254G1Params, Fq12, Fr, G1Affine as Bn254G1Affine, G2AffineElement, G2Element,
};
use bbrs_ecc::curves::bn254_pairing::reduced_ate_pairing;
use bbrs_ecc::groups::element::Element;
use bbrs_polynomials::polynomial::Polynomial;

use crate::batch_mul::batch_mul_native;
use crate::claim::{OpeningClaim, OpeningPair};
use crate::claim_batcher::{Batch, ClaimBatcher};
use crate::commitment_key::CommitmentKey;
use crate::gemini::{self, GeminiProver, GeminiVerifier, PolynomialBatcher};
use crate::kzg::KZG;

use super::shplemini::{ShpleminiProver, ShpleminiVerifier};
use super::{
    compute_shplonk_batching_challenge_powers, LinearCombinationOfClaims, ProverOpeningClaim,
    ShplonkProver, ShplonkVerifier,
};

use bbrs_transcript::NativeTranscript;

/// Copy an OpeningPair (Bn254FrParams doesn't implement Clone but Field is Copy).
fn copy_opening_pair(pair: &OpeningPair<Bn254FrParams>) -> OpeningPair<Bn254FrParams> {
    OpeningPair {
        challenge: pair.challenge,
        evaluation: pair.evaluation,
    }
}

const LOG_DEGREE: usize = 4;
const MAX_POLY_DEGREE: usize = 1 << LOG_DEGREE;

/// Generate a powers-of-tau SRS for testing.
fn create_test_srs(n: usize) -> (Vec<Bn254G1Affine>, G2AffineElement) {
    let tau = Fr::random_element();

    let g1 = Element::<Bn254G1Params>::one();
    let mut points = Vec::with_capacity(n);
    let mut tau_power = Fr::one();
    for _ in 0..n {
        points.push(g1.mul(&tau_power).to_affine());
        tau_power = tau_power * tau;
    }

    let g2 = G2Element::from_affine(&G2AffineElement::generator());
    let g2_x = g2.mul_scalar(&tau).to_affine();

    (points, g2_x)
}

/// Direct pairing check: e(P0, [1]_2) * e(P1, [x]_2) == 1_T.
fn direct_pairing_check(
    p0: &Bn254G1Affine,
    p1: &Bn254G1Affine,
    g2_x: &G2AffineElement,
) -> bool {
    let e0 = reduced_ate_pairing(p0, &G2AffineElement::generator());
    let e1 = reduced_ate_pairing(p1, g2_x);
    (e0 * e1) == Fq12::one()
}

/// Claim data for testing: polynomial, commitment, and opening pair.
struct ClaimData {
    polynomial: Polynomial<Bn254FrParams>,
    commitment: Bn254G1Affine,
    opening_pair: OpeningPair<Bn254FrParams>,
}

/// Generate random claim data for a given polynomial size.
fn generate_claim_data(
    ck: &CommitmentKey<Bn254G1Params>,
    poly_size: usize,
) -> ClaimData {
    let polynomial = Polynomial::random(poly_size, poly_size, 0);
    let challenge = Fr::random_element();
    let evaluation = polynomial.evaluate(&challenge);
    let commitment = ck.commit(&polynomial);
    ClaimData {
        polynomial,
        commitment,
        opening_pair: OpeningPair {
            challenge,
            evaluation,
        },
    }
}

/// Verify that p(r) = v for a prover opening claim.
fn verify_opening_pair(opening_pair: &OpeningPair<Bn254FrParams>, poly: &Polynomial<Bn254FrParams>) {
    let expected = poly.evaluate(&opening_pair.challenge);
    assert_eq!(
        opening_pair.evaluation, expected,
        "OpeningPair: evaluations mismatch"
    );
}

/// Verify that C = Commit(p) and p(r) = v.
fn verify_opening_claim(
    claim: &OpeningClaim<Bn254G1Params>,
    witness: &Polynomial<Bn254FrParams>,
    ck: &CommitmentKey<Bn254G1Params>,
) {
    let expected_eval = witness.evaluate(&claim.opening_pair.challenge);
    assert_eq!(
        claim.opening_pair.evaluation, expected_eval,
        "OpeningClaim: evaluations mismatch"
    );
    let expected_commitment = ck.commit(witness);
    assert_eq!(
        claim.commitment, expected_commitment,
        "OpeningClaim: commitment mismatch"
    );
}

// ── ShplonkSimple test ──────────────────────────────────────────────────────

/// Test of Shplonk prover/verifier for two polynomials of different sizes,
/// each opened at a single (different) point.
///
/// Port of C++ `TYPED_TEST(ShplonkTest, ShplonkSimple)` for BN254.
#[test]
fn shplonk_simple_bn254() {
    let (srs_points, _g2_x) = create_test_srs(MAX_POLY_DEGREE);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points.clone());

    let mut prover_transcript = NativeTranscript::prover_init_empty();

    // Generate two random polynomials of different sizes
    let setup = vec![
        generate_claim_data(&ck, MAX_POLY_DEGREE),
        generate_claim_data(&ck, MAX_POLY_DEGREE / 2),
    ];

    // Create prover opening claims
    let mut prover_opening_claims: Vec<ProverOpeningClaim<Bn254FrParams>> = setup
        .iter()
        .map(|cd| ProverOpeningClaim {
            polynomial: cd.polynomial.clone(),
            opening_pair: copy_opening_pair(&cd.opening_pair),
            gemini_fold: false,
        })
        .collect();

    // Execute the Shplonk prover
    let batched_opening_claim = ShplonkProver::prove(
        &ck,
        &mut prover_opening_claims,
        &mut prover_transcript,
        &mut [],
        &mut [],
        0,
    );

    // Intermediate check: verify the opening pair of the prover witness Q
    verify_opening_pair(
        &batched_opening_claim.opening_pair,
        &batched_opening_claim.polynomial,
    );

    // Initialize verifier transcript from prover transcript
    let mut verifier_transcript =
        NativeTranscript::verifier_init_empty(&mut prover_transcript);

    // Create verifier opening claims
    let verifier_claims: Vec<OpeningClaim<Bn254G1Params>> = setup
        .iter()
        .map(|cd| OpeningClaim {
            opening_pair: copy_opening_pair(&cd.opening_pair),
            commitment: cd.commitment,
        })
        .collect();

    // Get g1 identity
    let g1_identity = srs_points[0];

    // Execute the Shplonk verifier
    let batched_verifier_claim = ShplonkVerifier::<Bn254G1Params>::reduce_verification(
        &g1_identity,
        &verifier_claims,
        &mut verifier_transcript,
    );

    // Verify the claim
    verify_opening_claim(
        &batched_verifier_claim,
        &batched_opening_claim.polynomial,
        &ck,
    );
}

// ── ShplonkLinearlyDependent test ───────────────────────────────────────────

/// Test of Shplonk prover/verifier for polynomials that are linearly dependent.
///
/// Port of C++ `TYPED_TEST(ShplonkTest, ShplonkLinearlyDependent)` for BN254.
#[test]
fn shplonk_linearly_dependent_bn254() {
    let (srs_points, _g2_x) = create_test_srs(MAX_POLY_DEGREE);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points.clone());

    let mut prover_transcript = NativeTranscript::prover_init_empty();

    // Generate two random polynomials of different sizes
    let setup = vec![
        generate_claim_data(&ck, MAX_POLY_DEGREE),
        generate_claim_data(&ck, MAX_POLY_DEGREE / 2),
    ];

    // Generate random coefficients for linear combination
    let coefficients = vec![Fr::random_element(), Fr::random_element()];

    // Create the linear combination polynomial and evaluate it at a random point
    let lc_challenge = Fr::random_element();
    let mut lc_poly = Polynomial::<Bn254FrParams>::new(MAX_POLY_DEGREE, MAX_POLY_DEGREE, 0);
    let mut lc_eval = Fr::zero();

    for (coeff, cd) in coefficients.iter().zip(setup.iter()) {
        lc_poly.add_scaled(&cd.polynomial.as_span(), *coeff);
        lc_eval = lc_eval + *coeff * cd.polynomial.evaluate(&lc_challenge);
    }

    let lc_commitment = ck.commit(&lc_poly);

    // Extract commitments for the verifier
    let mut commitments: Vec<Bn254G1Affine> = setup.iter().map(|cd| cd.commitment).collect();

    // Create prover opening claims (original + linear combination)
    let mut prover_opening_claims: Vec<ProverOpeningClaim<Bn254FrParams>> = setup
        .iter()
        .map(|cd| ProverOpeningClaim {
            polynomial: cd.polynomial.clone(),
            opening_pair: copy_opening_pair(&cd.opening_pair),
            gemini_fold: false,
        })
        .collect();
    prover_opening_claims.push(ProverOpeningClaim {
        polynomial: lc_poly.clone(),
        opening_pair: OpeningPair {
            challenge: lc_challenge,
            evaluation: lc_eval,
        },
        gemini_fold: false,
    });

    // Execute the Shplonk prover
    let batched_opening_claim = ShplonkProver::prove(
        &ck,
        &mut prover_opening_claims,
        &mut prover_transcript,
        &mut [],
        &mut [],
        0,
    );

    // Intermediate check
    verify_opening_pair(
        &batched_opening_claim.opening_pair,
        &batched_opening_claim.polynomial,
    );

    // Create verifier opening claims
    let verifier_opening_claims: Vec<OpeningClaim<Bn254G1Params>> = vec![
        OpeningClaim {
            opening_pair: copy_opening_pair(&setup[0].opening_pair),
            commitment: setup[0].commitment,
        },
        OpeningClaim {
            opening_pair: copy_opening_pair(&setup[1].opening_pair),
            commitment: setup[1].commitment,
        },
        OpeningClaim {
            opening_pair: OpeningPair {
                challenge: lc_challenge,
                evaluation: lc_eval,
            },
            commitment: lc_commitment,
        },
    ];

    // Create update data for the verifier using LinearCombinationOfClaims
    let update_data = vec![
        LinearCombinationOfClaims {
            indices: vec![0],
            scalars: vec![Fr::one()],
            opening_pair: copy_opening_pair(&verifier_opening_claims[0].opening_pair),
        },
        LinearCombinationOfClaims {
            indices: vec![1],
            scalars: vec![Fr::one()],
            opening_pair: copy_opening_pair(&verifier_opening_claims[1].opening_pair),
        },
        LinearCombinationOfClaims {
            indices: vec![0, 1],
            scalars: coefficients.clone(),
            opening_pair: copy_opening_pair(&verifier_opening_claims[2].opening_pair),
        },
    ];

    let mut verifier_transcript =
        NativeTranscript::verifier_init_empty(&mut prover_transcript);

    let mut verifier = ShplonkVerifier::<Bn254G1Params>::new(
        &mut commitments,
        &mut verifier_transcript,
        verifier_opening_claims.len(),
    );

    let g1_identity = srs_points[0];

    // Execute the Shplonk verifier
    let batched_verifier_claim =
        verifier.reduce_verification_vector_claims(&g1_identity, &update_data);

    verify_opening_claim(
        &batched_verifier_claim,
        &batched_opening_claim.polynomial,
        &ck,
    );
}

// ── ExportBatchClaimAndVerify test ──────────────────────────────────────────

/// Test exporting batch claim from Shplonk verifier and verifying with KZG.
///
/// Port of C++ `TYPED_TEST(ShplonkTest, ExportBatchClaimAndVerify)` for BN254.
#[test]
fn export_batch_claim_and_verify_bn254() {
    let (srs_points, g2_x) = create_test_srs(MAX_POLY_DEGREE);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points.clone());

    let mut prover_transcript = NativeTranscript::prover_init_empty();

    // Generate claim data
    let setup = vec![
        generate_claim_data(&ck, MAX_POLY_DEGREE),
        generate_claim_data(&ck, MAX_POLY_DEGREE / 2),
    ];

    // Prover
    let mut prover_opening_claims: Vec<ProverOpeningClaim<Bn254FrParams>> = setup
        .iter()
        .map(|cd| ProverOpeningClaim {
            polynomial: cd.polynomial.clone(),
            opening_pair: copy_opening_pair(&cd.opening_pair),
            gemini_fold: false,
        })
        .collect();

    let batched_opening_claim = ShplonkProver::prove(
        &ck,
        &mut prover_opening_claims,
        &mut prover_transcript,
        &mut [],
        &mut [],
        0,
    );

    verify_opening_pair(
        &batched_opening_claim.opening_pair,
        &batched_opening_claim.polynomial,
    );

    // Compute KZG proof on the batched opening claim
    let kzg_claim = crate::claim::ProverOpeningClaim::<Bn254G1Params> {
        polynomial: batched_opening_claim.polynomial,
        opening_pair: copy_opening_pair(&batched_opening_claim.opening_pair),
        gemini_fold: false,
    };
    KZG::compute_opening_proof(&ck, kzg_claim, &mut prover_transcript);

    // Verifier
    let mut verifier_transcript =
        NativeTranscript::verifier_init_empty(&mut prover_transcript);

    let verifier_claims: Vec<OpeningClaim<Bn254G1Params>> = setup
        .iter()
        .map(|cd| OpeningClaim {
            opening_pair: copy_opening_pair(&cd.opening_pair),
            commitment: cd.commitment,
        })
        .collect();

    let mut verifier = ShplonkVerifier::<Bn254G1Params>::reduce_verification_no_finalize(
        &verifier_claims,
        &mut verifier_transcript,
    );

    let g1_identity = srs_points[0];

    // Export batch opening claim
    let batch_claim = verifier.export_batch_opening_claim(&g1_identity);

    // Verify with KZG
    let pairing_points =
        KZG::reduce_verify_batch_opening_claim(batch_claim, &mut verifier_transcript);

    assert!(
        direct_pairing_check(&pairing_points.p0, &pairing_points.p1, &g2_x),
        "Shplonk + KZG pairing check failed"
    );
}

// ── Unit tests for helpers ──────────────────────────────────────────────────

#[test]
fn compute_inverted_gemini_denominators_basic() {
    let z = Fr::from(7u64);
    let powers = vec![Fr::from(2u64), Fr::from(4u64)];

    let result =
        ShplonkVerifier::<Bn254G1Params>::compute_inverted_gemini_denominators(&z, &powers);

    assert_eq!(result.len(), 4);

    // Check: 1/(z-r) = 1/(7-2) = 1/5
    let expected_0 = (z - powers[0]).invert();
    assert_eq!(result[0], expected_0);

    // Check: 1/(z+r) = 1/(7+2) = 1/9
    let expected_1 = (z + powers[0]).invert();
    assert_eq!(result[1], expected_1);

    // Check: 1/(z-r^2) = 1/(7-4) = 1/3
    let expected_2 = (z - powers[1]).invert();
    assert_eq!(result[2], expected_2);

    // Check: 1/(z+r^2) = 1/(7+4) = 1/11
    let expected_3 = (z + powers[1]).invert();
    assert_eq!(result[3], expected_3);
}

#[test]
fn compute_shplonk_batching_challenge_powers_basic() {
    use super::compute_shplonk_batching_challenge_powers;

    let nu = Fr::from(3u64);
    let virtual_log_n = 2;

    let result = compute_shplonk_batching_challenge_powers(&nu, virtual_log_n, false, false);

    // num_powers = 2 * 2 + 2 = 6
    assert_eq!(result.len(), 6);
    assert_eq!(result[0], Fr::one());
    assert_eq!(result[1], nu);
    assert_eq!(result[2], nu * nu);
    assert_eq!(result[3], nu * nu * nu);
}

// ── Shplemini Tests (port of C++ shplemini.test.cpp) ────────────────────────

const SHPLEMINI_LOG_N: usize = 9;
const SHPLEMINI_N: usize = 1 << SHPLEMINI_LOG_N;
const NUM_POLYNOMIALS: usize = 7;
const NUM_SHIFTABLE: usize = 2;

/// Mock claim generator for Shplemini tests.
///
/// Port of C++ `MockClaimGenerator<Curve>`.
struct MockClaimGenerator {
    unshifted_polys: Vec<Polynomial<Bn254FrParams>>,
    unshifted_commitments: Vec<Bn254G1Affine>,
    unshifted_evals: Vec<Fr>,
    to_be_shifted_polys: Vec<Polynomial<Bn254FrParams>>,
    to_be_shifted_commitments: Vec<Bn254G1Affine>,
    to_be_shifted_evals: Vec<Fr>,
}

impl MockClaimGenerator {
    /// Construct claim data for random polynomials.
    ///
    /// `num_to_be_shifted` polynomials get both an unshifted and a shifted claim.
    fn new(
        poly_size: usize,
        num_polynomials: usize,
        num_to_be_shifted: usize,
        mle_opening_point: &[Fr],
        ck: &CommitmentKey<Bn254G1Params>,
    ) -> Self {
        assert!(num_polynomials >= num_to_be_shifted);
        let num_not_to_be_shifted = num_polynomials - num_to_be_shifted;

        let mut unshifted_polys = Vec::new();
        let mut unshifted_commitments = Vec::new();
        let mut unshifted_evals = Vec::new();
        let mut to_be_shifted_polys = Vec::new();
        let mut to_be_shifted_commitments = Vec::new();
        let mut to_be_shifted_evals = Vec::new();

        // Polynomials that are NOT to be shifted
        for _ in 0..num_not_to_be_shifted {
            let poly = Polynomial::random(poly_size, poly_size, 0);
            let commitment = ck.commit(&poly);
            let eval = poly.evaluate_mle(mle_opening_point, false);
            unshifted_polys.push(poly);
            unshifted_commitments.push(commitment);
            unshifted_evals.push(eval);
        }

        // Polynomials to be shifted: each gets both an unshifted and shifted claim
        for _ in 0..num_to_be_shifted {
            // Make shiftable: coeff[0] = 0
            let mut poly = Polynomial::random(poly_size, poly_size, 0);
            *poly.at_mut(0) = Fr::zero();
            let commitment = ck.commit(&poly);

            // Shifted evaluation: evaluate shifted(poly) = poly(X) / X at mle point
            let shifted_eval = poly.evaluate_mle(mle_opening_point, true);
            to_be_shifted_commitments.push(commitment);
            to_be_shifted_evals.push(shifted_eval);
            to_be_shifted_polys.push(poly.clone());

            // Unshifted counterpart
            let unshifted_eval = poly.evaluate_mle(mle_opening_point, false);
            unshifted_commitments.push(commitment);
            unshifted_evals.push(unshifted_eval);
            unshifted_polys.push(poly);
        }

        Self {
            unshifted_polys,
            unshifted_commitments,
            unshifted_evals,
            to_be_shifted_polys,
            to_be_shifted_commitments,
            to_be_shifted_evals,
        }
    }

    /// Build a `PolynomialBatcher` from the mock data.
    fn polynomial_batcher(&self) -> PolynomialBatcher<'_, Bn254FrParams> {
        let mut batcher = PolynomialBatcher::new(self.unshifted_polys[0].size());
        for poly in &self.unshifted_polys {
            batcher.unshifted.push(poly);
        }
        for poly in &self.to_be_shifted_polys {
            batcher.to_be_shifted_by_one.push(poly);
        }
        batcher
    }

    /// Build a `ClaimBatcher` from the mock data.
    fn claim_batcher(&self) -> ClaimBatcher<Bn254G1Params> {
        let mut cb = ClaimBatcher::new();
        cb.unshifted = Some(Batch {
            commitments: self.unshifted_commitments.clone(),
            evaluations: self.unshifted_evals.clone(),
            scalar: Fr::zero(),
        });
        if !self.to_be_shifted_polys.is_empty() {
            cb.shifted = Some(Batch {
                commitments: self.to_be_shifted_commitments.clone(),
                evaluations: self.to_be_shifted_evals.clone(),
                scalar: Fr::zero(),
            });
        }
        cb
    }

    /// Build a MockClaimGenerator with a single custom polynomial.
    fn with_custom_poly(
        _poly_size: usize,
        poly: Polynomial<Bn254FrParams>,
        claimed_eval: Fr,
        ck: &CommitmentKey<Bn254G1Params>,
    ) -> Self {
        let commitment = ck.commit(&poly);
        Self {
            unshifted_polys: vec![poly],
            unshifted_commitments: vec![commitment],
            unshifted_evals: vec![claimed_eval],
            to_be_shifted_polys: Vec::new(),
            to_be_shifted_commitments: Vec::new(),
            to_be_shifted_evals: Vec::new(),
        }
    }
}

/// Generate a random multilinear evaluation point.
fn random_evaluation_point(n: usize) -> Vec<Fr> {
    (0..n).map(|_| Fr::random_element()).collect()
}

/// Test that batch_multivariate_opening_claims operates correctly.
///
/// Port of C++ `TYPED_TEST(ShpleminiTest, CorrectnessOfMultivariateClaimBatching)` for BN254.
#[test]
fn shplemini_correctness_of_multivariate_claim_batching() {
    let (srs_points, _g2_x) = create_test_srs(SHPLEMINI_N);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points.clone());

    // Generate mock challenges
    let rho = Fr::random_element();
    let gemini_eval_challenge = Fr::random_element();
    let shplonk_batching_challenge = Fr::random_element();
    let shplonk_eval_challenge = Fr::random_element();

    // Generate multilinear polynomials and compute their commitments
    let mle_opening_point = random_evaluation_point(SHPLEMINI_LOG_N);

    let mock = MockClaimGenerator::new(SHPLEMINI_N, NUM_POLYNOMIALS, NUM_SHIFTABLE, &mle_opening_point, &ck);

    // Compute batched evaluation manually
    let mut rho_power = Fr::one();
    let mut batched_evaluation = Fr::zero();
    for eval in &mock.unshifted_evals {
        batched_evaluation = batched_evaluation + *eval * rho_power;
        rho_power = rho_power * rho;
    }
    for eval in &mock.to_be_shifted_evals {
        batched_evaluation = batched_evaluation + *eval * rho_power;
        rho_power = rho_power * rho;
    }

    // Compute batched commitments manually
    rho_power = Fr::one();
    let mut batched_commitment_unshifted = Element::<Bn254G1Params>::infinity();
    for comm in &mock.unshifted_commitments {
        batched_commitment_unshifted = batched_commitment_unshifted
            + Element::<Bn254G1Params>::from_affine(comm).mul(&rho_power);
        rho_power = rho_power * rho;
    }
    let mut batched_commitment_to_be_shifted = Element::<Bn254G1Params>::infinity();
    for comm in &mock.to_be_shifted_commitments {
        batched_commitment_to_be_shifted = batched_commitment_to_be_shifted
            + Element::<Bn254G1Params>::from_affine(comm).mul(&rho_power);
        rho_power = rho_power * rho;
    }

    // Compute expected result
    let to_be_shifted_contribution =
        batched_commitment_to_be_shifted.mul(&gemini_eval_challenge.invert());

    let commitment_to_univariate_pos = batched_commitment_unshifted + to_be_shifted_contribution;
    let commitment_to_univariate_neg = batched_commitment_unshifted - to_be_shifted_contribution;

    let expected_result = commitment_to_univariate_pos
        .mul(&(shplonk_eval_challenge - gemini_eval_challenge).invert())
        + commitment_to_univariate_neg.mul(
            &(shplonk_batching_challenge
                * (shplonk_eval_challenge + gemini_eval_challenge).invert()),
        );

    // Run the ShpleminiVerifier batching method
    let mut commitments = Vec::new();
    let mut scalars = Vec::new();
    let mut verifier_batched_evaluation = Fr::zero();

    let inverted_vanishing_eval_pos =
        (shplonk_eval_challenge - gemini_eval_challenge).invert();
    let inverted_vanishing_eval_neg =
        (shplonk_eval_challenge + gemini_eval_challenge).invert();

    let inverted_vanishing_evals = vec![inverted_vanishing_eval_pos, inverted_vanishing_eval_neg];

    let mut claim_batcher = mock.claim_batcher();
    claim_batcher.compute_scalars_for_each_batch(
        &inverted_vanishing_evals,
        &shplonk_batching_challenge,
        &gemini_eval_challenge,
    );

    claim_batcher.update_batch_mul_inputs_and_batched_evaluation(
        &mut commitments,
        &mut scalars,
        &mut verifier_batched_evaluation,
        &rho,
        Fr::zero(),
        Fr::zero(),
    );

    // Final check
    let shplemini_result = batch_mul_native::<Bn254G1Params>(&commitments, &scalars);

    assert_eq!(
        commitments.len(),
        mock.unshifted_commitments.len() + mock.to_be_shifted_commitments.len()
    );
    assert_eq!(batched_evaluation, verifier_batched_evaluation);
    assert_eq!((-expected_result).to_affine(), shplemini_result);
}

/// Test correctness of Gemini claim batching.
///
/// Port of C++ `TYPED_TEST(ShpleminiTest, CorrectnessOfGeminiClaimBatching)` for BN254.
#[test]
fn shplemini_correctness_of_gemini_claim_batching() {
    let (srs_points, _g2_x) = create_test_srs(SHPLEMINI_N);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points.clone());

    // Generate mock challenges
    let rho = Fr::random_element();
    let gemini_eval_challenge = Fr::random_element();
    let shplonk_batching_challenge = Fr::random_element();

    let shplonk_batching_challenge_powers =
        compute_shplonk_batching_challenge_powers(&shplonk_batching_challenge, SHPLEMINI_LOG_N, false, false);

    let shplonk_eval_challenge = Fr::random_element();

    let mle_opening_point = random_evaluation_point(SHPLEMINI_LOG_N);

    let mock = MockClaimGenerator::new(SHPLEMINI_N, NUM_POLYNOMIALS, NUM_SHIFTABLE, &mle_opening_point, &ck);

    let mut batcher = mock.polynomial_batcher();

    let batched = batcher.compute_batched(rho);

    // Compute fold polynomials
    let fold_polynomials =
        GeminiProver::compute_fold_polynomials(SHPLEMINI_LOG_N, &mle_opening_point, &batched);

    let mut prover_commitments = Vec::new();
    for l in 0..SHPLEMINI_LOG_N - 1 {
        let commitment = ck.commit(&fold_polynomials[l]);
        prover_commitments.push(commitment);
    }

    let (a_0_pos, a_0_neg) = batcher.compute_partially_evaluated_batch_polynomials(gemini_eval_challenge);

    let opening_claims = GeminiProver::construct_univariate_opening_claims(
        SHPLEMINI_LOG_N,
        a_0_pos,
        a_0_neg,
        fold_polynomials,
        gemini_eval_challenge,
    );

    let mut prover_evaluations = Vec::new();
    for l in 0..SHPLEMINI_LOG_N {
        prover_evaluations.push(opening_claims[l + 1].opening_pair.evaluation);
    }

    let r_squares = gemini::powers_of_evaluation_challenge(gemini_eval_challenge, SHPLEMINI_LOG_N);

    // Compute expected result
    let mut expected_result = Element::<Bn254G1Params>::infinity();
    let mut expected_inverse_vanishing_evals = Vec::with_capacity(2 * SHPLEMINI_LOG_N);
    for idx in 0..SHPLEMINI_LOG_N {
        expected_inverse_vanishing_evals.push((shplonk_eval_challenge - r_squares[idx]).invert());
        expected_inverse_vanishing_evals.push((shplonk_eval_challenge + r_squares[idx]).invert());
    }

    let mut current_challenge = shplonk_batching_challenge * shplonk_batching_challenge;
    for idx in 0..prover_commitments.len() {
        let c = Element::<Bn254G1Params>::from_affine(&prover_commitments[idx]);
        expected_result = expected_result
            - c.mul(&(current_challenge * expected_inverse_vanishing_evals[2 * idx + 2]));
        current_challenge = current_challenge * shplonk_batching_challenge;
        expected_result = expected_result
            - c.mul(&(current_challenge * expected_inverse_vanishing_evals[2 * idx + 3]));
        current_challenge = current_challenge * shplonk_batching_challenge;
    }

    // Run the ShpleminiVerifier batching method
    let inverse_vanishing_evals =
        ShplonkVerifier::<Bn254G1Params>::compute_inverted_gemini_denominators(
            &shplonk_eval_challenge,
            &r_squares,
        );

    let mut expected_constant_term_accumulator = Fr::zero();
    let padding_indicator_array = vec![Fr::one(); SHPLEMINI_LOG_N];

    let gemini_fold_pos_evaluations = GeminiVerifier::compute_fold_pos_evaluations(
        &padding_indicator_array,
        expected_constant_term_accumulator,
        &mle_opening_point,
        &r_squares,
        &prover_evaluations,
        expected_constant_term_accumulator,
    );

    let mut commitments = Vec::new();
    let mut scalars = Vec::new();

    ShpleminiVerifier::batch_gemini_claims_received_from_prover::<Bn254G1Params>(
        &padding_indicator_array,
        &prover_commitments,
        &prover_evaluations,
        &gemini_fold_pos_evaluations,
        &inverse_vanishing_evals,
        &shplonk_batching_challenge_powers,
        &mut commitments,
        &mut scalars,
        &mut expected_constant_term_accumulator,
    );

    // Compute the group element using the output of Shplemini method
    let shplemini_result = batch_mul_native::<Bn254G1Params>(&commitments, &scalars);

    assert_eq!(shplemini_result, expected_result.to_affine());
}

/// End-to-end Shplemini test: prove and verify with KZG.
///
/// Port of C++ `TYPED_TEST(ShpleminiTest, HighDegreeAttackAccept)` for BN254,
/// but using a correct polynomial (no attack).
#[test]
fn shplemini_end_to_end_bn254() {
    let (srs_points, g2_x) = create_test_srs(SHPLEMINI_N);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points.clone());

    let mle_opening_point = random_evaluation_point(SHPLEMINI_LOG_N);

    let mock = MockClaimGenerator::new(SHPLEMINI_N, NUM_POLYNOMIALS, NUM_SHIFTABLE, &mle_opening_point, &ck);

    let mut prover_transcript = NativeTranscript::prover_init_empty();

    let mut batcher = mock.polynomial_batcher();

    // Run Shplemini prover
    let opening_claim = ShpleminiProver::prove_without_zk(
        SHPLEMINI_N,
        &mut batcher,
        &mle_opening_point,
        &ck,
        &mut prover_transcript,
    );

    // Run KZG prover on the batched opening claim
    let kzg_claim = crate::claim::ProverOpeningClaim::<Bn254G1Params> {
        polynomial: opening_claim.polynomial,
        opening_pair: OpeningPair {
            challenge: opening_claim.opening_pair.challenge,
            evaluation: opening_claim.opening_pair.evaluation,
        },
        gemini_fold: false,
    };
    KZG::compute_opening_proof(&ck, kzg_claim, &mut prover_transcript);

    // Verifier side
    let mut verifier_transcript = NativeTranscript::verifier_init_empty(&mut prover_transcript);

    let padding_indicator_array = vec![Fr::one(); SHPLEMINI_LOG_N];
    let mut claim_batcher = mock.claim_batcher();

    let g1_identity = srs_points[0];

    let batch_opening_claim = ShpleminiVerifier::compute_batch_opening_claim_without_zk(
        &padding_indicator_array,
        &mut claim_batcher,
        &mle_opening_point,
        &g1_identity,
        &mut verifier_transcript,
    );

    // Verify with KZG
    let pairing_points =
        KZG::reduce_verify_batch_opening_claim(batch_opening_claim, &mut verifier_transcript);

    assert!(
        direct_pairing_check(&pairing_points.p0, &pairing_points.p1, &g2_x),
        "Shplemini + KZG pairing check failed"
    );
}

/// High degree attack test: prover commits to a crafted higher degree polynomial.
///
/// Port of C++ `TYPED_TEST(ShpleminiTest, HighDegreeAttackAccept)` for BN254.
#[test]
fn shplemini_high_degree_attack_accept() {
    let small_log_n: usize = 3;
    let (srs_points, g2_x) = create_test_srs(SHPLEMINI_N);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points.clone());

    let u = random_evaluation_point(small_log_n);
    let claimed_multilinear_eval = Fr::random_element();

    // Construct a polynomial that folds to the claimed eval after small_log_n rounds
    let mut poly = Polynomial::<Bn254FrParams>::new(SHPLEMINI_N, SHPLEMINI_N, 0);
    let tail = ((Fr::one() - u[0]) * (Fr::one() - u[1])).invert();
    *poly.at_mut(4) = claimed_multilinear_eval * tail * u[2].invert();
    *poly.at_mut(SHPLEMINI_N - 8) = tail;
    *poly.at_mut(SHPLEMINI_N - 4) = -tail * (Fr::one() - u[2]) * u[2].invert();

    let mock = MockClaimGenerator::with_custom_poly(SHPLEMINI_N, poly, claimed_multilinear_eval, &ck);

    let mut prover_transcript = NativeTranscript::prover_init_empty();
    let mut batcher = mock.polynomial_batcher();

    // circuit_size = SHPLEMINI_N: Gemini folds for log2(N) rounds over the full polynomial,
    // but only small_log_n challenges are meaningful. The polynomial is crafted so that
    // the Gemini folding produces the correct evaluation despite the high degree.
    let opening_claim = ShpleminiProver::prove_without_zk(
        SHPLEMINI_N,
        &mut batcher,
        &u,
        &ck,
        &mut prover_transcript,
    );

    // KZG prover
    let kzg_claim = crate::claim::ProverOpeningClaim::<Bn254G1Params> {
        polynomial: opening_claim.polynomial,
        opening_pair: OpeningPair {
            challenge: opening_claim.opening_pair.challenge,
            evaluation: opening_claim.opening_pair.evaluation,
        },
        gemini_fold: false,
    };
    KZG::compute_opening_proof(&ck, kzg_claim, &mut prover_transcript);

    // Verifier side
    let mut verifier_transcript = NativeTranscript::verifier_init_empty(&mut prover_transcript);
    let padding_indicator_array = vec![Fr::one(); small_log_n];
    let mut claim_batcher = mock.claim_batcher();
    let g1_identity = srs_points[0];

    let batch_opening_claim = ShpleminiVerifier::compute_batch_opening_claim_without_zk(
        &padding_indicator_array,
        &mut claim_batcher,
        &u,
        &g1_identity,
        &mut verifier_transcript,
    );

    let pairing_points =
        KZG::reduce_verify_batch_opening_claim(batch_opening_claim, &mut verifier_transcript);

    assert!(
        direct_pairing_check(&pairing_points.p0, &pairing_points.p1, &g2_x),
        "High degree attack (accept) pairing check should pass"
    );
}

/// High degree attack test: prover commits to a random higher degree polynomial.
///
/// Port of C++ `TYPED_TEST(ShpleminiTest, HighDegreeAttackReject)` for BN254.
#[test]
fn shplemini_high_degree_attack_reject() {
    let small_log_n: usize = 3;
    let big_n: usize = 1 << 12;
    let (srs_points, g2_x) = create_test_srs(big_n);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points.clone());

    let poly = Polynomial::random(big_n, big_n, 0);
    let u = random_evaluation_point(small_log_n);
    let claimed_multilinear_eval = Fr::random_element();

    let mock = MockClaimGenerator::with_custom_poly(big_n, poly, claimed_multilinear_eval, &ck);

    let mut prover_transcript = NativeTranscript::prover_init_empty();
    let mut batcher = mock.polynomial_batcher();

    // circuit_size = big_n: Gemini folds over the full polynomial size,
    // but only small_log_n challenges are meaningful. Random polynomial won't fold correctly.
    let opening_claim = ShpleminiProver::prove_without_zk(
        big_n,
        &mut batcher,
        &u,
        &ck,
        &mut prover_transcript,
    );

    // KZG prover
    let kzg_claim = crate::claim::ProverOpeningClaim::<Bn254G1Params> {
        polynomial: opening_claim.polynomial,
        opening_pair: OpeningPair {
            challenge: opening_claim.opening_pair.challenge,
            evaluation: opening_claim.opening_pair.evaluation,
        },
        gemini_fold: false,
    };
    KZG::compute_opening_proof(&ck, kzg_claim, &mut prover_transcript);

    // Verifier side
    let mut verifier_transcript = NativeTranscript::verifier_init_empty(&mut prover_transcript);
    let padding_indicator_array = vec![Fr::one(); small_log_n];
    let mut claim_batcher = mock.claim_batcher();
    let g1_identity = srs_points[0];

    let batch_opening_claim = ShpleminiVerifier::compute_batch_opening_claim_without_zk(
        &padding_indicator_array,
        &mut claim_batcher,
        &u,
        &g1_identity,
        &mut verifier_transcript,
    );

    let pairing_points =
        KZG::reduce_verify_batch_opening_claim(batch_opening_claim, &mut verifier_transcript);

    // Should FAIL - random polynomial doesn't fold correctly
    assert!(
        !direct_pairing_check(&pairing_points.p0, &pairing_points.p1, &g2_x),
        "High degree attack (reject) pairing check should fail"
    );
}
