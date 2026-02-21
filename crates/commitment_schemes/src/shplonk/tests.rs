//! Tests for Shplonk and Shplemini commitment schemes.
//!
//! Port of C++ shplonk.test.cpp.

use bbrs_ecc::curves::bn254::{
    Bn254FrParams, Bn254G1Params, Fq12, Fr, G1Affine as Bn254G1Affine, G2AffineElement, G2Element,
};
use bbrs_ecc::curves::bn254_pairing::reduced_ate_pairing;
use bbrs_ecc::groups::element::Element;
use bbrs_polynomials::polynomial::Polynomial;

use crate::claim::{OpeningClaim, OpeningPair};
use crate::commitment_key::CommitmentKey;
use crate::kzg::KZG;

use super::{LinearCombinationOfClaims, ProverOpeningClaim, ShplonkProver, ShplonkVerifier};

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
