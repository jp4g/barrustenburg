//! Tests for commitment schemes.

use bbrs_ecc::curves::bn254::{
    Bn254G1Params, Fq12, Fr, G1Affine as Bn254G1Affine, G2AffineElement, G2Element,
};
use bbrs_ecc::curves::bn254_pairing::reduced_ate_pairing;
use bbrs_ecc::groups::element::Element;
use bbrs_ecc::scalar_multiplication::naive_msm;
use bbrs_polynomials::polynomial::Polynomial;

use crate::batch_mul::batch_mul_native;
use crate::claim::{OpeningClaim, OpeningPair, ProverOpeningClaim};
use crate::commitment_key::CommitmentKey;
use crate::kzg::KZG;
use crate::verification_key::Bn254VerifierCommitmentKey;

use bbrs_transcript::NativeTranscript;

/// Generate n random BN254 G1 points for testing.
fn random_bn254_points(n: usize) -> Vec<Bn254G1Affine> {
    let generator = Element::<Bn254G1Params>::one();
    let mut points = Vec::with_capacity(n);
    let mut acc = generator;
    for _ in 0..n {
        points.push(acc.to_affine());
        let scalar = Fr::random_element();
        acc = acc.mul(&scalar);
    }
    points
}

#[test]
fn commitment_key_commit_matches_naive_msm() {
    let n = 8;
    let srs_points = random_bn254_points(n);

    // Use from_points to avoid global CRS singleton issues in parallel tests
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points.clone());

    // Create a random polynomial
    let coeffs: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let poly = Polynomial::from_coefficients(coeffs.clone(), n);

    // Commit using the commitment key
    let commitment = ck.commit(&poly);

    // Compute naive MSM for comparison
    let naive_result = naive_msm::<Bn254G1Params>(&coeffs, &srs_points[..n]).to_affine();

    assert_eq!(
        commitment, naive_result,
        "CommitmentKey::commit must match naive MSM"
    );
}

#[test]
fn commitment_key_commit_smaller_poly() {
    let srs_size = 16;
    let poly_size = 5;
    let srs_points = random_bn254_points(srs_size);

    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points.clone());

    let coeffs: Vec<Fr> = (0..poly_size).map(|_| Fr::random_element()).collect();
    let poly = Polynomial::from_coefficients(coeffs.clone(), srs_size);

    let commitment = ck.commit(&poly);
    let naive_result =
        naive_msm::<Bn254G1Params>(&coeffs, &srs_points[..poly_size]).to_affine();

    assert_eq!(
        commitment, naive_result,
        "Commitment of polynomial smaller than SRS must match naive MSM"
    );
}

#[test]
fn batch_mul_native_matches_naive_msm() {
    let n = 8;
    let points = random_bn254_points(n);
    let scalars: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();

    let batch_result = batch_mul_native::<Bn254G1Params>(&points, &scalars);
    let naive_result = naive_msm::<Bn254G1Params>(&scalars, &points).to_affine();

    assert_eq!(
        batch_result, naive_result,
        "batch_mul_native must match naive MSM"
    );
}

#[test]
fn opening_claim_verify_correct() {
    let n = 8;
    let srs_points = random_bn254_points(n);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points);

    // Create a random polynomial
    let coeffs: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let poly = Polynomial::from_coefficients(coeffs, n);

    // Choose a random evaluation point
    let challenge = Fr::random_element();
    let evaluation = poly.evaluate(&challenge);

    // Compute the commitment
    let commitment = ck.commit(&poly);

    // Create and verify the claim
    let claim = OpeningClaim::<Bn254G1Params> {
        opening_pair: OpeningPair {
            challenge,
            evaluation,
        },
        commitment,
    };

    assert!(
        claim.verify(&ck, &poly),
        "Valid opening claim must verify"
    );
}

#[test]
fn opening_claim_verify_wrong_evaluation() {
    let n = 8;
    let srs_points = random_bn254_points(n);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points);

    let coeffs: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let poly = Polynomial::from_coefficients(coeffs, n);

    let challenge = Fr::random_element();
    let commitment = ck.commit(&poly);

    // Wrong evaluation
    let claim = OpeningClaim::<Bn254G1Params> {
        opening_pair: OpeningPair {
            challenge,
            evaluation: Fr::random_element(),
        },
        commitment,
    };

    assert!(
        !claim.verify(&ck, &poly),
        "Opening claim with wrong evaluation must not verify"
    );
}

#[test]
fn verifier_commitment_key_pairing_check() {
    // The pairing check verifies: e(P0, [1]_2) * e(P1, [x]_2) == 1_T
    //
    // A trivial valid pair: P0 = tau * G1, P1 = -G1
    // Then e(tau*G1, G2) * e(-G1, tau*G2) = e(G1,G2)^tau * e(G1,G2)^{-tau} = 1
    let tau = Fr::from(42u64);

    // G2_x = tau * G2_gen
    use bbrs_ecc::curves::bn254::G2Element;
    let g2_gen = G2AffineElement::generator();
    let g2_proj = G2Element::from_affine(&g2_gen);
    let g2_x_proj = g2_proj.mul_scalar(&tau);
    let g2_x = g2_x_proj.to_affine();

    // P0 = tau * G1_gen
    let g1_gen = Element::<Bn254G1Params>::one();
    let p0 = g1_gen.mul(&tau);

    // P1 = -G1_gen
    let p1 = -g1_gen;

    // Initialize global CRS (needed by verifier key's initialize())
    let crs_points = random_bn254_points(4);
    bbrs_srs::global_crs::init_bn254_mem_crs_factory(&crs_points);

    let mut vkey = Bn254VerifierCommitmentKey::with_g2x(g2_x);

    assert!(
        vkey.pairing_check(&p0, &p1),
        "Valid pairing check must pass: e(tau*G1, G2) * e(-G1, tau*G2) == 1"
    );
}

// ── KZG tests ────────────────────────────────────────────────────────────────

const KZG_N: usize = 16;

/// Generate a powers-of-tau SRS for testing KZG.
///
/// Returns (G1 SRS points, G2_x = tau * G2_gen).
fn create_kzg_test_srs(n: usize) -> (Vec<Bn254G1Affine>, G2AffineElement) {
    let tau = Fr::random_element();

    // G1 SRS: [G, tau*G, tau^2*G, ..., tau^(n-1)*G]
    let g1 = Element::<Bn254G1Params>::one();
    let mut points = Vec::with_capacity(n);
    let mut tau_power = Fr::one();
    for _ in 0..n {
        points.push(g1.mul(&tau_power).to_affine());
        tau_power = tau_power * tau;
    }

    // G2 SRS: tau * G2_gen
    let g2 = G2Element::from_affine(&G2AffineElement::generator());
    let g2_x = g2.mul_scalar(&tau).to_affine();

    (points, g2_x)
}

/// Direct pairing check: e(P0, [1]_2) * e(P1, [x]_2) == 1_T.
///
/// Bypasses the global CRS to avoid OnceLock contention in tests.
fn direct_pairing_check(
    p0: &Bn254G1Affine,
    p1: &Bn254G1Affine,
    g2_x: &G2AffineElement,
) -> bool {
    let e0 = reduced_ate_pairing(p0, &G2AffineElement::generator());
    let e1 = reduced_ate_pairing(p1, g2_x);
    (e0 * e1) == Fq12::one()
}

/// Prove and verify a KZG opening claim.
fn kzg_prove_and_verify(
    ck: &CommitmentKey<Bn254G1Params>,
    g2_x: &G2AffineElement,
    challenge: Fr,
    evaluation: Fr,
    witness: &Polynomial<bbrs_ecc::curves::bn254::Bn254FrParams>,
) {
    let commitment = ck.commit(witness);
    let opening_claim = OpeningClaim::<Bn254G1Params> {
        opening_pair: OpeningPair {
            challenge,
            evaluation,
        },
        commitment,
    };

    let mut prover_transcript = NativeTranscript::prover_init_empty();

    KZG::compute_opening_proof(
        ck,
        ProverOpeningClaim::<Bn254G1Params> {
            polynomial: witness.clone(),
            opening_pair: OpeningPair {
                challenge,
                evaluation,
            },
            gemini_fold: false,
        },
        &mut prover_transcript,
    );

    let mut verifier_transcript =
        NativeTranscript::verifier_init_empty(&mut prover_transcript);
    let pairing_points = KZG::reduce_verify(&opening_claim, &mut verifier_transcript);

    assert!(
        direct_pairing_check(&pairing_points.p0, &pairing_points.p1, g2_x),
        "KZG pairing check failed"
    );
}

#[test]
fn kzg_single() {
    let (srs_points, g2_x) = create_kzg_test_srs(KZG_N);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points);

    let witness = Polynomial::random(KZG_N, KZG_N, 0);
    let challenge = Fr::random_element();
    let evaluation = witness.evaluate(&challenge);

    kzg_prove_and_verify(&ck, &g2_x, challenge, evaluation, &witness);
}

#[test]
fn kzg_zero_evaluation() {
    let (srs_points, g2_x) = create_kzg_test_srs(KZG_N);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points);

    let mut witness = Polynomial::random(KZG_N, KZG_N, 0);
    let challenge = Fr::random_element();
    let evaluation = witness.evaluate(&challenge);

    // Modify witness to achieve zero evaluation: p(r) - p(r) = 0
    *witness.at_mut(0) = witness.get(0) - evaluation;

    kzg_prove_and_verify(&ck, &g2_x, challenge, Fr::zero(), &witness);
}

#[test]
fn kzg_zero_polynomial() {
    let poly_size = 10;
    let (srs_points, g2_x) = create_kzg_test_srs(KZG_N);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points);

    let zero = Polynomial::new(poly_size, KZG_N, 0);
    assert!(zero.is_zero());

    let challenge = Fr::random_element();
    let evaluation = zero.evaluate(&challenge);

    kzg_prove_and_verify(&ck, &g2_x, challenge, evaluation, &zero);
}

#[test]
fn kzg_constant_polynomial() {
    let (srs_points, g2_x) = create_kzg_test_srs(KZG_N);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points);

    let constant = Polynomial::random(1, KZG_N, 0);
    let challenge = Fr::random_element();
    let evaluation = constant.evaluate(&challenge);

    kzg_prove_and_verify(&ck, &g2_x, challenge, evaluation, &constant);
}

#[test]
fn kzg_empty_polynomial() {
    let (srs_points, g2_x) = create_kzg_test_srs(KZG_N);
    let ck = CommitmentKey::<Bn254G1Params>::from_points(srs_points);

    let empty = Polynomial::new(0, 0, 0);
    let challenge = Fr::random_element();
    let evaluation = empty.evaluate(&challenge);

    kzg_prove_and_verify(&ck, &g2_x, challenge, evaluation, &empty);
}
