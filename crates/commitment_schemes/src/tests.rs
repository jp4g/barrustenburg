//! Tests for commitment schemes.

use bbrs_ecc::curves::bn254::{Bn254G1Params, Fr, G1Affine as Bn254G1Affine, G2AffineElement};
use bbrs_ecc::groups::element::Element;
use bbrs_ecc::scalar_multiplication::naive_msm;
use bbrs_polynomials::polynomial::Polynomial;

use crate::batch_mul::batch_mul_native;
use crate::claim::{OpeningClaim, OpeningPair};
use crate::commitment_key::CommitmentKey;
use crate::verification_key::Bn254VerifierCommitmentKey;

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
