use super::*;
use bbrs_ecc::curves::bn254::Bn254FrParams;
use bbrs_ecc::fields::field::Field;
use bbrs_polynomials::polynomial::Polynomial;

type Fr = Field<Bn254FrParams>;

// ── Unit tests for free functions ────────────────────────────────────────────

#[test]
fn test_powers_of_rho() {
    let rho = Fr::from(3u64);
    let powers = powers_of_rho(rho, 5);
    assert_eq!(powers.len(), 5);
    assert_eq!(powers[0], Fr::one());
    assert_eq!(powers[1], Fr::from(3u64));
    assert_eq!(powers[2], Fr::from(9u64));
    assert_eq!(powers[3], Fr::from(27u64));
    assert_eq!(powers[4], Fr::from(81u64));
}

#[test]
fn test_powers_of_rho_edge_cases() {
    assert_eq!(powers_of_rho(Fr::from(2u64), 0).len(), 0);
    let one = powers_of_rho(Fr::from(7u64), 1);
    assert_eq!(one.len(), 1);
    assert_eq!(one[0], Fr::one());
}

#[test]
fn test_powers_of_evaluation_challenge() {
    let r = Fr::from(5u64);
    let squares = powers_of_evaluation_challenge(r, 4);
    assert_eq!(squares.len(), 4);
    assert_eq!(squares[0], Fr::from(5u64)); // r
    assert_eq!(squares[1], Fr::from(25u64)); // r^2
    assert_eq!(squares[2], Fr::from(625u64)); // r^4
    assert_eq!(squares[3], Fr::from(390625u64)); // r^8
}

// ── PolynomialBatcher tests ──────────────────────────────────────────────────

#[test]
fn test_batcher_unshifted_only() {
    // Two polynomials: f0 = [1, 2, 3, 4], f1 = [5, 6, 7, 8]
    // With rho = 1, F = f0 + f1 = [6, 8, 10, 12]
    let n = 4;
    let f0 = Polynomial::<Bn254FrParams>::from_coefficients(
        vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)],
        n,
    );
    let f1 = Polynomial::<Bn254FrParams>::from_coefficients(
        vec![Fr::from(5u64), Fr::from(6u64), Fr::from(7u64), Fr::from(8u64)],
        n,
    );

    let mut batcher = PolynomialBatcher::new(n);
    batcher.unshifted = vec![&f0, &f1];

    let rho = Fr::one();
    let a_0 = batcher.compute_batched(rho);

    // With rho=1: running_scalar goes 1, 1, so F = 1*f0 + 1*f1
    assert_eq!(a_0.get(0), Fr::from(6u64));
    assert_eq!(a_0.get(1), Fr::from(8u64));
    assert_eq!(a_0.get(2), Fr::from(10u64));
    assert_eq!(a_0.get(3), Fr::from(12u64));
}

#[test]
fn test_batcher_with_shifted() {
    // f0 = [1, 2, 3, 4], g0 = [0, 20, 30, 40] (g0[0] = 0, shiftable)
    // With rho = 1: F = f0, G = g0
    // G/X = left-shift of G: [G[1], G[2], G[3], 0] = [20, 30, 40, 0]
    // A_0 = F + G/X = [1+20, 2+30, 3+40, 4+0] = [21, 32, 43, 4]
    let n = 4;
    let f0 = Polynomial::<Bn254FrParams>::from_coefficients(
        vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)],
        n,
    );
    let g0 = Polynomial::<Bn254FrParams>::from_coefficients(
        vec![Fr::zero(), Fr::from(20u64), Fr::from(30u64), Fr::from(40u64)],
        n,
    );

    let mut batcher = PolynomialBatcher::new(n);
    batcher.unshifted = vec![&f0];
    batcher.to_be_shifted_by_one = vec![&g0];

    let rho = Fr::one();
    let a_0 = batcher.compute_batched(rho);

    // G/X = [g0[1], g0[2], g0[3], 0] = [20, 30, 40, 0]
    // A_0 = F + G/X
    assert_eq!(a_0.get(0), Fr::from(21u64)); // 1 + 20
    assert_eq!(a_0.get(1), Fr::from(32u64)); // 2 + 30
    assert_eq!(a_0.get(2), Fr::from(43u64)); // 3 + 40
    assert_eq!(a_0.get(3), Fr::from(4u64)); // 4 + 0
}

#[test]
fn test_partially_evaluated_batch_polynomials() {
    let n = 4;
    let f0 = Polynomial::<Bn254FrParams>::from_coefficients(
        vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)],
        n,
    );
    // g0[0] = 0 (shiftable constraint)
    let g0 = Polynomial::<Bn254FrParams>::from_coefficients(
        vec![Fr::zero(), Fr::from(20u64), Fr::from(30u64), Fr::from(40u64)],
        n,
    );

    let mut batcher = PolynomialBatcher::new(n);
    batcher.unshifted = vec![&f0];
    batcher.to_be_shifted_by_one = vec![&g0];

    let rho = Fr::one();
    let _a_0 = batcher.compute_batched(rho);

    let r = Fr::from(7u64);
    let (a_0_pos, a_0_neg) = batcher.compute_partially_evaluated_batch_polynomials(r);

    // A_0+(X) = F(X) + G(X)/r
    // A_0-(X) = F(X) - G(X)/r
    // G here is the unshifted batched polynomial (not G/X)
    let r_inv = r.invert();

    for i in 0..n {
        let f_i = f0.get(i);
        let g_i = g0.get(i); // g0[0] = 0, so no contribution at index 0
        let expected_pos = f_i + g_i * r_inv;
        let expected_neg = f_i - g_i * r_inv;
        assert_eq!(a_0_pos.get(i), expected_pos, "A_0+ mismatch at index {}", i);
        assert_eq!(a_0_neg.get(i), expected_neg, "A_0- mismatch at index {}", i);
    }
}

// ── Fold polynomial tests ────────────────────────────────────────────────────

#[test]
fn test_compute_fold_polynomials_simple() {
    // log_n = 2, n = 4
    // A_0 = [a, b, c, d]
    // u = [u0, u1]
    // After one fold: A_1[j] = (1-u0)*A_0[2j] + u0*A_0[2j+1]
    //   A_1[0] = (1-u0)*a + u0*b = a + u0*(b-a)
    //   A_1[1] = (1-u0)*c + u0*d = c + u0*(d-c)
    let n = 4;
    let log_n = 2;
    let a = Fr::from(1u64);
    let b = Fr::from(2u64);
    let c = Fr::from(3u64);
    let d = Fr::from(4u64);
    let a_0 = Polynomial::<Bn254FrParams>::from_coefficients(vec![a, b, c, d], n);

    let u0 = Fr::from(5u64);
    let u1 = Fr::from(7u64);
    let challenge = vec![u0, u1];

    let folds = GeminiProver::compute_fold_polynomials(log_n, &challenge, &a_0);

    // log_n - 1 = 1 real fold (A_1) + 1 constant fold (A_2) = 2 total
    assert_eq!(folds.len(), 2);

    // A_1[0] = a + u0*(b - a) = 1 + 5*(2-1) = 6
    // A_1[1] = c + u0*(d - c) = 3 + 5*(4-3) = 8
    let a1_0 = a + u0 * (b - a);
    let a1_1 = c + u0 * (d - c);
    assert_eq!(folds[0].get(0), a1_0);
    assert_eq!(folds[0].get(1), a1_1);

    // Final eval: A_1[0] + u1*(A_1[1] - A_1[0])
    let final_eval = a1_0 + u1 * (a1_1 - a1_0);
    assert_eq!(folds[1].get(0), final_eval);
}

#[test]
fn test_compute_fold_polynomials_larger() {
    // log_n = 3, n = 8
    let n = 8;
    let log_n = 3;
    let coeffs: Vec<Fr> = (1..=8).map(|x| Fr::from(x as u64)).collect();
    let a_0 = Polynomial::<Bn254FrParams>::from_coefficients(coeffs.clone(), n);

    let u0 = Fr::from(2u64);
    let u1 = Fr::from(3u64);
    let u2 = Fr::from(5u64);
    let challenge = vec![u0, u1, u2];

    let folds = GeminiProver::compute_fold_polynomials(log_n, &challenge, &a_0);

    // log_n - 1 = 2 real folds + 1 const = 3
    assert_eq!(folds.len(), 3);

    // A_1 has 4 coefficients: A_1[j] = A_0[2j] + u0*(A_0[2j+1] - A_0[2j])
    assert_eq!(folds[0].size(), 4);
    for j in 0..4 {
        let even = coeffs[2 * j];
        let odd = coeffs[2 * j + 1];
        let expected = even + u0 * (odd - even);
        assert_eq!(folds[0].get(j), expected, "A_1[{}] mismatch", j);
    }

    // A_2 has 2 coefficients
    assert_eq!(folds[1].size(), 2);
    for j in 0..2 {
        let even = folds[0].get(2 * j);
        let odd = folds[0].get(2 * j + 1);
        let expected = even + u1 * (odd - even);
        assert_eq!(folds[1].get(j), expected, "A_2[{}] mismatch", j);
    }

    // FOLD_3 is the constant: A_2[0] + u2*(A_2[1] - A_2[0])
    let final_eval = folds[1].get(0) + u2 * (folds[1].get(1) - folds[1].get(0));
    assert_eq!(folds[2].get(0), final_eval);
}

// ── Helper: create a random "shiftable" polynomial (g[0] = 0) ───────────────

fn random_shiftable_polynomial(n: usize) -> Polynomial<Bn254FrParams> {
    let mut coeffs = vec![Fr::zero()]; // g[0] = 0
    for _ in 1..n {
        coeffs.push(Fr::random_element());
    }
    Polynomial::from_coefficients(coeffs, n)
}

// ── Full fold/unfold roundtrip ───────────────────────────────────────────────

/// The core Gemini fold/unfold roundtrip test.
///
/// 1. Create random multilinear polynomials
/// 2. Evaluate them at a random multilinear point u
/// 3. Run Gemini prover: batch, fold, compute opening claims
/// 4. Run Gemini verifier: reconstruct positive evaluations from negative ones
/// 5. Verify that the verifier's positive evaluations match the prover's claims
#[test]
fn test_gemini_fold_unfold_roundtrip() {
    let log_n = 4;
    let n: usize = 1 << log_n;

    // Random multilinear polynomial (just one, unshifted for simplicity)
    let poly = Polynomial::<Bn254FrParams>::random(n, n, 0);

    // Random multilinear evaluation point
    let u: Vec<Fr> = (0..log_n).map(|_| Fr::random_element()).collect();

    // Evaluate the multilinear polynomial at u
    let eval_at_u = poly.evaluate_mle(&u, false);

    // Challenges
    let rho = Fr::random_element();
    let r = Fr::random_element();

    // --- Prover side ---
    let mut batcher = PolynomialBatcher::new(n);
    batcher.unshifted = vec![&poly];

    // Compute A_0
    let a_0 = batcher.compute_batched(rho);

    // Compute fold polynomials
    let fold_polynomials = GeminiProver::compute_fold_polynomials(log_n, &u, &a_0);

    // Compute A_0+(X) and A_0-(X)
    let (a_0_pos, a_0_neg) =
        batcher.compute_partially_evaluated_batch_polynomials(r);

    // Construct claims
    let claims = GeminiProver::construct_univariate_opening_claims(
        log_n, a_0_pos, a_0_neg, fold_polynomials, r,
    );

    // --- Verifier side ---
    // batched_evaluation = rho^0 * eval_at_u = eval_at_u (single poly, rho^0 = 1)
    let batched_evaluation = eval_at_u;

    // challenge_powers: r, r^2, r^4, ...
    let challenge_powers = powers_of_evaluation_challenge(r, log_n);

    // fold_neg_evals: A_0(-r), A_1(-r^2), ..., A_{d-1}(-r^{2^{d-1}})
    let mut fold_neg_evals = Vec::with_capacity(log_n);
    fold_neg_evals.push(claims[1].opening_pair.evaluation); // A_0(-r)
    for i in 2..claims.len() {
        fold_neg_evals.push(claims[i].opening_pair.evaluation);
    }

    // Padding indicator: all 1s (no padding, virtual_log_n == log_n)
    let padding_indicator: Vec<Fr> = vec![Fr::one(); log_n];

    // Compute positive evaluations
    let fold_pos_evals = GeminiVerifier::compute_fold_pos_evaluations(
        &padding_indicator,
        batched_evaluation,
        &u,
        &challenge_powers,
        &fold_neg_evals,
        Fr::zero(),
    );

    // --- Verify roundtrip ---
    // fold_pos_evals[0] should equal A_0+(r) = claims[0].evaluation
    assert_eq!(
        fold_pos_evals[0],
        claims[0].opening_pair.evaluation,
        "A_0+(r) mismatch: verifier reconstruction != prover evaluation"
    );

    // fold_pos_evals[l] should be A_l(r^{2^l}) for l = 1, ..., log_n-1
    for l in 1..log_n {
        let r_sq_l = challenge_powers[l];
        let prover_eval_pos = claims[l + 1].polynomial.evaluate(&r_sq_l);
        assert_eq!(
            fold_pos_evals[l], prover_eval_pos,
            "A_{}(r^{{2^{}}}) mismatch at l={}",
            l, l, l
        );
    }
}

/// Roundtrip test with both unshifted and to-be-shifted polynomials.
#[test]
fn test_gemini_roundtrip_with_shifted() {
    let log_n = 3;
    let n: usize = 1 << log_n;

    let f = Polynomial::<Bn254FrParams>::random(n, n, 0);
    // g must have g[0] = 0 (shiftable constraint)
    let g = random_shiftable_polynomial(n);

    let u: Vec<Fr> = (0..log_n).map(|_| Fr::random_element()).collect();

    // Evaluations
    let f_eval = f.evaluate_mle(&u, false);
    let g_eval = g.evaluate_mle(&u, true); // shifted MLE evaluation

    let rho = Fr::random_element();
    let r = Fr::random_element();

    // Batched evaluation: rho^0 * f(u) + rho^1 * g_shift(u)
    let batched_evaluation = f_eval + rho * g_eval;

    // --- Prover ---
    let mut batcher = PolynomialBatcher::new(n);
    batcher.unshifted = vec![&f];
    batcher.to_be_shifted_by_one = vec![&g];

    let a_0 = batcher.compute_batched(rho);

    let fold_polynomials = GeminiProver::compute_fold_polynomials(log_n, &u, &a_0);

    let (a_0_pos, a_0_neg) =
        batcher.compute_partially_evaluated_batch_polynomials(r);

    let claims = GeminiProver::construct_univariate_opening_claims(
        log_n, a_0_pos, a_0_neg, fold_polynomials, r,
    );

    // --- Verifier ---
    let challenge_powers = powers_of_evaluation_challenge(r, log_n);

    let mut fold_neg_evals = Vec::with_capacity(log_n);
    fold_neg_evals.push(claims[1].opening_pair.evaluation);
    for i in 2..claims.len() {
        fold_neg_evals.push(claims[i].opening_pair.evaluation);
    }

    let padding_indicator: Vec<Fr> = vec![Fr::one(); log_n];

    let fold_pos_evals = GeminiVerifier::compute_fold_pos_evaluations(
        &padding_indicator,
        batched_evaluation,
        &u,
        &challenge_powers,
        &fold_neg_evals,
        Fr::zero(),
    );

    // Verify A_0+(r) matches
    assert_eq!(
        fold_pos_evals[0],
        claims[0].opening_pair.evaluation,
        "A_0+(r) mismatch with shifted polynomials"
    );

    // Verify fold evaluations match
    for l in 1..log_n {
        let r_sq_l = challenge_powers[l];
        let prover_eval_pos = claims[l + 1].polynomial.evaluate(&r_sq_l);
        assert_eq!(
            fold_pos_evals[l], prover_eval_pos,
            "A_{}(r^{{2^{}}}) mismatch (shifted test)",
            l, l
        );
    }
}

/// Test with multiple unshifted polynomials batched together.
#[test]
fn test_gemini_roundtrip_multiple_polys() {
    let log_n = 3;
    let n: usize = 1 << log_n;

    let f0 = Polynomial::<Bn254FrParams>::random(n, n, 0);
    let f1 = Polynomial::<Bn254FrParams>::random(n, n, 0);
    let f2 = Polynomial::<Bn254FrParams>::random(n, n, 0);

    let u: Vec<Fr> = (0..log_n).map(|_| Fr::random_element()).collect();

    let f0_eval = f0.evaluate_mle(&u, false);
    let f1_eval = f1.evaluate_mle(&u, false);
    let f2_eval = f2.evaluate_mle(&u, false);

    let rho = Fr::random_element();
    let r = Fr::random_element();

    // Batched evaluation: rho^0 * f0(u) + rho^1 * f1(u) + rho^2 * f2(u)
    let rho_powers = powers_of_rho(rho, 3);
    let batched_evaluation =
        rho_powers[0] * f0_eval + rho_powers[1] * f1_eval + rho_powers[2] * f2_eval;

    // --- Prover ---
    let mut batcher = PolynomialBatcher::new(n);
    batcher.unshifted = vec![&f0, &f1, &f2];

    let a_0 = batcher.compute_batched(rho);
    let fold_polynomials = GeminiProver::compute_fold_polynomials(log_n, &u, &a_0);
    let (a_0_pos, a_0_neg) =
        batcher.compute_partially_evaluated_batch_polynomials(r);
    let claims = GeminiProver::construct_univariate_opening_claims(
        log_n, a_0_pos, a_0_neg, fold_polynomials, r,
    );

    // --- Verifier ---
    let challenge_powers = powers_of_evaluation_challenge(r, log_n);

    let mut fold_neg_evals = Vec::with_capacity(log_n);
    fold_neg_evals.push(claims[1].opening_pair.evaluation);
    for i in 2..claims.len() {
        fold_neg_evals.push(claims[i].opening_pair.evaluation);
    }

    let padding_indicator: Vec<Fr> = vec![Fr::one(); log_n];

    let fold_pos_evals = GeminiVerifier::compute_fold_pos_evaluations(
        &padding_indicator,
        batched_evaluation,
        &u,
        &challenge_powers,
        &fold_neg_evals,
        Fr::zero(),
    );

    assert_eq!(
        fold_pos_evals[0],
        claims[0].opening_pair.evaluation,
        "A_0+(r) mismatch with multiple polys"
    );

    for l in 1..log_n {
        let r_sq_l = challenge_powers[l];
        let prover_eval_pos = claims[l + 1].polynomial.evaluate(&r_sq_l);
        assert_eq!(
            fold_pos_evals[l], prover_eval_pos,
            "A_{}(r^{{2^{}}}) mismatch (multi poly)",
            l, l
        );
    }
}

/// Test with multiple unshifted and shifted polynomials.
#[test]
fn test_gemini_roundtrip_mixed() {
    let log_n = 4;
    let n: usize = 1 << log_n;

    let f0 = Polynomial::<Bn254FrParams>::random(n, n, 0);
    let f1 = Polynomial::<Bn254FrParams>::random(n, n, 0);
    let g0 = random_shiftable_polynomial(n);
    let g1 = random_shiftable_polynomial(n);

    let u: Vec<Fr> = (0..log_n).map(|_| Fr::random_element()).collect();

    let rho = Fr::random_element();
    let r = Fr::random_element();

    // Batched evaluation: rho^0*f0(u) + rho^1*f1(u) + rho^2*g0_shift(u) + rho^3*g1_shift(u)
    let rho_powers = powers_of_rho(rho, 4);
    let batched_evaluation = rho_powers[0] * f0.evaluate_mle(&u, false)
        + rho_powers[1] * f1.evaluate_mle(&u, false)
        + rho_powers[2] * g0.evaluate_mle(&u, true)
        + rho_powers[3] * g1.evaluate_mle(&u, true);

    // --- Prover ---
    let mut batcher = PolynomialBatcher::new(n);
    batcher.unshifted = vec![&f0, &f1];
    batcher.to_be_shifted_by_one = vec![&g0, &g1];

    let a_0 = batcher.compute_batched(rho);
    let fold_polynomials = GeminiProver::compute_fold_polynomials(log_n, &u, &a_0);
    let (a_0_pos, a_0_neg) =
        batcher.compute_partially_evaluated_batch_polynomials(r);
    let claims = GeminiProver::construct_univariate_opening_claims(
        log_n, a_0_pos, a_0_neg, fold_polynomials, r,
    );

    // --- Verifier ---
    let challenge_powers = powers_of_evaluation_challenge(r, log_n);

    let mut fold_neg_evals = Vec::with_capacity(log_n);
    fold_neg_evals.push(claims[1].opening_pair.evaluation);
    for i in 2..claims.len() {
        fold_neg_evals.push(claims[i].opening_pair.evaluation);
    }

    let padding_indicator: Vec<Fr> = vec![Fr::one(); log_n];

    let fold_pos_evals = GeminiVerifier::compute_fold_pos_evaluations(
        &padding_indicator,
        batched_evaluation,
        &u,
        &challenge_powers,
        &fold_neg_evals,
        Fr::zero(),
    );

    assert_eq!(
        fold_pos_evals[0],
        claims[0].opening_pair.evaluation,
        "A_0+(r) mismatch with mixed polys"
    );

    for l in 1..log_n {
        let r_sq_l = challenge_powers[l];
        let prover_eval_pos = claims[l + 1].polynomial.evaluate(&r_sq_l);
        assert_eq!(
            fold_pos_evals[l], prover_eval_pos,
            "A_{}(r^{{2^{}}}) mismatch (mixed test)",
            l, l
        );
    }
}
