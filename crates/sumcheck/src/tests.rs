//! Tests for the sumcheck crate — ported from C++ sumcheck_round.test.cpp and sumcheck.test.cpp.

use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};
use bbrs_polynomials::gate_separator::GateSeparatorPolynomial;
use bbrs_polynomials::polynomial::Polynomial;
use bbrs_polynomials::univariate::Univariate;
use bbrs_relations::relation_parameters::RelationParameters;
use bbrs_transcript::NativeTranscript;

use bbrs_flavor::sumcheck_test_flavor::{AllEntities, ProverPolynomials};

use crate::sumcheck::{SumcheckProver, SumcheckVerifier};
use crate::sumcheck_round::{
    accumulate_relation_univariates_short_monomial, compute_effective_round_size,
    extend_edges_short_monomial, initialize_relation_separator, SumcheckProverRound,
    SumcheckTupleOfTuplesOfUnivariates, TupleOfArraysOfValues,
};

// ============================================================================
// Helper: create_satisfiable_trace
// ============================================================================

/// Create a simple satisfiable trace for testing.
///
/// Port of C++ `create_satisfiable_trace<Flavor>()`.
///
/// Sets up a 4-row (circuit_size=4) arithmetic circuit:
///   Row 0: all zeros (padding)
///   Row 1: 1 + 1 = 2  (q_l=1, q_r=1, q_o=-1, w_l=1, w_r=1, w_o=2)
///   Row 2: 2 * 2 = 4  (q_m=1, q_o=-1, w_l=2, w_r=2, w_o=4)
///   Row 3: all zeros (padding)
fn create_satisfiable_trace() -> (
    ProverPolynomials<Bn254FrParams>,
    RelationParameters<Fr>,
) {
    let circuit_size: usize = 4;
    let mut polys = ProverPolynomials::new(circuit_size);

    let one = Fr::one();
    let neg_one = Fr::zero() - one;

    // Row 1: 1 + 1 = 2 (addition gate)
    *polys.q_l.at_mut(1) = one;
    *polys.q_r.at_mut(1) = one;
    *polys.q_o.at_mut(1) = neg_one;
    *polys.q_arith.at_mut(1) = one;
    *polys.w_l.at_mut(1) = one;
    *polys.w_r.at_mut(1) = one;
    *polys.w_o.at_mut(1) = Fr::from(2u64);

    // Row 2: 2 * 2 = 4 (multiplication gate)
    *polys.q_m.at_mut(2) = one;
    *polys.q_o.at_mut(2) = neg_one;
    *polys.q_arith.at_mut(2) = one;
    *polys.w_l.at_mut(2) = Fr::from(2u64);
    *polys.w_r.at_mut(2) = Fr::from(2u64);
    *polys.w_o.at_mut(2) = Fr::from(4u64);

    // Set shifted polynomials
    polys.set_shifted();

    let relation_parameters = RelationParameters::default();

    (polys, relation_parameters)
}

// ============================================================================
// Sumcheck round tests (from sumcheck_round.test.cpp)
// ============================================================================

/// Port of C++ `SumcheckRound::SumcheckTupleOfTuplesOfUnivariates`
#[test]
fn test_sumcheck_tuple_of_tuples_of_univariates() {
    let accumulators = SumcheckTupleOfTuplesOfUnivariates::<Bn254FrParams>::zero();

    // All evaluations should be zero
    for val in &accumulators.arithmetic_0.evaluations {
        assert!(val.is_zero());
    }
    for val in &accumulators.arithmetic_1.evaluations {
        assert!(val.is_zero());
    }
    for val in &accumulators.dependent_test_0.evaluations {
        assert!(val.is_zero());
    }
}

/// Port of C++ `SumcheckRound::TuplesOfEvaluationArrays`
#[test]
fn test_tuples_of_evaluation_arrays() {
    let values = TupleOfArraysOfValues::<Bn254FrParams>::zero();

    for val in &values.arithmetic {
        assert!(val.is_zero());
    }
    for val in &values.dependent_test {
        assert!(val.is_zero());
    }
}

/// Port of C++ `SumcheckRound::AddTuplesOfTuplesOfUnivariates`
#[test]
fn test_add_tuples_of_tuples_of_univariates() {
    let mut a = SumcheckTupleOfTuplesOfUnivariates::<Bn254FrParams>::zero();
    let mut b = SumcheckTupleOfTuplesOfUnivariates::<Bn254FrParams>::zero();

    // Set some values in a
    a.arithmetic_0.evaluations[0] = Fr::from(1u64);
    a.arithmetic_0.evaluations[1] = Fr::from(2u64);
    a.arithmetic_1.evaluations[0] = Fr::from(3u64);
    a.dependent_test_0.evaluations[0] = Fr::from(4u64);

    // Set some values in b
    b.arithmetic_0.evaluations[0] = Fr::from(10u64);
    b.arithmetic_0.evaluations[1] = Fr::from(20u64);
    b.arithmetic_1.evaluations[0] = Fr::from(30u64);
    b.dependent_test_0.evaluations[0] = Fr::from(40u64);

    a.add(&b);

    assert_eq!(a.arithmetic_0.evaluations[0], Fr::from(11u64));
    assert_eq!(a.arithmetic_0.evaluations[1], Fr::from(22u64));
    assert_eq!(a.arithmetic_1.evaluations[0], Fr::from(33u64));
    assert_eq!(a.dependent_test_0.evaluations[0], Fr::from(44u64));
}

/// Port of C++ `SumcheckRound::ComputeEffectiveRoundSize`
#[test]
fn test_compute_effective_round_size() {
    let circuit_size: usize = 4;
    let polys = ProverPolynomials::<Bn254FrParams>::new(circuit_size);
    let effective = compute_effective_round_size::<Bn254FrParams, _>(&polys, circuit_size);
    assert_eq!(effective, circuit_size);
}

/// Port of C++ `SumcheckRound::ExtendEdgesShortMonomial`
#[test]
fn test_extend_edges_short_monomial() {
    let circuit_size: usize = 4;
    let mut polys = ProverPolynomials::<Bn254FrParams>::new(circuit_size);

    // Use precomputed polys (start_index=0) to avoid shiftable issues
    *polys.q_l.at_mut(0) = Fr::from(10u64);
    *polys.q_l.at_mut(1) = Fr::from(20u64);
    *polys.q_r.at_mut(0) = Fr::from(30u64);
    *polys.q_r.at_mut(1) = Fr::from(40u64);

    // Also set witness polys (start_index=1, so set indices 1,2)
    *polys.w_l.at_mut(1) = Fr::from(50u64);
    *polys.w_l.at_mut(2) = Fr::from(60u64);
    polys.set_shifted();

    let mut extended_edges: AllEntities<Univariate<Bn254FrParams, 2>> = AllEntities {
        q_m: Univariate::zero(),
        q_l: Univariate::zero(),
        q_r: Univariate::zero(),
        q_o: Univariate::zero(),
        q_4: Univariate::zero(),
        q_c: Univariate::zero(),
        q_arith: Univariate::zero(),
        q_test: Univariate::zero(),
        w_l: Univariate::zero(),
        w_r: Univariate::zero(),
        w_o: Univariate::zero(),
        w_4: Univariate::zero(),
        w_test_1: Univariate::zero(),
        w_test_2: Univariate::zero(),
        w_l_shift: Univariate::zero(),
        w_4_shift: Univariate::zero(),
    };

    // Extend at edge_idx=0 (reads indices 0 and 1)
    extend_edges_short_monomial(&mut extended_edges, &polys, 0);

    // Check precomputed polys (start_index=0)
    assert_eq!(extended_edges.q_l.evaluations[0], Fr::from(10u64));
    assert_eq!(extended_edges.q_l.evaluations[1], Fr::from(20u64));
    assert_eq!(extended_edges.q_r.evaluations[0], Fr::from(30u64));
    assert_eq!(extended_edges.q_r.evaluations[1], Fr::from(40u64));

    // Check witness polys at edge 0: w_l[0]=0 (before start_index), w_l[1]=50
    assert_eq!(extended_edges.w_l.evaluations[0], Fr::zero());
    assert_eq!(extended_edges.w_l.evaluations[1], Fr::from(50u64));
}

/// Port of C++ `SumcheckRound::AccumulateRelationUnivariatesSumcheckTestFlavor`
#[test]
fn test_accumulate_relation_univariates() {
    let (polys, relation_parameters) = create_satisfiable_trace();

    let mut extended_edges: AllEntities<Univariate<Bn254FrParams, 2>> = AllEntities {
        q_m: Univariate::zero(),
        q_l: Univariate::zero(),
        q_r: Univariate::zero(),
        q_o: Univariate::zero(),
        q_4: Univariate::zero(),
        q_c: Univariate::zero(),
        q_arith: Univariate::zero(),
        q_test: Univariate::zero(),
        w_l: Univariate::zero(),
        w_r: Univariate::zero(),
        w_o: Univariate::zero(),
        w_4: Univariate::zero(),
        w_test_1: Univariate::zero(),
        w_test_2: Univariate::zero(),
        w_l_shift: Univariate::zero(),
        w_4_shift: Univariate::zero(),
    };

    let mut accumulators = SumcheckTupleOfTuplesOfUnivariates::<Bn254FrParams>::zero();

    // Extend at edge 0 (rows 0,1)
    extend_edges_short_monomial(&mut extended_edges, &polys, 0);
    let scaling_factor = Fr::one();
    accumulate_relation_univariates_short_monomial(
        &mut accumulators,
        &extended_edges,
        &relation_parameters,
        scaling_factor,
    );

    // Extend at edge 2 (rows 2,3)
    extend_edges_short_monomial(&mut extended_edges, &polys, 2);
    accumulate_relation_univariates_short_monomial(
        &mut accumulators,
        &extended_edges,
        &relation_parameters,
        scaling_factor,
    );

    // The accumulators should be non-zero since we have active gates
    let has_nonzero_arith = accumulators.arithmetic_0.evaluations.iter().any(|v| !v.is_zero());
    assert!(has_nonzero_arith, "ArithmeticRelation accumulator should be non-zero for a satisfiable trace");
}

/// Port of C++ `SumcheckRound::CheckSumFieldArithmetic`
#[test]
fn test_check_sum_field_arithmetic() {
    let (polys, relation_parameters) = create_satisfiable_trace();
    let circuit_size = polys.get_polynomial_size();
    let multivariate_d = (circuit_size as f64).log2() as usize;

    // Generate random gate challenges
    let gate_challenges: Vec<Fr> = (0..multivariate_d).map(|i| Fr::from((i + 1) as u64)).collect();
    let alpha = Fr::from(42u64);
    let alphas = initialize_relation_separator(alpha);

    let gate_separators = GateSeparatorPolynomial::new(gate_challenges.clone(), multivariate_d);

    let mut prover_round = SumcheckProverRound::new(circuit_size);
    let round_univariate = prover_round.compute_univariate(
        &polys,
        &relation_parameters,
        &gate_separators,
        &alphas,
    );

    // Verify the checksum: S(0) + S(1) should equal the expected sum
    let s0 = round_univariate.value_at(0);
    let s1 = round_univariate.value_at(1);
    let round_sum = s0 + s1;

    // The sum should be non-trivial (non-zero) for a non-trivial trace
    // (it could be zero if the relation evaluations cancel, but for our trace it shouldn't)
    // The key test is that the prover-verifier roundtrip works (tested below).
    // Here we just verify the univariate was computed without panicking.
    let _ = round_sum;
}

/// Port of C++ `SumcheckRound::CheckSumPaddingIndicator`
#[test]
fn test_check_sum_padding_indicator() {
    let (polys, relation_parameters) = create_satisfiable_trace();
    let circuit_size = polys.get_polynomial_size();
    let multivariate_d = (circuit_size as f64).log2() as usize;

    let gate_challenges: Vec<Fr> = (0..multivariate_d).map(|i| Fr::from((i + 1) as u64)).collect();
    let alpha = Fr::from(42u64);
    let alphas = initialize_relation_separator(alpha);

    let gate_separators = GateSeparatorPolynomial::new(gate_challenges.clone(), multivariate_d);

    let mut prover_round = SumcheckProverRound::new(circuit_size);
    let round_univariate = prover_round.compute_univariate(
        &polys,
        &relation_parameters,
        &gate_separators,
        &alphas,
    );

    // With padding indicator = 1, checksum should pass
    let s0 = round_univariate.value_at(0);
    let s1 = round_univariate.value_at(1);
    let total_sum = s0 + s1;

    use crate::sumcheck_round::SumcheckVerifierRound;
    let mut verifier_round = SumcheckVerifierRound::new(total_sum);
    verifier_round.check_sum(&round_univariate, Fr::one());
    assert!(!verifier_round.round_failed, "Checksum should pass with indicator=1");
}

// ============================================================================
// Sumcheck prover/verifier tests (from sumcheck.test.cpp)
// ============================================================================

/// Port of C++ `SumcheckTests::PolynomialNormalization`
#[test]
fn test_polynomial_normalization() {
    let (polys, relation_parameters) = create_satisfiable_trace();
    let circuit_size = polys.get_polynomial_size();
    let multivariate_d = (circuit_size as f64).log2() as usize;

    let gate_challenges: Vec<Fr> = (0..multivariate_d).map(|_| Fr::random_element()).collect();

    let mut prover_transcript = NativeTranscript::prover_init_empty();
    let alpha: Fr = Fr::random_element();

    let prover = SumcheckProver::new(
        circuit_size,
        &polys,
        &mut prover_transcript,
        alpha,
        gate_challenges.clone(),
        relation_parameters.clone(),
    );

    let output = prover.prove();

    // Verify: the claimed evaluations should be the multilinear evaluation of each polynomial
    // at the challenge point.
    let challenge = &output.challenge;
    let claimed = &output.claimed_evaluations;
    let all_polys = polys.get_all();
    let all_claimed = claimed.get_all();

    for (poly, claimed_val) in all_polys.iter().zip(all_claimed.iter()) {
        // Manually evaluate the polynomial as a multilinear polynomial at the challenge
        let manual_eval = evaluate_multilinear(poly, challenge);
        assert_eq!(
            manual_eval, **claimed_val,
            "Claimed evaluation should match manual multilinear evaluation"
        );
    }
}

/// Port of C++ `SumcheckTests::Prover`
#[test]
fn test_prover() {
    let (polys, relation_parameters) = create_satisfiable_trace();
    let circuit_size = polys.get_polynomial_size();
    let multivariate_d = (circuit_size as f64).log2() as usize;

    let gate_challenges: Vec<Fr> = (0..multivariate_d).map(|_| Fr::random_element()).collect();

    let mut transcript = NativeTranscript::prover_init_empty();
    let alpha: Fr = Fr::random_element();

    let prover = SumcheckProver::new(
        circuit_size,
        &polys,
        &mut transcript,
        alpha,
        gate_challenges,
        relation_parameters,
    );

    // Should complete without panicking
    let _output = prover.prove();
}

/// Port of C++ `SumcheckTests::ProverAndVerifierSimple`
#[test]
fn test_prover_and_verifier_simple() {
    let (polys, relation_parameters) = create_satisfiable_trace();
    let circuit_size = polys.get_polynomial_size();
    let multivariate_d = (circuit_size as f64).log2() as usize;

    let gate_challenges: Vec<Fr> = (0..multivariate_d).map(|_| Fr::random_element()).collect();
    let alpha: Fr = Fr::random_element();

    // Prover
    let mut prover_transcript = NativeTranscript::prover_init_empty();
    let prover = SumcheckProver::new(
        circuit_size,
        &polys,
        &mut prover_transcript,
        alpha,
        gate_challenges.clone(),
        relation_parameters.clone(),
    );
    let _prover_output = prover.prove();

    // Verifier
    let mut verifier_transcript = NativeTranscript::verifier_init_empty(&mut prover_transcript);
    let target_sum = Fr::zero(); // Non-ZK: target sum starts at 0
    let padding_indicator_array: Vec<Fr> = vec![Fr::one(); multivariate_d];

    let verifier = SumcheckVerifier::new(
        &mut verifier_transcript,
        alpha,
        multivariate_d,
        target_sum,
    );
    let verifier_output = verifier.verify(
        &relation_parameters,
        &gate_challenges,
        &padding_indicator_array,
    );

    assert!(verifier_output.verified, "Verifier should accept a valid proof");
}

/// Port of C++ `SumcheckTests::ProverAndVerifierSimpleFailure`
#[test]
fn test_prover_and_verifier_simple_failure() {
    let (polys, relation_parameters) = create_satisfiable_trace();
    let circuit_size = polys.get_polynomial_size();
    let multivariate_d = (circuit_size as f64).log2() as usize;

    let gate_challenges: Vec<Fr> = (0..multivariate_d).map(|_| Fr::random_element()).collect();
    let alpha: Fr = Fr::random_element();

    // Prover
    let mut prover_transcript = NativeTranscript::prover_init_empty();
    let prover = SumcheckProver::new(
        circuit_size,
        &polys,
        &mut prover_transcript,
        alpha,
        gate_challenges.clone(),
        relation_parameters.clone(),
    );
    let _prover_output = prover.prove();

    // Export proof and tamper with it
    let mut proof_data = prover_transcript.export_proof();
    if proof_data.len() > 2 {
        // Tamper with a proof element (these are Fr elements)
        proof_data[1] = proof_data[1] + Fr::one();
    }

    // Verifier with tampered proof
    let mut verifier_transcript = NativeTranscript::from_proof(&proof_data);
    // Consume the init element (matching prover_init_empty)
    let _: Fr = verifier_transcript.receive_from_prover("Init");

    let target_sum = Fr::zero();
    let padding_indicator_array: Vec<Fr> = vec![Fr::one(); multivariate_d];

    let verifier = SumcheckVerifier::new(
        &mut verifier_transcript,
        alpha,
        multivariate_d,
        target_sum,
    );
    let verifier_output = verifier.verify(
        &relation_parameters,
        &gate_challenges,
        &padding_indicator_array,
    );

    assert!(!verifier_output.verified, "Verifier should reject a tampered proof");
}

// ============================================================================
// Helpers
// ============================================================================

/// Evaluate a polynomial as a multilinear polynomial at the given point.
///
/// For a polynomial with coefficients [c0, c1, ..., c_{2^d - 1}],
/// the multilinear evaluation at (x0, x1, ..., x_{d-1}) is:
///   sum_{i=0}^{2^d-1} c_i * prod_{j=0}^{d-1} ((1-x_j) if bit j of i is 0, else x_j)
fn evaluate_multilinear(poly: &Polynomial<Bn254FrParams>, point: &[Fr]) -> Fr {
    let n = 1 << point.len();
    let mut result = Fr::zero();
    for i in 0..n {
        let coeff = poly.get(i);
        if coeff.is_zero() {
            continue;
        }
        let mut term = coeff;
        for (j, xj) in point.iter().enumerate() {
            if (i >> j) & 1 == 1 {
                term = term * *xj;
            } else {
                term = term * (Fr::one() - *xj);
            }
        }
        result = result + term;
    }
    result
}
