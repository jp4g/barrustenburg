//! Port of `sumcheck_round.hpp` — SumcheckProverRound and SumcheckVerifierRound.
//!
//! Implements the core sumcheck round computations for the SumcheckTestFlavor.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_polynomials::gate_separator::GateSeparatorPolynomial;
use bbrs_polynomials::univariate::Univariate;
use bbrs_relations::relation_parameters::RelationParameters;

use bbrs_flavor::sumcheck_test_flavor::{
    AllEntities, AllValues, ExtendedEdges,
    BATCHED_RELATION_PARTIAL_LENGTH, MAX_PARTIAL_RELATION_LENGTH,
    NUM_SUBRELATIONS,
};

// ============================================================================
// Univariate Accumulators — replaces C++ tuple-of-tuples pattern
// ============================================================================

/// Per-relation univariate accumulators for SumcheckTestFlavor.
///
/// C++ uses `tuple<tuple<Univariate<FF,6>, Univariate<FF,5>>, tuple<Univariate<FF,2>>>`
/// We represent this as a struct with named fields.
pub struct SumcheckTupleOfTuplesOfUnivariates<P: FieldParams> {
    /// ArithmeticRelation subrelation 0: partial length 6
    pub arithmetic_0: Univariate<P, 6>,
    /// ArithmeticRelation subrelation 1: partial length 5
    pub arithmetic_1: Univariate<P, 5>,
    /// DependentTestRelation subrelation 0: partial length 2
    pub dependent_test_0: Univariate<P, 2>,
}

impl<P: FieldParams> SumcheckTupleOfTuplesOfUnivariates<P> {
    /// Create zero-initialized accumulators.
    pub fn zero() -> Self {
        Self {
            arithmetic_0: Univariate::zero(),
            arithmetic_1: Univariate::zero(),
            dependent_test_0: Univariate::zero(),
        }
    }

    /// Add another set of accumulators into this one.
    pub fn add(&mut self, other: &Self) {
        for i in 0..6 {
            self.arithmetic_0.evaluations[i] =
                self.arithmetic_0.evaluations[i] + other.arithmetic_0.evaluations[i];
        }
        for i in 0..5 {
            self.arithmetic_1.evaluations[i] =
                self.arithmetic_1.evaluations[i] + other.arithmetic_1.evaluations[i];
        }
        for i in 0..2 {
            self.dependent_test_0.evaluations[i] =
                self.dependent_test_0.evaluations[i] + other.dependent_test_0.evaluations[i];
        }
    }

    /// Scale accumulators by subrelation separator challenges.
    ///
    /// Port of C++ `RelationUtils::scale_univariates`.
    /// The first subrelation is not scaled; subsequent ones are scaled by challenge powers.
    pub fn scale(&mut self, alphas: &[Field<P>; NUM_SUBRELATIONS - 1]) {
        // arithmetic_0: subrelation 0, no scaling (implicitly *1)
        // arithmetic_1: subrelation 1, scale by alphas[0]
        for e in self.arithmetic_1.evaluations.iter_mut() {
            *e = *e * alphas[0];
        }
        // dependent_test_0: subrelation 2, scale by alphas[1]
        for e in self.dependent_test_0.evaluations.iter_mut() {
            *e = *e * alphas[1];
        }
    }
}

/// Array of subrelation evaluation values (for verifier).
///
/// C++ uses `tuple<array<FF, 2>, array<FF, 1>>`.
pub struct TupleOfArraysOfValues<P: FieldParams> {
    /// ArithmeticRelation: 2 subrelation evaluations.
    pub arithmetic: [Field<P>; 2],
    /// DependentTestRelation: 1 subrelation evaluation.
    pub dependent_test: [Field<P>; 1],
}

impl<P: FieldParams> TupleOfArraysOfValues<P> {
    /// Create zero-initialized values.
    pub fn zero() -> Self {
        Self {
            arithmetic: [Field::zero(); 2],
            dependent_test: [Field::zero(); 1],
        }
    }

    /// Flatten to an array of all subrelation values.
    pub fn to_flat_array(&self) -> [Field<P>; NUM_SUBRELATIONS] {
        [
            self.arithmetic[0],
            self.arithmetic[1],
            self.dependent_test[0],
        ]
    }
}

// ============================================================================
// SumcheckProverRound
// ============================================================================

/// Prover round for the sumcheck protocol.
///
/// Port of C++ `SumcheckProverRound<Flavor>` for SumcheckTestFlavor.
pub struct SumcheckProverRound<P: FieldParams> {
    /// In Round i = 0,...,d-1, equals 2^{d-i}.
    pub round_size: usize,
    /// Per-subrelation univariate accumulators.
    pub univariate_accumulators: SumcheckTupleOfTuplesOfUnivariates<P>,
}

/// Type alias for the subrelation separators array.
pub type SubrelationSeparators<P> = [Field<P>; NUM_SUBRELATIONS - 1];

/// The round univariate: Univariate of length BATCHED_RELATION_PARTIAL_LENGTH.
pub type SumcheckRoundUnivariate<P> = Univariate<P, BATCHED_RELATION_PARTIAL_LENGTH>;

impl<P: FieldParams> SumcheckProverRound<P> {
    /// Create a new prover round.
    pub fn new(initial_round_size: usize) -> Self {
        Self {
            round_size: initial_round_size,
            univariate_accumulators: SumcheckTupleOfTuplesOfUnivariates::zero(),
        }
    }

    /// Extend edges with full barycentric extension to MAX_PARTIAL_RELATION_LENGTH.
    ///
    /// Port of C++ `SumcheckProverRound::extend_edges` with USE_SHORT_MONOMIALS=false.
    pub fn extend_edges_full<T>(
        extended_edges: &mut ExtendedEdges<P, MAX_PARTIAL_RELATION_LENGTH>,
        multivariates: &AllEntities<T>,
        edge_idx: usize,
    ) where
        T: PolynomialAccess<P>,
    {
        let all_multivariates = multivariates.get_all();
        let all_edges = extended_edges.get_all_mut();
        for (edge, poly) in all_edges.into_iter().zip(all_multivariates.iter()) {
            if poly.end_index_value() < edge_idx {
                *edge = Univariate::zero();
            } else {
                let base = Univariate::<P, 2>::new([
                    poly.get_value(edge_idx),
                    poly.get_value(edge_idx + 1),
                ]);
                *edge = base.extend_to();
            }
        }
    }

    /// Compute the univariate for this round (short monomial, non-ZK).
    ///
    /// Port of C++ `compute_univariate_with_row_skipping` simplified for
    /// single-threaded, non-AVM, non-ZK case.
    pub fn compute_univariate<T>(
        &mut self,
        polynomials: &AllEntities<T>,
        relation_parameters: &RelationParameters<Field<P>>,
        gate_separators: &GateSeparatorPolynomial<P>,
        alphas: &SubrelationSeparators<P>,
    ) -> SumcheckRoundUnivariate<P>
    where
        T: PolynomialAccess<P> + HasEndIndex,
    {
        let effective_round_size = compute_effective_round_size::<P, T>(polynomials, self.round_size);

        let mut extended_edges: AllEntities<Univariate<P, 2>> = AllEntities {
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

        let mut edge_idx = 0;
        while edge_idx < effective_round_size {
            extend_edges_short_monomial(&mut extended_edges, polynomials, edge_idx);
            let scaling_factor = gate_separators.at(edge_idx);
            accumulate_relation_univariates_short_monomial(
                &mut self.univariate_accumulators,
                &extended_edges,
                relation_parameters,
                scaling_factor,
            );
            edge_idx += 2;
        }

        // Batch over relations
        self.batch_over_relations(alphas, gate_separators)
    }

    /// Batch the per-subrelation univariate accumulators into a single round univariate.
    ///
    /// Port of C++ `batch_over_relations`.
    pub fn batch_over_relations(
        &mut self,
        alphas: &SubrelationSeparators<P>,
        gate_separators: &GateSeparatorPolynomial<P>,
    ) -> SumcheckRoundUnivariate<P> {
        // Scale univariates by subrelation separators
        self.univariate_accumulators.scale(alphas);

        // Extend and batch
        let result = extend_and_batch_univariates(
            &self.univariate_accumulators,
            gate_separators,
        );

        // Reset accumulators for next round
        self.univariate_accumulators = SumcheckTupleOfTuplesOfUnivariates::zero();

        result
    }
}

// ============================================================================
// Free functions for prover round operations
// ============================================================================

/// Compute effective round size (non-ZK).
///
/// Find the maximum end_index across witness polynomials and round up to next even number.
pub fn compute_effective_round_size<P: FieldParams, T: HasEndIndex>(
    multivariates: &AllEntities<T>,
    round_size: usize,
) -> usize {
    let mut max_end_index: usize = 0;
    for witness_poly in multivariates.get_witness() {
        max_end_index = max_end_index.max(witness_poly.end_index());
    }
    let rounded = max_end_index + (max_end_index % 2);
    rounded.min(round_size)
}

/// Extend edges from multivariates at the given edge index (short monomial version).
///
/// For USE_SHORT_MONOMIALS=true, extended edges are just length-2 univariates.
pub fn extend_edges_short_monomial<P: FieldParams, T: PolynomialAccess<P>>(
    extended_edges: &mut AllEntities<Univariate<P, 2>>,
    multivariates: &AllEntities<T>,
    edge_idx: usize,
) {
    let all_multivariates = multivariates.get_all();
    let all_edges = extended_edges.get_all_mut();
    for (edge, poly) in all_edges.into_iter().zip(all_multivariates.iter()) {
        *edge = Univariate::new([poly.get_value(edge_idx), poly.get_value(edge_idx + 1)]);
    }
}

/// Accumulate relation univariates from extended edges (short monomial version).
pub fn accumulate_relation_univariates_short_monomial<P: FieldParams>(
    univariate_accumulators: &mut SumcheckTupleOfTuplesOfUnivariates<P>,
    extended_edges: &AllEntities<Univariate<P, 2>>,
    relation_parameters: &RelationParameters<Field<P>>,
    scaling_factor: Field<P>,
) {
    let q_arith_skip = extended_edges.q_arith.evaluations[0].is_zero()
        && extended_edges.q_arith.evaluations[1].is_zero();
    if !q_arith_skip {
        arithmetic_accumulate_short_monomial(
            &mut univariate_accumulators.arithmetic_0,
            &mut univariate_accumulators.arithmetic_1,
            extended_edges,
            relation_parameters,
            scaling_factor,
        );
    }

    let q_test_skip = extended_edges.q_test.evaluations[0].is_zero()
        && extended_edges.q_test.evaluations[1].is_zero();
    if !q_test_skip {
        dependent_test_accumulate_short_monomial(
            &mut univariate_accumulators.dependent_test_0,
            extended_edges,
        );
    }
}

// ============================================================================
// extend_and_batch_univariates
// ============================================================================

/// Extend sub-relation univariates to BATCHED_RELATION_PARTIAL_LENGTH and batch them.
///
/// Linearly independent subrelations are multiplied by the pow-polynomial:
///   result += extended * extended_random_polynomial * partial_evaluation_result
/// Linearly dependent subrelations (like DependentTestRelation) are added directly:
///   result += extended
///
/// Port of C++ `SumcheckProverRound::extend_and_batch_univariates`.
fn extend_and_batch_univariates<P: FieldParams>(
    accumulators: &SumcheckTupleOfTuplesOfUnivariates<P>,
    gate_separators: &GateSeparatorPolynomial<P>,
) -> SumcheckRoundUnivariate<P> {
    // Pow-factor: (1 - X) + X * beta_i
    let random_polynomial = Univariate::<P, 2>::new([
        Field::one(),
        gate_separators.current_element(),
    ]);
    let extended_random_polynomial: Univariate<P, BATCHED_RELATION_PARTIAL_LENGTH> =
        random_polynomial.extend_to();

    let mut result = Univariate::<P, BATCHED_RELATION_PARTIAL_LENGTH>::zero();

    // Subrelation 0: ArithmeticRelation[0], linearly independent
    {
        let extended: Univariate<P, BATCHED_RELATION_PARTIAL_LENGTH> =
            accumulators.arithmetic_0.extend_to();
        // Multiply by pow polynomial
        let contrib = univariate_mul(&extended, &extended_random_polynomial);
        let contrib = univariate_mul_scalar(&contrib, gate_separators.partial_evaluation_result);
        result = univariate_add(&result, &contrib);
    }

    // Subrelation 1: ArithmeticRelation[1], linearly independent
    {
        let extended: Univariate<P, BATCHED_RELATION_PARTIAL_LENGTH> =
            accumulators.arithmetic_1.extend_to();
        let contrib = univariate_mul(&extended, &extended_random_polynomial);
        let contrib = univariate_mul_scalar(&contrib, gate_separators.partial_evaluation_result);
        result = univariate_add(&result, &contrib);
    }

    // Subrelation 2: DependentTestRelation[0], linearly DEPENDENT → NOT multiplied by pow
    {
        let extended: Univariate<P, BATCHED_RELATION_PARTIAL_LENGTH> =
            accumulators.dependent_test_0.extend_to();
        result = univariate_add(&result, &extended);
    }

    result
}

// ============================================================================
// Arithmetic Relation accumulation (short monomial version)
// ============================================================================

/// Accumulate ArithmeticRelation contributions into univariate accumulators.
///
/// ArithmeticRelation subrelation 0 (partial length 6):
///   q_arith * (q_m * w_l * w_r + q_l * w_l + q_r * w_r + q_o * w_o + q_4 * w_4 + q_c)
///
/// ArithmeticRelation subrelation 1 (partial length 5):
///   q_arith * (q_arith - 1) * w_4
fn arithmetic_accumulate_short_monomial<P: FieldParams>(
    acc_0: &mut Univariate<P, 6>,
    acc_1: &mut Univariate<P, 5>,
    extended_edges: &AllEntities<Univariate<P, 2>>,
    _relation_parameters: &RelationParameters<Field<P>>,
    scaling_factor: Field<P>,
) {
    // For short monomials (length 2), we need to extend to the required partial lengths
    // to compute the product terms.
    let q_m: Univariate<P, 6> = extended_edges.q_m.extend_to();
    let q_l: Univariate<P, 6> = extended_edges.q_l.extend_to();
    let q_r: Univariate<P, 6> = extended_edges.q_r.extend_to();
    let q_o: Univariate<P, 6> = extended_edges.q_o.extend_to();
    let q_4: Univariate<P, 6> = extended_edges.q_4.extend_to();
    let q_c: Univariate<P, 6> = extended_edges.q_c.extend_to();
    let q_arith: Univariate<P, 6> = extended_edges.q_arith.extend_to();
    let w_l: Univariate<P, 6> = extended_edges.w_l.extend_to();
    let w_r: Univariate<P, 6> = extended_edges.w_r.extend_to();
    let w_o: Univariate<P, 6> = extended_edges.w_o.extend_to();
    let w_4: Univariate<P, 6> = extended_edges.w_4.extend_to();

    // Subrelation 0: q_arith * (q_m * w_l * w_r + q_l * w_l + q_r * w_r + q_o * w_o + q_4 * w_4 + q_c)
    {
        let mut tmp = Univariate::<P, 6>::zero();
        for i in 0..6 {
            tmp.evaluations[i] = q_m.evaluations[i] * w_l.evaluations[i] * w_r.evaluations[i]
                + q_l.evaluations[i] * w_l.evaluations[i]
                + q_r.evaluations[i] * w_r.evaluations[i]
                + q_o.evaluations[i] * w_o.evaluations[i]
                + q_4.evaluations[i] * w_4.evaluations[i]
                + q_c.evaluations[i];
            // Multiply by q_arith and scaling_factor
            tmp.evaluations[i] = tmp.evaluations[i] * q_arith.evaluations[i] * scaling_factor;
        }
        for i in 0..6 {
            acc_0.evaluations[i] = acc_0.evaluations[i] + tmp.evaluations[i];
        }
    }

    // Subrelation 1: q_arith * (q_arith - 1) * w_4
    {
        let q_arith5: Univariate<P, 5> = extended_edges.q_arith.extend_to();
        let w_45: Univariate<P, 5> = extended_edges.w_4.extend_to();
        for i in 0..5 {
            let contrib = q_arith5.evaluations[i]
                * (q_arith5.evaluations[i] - Field::one())
                * w_45.evaluations[i]
                * scaling_factor;
            acc_1.evaluations[i] = acc_1.evaluations[i] + contrib;
        }
    }
}

/// Accumulate DependentTestRelation (short monomial version).
///
/// DependentTestRelation subrelation 0 (partial length 2):
///   q_test * w_test_1
///
/// NOTE: Linearly dependent — NOT multiplied by scaling_factor!
fn dependent_test_accumulate_short_monomial<P: FieldParams>(
    acc: &mut Univariate<P, 2>,
    extended_edges: &AllEntities<Univariate<P, 2>>,
) {
    for i in 0..2 {
        acc.evaluations[i] = acc.evaluations[i]
            + extended_edges.q_test.evaluations[i] * extended_edges.w_test_1.evaluations[i];
    }
}

// ============================================================================
// SumcheckVerifierRound
// ============================================================================

/// Verifier round for the sumcheck protocol.
///
/// Port of C++ `SumcheckVerifierRound<Flavor>` (non-Grumpkin).
pub struct SumcheckVerifierRound<P: FieldParams> {
    /// Current target sum for verification.
    pub target_total_sum: Field<P>,
    /// Whether any round check has failed.
    pub round_failed: bool,
    /// Per-subrelation evaluation accumulators.
    pub relation_evaluations: TupleOfArraysOfValues<P>,
}

impl<P: FieldParams> SumcheckVerifierRound<P> {
    /// Create a new verifier round with the given target sum.
    pub fn new(target_total_sum: Field<P>) -> Self {
        Self {
            target_total_sum,
            round_failed: false,
            relation_evaluations: TupleOfArraysOfValues::zero(),
        }
    }

    /// Check that the round target sum is correct.
    ///
    /// total_sum = (1 - indicator) * target + indicator * (S(0) + S(1))
    /// Verifies target_total_sum == total_sum.
    ///
    /// Port of C++ `SumcheckVerifierRound::check_sum`.
    pub fn check_sum(
        &mut self,
        univariate: &Univariate<P, BATCHED_RELATION_PARTIAL_LENGTH>,
        indicator: Field<P>,
    ) {
        let total_sum = (Field::one() - indicator) * self.target_total_sum
            + indicator * (univariate.value_at(0) + univariate.value_at(1));
        let sumcheck_round_failed = self.target_total_sum != total_sum;
        self.round_failed = self.round_failed || sumcheck_round_failed;
    }

    /// Compute the next target sum: sigma_{i+1} = S^i(u_i).
    ///
    /// Port of C++ `SumcheckVerifierRound::compute_next_target_sum`.
    pub fn compute_next_target_sum(
        &mut self,
        univariate: &Univariate<P, BATCHED_RELATION_PARTIAL_LENGTH>,
        round_challenge: Field<P>,
        indicator: Field<P>,
    ) {
        self.target_total_sum = (Field::one() - indicator) * self.target_total_sum
            + indicator * univariate.evaluate(round_challenge);
    }

    /// Compute the full relation purported value at the final evaluation point.
    ///
    /// Port of C++ `SumcheckVerifierRound::compute_full_relation_purported_value`.
    pub fn compute_full_relation_purported_value(
        &mut self,
        purported_evaluations: &AllValues<P>,
        relation_parameters: &RelationParameters<Field<P>>,
        gate_separators: &GateSeparatorPolynomial<P>,
        alphas: &SubrelationSeparators<P>,
    ) -> Field<P> {
        // Accumulate relation evaluations (without skipping)
        accumulate_relation_evaluations(
            &mut self.relation_evaluations,
            purported_evaluations,
            relation_parameters,
            gate_separators.partial_evaluation_result,
        );

        // Scale and batch
        let flat = self.relation_evaluations.to_flat_array();
        scale_and_batch_with_independence(&flat, alphas)
    }

    /// Perform final verification.
    pub fn perform_final_verification(&self, full_honk_purported_value: Field<P>) -> bool {
        full_honk_purported_value == self.target_total_sum
    }
}

/// Accumulate relation evaluations at a single point (verifier side).
fn accumulate_relation_evaluations<P: FieldParams>(
    relation_evaluations: &mut TupleOfArraysOfValues<P>,
    purported_evaluations: &AllValues<P>,
    _relation_parameters: &RelationParameters<Field<P>>,
    scaling_factor: Field<P>,
) {
    // ArithmeticRelation
    // Subrelation 0: q_arith * (q_m * w_l * w_r + q_l * w_l + q_r * w_r + q_o * w_o + q_4 * w_4 + q_c) * scaling_factor
    let arith_inner = purported_evaluations.q_m * purported_evaluations.w_l * purported_evaluations.w_r
        + purported_evaluations.q_l * purported_evaluations.w_l
        + purported_evaluations.q_r * purported_evaluations.w_r
        + purported_evaluations.q_o * purported_evaluations.w_o
        + purported_evaluations.q_4 * purported_evaluations.w_4
        + purported_evaluations.q_c;
    relation_evaluations.arithmetic[0] = relation_evaluations.arithmetic[0]
        + purported_evaluations.q_arith * arith_inner * scaling_factor;

    // Subrelation 1: q_arith * (q_arith - 1) * w_4 * scaling_factor
    relation_evaluations.arithmetic[1] = relation_evaluations.arithmetic[1]
        + purported_evaluations.q_arith
            * (purported_evaluations.q_arith - Field::one())
            * purported_evaluations.w_4
            * scaling_factor;

    // DependentTestRelation
    // Subrelation 0: q_test * w_test_1 (linearly dependent — NOT multiplied by scaling_factor)
    relation_evaluations.dependent_test[0] =
        relation_evaluations.dependent_test[0] + purported_evaluations.q_test * purported_evaluations.w_test_1;
}

/// Scale and batch elements, respecting linear independence.
///
/// Linearly independent subrelations contribute: value * alpha^k * gate_separator
/// Linearly dependent subrelations contribute: value * alpha^k (no gate_separator already handled)
fn scale_and_batch_with_independence<P: FieldParams>(
    values: &[Field<P>; NUM_SUBRELATIONS],
    alphas: &SubrelationSeparators<P>,
) -> Field<P> {
    // Same as scale_and_batch_elements: first element not scaled, rest by alphas
    let mut result = values[0];
    for i in 1..NUM_SUBRELATIONS {
        result = result + values[i] * alphas[i - 1];
    }
    result
}

// ============================================================================
// Univariate arithmetic helpers
// ============================================================================

/// Pointwise multiply two Univariates.
fn univariate_mul<P: FieldParams, const N: usize>(
    a: &Univariate<P, N>,
    b: &Univariate<P, N>,
) -> Univariate<P, N> {
    let mut result = Univariate::zero();
    for i in 0..N {
        result.evaluations[i] = a.evaluations[i] * b.evaluations[i];
    }
    result
}

/// Pointwise multiply Univariate by scalar.
fn univariate_mul_scalar<P: FieldParams, const N: usize>(
    a: &Univariate<P, N>,
    scalar: Field<P>,
) -> Univariate<P, N> {
    let mut result = Univariate::zero();
    for i in 0..N {
        result.evaluations[i] = a.evaluations[i] * scalar;
    }
    result
}

/// Pointwise add two Univariates.
fn univariate_add<P: FieldParams, const N: usize>(
    a: &Univariate<P, N>,
    b: &Univariate<P, N>,
) -> Univariate<P, N> {
    let mut result = Univariate::zero();
    for i in 0..N {
        result.evaluations[i] = a.evaluations[i] + b.evaluations[i];
    }
    result
}

// ============================================================================
// Trait for polymorphic polynomial access
// ============================================================================

/// Trait for types that provide end_index (for compute_effective_round_size).
pub trait HasEndIndex {
    fn end_index(&self) -> usize;
}

/// Trait for types that support polynomial-like indexed access.
pub trait PolynomialAccess<P: FieldParams>: HasEndIndex {
    fn get_value(&self, idx: usize) -> Field<P>;
    fn end_index_value(&self) -> usize;
}

impl<P: FieldParams> HasEndIndex for bbrs_polynomials::polynomial::Polynomial<P> {
    fn end_index(&self) -> usize {
        self.end_index()
    }
}

impl<P: FieldParams> PolynomialAccess<P> for bbrs_polynomials::polynomial::Polynomial<P> {
    fn get_value(&self, idx: usize) -> Field<P> {
        self.get(idx)
    }

    fn end_index_value(&self) -> usize {
        self.end_index()
    }
}

// ============================================================================
// initialize_relation_separator
// ============================================================================

/// Create subrelation separators from a single alpha challenge.
///
/// alphas[0] = alpha, alphas[1] = alpha^2, ...
///
/// Port of C++ `initialize_relation_separator`.
pub fn initialize_relation_separator<P: FieldParams>(
    alpha: Field<P>,
) -> SubrelationSeparators<P> {
    let mut alphas = [Field::zero(); NUM_SUBRELATIONS - 1];
    alphas[0] = alpha;
    for i in 1..alphas.len() {
        alphas[i] = alphas[i - 1] * alpha;
    }
    alphas
}
