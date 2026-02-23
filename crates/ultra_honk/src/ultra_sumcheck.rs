//! Ultra-specific sumcheck prover and verifier.
//!
//! This implements the sumcheck protocol for the full Ultra flavor with all 9
//! relations and 28 subrelations. Uses a point-evaluation approach that reuses
//! the existing relation accumulate functions from the relations crate.

use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};
use bbrs_flavor::ultra_flavor::{
    AllEntities, AllValues, ProverPolynomials, NUM_ALL_ENTITIES,
    BATCHED_RELATION_PARTIAL_LENGTH, MAX_PARTIAL_RELATION_LENGTH, NUM_SUBRELATIONS,
};
use bbrs_polynomials::gate_separator::GateSeparatorPolynomial;
use bbrs_polynomials::polynomial::Polynomial;
use bbrs_polynomials::univariate::Univariate;
use bbrs_relations::relation_parameters::RelationParameters;
use bbrs_relations::ultra::input_elements::{InputElements, NUM_ELEMENTS};
use bbrs_relations::ultra::{
    arithmetic, delta_range, elliptic, logderiv_lookup, memory, non_native_field, permutation,
    poseidon2_external, poseidon2_internal,
};
use bbrs_transcript::codec::FieldSerializable;
use bbrs_transcript::NativeTranscript;

/// Subrelation separator challenges array.
pub type SubrelationSeparators = [Fr; NUM_SUBRELATIONS - 1];

/// The round univariate type.
pub type SumcheckRoundUnivariate = Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH>;

/// Sumcheck output containing the multivariate challenge and claimed evaluations.
pub struct UltraSumcheckOutput {
    /// The multivariate challenge vector u = (u_0, ..., u_{d-1}).
    pub challenge: Vec<Fr>,
    /// Evaluations of all polynomials at the challenge point.
    pub claimed_evaluations: AllValues<Bn254FrParams>,
    /// Whether verification passed (verifier only).
    pub verified: bool,
}

/// Create subrelation separators from a single alpha challenge.
pub fn initialize_relation_separator(alpha: Fr) -> SubrelationSeparators {
    let mut alphas = [Fr::zero(); NUM_SUBRELATIONS - 1];
    alphas[0] = alpha;
    for i in 1..alphas.len() {
        alphas[i] = alphas[i - 1] * alpha;
    }
    alphas
}

/// Wrapper for transcript serialization of Univariates.
struct UnivariateFrWrapper(SumcheckRoundUnivariate);

impl FieldSerializable for UnivariateFrWrapper {
    const NUM_FR: usize = BATCHED_RELATION_PARTIAL_LENGTH;

    fn serialize_to_frs(&self) -> Vec<Fr> {
        self.0.evaluations.to_vec()
    }

    fn deserialize_from_frs(frs: &[Fr]) -> Self {
        let mut evals = [Fr::zero(); BATCHED_RELATION_PARTIAL_LENGTH];
        evals.copy_from_slice(frs);
        Self(Univariate::new(evals))
    }
}

// ============================================================================
// AllValues <-> InputElements conversion
// ============================================================================

/// Convert AllValues (Ultra flavor, 41 fields) to InputElements (relations, 45 fields).
fn all_values_to_input_elements(av: &AllValues<Bn254FrParams>) -> InputElements<Bn254FrParams> {
    let mut ie = InputElements {
        data: [Fr::zero(); NUM_ELEMENTS],
    };
    ie.data[0] = av.q_c;
    ie.data[1] = av.q_l;
    ie.data[2] = av.q_r;
    ie.data[3] = av.q_o;
    ie.data[4] = av.q_4;
    ie.data[5] = av.q_m;
    ie.data[6] = av.q_arith;
    ie.data[7] = av.q_delta_range;
    ie.data[8] = av.q_elliptic;
    ie.data[9] = av.q_memory;
    ie.data[10] = av.q_nnf;
    ie.data[11] = av.q_lookup;
    ie.data[12] = av.q_poseidon2_external;
    ie.data[13] = av.q_poseidon2_internal;
    ie.data[14] = av.sigma_1;
    ie.data[15] = av.sigma_2;
    ie.data[16] = av.sigma_3;
    ie.data[17] = av.sigma_4;
    ie.data[18] = av.id_1;
    ie.data[19] = av.id_2;
    ie.data[20] = av.id_3;
    ie.data[21] = av.id_4;
    ie.data[22] = av.table_1;
    ie.data[23] = av.table_2;
    ie.data[24] = av.table_3;
    ie.data[25] = av.table_4;
    ie.data[26] = av.lagrange_first;
    ie.data[27] = av.lagrange_last;
    ie.data[28] = av.w_l;
    ie.data[29] = av.w_r;
    ie.data[30] = av.w_o;
    ie.data[31] = av.w_4;
    // 32: sorted_accum (unused)
    ie.data[33] = av.z_perm;
    // 34: z_lookup (unused)
    ie.data[35] = av.w_l_shift;
    ie.data[36] = av.w_r_shift;
    ie.data[37] = av.w_o_shift;
    ie.data[38] = av.w_4_shift;
    // 39: sorted_accum_shift (unused)
    ie.data[40] = av.z_perm_shift;
    // 41: z_lookup_shift (unused)
    ie.data[42] = av.lookup_read_counts;
    ie.data[43] = av.lookup_read_tags;
    ie.data[44] = av.lookup_inverses;
    ie
}

// ============================================================================
// Subrelation info
// ============================================================================

/// Partial lengths of all 28 Ultra subrelations.
const SUBRELATION_PARTIAL_LENGTHS: [usize; NUM_SUBRELATIONS] = [
    // Arithmetic: 2 subrelations
    6, 5,
    // Permutation: 2 subrelations
    6, 3,
    // LogDerivLookup: 3 subrelations
    5, 5, 3,
    // DeltaRange: 4 subrelations
    6, 6, 6, 6,
    // Elliptic: 2 subrelations
    6, 6,
    // Memory: 6 subrelations
    6, 6, 6, 6, 6, 6,
    // NonNativeField: 1 subrelation
    6,
    // Poseidon2External: 4 subrelations
    7, 7, 7, 7,
    // Poseidon2Internal: 4 subrelations
    7, 7, 7, 7,
];

/// Whether each subrelation is linearly independent.
/// LogDerivLookup subrelation index 1 (global index 4) is linearly dependent.
const SUBRELATION_LINEARLY_INDEPENDENT: [bool; NUM_SUBRELATIONS] = [
    // Arithmetic
    true, true,
    // Permutation
    true, true,
    // LogDerivLookup
    true, false, true,
    // DeltaRange
    true, true, true, true,
    // Elliptic
    true, true,
    // Memory
    true, true, true, true, true, true,
    // NonNativeField
    true,
    // Poseidon2External
    true, true, true, true,
    // Poseidon2Internal
    true, true, true, true,
];

/// Accumulate all 9 Ultra relations into a flat subrelation array.
fn accumulate_all_relations(
    evals: &mut [Fr; NUM_SUBRELATIONS],
    row: &InputElements<Bn254FrParams>,
    params: &RelationParameters<Fr>,
    scaling_factor: &Fr,
) {
    let mut offset = 0;

    // Arithmetic (2 subrelations)
    {
        let mut sub = [Fr::zero(); arithmetic::NUM_SUBRELATIONS];
        arithmetic::accumulate(&mut sub, row, params, scaling_factor);
        for i in 0..arithmetic::NUM_SUBRELATIONS {
            evals[offset + i] = evals[offset + i] + sub[i];
        }
        offset += arithmetic::NUM_SUBRELATIONS;
    }

    // Permutation (2 subrelations)
    {
        let mut sub = [Fr::zero(); permutation::NUM_SUBRELATIONS];
        permutation::accumulate(&mut sub, row, params, scaling_factor);
        for i in 0..permutation::NUM_SUBRELATIONS {
            evals[offset + i] = evals[offset + i] + sub[i];
        }
        offset += permutation::NUM_SUBRELATIONS;
    }

    // LogDerivLookup (3 subrelations)
    {
        let mut sub = [Fr::zero(); logderiv_lookup::NUM_SUBRELATIONS];
        logderiv_lookup::accumulate(&mut sub, row, params, scaling_factor);
        for i in 0..logderiv_lookup::NUM_SUBRELATIONS {
            evals[offset + i] = evals[offset + i] + sub[i];
        }
        offset += logderiv_lookup::NUM_SUBRELATIONS;
    }

    // DeltaRange (4 subrelations)
    {
        let mut sub = [Fr::zero(); delta_range::NUM_SUBRELATIONS];
        delta_range::accumulate(&mut sub, row, params, scaling_factor);
        for i in 0..delta_range::NUM_SUBRELATIONS {
            evals[offset + i] = evals[offset + i] + sub[i];
        }
        offset += delta_range::NUM_SUBRELATIONS;
    }

    // Elliptic (2 subrelations)
    {
        let mut sub = [Fr::zero(); elliptic::NUM_SUBRELATIONS];
        elliptic::accumulate(&mut sub, row, params, scaling_factor);
        for i in 0..elliptic::NUM_SUBRELATIONS {
            evals[offset + i] = evals[offset + i] + sub[i];
        }
        offset += elliptic::NUM_SUBRELATIONS;
    }

    // Memory (6 subrelations)
    {
        let mut sub = [Fr::zero(); memory::NUM_SUBRELATIONS];
        memory::accumulate(&mut sub, row, params, scaling_factor);
        for i in 0..memory::NUM_SUBRELATIONS {
            evals[offset + i] = evals[offset + i] + sub[i];
        }
        offset += memory::NUM_SUBRELATIONS;
    }

    // NonNativeField (1 subrelation)
    {
        let mut sub = [Fr::zero(); non_native_field::NUM_SUBRELATIONS];
        non_native_field::accumulate(&mut sub, row, params, scaling_factor);
        for i in 0..non_native_field::NUM_SUBRELATIONS {
            evals[offset + i] = evals[offset + i] + sub[i];
        }
        offset += non_native_field::NUM_SUBRELATIONS;
    }

    // Poseidon2External (4 subrelations)
    {
        let mut sub = [Fr::zero(); poseidon2_external::NUM_SUBRELATIONS];
        poseidon2_external::accumulate(&mut sub, row, params, scaling_factor);
        for i in 0..poseidon2_external::NUM_SUBRELATIONS {
            evals[offset + i] = evals[offset + i] + sub[i];
        }
        offset += poseidon2_external::NUM_SUBRELATIONS;
    }

    // Poseidon2Internal (4 subrelations)
    {
        let mut sub = [Fr::zero(); poseidon2_internal::NUM_SUBRELATIONS];
        poseidon2_internal::accumulate(&mut sub, row, params, scaling_factor);
        for i in 0..poseidon2_internal::NUM_SUBRELATIONS {
            evals[offset + i] = evals[offset + i] + sub[i];
        }
    }
}

// ============================================================================
// UltraSumcheckProver
// ============================================================================

/// Ultra sumcheck prover using point-evaluation approach.
pub struct UltraSumcheckProver<'a> {
    pub multivariate_n: usize,
    pub multivariate_d: usize,
    pub polynomials: &'a ProverPolynomials<Bn254FrParams>,
    pub alphas: SubrelationSeparators,
    pub gate_challenges: Vec<Fr>,
    pub relation_parameters: RelationParameters<Fr>,
}

impl<'a> UltraSumcheckProver<'a> {
    pub fn new(
        multivariate_n: usize,
        polynomials: &'a ProverPolynomials<Bn254FrParams>,
        alpha: Fr,
        gate_challenges: Vec<Fr>,
        relation_parameters: RelationParameters<Fr>,
    ) -> Self {
        let multivariate_d = multivariate_n.trailing_zeros() as usize;
        Self {
            multivariate_n,
            multivariate_d,
            polynomials,
            alphas: initialize_relation_separator(alpha),
            gate_challenges,
            relation_parameters,
        }
    }

    /// Run the sumcheck protocol.
    pub fn prove(self, transcript: &mut NativeTranscript) -> UltraSumcheckOutput {
        let mut gate_separators =
            GateSeparatorPolynomial::new(self.gate_challenges.clone(), self.multivariate_d);

        let mut multivariate_challenge = Vec::with_capacity(self.multivariate_d);

        // Create partially evaluated polynomials (initially a copy)
        let mut pep = create_partially_evaluated(self.polynomials, self.multivariate_n);

        let mut round_size = self.multivariate_n;

        for round_idx in 0..self.multivariate_d {
            // Compute the round univariate
            let round_univariate = self.compute_round_univariate(
                &pep,
                round_size,
                &gate_separators,
            );

            // Send to transcript
            transcript.send_to_verifier(
                &format!("Sumcheck:univariate_{}", round_idx),
                &UnivariateFrWrapper(round_univariate),
            );

            // Get challenge
            let round_challenge = transcript.get_challenge(
                &format!("Sumcheck:u_{}", round_idx),
            );
            multivariate_challenge.push(round_challenge);

            // Partially evaluate
            partially_evaluate_in_place(&mut pep, round_challenge, round_size);
            gate_separators.partially_evaluate(round_challenge);
            round_size >>= 1;
        }

        // Extract claimed evaluations
        let claimed_evaluations = extract_evaluations(&pep);

        // Send evaluations
        let eval_array = all_values_to_fr_array(&claimed_evaluations);
        transcript.send_to_verifier("Sumcheck:evaluations", &eval_array);

        UltraSumcheckOutput {
            challenge: multivariate_challenge,
            claimed_evaluations,
            verified: false,
        }
    }

    /// Compute the round univariate for the current round.
    fn compute_round_univariate(
        &self,
        pep: &AllEntities<Polynomial<Bn254FrParams>>,
        round_size: usize,
        gate_separators: &GateSeparatorPolynomial<Bn254FrParams>,
    ) -> SumcheckRoundUnivariate {
        // For each subrelation, accumulate contributions at each evaluation point
        let mut sub_accumulators = [[Fr::zero(); MAX_PARTIAL_RELATION_LENGTH]; NUM_SUBRELATIONS];

        let mut edge_idx = 0;
        while edge_idx < round_size {
            let scaling_factor = gate_separators.at(edge_idx);

            // For each evaluation point k in 0..MAX_PARTIAL_RELATION_LENGTH
            for k in 0..MAX_PARTIAL_RELATION_LENGTH {
                let k_fr = Fr::from(k as u64);

                // Evaluate all polynomials at point k by linear interpolation
                let row = interpolate_row(pep, edge_idx, k_fr);
                let ie = all_values_to_input_elements(&row);

                // Accumulate all relations
                let mut sub_evals = [Fr::zero(); NUM_SUBRELATIONS];
                accumulate_all_relations(
                    &mut sub_evals,
                    &ie,
                    &self.relation_parameters,
                    &scaling_factor,
                );

                for s in 0..NUM_SUBRELATIONS {
                    if k < SUBRELATION_PARTIAL_LENGTHS[s] {
                        sub_accumulators[s][k] = sub_accumulators[s][k] + sub_evals[s];
                    }
                }
            }

            edge_idx += 2;
        }

        // Batch: scale by alphas, extend to BATCHED_RELATION_PARTIAL_LENGTH, multiply by pow
        batch_subrelation_accumulators(
            &sub_accumulators,
            &self.alphas,
            gate_separators,
        )
    }
}

/// Extend evaluations from `partial_len` points to `BATCHED_RELATION_PARTIAL_LENGTH` points
/// using Newton forward differences (same algorithm as `Univariate::extend_to`).
fn extend_evaluations(
    src: &[Fr; MAX_PARTIAL_RELATION_LENGTH],
    partial_len: usize,
    dst: &mut [Fr; BATCHED_RELATION_PARTIAL_LENGTH],
) {
    // Copy known evaluations
    for i in 0..partial_len.min(BATCHED_RELATION_PARTIAL_LENGTH) {
        dst[i] = src[i];
    }

    if BATCHED_RELATION_PARTIAL_LENGTH <= partial_len {
        return;
    }

    if partial_len <= 1 {
        // Constant — already copied
        return;
    }

    if partial_len == 2 {
        // Linear extrapolation: f(x) = f(0) + x * delta
        let delta = src[1] - src[0];
        for i in 2..BATCHED_RELATION_PARTIAL_LENGTH {
            dst[i] = dst[i - 1] + delta;
        }
        return;
    }

    // General case: Newton forward differences
    let mut diffs = [Fr::zero(); MAX_PARTIAL_RELATION_LENGTH];
    for i in 0..partial_len {
        diffs[i] = src[i];
    }

    // Convert evaluations to forward-difference coefficients
    for i in 1..partial_len {
        for j in (i..partial_len).rev() {
            diffs[j] = diffs[j] - diffs[j - 1];
        }
    }

    // Shift origin iteratively and record extrapolated evaluations
    for x in 1..BATCHED_RELATION_PARTIAL_LENGTH {
        for j in 0..partial_len - 1 {
            diffs[j] = diffs[j] + diffs[j + 1];
        }
        if x >= partial_len {
            dst[x] = diffs[0];
        }
    }
}

/// Batch subrelation accumulators into a single round univariate.
fn batch_subrelation_accumulators(
    accumulators: &[[Fr; MAX_PARTIAL_RELATION_LENGTH]; NUM_SUBRELATIONS],
    alphas: &SubrelationSeparators,
    gate_separators: &GateSeparatorPolynomial<Bn254FrParams>,
) -> SumcheckRoundUnivariate {
    // Pow-factor univariate: (1 - X) + X * beta_i
    let random_polynomial = Univariate::<Bn254FrParams, 2>::new([
        Fr::one(),
        gate_separators.current_element(),
    ]);
    let extended_random_poly: Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH> =
        random_polynomial.extend_to();

    let mut result = Univariate::<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH>::zero();

    for s in 0..NUM_SUBRELATIONS {
        let partial_len = SUBRELATION_PARTIAL_LENGTHS[s];

        // Extend from partial_len evaluation points to BATCHED_RELATION_PARTIAL_LENGTH
        // using Newton forward differences (proper barycentric extension)
        let mut extended = [Fr::zero(); BATCHED_RELATION_PARTIAL_LENGTH];
        extend_evaluations(&accumulators[s], partial_len, &mut extended);
        let univariate = Univariate::<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH>::new(extended);

        // Scale by alpha (first subrelation not scaled)
        let scaled = if s == 0 {
            univariate
        } else {
            let mut u = univariate;
            for e in u.evaluations.iter_mut() {
                *e = *e * alphas[s - 1];
            }
            u
        };

        // Multiply linearly independent subrelations by pow polynomial
        if SUBRELATION_LINEARLY_INDEPENDENT[s] {
            let contrib = pointwise_mul(&scaled, &extended_random_poly);
            let contrib = pointwise_mul_scalar(&contrib, gate_separators.partial_evaluation_result);
            result = pointwise_add(&result, &contrib);
        } else {
            result = pointwise_add(&result, &scaled);
        }
    }

    result
}

/// Interpolate all polynomial values at a fractional point within an edge pair.
/// row[k] = poly[2i] + k * (poly[2i+1] - poly[2i])
fn interpolate_row(
    pep: &AllEntities<Polynomial<Bn254FrParams>>,
    edge_idx: usize,
    k: Fr,
) -> AllValues<Bn254FrParams> {
    let interp = |poly: &Polynomial<Bn254FrParams>| -> Fr {
        let v0 = poly.get(edge_idx);
        let v1 = poly.get(edge_idx + 1);
        v0 + k * (v1 - v0)
    };

    AllEntities {
        q_m: interp(&pep.q_m),
        q_c: interp(&pep.q_c),
        q_l: interp(&pep.q_l),
        q_r: interp(&pep.q_r),
        q_o: interp(&pep.q_o),
        q_4: interp(&pep.q_4),
        q_lookup: interp(&pep.q_lookup),
        q_arith: interp(&pep.q_arith),
        q_delta_range: interp(&pep.q_delta_range),
        q_elliptic: interp(&pep.q_elliptic),
        q_memory: interp(&pep.q_memory),
        q_nnf: interp(&pep.q_nnf),
        q_poseidon2_external: interp(&pep.q_poseidon2_external),
        q_poseidon2_internal: interp(&pep.q_poseidon2_internal),
        sigma_1: interp(&pep.sigma_1),
        sigma_2: interp(&pep.sigma_2),
        sigma_3: interp(&pep.sigma_3),
        sigma_4: interp(&pep.sigma_4),
        id_1: interp(&pep.id_1),
        id_2: interp(&pep.id_2),
        id_3: interp(&pep.id_3),
        id_4: interp(&pep.id_4),
        table_1: interp(&pep.table_1),
        table_2: interp(&pep.table_2),
        table_3: interp(&pep.table_3),
        table_4: interp(&pep.table_4),
        lagrange_first: interp(&pep.lagrange_first),
        lagrange_last: interp(&pep.lagrange_last),
        w_l: interp(&pep.w_l),
        w_r: interp(&pep.w_r),
        w_o: interp(&pep.w_o),
        w_4: interp(&pep.w_4),
        z_perm: interp(&pep.z_perm),
        lookup_inverses: interp(&pep.lookup_inverses),
        lookup_read_counts: interp(&pep.lookup_read_counts),
        lookup_read_tags: interp(&pep.lookup_read_tags),
        w_l_shift: interp(&pep.w_l_shift),
        w_r_shift: interp(&pep.w_r_shift),
        w_o_shift: interp(&pep.w_o_shift),
        w_4_shift: interp(&pep.w_4_shift),
        z_perm_shift: interp(&pep.z_perm_shift),
    }
}

/// Create partially evaluated polynomials from prover polynomials.
fn create_partially_evaluated(
    src: &ProverPolynomials<Bn254FrParams>,
    size: usize,
) -> AllEntities<Polynomial<Bn254FrParams>> {
    let copy = |p: &Polynomial<Bn254FrParams>| -> Polynomial<Bn254FrParams> {
        let mut new_p = Polynomial::new(size, size, 0);
        for i in 0..p.end_index().min(size) {
            *new_p.at_mut(i) = p.get(i);
        }
        new_p
    };

    AllEntities {
        q_m: copy(&src.q_m), q_c: copy(&src.q_c), q_l: copy(&src.q_l),
        q_r: copy(&src.q_r), q_o: copy(&src.q_o), q_4: copy(&src.q_4),
        q_lookup: copy(&src.q_lookup), q_arith: copy(&src.q_arith),
        q_delta_range: copy(&src.q_delta_range), q_elliptic: copy(&src.q_elliptic),
        q_memory: copy(&src.q_memory), q_nnf: copy(&src.q_nnf),
        q_poseidon2_external: copy(&src.q_poseidon2_external),
        q_poseidon2_internal: copy(&src.q_poseidon2_internal),
        sigma_1: copy(&src.sigma_1), sigma_2: copy(&src.sigma_2),
        sigma_3: copy(&src.sigma_3), sigma_4: copy(&src.sigma_4),
        id_1: copy(&src.id_1), id_2: copy(&src.id_2),
        id_3: copy(&src.id_3), id_4: copy(&src.id_4),
        table_1: copy(&src.table_1), table_2: copy(&src.table_2),
        table_3: copy(&src.table_3), table_4: copy(&src.table_4),
        lagrange_first: copy(&src.lagrange_first), lagrange_last: copy(&src.lagrange_last),
        w_l: copy(&src.w_l), w_r: copy(&src.w_r),
        w_o: copy(&src.w_o), w_4: copy(&src.w_4),
        z_perm: copy(&src.z_perm),
        lookup_inverses: copy(&src.lookup_inverses),
        lookup_read_counts: copy(&src.lookup_read_counts),
        lookup_read_tags: copy(&src.lookup_read_tags),
        w_l_shift: copy(&src.w_l_shift), w_r_shift: copy(&src.w_r_shift),
        w_o_shift: copy(&src.w_o_shift), w_4_shift: copy(&src.w_4_shift),
        z_perm_shift: copy(&src.z_perm_shift),
    }
}

/// Partially evaluate polynomials in place: pep[i/2] = pep[2i] + challenge * (pep[2i+1] - pep[2i]).
fn partially_evaluate_in_place(
    pep: &mut AllEntities<Polynomial<Bn254FrParams>>,
    challenge: Fr,
    round_size: usize,
) {
    let fold = |poly: &mut Polynomial<Bn254FrParams>| {
        let limit = poly.end_index().min(round_size);
        let mut i = 0;
        while i < limit {
            let v0 = poly.get(i);
            let v1 = poly.get(i + 1);
            *poly.at_mut(i >> 1) = v0 + challenge * (v1 - v0);
            i += 2;
        }
        let new_end = (limit / 2) + (limit % 2);
        poly.shrink_end_index(new_end);
    };

    for poly in pep.get_all_mut() {
        fold(poly);
    }
}

/// Extract evaluations from partially evaluated polynomials (all folded to a single value).
fn extract_evaluations(pep: &AllEntities<Polynomial<Bn254FrParams>>) -> AllValues<Bn254FrParams> {
    AllEntities {
        q_m: pep.q_m.get(0), q_c: pep.q_c.get(0), q_l: pep.q_l.get(0),
        q_r: pep.q_r.get(0), q_o: pep.q_o.get(0), q_4: pep.q_4.get(0),
        q_lookup: pep.q_lookup.get(0), q_arith: pep.q_arith.get(0),
        q_delta_range: pep.q_delta_range.get(0), q_elliptic: pep.q_elliptic.get(0),
        q_memory: pep.q_memory.get(0), q_nnf: pep.q_nnf.get(0),
        q_poseidon2_external: pep.q_poseidon2_external.get(0),
        q_poseidon2_internal: pep.q_poseidon2_internal.get(0),
        sigma_1: pep.sigma_1.get(0), sigma_2: pep.sigma_2.get(0),
        sigma_3: pep.sigma_3.get(0), sigma_4: pep.sigma_4.get(0),
        id_1: pep.id_1.get(0), id_2: pep.id_2.get(0),
        id_3: pep.id_3.get(0), id_4: pep.id_4.get(0),
        table_1: pep.table_1.get(0), table_2: pep.table_2.get(0),
        table_3: pep.table_3.get(0), table_4: pep.table_4.get(0),
        lagrange_first: pep.lagrange_first.get(0), lagrange_last: pep.lagrange_last.get(0),
        w_l: pep.w_l.get(0), w_r: pep.w_r.get(0),
        w_o: pep.w_o.get(0), w_4: pep.w_4.get(0),
        z_perm: pep.z_perm.get(0),
        lookup_inverses: pep.lookup_inverses.get(0),
        lookup_read_counts: pep.lookup_read_counts.get(0),
        lookup_read_tags: pep.lookup_read_tags.get(0),
        w_l_shift: pep.w_l_shift.get(0), w_r_shift: pep.w_r_shift.get(0),
        w_o_shift: pep.w_o_shift.get(0), w_4_shift: pep.w_4_shift.get(0),
        z_perm_shift: pep.z_perm_shift.get(0),
    }
}

/// Convert AllValues to a flat Fr array for transcript serialization.
fn all_values_to_fr_array(av: &AllValues<Bn254FrParams>) -> [Fr; NUM_ALL_ENTITIES] {
    let refs = av.get_all();
    let mut arr = [Fr::zero(); NUM_ALL_ENTITIES];
    for (i, val) in refs.iter().enumerate() {
        arr[i] = **val;
    }
    arr
}

/// Convert a flat Fr array from transcript back to AllValues.
fn fr_array_to_all_values(arr: &[Fr; NUM_ALL_ENTITIES]) -> AllValues<Bn254FrParams> {
    AllEntities {
        q_m: arr[0], q_c: arr[1], q_l: arr[2], q_r: arr[3], q_o: arr[4], q_4: arr[5],
        q_lookup: arr[6], q_arith: arr[7], q_delta_range: arr[8], q_elliptic: arr[9],
        q_memory: arr[10], q_nnf: arr[11], q_poseidon2_external: arr[12], q_poseidon2_internal: arr[13],
        sigma_1: arr[14], sigma_2: arr[15], sigma_3: arr[16], sigma_4: arr[17],
        id_1: arr[18], id_2: arr[19], id_3: arr[20], id_4: arr[21],
        table_1: arr[22], table_2: arr[23], table_3: arr[24], table_4: arr[25],
        lagrange_first: arr[26], lagrange_last: arr[27],
        w_l: arr[28], w_r: arr[29], w_o: arr[30], w_4: arr[31],
        z_perm: arr[32], lookup_inverses: arr[33], lookup_read_counts: arr[34], lookup_read_tags: arr[35],
        w_l_shift: arr[36], w_r_shift: arr[37], w_o_shift: arr[38], w_4_shift: arr[39], z_perm_shift: arr[40],
    }
}

// ============================================================================
// UltraSumcheckVerifier
// ============================================================================

/// Ultra sumcheck verifier.
pub struct UltraSumcheckVerifier {
    pub alphas: SubrelationSeparators,
    pub multivariate_d: usize,
}

impl UltraSumcheckVerifier {
    pub fn new(alpha: Fr, multivariate_d: usize) -> Self {
        Self {
            alphas: initialize_relation_separator(alpha),
            multivariate_d,
        }
    }

    /// Run the sumcheck verification protocol.
    pub fn verify(
        self,
        transcript: &mut NativeTranscript,
        relation_parameters: &RelationParameters<Fr>,
        gate_challenges: &[Fr],
    ) -> UltraSumcheckOutput {
        let mut gate_separators = GateSeparatorPolynomial::new_verifier(gate_challenges.to_vec());
        let mut target_total_sum = Fr::zero();
        let mut multivariate_challenge = Vec::with_capacity(self.multivariate_d);
        let mut verified = true;

        // Process each round
        for round_idx in 0..self.multivariate_d {
            let round_univariate: UnivariateFrWrapper = transcript
                .receive_from_prover(&format!("Sumcheck:univariate_{}", round_idx));
            let round_univariate = round_univariate.0;

            // Check sum: S(0) + S(1) should equal target from previous round
            // Round 0: target is 0 (total relation sum for a valid circuit)
            let sum = round_univariate.value_at(0) + round_univariate.value_at(1);
            let round_failed = target_total_sum != sum;
            #[cfg(test)]
            if round_failed {
                eprintln!("[Sumcheck Verifier] Round {} FAILED: target={:?}, sum={:?}", round_idx, target_total_sum, sum);
            }
            verified = verified && !round_failed;

            // Get challenge
            let round_challenge = transcript.get_challenge(
                &format!("Sumcheck:u_{}", round_idx),
            );
            multivariate_challenge.push(round_challenge);

            // Compute next target sum: S^i(u_i)
            target_total_sum = round_univariate.evaluate(round_challenge);
            gate_separators.partially_evaluate(round_challenge);
        }

        // Receive claimed evaluations
        let transcript_evaluations: [Fr; NUM_ALL_ENTITIES] =
            transcript.receive_from_prover("Sumcheck:evaluations");
        let purported_evaluations = fr_array_to_all_values(&transcript_evaluations);

        // Compute full relation purported value
        let ie = all_values_to_input_elements(&purported_evaluations);
        let mut sub_evals = [Fr::zero(); NUM_SUBRELATIONS];
        let scaling = gate_separators.partial_evaluation_result;
        accumulate_all_relations(&mut sub_evals, &ie, relation_parameters, &scaling);

        // Scale and batch
        let mut full_value = Fr::zero();
        for s in 0..NUM_SUBRELATIONS {
            let scaled = if s == 0 {
                sub_evals[s]
            } else {
                sub_evals[s] * self.alphas[s - 1]
            };

            if SUBRELATION_LINEARLY_INDEPENDENT[s] {
                full_value = full_value + scaled;
            } else {
                full_value = full_value + scaled;
            }
        }

        // Final check
        #[cfg(test)]
        {
            eprintln!("[Sumcheck Verifier] Final relation value: {:?}", full_value);
            eprintln!("[Sumcheck Verifier] Target total sum:     {:?}", target_total_sum);
            eprintln!("[Sumcheck Verifier] Match: {}", full_value == target_total_sum);
            eprintln!("[Sumcheck Verifier] Rounds verified: {}", verified);
        }
        verified = verified && (full_value == target_total_sum);

        UltraSumcheckOutput {
            challenge: multivariate_challenge,
            claimed_evaluations: purported_evaluations,
            verified,
        }
    }
}

// ============================================================================
// Univariate helpers
// ============================================================================

fn pointwise_mul(
    a: &Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH>,
    b: &Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH>,
) -> Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH> {
    let mut result = Univariate::zero();
    for i in 0..BATCHED_RELATION_PARTIAL_LENGTH {
        result.evaluations[i] = a.evaluations[i] * b.evaluations[i];
    }
    result
}

fn pointwise_mul_scalar(
    a: &Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH>,
    scalar: Fr,
) -> Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH> {
    let mut result = Univariate::zero();
    for i in 0..BATCHED_RELATION_PARTIAL_LENGTH {
        result.evaluations[i] = a.evaluations[i] * scalar;
    }
    result
}

fn pointwise_add(
    a: &Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH>,
    b: &Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH>,
) -> Univariate<Bn254FrParams, BATCHED_RELATION_PARTIAL_LENGTH> {
    let mut result = Univariate::zero();
    for i in 0..BATCHED_RELATION_PARTIAL_LENGTH {
        result.evaluations[i] = a.evaluations[i] + b.evaluations[i];
    }
    result
}
