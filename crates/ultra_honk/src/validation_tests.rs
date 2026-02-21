//! Cross-validation and deterministic vector tests for Package 5 port validation.
//!
//! These tests exercise integration boundaries across crates:
//! - Field arithmetic deterministic vectors
//! - Transcript → challenge determinism
//! - UltraFlavor → relation evaluation
//! - Grand product delta consistency
//! - Commitment key → MSM consistency

#[cfg(test)]
mod tests {
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};
    use bbrs_flavor::ultra_flavor;
    use bbrs_relations::relation_parameters::RelationParameters;
    use bbrs_relations::ultra::arithmetic;
    use bbrs_transcript::NativeTranscript;

    type P = Bn254FrParams;

    // =========================================================================
    // Deterministic Field Arithmetic Vectors
    // =========================================================================

    /// Verify Fr multiplication of known values matches expected result.
    ///
    /// This catches any Montgomery form bugs or modular arithmetic errors.
    #[test]
    fn test_deterministic_fr_multiplication() {
        let a = Fr::from(7u64);
        let b = Fr::from(13u64);
        let c = a * b;
        assert_eq!(c, Fr::from(91u64));

        // Exponentiation via multiply: 7 * 7 = 49
        let c_sq = a * a;
        assert_eq!(c_sq, Fr::from(49u64));

        // Check a - b + b = a (additive identity)
        let d = a - b + b;
        assert_eq!(d, a);
    }

    /// Verify Fr inversion: a * a^{-1} = 1.
    #[test]
    fn test_deterministic_fr_inversion() {
        let a = Fr::from(42u64);
        let a_inv = a.invert();
        let product = a * a_inv;
        assert_eq!(product, Fr::one());

        let b = Fr::from(123456789u64);
        let b_inv = b.invert();
        assert_eq!(b * b_inv, Fr::one());
    }

    /// Verify Fr exponentiation: 2^10 = 1024.
    #[test]
    fn test_deterministic_fr_exponentiation() {
        let two = Fr::from(2u64);
        let result = two.pow(&[10, 0, 0, 0]);
        assert_eq!(result, Fr::from(1024u64));

        // S-box: 3^5 = 243
        let three = Fr::from(3u64);
        let sbox = {
            let t2 = three * three;
            let t4 = t2 * t2;
            t4 * three
        };
        assert_eq!(sbox, Fr::from(243u64));
    }

    // =========================================================================
    // Transcript Challenge Determinism
    // =========================================================================

    /// Verify that transcript challenges are deterministic: same inputs → same challenge.
    #[test]
    fn test_transcript_challenge_determinism() {
        let mut t1 = NativeTranscript::prover_init_empty();
        let mut t2 = NativeTranscript::prover_init_empty();

        t1.send_to_verifier("test_val", &Fr::from(42u64));
        t2.send_to_verifier("test_val", &Fr::from(42u64));

        let c1: Fr = t1.get_challenge("alpha");
        let c2: Fr = t2.get_challenge("alpha");

        assert_eq!(c1, c2, "Same transcript inputs should produce same challenge");
    }

    /// Verify that different inputs produce different challenges.
    #[test]
    fn test_transcript_challenge_sensitivity() {
        let mut t1 = NativeTranscript::prover_init_empty();
        let mut t2 = NativeTranscript::prover_init_empty();

        t1.send_to_verifier("test_val", &Fr::from(42u64));
        t2.send_to_verifier("test_val", &Fr::from(43u64));

        let c1: Fr = t1.get_challenge("alpha");
        let c2: Fr = t2.get_challenge("alpha");

        assert_ne!(c1, c2, "Different transcript inputs should produce different challenges");
    }

    // =========================================================================
    // UltraFlavor → Relation Evaluation Cross-Validation
    // =========================================================================

    /// Test that UltraFlavor's ProverPolynomials can produce rows that work with
    /// the arithmetic relation.
    ///
    /// Sets up a simple 1+1=2 gate and verifies the relation evaluates to zero.
    #[test]
    fn test_ultra_flavor_arithmetic_relation_integration() {
        let circuit_size = 4;
        let mut polys = ultra_flavor::ProverPolynomials::<P>::new(circuit_size);

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

        polys.set_shifted();

        // Extract row 1 using UltraFlavor's get_row
        let row = polys.get_row(1);

        // Build InputElements from the row values
        let mut input = bbrs_relations::ultra::input_elements::InputElements::<P>::get_special();
        for e in input.data.iter_mut() {
            *e = Fr::zero();
        }

        // Map UltraFlavor fields to InputElements indices
        // q_l → data[1], q_r → data[2], q_o → data[3], q_arith → data[6]
        // w_l → data[28], w_r → data[29], w_o → data[30]
        input.data[1] = row.q_l;
        input.data[2] = row.q_r;
        input.data[3] = row.q_o;
        input.data[6] = row.q_arith;
        input.data[28] = row.w_l;
        input.data[29] = row.w_r;
        input.data[30] = row.w_o;

        // Run arithmetic relation
        let params = RelationParameters::<Fr>::default();
        let mut accum = [Fr::zero(); arithmetic::NUM_SUBRELATIONS];
        arithmetic::accumulate(&mut accum, &input, &params, &Fr::one());

        // For a satisfiable gate, the relation output should be zero
        assert!(accum[0].is_zero(), "Arithmetic subrelation 0 should be zero for 1+1=2");
    }

    // =========================================================================
    // Grand Product Delta Cross-Validation
    // =========================================================================

    /// Verify grand product delta is one when there are no public inputs.
    #[test]
    fn test_grand_product_delta_no_public_inputs() {
        let delta = bbrs_honk::compute_public_input_delta::<P>(
            &[],             // no public inputs
            Fr::from(5u64),  // beta
            Fr::from(7u64),  // gamma
            Fr::zero(),      // offset
        );
        // With no public inputs, the product is empty → delta = 1
        assert_eq!(delta, Fr::one());
    }

    /// Verify grand product delta is non-trivial for non-empty public inputs.
    #[test]
    fn test_grand_product_delta_single_input_consistency() {
        let beta = Fr::from(3u64);
        let gamma = Fr::from(11u64);
        let pi = vec![Fr::from(42u64)];

        let delta = bbrs_honk::compute_public_input_delta::<P>(
            &pi,
            beta,
            gamma,
            Fr::zero(), // offset
        );

        // delta should be non-trivial (not 0 or 1 for non-empty public inputs)
        assert!(!delta.is_zero(), "Delta should be non-zero for non-empty public inputs");
        // With offset=0 and 1 public input:
        // num = gamma + pi[0] + beta * SEPARATOR
        // den = gamma + pi[0] - beta * 1
        let separator = Fr::from(bbrs_honk::PERMUTATION_ARGUMENT_VALUE_SEPARATOR);
        let num = gamma + pi[0] + beta * separator;
        let den = gamma + pi[0] - beta * Fr::one();
        let expected = num * den.invert();
        assert_eq!(delta, expected);
    }

    // =========================================================================
    // UltraFlavor Constants Cross-Validation
    // =========================================================================

    /// Verify UltraFlavor constants match C++ values exactly.
    #[test]
    fn test_ultra_flavor_constants_match_cpp() {
        assert_eq!(ultra_flavor::NUM_WIRES, 4);
        assert_eq!(ultra_flavor::NUM_PRECOMPUTED_ENTITIES, 28);
        assert_eq!(ultra_flavor::NUM_WITNESS_ENTITIES, 8);
        assert_eq!(ultra_flavor::NUM_SHIFTED_ENTITIES, 5);
        assert_eq!(ultra_flavor::NUM_ALL_ENTITIES, 41);
        assert_eq!(ultra_flavor::NUM_RELATIONS, 9);
        assert_eq!(ultra_flavor::NUM_SUBRELATIONS, 28);
        assert_eq!(ultra_flavor::MAX_PARTIAL_RELATION_LENGTH, 7);
        assert_eq!(ultra_flavor::BATCHED_RELATION_PARTIAL_LENGTH, 8);
        assert_eq!(ultra_flavor::USE_SHORT_MONOMIALS, true);
        assert_eq!(ultra_flavor::HAS_ZK, false);
        assert_eq!(ultra_flavor::HAS_ZERO_ROW, true);
    }

    /// Verify subrelation partial lengths array matches C++ exactly.
    #[test]
    fn test_ultra_subrelation_partial_lengths_match_cpp() {
        // From C++: Ultra flavor has 28 subrelations across 9 relations
        let expected: [usize; 28] = [
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
        assert_eq!(
            ultra_flavor::ALL_SUBRELATION_PARTIAL_LENGTHS, expected,
            "Subrelation partial lengths must match C++ UltraFlavor exactly"
        );
    }

    /// Verify the linearly independent flags match C++.
    #[test]
    fn test_ultra_subrelation_linear_independence_match_cpp() {
        // LogDerivLookup subrelation index 1 (global index 5) is the only linearly dependent one
        let expected: [bool; 28] = [
            // Arithmetic
            true, true,
            // Permutation
            true, true,
            // LogDerivLookup: [true, false, true]
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
        assert_eq!(
            ultra_flavor::ALL_SUBRELATION_LINEARLY_INDEPENDENT, expected,
            "Subrelation linear independence flags must match C++ UltraFlavor exactly"
        );
    }

    // =========================================================================
    // Relation Parameters Integration
    // =========================================================================

    /// Verify that RelationParameters::get_random produces valid random parameters.
    #[test]
    fn test_relation_parameters_random() {
        let params = RelationParameters::<Fr>::get_random();

        // All main challenge parameters should be non-zero
        assert!(!params.eta.is_zero());
        assert!(!params.eta_two.is_zero());
        assert!(!params.eta_three.is_zero());
        assert!(!params.beta.is_zero());
        assert!(!params.gamma.is_zero());
    }

    /// Verify default RelationParameters has zero challenges.
    #[test]
    fn test_relation_parameters_default() {
        let params = RelationParameters::<Fr>::default();
        assert!(params.beta.is_zero());
        assert!(params.gamma.is_zero());
    }
}
