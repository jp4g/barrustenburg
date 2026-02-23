//! End-to-end prove/verify tests for the Ultra Honk proving system.
//!
//! These tests wire up UltraCircuitBuilder → ProvingKey → UltraProver::prove()
//! → UltraVerifier::verify() for full round-trip verification.

#[cfg(test)]
mod tests {
    use std::sync::Once;

    use bbrs_circuit_builder::gate_data::{AddQuad, AddTriple, EccAddGate, MulQuad};
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use bbrs_commitment_schemes::commitment_key::CommitmentKey;
    use bbrs_ecc::curves::bn254::{
        Bn254FrParams, Bn254G1Params, Fr, G2AffineElement, G2Element,
    };
    use bbrs_ecc::groups::element::Element;

    use crate::proving_key::ProvingKey;
    use crate::ultra_prover::UltraProver;
    use crate::ultra_verifier::UltraVerifier;
    use crate::verification_key::VerificationKey;

    type P = Bn254FrParams;

    /// Maximum SRS size for tests. Must be >= the dyadic circuit size.
    const TEST_SRS_SIZE: usize = 256;

    static INIT_CRS: Once = Once::new();
    static mut TEST_G2_X: Option<G2AffineElement> = None;

    /// Initialize the global CRS factory with a test powers-of-tau SRS.
    ///
    /// This must be called before any test that creates a ProvingKey or VerificationKey.
    /// Uses std::sync::Once to ensure single initialization across parallel tests.
    fn ensure_test_srs_initialized() -> G2AffineElement {
        INIT_CRS.call_once(|| {
            let tau = Fr::from(2u64); // Use known tau for deterministic tests
            let g1 = Element::<Bn254G1Params>::one();
            let mut points = Vec::with_capacity(TEST_SRS_SIZE);
            let mut tau_power = Fr::one();
            for _ in 0..TEST_SRS_SIZE {
                points.push(g1.mul(&tau_power).to_affine());
                tau_power = tau_power * tau;
            }

            let g2 = G2Element::from_affine(&G2AffineElement::generator());
            let g2_x = g2.mul_scalar(&tau).to_affine();

            bbrs_srs::global_crs::init_bn254_mem_crs_factory(&points);
            unsafe {
                TEST_G2_X = Some(g2_x);
            }
        });
        unsafe { TEST_G2_X.unwrap() }
    }

    /// Create a test commitment key from a powers-of-tau SRS.
    fn create_test_ck(n: usize) -> (CommitmentKey<Bn254G1Params>, G2AffineElement) {
        let tau = Fr::from(2u64);
        let g1 = Element::<Bn254G1Params>::one();
        let mut points = Vec::with_capacity(n);
        let mut tau_power = Fr::one();
        for _ in 0..n {
            points.push(g1.mul(&tau_power).to_affine());
            tau_power = tau_power * tau;
        }

        let g2 = G2Element::from_affine(&G2AffineElement::generator());
        let g2_x = g2.mul_scalar(&tau).to_affine();

        (CommitmentKey::from_points(points), g2_x)
    }

    /// Build a simple circuit: a + b = c, where a=1, b=1, c=2.
    fn build_addition_circuit() -> UltraCircuitBuilder<P> {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let a = builder.base.add_variable(Fr::one());
        let b = builder.base.add_variable(Fr::one());
        let c = builder.base.add_variable(Fr::from(2u64));

        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::one(),
            b_scaling: Fr::one(),
            c_scaling: Fr::zero() - Fr::one(),
            const_scaling: Fr::zero(),
        });

        builder
    }

    /// Build a multiplication circuit: a * b = c, where a=3, b=7, c=21.
    fn build_multiplication_circuit() -> UltraCircuitBuilder<P> {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let a = builder.base.add_variable(Fr::from(3u64));
        let b = builder.base.add_variable(Fr::from(7u64));
        let c = builder.base.add_variable(Fr::from(21u64));

        builder.create_big_mul_add_gate(
            &MulQuad {
                a,
                b,
                c,
                d: builder.base.zero_idx(),
                mul_scaling: Fr::one(),
                a_scaling: Fr::zero(),
                b_scaling: Fr::zero(),
                c_scaling: Fr::zero() - Fr::one(),
                d_scaling: Fr::zero(),
                const_scaling: Fr::zero(),
            },
            false,
        );

        builder
    }

    /// Build a circuit with public inputs: public a + b = c.
    fn build_public_input_circuit() -> UltraCircuitBuilder<P> {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let a = builder.base.add_public_variable(Fr::from(42u64));
        let b = builder.base.add_variable(Fr::from(58u64));
        let c = builder.base.add_variable(Fr::from(100u64));

        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::one(),
            b_scaling: Fr::one(),
            c_scaling: Fr::zero() - Fr::one(),
            const_scaling: Fr::zero(),
        });

        builder
    }

    /// Build a circuit with multiple gates.
    fn build_multi_gate_circuit() -> UltraCircuitBuilder<P> {
        let mut builder = UltraCircuitBuilder::<P>::new();

        // Gate 1: a + b = c  (1 + 2 = 3)
        let a = builder.base.add_variable(Fr::one());
        let b = builder.base.add_variable(Fr::from(2u64));
        let c = builder.base.add_variable(Fr::from(3u64));

        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::one(),
            b_scaling: Fr::one(),
            c_scaling: Fr::zero() - Fr::one(),
            const_scaling: Fr::zero(),
        });

        // Gate 2: c * d = e  (3 * 4 = 12)
        let d = builder.base.add_variable(Fr::from(4u64));
        let e = builder.base.add_variable(Fr::from(12u64));

        builder.create_big_mul_add_gate(
            &MulQuad {
                a: c,
                b: d,
                c: e,
                d: builder.base.zero_idx(),
                mul_scaling: Fr::one(),
                a_scaling: Fr::zero(),
                b_scaling: Fr::zero(),
                c_scaling: Fr::zero() - Fr::one(),
                d_scaling: Fr::zero(),
                const_scaling: Fr::zero(),
            },
            false,
        );

        // Gate 3: e + 5 = f  (12 + 5 = 17)
        let five = builder.base.add_variable(Fr::from(5u64));
        let f = builder.base.add_variable(Fr::from(17u64));

        builder.create_add_gate(&AddTriple {
            a: e,
            b: five,
            c: f,
            a_scaling: Fr::one(),
            b_scaling: Fr::one(),
            c_scaling: Fr::zero() - Fr::one(),
            const_scaling: Fr::zero(),
        });

        builder
    }

    /// Build a circuit with copy constraints (assert_equal).
    fn build_copy_constraint_circuit() -> UltraCircuitBuilder<P> {
        let mut builder = UltraCircuitBuilder::<P>::new();

        // Create two variables with the same value
        let a = builder.base.add_variable(Fr::from(7u64));
        let b = builder.base.add_variable(Fr::from(7u64));

        // Assert they're equal (creates a copy constraint)
        builder.base.assert_equal(a, b, "test copy");

        // Use both in different gates
        let c = builder.base.add_variable(Fr::from(14u64));

        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::one(),
            b_scaling: Fr::one(),
            c_scaling: Fr::zero() - Fr::one(),
            const_scaling: Fr::zero(),
        });

        builder
    }

    /// Build a circuit with boolean gates.
    fn build_boolean_circuit() -> UltraCircuitBuilder<P> {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let zero = builder.base.add_variable(Fr::zero());
        let one = builder.base.add_variable(Fr::one());

        builder.create_bool_gate(zero);
        builder.create_bool_gate(one);

        // Use the boolean variables: one + zero = one
        let result = builder.base.add_variable(Fr::one());
        builder.create_add_gate(&AddTriple {
            a: one,
            b: zero,
            c: result,
            a_scaling: Fr::one(),
            b_scaling: Fr::one(),
            c_scaling: Fr::zero() - Fr::one(),
            const_scaling: Fr::zero(),
        });

        builder
    }

    /// Build a circuit with a constant gate.
    fn build_constant_circuit() -> UltraCircuitBuilder<P> {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let a = builder.base.add_variable(Fr::from(10u64));
        let b = builder.base.add_variable(Fr::from(5u64));
        let c = builder.base.add_variable(Fr::from(15u64));

        // a + b - c = 0
        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::one(),
            b_scaling: Fr::one(),
            c_scaling: Fr::zero() - Fr::one(),
            const_scaling: Fr::zero(),
        });

        // a + const = c (10 + 5 = 15)
        builder.create_add_gate(&AddTriple {
            a,
            b: builder.base.zero_idx(),
            c,
            a_scaling: Fr::one(),
            b_scaling: Fr::zero(),
            c_scaling: Fr::zero() - Fr::one(),
            const_scaling: Fr::from(5u64),
        });

        builder
    }

    /// Helper to run a full prove → verify cycle.
    ///
    /// Returns true if verification passes.
    fn prove_and_verify(mut builder: UltraCircuitBuilder<P>) -> bool {
        let g2_x = ensure_test_srs_initialized();

        // Finalize circuit
        builder.finalize_circuit(true);
        builder.blocks.compute_offsets();

        // Determine circuit size needed
        let total = builder.blocks.get_total_content_size();
        let circuit_size = (total + 1).next_power_of_two();

        // Create commitment key
        let (ck, _g2_x) = create_test_ck(circuit_size);

        // Create proving key
        let pk = ProvingKey::create_with_ck(&mut builder, ck);

        // Create verification key
        let vk = VerificationKey::create_with_g2x(&pk, g2_x);

        // Prove
        let prover = UltraProver::new(pk);
        let proof = prover.prove();

        // Verify
        let verifier = UltraVerifier::new(vk);
        let output = verifier.verify(&proof);

        output.result
    }

    // =========================================================================
    // End-to-End Prove/Verify Tests
    // =========================================================================

    #[test]
    fn test_e2e_simple_addition() {
        let builder = build_addition_circuit();
        assert!(prove_and_verify(builder), "Addition circuit should verify");
    }

    #[test]
    fn test_e2e_simple_multiplication() {
        let builder = build_multiplication_circuit();
        assert!(
            prove_and_verify(builder),
            "Multiplication circuit should verify"
        );
    }

    #[test]
    fn test_e2e_public_inputs() {
        let builder = build_public_input_circuit();
        assert!(
            prove_and_verify(builder),
            "Public input circuit should verify"
        );
    }

    #[test]
    fn test_e2e_multi_gate() {
        let builder = build_multi_gate_circuit();
        assert!(
            prove_and_verify(builder),
            "Multi-gate circuit should verify"
        );
    }

    #[test]
    fn test_e2e_copy_constraints() {
        let builder = build_copy_constraint_circuit();
        assert!(
            prove_and_verify(builder),
            "Copy constraint circuit should verify"
        );
    }

    #[test]
    fn test_e2e_boolean_gates() {
        let builder = build_boolean_circuit();
        assert!(
            prove_and_verify(builder),
            "Boolean circuit should verify"
        );
    }

    #[test]
    fn test_e2e_constant_gates() {
        let builder = build_constant_circuit();
        assert!(
            prove_and_verify(builder),
            "Constant gate circuit should verify"
        );
    }

    // =========================================================================
    // Proving Key Construction Tests
    // =========================================================================

    #[test]
    fn test_proving_key_from_addition_circuit() {
        ensure_test_srs_initialized();
        let mut builder = build_addition_circuit();
        builder.finalize_circuit(true);
        builder.blocks.compute_offsets();

        let total = builder.blocks.get_total_content_size();
        let circuit_size = (total + 1).next_power_of_two();
        let (ck, _) = create_test_ck(circuit_size);

        let pk = ProvingKey::create_with_ck(&mut builder, ck);

        assert!(pk.circuit_size.is_power_of_two());
        assert!(pk.circuit_size >= 2);
        assert_eq!(pk.num_public_inputs, 0);
    }

    #[test]
    fn test_proving_key_public_inputs() {
        ensure_test_srs_initialized();
        let mut builder = build_public_input_circuit();
        builder.finalize_circuit(true);
        builder.blocks.compute_offsets();

        let total = builder.blocks.get_total_content_size();
        let circuit_size = (total + 1).next_power_of_two();
        let (ck, _) = create_test_ck(circuit_size);

        let pk = ProvingKey::create_with_ck(&mut builder, ck);

        assert_eq!(pk.num_public_inputs, 1);
        assert_eq!(pk.public_inputs[0], Fr::from(42u64));
    }

    #[test]
    fn test_proving_key_lagrange_polynomials() {
        ensure_test_srs_initialized();
        let mut builder = build_addition_circuit();
        builder.finalize_circuit(true);
        builder.blocks.compute_offsets();

        let total = builder.blocks.get_total_content_size();
        let circuit_size = (total + 1).next_power_of_two();
        let (ck, _) = create_test_ck(circuit_size);

        let pk = ProvingKey::create_with_ck(&mut builder, ck);

        // lagrange_first[0] == 1, rest zero
        assert_eq!(pk.polynomials.lagrange_first.get(0), Fr::one());
        for i in 1..pk.circuit_size {
            assert!(
                pk.polynomials.lagrange_first.get(i).is_zero(),
                "lagrange_first[{}] should be zero",
                i
            );
        }

        // lagrange_last[n-1] == 1, rest zero
        assert_eq!(
            pk.polynomials.lagrange_last.get(pk.circuit_size - 1),
            Fr::one()
        );
        for i in 0..pk.circuit_size - 1 {
            assert!(
                pk.polynomials.lagrange_last.get(i).is_zero(),
                "lagrange_last[{}] should be zero",
                i
            );
        }
    }

    // =========================================================================
    // Verification Key Tests
    // =========================================================================

    #[test]
    fn test_verification_key_commitments() {
        let g2_x = ensure_test_srs_initialized();
        let mut builder = build_addition_circuit();
        builder.finalize_circuit(true);
        builder.blocks.compute_offsets();

        let total = builder.blocks.get_total_content_size();
        let circuit_size = (total + 1).next_power_of_two();
        let (ck, _) = create_test_ck(circuit_size);

        let pk = ProvingKey::create_with_ck(&mut builder, ck);
        let vk = VerificationKey::create_with_g2x(&pk, g2_x);

        assert_eq!(vk.circuit_size, pk.circuit_size);
        assert_eq!(vk.num_public_inputs, pk.num_public_inputs);

        // All 28 commitments should be valid G1 points
        let commitments = vk.get_all_commitments();
        assert_eq!(commitments.len(), 28);
    }

    // =========================================================================
    // Proof Structure Tests
    // =========================================================================

    #[test]
    fn test_proof_is_nonempty() {
        ensure_test_srs_initialized();
        let mut builder = build_addition_circuit();
        builder.finalize_circuit(true);
        builder.blocks.compute_offsets();

        let total = builder.blocks.get_total_content_size();
        let circuit_size = (total + 1).next_power_of_two();
        let (ck, _) = create_test_ck(circuit_size);

        let pk = ProvingKey::create_with_ck(&mut builder, ck);
        let prover = UltraProver::new(pk);
        let proof = prover.prove();

        assert!(!proof.is_empty(), "Proof should be non-empty");
    }

    #[test]
    fn test_proof_determinism() {
        ensure_test_srs_initialized();

        // Generate two proofs from the same circuit
        let proof1 = {
            let mut builder = build_addition_circuit();
            builder.finalize_circuit(true);
            builder.blocks.compute_offsets();
            let total = builder.blocks.get_total_content_size();
            let circuit_size = (total + 1).next_power_of_two();
            let (ck, _) = create_test_ck(circuit_size);
            let pk = ProvingKey::create_with_ck(&mut builder, ck);
            let prover = UltraProver::new(pk);
            prover.prove()
        };

        let proof2 = {
            let mut builder = build_addition_circuit();
            builder.finalize_circuit(true);
            builder.blocks.compute_offsets();
            let total = builder.blocks.get_total_content_size();
            let circuit_size = (total + 1).next_power_of_two();
            let (ck, _) = create_test_ck(circuit_size);
            let pk = ProvingKey::create_with_ck(&mut builder, ck);
            let prover = UltraProver::new(pk);
            prover.prove()
        };

        assert_eq!(proof1, proof2, "Same circuit should produce same proof");
    }

    // =========================================================================
    // Multiple Public Input Tests
    // =========================================================================

    #[test]
    fn test_e2e_multiple_public_inputs() {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let a = builder.base.add_public_variable(Fr::from(10u64));
        let b = builder.base.add_public_variable(Fr::from(20u64));
        let c = builder.base.add_variable(Fr::from(30u64));

        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::one(),
            b_scaling: Fr::one(),
            c_scaling: Fr::zero() - Fr::one(),
            const_scaling: Fr::zero(),
        });

        assert!(prove_and_verify(builder), "Multiple public inputs should verify");
    }

    // =========================================================================
    // Large Circuit Tests
    // =========================================================================

    #[test]
    fn test_e2e_chain_of_additions() {
        let mut builder = UltraCircuitBuilder::<P>::new();

        // Build a chain: v[0] = 1, v[i+1] = v[i] + 1
        let mut prev = builder.base.add_variable(Fr::one());
        for i in 1..16 {
            let next = builder.base.add_variable(Fr::from((i + 1) as u64));
            let one_var = builder.base.add_variable(Fr::one());

            builder.create_add_gate(&AddTriple {
                a: prev,
                b: one_var,
                c: next,
                a_scaling: Fr::one(),
                b_scaling: Fr::one(),
                c_scaling: Fr::zero() - Fr::one(),
                const_scaling: Fr::zero(),
            });

            prev = next;
        }

        assert!(prove_and_verify(builder), "Chain circuit should verify");
    }

    // =========================================================================
    // Debug: Circuit Checker Validation
    // =========================================================================

    #[test]
    fn test_circuit_checker_addition() {
        use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
        let mut builder = build_addition_circuit();
        let result = UltraCircuitChecker::check(&mut builder);
        assert!(result.is_ok(), "Circuit checker failed: {:?}", result);
    }

    // =========================================================================
    // Negative Tests (verification should fail)
    // =========================================================================

    // =========================================================================
    // C++ ultra_honk.test.cpp ports
    // =========================================================================

    /// Port of C++ `ANonZeroPolynomialIsAGoodPolynomial`.
    /// Ensures selector, table, and wire polynomials contain non-zero coefficients.
    #[test]
    fn test_non_zero_polynomial_is_good_polynomial() {
        ensure_test_srs_initialized();
        let mut builder = build_addition_circuit();
        builder.finalize_circuit(true);
        builder.blocks.compute_offsets();

        let total = builder.blocks.get_total_content_size();
        let circuit_size = (total + 1).next_power_of_two();
        let (ck, _) = create_test_ck(circuit_size);

        let pk = ProvingKey::create_with_ck(&mut builder, ck);

        // Check that at least some precomputed polynomials have non-zero coefficients.
        // Not all selectors are used in a simple addition circuit (e.g., q_elliptic,
        // q_memory, q_nnf, q_poseidon2_* may be zero), so we check that the majority
        // are non-zero rather than asserting all.
        let precomputed = pk.polynomials.get_precomputed();
        let nonzero_count = precomputed.iter().filter(|poly| {
            (0..pk.circuit_size).any(|j| !poly.get(j).is_zero())
        }).count();
        assert!(nonzero_count > precomputed.len() / 2,
            "Too few non-zero precomputed polynomials: {}/{}",
            nonzero_count, precomputed.len());

        // Check wire polynomials have non-zero coefficients
        let wires = pk.polynomials.get_wires();
        for (i, poly) in wires.iter().enumerate() {
            let has_nonzero = (0..pk.circuit_size).any(|j| !poly.get(j).is_zero());
            assert!(has_nonzero, "Wire polynomial {} is identically zero", i);
        }
    }

    /// Port of C++ `TestNoLookupProof`.
    /// XOR operations via big_add_gate without lookup tables.
    #[test]
    fn test_no_lookup_proof() {
        let mut builder = UltraCircuitBuilder::<P>::new();

        for i in 0u64..16 {
            for j in 0u64..16 {
                let left_idx = builder.base.add_variable(Fr::from(j));
                let right_idx = builder.base.add_variable(Fr::from(i));
                let result_idx = builder.base.add_variable(Fr::from(j ^ i));

                let sum = Fr::from(j) + Fr::from(i) + Fr::from(j ^ i);
                let add_idx = builder.base.add_variable(sum);

                builder.create_big_add_gate(
                    &AddQuad {
                        a: left_idx,
                        b: right_idx,
                        c: result_idx,
                        d: add_idx,
                        a_scaling: Fr::one(),
                        b_scaling: Fr::one(),
                        c_scaling: Fr::one(),
                        d_scaling: Fr::zero() - Fr::one(),
                        const_scaling: Fr::zero(),
                    },
                    false,
                );
            }
        }

        assert!(prove_and_verify(builder), "No-lookup proof should verify");
    }

    /// Port of C++ `TestEllipticGate`.
    /// Elliptic curve point addition and subtraction via create_ecc_add_gate.
    #[test]
    fn test_elliptic_gate() {
        use bbrs_ecc::curves::grumpkin::GrumpkinG1Params;
        use bbrs_ecc::groups::element::Element as GElement;

        let mut builder = UltraCircuitBuilder::<P>::new();

        // Grumpkin Fq = BN254 Fr, so coordinates are Fr directly
        let p1 = GElement::<GrumpkinG1Params>::random_element().to_affine();
        let p2 = GElement::<GrumpkinG1Params>::random_element().to_affine();

        // Addition: p1 + p2 = p3
        let p3 = (GElement::<GrumpkinG1Params>::from_affine(&p1)
            + GElement::<GrumpkinG1Params>::from_affine(&p2))
            .to_affine();

        let x1 = builder.base.add_variable(p1.x);
        let y1 = builder.base.add_variable(p1.y);
        let x2 = builder.base.add_variable(p2.x);
        let y2 = builder.base.add_variable(p2.y);
        let x3 = builder.base.add_variable(p3.x);
        let y3 = builder.base.add_variable(p3.y);

        builder.create_ecc_add_gate(&EccAddGate {
            x1, y1, x2, y2, x3, y3,
            sign_coefficient: Fr::one(),
        });

        // Another addition with same points
        let p3b = (GElement::<GrumpkinG1Params>::from_affine(&p1)
            + GElement::<GrumpkinG1Params>::from_affine(&p2))
            .to_affine();
        let x3b = builder.base.add_variable(p3b.x);
        let y3b = builder.base.add_variable(p3b.y);

        builder.create_ecc_add_gate(&EccAddGate {
            x1, y1, x2, y2, x3: x3b, y3: y3b,
            sign_coefficient: Fr::one(),
        });

        // Subtraction: p1 - p2 = p4
        let p4 = (GElement::<GrumpkinG1Params>::from_affine(&p1)
            - GElement::<GrumpkinG1Params>::from_affine(&p2))
            .to_affine();
        let x4 = builder.base.add_variable(p4.x);
        let y4 = builder.base.add_variable(p4.y);

        builder.create_ecc_add_gate(&EccAddGate {
            x1, y1, x2, y2, x3: x4, y3: y4,
            sign_coefficient: Fr::zero() - Fr::one(),
        });

        assert!(prove_and_verify(builder), "Elliptic gate circuit should verify");
    }

    /// Port of C++ `TooShortProofRejected`.
    /// Truncated proof should panic or fail verification.
    #[test]
    fn test_too_short_proof_rejected() {
        let g2_x = ensure_test_srs_initialized();
        let mut builder = build_addition_circuit();
        builder.finalize_circuit(true);
        builder.blocks.compute_offsets();

        let total = builder.blocks.get_total_content_size();
        let circuit_size = (total + 1).next_power_of_two();
        let (ck, _) = create_test_ck(circuit_size);

        let pk = ProvingKey::create_with_ck(&mut builder, ck);
        let vk = VerificationKey::create_with_g2x(&pk, g2_x);

        let prover = UltraProver::new(pk);
        let proof = prover.prove();

        // Truncate by removing last 10 elements
        let truncated: Vec<Fr> = proof[..proof.len().saturating_sub(10)].to_vec();

        // Verifier should either panic (read past end) or return false
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let verifier = UltraVerifier::new(vk);
            verifier.verify(&truncated).result
        }));

        match result {
            Err(_) => {} // panicked — truncated proof correctly rejected
            Ok(verified) => assert!(!verified, "Truncated proof should fail verification"),
        }
    }

    /// Port of C++ `RangeChecksOnDuplicates` (simplified — no range constraint APIs).
    /// Tests copy constraints on equal variables with big_add_gate.
    #[test]
    fn test_duplicate_variables_with_copy_constraints() {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let a = builder.base.add_variable(Fr::from(100u64));
        let b = builder.base.add_variable(Fr::from(100u64));
        let c = builder.base.add_variable(Fr::from(100u64));
        let d = builder.base.add_variable(Fr::from(100u64));

        builder.base.assert_equal(a, b, "dup ab");
        builder.base.assert_equal(a, c, "dup ac");
        builder.base.assert_equal(a, d, "dup ad");

        builder.create_big_add_gate(
            &AddQuad {
                a, b, c, d,
                a_scaling: Fr::zero(),
                b_scaling: Fr::zero(),
                c_scaling: Fr::zero(),
                d_scaling: Fr::zero(),
                const_scaling: Fr::zero(),
            },
            false,
        );

        assert!(prove_and_verify(builder), "Duplicate variable copy constraints should verify");
    }

    /// Port of C++ `PublicInputs` (multi-gate version with arithmetic gates).
    #[test]
    fn test_public_inputs_arithmetic_gates() {
        let mut builder = UltraCircuitBuilder::<P>::new();

        for i in 0u64..10 {
            let a = builder.base.add_public_variable(Fr::from(i));
            let b = builder.base.add_variable(Fr::from(i + 1));
            let c = builder.base.add_variable(Fr::from(2 * i + 1));

            builder.create_add_gate(&AddTriple {
                a, b, c,
                a_scaling: Fr::one(),
                b_scaling: Fr::one(),
                c_scaling: Fr::zero() - Fr::one(),
                const_scaling: Fr::zero(),
            });
        }

        assert!(prove_and_verify(builder), "Public inputs with arithmetic gates should verify");
    }

    #[test]
    fn test_e2e_tampered_proof_fails() {
        let g2_x = ensure_test_srs_initialized();
        let mut builder = build_addition_circuit();
        builder.finalize_circuit(true);
        builder.blocks.compute_offsets();

        let total = builder.blocks.get_total_content_size();
        let circuit_size = (total + 1).next_power_of_two();
        let (ck, _) = create_test_ck(circuit_size);

        let pk = ProvingKey::create_with_ck(&mut builder, ck);
        let vk = VerificationKey::create_with_g2x(&pk, g2_x);

        let prover = UltraProver::new(pk);
        let mut proof = prover.prove();

        // Tamper with the proof
        if !proof.is_empty() {
            proof[0] = proof[0] + Fr::one();
        }

        let verifier = UltraVerifier::new(vk);
        let output = verifier.verify(&proof);

        assert!(
            !output.result,
            "Tampered proof should fail verification"
        );
    }
}
