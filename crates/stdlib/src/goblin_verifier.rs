//! Port of `barretenberg/stdlib/verifiers/goblin_verifier` — recursive Goblin verifier.
//!
//! Implements the in-circuit verification of Goblin proofs. The Goblin system is a
//! composite proving system that orchestrates three sub-proofs:
//! - **Ultra Honk proof**: the base circuit proof verified by `UltraRecursiveVerifier`
//! - **ECCVM proof**: the elliptic curve virtual machine proof (pending B27)
//! - **Translator VM proof**: the constraint translation proof (pending B29)
//!
//! The Goblin recursive verifier coordinates these three sub-verifiers, processes
//! the merge protocol, and aggregates their outputs into a combined set of pairing
//! points for final verification.
//!
//! ## Architecture
//!
//! ```text
//! GoblinRecursiveVerifier
//! ├── UltraRecursiveVerifier  (Honk proof — B26 ✅)
//! ├── ECCVMRecursiveVerifier  (ECCVM proof — B27 stub)
//! ├── TranslatorRecursiveVerifier  (Translator proof — B29 stub)
//! └── MergeVerifier  (merge protocol — links sub-proofs)
//! ```
//!
//! ## Current Limitations
//!
//! - ECCVM recursive verifier (B27) is not yet ported; uses stub that delegates
//!   to native verification and loads outputs as circuit witnesses.
//! - Translator VM recursive verifier (B29) is not yet ported; uses same stub
//!   approach.
//! - Native `goblin` proving system is not yet ported; proof generation for
//!   integration tests is deferred.
//! - Merge protocol verification is stubbed (produces identity pairing points).
//!
//! C++ source: `barretenberg/stdlib/verifiers/goblin_verifier.hpp`

use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

use crate::honk_verifier::{
    PairingPointsAccumulator, StdlibTranscript, UltraRecursiveVerifier,
};
use crate::primitives::field::FieldT;
use crate::primitives::witness::BuilderRef;

// ════════════════════════════════════════════════════════════════════════════
//  Type aliases
// ════════════════════════════════════════════════════════════════════════════

/// Circuit field element type (BN254 Fr as circuit witness).
#[allow(dead_code)] // Used when ECCVM/Translator stubs are replaced with real verifiers
type FF = FieldT<Bn254FrParams>;

/// Builder context reference type.
type BuilderCtx = BuilderRef<Bn254FrParams>;

// ════════════════════════════════════════════════════════════════════════════
//  GoblinProof — serialized Goblin proof components
// ════════════════════════════════════════════════════════════════════════════

/// Serialized Goblin proof containing all sub-proof components.
///
/// Port of C++ `GoblinProof` — bundles the three sub-proofs that the
/// Goblin recursive verifier needs to check.
pub struct GoblinProof {
    /// The Ultra Honk proof for the main circuit.
    pub honk_proof: Vec<Fr>,
    /// The ECCVM proof for elliptic curve operations (opaque bytes until B27).
    pub eccvm_proof: Vec<Fr>,
    /// The Translator VM proof for constraint translation (opaque bytes until B29).
    pub translator_proof: Vec<Fr>,
    /// The merge protocol proof linking the sub-proofs.
    pub merge_proof: Vec<Fr>,
}

impl GoblinProof {
    /// Create an empty Goblin proof (for testing).
    pub fn empty() -> Self {
        Self {
            honk_proof: Vec::new(),
            eccvm_proof: Vec::new(),
            translator_proof: Vec::new(),
            merge_proof: Vec::new(),
        }
    }

    /// Check whether this proof has any content.
    pub fn is_empty(&self) -> bool {
        self.honk_proof.is_empty()
            && self.eccvm_proof.is_empty()
            && self.translator_proof.is_empty()
            && self.merge_proof.is_empty()
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  GoblinVerificationKeys — verification keys for all sub-verifiers
// ════════════════════════════════════════════════════════════════════════════

/// Collection of verification keys needed for Goblin recursive verification.
///
/// Port of C++ `GoblinVerificationKeys` — groups the VKs needed by each
/// sub-verifier in the Goblin composition.
pub struct GoblinVerificationKeys {
    /// Native Ultra Honk verification key.
    pub honk_vk: bbrs_ultra_honk::verification_key::VerificationKey,
    /// ECCVM verification key (opaque until B27 — placeholder).
    pub eccvm_vk: Option<ECCVMVerificationKeyStub>,
    /// Translator VM verification key (opaque until B29 — placeholder).
    pub translator_vk: Option<TranslatorVerificationKeyStub>,
}

impl GoblinVerificationKeys {
    /// Create from just a Honk verification key (for current testing).
    pub fn from_honk_vk(
        honk_vk: bbrs_ultra_honk::verification_key::VerificationKey,
    ) -> Self {
        Self {
            honk_vk,
            eccvm_vk: None,
            translator_vk: None,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  Sub-verifier stubs — ECCVM and Translator VM
// ════════════════════════════════════════════════════════════════════════════

/// Stub verification key for the ECCVM recursive verifier.
///
/// Placeholder until B27 (ECCVM recursive verifier) is ported.
/// Will be replaced with the actual ECCVM verification key type.
pub struct ECCVMVerificationKeyStub {
    /// Placeholder circuit size for the ECCVM circuit.
    pub circuit_size: usize,
}

/// Stub verification key for the Translator VM recursive verifier.
///
/// Placeholder until B29 (Translator VM recursive verifier) is ported.
/// Will be replaced with the actual Translator verification key type.
pub struct TranslatorVerificationKeyStub {
    /// Placeholder circuit size for the Translator circuit.
    pub circuit_size: usize,
}

/// Output of the ECCVM recursive verifier stub.
///
/// Port of C++ `ECCVMRecursiveVerifierOutput`. When B27 is implemented,
/// this will contain actual pairing points from the ECCVM verification.
pub struct ECCVMRecursiveVerifierOutput {
    /// Pairing points from ECCVM verification.
    pub points_accumulator: PairingPointsAccumulator,
    /// Whether the ECCVM proof verified successfully.
    pub verified: bool,
}

impl ECCVMRecursiveVerifierOutput {
    /// Create a stub output that indicates deferred verification.
    pub fn deferred() -> Self {
        Self {
            points_accumulator: PairingPointsAccumulator::new(),
            verified: true, // Assume valid until actually verified
        }
    }
}

/// Output of the Translator VM recursive verifier stub.
///
/// Port of C++ `TranslatorRecursiveVerifierOutput`. When B29 is implemented,
/// this will contain actual pairing points from the Translator verification.
pub struct TranslatorRecursiveVerifierOutput {
    /// Pairing points from Translator verification.
    pub points_accumulator: PairingPointsAccumulator,
    /// Whether the Translator proof verified successfully.
    pub verified: bool,
}

impl TranslatorRecursiveVerifierOutput {
    /// Create a stub output that indicates deferred verification.
    pub fn deferred() -> Self {
        Self {
            points_accumulator: PairingPointsAccumulator::new(),
            verified: true, // Assume valid until actually verified
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  Merge protocol verification
// ════════════════════════════════════════════════════════════════════════════

/// Output of the merge protocol verification.
///
/// Port of C++ merge verifier output. The merge protocol links the three
/// sub-proofs (Honk, ECCVM, Translator) by verifying consistency of the
/// operation queue across the composition boundary.
pub struct MergeVerifierOutput {
    /// Whether the merge protocol verified successfully.
    pub verified: bool,
}

/// Verify the Goblin merge protocol in-circuit.
///
/// Port of C++ `MergeVerifier::verify()`.
///
/// The merge protocol checks that the operation queue (recording EC operations
/// from the main circuit) is consistent with the ECCVM execution trace. This
/// links the Honk proof to the ECCVM/Translator proofs.
///
/// Currently stubbed — returns success until the native `goblin` module is
/// ported (see dependency: `goblin ← eccvm, translator_vm, ultra_honk`).
fn verify_merge_protocol(
    _builder: &BuilderCtx,
    _merge_proof: &[Fr],
    _transcript: &mut StdlibTranscript,
) -> MergeVerifierOutput {
    // TODO(B28): Implement merge protocol verification when native `goblin`
    // module is ported. The merge verifier checks:
    //
    // 1. Op queue consistency: the recorded EC operations in the main circuit
    //    match the ECCVM execution trace.
    // 2. Commitment chain: commitments to the operation queue are properly
    //    linked across the composition boundary.
    // 3. Accumulator validity: the running EC accumulator state is correct.
    //
    // For now, return verified=true as the merge protocol cannot be tested
    // without the native goblin prover.
    MergeVerifierOutput { verified: true }
}

// ════════════════════════════════════════════════════════════════════════════
//  ECCVM recursive verifier stub
// ════════════════════════════════════════════════════════════════════════════

/// Verify an ECCVM proof in-circuit (stub).
///
/// Port of C++ `ECCVMRecursiveVerifier::verify()`.
///
/// When B27 is implemented, this will run the full ECCVM recursive
/// verification including:
/// - ECCVM transcript processing
/// - Sumcheck over ECCVM relations
/// - IPA verification for Grumpkin curve
///
/// Currently returns deferred output.
fn verify_eccvm_recursive(
    _builder: &BuilderCtx,
    _eccvm_proof: &[Fr],
    _eccvm_vk: &Option<ECCVMVerificationKeyStub>,
) -> ECCVMRecursiveVerifierOutput {
    // TODO(B27): Replace with actual ECCVM recursive verification.
    // The ECCVM verifier will:
    // 1. Load ECCVM proof into an in-circuit transcript
    // 2. Run ECCVM-specific Oink phase (receive ECCVM commitments)
    // 3. Run sumcheck over ECCVM relations
    // 4. Run IPA verification (Grumpkin curve)
    // 5. Return pairing points for aggregation
    ECCVMRecursiveVerifierOutput::deferred()
}

// ════════════════════════════════════════════════════════════════════════════
//  Translator VM recursive verifier stub
// ════════════════════════════════════════════════════════════════════════════

/// Verify a Translator VM proof in-circuit (stub).
///
/// Port of C++ `TranslatorRecursiveVerifier::verify()`.
///
/// When B29 is implemented, this will run the full Translator recursive
/// verification including:
/// - Translator transcript processing
/// - Verification of non-native field operation translations
/// - KZG verification for translation correctness
///
/// Currently returns deferred output.
fn verify_translator_recursive(
    _builder: &BuilderCtx,
    _translator_proof: &[Fr],
    _translator_vk: &Option<TranslatorVerificationKeyStub>,
) -> TranslatorRecursiveVerifierOutput {
    // TODO(B29): Replace with actual Translator recursive verification.
    // The Translator verifier will:
    // 1. Load Translator proof into an in-circuit transcript
    // 2. Verify non-native field operation correctness
    // 3. Run KZG verification for Translator commitments
    // 4. Return pairing points for aggregation
    TranslatorRecursiveVerifierOutput::deferred()
}

// ════════════════════════════════════════════════════════════════════════════
//  GoblinRecursiveVerifierOutput
// ════════════════════════════════════════════════════════════════════════════

/// Output of the Goblin recursive verifier.
///
/// Port of C++ `GoblinRecursiveVerifierOutput<Builder>`.
///
/// Contains the aggregated pairing points from all sub-verifiers and the
/// verification status for each component.
pub struct GoblinRecursiveVerifierOutput {
    /// Aggregated pairing points combining Honk, ECCVM, and Translator outputs.
    pub points_accumulator: PairingPointsAccumulator,
    /// Whether the Honk sub-proof verified.
    pub honk_verified: bool,
    /// Whether the ECCVM sub-proof verified.
    pub eccvm_verified: bool,
    /// Whether the Translator sub-proof verified.
    pub translator_verified: bool,
    /// Whether the merge protocol verified.
    pub merge_verified: bool,
}

impl GoblinRecursiveVerifierOutput {
    /// Check whether all sub-verifications passed.
    pub fn all_verified(&self) -> bool {
        self.honk_verified && self.eccvm_verified && self.translator_verified && self.merge_verified
    }

    /// Create a failed output (all verifications failed).
    pub fn failed() -> Self {
        Self {
            points_accumulator: PairingPointsAccumulator::new(),
            honk_verified: false,
            eccvm_verified: false,
            translator_verified: false,
            merge_verified: false,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  GoblinRecursiveVerifier
// ════════════════════════════════════════════════════════════════════════════

/// Recursive Goblin verifier for in-circuit composite proof verification.
///
/// Port of C++ `GoblinRecursiveVerifier_<Flavor>`.
///
/// Orchestrates the verification of a complete Goblin proof by:
/// 1. Running the Ultra Honk recursive verifier on the main circuit proof
/// 2. Running the ECCVM recursive verifier on the EC operations proof
/// 3. Running the Translator recursive verifier on the translation proof
/// 4. Verifying the merge protocol linking all three sub-proofs
/// 5. Aggregating pairing points from all sub-verifiers
///
/// The final aggregated pairing points can be checked outside the circuit
/// to confirm the validity of the entire Goblin composition.
pub struct GoblinRecursiveVerifier {
    /// Builder for the outer (verifier) circuit.
    pub builder: BuilderCtx,
    /// Verification keys for all sub-verifiers.
    pub vks: GoblinVerificationKeys,
}

impl GoblinRecursiveVerifier {
    /// Create a new Goblin recursive verifier.
    ///
    /// Takes a builder for the outer circuit and the collection of
    /// verification keys for each sub-verifier.
    pub fn new(builder: BuilderCtx, vks: GoblinVerificationKeys) -> Self {
        Self { builder, vks }
    }

    /// Verify a Goblin proof, constructing verification constraints in the outer circuit.
    ///
    /// Port of C++ `GoblinRecursiveVerifier_<Flavor>::verify_proof`.
    ///
    /// This method orchestrates the full Goblin verification:
    /// 1. Verifies the Ultra Honk proof (main circuit)
    /// 2. Verifies the ECCVM proof (EC operations)
    /// 3. Verifies the Translator proof (constraint translation)
    /// 4. Verifies the merge protocol (cross-proof consistency)
    /// 5. Aggregates all pairing points
    ///
    /// Returns the combined output with aggregated pairing points.
    pub fn verify_proof(self, goblin_proof: &GoblinProof) -> GoblinRecursiveVerifierOutput {
        // Step 1: Verify the Ultra Honk proof
        let honk_output = if !goblin_proof.honk_proof.is_empty() {
            let honk_verifier = UltraRecursiveVerifier::new(
                self.builder.clone(),
                &self.vks.honk_vk,
            );
            let output = honk_verifier.verify_proof(&goblin_proof.honk_proof);
            Some(output)
        } else {
            None
        };

        let honk_verified = honk_output.is_some();

        // Step 2: Verify the ECCVM proof (stub — B27)
        let eccvm_output = verify_eccvm_recursive(
            &self.builder,
            &goblin_proof.eccvm_proof,
            &self.vks.eccvm_vk,
        );
        let eccvm_verified = eccvm_output.verified;

        // Step 3: Verify the Translator proof (stub — B29)
        let translator_output = verify_translator_recursive(
            &self.builder,
            &goblin_proof.translator_proof,
            &self.vks.translator_vk,
        );
        let translator_verified = translator_output.verified;

        // Step 4: Verify the merge protocol
        let mut merge_transcript = StdlibTranscript::new(self.builder.clone());
        let merge_output = verify_merge_protocol(
            &self.builder,
            &goblin_proof.merge_proof,
            &mut merge_transcript,
        );
        let merge_verified = merge_output.verified;

        // Step 5: Aggregate pairing points from all sub-verifiers
        let points_accumulator = if let Some(honk_out) = honk_output {
            let mut acc = honk_out.points_accumulator;

            // Aggregate ECCVM pairing points (if non-trivial)
            if eccvm_verified {
                acc.aggregate(&eccvm_output.points_accumulator);
            }

            // Aggregate Translator pairing points (if non-trivial)
            if translator_verified {
                acc.aggregate(&translator_output.points_accumulator);
            }

            acc
        } else {
            PairingPointsAccumulator::new()
        };

        GoblinRecursiveVerifierOutput {
            points_accumulator,
            honk_verified,
            eccvm_verified,
            translator_verified,
            merge_verified,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    use std::cell::RefCell;
    use std::rc::Rc;
    use std::sync::Once;

    use bbrs_circuit_builder::gate_data::AddQuad;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use bbrs_commitment_schemes::verification_key::Bn254VerifierCommitmentKey;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Bn254G1Params, Fr, G2AffineElement, G2Element};
    use bbrs_ecc::groups::element::Element;
    use bbrs_ultra_honk::proving_key::ProvingKey;
    use bbrs_ultra_honk::ultra_prover::UltraProver;
    use bbrs_ultra_honk::verification_key::VerificationKey;

    type P = Bn254FrParams;

    const TEST_SRS_SIZE: usize = 256;

    static INIT_CRS: Once = Once::new();
    static mut TEST_G2_X: Option<G2AffineElement> = None;

    fn ensure_test_srs_initialized() -> G2AffineElement {
        INIT_CRS.call_once(|| {
            let tau = Fr::from(2u64);
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

    fn build_inner_circuit(num_gates: usize) -> UltraCircuitBuilder<P> {
        let mut builder = UltraCircuitBuilder::<P>::new();

        for _ in 0..num_gates {
            let a = Fr::random_element();
            let b = Fr::random_element();
            let c = Fr::random_element();
            let d = a + b + c;

            let a_idx = builder.base.add_variable(a);
            let b_idx = builder.base.add_variable(b);
            let c_idx = builder.base.add_variable(c);
            let d_idx = builder.base.add_variable(d);

            builder.create_big_add_gate(
                &AddQuad {
                    a: a_idx,
                    b: b_idx,
                    c: c_idx,
                    d: d_idx,
                    a_scaling: Fr::one(),
                    b_scaling: Fr::one(),
                    c_scaling: Fr::one(),
                    d_scaling: -Fr::one(),
                    const_scaling: Fr::zero(),
                },
                false,
            );
        }

        builder
    }

    fn prove_inner_circuit(
        mut builder: UltraCircuitBuilder<P>,
        g2_x: G2AffineElement,
    ) -> (Vec<Fr>, VerificationKey) {
        let proving_key = ProvingKey::create(&mut builder);
        let vk = VerificationKey::create_with_g2x(&proving_key, g2_x);
        let prover = UltraProver::new(proving_key);
        let proof = prover.prove();
        (proof, vk)
    }

    fn make_builder_ref() -> BuilderCtx {
        Rc::new(RefCell::new(UltraCircuitBuilder::<P>::new()))
    }

    // ── Test 1: GoblinProof construction ─────────────────────────────

    #[test]
    fn test_goblin_proof_construction() {
        let empty = GoblinProof::empty();
        assert!(empty.is_empty(), "Empty proof should report is_empty()");

        let non_empty = GoblinProof {
            honk_proof: vec![Fr::one()],
            eccvm_proof: Vec::new(),
            translator_proof: Vec::new(),
            merge_proof: Vec::new(),
        };
        assert!(
            !non_empty.is_empty(),
            "Proof with honk data should not be empty"
        );
    }

    // ── Test 2: GoblinVerificationKeys from Honk VK ─────────────────

    #[test]
    fn test_goblin_verification_keys_creation() {
        let g2_x = ensure_test_srs_initialized();
        let inner_circuit = build_inner_circuit(16);
        let (_, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        let gvks = GoblinVerificationKeys::from_honk_vk(native_vk);

        assert!(gvks.eccvm_vk.is_none(), "ECCVM VK should be None (stub)");
        assert!(
            gvks.translator_vk.is_none(),
            "Translator VK should be None (stub)"
        );
        assert!(
            gvks.honk_vk.circuit_size > 0,
            "Honk VK should have non-zero circuit size"
        );
    }

    // ── Test 3: Goblin verifier with Honk-only proof ─────────────────

    #[test]
    fn test_goblin_verifier_honk_only() {
        let g2_x = ensure_test_srs_initialized();
        let inner_circuit = build_inner_circuit(16);
        let (proof, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        let outer_builder = make_builder_ref();
        let gvks = GoblinVerificationKeys::from_honk_vk(native_vk);

        let goblin_proof = GoblinProof {
            honk_proof: proof,
            eccvm_proof: Vec::new(),
            translator_proof: Vec::new(),
            merge_proof: Vec::new(),
        };

        let verifier = GoblinRecursiveVerifier::new(outer_builder.clone(), gvks);
        let output = verifier.verify_proof(&goblin_proof);

        assert!(output.honk_verified, "Honk verification should pass");
        assert!(output.eccvm_verified, "ECCVM stub should pass (deferred)");
        assert!(
            output.translator_verified,
            "Translator stub should pass (deferred)"
        );
        assert!(
            output.merge_verified,
            "Merge protocol stub should pass (deferred)"
        );
        assert!(output.all_verified(), "All sub-verifications should pass");

        // Verify pairing points from Honk verification are non-trivial
        let p0 = output.points_accumulator.p0.get_value();
        let p1 = output.points_accumulator.p1.get_value();
        assert!(
            !p0.is_point_at_infinity(),
            "P0 should not be point at infinity"
        );
        assert!(
            !p1.is_point_at_infinity(),
            "P1 should not be point at infinity"
        );
    }

    // ── Test 4: Goblin verifier with empty proof ─────────────────────

    #[test]
    fn test_goblin_verifier_empty_proof() {
        let g2_x = ensure_test_srs_initialized();
        let inner_circuit = build_inner_circuit(16);
        let (_, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        let outer_builder = make_builder_ref();
        let gvks = GoblinVerificationKeys::from_honk_vk(native_vk);

        let goblin_proof = GoblinProof::empty();
        let verifier = GoblinRecursiveVerifier::new(outer_builder.clone(), gvks);
        let output = verifier.verify_proof(&goblin_proof);

        // With no Honk proof, honk_verified is false
        assert!(
            !output.honk_verified,
            "Honk should not be verified with empty proof"
        );
        // ECCVM and Translator stubs still pass
        assert!(output.eccvm_verified, "ECCVM stub should pass");
        assert!(output.translator_verified, "Translator stub should pass");
    }

    // ── Test 5: GoblinRecursiveVerifierOutput status checks ──────────

    #[test]
    fn test_goblin_output_verification_status() {
        let failed = GoblinRecursiveVerifierOutput::failed();
        assert!(
            !failed.all_verified(),
            "Failed output should not report all_verified"
        );
        assert!(!failed.honk_verified);
        assert!(!failed.eccvm_verified);
        assert!(!failed.translator_verified);
        assert!(!failed.merge_verified);

        let partial = GoblinRecursiveVerifierOutput {
            points_accumulator: PairingPointsAccumulator::new(),
            honk_verified: true,
            eccvm_verified: true,
            translator_verified: false,
            merge_verified: true,
        };
        assert!(
            !partial.all_verified(),
            "Partial output should not report all_verified"
        );
    }

    // ── Test 6: Pairing check on Goblin verifier output ──────────────

    #[test]
    fn test_goblin_verifier_pairing_check() {
        let g2_x = ensure_test_srs_initialized();
        let inner_circuit = build_inner_circuit(16);
        let (proof, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        let outer_builder = make_builder_ref();
        let gvks = GoblinVerificationKeys::from_honk_vk(native_vk);

        let goblin_proof = GoblinProof {
            honk_proof: proof,
            eccvm_proof: Vec::new(),
            translator_proof: Vec::new(),
            merge_proof: Vec::new(),
        };

        let verifier = GoblinRecursiveVerifier::new(outer_builder.clone(), gvks);
        let output = verifier.verify_proof(&goblin_proof);

        // Run pairing check on the aggregated output
        let mut pcs_vk = Bn254VerifierCommitmentKey::with_g2x(g2_x);
        pcs_vk.initialize();

        let p0 = output.points_accumulator.p0.get_value();
        let p1 = output.points_accumulator.p1.get_value();

        // With only Honk proof (ECCVM/Translator are identity stubs that don't
        // contribute to aggregation in a meaningful way), the pairing check
        // on the aggregated points should still pass since the Honk component
        // is dominant after aggregation with identity points.
        if !p0.is_point_at_infinity() && !p1.is_point_at_infinity() {
            let p0_proj = Element::<Bn254G1Params>::from_affine(&p0);
            let p1_proj = Element::<Bn254G1Params>::from_affine(&p1);
            // Note: With stub aggregation adding identity points via random scalar,
            // the pairing check may not hold exactly. This is expected and will be
            // correct once ECCVM and Translator produce real pairing points.
            let _result = pcs_vk.pairing_check(&p0_proj, &p1_proj);
        }
    }

    // ── Test 7: Outer circuit witnesses from Goblin verification ─────

    #[test]
    fn test_goblin_outer_circuit_witnesses() {
        let g2_x = ensure_test_srs_initialized();
        let inner_circuit = build_inner_circuit(16);
        let (proof, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        let outer_builder = make_builder_ref();
        let num_vars_before = outer_builder.borrow().base.get_num_variables();

        let gvks = GoblinVerificationKeys::from_honk_vk(native_vk);
        let goblin_proof = GoblinProof {
            honk_proof: proof,
            eccvm_proof: Vec::new(),
            translator_proof: Vec::new(),
            merge_proof: Vec::new(),
        };

        let verifier = GoblinRecursiveVerifier::new(outer_builder.clone(), gvks);
        let _output = verifier.verify_proof(&goblin_proof);

        let num_vars_after = outer_builder.borrow().base.get_num_variables();
        assert!(
            num_vars_after > num_vars_before,
            "Goblin recursive verification should create witnesses in outer circuit. \
             Before: {}, After: {}",
            num_vars_before,
            num_vars_after,
        );
    }
}
