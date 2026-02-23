//! Port of `barretenberg/stdlib/verifiers/chonk_verifier` — recursive Chonk (ClientIVC) verifier.
//!
//! Implements the in-circuit verification of Chonk proofs. The Chonk system is
//! the top-level IVC (Incrementally Verifiable Computation) composition that wraps
//! the Goblin recursive verifier and adds fold proof verification for accumulating
//! multiple IVC steps.
//!
//! ## Architecture
//!
//! ```text
//! ChonkRecursiveVerifier
//! ├── GoblinRecursiveVerifier  (composite proof — B28 ✅)
//! │   ├── UltraRecursiveVerifier  (Honk proof — B26 ✅)
//! │   ├── ECCVMRecursiveVerifier  (ECCVM proof — B27 stub)
//! │   ├── TranslatorRecursiveVerifier  (Translator proof — B29 stub)
//! │   └── MergeVerifier  (merge protocol)
//! └── FoldVerifier  (HyperNova fold proof — stub, depends on native chonk)
//! ```
//!
//! ## Current Limitations
//!
//! - Native `chonk` module (ClientIVC) is not yet ported; depends on `dsl`,
//!   `hypernova`, and `multilinear_batching`.
//! - Fold proof verification is stubbed — returns success until HyperNova
//!   folding is implemented.
//! - IVC accumulator state management is deferred.
//!
//! C++ source: `barretenberg/stdlib/verifiers/chonk_verifier.hpp`

use bbrs_ecc::curves::bn254::Fr;

use crate::goblin_verifier::{
    GoblinProof, GoblinRecursiveVerifier, GoblinRecursiveVerifierOutput, GoblinVerificationKeys,
};
use crate::honk_verifier::PairingPointsAccumulator;
use crate::primitives::witness::BuilderRef;

// ════════════════════════════════════════════════════════════════════════════
//  Type aliases
// ════════════════════════════════════════════════════════════════════════════

/// Builder context reference type.
type BuilderCtx = BuilderRef<bbrs_ecc::curves::bn254::Bn254FrParams>;

// ════════════════════════════════════════════════════════════════════════════
//  ChonkProof — serialized Chonk (ClientIVC) proof components
// ════════════════════════════════════════════════════════════════════════════

/// Serialized Chonk proof containing the Goblin proof and IVC fold proof.
///
/// Port of C++ `ClientIVC::Proof` — bundles the Goblin composite proof with
/// the HyperNova fold proof that links consecutive IVC steps.
pub struct ChonkProof {
    /// The Goblin composite proof (Honk + ECCVM + Translator + merge).
    pub goblin_proof: GoblinProof,
    /// The fold proof from HyperNova accumulation (opaque until native chonk is ported).
    pub fold_proof: Vec<Fr>,
}

impl ChonkProof {
    /// Create an empty Chonk proof (for testing).
    pub fn empty() -> Self {
        Self {
            goblin_proof: GoblinProof::empty(),
            fold_proof: Vec::new(),
        }
    }

    /// Check whether this proof has any content.
    pub fn is_empty(&self) -> bool {
        self.goblin_proof.is_empty() && self.fold_proof.is_empty()
    }

    /// Create a Chonk proof from a Goblin proof (no fold proof).
    ///
    /// Useful for testing when only the Goblin layer is active and
    /// fold verification is stubbed.
    pub fn from_goblin_proof(goblin_proof: GoblinProof) -> Self {
        Self {
            goblin_proof,
            fold_proof: Vec::new(),
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  ChonkVerificationKeys
// ════════════════════════════════════════════════════════════════════════════

/// Collection of verification keys needed for Chonk recursive verification.
///
/// Port of C++ `ClientIVC::VerificationKeys` — groups the Goblin VKs with
/// the fold verification key for HyperNova accumulation.
pub struct ChonkVerificationKeys {
    /// Verification keys for the Goblin sub-verifiers.
    pub goblin_vks: GoblinVerificationKeys,
    /// Fold verification key for HyperNova (placeholder until native chonk is ported).
    pub fold_vk: Option<FoldVerificationKeyStub>,
}

impl ChonkVerificationKeys {
    /// Create from Goblin verification keys only (for current testing).
    pub fn from_goblin_vks(goblin_vks: GoblinVerificationKeys) -> Self {
        Self {
            goblin_vks,
            fold_vk: None,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  Fold verification stub
// ════════════════════════════════════════════════════════════════════════════

/// Stub verification key for the HyperNova fold verifier.
///
/// Placeholder until the native `chonk` module (which depends on `hypernova`
/// and `multilinear_batching`) is ported.
pub struct FoldVerificationKeyStub {
    /// Placeholder accumulator size for the folding circuit.
    pub accumulator_size: usize,
}

/// Output of the fold proof verification.
///
/// Port of the HyperNova fold verifier output. When the native `chonk`
/// module is implemented, this will contain the verified IVC accumulator
/// state and any pairing points from the fold verification.
pub struct FoldVerifierOutput {
    /// Whether the fold proof verified successfully.
    pub verified: bool,
    /// Pairing points from fold verification (if any).
    pub points_accumulator: PairingPointsAccumulator,
}

impl FoldVerifierOutput {
    /// Create a stub output that indicates deferred verification.
    fn deferred() -> Self {
        Self {
            verified: true,
            points_accumulator: PairingPointsAccumulator::new(),
        }
    }
}

/// Verify the HyperNova fold proof in-circuit (stub).
///
/// Port of C++ `ClientIVC::FoldVerifier::verify()`.
///
/// When the native `chonk` module is implemented, this will:
/// 1. Load the fold proof into an in-circuit transcript
/// 2. Verify the HyperNova accumulation step
/// 3. Check the IVC accumulator state transition
/// 4. Return pairing points for aggregation
///
/// Currently stubbed — returns success until HyperNova is ported.
fn verify_fold_proof(
    _builder: &BuilderCtx,
    _fold_proof: &[Fr],
    _fold_vk: &Option<FoldVerificationKeyStub>,
) -> FoldVerifierOutput {
    // TODO(chonk): Implement fold proof verification when native `chonk`,
    // `hypernova`, and `multilinear_batching` modules are ported.
    //
    // The fold verifier checks:
    // 1. Accumulator consistency: the running IVC accumulator from the
    //    previous step is correctly updated.
    // 2. Relation check: the folded relation holds for the new accumulator.
    // 3. Commitment chain: commitments to the accumulator state are
    //    properly linked across IVC steps.
    FoldVerifierOutput::deferred()
}

// ════════════════════════════════════════════════════════════════════════════
//  ChonkRecursiveVerifierOutput
// ════════════════════════════════════════════════════════════════════════════

/// Output of the Chonk recursive verifier.
///
/// Port of C++ `ClientIVC::RecursiveVerifierOutput`.
///
/// Contains the aggregated pairing points from the Goblin verifier and fold
/// verifier, along with verification status for each component.
pub struct ChonkRecursiveVerifierOutput {
    /// Aggregated pairing points combining Goblin and fold outputs.
    pub points_accumulator: PairingPointsAccumulator,
    /// Output from the Goblin recursive verifier.
    pub goblin_output: GoblinRecursiveVerifierOutput,
    /// Whether the fold proof verified successfully.
    pub fold_verified: bool,
}

impl ChonkRecursiveVerifierOutput {
    /// Check whether all verifications passed (Goblin + fold).
    pub fn all_verified(&self) -> bool {
        self.goblin_output.all_verified() && self.fold_verified
    }

    /// Create a failed output (all verifications failed).
    pub fn failed() -> Self {
        Self {
            points_accumulator: PairingPointsAccumulator::new(),
            goblin_output: GoblinRecursiveVerifierOutput::failed(),
            fold_verified: false,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  ChonkRecursiveVerifier
// ════════════════════════════════════════════════════════════════════════════

/// Recursive Chonk (ClientIVC) verifier for in-circuit IVC proof verification.
///
/// Port of C++ `ClientIVC::RecursiveVerifier`.
///
/// Orchestrates the verification of a complete Chonk proof by:
/// 1. Running the Goblin recursive verifier on the composite proof
/// 2. Running the fold verifier on the HyperNova accumulation proof
/// 3. Aggregating pairing points from both verifiers
///
/// The final aggregated pairing points can be checked outside the circuit
/// to confirm the validity of the entire IVC composition.
pub struct ChonkRecursiveVerifier {
    /// Builder for the outer (verifier) circuit.
    pub builder: BuilderCtx,
    /// Verification keys for all sub-verifiers.
    pub vks: ChonkVerificationKeys,
}

impl ChonkRecursiveVerifier {
    /// Create a new Chonk recursive verifier.
    pub fn new(builder: BuilderCtx, vks: ChonkVerificationKeys) -> Self {
        Self { builder, vks }
    }

    /// Verify a Chonk proof, constructing verification constraints in the outer circuit.
    ///
    /// Port of C++ `ClientIVC::RecursiveVerifier::verify_proof`.
    ///
    /// This method orchestrates the full Chonk verification:
    /// 1. Verifies the Goblin composite proof (Honk + ECCVM + Translator + merge)
    /// 2. Verifies the fold proof (HyperNova accumulation step)
    /// 3. Aggregates all pairing points
    ///
    /// Returns the combined output with aggregated pairing points.
    pub fn verify_proof(self, chonk_proof: &ChonkProof) -> ChonkRecursiveVerifierOutput {
        // Step 1: Verify the Goblin composite proof
        let goblin_verifier = GoblinRecursiveVerifier::new(
            self.builder.clone(),
            self.vks.goblin_vks,
        );
        let goblin_output = goblin_verifier.verify_proof(&chonk_proof.goblin_proof);

        // Destructure Goblin output to move the accumulator out
        let honk_verified = goblin_output.honk_verified;
        let eccvm_verified = goblin_output.eccvm_verified;
        let translator_verified = goblin_output.translator_verified;
        let merge_verified = goblin_output.merge_verified;
        let mut points_accumulator = goblin_output.points_accumulator;

        // Step 2: Verify the fold proof (stub — depends on native chonk)
        let fold_output = verify_fold_proof(
            &self.builder,
            &chonk_proof.fold_proof,
            &self.vks.fold_vk,
        );
        let fold_verified = fold_output.verified;

        // Step 3: Aggregate fold pairing points into the accumulator
        if fold_verified {
            points_accumulator.aggregate(&fold_output.points_accumulator);
        }

        ChonkRecursiveVerifierOutput {
            points_accumulator,
            goblin_output: GoblinRecursiveVerifierOutput {
                points_accumulator: PairingPointsAccumulator::new(),
                honk_verified,
                eccvm_verified,
                translator_verified,
                merge_verified,
            },
            fold_verified,
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

    // ── Test 1: Chonk verification with Honk-only Goblin proof ──────

    #[test]
    fn test_chonk_verifier_honk_only() {
        let g2_x = ensure_test_srs_initialized();
        let inner_circuit = build_inner_circuit(16);
        let (proof, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        let outer_builder = make_builder_ref();
        let goblin_vks = GoblinVerificationKeys::from_honk_vk(native_vk);
        let chonk_vks = ChonkVerificationKeys::from_goblin_vks(goblin_vks);

        let chonk_proof = ChonkProof::from_goblin_proof(GoblinProof {
            honk_proof: proof,
            eccvm_proof: Vec::new(),
            translator_proof: Vec::new(),
            merge_proof: Vec::new(),
        });

        let verifier = ChonkRecursiveVerifier::new(outer_builder.clone(), chonk_vks);
        let output = verifier.verify_proof(&chonk_proof);

        assert!(
            output.goblin_output.honk_verified,
            "Honk verification should pass through Chonk"
        );
        assert!(
            output.goblin_output.eccvm_verified,
            "ECCVM stub should pass (deferred)"
        );
        assert!(
            output.goblin_output.translator_verified,
            "Translator stub should pass (deferred)"
        );
        assert!(
            output.goblin_output.merge_verified,
            "Merge protocol stub should pass (deferred)"
        );
        assert!(
            output.fold_verified,
            "Fold verification stub should pass (deferred)"
        );
        assert!(
            output.all_verified(),
            "All sub-verifications should pass"
        );

        // Verify aggregated pairing points are non-trivial
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

        // Verify outer circuit has witnesses from verification
        let num_vars = outer_builder.borrow().base.get_num_variables();
        assert!(
            num_vars > 1,
            "Chonk verification should create witnesses in outer circuit"
        );
    }

    // ── Test 2: ChonkRecursiveVerifierOutput status checks ──────────

    #[test]
    fn test_chonk_output_verification_status() {
        let failed = ChonkRecursiveVerifierOutput::failed();
        assert!(
            !failed.all_verified(),
            "Failed output should not report all_verified"
        );
        assert!(!failed.goblin_output.honk_verified);
        assert!(!failed.goblin_output.eccvm_verified);
        assert!(!failed.goblin_output.translator_verified);
        assert!(!failed.goblin_output.merge_verified);
        assert!(!failed.fold_verified);

        // Partial: Goblin passes but fold fails
        let partial_fold_fail = ChonkRecursiveVerifierOutput {
            points_accumulator: PairingPointsAccumulator::new(),
            goblin_output: GoblinRecursiveVerifierOutput {
                points_accumulator: PairingPointsAccumulator::new(),
                honk_verified: true,
                eccvm_verified: true,
                translator_verified: true,
                merge_verified: true,
            },
            fold_verified: false,
        };
        assert!(
            !partial_fold_fail.all_verified(),
            "Should not be all_verified when fold fails"
        );

        // Partial: fold passes but Goblin partially fails
        let partial_goblin_fail = ChonkRecursiveVerifierOutput {
            points_accumulator: PairingPointsAccumulator::new(),
            goblin_output: GoblinRecursiveVerifierOutput {
                points_accumulator: PairingPointsAccumulator::new(),
                honk_verified: true,
                eccvm_verified: false,
                translator_verified: true,
                merge_verified: true,
            },
            fold_verified: true,
        };
        assert!(
            !partial_goblin_fail.all_verified(),
            "Should not be all_verified when Goblin sub-verifier fails"
        );

        // Full success
        let success = ChonkRecursiveVerifierOutput {
            points_accumulator: PairingPointsAccumulator::new(),
            goblin_output: GoblinRecursiveVerifierOutput {
                points_accumulator: PairingPointsAccumulator::new(),
                honk_verified: true,
                eccvm_verified: true,
                translator_verified: true,
                merge_verified: true,
            },
            fold_verified: true,
        };
        assert!(
            success.all_verified(),
            "Full success should report all_verified"
        );

        // ChonkProof construction helpers
        let empty = ChonkProof::empty();
        assert!(empty.is_empty(), "Empty proof should report is_empty()");

        let non_empty = ChonkProof::from_goblin_proof(GoblinProof {
            honk_proof: vec![Fr::one()],
            eccvm_proof: Vec::new(),
            translator_proof: Vec::new(),
            merge_proof: Vec::new(),
        });
        assert!(
            !non_empty.is_empty(),
            "Proof with Goblin data should not be empty"
        );
    }
}
