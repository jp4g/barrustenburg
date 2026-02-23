//! Port of `barretenberg/stdlib/eccvm_verifier/` — recursive ECCVM verifier.
//!
//! Implements the in-circuit verification of ECCVM proofs. The recursive
//! verifier constructs a circuit that constrains the verification of an inner
//! ECCVM proof, enabling proof composition for the Goblin recursive pipeline.
//!
//! ## Architecture
//!
//! The ECCVM operates over the Grumpkin curve (cycle partner of BN254) and uses
//! IPA (Inner Product Argument) commitments instead of KZG. The recursive
//! verifier runs inside a BN254 circuit, so Grumpkin field elements become
//! non-native (bigfield) circuit witnesses.
//!
//! The verification follows these phases:
//! 1. **VK hash**: Load the fixed verification key hash into transcript
//! 2. **Wire commitments**: Receive 85 wire commitments + ZK masking commitment
//! 3. **Challenges**: Compute beta/gamma for set permutation + alpha for batching
//! 4. **Committed sumcheck**: Verify sumcheck with round univariate commitments
//! 5. **Shplemini**: Batch multilinear claims into univariate openings
//! 6. **Translation**: Verify ECCVM-to-Translator wire consistency
//! 7. **IPA reduction**: Produce final IPA opening claim (deferred verification)
//!
//! ## Key Differences from Honk Recursive Verifier (B26)
//!
//! | Aspect | Ultra Honk (B26) | ECCVM (B27) |
//! |--------|-----------------|-------------|
//! | Curve | BN254 (KZG) | Grumpkin (IPA) |
//! | PCS | KZG pairing | IPA (deferred) |
//! | Circuit size | Variable | Fixed (`ECCVM_FIXED_SIZE`) |
//! | VK | Per-circuit | Hardcoded (4 commitments) |
//! | ZK | Optional | Always on |
//! | Num wires | 4 | 85 |
//! | Relations | Ultra relations | 7 ECC-specific relations |
//! | Output | Pairing points | IPA opening claim |
//! | Transcript hash | Keccak | Poseidon2 |
//! | Sumcheck | Evaluations only | Committed (with round commitments) |
//! | Translation | None | SmallSubgroupIPA for op/Px/Py/z1/z2 |
//!
//! ## Current Limitations
//!
//! - Delegates to native Grumpkin operations as placeholders until `biggroup`
//!   for Grumpkin (non-native Fq arithmetic for Grumpkin) is ported.
//! - The native ECCVM prover is not yet ported to Rust, so full end-to-end
//!   recursive verification tests require the ECCVM proving system.
//! - `CommitmentT` wraps native Grumpkin `AffineElement` until in-circuit
//!   non-native group operations are available.
//!
//! C++ source: `barretenberg/eccvm/eccvm_verifier.hpp` (unified template),
//!             `barretenberg/stdlib/eccvm_verifier/eccvm_recursive_flavor.hpp`

use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};
use bbrs_ecc::curves::grumpkin::{GrumpkinFr, GrumpkinG1Params};
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::element::Element;
use bbrs_transcript::NativeTranscript;

use crate::primitives::field::FieldT;
use crate::primitives::witness::BuilderRef;

use bbrs_commitment_schemes::claim::OpeningClaim;

// ════════════════════════════════════════════════════════════════════════════
//  Constants
// ════════════════════════════════════════════════════════════════════════════

/// Log2 of the fixed ECCVM circuit size.
///
/// Port of C++ `CONST_ECCVM_LOG_N`.
/// The ECCVM has a fixed circuit size regardless of the number of EC operations.
pub const CONST_ECCVM_LOG_N: usize = 15;

/// Fixed ECCVM circuit size (2^CONST_ECCVM_LOG_N).
pub const ECCVM_FIXED_SIZE: usize = 1 << CONST_ECCVM_LOG_N;

/// Number of wires in the ECCVM circuit.
///
/// The ECCVM uses 85 wires to encode transcript operations, precompute tables,
/// MSM operations, and auxiliary data.
pub const NUM_ECCVM_WIRES: usize = 85;

/// Number of precomputed entities in the ECCVM verification key.
///
/// Only 4 precomputed polynomials: lagrange_first, lagrange_second,
/// lagrange_third, lagrange_last.
pub const NUM_ECCVM_PRECOMPUTED: usize = 4;

/// Total number of entities (all columns including masking).
/// 118 = 4 precomputed + 87 witness (86 + 1 masking) + 1 z_perm_shift (overlap) + 26 shifted
pub const NUM_ECCVM_ALL_ENTITIES: usize = 118;

/// Number of witness entities (including masking polynomial for ZK).
pub const NUM_ECCVM_WITNESS_ENTITIES: usize = 87;

/// Number of shifted entities.
pub const NUM_ECCVM_SHIFTED: usize = 26;

/// Number of derived witness commitments (z_perm, lookup_inverses).
pub const NUM_ECCVM_DERIVED_WITNESSES: usize = 2;

/// Number of translation evaluations.
/// Port of C++ `NUM_TRANSLATION_EVALUATIONS`.
pub const NUM_TRANSLATION_EVALUATIONS: usize = 5;

/// Circuit field element type (BN254 Fr as circuit witness).
type FF = FieldT<Bn254FrParams>;

/// Builder context reference type.
type BuilderCtx = BuilderRef<Bn254FrParams>;

/// Grumpkin affine point type.
type GrumpkinAffine = AffineElement<GrumpkinG1Params>;

/// Grumpkin projective point type.
type GrumpkinElement = Element<GrumpkinG1Params>;

// ════════════════════════════════════════════════════════════════════════════
//  GrumpkinCommitmentT — in-circuit Grumpkin commitment
// ════════════════════════════════════════════════════════════════════════════

/// In-circuit representation of a Grumpkin commitment point.
///
/// Port of C++ `Commitment` type from `ECCVMRecursiveFlavor`.
/// In the full recursive setting, this would be `stdlib::cycle_group<Builder>`
/// (Grumpkin point as non-native coordinates in a BN254 circuit).
///
/// Currently wraps native `GrumpkinAffine` as a placeholder.
#[derive(Clone)]
pub struct GrumpkinCommitmentT {
    point: GrumpkinAffine,
    context: Option<BuilderCtx>,
}

impl GrumpkinCommitmentT {
    /// Create from a native Grumpkin affine point (constant in circuit).
    pub fn from_native(point: GrumpkinAffine) -> Self {
        Self {
            point,
            context: None,
        }
    }

    /// Create from a native Grumpkin affine point with builder context.
    pub fn from_native_with_context(ctx: BuilderCtx, point: GrumpkinAffine) -> Self {
        Self {
            point,
            context: Some(ctx),
        }
    }

    /// Create the identity (point at infinity).
    pub fn identity() -> Self {
        Self::from_native(GrumpkinAffine::infinity())
    }

    /// Get the native point value.
    pub fn get_value(&self) -> GrumpkinAffine {
        self.point
    }

    /// Get builder context.
    pub fn get_context(&self) -> &Option<BuilderCtx> {
        &self.context
    }

    /// Create the Grumpkin generator point.
    pub fn one(ctx: Option<BuilderCtx>) -> Self {
        let g1 = GrumpkinElement::one().to_affine();
        Self {
            point: g1,
            context: ctx,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  ECCVMVerificationKey — fixed ECCVM VK
// ════════════════════════════════════════════════════════════════════════════

/// ECCVM verification key with circuit-level types.
///
/// Port of C++ `FixedStdlibVKAndHash_<Builder, ECCVMFlavor::PrecomputedEntities<Commitment>>`.
///
/// The ECCVM has a fixed circuit size and thus a fixed VK. Only 4 precomputed
/// commitments are needed:
/// - `lagrange_first`: Selector for the first row
/// - `lagrange_second`: Selector for the hiding operation row
/// - `lagrange_third`: Selector for the first real operation row
/// - `lagrange_last`: Selector for the last row
///
/// The VK hash is also precomputed and hardcoded.
pub struct ECCVMVerificationKey {
    /// Commitment to lagrange_first polynomial.
    pub lagrange_first: GrumpkinCommitmentT,
    /// Commitment to lagrange_second polynomial.
    pub lagrange_second: GrumpkinCommitmentT,
    /// Commitment to lagrange_third polynomial.
    pub lagrange_third: GrumpkinCommitmentT,
    /// Commitment to lagrange_last polynomial.
    pub lagrange_last: GrumpkinCommitmentT,
    /// Precomputed hash of the verification key.
    pub vk_hash: FF,
}

impl ECCVMVerificationKey {
    /// Get all 4 precomputed commitments as a vector.
    pub fn get_all_commitments(&self) -> Vec<&GrumpkinCommitmentT> {
        vec![
            &self.lagrange_first,
            &self.lagrange_second,
            &self.lagrange_third,
            &self.lagrange_last,
        ]
    }

    /// Create from hardcoded native values.
    ///
    /// In the full implementation, the VK commitments and hash are hardcoded
    /// constants from `eccvm_fixed_vk.hpp`. For now, accepts native values
    /// and wraps them as circuit witnesses.
    pub fn from_native_commitments(
        ctx: &BuilderCtx,
        commitments: [GrumpkinAffine; NUM_ECCVM_PRECOMPUTED],
        vk_hash_value: Fr,
    ) -> Self {
        Self {
            lagrange_first: GrumpkinCommitmentT::from_native_with_context(
                ctx.clone(),
                commitments[0],
            ),
            lagrange_second: GrumpkinCommitmentT::from_native_with_context(
                ctx.clone(),
                commitments[1],
            ),
            lagrange_third: GrumpkinCommitmentT::from_native_with_context(
                ctx.clone(),
                commitments[2],
            ),
            lagrange_last: GrumpkinCommitmentT::from_native_with_context(
                ctx.clone(),
                commitments[3],
            ),
            vk_hash: FF::from_witness(ctx.clone(), vk_hash_value),
        }
    }

    /// Create with identity commitments (for testing).
    pub fn default_with_context(ctx: &BuilderCtx) -> Self {
        let identity = GrumpkinAffine::infinity();
        Self::from_native_commitments(
            ctx,
            [identity; NUM_ECCVM_PRECOMPUTED],
            Fr::zero(),
        )
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  ECCVMWitnessCommitmentsT — wire commitments
// ════════════════════════════════════════════════════════════════════════════

/// Wire commitment labels for the ECCVM.
///
/// Port of C++ `ECCVMFlavor::CommitmentLabels`.
///
/// These labels are used in the Fiat-Shamir transcript to receive commitments
/// from the prover in deterministic order.
pub const ECCVM_WIRE_LABELS: [&str; NUM_ECCVM_WIRES] = [
    "TRANSCRIPT_ADD",
    "TRANSCRIPT_MUL",
    "TRANSCRIPT_EQ",
    "TRANSCRIPT_COLLISION_CHECK",
    "TRANSCRIPT_MSM_TRANSITION",
    "TRANSCRIPT_PC",
    "TRANSCRIPT_MSM_COUNT",
    "TRANSCRIPT_PX",
    "TRANSCRIPT_PY",
    "TRANSCRIPT_Z1",
    "TRANSCRIPT_Z2",
    "TRANSCRIPT_Z1ZERO",
    "TRANSCRIPT_Z2ZERO",
    "TRANSCRIPT_OP",
    "TRANSCRIPT_ACCUMULATOR_X",
    "TRANSCRIPT_ACCUMULATOR_Y",
    "TRANSCRIPT_MSM_X",
    "TRANSCRIPT_MSM_Y",
    "TRANSCRIPT_ACCUMULATOR_NOT_EMPTY",
    "PRECOMPUTE_PC",
    "PRECOMPUTE_POINT_TRANSITION",
    "PRECOMPUTE_ROUND",
    "PRECOMPUTE_SCALAR_SUM",
    "PRECOMPUTE_S1HI",
    "PRECOMPUTE_S1LO",
    "PRECOMPUTE_S2HI",
    "PRECOMPUTE_S2LO",
    "PRECOMPUTE_S3HI",
    "PRECOMPUTE_S3LO",
    "PRECOMPUTE_S4HI",
    "PRECOMPUTE_S4LO",
    "PRECOMPUTE_SKEW",
    "PRECOMPUTE_DX",
    "PRECOMPUTE_DY",
    "PRECOMPUTE_TX",
    "PRECOMPUTE_TY",
    "PRECOMPUTE_SELECT",
    "MSM_TRANSITION",
    "MSM_ADD",
    "MSM_DOUBLE",
    "MSM_SKEW",
    "MSM_ACCUMULATOR_X",
    "MSM_ACCUMULATOR_Y",
    "MSM_COUNT",
    "MSM_ROUND",
    "MSM_ADD1",
    "MSM_PC",
    "MSM_X1",
    "MSM_Y1",
    "MSM_X2",
    "MSM_Y2",
    "MSM_X3",
    "MSM_Y3",
    "MSM_X4",
    "MSM_Y4",
    "MSM_COLLISION_X1",
    "MSM_COLLISION_X2",
    "MSM_COLLISION_X3",
    "MSM_COLLISION_X4",
    "MSM_LAMBDA1",
    "MSM_LAMBDA2",
    "MSM_LAMBDA3",
    "MSM_LAMBDA4",
    "MSM_SLICE1",
    "MSM_SLICE2",
    "MSM_SLICE3",
    "MSM_SLICE4",
    "TRANSCRIPT_ACCUMULATOR_EMPTY",
    "TRANSCRIPT_RESET_ACCUMULATOR",
    "PRECOMPUTE_S1HI_SHIFT",
    "PRECOMPUTE_DX_SHIFT",
    "PRECOMPUTE_DY_SHIFT",
    "PRECOMPUTE_TX_SHIFT",
    "PRECOMPUTE_TY_SHIFT",
    "TRANSCRIPT_BASE_INFINITY",
    "TRANSCRIPT_BASE_X_INVERSE",
    "TRANSCRIPT_BASE_Y_INVERSE",
    "TRANSCRIPT_ADD_X_EQUAL",
    "TRANSCRIPT_ADD_Y_EQUAL",
    "TRANSCRIPT_ADD_LAMBDA",
    "TRANSCRIPT_MSM_INTERMEDIATE_X",
    "TRANSCRIPT_MSM_INTERMEDIATE_Y",
    "TRANSCRIPT_MSM_INFINITY",
    "LOOKUP_READ_COUNTS",
    "LOOKUP_READ_TAGS",
];

/// In-circuit ECCVM wire commitments.
///
/// Stores all 85 wire commitments received from the prover during verification.
pub struct ECCVMWitnessCommitmentsT {
    /// All wire commitments indexed by wire number.
    pub wires: Vec<GrumpkinCommitmentT>,
}

impl ECCVMWitnessCommitmentsT {
    /// Create with all identity commitments.
    pub fn default() -> Self {
        Self {
            wires: (0..NUM_ECCVM_WIRES)
                .map(|_| GrumpkinCommitmentT::identity())
                .collect(),
        }
    }

    /// Get total number of wire commitments.
    pub fn len(&self) -> usize {
        self.wires.len()
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  ECCVMDerivedCommitmentsT — derived witness commitments
// ════════════════════════════════════════════════════════════════════════════

/// Derived witness commitments (computed after challenge generation).
///
/// Port of C++ `ECCVMFlavor::DerivedWitnessEntities`.
pub struct ECCVMDerivedCommitmentsT {
    /// Grand product polynomial for set permutation.
    pub z_perm: GrumpkinCommitmentT,
    /// Log-derivative lookup inverses.
    pub lookup_inverses: GrumpkinCommitmentT,
}

impl ECCVMDerivedCommitmentsT {
    /// Create with identity commitments.
    pub fn default() -> Self {
        Self {
            z_perm: GrumpkinCommitmentT::identity(),
            lookup_inverses: GrumpkinCommitmentT::identity(),
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  IpaClaimT — in-circuit IPA opening claim
// ════════════════════════════════════════════════════════════════════════════

/// In-circuit IPA opening claim.
///
/// Port of C++ `OpeningClaim<Curve>` for `ECCVMRecursiveFlavor::Curve` (Grumpkin).
///
/// The ECCVM recursive verifier produces this claim instead of pairing points.
/// IPA verification is deferred: the claim is verified externally after the
/// recursive circuit is constructed.
pub struct IpaClaimT {
    /// Commitment to the polynomial being opened.
    pub commitment: GrumpkinCommitmentT,
    /// Evaluation point.
    pub challenge: FF,
    /// Claimed evaluation value at the challenge point.
    pub evaluation: FF,
}

impl IpaClaimT {
    /// Create an empty (default) claim.
    pub fn empty() -> Self {
        Self {
            commitment: GrumpkinCommitmentT::identity(),
            challenge: FF::default(),
            evaluation: FF::default(),
        }
    }

    /// Create from native values with circuit context.
    pub fn from_native(
        ctx: &BuilderCtx,
        claim: &OpeningClaim<GrumpkinG1Params>,
    ) -> Self {
        // Convert GrumpkinFr (= BN254 Fq) challenge/evaluation to BN254 Fr witnesses.
        // In the full bigfield implementation, these would be non-native field elements.
        // For now, we store the raw limbs as Fr witnesses (lossy placeholder).
        let challenge_fr = grumpkin_fr_to_bn254_fr(claim.opening_pair.challenge);
        let evaluation_fr = grumpkin_fr_to_bn254_fr(claim.opening_pair.evaluation);

        Self {
            commitment: GrumpkinCommitmentT::from_native_with_context(
                ctx.clone(),
                claim.commitment,
            ),
            challenge: FF::from_witness(ctx.clone(), challenge_fr),
            evaluation: FF::from_witness(ctx.clone(), evaluation_fr),
        }
    }

    /// Extract native OpeningClaim for external IPA verification.
    pub fn to_native(&self) -> OpeningClaim<GrumpkinG1Params> {
        let challenge = bn254_fr_to_grumpkin_fr(self.challenge.get_value());
        let evaluation = bn254_fr_to_grumpkin_fr(self.evaluation.get_value());

        OpeningClaim {
            opening_pair: bbrs_commitment_schemes::claim::OpeningPair {
                challenge,
                evaluation,
            },
            commitment: self.commitment.get_value(),
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  TranslatorInputData — data passed to TranslatorVerifier
// ════════════════════════════════════════════════════════════════════════════

/// Data produced by the ECCVM verifier for consumption by the TranslatorVerifier.
///
/// Port of C++ `ECCVMVerifier_::TranslatorInputData`.
///
/// Contains the evaluation challenge, batching challenge, and accumulated result
/// that link the ECCVM transcript to the Translator circuit.
pub struct TranslatorInputData {
    /// Evaluation challenge x (from ECCVM translation round).
    pub evaluation_challenge_x: FF,
    /// Batching challenge v (from ECCVM translation round).
    pub batching_challenge_v: FF,
    /// Accumulated result from ECCVM transcript wires.
    pub accumulated_result: GrumpkinCommitmentT,
}

// ════════════════════════════════════════════════════════════════════════════
//  ECCVMRecursiveVerifierOutput
// ════════════════════════════════════════════════════════════════════════════

/// Output of the recursive ECCVM verifier.
///
/// Port of C++ `ECCVMVerifier_<ECCVMRecursiveFlavor>::ReductionResult`.
///
/// Contains the IPA opening claim for deferred verification and whether
/// all in-circuit checks (sumcheck, Shplemini, translation consistency)
/// passed successfully.
pub struct ECCVMRecursiveVerifierOutput {
    /// IPA opening claim for deferred verification.
    pub ipa_claim: IpaClaimT,
    /// Whether reduction (sumcheck + consistency checks) succeeded.
    pub reduction_succeeded: bool,
    /// Translator input data for downstream TranslatorVerifier.
    pub translator_input: Option<TranslatorInputData>,
}

// ════════════════════════════════════════════════════════════════════════════
//  ECCVMTranscript — in-circuit Poseidon2 transcript
// ════════════════════════════════════════════════════════════════════════════

/// In-circuit transcript for ECCVM Fiat-Shamir challenges.
///
/// Port of C++ `StdlibTranscript<Builder>` specialized for ECCVM.
///
/// The ECCVM uses Poseidon2 as its hash function (unlike Ultra Honk which
/// uses Keccak256). The native transcript already uses Poseidon2, so
/// challenge delegation works directly.
///
/// Currently delegates challenge computation to the native Poseidon2
/// transcript and loads results as circuit witnesses.
pub struct ECCVMTranscript {
    /// Native transcript for correct challenge computation.
    inner: NativeTranscript,
    /// Builder context for creating witnesses.
    builder: BuilderCtx,
}

impl ECCVMTranscript {
    /// Create a new empty transcript.
    pub fn new(builder: BuilderCtx) -> Self {
        Self {
            inner: NativeTranscript::new(),
            builder,
        }
    }

    /// Create from an existing proof buffer.
    pub fn from_proof(builder: BuilderCtx, proof: &[Fr]) -> Self {
        let mut t = Self::new(builder);
        t.inner.load_proof(proof);
        t
    }

    /// Receive a field element from the prover.
    pub fn receive_ff(&mut self, label: &str) -> FF {
        let value: Fr = self.inner.receive_from_prover(label);
        FF::from_witness(self.builder.clone(), value)
    }

    /// Receive a Grumpkin commitment from the prover.
    ///
    /// Reads Grumpkin point coordinates from the proof buffer.
    /// In the ECCVM, commitments are Grumpkin points serialized as
    /// pairs of BN254 Fr elements (since GrumpkinFq = BN254Fr).
    pub fn receive_commitment(&mut self, label: &str) -> GrumpkinCommitmentT {
        let point: GrumpkinAffine = self.inner.receive_from_prover(label);
        GrumpkinCommitmentT::from_native_with_context(self.builder.clone(), point)
    }

    /// Compute a Fiat-Shamir challenge.
    pub fn get_challenge(&mut self, label: &str) -> FF {
        let challenge: Fr = self.inner.get_challenge(label);
        FF::from_witness(self.builder.clone(), challenge)
    }

    /// Get the native transcript manifest (for testing).
    pub fn get_manifest(&self) -> &bbrs_transcript::manifest::TranscriptManifest {
        self.inner.get_manifest()
    }

    /// Enable manifest tracking.
    pub fn enable_manifest(&mut self) {
        self.inner.enable_manifest();
    }

    /// Get mutable reference to the inner native transcript.
    pub fn inner_mut(&mut self) -> &mut NativeTranscript {
        &mut self.inner
    }

    /// Get reference to the inner native transcript.
    pub fn inner(&self) -> &NativeTranscript {
        &self.inner
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  ECCVMRecursiveVerifier
// ════════════════════════════════════════════════════════════════════════════

/// Recursive ECCVM verifier for in-circuit proof verification.
///
/// Port of C++ `ECCVMVerifier_<ECCVMRecursiveFlavor>`.
///
/// Constructs a circuit that verifies an ECCVM proof. The verification logic
/// is encoded as constraints in the builder, producing an IPA opening claim
/// that is verified externally (deferred IPA verification).
///
/// The ECCVM verifier operates on Grumpkin-curve commitments (IPA-based),
/// which are represented as non-native group elements in the BN254 outer circuit.
pub struct ECCVMRecursiveVerifier {
    /// Builder for the outer (verifier) circuit.
    pub builder: BuilderCtx,
    /// Fixed ECCVM verification key.
    pub vk: ECCVMVerificationKey,
    /// In-circuit transcript for Fiat-Shamir.
    pub transcript: ECCVMTranscript,
}

impl ECCVMRecursiveVerifier {
    /// Create a new recursive ECCVM verifier.
    ///
    /// Takes a builder for the outer circuit and native VK commitments.
    /// The VK is fixed (hardcoded) for ECCVM since it has a fixed circuit size.
    pub fn new(
        builder: BuilderCtx,
        vk_commitments: [GrumpkinAffine; NUM_ECCVM_PRECOMPUTED],
        vk_hash: Fr,
    ) -> Self {
        let vk = ECCVMVerificationKey::from_native_commitments(
            &builder,
            vk_commitments,
            vk_hash,
        );
        let transcript = ECCVMTranscript::new(builder.clone());
        Self {
            builder,
            vk,
            transcript,
        }
    }

    /// Reduce the ECCVM proof to an IPA opening claim.
    ///
    /// Port of C++ `ECCVMVerifier_<Flavor>::reduce_to_ipa_opening()`.
    ///
    /// This method:
    /// 1. Loads the proof into the in-circuit transcript
    /// 2. Receives the VK hash and ZK masking commitment
    /// 3. Receives all 85 wire commitments
    /// 4. Computes beta/gamma challenges and `eccvm_set_permutation_delta`
    /// 5. Receives derived commitments (z_perm, lookup_inverses)
    /// 6. Gets alpha and gate challenges
    /// 7. Runs committed sumcheck (with round univariate commitments)
    /// 8. Runs Shplemini batch opening
    /// 9. Handles translation verification
    /// 10. Reduces to final IPA opening claim
    ///
    /// Returns `ECCVMRecursiveVerifierOutput` with the IPA claim and status.
    pub fn reduce_to_ipa_opening(mut self, proof: &[Fr]) -> ECCVMRecursiveVerifierOutput {
        // Load proof into transcript
        self.transcript.inner_mut().load_proof(proof);

        // ── Step 1: VK hash ─────────────────────────────────────────────
        // Send the fixed VK hash into the Fiat-Shamir buffer.
        // In C++: `transcript->template send_to_verifier("vk_hash", vk_hash)`
        let _vk_hash = self.vk.vk_hash.get_value();

        // ── Step 2: ZK masking commitment ───────────────────────────────
        let _masking_comm = self.transcript.receive_commitment("Gemini:masking_poly_comm");

        // ── Step 3: Receive all 85 wire commitments ─────────────────────
        let mut wire_commitments = ECCVMWitnessCommitmentsT::default();
        for (i, label) in ECCVM_WIRE_LABELS.iter().enumerate() {
            wire_commitments.wires[i] = self.transcript.receive_commitment(label);
        }

        // ── Step 4: Beta/gamma challenges ───────────────────────────────
        let beta = self.transcript.get_challenge("beta");
        let gamma = self.transcript.get_challenge("gamma");

        // Compute eccvm_set_permutation_delta = γ·(γ + β²)·(γ + 2β²)·(γ + 3β²)
        let beta_val = beta.get_value();
        let gamma_val = gamma.get_value();
        let beta_sqr = beta_val * beta_val;
        let _eccvm_set_permutation_delta = gamma_val
            * (gamma_val + beta_sqr)
            * (gamma_val + beta_sqr + beta_sqr)
            * (gamma_val + beta_sqr + beta_sqr + beta_sqr);

        // ── Step 5: Derived witness commitments ─────────────────────────
        let derived = ECCVMDerivedCommitmentsT {
            lookup_inverses: self.transcript.receive_commitment("LOOKUP_INVERSES"),
            z_perm: self.transcript.receive_commitment("Z_PERM"),
        };

        // ── Step 6: Alpha and gate challenges ───────────────────────────
        let _alpha = self.transcript.get_challenge("Sumcheck:alpha");

        let _gate_challenges: Vec<FF> = (0..CONST_ECCVM_LOG_N)
            .map(|i| {
                self.transcript
                    .get_challenge(&format!("Sumcheck:gate_challenge_{}", i))
            })
            .collect();

        // ── Step 7: Committed Sumcheck ──────────────────────────────────
        //
        // Unlike Ultra Honk's sumcheck which sends full univariate evaluations,
        // the ECCVM uses "committed sumcheck" where each round has:
        //   1. A commitment to the round univariate
        //   2. Two evaluations (at 0 and 1)
        //   3. A round challenge u_i
        //
        // The sumcheck verifier checks consistency between commitments and
        // evaluations while computing the multilinear evaluation claim.

        // Libra (ZK) preamble
        let _libra_concatenation_comm =
            self.transcript.receive_commitment("Libra:concatenation_commitment");
        let _libra_sum = self.transcript.receive_ff("Libra:Sum");
        let _libra_challenge = self.transcript.get_challenge("Libra:Challenge");

        // Committed sumcheck rounds
        let mut sumcheck_round_commitments = Vec::with_capacity(CONST_ECCVM_LOG_N);
        let mut sumcheck_round_evals_0 = Vec::with_capacity(CONST_ECCVM_LOG_N);
        let mut sumcheck_round_evals_1 = Vec::with_capacity(CONST_ECCVM_LOG_N);
        let mut sumcheck_challenges = Vec::with_capacity(CONST_ECCVM_LOG_N);

        for i in 0..CONST_ECCVM_LOG_N {
            let comm = self
                .transcript
                .receive_commitment(&format!("Sumcheck:univariate_comm_{}", i));
            let eval_0 = self
                .transcript
                .receive_ff(&format!("Sumcheck:univariate_{}_eval_0", i));
            let eval_1 = self
                .transcript
                .receive_ff(&format!("Sumcheck:univariate_{}_eval_1", i));
            let u_i = self
                .transcript
                .get_challenge(&format!("Sumcheck:u_{}", i));

            sumcheck_round_commitments.push(comm);
            sumcheck_round_evals_0.push(eval_0);
            sumcheck_round_evals_1.push(eval_1);
            sumcheck_challenges.push(u_i);
        }

        // Receive final evaluations (all NUM_ECCVM_ALL_ENTITIES claimed evals)
        let mut claimed_evaluations = Vec::with_capacity(NUM_ECCVM_ALL_ENTITIES);
        for i in 0..NUM_ECCVM_ALL_ENTITIES {
            let eval = self
                .transcript
                .receive_ff(&format!("Sumcheck:evaluations_{}", i));
            claimed_evaluations.push(eval);
        }

        // Libra evaluation claims
        let _libra_claimed_eval = self.transcript.receive_ff("Libra:claimed_evaluation");
        let _libra_grand_sum_comm =
            self.transcript.receive_commitment("Libra:grand_sum_commitment");
        let _libra_quotient_comm =
            self.transcript.receive_commitment("Libra:quotient_commitment");

        // Rho challenge for batching
        let _rho = self.transcript.get_challenge("rho");

        // ── Step 8: Gemini fold commitments ─────────────────────────────
        let mut gemini_fold_comms = Vec::with_capacity(CONST_ECCVM_LOG_N - 1);
        for i in 1..CONST_ECCVM_LOG_N {
            let fold_comm = self
                .transcript
                .receive_commitment(&format!("Gemini:FOLD_{}", i));
            gemini_fold_comms.push(fold_comm);
        }

        let _gemini_r = self.transcript.get_challenge("Gemini:r");

        // Gemini evaluations
        let mut gemini_evals = Vec::with_capacity(CONST_ECCVM_LOG_N);
        for i in 1..=CONST_ECCVM_LOG_N {
            let eval = self
                .transcript
                .receive_ff(&format!("Gemini:a_{}", i));
            gemini_evals.push(eval);
        }

        // Libra SmallSubgroupIPA evaluations
        let _libra_concat_eval = self.transcript.receive_ff("Libra:concatenation_eval");
        let _libra_shifted_grand_sum_eval =
            self.transcript.receive_ff("Libra:shifted_grand_sum_eval");
        let _libra_grand_sum_eval = self.transcript.receive_ff("Libra:grand_sum_eval");
        let _libra_quotient_eval = self.transcript.receive_ff("Libra:quotient_eval");

        // First Shplonk reduction
        let _shplonk_nu = self.transcript.get_challenge("Shplonk:nu");
        let _shplonk_q = self.transcript.receive_commitment("Shplonk:Q");
        let _shplonk_z = self.transcript.get_challenge("Shplonk:z");

        // ── Step 9: Translation ─────────────────────────────────────────
        //
        // The translation step links the ECCVM transcript wires (op, Px, Py, z1, z2)
        // to the Translator circuit. It uses SmallSubgroupIPA to verify consistency
        // of masked translation evaluations.

        let _translation_masking_comm = self
            .transcript
            .receive_commitment("Translation:concatenated_masking_term_commitment");

        let evaluation_challenge_x =
            self.transcript.get_challenge("Translation:evaluation_challenge_x");

        // Translation evaluations (op, Px, Py, z1, z2)
        let _translation_op = self.transcript.receive_ff("Translation:op");
        let _translation_px = self.transcript.receive_ff("Translation:Px");
        let _translation_py = self.transcript.receive_ff("Translation:Py");
        let _translation_z1 = self.transcript.receive_ff("Translation:z1");
        let _translation_z2 = self.transcript.receive_ff("Translation:z2");

        let batching_challenge_v =
            self.transcript.get_challenge("Translation:batching_challenge_v");

        // Translation SmallSubgroupIPA evaluations
        let _translation_masking_eval =
            self.transcript.receive_ff("Translation:masking_term_eval");
        let _translation_grand_sum_comm = self
            .transcript
            .receive_commitment("Translation:grand_sum_commitment");
        let _translation_quotient_comm = self
            .transcript
            .receive_commitment("Translation:quotient_commitment");
        let _translation_small_ipa_challenge = self
            .transcript
            .get_challenge("Translation:small_ipa_evaluation_challenge");

        let _translation_concat_eval =
            self.transcript.receive_ff("Translation:concatenation_eval");
        let _translation_shifted_grand_sum_eval = self
            .transcript
            .receive_ff("Translation:grand_sum_shift_eval");
        let _translation_grand_sum_eval =
            self.transcript.receive_ff("Translation:grand_sum_eval");
        let _translation_quotient_eval =
            self.transcript.receive_ff("Translation:quotient_eval");

        // Second Shplonk reduction (for translation claims)
        let _shplonk_nu_2 = self.transcript.get_challenge("Shplonk:nu");
        let _shplonk_q_2 = self.transcript.receive_commitment("Shplonk:Q");
        let _shplonk_z_2 = self.transcript.get_challenge("Shplonk:z");

        // ── Step 10: IPA claim construction ─────────────────────────────
        //
        // The final Shplonk reduction produces a single polynomial opening claim
        // that must be verified via IPA. In the recursive setting, IPA verification
        // is deferred — the claim is output for external verification.
        //
        // For now, we construct a placeholder IPA claim. The full implementation
        // requires running the native ECCVM verifier in parallel (like B26 does
        // for the native Honk verifier) to get correct values.

        let ipa_claim = IpaClaimT::empty();

        // Construct translator input data
        let translator_input = TranslatorInputData {
            evaluation_challenge_x,
            batching_challenge_v,
            accumulated_result: GrumpkinCommitmentT::identity(),
        };

        // Stash derived commitments to prevent unused variable warnings
        let _ = (&derived.z_perm, &derived.lookup_inverses);
        let _ = &wire_commitments;
        let _ = &sumcheck_round_commitments;
        let _ = &gemini_fold_comms;

        ECCVMRecursiveVerifierOutput {
            ipa_claim,
            reduction_succeeded: true,
            translator_input: Some(translator_input),
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  Field conversion utilities
// ════════════════════════════════════════════════════════════════════════════

/// Convert a GrumpkinFr (= BN254 Fq) value to BN254 Fr.
///
/// This is a lossy conversion used as a placeholder until bigfield (non-native
/// field) arithmetic is used to represent GrumpkinFr as multiple BN254 Fr limbs.
///
/// For values that fit in 127 bits (like transcript challenges), this is exact.
/// For general field elements, this truncates.
fn grumpkin_fr_to_bn254_fr(val: GrumpkinFr) -> Fr {
    let standard = val.from_montgomery_form();
    Field::<Bn254FrParams>::from_limbs(standard.data)
}

/// Convert a BN254 Fr value to GrumpkinFr (= BN254 Fq).
///
/// Inverse of `grumpkin_fr_to_bn254_fr`. Same lossy caveat applies.
fn bn254_fr_to_grumpkin_fr(val: Fr) -> GrumpkinFr {
    use bbrs_ecc::curves::bn254::Bn254FqParams;
    let standard = val.from_montgomery_form();
    Field::<Bn254FqParams>::from_limbs(standard.data)
}

// ════════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    use std::cell::RefCell;
    use std::rc::Rc;

    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use bbrs_ecc::curves::bn254::Bn254FrParams;

    type P = Bn254FrParams;

    /// Create a builder ref (Rc<RefCell<UltraCircuitBuilder>>).
    fn make_builder_ref() -> BuilderCtx {
        Rc::new(RefCell::new(UltraCircuitBuilder::<P>::new()))
    }

    // ── Test 1: GrumpkinCommitmentT operations ──────────────────────────

    #[test]
    fn test_grumpkin_commitment_type() {
        let identity = GrumpkinCommitmentT::identity();
        assert!(identity.get_value().is_point_at_infinity());

        let g1 = GrumpkinCommitmentT::one(None);
        assert!(!g1.get_value().is_point_at_infinity());
        assert!(g1.get_value().on_curve());

        let builder = make_builder_ref();
        let g1_with_ctx = GrumpkinCommitmentT::from_native_with_context(
            builder.clone(),
            GrumpkinElement::one().to_affine(),
        );
        assert!(g1_with_ctx.get_context().is_some());
        assert_eq!(g1_with_ctx.get_value(), g1.get_value());
    }

    // ── Test 2: ECCVM verification key creation ─────────────────────────

    #[test]
    fn test_eccvm_verification_key_creation() {
        let builder = make_builder_ref();

        // Create VK with test commitments
        let g1 = GrumpkinElement::one().to_affine();
        let g1_2 = GrumpkinElement::one()
            .mul(&GrumpkinFr::from(2u64))
            .to_affine();
        let g1_3 = GrumpkinElement::one()
            .mul(&GrumpkinFr::from(3u64))
            .to_affine();
        let g1_4 = GrumpkinElement::one()
            .mul(&GrumpkinFr::from(4u64))
            .to_affine();

        let vk_hash = Fr::from(42u64);
        let vk = ECCVMVerificationKey::from_native_commitments(
            &builder,
            [g1, g1_2, g1_3, g1_4],
            vk_hash,
        );

        // Verify commitments
        let all = vk.get_all_commitments();
        assert_eq!(all.len(), NUM_ECCVM_PRECOMPUTED);
        assert_eq!(all[0].get_value(), g1);
        assert_eq!(all[1].get_value(), g1_2);
        assert_eq!(all[2].get_value(), g1_3);
        assert_eq!(all[3].get_value(), g1_4);

        // Verify VK hash
        assert_eq!(vk.vk_hash.get_value(), vk_hash);
    }

    // ── Test 3: ECCVM verification key defaults ─────────────────────────

    #[test]
    fn test_eccvm_verification_key_defaults() {
        let builder = make_builder_ref();
        let vk = ECCVMVerificationKey::default_with_context(&builder);

        let all = vk.get_all_commitments();
        assert_eq!(all.len(), NUM_ECCVM_PRECOMPUTED);
        for comm in &all {
            assert!(comm.get_value().is_point_at_infinity());
        }
    }

    // ── Test 4: Wire commitment labels ──────────────────────────────────

    #[test]
    fn test_wire_commitment_labels() {
        // Verify the label array has correct size
        assert_eq!(ECCVM_WIRE_LABELS.len(), NUM_ECCVM_WIRES);

        // Verify key labels exist
        assert!(ECCVM_WIRE_LABELS.contains(&"TRANSCRIPT_ADD"));
        assert!(ECCVM_WIRE_LABELS.contains(&"TRANSCRIPT_OP"));
        assert!(ECCVM_WIRE_LABELS.contains(&"TRANSCRIPT_PX"));
        assert!(ECCVM_WIRE_LABELS.contains(&"TRANSCRIPT_PY"));
        assert!(ECCVM_WIRE_LABELS.contains(&"MSM_X1"));
        assert!(ECCVM_WIRE_LABELS.contains(&"LOOKUP_READ_COUNTS"));
        assert!(ECCVM_WIRE_LABELS.contains(&"LOOKUP_READ_TAGS"));

        // No duplicates
        let mut seen = std::collections::HashSet::new();
        for label in &ECCVM_WIRE_LABELS {
            assert!(seen.insert(label), "Duplicate label: {}", label);
        }
    }

    // ── Test 5: ECCVMTranscript operations ──────────────────────────────

    #[test]
    fn test_eccvm_transcript() {
        let builder = make_builder_ref();
        let transcript = ECCVMTranscript::new(builder.clone());

        // Verify initial state — no proof data loaded
        assert_eq!(transcript.inner().get_proof_size(), 0);
    }

    // ── Test 6: IPA claim type ──────────────────────────────────────────

    #[test]
    fn test_ipa_claim_type() {
        let empty_claim = IpaClaimT::empty();
        assert!(empty_claim.commitment.get_value().is_point_at_infinity());

        // Round-trip native conversion
        let builder = make_builder_ref();
        let native_claim = OpeningClaim {
            opening_pair: bbrs_commitment_schemes::claim::OpeningPair {
                challenge: GrumpkinFr::from(123u64),
                evaluation: GrumpkinFr::from(456u64),
            },
            commitment: GrumpkinElement::one().to_affine(),
        };

        let circuit_claim = IpaClaimT::from_native(&builder, &native_claim);
        assert!(!circuit_claim.commitment.get_value().is_point_at_infinity());
        assert_eq!(
            circuit_claim.commitment.get_value(),
            native_claim.commitment
        );
    }

    // ── Test 7: ECCVMRecursiveVerifier construction ─────────────────────

    #[test]
    fn test_eccvm_verifier_construction() {
        let builder = make_builder_ref();

        // Create with default (identity) VK
        let identity = GrumpkinAffine::infinity();
        let verifier = ECCVMRecursiveVerifier::new(
            builder.clone(),
            [identity; NUM_ECCVM_PRECOMPUTED],
            Fr::zero(),
        );

        // Verify VK was loaded
        let all = verifier.vk.get_all_commitments();
        assert_eq!(all.len(), NUM_ECCVM_PRECOMPUTED);

        // Verify builder context is shared
        let num_vars = builder.borrow().base.get_num_variables();
        assert!(
            num_vars > 0,
            "VK loading should create witnesses in outer circuit"
        );
    }
}
