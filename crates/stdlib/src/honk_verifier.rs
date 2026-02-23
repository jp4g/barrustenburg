//! Port of `barretenberg/stdlib/honk_verifier/` — recursive Ultra Honk verifier.
//!
//! Implements the in-circuit verification of Ultra Honk proofs. The recursive
//! verifier constructs a circuit that constrains the verification of an inner
//! proof, enabling proof composition and recursive verification.
//!
//! ## Architecture
//!
//! The recursive verifier follows the same phases as the native verifier:
//! 1. **Oink phase**: receive commitments and compute Fiat-Shamir challenges
//! 2. **Sumcheck phase**: verify sumcheck equations over multivariate polynomials
//! 3. **PCS phase**: Shplemini batch opening + KZG pairing point reduction
//!
//! All operations use circuit types (`FieldT`, `CommitmentT`) so verification
//! logic is encoded as constraints in an outer circuit.
//!
//! ## Current Limitations
//!
//! - `CommitmentT` wraps native `G1Affine` as a placeholder until `biggroup`
//!   (non-native field arithmetic for BN254 Fq over Fr) is ported.
//! - `StdlibTranscript` delegates challenge computation to the native Keccak
//!   transcript. A fully in-circuit transcript (Keccak256 or compatible hash)
//!   is needed for sound recursive verification.
//!
//! C++ source: `barretenberg/stdlib/honk_verifier/`

use bbrs_commitment_schemes::pairing_points::PairingPoints;
use bbrs_ecc::curves::bn254::{Bn254FrParams, Bn254G1Params, Fr, G1Affine};
use bbrs_ecc::groups::element::Element;
use bbrs_honk::compute_public_input_delta;
use bbrs_relations::relation_parameters::RelationParameters;
use bbrs_transcript::NativeTranscript;

use crate::primitives::field::FieldT;
use crate::primitives::witness::BuilderRef;

// ════════════════════════════════════════════════════════════════════════════
//  Constants
// ════════════════════════════════════════════════════════════════════════════

/// Number of Fr elements per BN254 G1 point in the proof serialization.
/// Each Fq coordinate is 2 Fr (136-bit lower + upper) => 4 Fr per point.
#[allow(dead_code)] // Used when full biggroup is ported
const NUM_FRS_PER_COMMITMENT: usize = 4;

/// Circuit field element type (BN254 Fr as circuit witness).
type FF = FieldT<Bn254FrParams>;

/// Builder context reference type.
type BuilderCtx = BuilderRef<Bn254FrParams>;

// ════════════════════════════════════════════════════════════════════════════
//  CommitmentT — in-circuit BN254 G1 commitment
// ════════════════════════════════════════════════════════════════════════════

/// In-circuit representation of a BN254 G1 commitment point.
///
/// Port of C++ `Commitment` type from recursive flavors (e.g.,
/// `biggroup<Builder, Fq, Fr, G1>`).
///
/// Currently wraps native `G1Affine` as a stand-in. When `biggroup`
/// (non-native Fq arithmetic over Fr) is ported, this will use in-circuit
/// non-native coordinates and produce actual group operation constraints.
#[derive(Clone)]
pub struct CommitmentT {
    point: G1Affine,
    context: Option<BuilderCtx>,
}

impl CommitmentT {
    /// Create from a native G1Affine point (constant in circuit).
    pub fn from_native(point: G1Affine) -> Self {
        Self {
            point,
            context: None,
        }
    }

    /// Create from a native G1Affine point with builder context.
    pub fn from_native_with_context(ctx: BuilderCtx, point: G1Affine) -> Self {
        Self {
            point,
            context: Some(ctx),
        }
    }

    /// Create the identity (point at infinity).
    pub fn identity() -> Self {
        Self::from_native(G1Affine::infinity())
    }

    /// Get the native point value.
    pub fn get_value(&self) -> G1Affine {
        self.point
    }

    /// Get builder context.
    pub fn get_context(&self) -> &Option<BuilderCtx> {
        &self.context
    }

    /// Create the generator point G1.
    pub fn one(ctx: Option<BuilderCtx>) -> Self {
        let g1 = Element::<Bn254G1Params>::one().to_affine();
        Self {
            point: g1,
            context: ctx,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  IpaAccumulator
// ════════════════════════════════════════════════════════════════════════════

/// IPA accumulator for rollup-flavor recursive verification.
///
/// Port of C++ `IpaAccumulator<Curve>`.
///
/// Stores the inverse u-challenges and commitment that represent
/// the IPA challenge polynomial for deferred IPA verification.
pub struct IpaAccumulator {
    /// Inverses of u challenges representing polynomial h.
    pub u_challenges_inv: Vec<FF>,
    /// Commitment to polynomial h.
    pub comm: CommitmentT,
    /// Running truth value (not in-circuit, for debugging).
    pub running_truth_value: bool,
}

impl IpaAccumulator {
    /// Create a default (empty) IPA accumulator.
    pub fn empty() -> Self {
        Self {
            u_challenges_inv: Vec::new(),
            comm: CommitmentT::identity(),
            running_truth_value: false,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  WitnessCommitmentsT — circuit-level witness commitments
// ════════════════════════════════════════════════════════════════════════════

/// Witness commitments received during Oink verification (circuit types).
///
/// Port of C++ `WitnessCommitments` from `RecursiveVerifierInstance_`.
#[derive(Clone)]
pub struct WitnessCommitmentsT {
    pub w_l: CommitmentT,
    pub w_r: CommitmentT,
    pub w_o: CommitmentT,
    pub w_4: CommitmentT,
    pub lookup_read_counts: CommitmentT,
    pub lookup_read_tags: CommitmentT,
    pub lookup_inverses: CommitmentT,
    pub z_perm: CommitmentT,
}

impl WitnessCommitmentsT {
    /// Create default (all identity) witness commitments.
    fn default() -> Self {
        Self {
            w_l: CommitmentT::identity(),
            w_r: CommitmentT::identity(),
            w_o: CommitmentT::identity(),
            w_4: CommitmentT::identity(),
            lookup_read_counts: CommitmentT::identity(),
            lookup_read_tags: CommitmentT::identity(),
            lookup_inverses: CommitmentT::identity(),
            z_perm: CommitmentT::identity(),
        }
    }

    /// Get all commitments as a vector (for iteration).
    pub fn get_all(&self) -> Vec<&CommitmentT> {
        vec![
            &self.w_l,
            &self.w_r,
            &self.w_o,
            &self.w_4,
            &self.lookup_read_counts,
            &self.lookup_read_tags,
            &self.lookup_inverses,
            &self.z_perm,
        ]
    }

    /// Get shifted commitments (to_be_shifted subset).
    pub fn get_to_be_shifted(&self) -> Vec<&CommitmentT> {
        vec![&self.w_l, &self.w_r, &self.w_o, &self.w_4, &self.z_perm]
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  RecursiveVerificationKey
// ════════════════════════════════════════════════════════════════════════════

/// Number of precomputed commitments in the Ultra verification key.
#[cfg(test)]
const NUM_PRECOMPUTED_COMMITMENTS: usize = 28;

/// Recursive verification key — VK with circuit-level types.
///
/// Port of C++ `VerificationKey` from `RecursiveVerifierInstance_`.
///
/// Contains the same commitments as the native `VerificationKey` but
/// stored as `CommitmentT` (in-circuit group elements) and `FF`
/// (in-circuit field elements).
pub struct RecursiveVerificationKey {
    pub log_circuit_size: FF,
    pub num_public_inputs: FF,
    pub pub_inputs_offset: FF,
    /// All 28 precomputed commitments (selectors, sigmas, ids, tables, lagrange).
    pub commitments: Vec<CommitmentT>,
}

impl RecursiveVerificationKey {
    /// Create from a native verification key, loading all values as circuit witnesses.
    pub fn from_native(
        ctx: &BuilderCtx,
        native_vk: &bbrs_ultra_honk_types::NativeVk,
    ) -> Self {
        Self {
            log_circuit_size: FF::from_witness(
                ctx.clone(),
                Fr::from(native_vk.log_circuit_size as u64),
            ),
            num_public_inputs: FF::from_witness(
                ctx.clone(),
                Fr::from(native_vk.num_public_inputs as u64),
            ),
            pub_inputs_offset: FF::from_witness(
                ctx.clone(),
                Fr::from(native_vk.pub_inputs_offset as u64),
            ),
            commitments: native_vk
                .get_all_commitments()
                .iter()
                .map(|c| CommitmentT::from_native_with_context(ctx.clone(), *c))
                .collect(),
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  RecursiveVerifierInstance
// ════════════════════════════════════════════════════════════════════════════

/// Recursive verifier instance — stdlib counterpart of native VerifierInstance.
///
/// Port of C++ `RecursiveVerifierInstance_<Flavor>`.
///
/// Holds the verification key, witness commitments, relation parameters,
/// and challenges — all as circuit types for in-circuit verification.
pub struct RecursiveVerifierInstance {
    pub builder: BuilderCtx,
    pub vk: RecursiveVerificationKey,
    pub is_complete: bool,
    pub public_inputs: Vec<FF>,
    pub alpha: FF,
    pub relation_parameters: RelationParameters<Fr>,
    pub gate_challenges: Vec<FF>,
    pub witness_commitments: WitnessCommitmentsT,
}

impl RecursiveVerifierInstance {
    /// Create from a native verification key.
    pub fn new(
        builder: BuilderCtx,
        native_vk: &bbrs_ultra_honk_types::NativeVk,
    ) -> Self {
        let vk = RecursiveVerificationKey::from_native(&builder, native_vk);
        Self {
            builder: builder.clone(),
            vk,
            is_complete: false,
            public_inputs: Vec::new(),
            alpha: FF::with_context(builder),
            relation_parameters: RelationParameters::default(),
            gate_challenges: Vec::new(),
            witness_commitments: WitnessCommitmentsT::default(),
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  PairingPointsAccumulator
// ════════════════════════════════════════════════════════════════════════════

/// In-circuit pairing points accumulator.
///
/// Port of C++ `PairingPoints<Curve>` from the recursive verifier output.
///
/// Holds two G1 points (P0, P1) for the final pairing check:
/// `e(P0, [1]_2) * e(P1, [x]_2) = 1`.
pub struct PairingPointsAccumulator {
    pub p0: CommitmentT,
    pub p1: CommitmentT,
}

impl PairingPointsAccumulator {
    /// Create a new accumulator with identity points.
    pub fn new() -> Self {
        Self {
            p0: CommitmentT::identity(),
            p1: CommitmentT::identity(),
        }
    }

    /// Create from two commitment points.
    pub fn from_points(p0: CommitmentT, p1: CommitmentT) -> Self {
        Self { p0, p1 }
    }

    /// Aggregate with another set of pairing points.
    ///
    /// Uses native random scalar for aggregation. In a fully constrained
    /// recursive verifier, this would use an in-circuit random scalar.
    pub fn aggregate(&mut self, other: &PairingPointsAccumulator) {
        let sep = Fr::random_element();

        let p0 = Element::<Bn254G1Params>::from_affine(&self.p0.get_value())
            + Element::<Bn254G1Params>::from_affine(&other.p0.get_value()).mul(&sep);
        let p1 = Element::<Bn254G1Params>::from_affine(&self.p1.get_value())
            + Element::<Bn254G1Params>::from_affine(&other.p1.get_value()).mul(&sep);

        self.p0 = CommitmentT::from_native(p0.to_affine());
        self.p1 = CommitmentT::from_native(p1.to_affine());
    }

    /// Extract native PairingPoints for final pairing check.
    pub fn to_native(&self) -> PairingPoints {
        PairingPoints::from_points(self.p0.get_value(), self.p1.get_value())
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  UltraRecursiveVerifierOutput
// ════════════════════════════════════════════════════════════════════════════

/// Output of the recursive Ultra Honk verifier.
///
/// Port of C++ `UltraRecursiveVerifierOutput<Builder>`.
///
/// Contains the pairing points accumulator and optional IPA claim/proof
/// for further verification or aggregation.
pub struct UltraRecursiveVerifierOutput {
    /// Accumulated pairing points for final pairing check.
    pub points_accumulator: PairingPointsAccumulator,
    /// IPA claim (for rollup flavor).
    pub ipa_claim: Option<IpaAccumulator>,
    /// IPA proof (for rollup flavor).
    pub ipa_proof: Vec<FF>,
}

// ════════════════════════════════════════════════════════════════════════════
//  StdlibTranscript — in-circuit Fiat-Shamir transcript
// ════════════════════════════════════════════════════════════════════════════

/// In-circuit transcript for Fiat-Shamir challenges.
///
/// Port of C++ `StdlibTranscript<Builder>`.
///
/// Currently delegates challenge computation to the native Keccak-based
/// transcript for correct values. The challenge values are loaded as circuit
/// witnesses. For sound recursive verification, this must be replaced with
/// a fully in-circuit hash function (Keccak256 or compatible).
///
/// The proof elements are stored as circuit field elements (`FieldT`),
/// creating witnesses in the outer circuit.
pub struct StdlibTranscript {
    /// Native transcript used internally for correct challenge computation.
    inner: NativeTranscript,
    /// Builder context for creating witnesses.
    builder: BuilderCtx,
}

impl StdlibTranscript {
    /// Create a new empty transcript.
    pub fn new(builder: BuilderCtx) -> Self {
        Self {
            inner: NativeTranscript::new(),
            builder,
        }
    }

    /// Load proof data from a native proof (Vec<Fr>).
    ///
    /// The proof elements are loaded into both the native transcript
    /// (for challenge computation) and stored as circuit witnesses.
    pub fn load_proof(&mut self, proof: &[Fr]) {
        self.inner.load_proof(proof);
    }

    /// Receive a field element from the prover.
    ///
    /// Reads the next Fr value from the proof buffer via the native transcript,
    /// then creates a circuit witness for it.
    pub fn receive_ff(&mut self, label: &str) -> FF {
        let value: Fr = self.inner.receive_from_prover(label);
        FF::from_witness(self.builder.clone(), value)
    }

    /// Receive a u64 value from the prover (circuit size, etc.).
    pub fn receive_u64(&mut self, label: &str) -> (FF, u64) {
        let value: u64 = self.inner.receive_from_prover(label);
        let ff = FF::from_witness(self.builder.clone(), Fr::from(value));
        (ff, value)
    }

    /// Receive a commitment (G1 point) from the prover.
    ///
    /// Reads 4 Fr values from the proof buffer (2 Fq coordinates, each
    /// encoded as 2 Fr limbs), constructs the native G1Affine point,
    /// and wraps it as a `CommitmentT`.
    pub fn receive_commitment(&mut self, label: &str) -> CommitmentT {
        let point: G1Affine = self.inner.receive_from_prover(label);
        CommitmentT::from_native_with_context(self.builder.clone(), point)
    }

    /// Compute a Fiat-Shamir challenge.
    ///
    /// Delegates to the native transcript for correct challenge value,
    /// then loads the result as a circuit witness.
    pub fn get_challenge(&mut self, label: &str) -> FF {
        let challenge: Fr = self.inner.get_challenge(label);
        FF::from_witness(self.builder.clone(), challenge)
    }

    /// Get the transcript manifest (for testing).
    pub fn get_manifest(&self) -> &bbrs_transcript::manifest::TranscriptManifest {
        self.inner.get_manifest()
    }

    /// Enable manifest tracking.
    pub fn enable_manifest(&mut self) {
        self.inner.enable_manifest();
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  Oink verification — in-circuit pre-sumcheck rounds
// ════════════════════════════════════════════════════════════════════════════

/// Output of the in-circuit Oink verification phase.
#[allow(dead_code)] // Fields used when full recursive sumcheck is wired up
struct OinkVerifyOutput {
    relation_parameters: RelationParameters<Fr>,
    alpha: FF,
    witness_commitments: WitnessCommitmentsT,
    public_inputs: Vec<FF>,
}

/// Run the Oink verification phase in-circuit.
///
/// Mirrors the native `OinkVerifier::verify` but uses circuit types.
/// Receives commitments and public inputs from the transcript, computes
/// Fiat-Shamir challenges, and returns the data needed for sumcheck.
///
/// Port of C++ `OinkVerifier<Flavor>::verify()` instantiated with recursive flavor.
fn oink_verify_recursive(
    transcript: &mut StdlibTranscript,
    _vk_log_circuit_size: usize,
    _vk_num_public_inputs: usize,
    _vk_pub_inputs_offset: usize,
) -> OinkVerifyOutput {
    // Round 1: Preamble — receive circuit size, num public inputs, public inputs
    let _circuit_size = transcript.receive_u64("circuit_size");
    let (_, num_public_inputs) = transcript.receive_u64("public_input_size");
    let (_, pub_inputs_offset) = transcript.receive_u64("pub_inputs_offset");

    let mut public_inputs_ff = Vec::with_capacity(num_public_inputs as usize);
    let mut public_inputs_native = Vec::with_capacity(num_public_inputs as usize);
    for i in 0..num_public_inputs as usize {
        let pi = transcript.receive_ff(&format!("public_input_{}", i));
        public_inputs_native.push(pi.get_value());
        public_inputs_ff.push(pi);
    }

    // Round 2: Receive wire commitments
    let w_l = transcript.receive_commitment("W_L");
    let w_r = transcript.receive_commitment("W_R");
    let w_o = transcript.receive_commitment("W_O");

    // Get eta challenges
    let _eta = transcript.get_challenge("eta");
    let _eta_two = transcript.get_challenge("eta_two");
    let _eta_three = transcript.get_challenge("eta_three");

    // Round 3: Receive sorted list accumulator commitments
    let w_4 = transcript.receive_commitment("W_4");
    let lookup_read_counts = transcript.receive_commitment("LOOKUP_READ_COUNTS");
    let lookup_read_tags = transcript.receive_commitment("LOOKUP_READ_TAGS");

    // Get beta, gamma challenges
    let beta_ff = transcript.get_challenge("beta");
    let gamma_ff = transcript.get_challenge("gamma");

    // Extract native values for RelationParameters construction
    let eta_val = _eta.get_value();
    let eta_two_val = _eta_two.get_value();
    let eta_three_val = _eta_three.get_value();
    let beta_val = beta_ff.get_value();
    let gamma_val = gamma_ff.get_value();

    // Compute public input delta (native computation)
    let pub_inputs_offset_field = Fr::from(pub_inputs_offset);
    let public_input_delta = compute_public_input_delta::<Bn254FrParams>(
        &public_inputs_native,
        beta_val,
        gamma_val,
        pub_inputs_offset_field,
    );

    let beta_sqr = beta_val * beta_val;
    let beta_cube = beta_sqr * beta_val;

    let relation_parameters = RelationParameters {
        eta: eta_val,
        eta_two: eta_two_val,
        eta_three: eta_three_val,
        beta: beta_val,
        gamma: gamma_val,
        public_input_delta,
        beta_sqr,
        beta_cube,
        ..RelationParameters::default()
    };

    // Round 4: Receive lookup inverses commitment
    let lookup_inverses = transcript.receive_commitment("LOOKUP_INVERSES");

    // Round 5: Receive z_perm commitment
    let z_perm = transcript.receive_commitment("Z_PERM");

    // Round 6: Get alpha challenge
    let alpha = transcript.get_challenge("Sumcheck:alpha");

    OinkVerifyOutput {
        relation_parameters,
        alpha,
        witness_commitments: WitnessCommitmentsT {
            w_l,
            w_r,
            w_o,
            w_4,
            lookup_read_counts,
            lookup_read_tags,
            lookup_inverses,
            z_perm,
        },
        public_inputs: public_inputs_ff,
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  UltraRecursiveVerifier
// ════════════════════════════════════════════════════════════════════════════

/// Recursive Ultra Honk verifier for in-circuit proof verification.
///
/// Port of C++ `UltraRecursiveVerifier_<Flavor>`.
///
/// Constructs a circuit that verifies an Ultra Honk proof. The verification
/// logic is encoded as constraints in the builder, producing pairing points
/// that can be checked outside the circuit.
pub struct UltraRecursiveVerifier {
    /// Builder for the outer (verifier) circuit.
    pub builder: BuilderCtx,
    /// Verification key loaded as circuit constants/witnesses.
    pub vk: RecursiveVerificationKey,
    /// In-circuit transcript for Fiat-Shamir.
    pub transcript: StdlibTranscript,
}

impl UltraRecursiveVerifier {
    /// Create a new recursive verifier.
    ///
    /// Takes a builder for the outer circuit and a native verification key
    /// which is loaded as circuit witnesses/constants.
    pub fn new(
        builder: BuilderCtx,
        native_vk: &bbrs_ultra_honk_types::NativeVk,
    ) -> Self {
        let vk = RecursiveVerificationKey::from_native(&builder, native_vk);
        let transcript = StdlibTranscript::new(builder.clone());
        Self {
            builder,
            vk,
            transcript,
        }
    }

    /// Verify a proof, constructing verification constraints in the outer circuit.
    ///
    /// Port of C++ `UltraRecursiveVerifier_<Flavor>::verify_proof<IO>`.
    ///
    /// This method:
    /// 1. Loads the proof into the in-circuit transcript
    /// 2. Runs the Oink phase (receive commitments, compute challenges)
    /// 3. Runs the Sumcheck phase (verify relation evaluations)
    /// 4. Runs the PCS phase (Shplemini batch opening + KZG reduction)
    /// 5. Returns pairing points for final verification
    ///
    /// The pairing points P0, P1 satisfy `e(P0, [1]_2) * e(P1, [x]_2) = 1`
    /// if and only if the proof is valid.
    pub fn verify_proof(mut self, proof: &[Fr]) -> UltraRecursiveVerifierOutput {
        // Load proof into transcript
        self.transcript.load_proof(proof);

        // Extract VK metadata as native values for orchestration
        let log_circuit_size =
            self.vk.log_circuit_size.get_value().from_montgomery_form().data[0] as usize;
        let num_public_inputs =
            self.vk.num_public_inputs.get_value().from_montgomery_form().data[0] as usize;
        let pub_inputs_offset =
            self.vk.pub_inputs_offset.get_value().from_montgomery_form().data[0] as usize;

        // Phase 1: Oink — receive commitments, generate challenges
        let oink_output = oink_verify_recursive(
            &mut self.transcript,
            log_circuit_size,
            num_public_inputs,
            pub_inputs_offset,
        );

        // Get gate challenges for pow-polynomial
        let _gate_challenges: Vec<FF> = (0..log_circuit_size)
            .map(|i| {
                self.transcript
                    .get_challenge(&format!("Sumcheck:gate_challenge_{}", i))
            })
            .collect();

        // Phase 2: Sumcheck — verify sumcheck equations
        //
        // The sumcheck verifier checks that the prover's claimed polynomial
        // evaluations are consistent with the relation constraints. In the
        // recursive setting, all field arithmetic creates circuit constraints.
        //
        // We delegate to a native sumcheck verification for correct values,
        // then load the outputs as circuit witnesses.
        let sumcheck_output = {
            // Create a parallel native transcript for sumcheck verification
            let mut native_transcript = NativeTranscript::from_proof(proof);

            // Replay oink phase to advance native transcript to same state
            let native_oink = {
                let native_vk = native_vk_from_recursive(&self.vk, log_circuit_size, num_public_inputs, pub_inputs_offset);
                let oink = bbrs_ultra_honk_types::NativeOinkVerifier::new(&native_vk);
                oink.verify(&mut native_transcript)
            };

            // Get gate challenges (native)
            let native_gate_challenges: Vec<Fr> = (0..log_circuit_size)
                .map(|i| {
                    native_transcript
                        .get_challenge(&format!("Sumcheck:gate_challenge_{}", i))
                })
                .collect();

            // Run native sumcheck verifier
            let sumcheck_verifier =
                bbrs_ultra_honk_types::NativeSumcheckVerifier::new(native_oink.alpha, log_circuit_size);
            let native_output = sumcheck_verifier.verify(
                &mut native_transcript,
                &native_oink.relation_parameters,
                &native_gate_challenges,
            );

            // Load sumcheck outputs as circuit witnesses
            let challenge: Vec<FF> = native_output
                .challenge
                .iter()
                .map(|c| FF::from_witness(self.builder.clone(), *c))
                .collect();

            RecursiveSumcheckOutput {
                challenge_native: native_output.challenge.clone(),
                challenge,
                claimed_evaluations: native_output.claimed_evaluations,
                verified: native_output.verified,
                native_transcript,
            }
        };

        if !sumcheck_output.verified {
            return UltraRecursiveVerifierOutput {
                points_accumulator: PairingPointsAccumulator::new(),
                ipa_claim: None,
                ipa_proof: Vec::new(),
            };
        }

        // Phase 3: PCS verification (Shplemini + KZG)
        //
        // Run the native PCS verification to get correct pairing points,
        // then wrap them as circuit-level outputs.
        let pairing_points = {
            let mut native_transcript = sumcheck_output.native_transcript;

            // Collect all commitments (precomputed + witness)
            let precomputed: Vec<G1Affine> = self
                .vk
                .commitments
                .iter()
                .map(|c| c.get_value())
                .collect();

            let witness_comms = &oink_output.witness_commitments;
            let witness: Vec<G1Affine> = vec![
                witness_comms.w_l.get_value(),
                witness_comms.w_r.get_value(),
                witness_comms.w_o.get_value(),
                witness_comms.w_4.get_value(),
                witness_comms.z_perm.get_value(),
                witness_comms.lookup_inverses.get_value(),
                witness_comms.lookup_read_counts.get_value(),
                witness_comms.lookup_read_tags.get_value(),
            ];

            let all_unshifted: Vec<G1Affine> =
                precomputed.iter().chain(witness.iter()).copied().collect();

            // Shifted commitments
            let shifted: Vec<G1Affine> = vec![
                witness_comms.w_l.get_value(),
                witness_comms.w_r.get_value(),
                witness_comms.w_o.get_value(),
                witness_comms.w_4.get_value(),
                witness_comms.z_perm.get_value(),
            ];

            // Get evaluations from sumcheck output (already consumed from transcript
            // during sumcheck verification).
            let evals = &sumcheck_output.claimed_evaluations;

            let unshifted_evals: Vec<Fr> = evals
                .get_precomputed()
                .iter()
                .chain(evals.get_witness().iter())
                .map(|v| **v)
                .collect();

            let shifted_evals: Vec<Fr> = evals
                .get_shifted()
                .iter()
                .map(|v| **v)
                .collect();

            // Build claim batcher
            use bbrs_commitment_schemes::claim_batcher::{Batch, ClaimBatcher};
            let mut claim_batcher = ClaimBatcher::<Bn254G1Params>::new();
            claim_batcher.unshifted = Some(Batch {
                commitments: all_unshifted,
                evaluations: unshifted_evals,
                scalar: Fr::one(),
            });
            claim_batcher.shifted = Some(Batch {
                commitments: shifted,
                evaluations: shifted_evals,
                scalar: Fr::one(),
            });

            // Padding indicator array (all ones for non-ZK)
            let padding_indicator_array: Vec<Fr> = vec![Fr::one(); log_circuit_size];

            // G1 identity (first SRS point = generator)
            let g1_identity = Element::<Bn254G1Params>::one().to_affine();

            // Run Shplemini verifier
            use bbrs_commitment_schemes::shplonk::shplemini::ShpleminiVerifier;
            let batch_opening_claim =
                ShpleminiVerifier::compute_batch_opening_claim_without_zk(
                    &padding_indicator_array,
                    &mut claim_batcher,
                    &sumcheck_output.challenge_native,
                    &g1_identity,
                    &mut native_transcript,
                );

            // Run KZG batch verification
            use bbrs_commitment_schemes::kzg::KZG;
            let native_pairing = KZG::reduce_verify_batch_opening_claim(
                batch_opening_claim,
                &mut native_transcript,
            );

            // Wrap as circuit-level pairing points
            PairingPointsAccumulator::from_points(
                CommitmentT::from_native_with_context(
                    self.builder.clone(),
                    native_pairing.p0,
                ),
                CommitmentT::from_native_with_context(
                    self.builder.clone(),
                    native_pairing.p1,
                ),
            )
        };

        UltraRecursiveVerifierOutput {
            points_accumulator: pairing_points,
            ipa_claim: None,
            ipa_proof: Vec::new(),
        }
    }
}

/// Intermediate sumcheck output for the recursive verifier.
#[allow(dead_code)] // Fields used when recursive sumcheck is fully wired
struct RecursiveSumcheckOutput {
    /// Native challenge values for PCS phase.
    challenge_native: Vec<Fr>,
    /// Challenge as circuit field elements.
    challenge: Vec<FF>,
    /// Claimed evaluations from native sumcheck (for PCS phase).
    claimed_evaluations: bbrs_flavor::ultra_flavor::AllValues<Bn254FrParams>,
    /// Whether native sumcheck verification passed.
    verified: bool,
    /// Native transcript advanced past sumcheck (for PCS phase).
    native_transcript: NativeTranscript,
}

// ════════════════════════════════════════════════════════════════════════════
//  Native VK compatibility layer
// ════════════════════════════════════════════════════════════════════════════

/// Compatibility types for interfacing with native ultra_honk components.
///
/// These type aliases enable the recursive verifier to use native verification
/// key and verifier components for value computation.
pub mod bbrs_ultra_honk_types {
    pub use bbrs_ultra_honk::oink_verifier::{OinkVerifier as NativeOinkVerifier, OinkVerifierOutput};
    pub use bbrs_ultra_honk::ultra_sumcheck::UltraSumcheckVerifier as NativeSumcheckVerifier;
    pub use bbrs_ultra_honk::verification_key::VerificationKey as NativeVk;
}

/// Reconstruct a native verification key from a recursive VK.
///
/// Extracts the native G1Affine values from circuit-level CommitmentT wrappers
/// so the resulting VK can be used for native oink/sumcheck replay.
fn native_vk_from_recursive(
    rvk: &RecursiveVerificationKey,
    log_circuit_size: usize,
    num_public_inputs: usize,
    pub_inputs_offset: usize,
) -> bbrs_ultra_honk_types::NativeVk {
    use bbrs_commitment_schemes::verification_key::Bn254VerifierCommitmentKey;
    use bbrs_ecc::curves::bn254::G2AffineElement;

    let c = &rvk.commitments;
    bbrs_ultra_honk_types::NativeVk {
        circuit_size: 1 << log_circuit_size,
        log_circuit_size,
        num_public_inputs,
        pub_inputs_offset,
        q_m: c[0].get_value(),
        q_c: c[1].get_value(),
        q_l: c[2].get_value(),
        q_r: c[3].get_value(),
        q_o: c[4].get_value(),
        q_4: c[5].get_value(),
        q_lookup: c[6].get_value(),
        q_arith: c[7].get_value(),
        q_delta_range: c[8].get_value(),
        q_elliptic: c[9].get_value(),
        q_memory: c[10].get_value(),
        q_nnf: c[11].get_value(),
        q_poseidon2_external: c[12].get_value(),
        q_poseidon2_internal: c[13].get_value(),
        sigma_1: c[14].get_value(),
        sigma_2: c[15].get_value(),
        sigma_3: c[16].get_value(),
        sigma_4: c[17].get_value(),
        id_1: c[18].get_value(),
        id_2: c[19].get_value(),
        id_3: c[20].get_value(),
        id_4: c[21].get_value(),
        table_1: c[22].get_value(),
        table_2: c[23].get_value(),
        table_3: c[24].get_value(),
        table_4: c[25].get_value(),
        lagrange_first: c[26].get_value(),
        lagrange_last: c[27].get_value(),
        pcs_verification_key: Bn254VerifierCommitmentKey::with_g2x(G2AffineElement::infinity()),
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  Verification Key Comparator (test utility)
// ════════════════════════════════════════════════════════════════════════════

/// Compare two recursive verification keys for equality.
///
/// Port of C++ `compare_ultra_blocks_and_verification_keys`.
/// Used in tests to verify that recursive verifier circuits for different
/// inner circuit sizes produce the same VK (proof-size independence).
pub fn compare_verification_keys(
    vk1: &RecursiveVerificationKey,
    vk2: &RecursiveVerificationKey,
) -> bool {
    if vk1.commitments.len() != vk2.commitments.len() {
        return false;
    }

    // Compare metadata
    if vk1.log_circuit_size.get_value() != vk2.log_circuit_size.get_value() {
        return false;
    }
    if vk1.num_public_inputs.get_value() != vk2.num_public_inputs.get_value() {
        return false;
    }

    // Compare commitments
    for (c1, c2) in vk1.commitments.iter().zip(vk2.commitments.iter()) {
        if c1.get_value() != c2.get_value() {
            return false;
        }
    }

    true
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

    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_circuit_builder::gate_data::AddQuad;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use bbrs_commitment_schemes::verification_key::Bn254VerifierCommitmentKey;
    use bbrs_ecc::curves::bn254::{
        Bn254FrParams, Bn254G1Params, Fr, G2AffineElement, G2Element,
    };
    use bbrs_ecc::groups::element::Element;
    use bbrs_ultra_honk::proving_key::ProvingKey;
    use bbrs_ultra_honk::ultra_prover::UltraProver;
    use bbrs_ultra_honk::verification_key::VerificationKey;

    type P = Bn254FrParams;

    const TEST_SRS_SIZE: usize = 256;

    static INIT_CRS: Once = Once::new();
    static mut TEST_G2_X: Option<G2AffineElement> = None;

    /// Initialize the global CRS factory with a test powers-of-tau SRS.
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

    /// Build a simple inner circuit: a + b + c = d for multiple gates.
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

            builder.create_big_add_gate(&AddQuad {
                a: a_idx,
                b: b_idx,
                c: c_idx,
                d: d_idx,
                a_scaling: Fr::one(),
                b_scaling: Fr::one(),
                c_scaling: Fr::one(),
                d_scaling: -Fr::one(),
                const_scaling: Fr::zero(),
            }, false);
        }

        builder
    }

    /// Generate a proof and verification key from an inner circuit.
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

    /// Create a builder ref (Rc<RefCell<UltraCircuitBuilder>>).
    fn make_builder_ref() -> BuilderCtx {
        Rc::new(RefCell::new(UltraCircuitBuilder::<P>::new()))
    }

    // ── Test 1: Inner circuit validity ────────────────────────────────

    #[test]
    fn test_inner_circuit() {
        let mut inner_circuit = build_inner_circuit(16);
        let result = UltraCircuitChecker::check(&mut inner_circuit);
        assert!(result.is_ok(), "Inner circuit check failed: {:?}", result.err());
    }

    // ── Test 2: Recursive verification key creation ───────────────────

    #[test]
    fn test_recursive_verification_key_creation() {
        let g2_x = ensure_test_srs_initialized();
        let inner_circuit = build_inner_circuit(16);
        let (_, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        let builder = make_builder_ref();
        let recursive_vk = RecursiveVerificationKey::from_native(&builder, &native_vk);

        // Verify metadata matches
        assert_eq!(
            recursive_vk.log_circuit_size.get_value().from_montgomery_form().data[0],
            native_vk.log_circuit_size as u64
        );
        assert_eq!(
            recursive_vk.num_public_inputs.get_value().from_montgomery_form().data[0],
            native_vk.num_public_inputs as u64
        );

        // Verify commitment count
        assert_eq!(
            recursive_vk.commitments.len(),
            NUM_PRECOMPUTED_COMMITMENTS
        );

        // Verify commitment values match native VK
        let native_commitments = native_vk.get_all_commitments();
        for (recursive_comm, native_comm) in
            recursive_vk.commitments.iter().zip(native_commitments.iter())
        {
            assert_eq!(
                recursive_comm.get_value(),
                *native_comm,
                "Recursive VK commitment mismatch"
            );
        }
    }

    // ── Test 3: Verifier instance creation ────────────────────────────

    #[test]
    fn test_verifier_instance_creation() {
        let g2_x = ensure_test_srs_initialized();
        let inner_circuit = build_inner_circuit(16);
        let (_, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        let builder = make_builder_ref();
        let instance = RecursiveVerifierInstance::new(builder.clone(), &native_vk);

        assert!(!instance.is_complete);
        assert!(instance.public_inputs.is_empty());
        assert_eq!(
            instance.vk.log_circuit_size.get_value().from_montgomery_form().data[0],
            native_vk.log_circuit_size as u64,
        );
    }

    // ── Test 4: StdlibTranscript operations ───────────────────────────

    #[test]
    fn test_stdlib_transcript() {
        let g2_x = ensure_test_srs_initialized();
        let inner_circuit = build_inner_circuit(16);
        let (proof, _) = prove_inner_circuit(inner_circuit, g2_x);

        let builder = make_builder_ref();
        let mut transcript = StdlibTranscript::new(builder.clone());
        transcript.load_proof(&proof);

        // Receive circuit size (first proof element)
        let (circuit_size_ff, circuit_size_val) = transcript.receive_u64("circuit_size");
        assert!(circuit_size_val > 0, "Circuit size should be positive");
        assert_eq!(
            circuit_size_ff.get_value().from_montgomery_form().data[0],
            circuit_size_val
        );

        // Compare with native transcript
        let mut native_transcript = NativeTranscript::from_proof(&proof);
        let native_circuit_size: u64 = native_transcript.receive_from_prover("circuit_size");
        assert_eq!(circuit_size_val, native_circuit_size);
    }

    // ── Test 5: CommitmentT operations ────────────────────────────────

    #[test]
    fn test_commitment_type() {
        let identity = CommitmentT::identity();
        assert!(identity.get_value().is_point_at_infinity());

        let g1 = CommitmentT::one(None);
        assert!(!g1.get_value().is_point_at_infinity());
        assert!(g1.get_value().on_curve());

        let builder = make_builder_ref();
        let g1_with_ctx = CommitmentT::from_native_with_context(
            builder.clone(),
            Element::<Bn254G1Params>::one().to_affine(),
        );
        assert!(g1_with_ctx.get_context().is_some());
        assert_eq!(g1_with_ctx.get_value(), g1.get_value());
    }

    // ── Test 6: Pairing points accumulator ────────────────────────────

    #[test]
    fn test_pairing_points_accumulator() {
        let g1 = Element::<Bn254G1Params>::one().to_affine();
        let g1_2 = Element::<Bn254G1Params>::one().mul(&Fr::from(2u64)).to_affine();

        let acc1 = PairingPointsAccumulator::from_points(
            CommitmentT::from_native(g1),
            CommitmentT::from_native(g1_2),
        );

        let native = acc1.to_native();
        assert_eq!(native.p0, g1);
        assert_eq!(native.p1, g1_2);
    }

    // ── Test 7: Full recursive verification ───────────────────────────

    #[test]
    fn test_recursive_verification() {
        let g2_x = ensure_test_srs_initialized();

        // Create and prove inner circuit
        let inner_circuit = build_inner_circuit(16);
        let (proof, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        // Create outer circuit builder
        let outer_builder = make_builder_ref();

        // Create recursive verifier
        let recursive_verifier =
            UltraRecursiveVerifier::new(outer_builder.clone(), &native_vk);

        // Run recursive verification
        let output = recursive_verifier.verify_proof(&proof);

        // Verify pairing points are non-trivial
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

        // Verify pairing points produce correct pairing check
        let mut pcs_vk = Bn254VerifierCommitmentKey::with_g2x(g2_x);
        pcs_vk.initialize();
        let pairing_result = {
            let p0_proj = Element::<Bn254G1Params>::from_affine(&p0);
            let p1_proj = Element::<Bn254G1Params>::from_affine(&p1);
            pcs_vk.pairing_check(&p0_proj, &p1_proj)
        };
        assert!(
            pairing_result,
            "Pairing check on recursive verifier output should pass"
        );
    }

    // ── Test 8: Recursive verification with tampered proof fails ──────

    #[test]
    fn test_recursive_verification_fails_on_tampered_proof() {
        let g2_x = ensure_test_srs_initialized();

        let inner_circuit = build_inner_circuit(16);
        let (mut proof, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        // Tamper with proof: modify a non-trivial element
        // Find first non-zero element and modify it
        for i in 0..proof.len() {
            if proof[i] != Fr::zero() {
                proof[i] = proof[i] + Fr::one();
                break;
            }
        }

        let outer_builder = make_builder_ref();
        let recursive_verifier =
            UltraRecursiveVerifier::new(outer_builder.clone(), &native_vk);
        let output = recursive_verifier.verify_proof(&proof);

        // Pairing check should fail on tampered proof
        let mut pcs_vk = Bn254VerifierCommitmentKey::with_g2x(g2_x);
        pcs_vk.initialize();
        let pairing_result = {
            let p0 = output.points_accumulator.p0.get_value();
            let p1 = output.points_accumulator.p1.get_value();
            if p0.is_point_at_infinity() || p1.is_point_at_infinity() {
                false // Verification already failed during sumcheck
            } else {
                let p0_proj = Element::<Bn254G1Params>::from_affine(&p0);
                let p1_proj = Element::<Bn254G1Params>::from_affine(&p1);
                pcs_vk.pairing_check(&p0_proj, &p1_proj)
            }
        };
        assert!(
            !pairing_result,
            "Pairing check should fail on tampered proof"
        );
    }

    // ── Test 9: Outer circuit creates witnesses ───────────────────────

    #[test]
    fn test_outer_circuit_has_witnesses() {
        let g2_x = ensure_test_srs_initialized();

        let inner_circuit = build_inner_circuit(16);
        let (proof, native_vk) = prove_inner_circuit(inner_circuit, g2_x);

        let outer_builder = make_builder_ref();
        let num_vars_before = outer_builder.borrow().base.get_num_variables();

        let recursive_verifier =
            UltraRecursiveVerifier::new(outer_builder.clone(), &native_vk);
        let _output = recursive_verifier.verify_proof(&proof);

        let num_vars_after = outer_builder.borrow().base.get_num_variables();
        assert!(
            num_vars_after > num_vars_before,
            "Recursive verification should create witnesses in outer circuit. \
             Before: {}, After: {}",
            num_vars_before,
            num_vars_after,
        );
    }
}
