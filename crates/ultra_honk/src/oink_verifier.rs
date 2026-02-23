//! Port of `oink_verifier.hpp`/`.cpp` — Oink Verifier for Ultra Honk.
//!
//! The Oink verifier mirrors the prover's pre-sumcheck rounds:
//! 1. Preamble: receive circuit size, public inputs
//! 2. Wire commitments: receive w_l, w_r, w_o commitments
//! 3. Sorted list accumulator: receive w_4, lookup_read_counts/tags commitments
//! 4. Log-derivative inverse: receive lookup_inverses commitment
//! 5. Grand product: receive z_perm commitment
//! 6. Generate alpha challenge for sumcheck

use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr, G1Affine};
use bbrs_honk::compute_public_input_delta;
use bbrs_relations::relation_parameters::RelationParameters;
use bbrs_transcript::NativeTranscript;

use crate::verification_key::VerificationKey;

/// Output data from the Oink verifier rounds.
pub struct OinkVerifierOutput {
    /// Relation parameters with beta, gamma, eta challenges filled in.
    pub relation_parameters: RelationParameters<Fr>,
    /// The alpha challenge for sumcheck.
    pub alpha: Fr,
    /// Witness commitments received from prover.
    pub witness_commitments: WitnessCommitments,
}

/// Commitments to witness polynomials received during Oink.
pub struct WitnessCommitments {
    pub w_l: G1Affine,
    pub w_r: G1Affine,
    pub w_o: G1Affine,
    pub w_4: G1Affine,
    pub lookup_read_counts: G1Affine,
    pub lookup_read_tags: G1Affine,
    pub lookup_inverses: G1Affine,
    pub z_perm: G1Affine,
}

/// Oink verifier for Ultra Honk.
pub struct OinkVerifier<'a> {
    vk: &'a VerificationKey,
}

impl<'a> OinkVerifier<'a> {
    pub fn new(vk: &'a VerificationKey) -> Self {
        Self { vk }
    }

    /// Execute the Oink verification rounds.
    pub fn verify(self, transcript: &mut NativeTranscript) -> OinkVerifierOutput {
        let _vk = self.vk;

        // Round 1: Preamble — receive circuit size, num public inputs, public inputs
        let _circuit_size: u64 = transcript.receive_from_prover("circuit_size");
        let num_public_inputs: u64 = transcript.receive_from_prover("public_input_size");
        let pub_inputs_offset: u64 = transcript.receive_from_prover("pub_inputs_offset");

        let mut public_inputs = Vec::with_capacity(num_public_inputs as usize);
        for i in 0..num_public_inputs as usize {
            let pi: Fr = transcript.receive_from_prover(&format!("public_input_{}", i));
            public_inputs.push(pi);
        }

        // Round 2: Receive wire commitments
        let w_l: G1Affine = transcript.receive_from_prover("W_L");
        let w_r: G1Affine = transcript.receive_from_prover("W_R");
        let w_o: G1Affine = transcript.receive_from_prover("W_O");

        // Get eta challenges
        let eta = transcript.get_challenge("eta");
        let eta_two = transcript.get_challenge("eta_two");
        let eta_three = transcript.get_challenge("eta_three");

        // Round 3: Receive sorted list accumulator commitments
        let w_4: G1Affine = transcript.receive_from_prover("W_4");
        let lookup_read_counts: G1Affine =
            transcript.receive_from_prover("LOOKUP_READ_COUNTS");
        let lookup_read_tags: G1Affine =
            transcript.receive_from_prover("LOOKUP_READ_TAGS");

        // Get beta, gamma challenges
        let beta = transcript.get_challenge("beta");
        let gamma = transcript.get_challenge("gamma");

        // Compute public input delta
        let pub_inputs_offset_field = Fr::from(pub_inputs_offset);
        let public_input_delta = compute_public_input_delta::<Bn254FrParams>(
            &public_inputs,
            beta,
            gamma,
            pub_inputs_offset_field,
        );

        let beta_sqr = beta * beta;
        let beta_cube = beta_sqr * beta;

        let relation_parameters = RelationParameters {
            eta,
            eta_two,
            eta_three,
            beta,
            gamma,
            public_input_delta,
            beta_sqr,
            beta_cube,
            ..RelationParameters::default()
        };

        // Round 4: Receive lookup inverses commitment
        let lookup_inverses: G1Affine =
            transcript.receive_from_prover("LOOKUP_INVERSES");

        // Round 5: Receive z_perm commitment
        let z_perm: G1Affine = transcript.receive_from_prover("Z_PERM");

        // Round 6: Get alpha challenge
        let alpha = transcript.get_challenge("Sumcheck:alpha");

        OinkVerifierOutput {
            relation_parameters,
            alpha,
            witness_commitments: WitnessCommitments {
                w_l,
                w_r,
                w_o,
                w_4,
                lookup_read_counts,
                lookup_read_tags,
                lookup_inverses,
                z_perm,
            },
        }
    }
}
