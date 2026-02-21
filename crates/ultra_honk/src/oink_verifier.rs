//! Port of `oink_verifier.hpp`/`.cpp` — Oink Verifier for Ultra Honk.
//!
//! The Oink verifier mirrors the prover's pre-sumcheck rounds:
//! 1. Preamble: receive circuit size, public inputs
//! 2. Wire commitments: receive w_l, w_r, w_o commitments
//! 3. Sorted list accumulator: receive w_4, lookup read counts/tags commitments
//! 4. Log-derivative inverse: receive lookup_inverses commitment
//! 5. Grand product: receive z_perm commitment
//! 6. Generate alpha challenge for sumcheck
//!
//! NOTE: Full implementation requires VerificationKey and commitment point types.
//! The structural types are defined here.

use bbrs_ecc::curves::bn254::Fr;
use bbrs_relations::relation_parameters::RelationParameters;

/// Output data from the Oink verifier rounds.
pub struct OinkVerifierOutput {
    /// Relation parameters with beta, gamma, eta challenges filled in.
    pub relation_parameters: RelationParameters<Fr>,
    /// The alpha challenge for sumcheck.
    pub alpha: Fr,
}
