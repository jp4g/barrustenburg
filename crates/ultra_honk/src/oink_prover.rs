//! Port of `oink_prover.hpp`/`.cpp` — Oink Prover for Ultra Honk.
//!
//! The Oink prover executes the pre-sumcheck rounds:
//! 1. Preamble: send circuit size, public inputs
//! 2. Wire commitments: commit to w_l, w_r, w_o
//! 3. Sorted list accumulator: compute and commit to w_4, lookup read counts/tags
//! 4. Log-derivative inverse: compute and commit to lookup_inverses
//! 5. Grand product: compute and commit to z_perm
//! 6. Generate alpha challenge for sumcheck
//!
//! NOTE: Full implementation requires ProverInstance (with circuit builder / witness
//! computation infrastructure from Package 6). The types and orchestration are defined
//! here; the actual prove() method will be completed when Package 6 is available.

use bbrs_ecc::curves::bn254::Fr;
use bbrs_flavor::ultra_flavor;
use bbrs_relations::relation_parameters::RelationParameters;

/// Output data from the Oink prover rounds.
///
/// Contains the relation parameters (with challenges) and the alpha separator
/// needed for sumcheck.
pub struct OinkOutput {
    /// Relation parameters with beta, gamma, eta challenges filled in.
    pub relation_parameters: RelationParameters<Fr>,
    /// The alpha challenge used to batch subrelations in sumcheck.
    pub alpha: Fr,
}

/// Oink prover for Ultra Honk.
///
/// Port of C++ `OinkProver<UltraFlavor>`.
///
/// The Oink phase computes and commits to witness polynomials before sumcheck.
/// Full implementation requires ProverInstance from Package 6.
pub struct OinkProver {
    /// The prover polynomials (filled by witness computation).
    pub polynomials: ultra_flavor::ProverPolynomials<bbrs_ecc::curves::bn254::Bn254FrParams>,
    /// Public inputs for the circuit.
    pub public_inputs: Vec<Fr>,
    /// Relation parameters (will be filled with challenges during prove()).
    pub relation_parameters: RelationParameters<Fr>,
}

impl OinkProver {
    /// Create a new OinkProver.
    ///
    /// NOTE: In the full implementation, this takes a ProverInstance.
    /// For now, it takes pre-constructed polynomials.
    pub fn new(
        polynomials: ultra_flavor::ProverPolynomials<bbrs_ecc::curves::bn254::Bn254FrParams>,
        public_inputs: Vec<Fr>,
    ) -> Self {
        Self {
            polynomials,
            public_inputs,
            relation_parameters: RelationParameters::default(),
        }
    }
}
