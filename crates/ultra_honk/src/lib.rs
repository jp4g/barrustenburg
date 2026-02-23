//! Port of the `ultra_honk/` directory — Ultra Honk Prover and Verifier.
//!
//! This crate provides the Ultra Honk proving system:
//! - **Proving key**: Construction from UltraCircuitBuilder
//! - **Verification key**: Commitments to precomputed polynomials
//! - **Oink prover/verifier**: Pre-sumcheck commitment rounds
//! - **Ultra sumcheck**: Full sumcheck protocol for Ultra flavor
//! - **Ultra prover**: Full proof generation (oink + sumcheck + PCS)
//! - **Ultra verifier**: Full proof verification

pub mod oink_prover;
pub mod oink_verifier;
pub mod proving_key;
pub mod ultra_prover;
pub mod ultra_sumcheck;
pub mod ultra_verifier;
pub mod verification_key;

#[cfg(test)]
mod validation_tests;
#[cfg(test)]
mod e2e_tests;

pub use oink_prover::{OinkOutput, OinkProver};
pub use oink_verifier::OinkVerifierOutput;
pub use proving_key::ProvingKey;
pub use ultra_prover::UltraProver;
pub use ultra_verifier::{UltraVerifier, UltraVerifierOutput};
pub use verification_key::VerificationKey;
