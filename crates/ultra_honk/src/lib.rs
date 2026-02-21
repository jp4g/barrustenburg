//! Port of the `ultra_honk/` directory — Ultra Honk Prover and Verifier.
//!
//! This crate provides the Ultra Honk proving system:
//! - **Oink prover/verifier**: Pre-sumcheck commitment rounds
//! - **Ultra prover**: Full proof generation (oink + sumcheck + PCS)
//! - **Ultra verifier**: Full proof verification
//!
//! NOTE: Most implementations are structural stubs until Package 6 (Circuit Builder)
//! provides ProverInstance and witness computation infrastructure.

pub mod oink_prover;
pub mod oink_verifier;
pub mod ultra_prover;
pub mod ultra_verifier;
mod validation_tests;

pub use oink_prover::{OinkOutput, OinkProver};
pub use oink_verifier::OinkVerifierOutput;
pub use ultra_prover::UltraProver;
pub use ultra_verifier::{UltraVerifier, UltraVerifierOutput};
