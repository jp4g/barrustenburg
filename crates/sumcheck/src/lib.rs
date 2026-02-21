//! Port of barretenberg's sumcheck protocol.
//!
//! Implements the sumcheck prover and verifier for SumcheckTestFlavor.

pub mod sumcheck;
pub mod sumcheck_output;
pub mod sumcheck_round;

#[cfg(test)]
mod tests;
