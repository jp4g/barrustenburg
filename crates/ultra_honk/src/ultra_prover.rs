//! Port of `ultra_prover.hpp`/`.cpp` — Ultra Honk Prover.
//!
//! The Ultra prover orchestrates the full proof generation:
//! 1. Oink phase: commit to witness polynomials, compute challenges
//! 2. Sumcheck phase: run sumcheck protocol
//! 3. PCS phase: run Shplemini/SHPLONK polynomial commitment opening
//!
//! NOTE: Full implementation requires ProverInstance, CommitmentKey, and the
//! complete Oink prover. The structural types are defined here.

/// Ultra Honk prover.
///
/// Port of C++ `UltraProver<UltraFlavor>`.
pub struct UltraProver;

// Full implementation deferred to Package 6 (Circuit Builder).
// The prove() method will orchestrate:
//   1. oink_prover.prove()
//   2. sumcheck_prover.prove()
//   3. shplemini_prover.prove()
