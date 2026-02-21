//! Port of `ultra_verifier.hpp`/`.cpp` — Ultra Honk Verifier.
//!
//! The Ultra verifier orchestrates the full proof verification:
//! 1. Oink phase: receive commitments, generate challenges
//! 2. Sumcheck phase: run sumcheck verifier
//! 3. PCS phase: run Shplemini/SHPLONK verification + pairing check
//!
//! NOTE: Full implementation requires VerificationKey, VerifierCommitments,
//! and the complete Oink verifier. The structural types are defined here.

/// Output of the Ultra Honk verifier.
///
/// Port of C++ `UltraVerifier_::UltraVerifierOutput`.
pub struct UltraVerifierOutput {
    /// Whether the proof verified successfully.
    pub result: bool,
}

/// Ultra Honk verifier.
///
/// Port of C++ `UltraVerifier_<UltraFlavor>`.
pub struct UltraVerifier;

// Full implementation deferred to Package 6 (Circuit Builder).
// The verify_proof() method will orchestrate:
//   1. oink_verifier.verify()
//   2. sumcheck_verifier.verify()
//   3. shplemini_verifier.compute_batch_opening_claim()
//   4. PCS::reduce_verify_batch_opening_claim()
//   5. pairing_points.check()
