//! Barretenberg commitment schemes crate.
//!
//! Port of C++ `barretenberg/commitment_schemes/`.
//! Provides polynomial commitment scheme (PCS) types for KZG and IPA.

pub mod claim;
pub mod commitment_key;
pub mod verification_key;
pub mod pairing_points;
pub mod batch_mul;
pub mod claim_batcher;

// Stub modules — not yet implemented
pub mod kzg {}
pub mod ipa {}
pub mod gemini {}
pub mod shplonk {}
pub mod small_subgroup_ipa {}

#[cfg(test)]
mod tests;
