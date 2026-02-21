// Group utilities (derive_generators, hash-to-curve)
//
// The actual implementation is in crypto::generators.
// This module re-exports for convenience.

use crate::ecc::curves::grumpkin;

type GrumpkinAffine = grumpkin::G1Affine;

/// Derive `num_generators` independent Grumpkin curve points via BLAKE3 hash-to-curve.
///
/// Delegates to `crypto::generators::derive_generators`.
pub fn derive_generators(
    num_generators: usize,
    domain_separator: &[u8],
) -> Vec<GrumpkinAffine> {
    crate::crypto::generators::derive_generators(domain_separator, num_generators, 0)
}
