// STUB: Batched affine point addition
//
// C++ source: ecc/batched_affine_addition/batched_affine_addition.hpp/cpp
//
// Provides:
// - BatchedAffineAddition<Curve>: reduce sequences of affine points via batched addition
// - Montgomery batch inversion trick: compute all 1/(x2-x1) denominators in one batch
// - Thread-parallel sequence reduction
//
// Required by: Pippenger MSM (for bucket accumulation in affine coordinates)
// Depends on: basic field arithmetic (already ported)
// Priority: Performance optimization (medium — only needed for MSM perf)

use crate::ecc::groups::affine_element::AffineElement;
use crate::ecc::groups::curve_params::CurveParams;

/// Batch-add a sequence of affine points using Montgomery's batch inversion trick.
///
/// Given points [P0, P1, P2, ...], computes P0 + P1 + P2 + ... efficiently
/// by batching all the modular inversions needed for affine addition.
///
/// TODO: Port from `ecc/batched_affine_addition/batched_affine_addition.hpp`.
pub fn batched_affine_add<C: CurveParams>(
    _points: &[AffineElement<C>],
) -> AffineElement<C> {
    todo!("port batched affine addition from C++")
}
