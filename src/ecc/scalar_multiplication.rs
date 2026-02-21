// STUB: Multi-Scalar Multiplication (Pippenger algorithm)
//
// C++ source: ecc/scalar_multiplication/scalar_multiplication.hpp/cpp
//
// Provides:
// - MSM<Curve>::batched_mul(): Pippenger multi-scalar multiplication
// - Dynamic bucket width selection based on point count
// - Thread-parallel bucket accumulation with affine trick
// - Small-point optimization for < ~100 points
//
// Components to port:
// - scalar_multiplication.hpp — MSM class, work unit distribution
// - process_buckets.hpp — radix sort and bucket ordering
// - bitvector.hpp — bit-packed vector for bucket existence tracking
//
// Required by: polynomial commitments (KZG, IPA), prover inner loops
// Depends on: WNAF (groups/wnaf.rs), batched affine addition
// Priority: Performance optimization (critical for prover speed)

use crate::ecc::groups::affine_element::AffineElement;
use crate::ecc::groups::element::Element;
use crate::ecc::groups::curve_params::CurveParams;
use crate::ecc::fields::field::Field;

/// Compute sum(scalars[i] * points[i]) using Pippenger's bucket method.
///
/// TODO: Port from `ecc/scalar_multiplication/scalar_multiplication.hpp`.
/// Current fallback: naive double-and-add per point.
pub fn pippenger_msm<C: CurveParams>(
    _scalars: &[Field<C::ScalarFieldParams>],
    _points: &[AffineElement<C>],
) -> Element<C> {
    todo!("port Pippenger MSM from C++")
}

/// Naive MSM fallback: sum of individual scalar multiplications.
pub fn naive_msm<C: CurveParams>(
    scalars: &[Field<C::ScalarFieldParams>],
    points: &[AffineElement<C>],
) -> Element<C> {
    assert_eq!(scalars.len(), points.len());
    let mut acc = Element::<C>::infinity();
    for (s, p) in scalars.iter().zip(points.iter()) {
        let proj = Element::from_affine(p);
        acc += proj.mul_without_endomorphism(s);
    }
    acc
}
