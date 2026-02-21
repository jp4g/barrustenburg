//! Batch scalar multiplication utility.
//!
//! Port of C++ `commitment_schemes/utils/batch_mul_native.hpp`.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;
use bbrs_ecc::scalar_multiplication::pippenger_msm;

/// Wrapper around Pippenger MSM for native batch multiplication.
///
/// Computes sum_i scalars[i] * points[i] and returns the result as an
/// affine element.
///
/// Port of C++ `batch_mul_native<Curve>()`.
pub fn batch_mul_native<C: CurveParams>(
    points: &[AffineElement<C>],
    scalars: &[Field<C::ScalarFieldParams>],
) -> AffineElement<C> {
    assert_eq!(
        points.len(),
        scalars.len(),
        "batch_mul_native: points and scalars must have the same length"
    );

    if points.is_empty() {
        return AffineElement::infinity();
    }

    // Copy scalars since pippenger_msm may mutate them internally
    let scalars_copy: Vec<_> = scalars.to_vec();
    pippenger_msm(&scalars_copy, points).to_affine()
}
