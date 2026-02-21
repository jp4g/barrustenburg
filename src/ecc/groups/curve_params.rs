use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;

/// Trait defining an elliptic curve's parameters.
///
/// A curve is defined by y^2 = x^3 + a*x + b over the base field Fq,
/// with scalar field Fr for the prime-order subgroup.
///
/// Constants marked "Montgomery form" are stored as raw Montgomery-form limbs
/// (matching C++ brace-init `fq{...}`). Use `Field::from_raw()` to construct.
///
/// Constants marked "standard form" are stored as raw integer limbs.
/// Use `Field::from_limbs()` to construct (which calls `to_montgomery_form()`).
pub trait CurveParams: 'static + Send + Sync + Sized {
    type BaseFieldParams: FieldParams;
    type ScalarFieldParams: FieldParams;

    /// Whether the curve has a non-zero `a` coefficient (e.g., secp256r1 has a = -3).
    const HAS_A: bool;

    /// Curve coefficient `a` in Montgomery form. Zero if HAS_A is false.
    const A: [u64; 4];

    /// Curve coefficient `b` in Montgomery form.
    const B: [u64; 4];

    /// Generator point x-coordinate in Montgomery form.
    const GENERATOR_X: [u64; 4];

    /// Generator point y-coordinate in Montgomery form.
    const GENERATOR_Y: [u64; 4];

    /// Whether to use the GLV endomorphism for scalar multiplication.
    const USE_ENDOMORPHISM: bool;

    /// Construct the generator point's x-coordinate as a field element.
    /// Default: from_raw (assumes Montgomery form). Override for runtime conversion.
    fn generator_x() -> Field<Self::BaseFieldParams> {
        Field::from_raw(Self::GENERATOR_X)
    }

    /// Construct the generator point's y-coordinate as a field element.
    fn generator_y() -> Field<Self::BaseFieldParams> {
        Field::from_raw(Self::GENERATOR_Y)
    }

    /// Construct curve coefficient a as a field element.
    fn coeff_a() -> Field<Self::BaseFieldParams> {
        Field::from_raw(Self::A)
    }

    /// Construct curve coefficient b as a field element.
    fn coeff_b() -> Field<Self::BaseFieldParams> {
        Field::from_raw(Self::B)
    }
}

// Convenience type aliases
pub type BaseField<C> = Field<<C as CurveParams>::BaseFieldParams>;
pub type ScalarField<C> = Field<<C as CurveParams>::ScalarFieldParams>;
