use std::marker::PhantomData;

use crate::ecc::fields::field_params::FieldParams;
#[allow(unused_imports)]
use crate::ecc::fields::field::Field;
use crate::ecc::groups::curve_params::{BaseField, CurveParams};

/// An elliptic curve point in affine coordinates (x, y).
///
/// Point at infinity is represented differently depending on the field modulus size:
/// - Big modulus (>= 2^254): x is set to the modulus value
/// - Small modulus (< 2^254): MSB of x.data[3] is set, all other limbs zero
pub struct AffineElement<C: CurveParams> {
    pub x: BaseField<C>,
    pub y: BaseField<C>,
    _phantom: PhantomData<C>,
}

impl<C: CurveParams> Clone for AffineElement<C> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            x: self.x,
            y: self.y,
            _phantom: PhantomData,
        }
    }
}

impl<C: CurveParams> Copy for AffineElement<C> {}

impl<C: CurveParams> std::fmt::Debug for AffineElement<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_point_at_infinity() {
            write!(f, "AffineElement(infinity)")
        } else {
            write!(f, "AffineElement({:?}, {:?})", self.x, self.y)
        }
    }
}

impl<C: CurveParams> AffineElement<C> {
    /// Construct from coordinates.
    #[inline]
    pub fn new(x: BaseField<C>, y: BaseField<C>) -> Self {
        Self {
            x,
            y,
            _phantom: PhantomData,
        }
    }

    /// The generator point (affine one).
    #[inline]
    pub fn one() -> Self {
        Self::new(C::generator_x(), C::generator_y())
    }

    /// The point at infinity.
    #[inline]
    pub fn infinity() -> Self {
        let mut result = Self::new(BaseField::<C>::zero(), BaseField::<C>::zero());
        result.self_set_infinity();
        result
    }

    /// Check if this point represents the point at infinity.
    #[inline]
    pub fn is_point_at_infinity(&self) -> bool {
        if C::BaseFieldParams::MODULUS_IS_BIG {
            let m = C::BaseFieldParams::MODULUS;
            ((self.x.data[0] ^ m[0])
                | (self.x.data[1] ^ m[1])
                | (self.x.data[2] ^ m[2])
                | (self.x.data[3] ^ m[3]))
                == 0
        } else {
            self.x.is_msb_set()
        }
    }

    /// Set this point to the point at infinity (in-place).
    #[inline]
    pub fn self_set_infinity(&mut self) {
        if C::BaseFieldParams::MODULUS_IS_BIG {
            let m = C::BaseFieldParams::MODULUS;
            self.x.data[0] = m[0];
            self.x.data[1] = m[1];
            self.x.data[2] = m[2];
            self.x.data[3] = m[3];
        } else {
            self.x = BaseField::<C>::zero();
            self.y = BaseField::<C>::zero();
            self.x.self_set_msb();
        }
    }

    /// Return a copy with infinity flag set.
    #[inline]
    pub fn set_infinity(&self) -> Self {
        let mut result = *self;
        result.self_set_infinity();
        result
    }

    /// Try to construct a point from an x-coordinate and a sign bit.
    ///
    /// Computes y^2 = x^3 + a*x + b, attempts sqrt. Returns `None` if
    /// the value is not a quadratic residue (no valid point).
    /// The sign_bit selects which of the two square roots to use:
    /// if the LSB of reduced y doesn't match sign_bit, y is negated.
    pub fn from_x_coordinate(x: BaseField<C>, sign_bit: bool) -> Option<Self> {
        let mut yy = x.sqr() * x + C::coeff_b();
        if C::HAS_A {
            yy = yy + (x * C::coeff_a());
        }
        let (found_root, y) = yy.sqrt();
        if found_root {
            let y_reduced = y.from_montgomery_form();
            let y_parity = y_reduced.data[0] & 1 == 1;
            let y = if y_parity != sign_bit { -y } else { y };
            Some(Self::new(x, y))
        } else {
            None
        }
    }

    /// Check if the point lies on the curve: y^2 == x^3 + a*x + b.
    pub fn on_curve(&self) -> bool {
        if self.is_point_at_infinity() {
            return true;
        }
        let mut xxx = self.x.sqr() * self.x + C::coeff_b();
        if C::HAS_A {
            xxx = xxx + (self.x * C::coeff_a());
        }
        let yy = self.y.sqr();
        xxx == yy
    }
}

// ---------------------------------------------------------------------------
// Operator impls
// ---------------------------------------------------------------------------

impl<C: CurveParams> PartialEq for AffineElement<C> {
    fn eq(&self, other: &Self) -> bool {
        let this_inf = self.is_point_at_infinity();
        let other_inf = other.is_point_at_infinity();
        let both_inf = this_inf && other_inf;
        let only_one_inf = this_inf != other_inf;
        !only_one_inf && (both_inf || (self.x == other.x && self.y == other.y))
    }
}

impl<C: CurveParams> Eq for AffineElement<C> {}

impl<C: CurveParams> From<crate::ecc::groups::element::Element<C>> for AffineElement<C> {
    fn from(e: crate::ecc::groups::element::Element<C>) -> Self {
        e.to_affine()
    }
}

impl<C: CurveParams> std::ops::Neg for AffineElement<C> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(self.x, -self.y)
    }
}
