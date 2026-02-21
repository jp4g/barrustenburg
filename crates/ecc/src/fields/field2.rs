// Quadratic extension field F_{p^2} = F_p[u] / (u^2 + 1)
//
// C++ source: ecc/fields/field2.hpp, field2_declarations.hpp
//
// Elements are pairs (c0, c1) representing c0 + c1*u where u^2 = -1.
// Requires p ≡ 3 (mod 4) so that -1 is not a quadratic residue.

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::field::Field;
use super::field_params::FieldParams;

pub struct Field2<P: FieldParams> {
    pub c0: Field<P>,
    pub c1: Field<P>,
}

impl<P: FieldParams> Clone for Field2<P> {
    #[inline]
    fn clone(&self) -> Self {
        Self { c0: self.c0, c1: self.c1 }
    }
}

impl<P: FieldParams> Copy for Field2<P> {}

impl<P: FieldParams> std::fmt::Debug for Field2<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Field2({:?}, {:?})", self.c0, self.c1)
    }
}

impl<P: FieldParams> Field2<P> {
    #[inline]
    pub fn new(c0: Field<P>, c1: Field<P>) -> Self {
        Self { c0, c1 }
    }

    #[inline]
    pub fn zero() -> Self {
        Self {
            c0: Field::zero(),
            c1: Field::zero(),
        }
    }

    #[inline]
    pub fn one() -> Self {
        Self {
            c0: Field::one(),
            c1: Field::zero(),
        }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// Multiply each component by a base field element.
    #[inline]
    pub fn mul_by_fq(&self, a: Field<P>) -> Self {
        Self {
            c0: a * self.c0,
            c1: a * self.c1,
        }
    }

    /// Squaring: (c0 + c1*u)^2 = (c0+c1)(c0-c1) + 2*c0*c1*u
    #[inline]
    pub fn sqr(&self) -> Self {
        let t1 = self.c0 * self.c1;
        Self {
            c0: (self.c0 + self.c1) * (self.c0 - self.c1),
            c1: t1 + t1,
        }
    }

    /// Inversion: 1/(c0 + c1*u) = (c0 - c1*u) / (c0^2 + c1^2)
    pub fn invert(&self) -> Self {
        let t3 = (self.c0.sqr() + self.c1.sqr()).invert();
        Self {
            c0: self.c0 * t3,
            c1: -(self.c1 * t3),
        }
    }

    /// Frobenius map: conjugation (c0, -c1)
    #[inline]
    pub fn frobenius_map(&self) -> Self {
        Self {
            c0: self.c0,
            c1: -self.c1,
        }
    }

    #[inline]
    pub fn to_montgomery_form(&self) -> Self {
        Self {
            c0: self.c0.to_montgomery_form(),
            c1: self.c1.to_montgomery_form(),
        }
    }

    #[inline]
    pub fn from_montgomery_form(&self) -> Self {
        Self {
            c0: self.c0.from_montgomery_form(),
            c1: self.c1.from_montgomery_form(),
        }
    }
}

// ---------------------------------------------------------------------------
// Operator impls
// ---------------------------------------------------------------------------

impl<P: FieldParams> Add for Field2<P> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}

impl<P: FieldParams> AddAssign for Field2<P> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<P: FieldParams> Sub for Field2<P> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}

impl<P: FieldParams> SubAssign for Field2<P> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<P: FieldParams> Mul for Field2<P> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        // Karatsuba: (c0 + c1*u)(d0 + d1*u)
        // = (c0*d0 - c1*d1) + ((c0+c1)*(d0+d1) - c0*d0 - c1*d1)*u
        let t1 = self.c0 * rhs.c0;
        let t2 = self.c1 * rhs.c1;
        let t3 = self.c0 + self.c1;
        let t4 = rhs.c0 + rhs.c1;
        Self {
            c0: t1 - t2,
            c1: t3 * t4 - (t1 + t2),
        }
    }
}

impl<P: FieldParams> MulAssign for Field2<P> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<P: FieldParams> Neg for Field2<P> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}

impl<P: FieldParams> PartialEq for Field2<P> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<P: FieldParams> Eq for Field2<P> {}
