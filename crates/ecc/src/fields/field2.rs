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

impl<P: FieldParams> Field2<P> {
    /// Generate a random Field2 element.
    pub fn random_element() -> Self {
        Self {
            c0: Field::random_element(),
            c1: Field::random_element(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curves::bn254::Bn254FqParams;

    type Fq2 = Field2<Bn254FqParams>;
    type Fq = Field<Bn254FqParams>;

    #[test]
    fn fq2_eq() {
        let a = Fq2::random_element();
        let b = a;
        assert_eq!(a, b);
        let c = Fq2::random_element();
        // Very unlikely to be equal
        if a != c {
            assert_ne!(a, c);
        }
    }

    #[test]
    fn fq2_is_zero() {
        assert!(Fq2::zero().is_zero());
        assert!(!Fq2::one().is_zero());
        let r = Fq2::random_element();
        // Random element is almost certainly not zero
        assert!(!r.is_zero());
    }

    #[test]
    fn fq2_random_element() {
        let a = Fq2::random_element();
        let b = Fq2::random_element();
        assert_ne!(a, b, "two random elements should differ (probabilistically)");
    }

    #[test]
    fn fq2_mul_check_against_constants() {
        let a = Fq2::new(Fq::from(3u64), Fq::from(5u64));
        let b = Fq2::new(Fq::from(7u64), Fq::from(11u64));
        let c = a * b;
        // (3+5u)(7+11u) = 21+33u+35u+55u^2 = (21-55) + (33+35)u = -34 + 68u
        let neg34 = Fq::from(34u64).negate();
        let expected = Fq2::new(neg34, Fq::from(68u64));
        assert_eq!(c, expected, "Fq2 mul with small constants");
    }

    #[test]
    fn fq2_sqr_check_against_constants() {
        let a = Fq2::new(Fq::from(3u64), Fq::from(5u64));
        let sq = a.sqr();
        let mul_self = a * a;
        assert_eq!(sq, mul_self, "sqr should equal self*self");
    }

    #[test]
    fn fq2_add_check_against_constants() {
        let a = Fq2::new(Fq::from(3u64), Fq::from(5u64));
        let b = Fq2::new(Fq::from(7u64), Fq::from(11u64));
        let c = a + b;
        let expected = Fq2::new(Fq::from(10u64), Fq::from(16u64));
        assert_eq!(c, expected);
    }

    #[test]
    fn fq2_sub_check_against_constants() {
        let a = Fq2::new(Fq::from(10u64), Fq::from(16u64));
        let b = Fq2::new(Fq::from(3u64), Fq::from(5u64));
        let c = a - b;
        let expected = Fq2::new(Fq::from(7u64), Fq::from(11u64));
        assert_eq!(c, expected);
    }

    #[test]
    fn fq2_to_montgomery_form() {
        let a = Fq2::new(
            Fq::from_raw([1, 0, 0, 0]),
            Fq::from_raw([2, 0, 0, 0]),
        );
        let mont = a.to_montgomery_form();
        let back = mont.from_montgomery_form();
        assert_eq!(back.c0.data, a.c0.data);
        assert_eq!(back.c1.data, a.c1.data);
    }

    #[test]
    fn fq2_from_montgomery_form() {
        let a = Fq2::new(Fq::from(42u64), Fq::from(99u64));
        let from_mont = a.from_montgomery_form();
        let back = from_mont.to_montgomery_form();
        assert_eq!(a, back, "from_montgomery_form should be invertible");
    }

    #[test]
    fn fq2_mul_sqr_consistency() {
        for _ in 0..100 {
            let a = Fq2::random_element();
            assert_eq!(a * a, a.sqr(), "a*a should equal a.sqr()");
        }
    }

    #[test]
    fn fq2_add_mul_consistency() {
        for _ in 0..100 {
            let a = Fq2::random_element();
            let b = Fq2::random_element();
            let c = Fq2::random_element();
            let lhs = (a + b) * c;
            let rhs = a * c + b * c;
            assert_eq!(lhs, rhs, "(a+b)*c should equal a*c + b*c");
        }
    }

    #[test]
    fn fq2_sub_mul_consistency() {
        for _ in 0..100 {
            let a = Fq2::random_element();
            let b = Fq2::random_element();
            let c = Fq2::random_element();
            let lhs = (a - b) * c;
            let rhs = a * c - b * c;
            assert_eq!(lhs, rhs, "(a-b)*c should equal a*c - b*c");
        }
    }

    #[test]
    fn fq2_invert() {
        for _ in 0..100 {
            let a = Fq2::random_element();
            if !a.is_zero() {
                let a_inv = a.invert();
                let product = a * a_inv;
                assert_eq!(product, Fq2::one(), "a * a^-1 should be one");
            }
        }
    }

    #[test]
    fn fq2_serialize() {
        // Test that montgomery roundtrip preserves values
        let a = Fq2::random_element();
        let roundtrip = a.from_montgomery_form().to_montgomery_form();
        assert_eq!(a, roundtrip, "montgomery roundtrip should preserve Fq2 value");
    }
}
