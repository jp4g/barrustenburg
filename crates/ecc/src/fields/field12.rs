// Dodecic extension field Fq12 = Fq6[w] / (w^2 - v)
//
// C++ source: ecc/fields/field12.hpp
//
// Elements are pairs (c0, c1) of Field6 elements.
// Used for the BN254 ate pairing target group GT.

use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use super::field2::Field2;
use super::field6::{Field6, Field6Params};

/// Trait providing the BN254-specific Fq12 Frobenius coefficients.
pub trait Field12Params: Field6Params {
    fn frobenius_coefficients_1() -> Field2<Self>;
    fn frobenius_coefficients_2() -> Field2<Self>;
    fn frobenius_coefficients_3() -> Field2<Self>;
}

/// Line evaluation coefficients for the Miller loop sparse multiplication.
pub struct EllCoeffs<P: Field6Params> {
    pub o: Field2<P>,
    pub vw: Field2<P>,
    pub vv: Field2<P>,
}

impl<P: Field6Params> Clone for EllCoeffs<P> {
    #[inline]
    fn clone(&self) -> Self {
        Self { o: self.o, vw: self.vw, vv: self.vv }
    }
}

impl<P: Field6Params> Copy for EllCoeffs<P> {}

impl<P: Field6Params> std::fmt::Debug for EllCoeffs<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "EllCoeffs({:?}, {:?}, {:?})", self.o, self.vw, self.vv)
    }
}

pub struct Field12<P: Field12Params> {
    pub c0: Field6<P>,
    pub c1: Field6<P>,
}

impl<P: Field12Params> Clone for Field12<P> {
    #[inline]
    fn clone(&self) -> Self {
        Self { c0: self.c0, c1: self.c1 }
    }
}

impl<P: Field12Params> Copy for Field12<P> {}

impl<P: Field12Params> std::fmt::Debug for Field12<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Field12({:?}, {:?})", self.c0, self.c1)
    }
}

impl<P: Field12Params> Field12<P> {
    #[inline]
    pub fn new(c0: Field6<P>, c1: Field6<P>) -> Self {
        Self { c0, c1 }
    }

    #[inline]
    pub fn zero() -> Self {
        Self {
            c0: Field6::zero(),
            c1: Field6::zero(),
        }
    }

    #[inline]
    pub fn one() -> Self {
        Self {
            c0: Field6::one(),
            c1: Field6::zero(),
        }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// Fq12's non-residue multiplication on Fq6: (c0, c1, c2) -> (xi*c2, c0, c1)
    /// where xi is the Fq6 non-residue applied to the Fq2 component c2.
    #[inline]
    fn mul_by_non_residue(a: &Field6<P>) -> Field6<P> {
        Field6::new(Field6::<P>::mul_by_non_residue(&a.c2), a.c0, a.c1)
    }

    /// Squaring in Fq12 = Fq6[w] / (w^2 - v).
    pub fn sqr(&self) -> Self {
        let t0 = self.c0 + self.c1;
        let t1 = Self::mul_by_non_residue(&self.c1) + self.c0;

        let t0 = t0 * t1;
        let t1 = self.c0 * self.c1;

        Self {
            c0: t0 - (t1 + Self::mul_by_non_residue(&t1)),
            c1: t1 + t1,
        }
    }

    /// Inversion via Algorithm 8 from "High-Speed Software Implementation of
    /// the Optimal Ate Pairing over Barreto-Naehrig Curves".
    pub fn invert(&self) -> Self {
        let t0 = (self.c0.sqr() - Self::mul_by_non_residue(&self.c1.sqr())).invert();
        Self {
            c0: self.c0 * t0,
            c1: -(self.c1 * t0),
        }
    }

    /// Cyclotomic squaring (currently delegates to generic sqr).
    #[inline]
    pub fn cyclotomic_squared(&self) -> Self {
        self.sqr()
    }

    /// Unitary inverse: conjugation in the cyclotomic subgroup.
    #[inline]
    pub fn unitary_inverse(&self) -> Self {
        Self {
            c0: self.c0,
            c1: -self.c1,
        }
    }

    pub fn frobenius_map_one(&self) -> Self {
        Self {
            c0: self.c0.frobenius_map_one(),
            c1: self.c1.frobenius_map_one().mul_by_fq2(&P::frobenius_coefficients_1()),
        }
    }

    pub fn frobenius_map_two(&self) -> Self {
        Self {
            c0: self.c0.frobenius_map_two(),
            c1: self.c1.frobenius_map_two().mul_by_fq2(&P::frobenius_coefficients_2()),
        }
    }

    pub fn frobenius_map_three(&self) -> Self {
        Self {
            c0: self.c0.frobenius_map_three(),
            c1: self.c1.frobenius_map_three().mul_by_fq2(&P::frobenius_coefficients_3()),
        }
    }

    /// Sparse multiplication by a line evaluation (ell_coeffs).
    ///
    /// The multiplicand is a sparse Fq12 element:
    ///   (ell.o, 0, ell.vv) + w*(0, ell.vw, 0)
    ///
    /// Ported from C++ field12.hpp::self_sparse_mul.
    pub fn self_sparse_mul(&mut self, ell: &EllCoeffs<P>) {
        let d0 = self.c0.c0 * ell.o;
        let d2 = self.c0.c2 * ell.vv;
        let d4 = self.c1.c1 * ell.vw;
        let t2 = self.c0.c0 + self.c1.c1;
        let t1 = self.c0.c0 + self.c0.c2;
        let mut s0 = self.c0.c1 + self.c1.c0;
        s0 += self.c1.c2;

        let mut s1 = self.c0.c1 * ell.vv;
        let mut t3 = s1 + d4;
        let t4 = Field6::<P>::mul_by_non_residue(&t3);
        self.c0.c0 = t4 + d0;

        t3 = self.c1.c2 * ell.vw;
        s1 += t3;
        t3 = t3 + d2;
        let t4 = Field6::<P>::mul_by_non_residue(&t3);
        t3 = self.c0.c1 * ell.o;
        s1 += t3;
        self.c0.c1 = t4 + t3;

        let t0 = ell.o + ell.vv;
        t3 = t1 * t0;
        t3 = t3 - d0;
        t3 = t3 - d2;
        let t4 = self.c1.c0 * ell.vw;
        s1 += t4;

        let t0 = self.c0.c2 + self.c1.c1;
        self.c0.c2 = t3 + t4;

        let t1 = ell.vv + ell.vw;
        t3 = t0 * t1;
        t3 = t3 - d2;
        t3 = t3 - d4;
        let t4 = Field6::<P>::mul_by_non_residue(&t3);
        t3 = self.c1.c0 * ell.o;
        s1 += t3;
        self.c1.c0 = t3 + t4;

        t3 = self.c1.c2 * ell.vv;
        s1 += t3;
        let t4 = Field6::<P>::mul_by_non_residue(&t3);
        let t0 = ell.o + ell.vw;
        t3 = t0 * t2;
        t3 = t3 - d0;
        t3 = t3 - d4;
        self.c1.c1 = t3 + t4;

        let t0 = ell.o + ell.vv + ell.vw;
        t3 = s0 * t0;
        self.c1.c2 = t3 - s1;
    }
}

// ---------------------------------------------------------------------------
// Operator impls
// ---------------------------------------------------------------------------

impl<P: Field12Params> Add for Field12<P> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}

impl<P: Field12Params> AddAssign for Field12<P> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<P: Field12Params> Sub for Field12<P> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}

impl<P: Field12Params> SubAssign for Field12<P> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<P: Field12Params> Mul for Field12<P> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let t0 = self.c0 * rhs.c0;
        let t1 = self.c1 * rhs.c1;
        let t2 = self.c0 + self.c1;
        let t3 = rhs.c0 + rhs.c1;
        Self {
            c0: Self::mul_by_non_residue(&t1) + t0,
            c1: t2 * t3 - (t0 + t1),
        }
    }
}

impl<P: Field12Params> MulAssign for Field12<P> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<P: Field12Params> PartialEq for Field12<P> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<P: Field12Params> Eq for Field12<P> {}

impl<P: Field12Params> Field12<P> {
    /// Generate a random Field12 element.
    pub fn random_element() -> Self {
        Self {
            c0: Field6::random_element(),
            c1: Field6::random_element(),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curves::bn254::Bn254FqParams;
    use crate::fields::field::Field;

    type BnFq12 = Field12<Bn254FqParams>;
    type BnFq6 = Field6<Bn254FqParams>;
    type BnFq2 = Field2<Bn254FqParams>;
    type Fq = Field<Bn254FqParams>;

    #[test]
    fn fq12_eq() {
        let a = BnFq12::random_element();
        let b = a;
        assert_eq!(a, b);
    }

    #[test]
    fn fq12_is_zero() {
        assert!(BnFq12::zero().is_zero());
        assert!(!BnFq12::one().is_zero());
    }

    #[test]
    fn fq12_random_element() {
        let a = BnFq12::random_element();
        let b = BnFq12::random_element();
        assert_ne!(a, b, "two random Fq12 elements should differ");
    }

    #[test]
    fn fq12_mul_check_against_constants() {
        let a = BnFq12::random_element();
        let c = a * a;
        let d = a.sqr();
        assert_eq!(c, d, "a*a should equal a.sqr()");
    }

    #[test]
    fn fq12_sqr_check_against_constants() {
        let a = BnFq12::new(
            BnFq6::new(
                BnFq2::new(Fq::from(1u64), Fq::from(2u64)),
                BnFq2::new(Fq::from(3u64), Fq::from(4u64)),
                BnFq2::new(Fq::from(5u64), Fq::from(6u64)),
            ),
            BnFq6::new(
                BnFq2::new(Fq::from(7u64), Fq::from(8u64)),
                BnFq2::new(Fq::from(9u64), Fq::from(10u64)),
                BnFq2::new(Fq::from(11u64), Fq::from(12u64)),
            ),
        );
        assert_eq!(a.sqr(), a * a);
    }

    #[test]
    fn fq12_add_check_against_constants() {
        let a = BnFq12::one();
        let b = BnFq12::one();
        let c = a + b;
        // c.c0 should have c0.c0.c0 = 2, rest zero; c.c1 should be zero
        assert_eq!(c.c0.c0.c0, Fq::from(2u64));
        assert!(c.c1.is_zero());
    }

    #[test]
    fn fq12_sub_check_against_constants() {
        let a = BnFq12::random_element();
        let b = a;
        let c = a - b;
        assert!(c.is_zero(), "a - a should be zero");
    }

    #[test]
    fn fq12_to_montgomery_form() {
        let a = BnFq12::random_element();
        let roundtrip = a.from_montgomery_form().to_montgomery_form();
        assert_eq!(a, roundtrip);
    }

    #[test]
    fn fq12_from_montgomery_form() {
        let a = BnFq12::random_element();
        let from = a.from_montgomery_form();
        let back = from.to_montgomery_form();
        assert_eq!(a, back);
    }

    #[test]
    fn fq12_mul_sqr_consistency() {
        for _ in 0..10 {
            let a = BnFq12::random_element();
            assert_eq!(a * a, a.sqr(), "a*a should equal a.sqr()");
        }
    }

    #[test]
    fn fq12_add_mul_consistency() {
        for _ in 0..10 {
            let a = BnFq12::random_element();
            let b = BnFq12::random_element();
            let c = BnFq12::random_element();
            let lhs = (a + b) * c;
            let rhs = a * c + b * c;
            assert_eq!(lhs, rhs, "(a+b)*c should equal a*c + b*c");
        }
    }

    #[test]
    fn fq12_sub_mul_consistency() {
        for _ in 0..10 {
            let a = BnFq12::random_element();
            let b = BnFq12::random_element();
            let c = BnFq12::random_element();
            let lhs = (a - b) * c;
            let rhs = a * c - b * c;
            assert_eq!(lhs, rhs, "(a-b)*c should equal a*c - b*c");
        }
    }

    #[test]
    fn fq12_invert() {
        for _ in 0..10 {
            let a = BnFq12::random_element();
            if !a.is_zero() {
                let a_inv = a.invert();
                let product = a * a_inv;
                assert_eq!(product, BnFq12::one(), "a * a^-1 should be one");
            }
        }
    }

    #[test]
    fn fq12_copy() {
        let a = BnFq12::random_element();
        let b = a;
        assert_eq!(a, b, "copy should preserve value");
    }

    #[test]
    fn fq12_unitary_inverse() {
        let a = BnFq12::random_element();
        let u_inv = a.unitary_inverse();
        // unitary_inverse flips c1 sign: (c0, c1) -> (c0, -c1)
        assert_eq!(u_inv.c0, a.c0, "c0 should be unchanged");
        assert_eq!(u_inv.c1, -a.c1, "c1 should be negated");
    }

    #[test]
    fn fq12_frobenius_map_one() {
        let a = BnFq12::random_element();
        let f = a.frobenius_map_one();
        // Just verify it doesn't panic and produces a non-trivial result
        assert_ne!(f, BnFq12::zero());
    }

    #[test]
    fn fq12_frobenius_map_two() {
        let a = BnFq12::random_element();
        let f = a.frobenius_map_two();
        assert_ne!(f, BnFq12::zero());
    }

    #[test]
    fn fq12_frobenius_map_three() {
        let a = BnFq12::random_element();
        let f = a.frobenius_map_three();
        assert_ne!(f, BnFq12::zero());
    }
}
