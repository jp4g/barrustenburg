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
