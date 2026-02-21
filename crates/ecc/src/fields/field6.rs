// Sextic extension field Fq6 = Fq2[v] / (v^3 - xi)
//
// C++ source: ecc/fields/field6.hpp
//
// Elements are triples (c0, c1, c2) of Field2 elements.
// xi is a non-residue in Fq2, provided by the Field6Params trait.

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::field2::Field2;
use super::field_params::FieldParams;

/// Trait providing the BN254-specific Fq6 parameters.
/// Extends FieldParams (the base prime field) with extension-field constants.
pub trait Field6Params: FieldParams {
    /// Multiply an Fq2 element by the non-residue xi in Fq2.
    fn mul_by_non_residue(a: &Field2<Self>) -> Field2<Self>;

    fn frobenius_coeffs_c1_1() -> Field2<Self>;
    fn frobenius_coeffs_c1_2() -> Field2<Self>;
    fn frobenius_coeffs_c1_3() -> Field2<Self>;
    fn frobenius_coeffs_c2_1() -> Field2<Self>;
    fn frobenius_coeffs_c2_2() -> Field2<Self>;
    fn frobenius_coeffs_c2_3() -> Field2<Self>;
}

pub struct Field6<P: Field6Params> {
    pub c0: Field2<P>,
    pub c1: Field2<P>,
    pub c2: Field2<P>,
}

impl<P: Field6Params> Clone for Field6<P> {
    #[inline]
    fn clone(&self) -> Self {
        Self { c0: self.c0, c1: self.c1, c2: self.c2 }
    }
}

impl<P: Field6Params> Copy for Field6<P> {}

impl<P: Field6Params> std::fmt::Debug for Field6<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Field6({:?}, {:?}, {:?})", self.c0, self.c1, self.c2)
    }
}

impl<P: Field6Params> Field6<P> {
    #[inline]
    pub fn new(c0: Field2<P>, c1: Field2<P>, c2: Field2<P>) -> Self {
        Self { c0, c1, c2 }
    }

    #[inline]
    pub fn zero() -> Self {
        Self {
            c0: Field2::zero(),
            c1: Field2::zero(),
            c2: Field2::zero(),
        }
    }

    #[inline]
    pub fn one() -> Self {
        Self {
            c0: Field2::one(),
            c1: Field2::zero(),
            c2: Field2::zero(),
        }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }

    /// Multiply an Fq2 element by the Fq6 non-residue (delegates to params).
    #[inline]
    pub fn mul_by_non_residue(a: &Field2<P>) -> Field2<P> {
        P::mul_by_non_residue(a)
    }

    /// Multiply each Fq2 component by an Fq2 scalar.
    #[inline]
    pub fn mul_by_fq2(&self, other: &Field2<P>) -> Self {
        Self {
            c0: *other * self.c0,
            c1: *other * self.c1,
            c2: *other * self.c2,
        }
    }

    /// Squaring via CH-SQR2 (Devegili et al.).
    pub fn sqr(&self) -> Self {
        let s0 = self.c0.sqr();
        let mut s1 = self.c0 * self.c1;
        s1 = s1 + s1;
        let s2 = (self.c0 + self.c2 - self.c1).sqr();
        let mut s3 = self.c1 * self.c2;
        s3 = s3 + s3;
        let s4 = self.c2.sqr();
        Self {
            c0: Self::mul_by_non_residue(&s3) + s0,
            c1: Self::mul_by_non_residue(&s4) + s1,
            c2: s1 + s2 + s3 - s0 - s4,
        }
    }

    /// Inversion via Algorithm 17 from "High-Speed Software Implementation of
    /// the Optimal Ate Pairing over Barreto-Naehrig Curves".
    pub fn invert(&self) -> Self {
        let c0_sq = self.c0.sqr();
        let c1_c2 = self.c1 * self.c2;
        let cap_c0 = c0_sq - Self::mul_by_non_residue(&c1_c2);

        let c2_sq = self.c2.sqr();
        let c0_c1 = self.c0 * self.c1;
        let cap_c1 = Self::mul_by_non_residue(&c2_sq) - c0_c1;

        let c1_sq = self.c1.sqr();
        let c0_c2 = self.c0 * self.c2;
        let cap_c2 = c1_sq - c0_c2;

        let t0 = (self.c0 * cap_c0
            + Self::mul_by_non_residue(&(self.c2 * cap_c1 + self.c1 * cap_c2)))
        .invert();

        Self {
            c0: t0 * cap_c0,
            c1: t0 * cap_c1,
            c2: t0 * cap_c2,
        }
    }

    pub fn frobenius_map_one(&self) -> Self {
        Self {
            c0: self.c0.frobenius_map(),
            c1: P::frobenius_coeffs_c1_1() * self.c1.frobenius_map(),
            c2: P::frobenius_coeffs_c2_1() * self.c2.frobenius_map(),
        }
    }

    pub fn frobenius_map_two(&self) -> Self {
        Self {
            c0: self.c0,
            c1: P::frobenius_coeffs_c1_2() * self.c1,
            c2: P::frobenius_coeffs_c2_2() * self.c2,
        }
    }

    pub fn frobenius_map_three(&self) -> Self {
        Self {
            c0: self.c0.frobenius_map(),
            c1: P::frobenius_coeffs_c1_3() * self.c1.frobenius_map(),
            c2: P::frobenius_coeffs_c2_3() * self.c2.frobenius_map(),
        }
    }

    #[inline]
    pub fn to_montgomery_form(&self) -> Self {
        Self {
            c0: self.c0.to_montgomery_form(),
            c1: self.c1.to_montgomery_form(),
            c2: self.c2.to_montgomery_form(),
        }
    }

    #[inline]
    pub fn from_montgomery_form(&self) -> Self {
        Self {
            c0: self.c0.from_montgomery_form(),
            c1: self.c1.from_montgomery_form(),
            c2: self.c2.from_montgomery_form(),
        }
    }
}

// ---------------------------------------------------------------------------
// Operator impls
// ---------------------------------------------------------------------------

impl<P: Field6Params> Add for Field6<P> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
            c2: self.c2 + rhs.c2,
        }
    }
}

impl<P: Field6Params> AddAssign for Field6<P> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<P: Field6Params> Sub for Field6<P> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
            c2: self.c2 - rhs.c2,
        }
    }
}

impl<P: Field6Params> SubAssign for Field6<P> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<P: Field6Params> Mul for Field6<P> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        // Karatsuba (Devegili et al., Section 4)
        let t0 = self.c0 * rhs.c0;
        let t1 = self.c1 * rhs.c1;
        let t2 = self.c2 * rhs.c2;

        let t3 = (self.c0 + self.c2) * (rhs.c0 + rhs.c2);
        let t4 = (self.c0 + self.c1) * (rhs.c0 + rhs.c1);
        let t5 = (self.c1 + self.c2) * (rhs.c1 + rhs.c2);

        Self {
            c0: t0 + Self::mul_by_non_residue(&(t5 - (t1 + t2))),
            c1: t4 - (t0 + t1) + Self::mul_by_non_residue(&t2),
            c2: t3 + t1 - (t0 + t2),
        }
    }
}

impl<P: Field6Params> MulAssign for Field6<P> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<P: Field6Params> Neg for Field6<P> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
        }
    }
}

impl<P: Field6Params> PartialEq for Field6<P> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1 && self.c2 == other.c2
    }
}

impl<P: Field6Params> Eq for Field6<P> {}
