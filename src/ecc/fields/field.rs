use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::field_params::FieldParams;

// ---------------------------------------------------------------------------
// Helper functions (from field_impl_generic.hpp, __SIZEOF_INT128__ path)
// ---------------------------------------------------------------------------

/// 64x64 -> 128-bit wide multiply, returns (lo, hi).
#[inline(always)]
const fn mul_wide(a: u64, b: u64) -> (u64, u64) {
    let res = a as u128 * b as u128;
    (res as u64, (res >> 64) as u64)
}

/// Multiply-accumulate: a + b*c + carry_in -> (result, carry_out).
#[inline(always)]
const fn mac(a: u64, b: u64, c: u64, carry_in: u64) -> (u64, u64) {
    let res = a as u128 + (b as u128 * c as u128) + carry_in as u128;
    (res as u64, (res >> 64) as u64)
}

/// Multiply-accumulate without carry_in: a + b*c -> (result, carry_out).
#[inline(always)]
const fn mac_mini(a: u64, b: u64, c: u64) -> (u64, u64) {
    let res = a as u128 + (b as u128 * c as u128);
    (res as u64, (res >> 64) as u64)
}

/// Multiply-accumulate, discard low 64 bits: returns only high word of a + b*c.
#[inline(always)]
const fn mac_discard_lo(a: u64, b: u64, c: u64) -> u64 {
    let res = a as u128 + (b as u128 * c as u128);
    (res >> 64) as u64
}

/// Add with carry: a + b + carry_in -> (result, carry_out).
#[inline(always)]
const fn addc(a: u64, b: u64, carry_in: u64) -> (u64, u64) {
    let res = a as u128 + b as u128 + carry_in as u128;
    (res as u64, (res >> 64) as u64)
}

/// Subtract with borrow: a - b - (borrow_in >> 63) -> (result, borrow_out).
/// borrow_out has its MSB set (i.e. is all-ones or zero in the high bits) to indicate underflow.
/// Uses wrapping arithmetic to match C++ unsigned overflow semantics.
#[inline(always)]
const fn sbb(a: u64, b: u64, borrow_in: u64) -> (u64, u64) {
    let res = (a as u128).wrapping_sub(b as u128 + (borrow_in >> 63) as u128);
    (res as u64, (res >> 64) as u64)
}

/// Square-accumulate helper for optimized squaring.
/// Computes: out = a + 2*b*c + carry_in_lo, with carries propagated through carry_in_hi.
#[inline(always)]
const fn square_accumulate(
    a: u64,
    b: u64,
    c: u64,
    carry_in_lo: u64,
    carry_in_hi: u64,
) -> (u64, u64, u64) {
    let product = b as u128 * c as u128;
    let r0 = product as u64;
    let r1 = (product >> 64) as u64;

    let mut out = r0.wrapping_add(r0);
    let mut carry_lo: u64 = if out < r0 { 1 } else { 0 };
    out = out.wrapping_add(a);
    carry_lo += if out < a { 1 } else { 0 };
    out = out.wrapping_add(carry_in_lo);
    carry_lo += if out < carry_in_lo { 1 } else { 0 };
    carry_lo = carry_lo.wrapping_add(r1);
    let mut carry_hi: u64 = if carry_lo < r1 { 1 } else { 0 };
    carry_lo = carry_lo.wrapping_add(r1);
    carry_hi += if carry_lo < r1 { 1 } else { 0 };
    carry_lo = carry_lo.wrapping_add(carry_in_hi);
    carry_hi += if carry_lo < carry_in_hi { 1 } else { 0 };
    (out, carry_lo, carry_hi)
}

// ---------------------------------------------------------------------------
// Field<P> struct
// ---------------------------------------------------------------------------

/// A prime field element in Montgomery form, generic over parameters `P`.
///
/// Internally stores 4 x u64 limbs (little-endian). Values are kept in
/// the range [0, 2p) ("coarse" form) for most operations; full reduction
/// to [0, p) happens only on `reduce()`, `from_montgomery_form()`, and comparisons.
#[repr(C, align(32))]
pub struct Field<P: FieldParams> {
    pub data: [u64; 4],
    _phantom: PhantomData<P>,
}

// Manual Clone/Copy because PhantomData<P> doesn't require P: Copy
impl<P: FieldParams> Clone for Field<P> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            data: self.data,
            _phantom: PhantomData,
        }
    }
}

impl<P: FieldParams> Copy for Field<P> {}

impl<P: FieldParams> std::fmt::Debug for Field<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let r = self.reduce();
        write!(
            f,
            "Field(0x{:016x}{:016x}{:016x}{:016x})",
            r.data[3], r.data[2], r.data[1], r.data[0]
        )
    }
}

// Derived constants computed from FieldParams
impl<P: FieldParams> Field<P> {
    const MODULUS: [u64; 4] = P::MODULUS;

    /// -modulus mod 2^256 = 2^256 - modulus, used for branchless reduction.
    /// This is two's complement (!modulus + 1), NOT bitwise complement.
    const NOT_MODULUS: [u64; 4] = {
        let m = P::MODULUS;
        // Two's complement: -m = !m + 1
        let r0 = (!m[0]).wrapping_add(1);
        let c0 = (r0 < 1) as u64;
        let r1 = (!m[1]).wrapping_add(c0);
        let c1 = (r1 < c0) as u64;
        let r2 = (!m[2]).wrapping_add(c1);
        let c2 = (r2 < c1) as u64;
        let r3 = (!m[3]).wrapping_add(c2);
        [r0, r1, r2, r3]
    };

    /// 2 * modulus
    const TWICE_MODULUS: [u64; 4] = {
        let m = P::MODULUS;
        // Safe because modulus < 2^255 for all our curves
        let (r0, c) = (m[0] << 1, m[0] >> 63);
        let (r1, c) = ((m[1] << 1) | c, m[1] >> 63);
        let (r2, c) = ((m[2] << 1) | c, m[2] >> 63);
        let r3 = (m[3] << 1) | c;
        [r0, r1, r2, r3]
    };

    /// -(2*modulus) mod 2^256 = 2^256 - 2*modulus. Two's complement of twice_modulus.
    const TWICE_NOT_MODULUS: [u64; 4] = {
        let tm = Self::TWICE_MODULUS;
        let r0 = (!tm[0]).wrapping_add(1);
        let c0 = (r0 < 1) as u64;
        let r1 = (!tm[1]).wrapping_add(c0);
        let c1 = (r1 < c0) as u64;
        let r2 = (!tm[2]).wrapping_add(c1);
        let c2 = (r2 < c1) as u64;
        let r3 = (!tm[3]).wrapping_add(c2);
        [r0, r1, r2, r3]
    };

    /// modulus - 2, used for Fermat inversion
    const MODULUS_MINUS_TWO: [u64; 4] = {
        let m = P::MODULUS;
        [m[0].wrapping_sub(2), m[1], m[2], m[3]]
    };
}

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

impl<P: FieldParams> Field<P> {
    /// Zero element (additive identity). Not in Montgomery form since 0*R = 0.
    #[inline]
    pub const fn zero() -> Self {
        Self {
            data: [0, 0, 0, 0],
            _phantom: PhantomData,
        }
    }

    /// One element (multiplicative identity) in Montgomery form.
    #[inline]
    pub fn one() -> Self {
        Self::from(1u64)
    }

    /// Construct from a u64, auto-converting to Montgomery form.
    /// Matches C++ `field(uint64_t)`.
    #[inline]
    pub fn from(val: u64) -> Self {
        let mut f = Self {
            data: [val, 0, 0, 0],
            _phantom: PhantomData,
        };
        f = f.to_montgomery_form();
        f
    }

    /// Construct from raw limbs already in Montgomery form.
    /// Matches C++ `field{a, b, c, d}` brace-init.
    #[inline]
    pub const fn from_raw(data: [u64; 4]) -> Self {
        Self {
            data,
            _phantom: PhantomData,
        }
    }

    /// Construct from non-Montgomery limbs (standard integer representation)
    /// and convert to Montgomery form.
    #[inline]
    pub fn from_limbs(data: [u64; 4]) -> Self {
        let f = Self {
            data,
            _phantom: PhantomData,
        };
        f.to_montgomery_form()
    }
}

// ---------------------------------------------------------------------------
// Core arithmetic
// ---------------------------------------------------------------------------

impl<P: FieldParams> Field<P> {
    /// Full reduction from [0, 2p) to [0, p).
    #[inline]
    pub fn reduce(&self) -> Self {
        if P::MODULUS_IS_BIG {
            // Compare with modulus directly
            let ge = self.ge_modulus();
            if ge {
                let mut borrow = 0u64;
                let (r0, b) = sbb(self.data[0], Self::MODULUS[0], borrow);
                borrow = b;
                let (r1, b) = sbb(self.data[1], Self::MODULUS[1], borrow);
                borrow = b;
                let (r2, b) = sbb(self.data[2], Self::MODULUS[2], borrow);
                borrow = b;
                let (r3, _) = sbb(self.data[3], Self::MODULUS[3], borrow);
                return Self::from_raw([r0, r1, r2, r3]);
            }
            *self
        } else {
            // Branchless: add not_modulus, check carry
            let t0 = self.data[0].wrapping_add(Self::NOT_MODULUS[0]);
            let c = if t0 < self.data[0] { 1u64 } else { 0 };
            let (t1, c) = addc(self.data[1], Self::NOT_MODULUS[1], c);
            let (t2, c) = addc(self.data[2], Self::NOT_MODULUS[2], c);
            let (t3, c) = addc(self.data[3], Self::NOT_MODULUS[3], c);
            let mask = 0u64.wrapping_sub(c);
            let inv_mask = !mask;
            Self::from_raw([
                (self.data[0] & inv_mask) | (t0 & mask),
                (self.data[1] & inv_mask) | (t1 & mask),
                (self.data[2] & inv_mask) | (t2 & mask),
                (self.data[3] & inv_mask) | (t3 & mask),
            ])
        }
    }

    /// Check if self >= modulus (used for big modulus reduction).
    #[inline]
    fn ge_modulus(&self) -> bool {
        if self.data[3] > Self::MODULUS[3] {
            return true;
        }
        if self.data[3] < Self::MODULUS[3] {
            return false;
        }
        if self.data[2] > Self::MODULUS[2] {
            return true;
        }
        if self.data[2] < Self::MODULUS[2] {
            return false;
        }
        if self.data[1] > Self::MODULUS[1] {
            return true;
        }
        if self.data[1] < Self::MODULUS[1] {
            return false;
        }
        self.data[0] >= Self::MODULUS[0]
    }

    /// Modular addition. Result stays in [0, 2p).
    #[inline]
    pub fn add(&self, other: &Self) -> Self {
        if P::MODULUS_IS_BIG {
            // Add limbs
            let r0 = self.data[0].wrapping_add(other.data[0]);
            let c = if r0 < self.data[0] { 1u64 } else { 0 };
            let (r1, c) = addc(self.data[1], other.data[1], c);
            let (r2, c) = addc(self.data[2], other.data[2], c);
            let (r3, c) = addc(self.data[3], other.data[3], c);

            if c != 0 {
                // Overflow: subtract modulus (possibly twice)
                let mut borrow = 0u64;
                let (mut r0, b) = sbb(r0, Self::MODULUS[0], borrow);
                borrow = b;
                let (mut r1, b) = sbb(r1, Self::MODULUS[1], borrow);
                borrow = b;
                let (mut r2, b) = sbb(r2, Self::MODULUS[2], borrow);
                borrow = b;
                let (mut r3, b) = sbb(r3, Self::MODULUS[3], borrow);

                // If no borrow (still >= p), subtract again
                if b == 0 {
                    borrow = 0;
                    let (s0, b2) = sbb(r0, Self::MODULUS[0], borrow);
                    borrow = b2;
                    let (s1, b2) = sbb(r1, Self::MODULUS[1], borrow);
                    borrow = b2;
                    let (s2, b2) = sbb(r2, Self::MODULUS[2], borrow);
                    borrow = b2;
                    let (s3, _) = sbb(r3, Self::MODULUS[3], borrow);
                    r0 = s0;
                    r1 = s1;
                    r2 = s2;
                    r3 = s3;
                }
                return Self::from_raw([r0, r1, r2, r3]);
            }
            Self::from_raw([r0, r1, r2, r3])
        } else {
            // Add limbs
            let r0 = self.data[0].wrapping_add(other.data[0]);
            let c = if r0 < self.data[0] { 1u64 } else { 0 };
            let (r1, c) = addc(self.data[1], other.data[1], c);
            let (r2, c) = addc(self.data[2], other.data[2], c);
            let r3 = self.data[3].wrapping_add(other.data[3]).wrapping_add(c);

            // Branchless: add twice_not_modulus, use carry to select
            let t0 = r0.wrapping_add(Self::TWICE_NOT_MODULUS[0]);
            let c2 = if t0 < Self::TWICE_NOT_MODULUS[0] { 1u64 } else { 0 };
            let (t1, c2) = addc(r1, Self::TWICE_NOT_MODULUS[1], c2);
            let (t2, c2) = addc(r2, Self::TWICE_NOT_MODULUS[2], c2);
            let (t3, c2) = addc(r3, Self::TWICE_NOT_MODULUS[3], c2);
            let mask = 0u64.wrapping_sub(c2);
            let inv_mask = !mask;
            Self::from_raw([
                (r0 & inv_mask) | (t0 & mask),
                (r1 & inv_mask) | (t1 & mask),
                (r2 & inv_mask) | (t2 & mask),
                (r3 & inv_mask) | (t3 & mask),
            ])
        }
    }

    /// Modular subtraction. Result in [0, 2p).
    #[inline]
    pub fn subtract(&self, other: &Self) -> Self {
        let mut borrow = 0u64;
        let (mut r0, b) = sbb(self.data[0], other.data[0], borrow);
        borrow = b;
        let (mut r1, b) = sbb(self.data[1], other.data[1], borrow);
        borrow = b;
        let (mut r2, b) = sbb(self.data[2], other.data[2], borrow);
        borrow = b;
        let (mut r3, b) = sbb(self.data[3], other.data[3], borrow);
        borrow = b;

        // If borrow, add modulus (first correction)
        r0 = r0.wrapping_add(Self::MODULUS[0] & borrow);
        let mut carry = if r0 < (Self::MODULUS[0] & borrow) { 1u64 } else { 0 };
        let (v1, c) = addc(r1, Self::MODULUS[1] & borrow, carry);
        r1 = v1;
        carry = c;
        let (v2, c) = addc(r2, Self::MODULUS[2] & borrow, carry);
        r2 = v2;
        carry = c;
        // Track carry out of r3 to detect if we need a second add
        let r3_wide = r3 as u128 + (Self::MODULUS[3] & borrow) as u128 + carry as u128;
        r3 = r3_wide as u64;
        let carry_out = (r3_wide >> 64) as u64;

        // If still negative (first add didn't produce 256-bit overflow), add modulus again
        if carry_out == 0 && borrow != 0 {
            let old_r0 = r0;
            r0 = r0.wrapping_add(Self::MODULUS[0] & borrow);
            carry = if r0 < old_r0 { 1 } else { 0 };
            let (v1, c) = addc(r1, Self::MODULUS[1] & borrow, carry);
            r1 = v1;
            carry = c;
            let (v2, c) = addc(r2, Self::MODULUS[2] & borrow, carry);
            r2 = v2;
            carry = c;
            r3 = r3.wrapping_add((Self::MODULUS[3] & borrow).wrapping_add(carry));
        }

        Self::from_raw([r0, r1, r2, r3])
    }

    /// Coarse subtraction: on underflow adds 2p (stays in [0, 2p)).
    /// For big modulus, falls through to regular subtract.
    #[inline]
    pub fn subtract_coarse(&self, other: &Self) -> Self {
        if P::MODULUS_IS_BIG {
            return self.subtract(other);
        }
        let mut borrow = 0u64;
        let (r0, b) = sbb(self.data[0], other.data[0], borrow);
        borrow = b;
        let (r1, b) = sbb(self.data[1], other.data[1], borrow);
        borrow = b;
        let (r2, b) = sbb(self.data[2], other.data[2], borrow);
        borrow = b;
        let (r3, b) = sbb(self.data[3], other.data[3], borrow);
        borrow = b;

        let out0 = r0.wrapping_add(Self::TWICE_MODULUS[0] & borrow);
        let carry = if out0 < (Self::TWICE_MODULUS[0] & borrow) {
            1u64
        } else {
            0
        };
        let (out1, carry) = addc(r1, Self::TWICE_MODULUS[1] & borrow, carry);
        let (out2, carry) = addc(r2, Self::TWICE_MODULUS[2] & borrow, carry);
        let out3 = r3
            .wrapping_add(Self::TWICE_MODULUS[3] & borrow)
            .wrapping_add(carry);

        Self::from_raw([out0, out1, out2, out3])
    }

    /// Montgomery multiplication for big modulus (>= 2^254).
    #[inline]
    fn montgomery_mul_big(&self, other: &Self) -> Self {
        let modulus = Self::MODULUS;
        let r_inv = P::R_INV;

        let mut c: u64;
        let mut t0: u64 = 0;
        let mut t1: u64 = 0;
        let mut t2: u64 = 0;
        let mut t3: u64 = 0;
        let mut t4: u64 = 0;
        let mut t5: u64;

        for &element in &self.data {
            c = 0;
            let (v, co) = mac(t0, element, other.data[0], c);
            t0 = v;
            c = co;
            let (v, co) = mac(t1, element, other.data[1], c);
            t1 = v;
            c = co;
            let (v, co) = mac(t2, element, other.data[2], c);
            t2 = v;
            c = co;
            let (v, co) = mac(t3, element, other.data[3], c);
            t3 = v;
            c = co;
            let (v, ts) = addc(t4, c, 0);
            t4 = v;
            t5 = ts;

            let k = t0.wrapping_mul(r_inv);
            c = mac_discard_lo(t0, k, modulus[0]);
            let (v, co) = mac(t1, k, modulus[1], c);
            t0 = v;
            c = co;
            let (v, co) = mac(t2, k, modulus[2], c);
            t1 = v;
            c = co;
            let (v, co) = mac(t3, k, modulus[3], c);
            t2 = v;
            c = co;
            let (v, co2) = addc(c, t4, 0);
            t3 = v;
            t4 = t5 + co2;
        }

        // Final reduction: subtract modulus and conditionally add back
        let mut borrow = 0u64;
        let (r0, b) = sbb(t0, modulus[0], borrow);
        borrow = b;
        let (r1, b) = sbb(t1, modulus[1], borrow);
        borrow = b;
        let (r2, b) = sbb(t2, modulus[2], borrow);
        borrow = b;
        let (r3, b) = sbb(t3, modulus[3], borrow);
        borrow = b;

        // borrow ^ (0 - t4): if t4 != 0, we need to add back
        let borrow = borrow ^ (0u64.wrapping_sub(t4));

        let out0 = r0.wrapping_add(modulus[0] & borrow);
        let carry = if out0 < (modulus[0] & borrow) {
            1u64
        } else {
            0
        };
        let (out1, carry) = addc(r1, modulus[1] & borrow, carry);
        let (out2, carry) = addc(r2, modulus[2] & borrow, carry);
        let out3 = r3
            .wrapping_add(modulus[3] & borrow)
            .wrapping_add(carry);

        Self::from_raw([out0, out1, out2, out3])
    }

    /// Montgomery multiplication (small modulus, < 2^254).
    /// Unrolled 4-round interleaved multiply-reduce.
    #[inline]
    fn montgomery_mul_small(&self, other: &Self) -> Self {
        let modulus = Self::MODULUS;
        let r_inv = P::R_INV;

        // Round 0
        let (t0, c) = mul_wide(self.data[0], other.data[0]);
        let k = t0.wrapping_mul(r_inv);
        let a = mac_discard_lo(t0, k, modulus[0]);

        let (t1, a2) = mac_mini(a, self.data[0], other.data[1]);
        let (t0, c) = mac(t1, k, modulus[1], c);
        let (t2, a2b) = mac_mini(a2, self.data[0], other.data[2]);
        let (t1, c) = mac(t2, k, modulus[2], c);
        let (t3, a3) = mac_mini(a2b, self.data[0], other.data[3]);
        let (t2, c) = mac(t3, k, modulus[3], c);
        let t3 = c.wrapping_add(a3);

        // Round 1
        let (t0_new, a) = mac_mini(t0, self.data[1], other.data[0]);
        let k = t0_new.wrapping_mul(r_inv);
        let c = mac_discard_lo(t0_new, k, modulus[0]);
        let (t1_tmp, a2) = mac(t1, self.data[1], other.data[1], a);
        let (t0, c) = mac(t1_tmp, k, modulus[1], c);
        let (t2_tmp, a2b) = mac(t2, self.data[1], other.data[2], a2);
        let (t1, c) = mac(t2_tmp, k, modulus[2], c);
        let (t3_tmp, a3) = mac(t3, self.data[1], other.data[3], a2b);
        let (t2, c) = mac(t3_tmp, k, modulus[3], c);
        let t3 = c.wrapping_add(a3);

        // Round 2
        let (t0_new, a) = mac_mini(t0, self.data[2], other.data[0]);
        let k = t0_new.wrapping_mul(r_inv);
        let c = mac_discard_lo(t0_new, k, modulus[0]);
        let (t1_tmp, a2) = mac(t1, self.data[2], other.data[1], a);
        let (t0, c) = mac(t1_tmp, k, modulus[1], c);
        let (t2_tmp, a2b) = mac(t2, self.data[2], other.data[2], a2);
        let (t1, c) = mac(t2_tmp, k, modulus[2], c);
        let (t3_tmp, a3) = mac(t3, self.data[2], other.data[3], a2b);
        let (t2, c) = mac(t3_tmp, k, modulus[3], c);
        let t3 = c.wrapping_add(a3);

        // Round 3
        let (t0_new, a) = mac_mini(t0, self.data[3], other.data[0]);
        let k = t0_new.wrapping_mul(r_inv);
        let c = mac_discard_lo(t0_new, k, modulus[0]);
        let (t1_tmp, a2) = mac(t1, self.data[3], other.data[1], a);
        let (t0, c) = mac(t1_tmp, k, modulus[1], c);
        let (t2_tmp, a2b) = mac(t2, self.data[3], other.data[2], a2);
        let (t1, c) = mac(t2_tmp, k, modulus[2], c);
        let (t3_tmp, a3) = mac(t3, self.data[3], other.data[3], a2b);
        let (t2, c) = mac(t3_tmp, k, modulus[3], c);
        let t3 = c.wrapping_add(a3);

        Self::from_raw([t0, t1, t2, t3])
    }

    /// Montgomery multiplication dispatching to big or small path.
    #[inline]
    pub fn montgomery_mul(&self, other: &Self) -> Self {
        if P::MODULUS_IS_BIG {
            self.montgomery_mul_big(other)
        } else {
            self.montgomery_mul_small(other)
        }
    }

    /// Montgomery squaring (optimized for small modulus, delegates to mul_big for big).
    #[inline]
    pub fn montgomery_square(&self) -> Self {
        if P::MODULUS_IS_BIG {
            return self.montgomery_mul_big(self);
        }

        let modulus = Self::MODULUS;
        let r_inv = P::R_INV;

        // Round 0: diagonal + off-diagonal cross terms with data[0]
        let (t0, carry_lo) = mul_wide(self.data[0], self.data[0]);
        let carry_hi = 0u64;
        let (t1, carry_lo, carry_hi) =
            square_accumulate(0, self.data[1], self.data[0], carry_lo, carry_hi);
        let (t2, carry_lo, carry_hi) =
            square_accumulate(0, self.data[2], self.data[0], carry_lo, carry_hi);
        let (t3, carry_lo, _carry_hi) =
            square_accumulate(0, self.data[3], self.data[0], carry_lo, carry_hi);

        let round_carry = carry_lo;
        let k = t0.wrapping_mul(r_inv);
        let c = mac_discard_lo(t0, k, modulus[0]);
        let (t0, c) = mac(t1, k, modulus[1], c);
        let (t1, c) = mac(t2, k, modulus[2], c);
        let (t2, c) = mac(t3, k, modulus[3], c);
        let t3 = c.wrapping_add(round_carry);

        // Round 1
        let (t1_new, carry_lo) = mac_mini(t1, self.data[1], self.data[1]);
        let carry_hi = 0u64;
        let (t2_new, carry_lo, carry_hi) =
            square_accumulate(t2, self.data[2], self.data[1], carry_lo, carry_hi);
        let (t3_new, carry_lo, _carry_hi) =
            square_accumulate(t3, self.data[3], self.data[1], carry_lo, carry_hi);
        let round_carry = carry_lo;
        let k = t0.wrapping_mul(r_inv);
        let c = mac_discard_lo(t0, k, modulus[0]);
        let (t0, c) = mac(t1_new, k, modulus[1], c);
        let (t1, c) = mac(t2_new, k, modulus[2], c);
        let (t2, c) = mac(t3_new, k, modulus[3], c);
        let t3 = c.wrapping_add(round_carry);

        // Round 2
        let (t2_new, carry_lo) = mac_mini(t2, self.data[2], self.data[2]);
        let carry_hi = 0u64;
        let (t3_new, carry_lo, _carry_hi) =
            square_accumulate(t3, self.data[3], self.data[2], carry_lo, carry_hi);
        let round_carry = carry_lo;
        let k = t0.wrapping_mul(r_inv);
        let c = mac_discard_lo(t0, k, modulus[0]);
        let (t0, c) = mac(t1, k, modulus[1], c);
        let (t1, c) = mac(t2_new, k, modulus[2], c);
        let (t2, c) = mac(t3_new, k, modulus[3], c);
        let t3 = c.wrapping_add(round_carry);

        // Round 3
        let (t3_new, carry_lo) = mac_mini(t3, self.data[3], self.data[3]);
        let k = t0.wrapping_mul(r_inv);
        let round_carry = carry_lo;
        let c = mac_discard_lo(t0, k, modulus[0]);
        let (t0, c) = mac(t1, k, modulus[1], c);
        let (t1, c) = mac(t2, k, modulus[2], c);
        let (t2, c) = mac(t3_new, k, modulus[3], c);
        let t3 = c.wrapping_add(round_carry);

        Self::from_raw([t0, t1, t2, t3])
    }

    /// Alias for montgomery_square.
    #[inline]
    pub fn sqr(&self) -> Self {
        self.montgomery_square()
    }

    /// Convert from standard form to Montgomery form: self * R^2 mod p.
    #[inline]
    pub fn to_montgomery_form(&self) -> Self {
        let r_squared = Self::from_raw(P::R_SQUARED);
        // C++ does self_reduce_once() x3 before multiplying
        let mut tmp = self.reduce();
        tmp = tmp.reduce();
        tmp = tmp.reduce();
        let result = tmp.montgomery_mul(&r_squared);
        result.reduce()
    }

    /// Convert from Montgomery form to standard form: self * 1 mod p.
    #[inline]
    pub fn from_montgomery_form(&self) -> Self {
        let one_raw = Self::from_raw([1, 0, 0, 0]);
        let result = self.montgomery_mul(&one_raw);
        result.reduce()
    }

    /// Negate: returns -self mod p.
    #[inline]
    pub fn negate(&self) -> Self {
        if P::MODULUS_IS_BIG {
            let p = Self::from_raw(Self::MODULUS);
            p.subtract(self)
        } else {
            let p2 = Self::from_raw(Self::TWICE_MODULUS);
            p2.subtract_coarse(self).reduce()
        }
    }

    /// Check if zero (either 0 or p in limb representation).
    #[inline]
    pub fn is_zero(&self) -> bool {
        ((self.data[0] | self.data[1] | self.data[2] | self.data[3]) == 0)
            || (self.data[0] == P::MODULUS[0]
                && self.data[1] == P::MODULUS[1]
                && self.data[2] == P::MODULUS[2]
                && self.data[3] == P::MODULUS[3])
    }

    /// Equality: reduce both, compare limbs.
    #[inline]
    pub fn eq_field(&self, other: &Self) -> bool {
        let a = self.reduce();
        let b = other.reduce();
        a.data[0] == b.data[0]
            && a.data[1] == b.data[1]
            && a.data[2] == b.data[2]
            && a.data[3] == b.data[3]
    }

    /// Exponentiation via square-and-multiply.
    pub fn pow(&self, exp: &[u64; 4]) -> Self {
        // Find MSB
        let mut msb = 0u32;
        for i in (0..4).rev() {
            if exp[i] != 0 {
                msb = (i as u32) * 64 + (63 - exp[i].leading_zeros());
                break;
            }
        }
        // Check if exponent is zero
        if exp[0] == 0 && exp[1] == 0 && exp[2] == 0 && exp[3] == 0 {
            return Self::one();
        }
        if self.is_zero() {
            return Self::zero();
        }

        let mut accumulator = *self;
        let to_mul = *self;
        for i in (0..msb).rev() {
            accumulator = accumulator.sqr();
            let limb_idx = (i / 64) as usize;
            let bit_idx = i % 64;
            if (exp[limb_idx] >> bit_idx) & 1 == 1 {
                accumulator = accumulator.montgomery_mul(&to_mul);
            }
        }
        accumulator
    }

    /// Modular inverse via Fermat's little theorem: self^(p-2) mod p.
    pub fn invert(&self) -> Self {
        debug_assert!(!self.is_zero(), "cannot invert zero");
        self.pow(&Self::MODULUS_MINUS_TWO)
    }

    /// Square root. Returns (true, root) if QR, (false, zero) otherwise.
    /// Uses (p+1)/4 exponentiation when p ≡ 3 (mod 4).
    pub fn sqrt(&self) -> (bool, Self) {
        if P::MODULUS[0] & 0x3 == 0x3 {
            // p ≡ 3 mod 4: sqrt = self^((p+1)/4)
            let exp = {
                let m = P::MODULUS;
                // (modulus + 1) >> 2
                let (a0, c) = addc(m[0], 1, 0);
                let (a1, c) = addc(m[1], 0, c);
                let (a2, c) = addc(m[2], 0, c);
                let a3 = m[3].wrapping_add(c);
                // >> 2
                [
                    (a0 >> 2) | (a1 << 62),
                    (a1 >> 2) | (a2 << 62),
                    (a2 >> 2) | (a3 << 62),
                    a3 >> 2,
                ]
            };
            let root = self.pow(&exp);
            let check = root.sqr();
            if check.eq_field(self) {
                (true, root)
            } else {
                (false, Self::zero())
            }
        } else {
            // Tonelli-Shanks would go here for other moduli
            // For now, not needed for our curves (BN254 Fr has modulus_0 & 3 == 1)
            // Fall back to generic Tonelli-Shanks
            self.tonelli_shanks_sqrt()
        }
    }

    /// Tonelli-Shanks square root (simplified version).
    fn tonelli_shanks_sqrt(&self) -> (bool, Self) {
        // Factor p-1 = Q * 2^S
        let mut q = {
            let m = P::MODULUS;
            // p - 1
            [m[0].wrapping_sub(1), m[1], m[2], m[3]]
        };
        let mut s = 0u32;
        while {
            let limb = (s / 64) as usize;
            let bit = s % 64;
            limb < 4 && (q[limb] >> bit) & 1 == 0
        } {
            s += 1;
        }
        // q = (p-1) >> s
        let shift = s;
        q = Self::shr_limbs(&q, shift);

        // Find a non-residue z (try 2, 3, 4, ...)
        let p_minus_1_over_2 = {
            let m = P::MODULUS;
            let pm1 = [m[0].wrapping_sub(1), m[1], m[2], m[3]];
            Self::shr_limbs(&pm1, 1)
        };
        let mut z_val = 2u64;
        let neg_one = Self::from(1u64).negate();
        let z = loop {
            let z_field = Self::from(z_val);
            let check = z_field.pow(&p_minus_1_over_2);
            if check.eq_field(&neg_one) {
                break z_field;
            }
            z_val += 1;
        };

        let mut m_val = s;
        let mut c = z.pow(&q);
        let mut t = self.pow(&q);
        // r = self^((q+1)/2)
        let q_plus_1_over_2 = {
            let (a0, c) = addc(q[0], 1, 0);
            let (a1, c) = addc(q[1], 0, c);
            let (a2, c) = addc(q[2], 0, c);
            let a3 = q[3].wrapping_add(c);
            Self::shr_limbs(&[a0, a1, a2, a3], 1)
        };
        let mut r = self.pow(&q_plus_1_over_2);

        loop {
            if t.eq_field(&Self::one()) {
                return (true, r);
            }
            if t.is_zero() {
                return (true, Self::zero());
            }
            // Find least i such that t^(2^i) = 1
            let mut i = 1u32;
            let mut tmp = t.sqr();
            while !tmp.eq_field(&Self::one()) {
                tmp = tmp.sqr();
                i += 1;
                if i >= m_val {
                    return (false, Self::zero());
                }
            }
            // b = c^(2^(m-i-1))
            let mut b = c;
            for _ in 0..(m_val - i - 1) {
                b = b.sqr();
            }
            m_val = i;
            c = b.sqr();
            t = t.montgomery_mul(&c);
            r = r.montgomery_mul(&b);
        }
    }

    /// Right-shift a 4-limb number by `shift` bits.
    fn shr_limbs(val: &[u64; 4], shift: u32) -> [u64; 4] {
        if shift == 0 {
            return *val;
        }
        if shift >= 256 {
            return [0; 4];
        }
        let limb_shift = (shift / 64) as usize;
        let bit_shift = shift % 64;
        let mut result = [0u64; 4];
        for i in 0..4 {
            let src = i + limb_shift;
            if src < 4 {
                result[i] = val[src] >> bit_shift;
                if bit_shift > 0 && src + 1 < 4 {
                    result[i] |= val[src + 1] << (64 - bit_shift);
                }
            }
        }
        result
    }

    // --- Infinity flag helpers (used by group operations) ---

    /// Set MSB of data[3] to flag infinity (small modulus convention).
    #[inline]
    pub fn self_set_msb(&mut self) {
        self.data[3] = 1u64 << 63;
    }

    /// Check if MSB of data[3] is set.
    #[inline]
    pub fn is_msb_set(&self) -> bool {
        (self.data[3] >> 63) == 1
    }

    /// Conditional negate (branchless via predicate).
    #[inline]
    pub fn conditional_negate(&self, predicate: bool) -> Self {
        if predicate {
            self.negate()
        } else {
            *self
        }
    }

    /// Cube root of unity from field params.
    #[inline]
    pub fn cube_root_of_unity() -> Self {
        Self::from_raw(P::CUBE_ROOT)
    }

    /// Serialize to 32 big-endian bytes (from Montgomery form).
    ///
    /// Converts from Montgomery form to standard integer representation,
    /// then writes data[3] (MSB) first as big-endian u64s.
    pub fn to_be_bytes(&self) -> [u8; 32] {
        let reduced = self.from_montgomery_form();
        let mut bytes = [0u8; 32];
        bytes[0..8].copy_from_slice(&reduced.data[3].to_be_bytes());
        bytes[8..16].copy_from_slice(&reduced.data[2].to_be_bytes());
        bytes[16..24].copy_from_slice(&reduced.data[1].to_be_bytes());
        bytes[24..32].copy_from_slice(&reduced.data[0].to_be_bytes());
        bytes
    }

    /// Deserialize from 32 big-endian bytes (into Montgomery form).
    ///
    /// Reads 4 big-endian u64s, constructs limbs, and converts to Montgomery form.
    /// Performs modular reduction if value >= modulus.
    pub fn from_be_bytes(bytes: &[u8; 32]) -> Self {
        let data3 = u64::from_be_bytes(bytes[0..8].try_into().unwrap());
        let data2 = u64::from_be_bytes(bytes[8..16].try_into().unwrap());
        let data1 = u64::from_be_bytes(bytes[16..24].try_into().unwrap());
        let data0 = u64::from_be_bytes(bytes[24..32].try_into().unwrap());
        Self::from_limbs([data0, data1, data2, data3])
    }

    /// Extract bit at position `idx` from the reduced (non-Montgomery) form.
    pub fn get_bit(&self, idx: usize) -> bool {
        let reduced = self.from_montgomery_form();
        let limb_idx = idx / 64;
        let bit_idx = idx % 64;
        if limb_idx >= 4 {
            return false;
        }
        (reduced.data[limb_idx] >> bit_idx) & 1 == 1
    }

    /// Generate a uniformly random field element.
    ///
    /// Generates 512 random bits and reduces mod p for uniform distribution.
    pub fn random_element() -> Self {
        use rand::Rng;
        let mut rng = rand::rng();
        let lo = [
            rng.random::<u64>(),
            rng.random::<u64>(),
            rng.random::<u64>(),
            rng.random::<u64>(),
        ];
        let hi = [
            rng.random::<u64>(),
            rng.random::<u64>(),
            rng.random::<u64>(),
            rng.random::<u64>(),
        ];
        Self::from_u512(lo, hi)
    }

    /// Reduce a 512-bit value (lo || hi) modulo the field modulus.
    ///
    /// Used by hash-to-curve to map a 512-bit hash output into a field element,
    /// ensuring negligible bias.
    pub fn from_u512(lo: [u64; 4], hi: [u64; 4]) -> Self {
        use crypto_bigint::{NonZero, U256, U512};

        // Build the 512-bit value: val = (hi << 256) | lo
        // Tuple order: (lower_half, upper_half)
        let lo_256 = U256::from_words(lo);
        let hi_256 = U256::from_words(hi);
        let val = U512::from((lo_256, hi_256));

        // Build modulus as U512 (modulus in lower half, zero in upper)
        let modulus = U256::from_words(P::MODULUS);
        let modulus_wide = U512::from((modulus, U256::ZERO));
        let nz_mod = NonZero::new(modulus_wide).expect("modulus is nonzero");

        // Compute val % modulus
        let (_, remainder) = val.div_rem(&nz_mod);

        // Extract low 256 bits (remainder < modulus < 2^256, so upper half is zero)
        let words: [u64; 8] = remainder.to_words();
        let limbs = [words[0], words[1], words[2], words[3]];

        Self::from_limbs(limbs)
    }
}

// ---------------------------------------------------------------------------
// Operator impls
// ---------------------------------------------------------------------------

impl<P: FieldParams> Add for Field<P> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Field::add(&self, &rhs)
    }
}

impl<P: FieldParams> AddAssign for Field<P> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = Field::add(self, &rhs);
    }
}

impl<P: FieldParams> Sub for Field<P> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Field::subtract_coarse(&self, &rhs)
    }
}

impl<P: FieldParams> SubAssign for Field<P> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = Field::subtract_coarse(self, &rhs);
    }
}

impl<P: FieldParams> Mul for Field<P> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Field::montgomery_mul(&self, &rhs)
    }
}

impl<P: FieldParams> MulAssign for Field<P> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = Field::montgomery_mul(self, &rhs);
    }
}

impl<P: FieldParams> Neg for Field<P> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Field::negate(&self)
    }
}

impl<P: FieldParams> PartialEq for Field<P> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.eq_field(other)
    }
}

impl<P: FieldParams> Eq for Field<P> {}
