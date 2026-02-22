// Extended-width unsigned integer types (512-bit, 1024-bit).
//
// C++ source: barretenberg/cpp/src/barretenberg/numeric/uintx/uintx.hpp
//
// BB's `uintx<T>` is a generic double-width type composed of (lo, hi) halves.
// We map this to crypto-bigint's wider Uint types:
//   uint512_t  → U512  = Uint<8>
//   uint1024_t → U1024 = Uint<16>
//
// crypto-bigint provides concat/split for converting between widths,
// which mirrors BB's `.lo` / `.hi` access pattern.

use crypto_bigint::Uint;

/// 512-bit unsigned integer.
pub type U512 = Uint<8>;

/// 1024-bit unsigned integer.
pub type U1024 = Uint<16>;

/// BB-compatible extension methods for U512.
///
/// Mirrors the `.lo` / `.hi` access pattern from BB's `uintx<uint256_t>`.
pub trait U512Ext {
    /// Extract the low 256 bits. Mirrors BB's `uint512_t.lo`.
    fn lo(&self) -> super::U256;

    /// Extract the high 256 bits. Mirrors BB's `uint512_t.hi`.
    fn hi(&self) -> super::U256;

    /// Construct from (lo, hi) pair. Mirrors BB's `uint512_t(lo, hi)`.
    fn from_lo_hi(lo: super::U256, hi: super::U256) -> Self;
}

impl U512Ext for U512 {
    fn lo(&self) -> super::U256 {
        let (lo, _hi) = self.split();
        lo
    }

    fn hi(&self) -> super::U256 {
        let (_lo, hi) = self.split();
        hi
    }

    fn from_lo_hi(lo: super::U256, hi: super::U256) -> Self {
        lo.concat(&hi)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::uint256::U256Ext;
    use crate::{U256, U1024};
    use crypto_bigint::Zero;

    #[test]
    fn lo_hi_roundtrip() {
        let lo = U256::from_limbs([1, 2, 3, 4]);
        let hi = U256::from_limbs([5, 6, 7, 8]);
        let wide = U512::from_lo_hi(lo, hi);
        assert_eq!(wide.lo(), lo);
        assert_eq!(wide.hi(), hi);
    }

    #[test]
    fn u256_widening_mul_to_u512() {
        let a = U256::from_limbs([u64::MAX, 0, 0, 0]);
        let b = U256::from_limbs([u64::MAX, 0, 0, 0]);
        let wide = a.widening_mul(&b);
        let (lo, hi) = wide.split();
        // (2^64 - 1)^2 = 2^128 - 2^65 + 1
        assert_eq!(lo, U256::from_limbs([1, 0xFFFF_FFFF_FFFF_FFFE, 0, 0]));
        assert_eq!(hi, U256::ZERO);
    }

    #[test]
    fn u512_mod_reduction() {
        // Test the pattern BB uses: (uint512_t(...) % modulus).lo
        let a = U256::from_limbs([100, 0, 0, 0]);
        let b = U256::from_limbs([7, 0, 0, 0]);
        let wide_a = U512::from_lo_hi(a, U256::ZERO);
        let wide_b = U512::from_lo_hi(b, U256::ZERO);
        let (_q, r) = wide_a.div_rem(&wide_b.to_nz().unwrap());
        let result = r.lo();
        assert_eq!(result, U256::from_limbs([2, 0, 0, 0]));
    }

    #[test]
    fn u512_add() {
        let a = U512::from_lo_hi(U256::from_limbs([100, 0, 0, 0]), U256::ZERO);
        let b = U512::from_lo_hi(U256::from_limbs([200, 0, 0, 0]), U256::ZERO);
        let c = a.wrapping_add(&b);
        assert_eq!(c.lo(), U256::from_limbs([300, 0, 0, 0]));
        assert_eq!(c.hi(), U256::ZERO);
    }

    #[test]
    fn u512_sub() {
        let a = U512::from_lo_hi(U256::from_limbs([300, 0, 0, 0]), U256::ZERO);
        let b = U512::from_lo_hi(U256::from_limbs([100, 0, 0, 0]), U256::ZERO);
        let c = a.wrapping_sub(&b);
        assert_eq!(c.lo(), U256::from_limbs([200, 0, 0, 0]));
    }

    #[test]
    fn u512_mul() {
        let a = U512::from_lo_hi(U256::from_limbs([7, 0, 0, 0]), U256::ZERO);
        let b = U512::from_lo_hi(U256::from_limbs([6, 0, 0, 0]), U256::ZERO);
        let c = a.wrapping_mul(&b);
        assert_eq!(c.lo(), U256::from_limbs([42, 0, 0, 0]));
    }

    #[test]
    fn u512_div_and_mod() {
        let a = U512::from_lo_hi(U256::from_limbs([100, 0, 0, 0]), U256::ZERO);
        let b = U512::from_lo_hi(U256::from_limbs([7, 0, 0, 0]), U256::ZERO);
        let (q, r) = a.div_rem(&b.to_nz().unwrap());
        assert_eq!(q.lo(), U256::from_limbs([14, 0, 0, 0]));
        assert_eq!(r.lo(), U256::from_limbs([2, 0, 0, 0]));
    }

    #[test]
    fn u512_get_bit() {
        let val = U512::from_lo_hi(U256::from_limbs([0b1010, 0, 0, 0]), U256::ZERO);
        assert!(val.bit_vartime(1));
        assert!(!val.bit_vartime(2));
        assert!(val.bit_vartime(3));
    }

    #[test]
    fn u512_and() {
        let a = U512::from_lo_hi(U256::from_limbs([0xFF00, 0, 0, 0]), U256::ZERO);
        let b = U512::from_lo_hi(U256::from_limbs([0x0FF0, 0, 0, 0]), U256::ZERO);
        let c = a.bitand(&b);
        assert_eq!(c.lo(), U256::from_limbs([0x0F00, 0, 0, 0]));
    }

    #[test]
    fn u512_or() {
        let a = U512::from_lo_hi(U256::from_limbs([0xFF00, 0, 0, 0]), U256::ZERO);
        let b = U512::from_lo_hi(U256::from_limbs([0x00FF, 0, 0, 0]), U256::ZERO);
        let c = a.bitor(&b);
        assert_eq!(c.lo(), U256::from_limbs([0xFFFF, 0, 0, 0]));
    }

    #[test]
    fn u512_xor() {
        let a = U512::from_lo_hi(U256::from_limbs([0xFF00, 0, 0, 0]), U256::ZERO);
        let b = U512::from_lo_hi(U256::from_limbs([0x0FF0, 0, 0, 0]), U256::ZERO);
        let c = a.bitxor(&b);
        assert_eq!(c.lo(), U256::from_limbs([0xF0F0, 0, 0, 0]));
    }

    #[test]
    fn u512_not_equal() {
        let a = U512::from_lo_hi(U256::from_limbs([42, 0, 0, 0]), U256::ZERO);
        let b = U512::from_lo_hi(U256::from_limbs([43, 0, 0, 0]), U256::ZERO);
        assert_ne!(a, b);
    }

    #[test]
    fn u512_equality() {
        let a = U512::from_lo_hi(U256::from_limbs([42, 0, 0, 0]), U256::ZERO);
        let b = U512::from_lo_hi(U256::from_limbs([42, 0, 0, 0]), U256::ZERO);
        assert_eq!(a, b);
    }

    #[test]
    fn u512_cross_boundary() {
        // Test addition that crosses the lo/hi boundary
        let a = U512::from_lo_hi(
            U256::from_limbs([u64::MAX, u64::MAX, u64::MAX, u64::MAX]),
            U256::ZERO,
        );
        let b = U512::from_lo_hi(U256::ONE, U256::ZERO);
        let c = a.wrapping_add(&b);
        assert_eq!(c.lo(), U256::ZERO);
        assert_eq!(c.hi(), U256::ONE);
    }

    #[test]
    fn u1024_mod_reduction() {
        // Create U1024 from two U512 halves, reduce mod known modulus
        let lo = U512::from_lo_hi(U256::from_limbs([100, 0, 0, 0]), U256::ZERO);
        let hi = U512::from_lo_hi(U256::ZERO, U256::ZERO);
        let wide: U1024 = lo.concat(&hi);
        let modulus_lo = U512::from_lo_hi(U256::from_limbs([7, 0, 0, 0]), U256::ZERO);
        let modulus: U1024 = modulus_lo.concat(&U512::from_lo_hi(U256::ZERO, U256::ZERO));
        let (_q, r) = wide.div_rem(&modulus.to_nz().unwrap());
        let (r_lo, _r_hi) = r.split();
        assert_eq!(r_lo.lo(), U256::from_limbs([2, 0, 0, 0]));
    }

    #[test]
    fn u512_bit_not() {
        let a = U512::from_lo_hi(U256::from_limbs([0xFF00, 0, 0, 0]), U256::ZERO);
        let not_a = a.not();
        let lo = not_a.lo();
        assert_eq!(lo.limbs()[0], !0xFF00u64);
        assert_eq!(lo.limbs()[1], u64::MAX);
        assert_eq!(not_a.hi(), U256::from_limbs([u64::MAX; 4]));
        assert_eq!(not_a.not(), a);
    }

    #[test]
    fn u512_logic_not() {
        assert!(bool::from(U512::ZERO.is_zero()));
        assert!(!bool::from(U512::ONE.is_zero()));
        let nonzero = U512::from_lo_hi(U256::ZERO, U256::ONE);
        assert!(!bool::from(nonzero.is_zero()));
    }

    #[test]
    fn u512_invmod_regression() {
        // Verify a * a_inv mod m == 1 using widening mul + div_rem
        let a = U512::from_lo_hi(U256::from_limbs([7, 0, 0, 0]), U256::ZERO);
        let m = U512::from_lo_hi(U256::from_limbs([11, 0, 0, 0]), U256::ZERO);
        // 7 * 8 = 56, 56 mod 11 = 1 (so 8 is the inverse of 7 mod 11)
        let a_inv = U512::from_lo_hi(U256::from_limbs([8, 0, 0, 0]), U256::ZERO);
        let product = a.wrapping_mul(&a_inv);
        let (_q, r) = product.div_rem(&m.to_nz().unwrap());
        assert_eq!(r, U512::from_lo_hi(U256::from_limbs([1, 0, 0, 0]), U256::ZERO));
    }

    #[test]
    fn u512_barrett_reduction_regression() {
        // Reduce specific 512-bit value mod known modulus, verify exact result
        // 1000000 mod 997 = 9
        let a = U512::from_lo_hi(U256::from_limbs([1_000_000, 0, 0, 0]), U256::ZERO);
        let m = U512::from_lo_hi(U256::from_limbs([997, 0, 0, 0]), U256::ZERO);
        let (_q, r) = a.div_rem(&m.to_nz().unwrap());
        assert_eq!(r.lo(), U256::from_limbs([9, 0, 0, 0]));

        // Larger: verify cross-limb reduction
        let big = U512::from_lo_hi(
            U256::from_limbs([u64::MAX, u64::MAX, 0, 0]),
            U256::ZERO,
        );
        let p = U512::from_lo_hi(U256::from_limbs([0xFFFFFFFFFFFFFFFD, u64::MAX, 0, 0]), U256::ZERO);
        let (_q2, r2) = big.div_rem(&p.to_nz().unwrap());
        assert_eq!(r2, U512::from_lo_hi(U256::from_limbs([2, 0, 0, 0]), U256::ZERO));
    }
}
