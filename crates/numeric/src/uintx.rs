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
    use crate::U256;

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
}
