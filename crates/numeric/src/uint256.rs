// 256-bit unsigned integer type.
//
// Wraps `crypto_bigint::U256` and adds BB-compatible convenience methods.
//
// C++ source: barretenberg/cpp/src/barretenberg/numeric/uint256/uint256.hpp
//
// BB's uint256_t is stored as 4 x u64 limbs (data[0..3]) in little-endian
// limb order (data[0] is least significant). crypto-bigint uses the same
// internal layout.

use crypto_bigint::Uint;

/// 256-bit unsigned integer, backed by `crypto_bigint::U256`.
pub type U256 = Uint<4>;

/// BB-compatible extension methods for U256.
///
/// These mirror the methods on BB's `uint256_t` that don't have direct
/// equivalents in crypto-bigint's API.
pub trait U256Ext {
    /// Get the position of the most significant bit (0-indexed).
    /// Returns 0 for zero input (matching BB behavior).
    ///
    /// Mirrors BB's `uint256_t::get_msb()`.
    fn get_msb(&self) -> u32;

    /// Extract a single bit at the given index.
    ///
    /// Mirrors BB's `uint256_t::get_bit(index)`.
    fn get_bit(&self, index: u32) -> bool;

    /// Extract a bit-range [start, end) as a u64.
    ///
    /// Mirrors BB's `uint256_t::slice(start, end)`. The range must fit in 64
    /// bits (end - start <= 64).
    fn slice(&self, start: u32, end: u32) -> u64;

    /// Construct from 4 x u64 limbs in little-endian limb order.
    ///
    /// Matches BB's `uint256_t(lo, mid_lo, mid_hi, hi)` constructor where
    /// `data[0] = lo` is least significant.
    fn from_limbs(limbs: [u64; 4]) -> Self;

    /// Access the raw u64 limbs in little-endian limb order.
    fn limbs(&self) -> [u64; 4];
}

impl U256Ext for U256 {
    fn get_msb(&self) -> u32 {
        let bits = self.bits_vartime();
        if bits == 0 { 0 } else { bits - 1 }
    }

    fn get_bit(&self, index: u32) -> bool {
        self.bit_vartime(index)
    }

    fn slice(&self, start: u32, end: u32) -> u64 {
        assert!(end > start, "end must be greater than start");
        assert!(end - start <= 64, "slice range must fit in u64");

        let shifted = self.wrapping_shr_vartime(start);
        let mask = if end - start == 64 {
            u64::MAX
        } else {
            (1u64 << (end - start)) - 1
        };
        shifted.as_words()[0] & mask
    }

    fn from_limbs(limbs: [u64; 4]) -> Self {
        U256::from_words(limbs)
    }

    fn limbs(&self) -> [u64; 4] {
        *self.as_words()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::{Encoding, U256, Zero};

    #[test]
    fn from_limbs_roundtrip() {
        let limbs = [0x1111_2222_3333_4444u64, 0x5555_6666_7777_8888, 0x9999_aaaa_bbbb_cccc, 0xdddd_eeee_ffff_0000];
        let val = U256::from_limbs(limbs);
        assert_eq!(val.limbs(), limbs);
    }

    #[test]
    fn get_msb_basic() {
        assert_eq!(U256::ZERO.get_msb(), 0);
        assert_eq!(U256::ONE.get_msb(), 0);
        assert_eq!(U256::from_limbs([0, 0, 0, 1]).get_msb(), 192);
        assert_eq!(U256::from_limbs([0, 0, 0, 1 << 63]).get_msb(), 255);
    }

    #[test]
    fn get_bit_basic() {
        let val = U256::from_limbs([0b1010, 0, 0, 0]);
        assert!(val.get_bit(1));
        assert!(!val.get_bit(2));
        assert!(val.get_bit(3));
        assert!(!val.get_bit(4));
    }

    #[test]
    fn slice_basic() {
        let val = U256::from_limbs([0xABCD_EF01_2345_6789, 0, 0, 0]);
        // Extract bottom 16 bits
        assert_eq!(val.slice(0, 16), 0x6789);
        // Extract bits [16, 32)
        assert_eq!(val.slice(16, 32), 0x2345);
    }

    #[test]
    fn slice_cross_limb() {
        let val = U256::from_limbs([u64::MAX, 0x00FF, 0, 0]);
        // Extract bits [60, 72) — crosses limb boundary
        assert_eq!(val.slice(60, 72), 0xFF_F);
    }

    #[test]
    fn big_endian_roundtrip() {
        let limbs = [1u64, 2, 3, 4];
        let val = U256::from_limbs(limbs);
        let bytes = val.to_be_bytes();
        let recovered = U256::from_be_bytes(bytes);
        assert_eq!(val, recovered);
    }

    #[test]
    fn div_rem_basic() {
        let a = U256::from_limbs([100, 0, 0, 0]);
        let b = U256::from_limbs([7, 0, 0, 0]);
        let (q, r) = a.div_rem(&b.to_nz().unwrap());
        assert_eq!(q, U256::from_limbs([14, 0, 0, 0]));
        assert_eq!(r, U256::from_limbs([2, 0, 0, 0]));
    }

    #[test]
    fn widening_mul_basic() {
        let a = U256::from_limbs([u64::MAX, u64::MAX, u64::MAX, u64::MAX]);
        let b = U256::from_limbs([2, 0, 0, 0]);
        let wide = a.widening_mul(&b);
        // (2^256 - 1) * 2 = 2^257 - 2
        // lo should be 0xFFFF...FFFE, hi should be 1
        let (lo, hi) = wide.split();
        assert_eq!(hi, U256::ONE);
        assert_eq!(lo, U256::from_limbs([u64::MAX - 1, u64::MAX, u64::MAX, u64::MAX]));
    }

    #[test]
    fn u256_add() {
        let a = U256::from_limbs([100, 0, 0, 0]);
        let b = U256::from_limbs([200, 0, 0, 0]);
        let c = a.wrapping_add(&b);
        assert_eq!(c, U256::from_limbs([300, 0, 0, 0]));
    }

    #[test]
    fn u256_add_overflow() {
        let a = U256::from_limbs([u64::MAX, 0, 0, 0]);
        let b = U256::from_limbs([1, 0, 0, 0]);
        let c = a.wrapping_add(&b);
        assert_eq!(c, U256::from_limbs([0, 1, 0, 0]));
    }

    #[test]
    fn u256_sub() {
        let a = U256::from_limbs([300, 0, 0, 0]);
        let b = U256::from_limbs([100, 0, 0, 0]);
        let c = a.wrapping_sub(&b);
        assert_eq!(c, U256::from_limbs([200, 0, 0, 0]));
    }

    #[test]
    fn u256_mul() {
        let a = U256::from_limbs([7, 0, 0, 0]);
        let b = U256::from_limbs([6, 0, 0, 0]);
        let c = a.wrapping_mul(&b);
        assert_eq!(c, U256::from_limbs([42, 0, 0, 0]));
    }

    #[test]
    fn u256_right_shift() {
        let a = U256::from_limbs([0x100, 0, 0, 0]);
        let b = a.wrapping_shr_vartime(4);
        assert_eq!(b, U256::from_limbs([0x10, 0, 0, 0]));
    }

    #[test]
    fn u256_left_shift() {
        let a = U256::from_limbs([0x10, 0, 0, 0]);
        let b = a.wrapping_shl_vartime(4);
        assert_eq!(b, U256::from_limbs([0x100, 0, 0, 0]));
    }

    #[test]
    fn u256_and() {
        let a = U256::from_limbs([0xFF00, 0, 0, 0]);
        let b = U256::from_limbs([0x0FF0, 0, 0, 0]);
        let c = a.bitand(&b);
        assert_eq!(c, U256::from_limbs([0x0F00, 0, 0, 0]));
    }

    #[test]
    fn u256_or() {
        let a = U256::from_limbs([0xFF00, 0, 0, 0]);
        let b = U256::from_limbs([0x00FF, 0, 0, 0]);
        let c = a.bitor(&b);
        assert_eq!(c, U256::from_limbs([0xFFFF, 0, 0, 0]));
    }

    #[test]
    fn u256_xor() {
        let a = U256::from_limbs([0xFF00, 0, 0, 0]);
        let b = U256::from_limbs([0x0FF0, 0, 0, 0]);
        let c = a.bitxor(&b);
        assert_eq!(c, U256::from_limbs([0xF0F0, 0, 0, 0]));
    }

    #[test]
    fn u256_equality() {
        let a = U256::from_limbs([42, 0, 0, 0]);
        let b = U256::from_limbs([42, 0, 0, 0]);
        let c = U256::from_limbs([43, 0, 0, 0]);
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn u256_ordering() {
        let a = U256::from_limbs([42, 0, 0, 0]);
        let b = U256::from_limbs([43, 0, 0, 0]);
        assert!(a < b);
        assert!(b > a);
        assert!(a <= a);
        assert!(a >= a);
    }

    #[test]
    fn u256_bit_not() {
        let a = U256::from_limbs([0xFF00, 0, 0, 0]);
        let not_a = a.not();
        assert_eq!(not_a.limbs()[0], !0xFF00u64);
        assert_eq!(not_a.limbs()[1], u64::MAX);
        assert_eq!(not_a.limbs()[2], u64::MAX);
        assert_eq!(not_a.limbs()[3], u64::MAX);
        // Double NOT should return original
        assert_eq!(not_a.not(), a);
    }

    #[test]
    fn u256_logic_not() {
        assert!(bool::from(U256::ZERO.is_zero()));
        assert!(!bool::from(U256::ONE.is_zero()));
        assert!(!bool::from(U256::from_limbs([0, 0, 0, 1]).is_zero()));
    }

    #[test]
    fn u256_from_string() {
        // Construct from big-endian hex string and verify matches from_limbs
        let hex_str = "000000000000000000000000000000000000000000000000000000000000002A";
        let from_hex = U256::from_be_hex(hex_str);
        let from_limbs = U256::from_limbs([42, 0, 0, 0]);
        assert_eq!(from_hex, from_limbs);

        // Larger value: 0x1234...
        let hex_str2 = "00000000000000010000000000000002000000000000000300000000000000FF";
        let from_hex2 = U256::from_be_hex(hex_str2);
        let from_limbs2 = U256::from_limbs([0x00000000000000FF, 0x0000000000000003, 0x0000000000000002, 0x0000000000000001]);
        assert_eq!(from_hex2, from_limbs2);
    }
}
