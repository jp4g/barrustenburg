// Bit manipulation utilities.
//
// C++ source: barretenberg/cpp/src/barretenberg/numeric/bitop/
//
// Most of BB's bitop functions map directly to Rust built-in methods on
// primitive types. We provide standalone functions here to match BB's
// `numeric::get_msb()` calling convention used by downstream modules
// (e.g., wnaf.hpp calls `numeric::get_msb(uint64_t)`).

/// Position of the most significant set bit (0-indexed).
/// Returns 0 for input 0 (matching BB's `get_msb` behavior).
///
/// Mirrors BB's `numeric::get_msb<T>()` from `get_msb.hpp`.
#[inline]
pub fn get_msb64(val: u64) -> u32 {
    if val == 0 { 0 } else { 63 - val.leading_zeros() }
}

/// Position of the most significant set bit for u32.
#[inline]
pub fn get_msb32(val: u32) -> u32 {
    if val == 0 { 0 } else { 31 - val.leading_zeros() }
}

/// Count leading zeros.
///
/// Mirrors BB's `numeric::count_leading_zeros<T>()` from `count_leading_zeros.hpp`.
#[inline]
pub fn count_leading_zeros64(val: u64) -> u32 {
    val.leading_zeros()
}

/// Keep only the `n` least significant bits, zeroing all higher bits.
///
/// Mirrors BB's `numeric::keep_n_lsb()` from `keep_n_lsb.hpp`.
#[inline]
pub fn keep_n_lsb(val: u64, n: u32) -> u64 {
    if n >= 64 { val } else { val & ((1u64 << n) - 1) }
}

/// Integer exponentiation via binary method.
///
/// Mirrors BB's `numeric::pow64()` from `pow.hpp`.
#[inline]
pub fn pow64(base: u64, exp: u64) -> u64 {
    let mut result = 1u64;
    let mut b = base;
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            result = result.wrapping_mul(b);
        }
        b = b.wrapping_mul(b);
        e >>= 1;
    }
    result
}

/// Check if a value is a power of two.
#[inline]
pub fn is_power_of_two(val: u64) -> bool {
    val != 0 && (val & (val - 1)) == 0
}

/// Round up to the next power of two.
///
/// Mirrors BB's `numeric::round_up_power_2()`.
#[inline]
pub fn round_up_power_2(val: u64) -> u64 {
    val.next_power_of_two()
}

/// Ceiling division: ceil(numerator / denominator).
///
/// Mirrors BB's `numeric::ceil_div()` from `general/general.hpp`.
#[inline]
pub fn ceil_div(numerator: u64, denominator: u64) -> u64 {
    (numerator + denominator - 1) / denominator
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn get_msb64_cases() {
        assert_eq!(get_msb64(0), 0);
        assert_eq!(get_msb64(1), 0);
        assert_eq!(get_msb64(2), 1);
        assert_eq!(get_msb64(0xFF), 7);
        assert_eq!(get_msb64(1 << 63), 63);
    }

    #[test]
    fn keep_n_lsb_cases() {
        assert_eq!(keep_n_lsb(0xFF, 4), 0x0F);
        assert_eq!(keep_n_lsb(u64::MAX, 1), 1);
        assert_eq!(keep_n_lsb(u64::MAX, 64), u64::MAX);
    }

    #[test]
    fn pow64_cases() {
        assert_eq!(pow64(2, 10), 1024);
        assert_eq!(pow64(3, 0), 1);
        assert_eq!(pow64(0, 5), 0);
    }

    #[test]
    fn ceil_div_cases() {
        assert_eq!(ceil_div(10, 3), 4);
        assert_eq!(ceil_div(9, 3), 3);
        assert_eq!(ceil_div(1, 1), 1);
    }

    #[test]
    fn power_of_two_checks() {
        assert!(is_power_of_two(1));
        assert!(is_power_of_two(256));
        assert!(!is_power_of_two(0));
        assert!(!is_power_of_two(3));
    }
}
