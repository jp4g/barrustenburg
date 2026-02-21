// Windowed Non-Adjacent Form (WNAF) encoding for scalar multiplication.
//
// Ported from C++ `ecc/groups/wnaf.hpp`.
// Used by endomorphism-based scalar multiplication and Pippenger MSM.

/// Number of bits in the half-scalar (after endomorphism splitting).
pub const SCALAR_BITS: usize = 127;

/// Extract `bits` bits from a 128-bit scalar starting at `bit_position`.
///
/// The scalar is stored as [u64; 2] in little-endian limb order.
#[inline]
pub fn get_wnaf_bits(scalar: &[u64; 2], bits: u64, bit_position: u64) -> u64 {
    let lo_limb_idx = (bit_position >> 6) as usize;
    let hi_limb_idx = ((bit_position + bits - 1) >> 6) as usize;
    let lo_shift = bit_position & 63;
    let bit_mask = (1u64 << bits) - 1;

    let lo = scalar[lo_limb_idx] >> lo_shift;
    let hi_shift = if lo_shift != 0 { 64 - lo_shift } else { 0 };
    // Use wrapping_shl to avoid panic when hi_shift == 64 (result is masked out anyway)
    let hi = scalar[hi_limb_idx].wrapping_shl(hi_shift as u32);
    let hi_mask = bit_mask & (0u64.wrapping_sub((lo_limb_idx != hi_limb_idx) as u64));

    (lo & bit_mask) | (hi & hi_mask)
}

/// Compute fixed-window WNAF representation for a 128-bit scalar.
///
/// - `scalar`: 128-bit half-scalar as [u64; 2]
/// - `wnaf`: output slice, must have room for at least `wnaf_entries * num_points` entries
/// - `skew`: set to true if scalar is even (requires skew correction)
/// - `point_index`: index of this point (encoded in upper bits of WNAF entry)
/// - `num_points`: stride between WNAF entries (for interleaved storage)
/// - `wnaf_bits`: window width in bits
///
/// WNAF entry format: bits 0-30 = abs value, bit 31 = sign flag, bits 32-63 = point_index.
pub fn fixed_wnaf(
    scalar: &[u64; 2],
    wnaf: &mut [u64],
    skew: &mut bool,
    point_index: u64,
    num_points: u64,
    wnaf_bits: usize,
) {
    *skew = (scalar[0] & 1) == 0;
    let mut previous = get_wnaf_bits(scalar, wnaf_bits as u64, 0) + (*skew as u64);
    let wnaf_entries = (SCALAR_BITS + wnaf_bits - 1) / wnaf_bits;

    for round_i in 1..(wnaf_entries - 1) {
        let slice = get_wnaf_bits(scalar, wnaf_bits as u64, (round_i * wnaf_bits) as u64);
        let predicate = ((slice & 1) == 0) as u64;
        let idx = (wnaf_entries - round_i) * (num_points as usize);
        wnaf[idx] = ((((previous.wrapping_sub(predicate << wnaf_bits))
            ^ (0u64.wrapping_sub(predicate)))
            >> 1)
            | (predicate << 31))
            | point_index;
        previous = slice + predicate;
    }

    let final_bits = SCALAR_BITS - (wnaf_bits * (wnaf_entries - 1));
    let slice = get_wnaf_bits(
        scalar,
        final_bits as u64,
        ((wnaf_entries - 1) * wnaf_bits) as u64,
    );
    let predicate = ((slice & 1) == 0) as u64;

    wnaf[num_points as usize] = ((((previous.wrapping_sub(predicate << wnaf_bits))
        ^ (0u64.wrapping_sub(predicate)))
        >> 1)
        | (predicate << 31))
        | point_index;
    wnaf[0] = ((slice + predicate) >> 1) | point_index;
}

/// Determine optimal bucket width for Pippenger based on point count.
pub fn get_optimal_bucket_width(num_points: usize) -> usize {
    if num_points >= 14617149 {
        return 21;
    }
    if num_points >= 1139094 {
        return 18;
    }
    if num_points >= 155975 {
        return 15;
    }
    if num_points >= 144834 {
        return 14;
    }
    if num_points >= 25067 {
        return 12;
    }
    if num_points >= 13926 {
        return 11;
    }
    if num_points >= 7659 {
        return 10;
    }
    if num_points >= 2436 {
        return 9;
    }
    if num_points >= 376 {
        return 7;
    }
    if num_points >= 231 {
        return 6;
    }
    if num_points >= 97 {
        return 5;
    }
    if num_points >= 35 {
        return 4;
    }
    if num_points >= 10 {
        return 3;
    }
    if num_points >= 2 {
        return 2;
    }
    1
}

/// Get the number of buckets for Pippenger MSM.
pub fn get_num_buckets(num_points: usize) -> usize {
    let bits_per_bucket = get_optimal_bucket_width(num_points / 2);
    1 << bits_per_bucket
}

/// Get the number of rounds for Pippenger MSM.
pub fn get_num_rounds(num_points: usize) -> usize {
    let bits_per_bucket = get_optimal_bucket_width(num_points / 2);
    (SCALAR_BITS + bits_per_bucket) / (bits_per_bucket + 1)
}
