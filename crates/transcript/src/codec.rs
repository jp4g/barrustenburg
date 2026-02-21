//! Field codec for transcript serialization.
//!
//! Port of C++ `FrCodec` from `ecc/fields/field_conversion.hpp`.
//! Handles serialization/deserialization of various types to/from BN254 Fr elements.

use bbrs_ecc::curves::bn254::{Bn254FqParams, Fr};
use bbrs_ecc::curves::grumpkin::GrumpkinFr;
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;

// Number of limb bits in stdlib field simulation (matches C++ NUM_LIMB_BITS_IN_FIELD_SIMULATION)
const NUM_LIMB_BITS: u32 = 68;

/// Trait for types that can be serialized to/from Fr field elements for the transcript.
pub trait FieldSerializable: Sized {
    /// Number of Fr elements needed to represent this type.
    const NUM_FR: usize;

    /// Serialize to a vector of Fr elements.
    fn serialize_to_frs(&self) -> Vec<Fr>;

    /// Deserialize from a slice of Fr elements.
    fn deserialize_from_frs(frs: &[Fr]) -> Self;
}

// --- Primitive types ---

impl FieldSerializable for bool {
    const NUM_FR: usize = 1;

    fn serialize_to_frs(&self) -> Vec<Fr> {
        vec![Fr::from(u64::from(*self))]
    }

    fn deserialize_from_frs(frs: &[Fr]) -> Self {
        assert_eq!(frs.len(), 1);
        let reduced = frs[0].from_montgomery_form();
        reduced.data[0] != 0
    }
}

impl FieldSerializable for u32 {
    const NUM_FR: usize = 1;

    fn serialize_to_frs(&self) -> Vec<Fr> {
        vec![Fr::from(u64::from(*self))]
    }

    fn deserialize_from_frs(frs: &[Fr]) -> Self {
        assert_eq!(frs.len(), 1);
        let reduced = frs[0].from_montgomery_form();
        reduced.data[0] as u32
    }
}

impl FieldSerializable for u64 {
    const NUM_FR: usize = 1;

    fn serialize_to_frs(&self) -> Vec<Fr> {
        vec![Fr::from(*self)]
    }

    fn deserialize_from_frs(frs: &[Fr]) -> Self {
        assert_eq!(frs.len(), 1);
        let reduced = frs[0].from_montgomery_form();
        reduced.data[0]
    }
}

// --- BN254 Fr (1 field element) ---

impl FieldSerializable for Fr {
    const NUM_FR: usize = 1;

    fn serialize_to_frs(&self) -> Vec<Fr> {
        vec![*self]
    }

    fn deserialize_from_frs(frs: &[Fr]) -> Self {
        assert_eq!(frs.len(), 1);
        frs[0]
    }
}

// --- Fixed-size arrays of Fr ---

impl<const N: usize> FieldSerializable for [Fr; N] {
    const NUM_FR: usize = N;

    fn serialize_to_frs(&self) -> Vec<Fr> {
        self.to_vec()
    }

    fn deserialize_from_frs(frs: &[Fr]) -> Self {
        assert_eq!(frs.len(), N);
        let mut arr = [Fr::zero(); N];
        arr.copy_from_slice(frs);
        arr
    }
}

// --- Grumpkin Fr = BN254 Fq (2 field elements) ---
//
// Splits a Grumpkin Fr (BN254 Fq) into two BN254 Fr elements:
//   - lo: lower 136 bits (2 * NUM_LIMB_BITS)
//   - hi: upper 118 bits (254 - 136)
// Mirrors stdlib bigfield limbs (68-bit each → 136-bit lower chunk).

impl FieldSerializable for GrumpkinFr {
    const NUM_FR: usize = 2;

    fn serialize_to_frs(&self) -> Vec<Fr> {
        convert_grumpkin_fr_to_bn254_frs(self)
    }

    fn deserialize_from_frs(frs: &[Fr]) -> Self {
        assert_eq!(frs.len(), 2);
        convert_grumpkin_fr_from_bn254_frs(frs)
    }
}

// --- Affine curve points ---

/// Implementation for AffineElement where the base field IS Fr (e.g., Grumpkin G1).
/// Each coordinate is 1 Fr → total 2 Fr.
impl<C> FieldSerializable for AffineElement<C>
where
    C: CurveParams,
    Field<C::BaseFieldParams>: FieldSerializable,
{
    const NUM_FR: usize = 2 * <Field<C::BaseFieldParams> as FieldSerializable>::NUM_FR;

    fn serialize_to_frs(&self) -> Vec<Fr> {
        if self.is_point_at_infinity() {
            let zero = Field::<C::BaseFieldParams>::zero();
            let mut frs = zero.serialize_to_frs();
            frs.extend(zero.serialize_to_frs());
            return frs;
        }
        let mut frs = self.x.serialize_to_frs();
        frs.extend(self.y.serialize_to_frs());
        frs
    }

    fn deserialize_from_frs(frs: &[Fr]) -> Self {
        let base_size = <Field<C::BaseFieldParams> as FieldSerializable>::NUM_FR;
        assert_eq!(frs.len(), 2 * base_size);

        let x = Field::<C::BaseFieldParams>::deserialize_from_frs(&frs[..base_size]);
        let y = Field::<C::BaseFieldParams>::deserialize_from_frs(&frs[base_size..]);

        if x == Field::<C::BaseFieldParams>::zero() && y == Field::<C::BaseFieldParams>::zero() {
            let mut pt = AffineElement::new(x, y);
            pt.self_set_infinity();
            return pt;
        }

        let pt = AffineElement::new(x, y);
        assert!(pt.on_curve(), "Deserialized point not on curve");
        pt
    }
}

// --- Grumpkin Fr <-> BN254 Fr conversion helpers ---

/// Convert 2 BN254 Fr elements to a Grumpkin Fr (BN254 Fq).
///
/// Inverse of `convert_grumpkin_fr_to_bn254_frs`.
/// lo covers lower 136 bits, hi covers upper bits.
fn convert_grumpkin_fr_from_bn254_frs(frs: &[Fr]) -> GrumpkinFr {
    assert_eq!(frs.len(), 2);

    let lower_bits = 2 * NUM_LIMB_BITS; // 136
    let total_bits: u32 = 254;

    // Extract lo and hi in standard form
    let lo = frs[0].from_montgomery_form();
    let hi = frs[1].from_montgomery_form();

    // Validate bounds
    // lo must fit in lower_bits (136 bits)
    assert!(
        limbs_msb(&lo.data) < lower_bits,
        "Conversion error: lo exceeds 136 bits"
    );
    // hi must fit in (total_bits - lower_bits) = 118 bits
    assert!(
        limbs_msb(&hi.data) < (total_bits - lower_bits),
        "Conversion error: hi exceeds 118 bits"
    );

    // Reconstruct: value = lo + (hi << 136)
    let combined = limbs_add_shifted(&lo.data, &hi.data, lower_bits);

    // Convert to GrumpkinFr (BN254 Fq)
    Field::<Bn254FqParams>::from_limbs(combined)
}

/// Convert a Grumpkin Fr (BN254 Fq) to 2 BN254 Fr elements.
///
/// Splits into 136-bit lower chunk and remaining upper chunk.
fn convert_grumpkin_fr_to_bn254_frs(val: &GrumpkinFr) -> Vec<Fr> {
    let lower_bits = 2 * NUM_LIMB_BITS; // 136

    // Get value in standard form
    let standard = val.from_montgomery_form();

    // Extract lo (lower 136 bits) and hi (upper bits)
    let lo = limbs_mask(&standard.data, lower_bits);
    let hi = limbs_shr(&standard.data, lower_bits);

    vec![Fr::from_limbs(lo), Fr::from_limbs(hi)]
}

/// FrCodec: Split a challenge Fr element into two 127-bit halves.
///
/// Both `lo` and `hi` are 127 bits each (254/2), providing >100-bit security.
pub fn split_challenge(challenge: &Fr) -> [Fr; 2] {
    const TOTAL_BITS: u32 = 254; // BN254 Fr modulus MSB + 1
    const LO_BITS: u32 = TOTAL_BITS / 2; // 127

    let standard = challenge.from_montgomery_form();
    let lo = limbs_mask(&standard.data, LO_BITS);
    let hi = limbs_shr(&standard.data, LO_BITS);

    [Fr::from_limbs(lo), Fr::from_limbs(hi)]
}

// --- Low-level limb arithmetic helpers ---

/// Get the MSB position of a 256-bit number stored in 4 limbs (0-indexed, returns 0 for zero).
fn limbs_msb(data: &[u64; 4]) -> u32 {
    for i in (0..4).rev() {
        if data[i] != 0 {
            return (i as u32) * 64 + (63 - data[i].leading_zeros());
        }
    }
    0
}

/// Mask a 256-bit limb array to keep only the lower `bits` bits.
fn limbs_mask(data: &[u64; 4], bits: u32) -> [u64; 4] {
    let mut result = [0u64; 4];
    let full_limbs = (bits / 64) as usize;
    let remaining = bits % 64;

    for i in 0..full_limbs.min(4) {
        result[i] = data[i];
    }
    if full_limbs < 4 && remaining > 0 {
        result[full_limbs] = data[full_limbs] & ((1u64 << remaining) - 1);
    }
    result
}

/// Shift a 256-bit limb array right by `shift` bits.
fn limbs_shr(data: &[u64; 4], shift: u32) -> [u64; 4] {
    let mut result = [0u64; 4];
    let limb_shift = (shift / 64) as usize;
    let bit_shift = shift % 64;

    if bit_shift == 0 {
        for i in 0..(4 - limb_shift) {
            result[i] = data[i + limb_shift];
        }
    } else {
        for i in 0..(4 - limb_shift) {
            result[i] = data[i + limb_shift] >> bit_shift;
            if i + limb_shift + 1 < 4 {
                result[i] |= data[i + limb_shift + 1] << (64 - bit_shift);
            }
        }
    }
    result
}

/// Add hi << shift to lo (256-bit limb arithmetic, no overflow check).
fn limbs_add_shifted(lo: &[u64; 4], hi: &[u64; 4], shift: u32) -> [u64; 4] {
    let mut result = *lo;
    let limb_offset = (shift / 64) as usize;
    let bit_offset = shift % 64;

    if bit_offset == 0 {
        for i in 0..4 {
            if i + limb_offset < 4 {
                let (sum, carry) = result[i + limb_offset].overflowing_add(hi[i]);
                result[i + limb_offset] = sum;
                // Propagate carry
                if carry {
                    let mut j = i + limb_offset + 1;
                    while j < 4 {
                        let (s, c) = result[j].overflowing_add(1);
                        result[j] = s;
                        if !c {
                            break;
                        }
                        j += 1;
                    }
                }
            }
        }
    } else {
        for i in 0..4 {
            let shifted_lo = hi[i] << bit_offset;
            let carry_in = if i > 0 {
                hi[i - 1] >> (64 - bit_offset)
            } else {
                0
            };
            let val = shifted_lo | carry_in;

            if i + limb_offset < 4 {
                let (sum, carry) = result[i + limb_offset].overflowing_add(val);
                result[i + limb_offset] = sum;
                if carry {
                    let mut j = i + limb_offset + 1;
                    while j < 4 {
                        let (s, c) = result[j].overflowing_add(1);
                        result[j] = s;
                        if !c {
                            break;
                        }
                        j += 1;
                    }
                }
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fr_roundtrip() {
        let val = Fr::from(42u64);
        let frs = val.serialize_to_frs();
        assert_eq!(frs.len(), 1);
        let recovered = Fr::deserialize_from_frs(&frs);
        assert_eq!(val, recovered);
    }

    #[test]
    fn test_u32_roundtrip() {
        let val: u32 = 12345;
        let frs = val.serialize_to_frs();
        assert_eq!(frs.len(), 1);
        let recovered = u32::deserialize_from_frs(&frs);
        assert_eq!(val, recovered);
    }

    #[test]
    fn test_bool_roundtrip() {
        for val in [true, false] {
            let frs = val.serialize_to_frs();
            assert_eq!(frs.len(), 1);
            let recovered = bool::deserialize_from_frs(&frs);
            assert_eq!(val, recovered);
        }
    }

    #[test]
    fn test_grumpkin_fr_roundtrip() {
        // Create a GrumpkinFr value and verify roundtrip through 2 Fr elements
        let val = GrumpkinFr::from(123456789u64);
        let frs = val.serialize_to_frs();
        assert_eq!(frs.len(), 2);
        let recovered = GrumpkinFr::deserialize_from_frs(&frs);
        assert_eq!(val, recovered);
    }

    #[test]
    fn test_split_challenge() {
        let challenge = Fr::from(0xDEADBEEFu64);
        let [lo, hi] = split_challenge(&challenge);
        // For small values, hi should be zero and lo should equal the original
        assert_eq!(hi, Fr::zero());
        assert_eq!(lo, challenge);
    }

    #[test]
    fn test_split_challenge_large() {
        // Create a value with bits in both halves
        let val = Fr::from_limbs([0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x7F, 0]);
        let [lo, hi] = split_challenge(&val);
        // lo should have 127 bits, hi should have the rest
        let lo_std = lo.from_montgomery_form();
        let hi_std = hi.from_montgomery_form();
        assert!(limbs_msb(&lo_std.data) < 127);
        assert!(limbs_msb(&hi_std.data) < 127);
    }

    #[test]
    fn test_limbs_mask() {
        let data = [0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF];
        let masked = limbs_mask(&data, 136);
        assert_eq!(masked[0], 0xFFFFFFFFFFFFFFFF);
        assert_eq!(masked[1], 0xFFFFFFFFFFFFFFFF);
        assert_eq!(masked[2], 0xFF); // only lower 8 bits of limb 2
        assert_eq!(masked[3], 0);
    }

    #[test]
    fn test_limbs_shr() {
        let data = [0, 0, 0xFF00, 0];
        let shifted = limbs_shr(&data, 136);
        // 136 = 2*64 + 8, so limb[2] >> 8 = 0xFF
        assert_eq!(shifted[0], 0xFF);
        assert_eq!(shifted[1], 0);
    }
}
