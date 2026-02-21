// Keccak-256 (Ethereum variant) with field element hashing glue.
//
// The core Keccak-256 primitive is provided by `tiny-keccak`. This module adds
// the field-element serialization layer that Barretenberg uses to hash BN254
// field elements represented as 4 x u64 limbs.
//
// C++ source: barretenberg/cpp/src/barretenberg/crypto/keccak/keccak.cpp
// C++ header: barretenberg/cpp/src/barretenberg/crypto/keccak/keccak.hpp
//
// The core `ethash_keccak256` function is standard Ethereum Keccak-256 (0x01
// padding, not SHA-3's 0x06), so we delegate directly to `tiny-keccak`.

use tiny_keccak::{Hasher, Keccak};

/// 256-bit Keccak hash output, stored as 4 x u64 words.
///
/// Mirrors C++ `struct keccak256 { uint64_t word64s[4]; }` from hash_types.hpp.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Keccak256([u64; 4]);

impl Keccak256 {
    pub fn word64s(&self) -> &[u64; 4] {
        &self.0
    }

    pub fn to_bytes(&self) -> [u8; 32] {
        let mut out = [0u8; 32];
        for (i, word) in self.0.iter().enumerate() {
            out[i * 8..(i + 1) * 8].copy_from_slice(&word.to_le_bytes());
        }
        out
    }
}

/// Standard Ethereum Keccak-256 hash of arbitrary bytes.
///
/// Equivalent to C++ `ethash_keccak256(data, size)`.
pub fn keccak256(data: &[u8]) -> Keccak256 {
    let mut hasher = Keccak::v256();
    hasher.update(data);
    let mut output = [0u8; 32];
    hasher.finalize(&mut output);

    // Pack bytes into u64 words (little-endian), matching C++ output layout
    let mut word64s = [0u64; 4];
    for i in 0..4 {
        word64s[i] = u64::from_le_bytes(output[i * 8..(i + 1) * 8].try_into().unwrap());
    }
    Keccak256(word64s)
}

/// Hash field elements represented as 4 x u64 limbs each.
///
/// Each element's limbs are serialized as big-endian bytes before hashing.
/// This matches the C++ `hash_field_elements(limbs, num_elements)` function
/// in barretenberg/cpp/src/barretenberg/crypto/keccak/keccak.cpp lines 114-134.
///
/// `limbs` is a flat slice where each consecutive 4 entries represent one
/// 256-bit field element (e.g. BN254 Fr).
pub fn hash_field_elements(limbs: &[u64]) -> Keccak256 {
    assert!(limbs.len() % 4 == 0, "limbs length must be a multiple of 4");
    let num_elements = limbs.len() / 4;
    let mut buffer = vec![0u8; num_elements * 32];

    for i in 0..num_elements {
        for j in 0..4 {
            let word = limbs[i * 4 + j];
            let idx = i * 32 + j * 8;
            buffer[idx..idx + 8].copy_from_slice(&word.to_be_bytes());
        }
    }

    keccak256(&buffer)
}

/// Hash a single field element (4 x u64 limbs).
///
/// Equivalent to C++ `hash_field_element(limb)`.
pub fn hash_field_element(limbs: &[u64; 4]) -> Keccak256 {
    hash_field_elements(limbs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_hash_matches_known_keccak256() {
        // Keccak-256 of empty input is a well-known constant
        let result = keccak256(&[]);
        let bytes = result.to_bytes();
        let expected: [u8; 32] = [
            0xc5, 0xd2, 0x46, 0x01, 0x86, 0xf7, 0x23, 0x3c, 0x92, 0x7e, 0x7d, 0xb2, 0xdc, 0xc7,
            0x03, 0xc0, 0xe5, 0x00, 0xb6, 0x53, 0xca, 0x82, 0x27, 0x3b, 0x7b, 0xfa, 0xd8, 0x04,
            0x5d, 0x85, 0xa4, 0x70,
        ];
        assert_eq!(bytes, expected);
    }

    #[test]
    fn field_element_serialization_is_big_endian() {
        // Verify that limbs are serialized in big-endian byte order, matching
        // the C++ implementation's manual byte extraction.
        let limbs: [u64; 4] = [
            0x0102030405060708,
            0x090a0b0c0d0e0f10,
            0x1112131415161718,
            0x191a1b1c1d1e1f20,
        ];
        let result = hash_field_element(&limbs);

        // Construct the expected input buffer manually
        let mut expected_input = [0u8; 32];
        expected_input[0..8].copy_from_slice(&limbs[0].to_be_bytes());
        expected_input[8..16].copy_from_slice(&limbs[1].to_be_bytes());
        expected_input[16..24].copy_from_slice(&limbs[2].to_be_bytes());
        expected_input[24..32].copy_from_slice(&limbs[3].to_be_bytes());

        let expected = keccak256(&expected_input);
        assert_eq!(result, expected);
    }

    #[test]
    fn hash_multiple_field_elements() {
        let limbs: [u64; 8] = [1, 2, 3, 4, 5, 6, 7, 8];
        let result = hash_field_elements(&limbs);

        // Build expected buffer: two field elements, each 4 limbs big-endian
        let mut expected_input = [0u8; 64];
        for i in 0..8 {
            let offset = i * 8;
            expected_input[offset..offset + 8].copy_from_slice(&limbs[i].to_be_bytes());
        }
        let expected = keccak256(&expected_input);
        assert_eq!(result, expected);
    }

    #[test]
    #[should_panic(expected = "limbs length must be a multiple of 4")]
    fn rejects_invalid_limb_count() {
        hash_field_elements(&[1, 2, 3]);
    }
}
