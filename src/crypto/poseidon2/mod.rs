pub mod params;
pub mod permutation;
pub mod sponge;

use crate::ecc::curves::bn254::Fr;

/// Poseidon2 hash function over BN254 Fr.
pub struct Poseidon2;

impl Poseidon2 {
    /// Hash a slice of field elements into a single field element.
    pub fn hash(input: &[Fr]) -> Fr {
        sponge::hash(input)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecc::curves::bn254::Fr;

    #[test]
    fn test_permutation_vector() {
        let input = *params::TEST_VECTOR_INPUT;
        let expected = *params::TEST_VECTOR_OUTPUT;
        let result = permutation::permutation(&input);

        for i in 0..4 {
            assert_eq!(result[i], expected[i], "permutation output mismatch at index {}", i);
        }
    }

    #[test]
    fn test_hash_consistency() {
        // Input: 4 copies of 0x9a807b615c4d3e2fa0b1c2d3e4f56789fedcba9876543210abcdef0123456789
        let val = Fr::from_limbs([
            0xabcdef0123456789,
            0xfedcba9876543210,
            0xa0b1c2d3e4f56789,
            0x9a807b615c4d3e2f,
        ]);
        // Note: In C++ test, a and b lack "0x" prefix but c and d have it.
        // BB's fr(string) treats both the same way - as hex. So all 4 are identical.
        let input = [val, val, val, val];
        let result = Poseidon2::hash(&input);

        let expected = Fr::from_limbs([
            0x41930a25a0bcec11,
            0x47637802a579fa98,
            0xfc839dea0ecec749,
            0x2f43a0f83b51a6f5,
        ]);
        assert_eq!(result, expected, "hash consistency check failed");
    }

    #[test]
    fn test_hash_determinism() {
        let a = Fr::from(42u64);
        let b = Fr::from(99u64);
        let c = Fr::from(1337u64);
        let d = Fr::from(7u64);

        let input1 = [a, b, c, d];
        let input2 = [d, c, b, a];

        let r0 = Poseidon2::hash(&input1);
        let r1 = Poseidon2::hash(&input1);
        let r2 = Poseidon2::hash(&input2);

        assert_eq!(r0, r1, "same input must produce same output");
        assert_ne!(r0, r2, "different input order must produce different output");
    }

    #[test]
    fn test_hash_sensitivity() {
        let a = Fr::from(1u64);
        let b = Fr::from(2u64);
        let input1 = [a, b];
        let input2 = [a, b + Fr::from(1u64)];

        let r0 = Poseidon2::hash(&input1);
        let r1 = Poseidon2::hash(&input2);
        assert_ne!(r0, r1, "small input change must produce different output");
    }
}
