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
        let val = Fr::from_limbs([
            0xabcdef0123456789,
            0xfedcba9876543210,
            0xa0b1c2d3e4f56789,
            0x9a807b615c4d3e2f,
        ]);
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

    #[test]
    fn test_permutation_determinism_and_sensitivity() {
        let a = Fr::random_element();
        let b = Fr::random_element();
        let c = Fr::random_element();
        let d = Fr::random_element();

        let input1 = [a, b, c, d];
        let input2 = [d, c, b, a];

        let r0 = permutation::permutation(&input1);
        let r1 = permutation::permutation(&input1);
        let r2 = permutation::permutation(&input2);

        assert_eq!(r0, r1, "same input must produce same permutation output");
        assert_ne!(r0, r2, "different input order must produce different permutation output");
    }

    #[test]
    fn test_permutation_consistency_check() {
        // C++ ConsistencyCheck: all 4 inputs = 0x9a807b615c4d3e2fa0b1c2d3e4f56789fedcba9876543210abcdef0123456789
        let val = Fr::from_limbs([
            0xabcdef0123456789,
            0xfedcba9876543210,
            0xa0b1c2d3e4f56789,
            0x9a807b615c4d3e2f,
        ]);
        let input = [val, val, val, val];
        let result = permutation::permutation(&input);

        // C++ expected values from hex strings (standard form, will be converted to Montgomery by from_limbs)
        let expected = [
            // 0x2bf1eaf87f7d27e8dc4056e9af975985bccc89077a21891d6c7b6ccce0631f95
            Fr::from_limbs([0x6c7b6ccce0631f95, 0xbccc89077a21891d, 0xdc4056e9af975985, 0x2bf1eaf87f7d27e8]),
            // 0x0c01fa1b8d0748becafbe452c0cb0231c38224ea824554c9362518eebdd5701f
            Fr::from_limbs([0x362518eebdd5701f, 0xc38224ea824554c9, 0xcafbe452c0cb0231, 0x0c01fa1b8d0748be]),
            // 0x018555a8eb50cf07f64b019ebaf3af3c925c93e631f3ecd455db07bbb52bbdd3
            Fr::from_limbs([0x55db07bbb52bbdd3, 0x925c93e631f3ecd4, 0xf64b019ebaf3af3c, 0x018555a8eb50cf07]),
            // 0x0cbea457c91c22c6c31fd89afd2541efc2edf31736b9f721e823b2165c90fd41
            Fr::from_limbs([0xe823b2165c90fd41, 0xc2edf31736b9f721, 0xc31fd89afd2541ef, 0x0cbea457c91c22c6]),
        ];

        for i in 0..4 {
            assert_eq!(result[i], expected[i], "permutation consistency mismatch at index {}", i);
        }
    }
}
