use crate::ecc::curves::bn254::Fr;

use super::permutation;

const RATE: usize = 3;
const T: usize = 4; // rate + capacity

/// Poseidon2 sponge hash over BN254 Fr.
///
/// Implements the sponge specification with rate=3, capacity=1.
/// The capacity element (state[3]) is initialized with input_length << 64.
pub fn hash(input: &[Fr]) -> Fr {
    // IV = input_length << 64
    let iv = {
        let len = input.len() as u64;
        // The IV is (len << 64) as a field element.
        // In [u64;4] limbs (little-endian): [0, len, 0, 0]
        Fr::from_limbs([0, len, 0, 0])
    };

    hash_with_iv(input, iv)
}

/// Hash with a custom initial value (IV).
fn hash_with_iv(input: &[Fr], iv: Fr) -> Fr {
    let mut state = [Fr::zero(); T];
    state[RATE] = iv; // capacity element

    let mut cache = [Fr::zero(); RATE];
    let mut cache_size = 0usize;

    // Absorb all input elements
    for &elem in input {
        if cache_size == RATE {
            // Cache full: duplex
            for i in 0..RATE {
                state[i] = state[i] + cache[i];
            }
            state = permutation::permutation(&state);
            cache = [Fr::zero(); RATE];
            cache[0] = elem;
            cache_size = 1;
        } else {
            cache[cache_size] = elem;
            cache_size += 1;
        }
    }

    // Squeeze: final duplex
    for i in 0..RATE {
        state[i] = state[i] + cache[i];
    }
    state = permutation::permutation(&state);

    state[0]
}
