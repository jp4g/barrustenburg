// Random number generation.
//
// C++ source: barretenberg/cpp/src/barretenberg/numeric/random/engine.{hpp,cpp}
//
// BB provides two RNG implementations:
// - RandomEngine: CSPRNG backed by OS entropy (getrandom/getentropy)
// - DebugEngine: Deterministic MT19937 seeded for reproducibility
//
// We delegate to the `rand` crate which provides equivalent functionality.

use rand::{rngs::StdRng, Rng, SeedableRng};

use crate::U256;
use crate::uint256::U256Ext;

/// Get a random u64 from OS entropy.
pub fn get_random_u64() -> u64 {
    rand::rng().random()
}

/// Get a random U256 from OS entropy.
pub fn get_random_u256() -> U256 {
    let mut rng = rand::rng();
    U256::from_limbs([rng.random(), rng.random(), rng.random(), rng.random()])
}

/// Deterministic RNG for testing, seeded from a u64.
///
/// Mirrors BB's `DebugEngine` backed by MT19937.
pub struct DebugRng {
    inner: StdRng,
}

impl DebugRng {
    pub fn new(seed: u64) -> Self {
        Self { inner: StdRng::seed_from_u64(seed) }
    }

    pub fn get_random_u64(&mut self) -> u64 {
        self.inner.random()
    }

    pub fn get_random_u256(&mut self) -> U256 {
        U256::from_limbs([
            self.inner.random(),
            self.inner.random(),
            self.inner.random(),
            self.inner.random(),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn debug_rng_is_deterministic() {
        let mut rng1 = DebugRng::new(42);
        let mut rng2 = DebugRng::new(42);
        for _ in 0..10 {
            assert_eq!(rng1.get_random_u64(), rng2.get_random_u64());
        }
    }

    #[test]
    fn debug_rng_different_seeds_differ() {
        let mut rng1 = DebugRng::new(1);
        let mut rng2 = DebugRng::new(2);
        // Overwhelmingly likely to differ
        assert_ne!(rng1.get_random_u64(), rng2.get_random_u64());
    }

    #[test]
    fn os_rng_produces_values() {
        let a = get_random_u64();
        let b = get_random_u64();
        // Extremely unlikely to be equal
        assert_ne!(a, b);
    }
}
