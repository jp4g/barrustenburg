pub mod precomputed;

use bbrs_ecc::curves::grumpkin::{self, GrumpkinFq};
use bbrs_ecc::groups::affine_element::AffineElement;

type GrumpkinAffine = grumpkin::G1Affine;

/// Generator context for Pedersen operations.
pub struct GeneratorContext {
    pub offset: usize,
    pub domain_separator: String,
}

impl Default for GeneratorContext {
    fn default() -> Self {
        Self {
            offset: 0,
            domain_separator: "DEFAULT_DOMAIN_SEPARATOR".to_string(),
        }
    }
}

impl GeneratorContext {
    pub fn new(offset: usize, domain_separator: &str) -> Self {
        Self {
            offset,
            domain_separator: domain_separator.to_string(),
        }
    }

    pub fn with_offset(offset: usize) -> Self {
        Self {
            offset,
            ..Default::default()
        }
    }

    /// Get generators for this context.
    pub fn get_generators(&self, count: usize) -> Vec<GrumpkinAffine> {
        let is_default = self.domain_separator == "DEFAULT_DOMAIN_SEPARATOR";
        let total_needed = count + self.offset;

        // Use precomputed generators if possible
        if is_default && total_needed <= 8 {
            let precomputed = &*precomputed::DEFAULT_GENERATORS;
            return precomputed[self.offset..self.offset + count].to_vec();
        }

        // Otherwise, derive generators
        let derived = derive_generators(
            self.domain_separator.as_bytes(),
            total_needed,
            0,
        );
        derived[self.offset..self.offset + count].to_vec()
    }
}

/// Derive `num_generators` independent Grumpkin curve points via BLAKE3 hash-to-curve.
///
/// Algorithm (matching C++ `group::derive_generators`):
/// 1. Hash domain separator with BLAKE3 to get 32 bytes
/// 2. Build a 64-byte preimage: [hash(32 bytes), generator_index(4 bytes big-endian), zeros(28 bytes)]
/// 3. For each generator, call hash_to_curve with the preimage
pub fn derive_generators(
    domain_separator: &[u8],
    num_generators: usize,
    starting_index: usize,
) -> Vec<GrumpkinAffine> {
    let domain_hash = blake3::hash(domain_separator);
    let domain_hash_bytes = domain_hash.as_bytes();

    let mut generator_preimage = vec![0u8; 64];
    generator_preimage[..32].copy_from_slice(domain_hash_bytes);

    let mut result = Vec::with_capacity(num_generators);

    for i in starting_index..starting_index + num_generators {
        let idx = i as u32;
        // Big-endian u32 at bytes 32..36
        generator_preimage[32] = (idx >> 24) as u8;
        generator_preimage[33] = ((idx >> 16) & 0xff) as u8;
        generator_preimage[34] = ((idx >> 8) & 0xff) as u8;
        generator_preimage[35] = (idx & 0xff) as u8;
        // Bytes 36..64 remain zero (reset not needed since we initialized to 0)

        result.push(hash_to_curve(&generator_preimage));
    }

    result
}

/// Hash a seed buffer into a Grumpkin curve point.
///
/// Algorithm (matching C++ `affine_element::hash_to_curve`):
/// 1. Append 2 bytes to seed: [attempt_count, 0x00]
/// 2. BLAKE3(seed || attempt_count || 0x00) -> hash_hi (32 bytes)
/// 3. BLAKE3(seed || attempt_count || 0x01) -> hash_lo (32 bytes)
/// 4. Interpret as 512-bit integer (hi || lo), reduce mod Fq -> x-coordinate
/// 5. Try to derive y from curve equation
/// 6. If not a quadratic residue, increment attempt_count and retry
/// 7. Use MSB of hash_hi as sign bit for y-parity
fn hash_to_curve(seed: &[u8]) -> GrumpkinAffine {
    let seed_size = seed.len();
    let mut target_seed = Vec::with_capacity(seed_size + 2);
    target_seed.extend_from_slice(seed);
    target_seed.push(0); // attempt_count
    target_seed.push(0); // hash selector

    for attempt_count in 0u8..=255 {
        target_seed[seed_size] = attempt_count;

        // Hash with selector = 0 -> hi
        target_seed[seed_size + 1] = 0;
        let hash_hi = blake3::hash(&target_seed);
        let hash_hi_bytes = hash_hi.as_bytes();

        // Hash with selector = 1 -> lo
        target_seed[seed_size + 1] = 1;
        let hash_lo = blake3::hash(&target_seed);
        let hash_lo_bytes = hash_lo.as_bytes();

        // Read each 32-byte hash as a big-endian uint256
        let hi_limbs = read_uint256_be(hash_hi_bytes);
        let lo_limbs = read_uint256_be(hash_lo_bytes);

        // Reduce 512-bit (hi || lo) mod Fq to get x-coordinate
        // lo is the lower 256 bits, hi is the upper 256 bits
        let x = GrumpkinFq::from_u512(lo_limbs, hi_limbs);

        // Sign bit = MSB of hash_hi
        let sign_bit = hash_hi_bytes[0] > 127;

        if let Some(point) = AffineElement::from_x_coordinate(x, sign_bit) {
            return point;
        }
    }

    panic!("hash_to_curve: failed to find a valid point after 256 attempts");
}

/// Read 32 bytes (big-endian) as [u64; 4] in little-endian limb order.
fn read_uint256_be(bytes: &[u8; 32]) -> [u64; 4] {
    let read_limb = |offset: usize| -> u64 {
        let mut val = 0u64;
        for i in 0..8 {
            val |= (bytes[offset + i] as u64) << ((7 - i) * 8);
        }
        val
    };

    // bytes[0..8] = most significant limb (data[3])
    // bytes[8..16] = data[2]
    // bytes[16..24] = data[1]
    // bytes[24..32] = data[0] (least significant)
    [read_limb(24), read_limb(16), read_limb(8), read_limb(0)]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_precomputed_generators_on_curve() {
        let generators = &*precomputed::DEFAULT_GENERATORS;
        for (i, g) in generators.iter().enumerate() {
            assert!(g.on_curve(), "precomputed generator {} is not on curve", i);
        }
    }

    #[test]
    fn test_length_generator_on_curve() {
        let g = *precomputed::LENGTH_GENERATOR;
        assert!(g.on_curve(), "length generator is not on curve");
    }

    #[test]
    fn test_length_generator_matches() {
        let expected = *precomputed::LENGTH_GENERATOR;
        let derived = derive_generators(b"pedersen_hash_length", 1, 0);
        assert_eq!(derived[0], expected, "length generator does not match derive_generators");
    }

    #[test]
    fn test_default_generators_match_derived() {
        let precomputed = &*precomputed::DEFAULT_GENERATORS;
        let derived = derive_generators(b"DEFAULT_DOMAIN_SEPARATOR", 8, 0);

        for i in 0..8 {
            assert_eq!(
                precomputed[i], derived[i],
                "precomputed generator {} does not match derived", i
            );
        }
    }

    #[test]
    fn test_derived_generators_on_curve() {
        let derived = derive_generators(b"test_domain", 4, 0);
        for (i, g) in derived.iter().enumerate() {
            assert!(g.on_curve(), "derived generator {} is not on curve", i);
        }
    }
}
