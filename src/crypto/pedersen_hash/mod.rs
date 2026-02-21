use crate::crypto::generators::precomputed::LENGTH_GENERATOR;
use crate::crypto::generators::GeneratorContext;
use crate::crypto::pedersen_commitment;
use crate::ecc::curves::grumpkin::{self, GrumpkinFq, GrumpkinFr};

type GrumpkinElement = grumpkin::G1Element;

/// Pedersen hash over Grumpkin.
///
/// hash(inputs) = (len * length_generator + commit(inputs)).to_affine().x
///
/// The length inclusion prevents length-extension attacks.
/// Output is the x-coordinate only.
pub fn hash(inputs: &[GrumpkinFq], ctx: &GeneratorContext) -> GrumpkinFq {
    // len * length_generator
    // For Grumpkin: ScalarField = BN254 Fq = GrumpkinFr
    let len_scalar = GrumpkinFr::from(inputs.len() as u64);
    let length_gen = GrumpkinElement::from_affine(&*LENGTH_GENERATOR);
    let result = length_gen.mul_without_endomorphism(&len_scalar);

    // Add commitment
    let commitment = pedersen_commitment::commit_native(inputs, ctx);
    let combined = result + commitment;

    combined.to_affine().x
}

/// Hash a byte buffer using Pedersen.
///
/// Splits input into 31-byte big-endian field elements, then iteratively
/// hashes pairs: hash([x0, x1]), then hash([result, x2]), etc.
pub fn hash_buffer(input: &[u8], ctx: &GeneratorContext) -> GrumpkinFq {
    let converted = convert_buffer(input);

    if converted.len() < 2 {
        return hash(&converted, ctx);
    }

    let mut result = hash(&[converted[0], converted[1]], ctx);
    for i in 2..converted.len() {
        result = hash(&[result, converted[i]], ctx);
    }
    result
}

/// Convert a byte buffer into 31-byte-chunk field elements (big-endian).
fn convert_buffer(input: &[u8]) -> Vec<GrumpkinFq> {
    let num_bytes = input.len();
    let bytes_per_element = 31;
    let num_elements = (num_bytes / bytes_per_element)
        + if num_bytes % bytes_per_element != 0 {
            1
        } else {
            0
        };

    let slice = |data: &[u8], start: usize, size: usize| -> GrumpkinFq {
        // Read big-endian bytes into [u64; 4] limbs
        let mut limbs = [0u64; 4];
        for i in 0..size {
            let byte_pos = size - 1 - i;
            let limb_idx = byte_pos / 8;
            let bit_offset = (byte_pos % 8) * 8;
            limbs[limb_idx] |= (data[start + i] as u64) << bit_offset;
        }
        GrumpkinFq::from_limbs(limbs)
    };

    let mut elements = Vec::with_capacity(num_elements);
    for i in 0..num_elements.saturating_sub(1) {
        elements.push(slice(input, i * bytes_per_element, bytes_per_element));
    }
    let last_start = (num_elements - 1) * bytes_per_element;
    let bytes_remaining = num_bytes - last_start;
    elements.push(slice(input, last_start, bytes_remaining));

    elements
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash_one_one() {
        let x = GrumpkinFq::one();
        let result = hash(&[x, x], &GeneratorContext::default());

        let expected = GrumpkinFq::from_limbs([
            0x378ab8845463297b,
            0xa923013ddbbcbdc3,
            0x6cd6dca13d4bb9d1,
            0x07ebfbf4df29888c,
        ]);

        assert_eq!(result, expected, "hash([1, 1]) mismatch");
    }

    #[test]
    fn test_hash_with_offset() {
        let x = GrumpkinFq::one();
        let ctx = GeneratorContext::with_offset(5);
        let result = hash(&[x, x], &ctx);

        let expected = GrumpkinFq::from_limbs([
            0x17aab70d7861daa6,
            0x6df0cec333fad876,
            0xcda124524e6b03f3,
            0x1c446df60816b897,
        ]);

        assert_eq!(result, expected, "hash([1, 1], offset=5) mismatch");
    }
}
