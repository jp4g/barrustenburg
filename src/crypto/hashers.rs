/// Trait for hash functions used by HMAC, ECDSA, and Schnorr.
///
/// Mirrors the C++ Hash concept: BLOCK_SIZE, OUTPUT_SIZE, and a hash function.
pub trait Hasher {
    const BLOCK_SIZE: usize;
    const OUTPUT_SIZE: usize;
    fn hash(message: &[u8]) -> Vec<u8>;
}

/// SHA-256 hasher.
pub struct Sha256Hasher;

impl Hasher for Sha256Hasher {
    const BLOCK_SIZE: usize = 64;
    const OUTPUT_SIZE: usize = 32;

    fn hash(message: &[u8]) -> Vec<u8> {
        crate::crypto::sha256::sha256(message).to_vec()
    }
}

/// BLAKE2s-256 hasher.
pub struct Blake2sHasher;

impl Hasher for Blake2sHasher {
    const BLOCK_SIZE: usize = 64;
    const OUTPUT_SIZE: usize = 32;

    fn hash(message: &[u8]) -> Vec<u8> {
        crate::crypto::blake2s::blake2s(message).to_vec()
    }
}

/// Keccak-256 hasher (Ethereum variant).
pub struct KeccakHasher;

impl Hasher for KeccakHasher {
    const BLOCK_SIZE: usize = 64;
    const OUTPUT_SIZE: usize = 32;

    fn hash(message: &[u8]) -> Vec<u8> {
        crate::crypto::keccak::keccak256(message).to_bytes().to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_hasher_known_output() {
        let result = Sha256Hasher::hash(b"abc");
        assert_eq!(result.len(), 32);
        assert_eq!(result[0], 0xBA);
        assert_eq!(result[1], 0x78);
    }

    #[test]
    fn blake2s_hasher_known_output() {
        let result = Blake2sHasher::hash(b"abc");
        assert_eq!(result.len(), 32);
        assert_eq!(result[0], 0x50);
        assert_eq!(result[1], 0x8C);
    }

    #[test]
    fn keccak_hasher_known_output() {
        let result = KeccakHasher::hash(b"");
        assert_eq!(result.len(), 32);
        assert_eq!(result[0], 0xc5);
        assert_eq!(result[1], 0xd2);
    }
}
