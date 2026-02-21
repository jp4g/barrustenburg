use crate::crypto::hashers::Hasher;
use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;

/// HMAC (RFC 2104) generic over a hash function.
///
/// Matches C++ `crypto::hmac<Hash>(message, key)`.
pub fn hmac<H: Hasher>(message: &[u8], key: &[u8]) -> Vec<u8> {
    // If key > block size, hash it first
    let k_prime = if key.len() > H::BLOCK_SIZE {
        H::hash(key)
    } else {
        key.to_vec()
    };

    // Pad key to block size with zeros
    let mut k_padded = vec![0u8; H::BLOCK_SIZE];
    k_padded[..k_prime.len()].copy_from_slice(&k_prime);

    // ipad = k_padded XOR 0x36
    let mut ipad = vec![0u8; H::BLOCK_SIZE];
    for i in 0..H::BLOCK_SIZE {
        ipad[i] = k_padded[i] ^ 0x36;
    }

    // opad = k_padded XOR 0x5c
    let mut opad = vec![0u8; H::BLOCK_SIZE];
    for i in 0..H::BLOCK_SIZE {
        opad[i] = k_padded[i] ^ 0x5c;
    }

    // inner = H(ipad || message)
    let mut inner_input = Vec::with_capacity(H::BLOCK_SIZE + message.len());
    inner_input.extend_from_slice(&ipad);
    inner_input.extend_from_slice(message);
    let inner = H::hash(&inner_input);

    // outer = H(opad || inner)
    let mut outer_input = Vec::with_capacity(H::BLOCK_SIZE + H::OUTPUT_SIZE);
    outer_input.extend_from_slice(&opad);
    outer_input.extend_from_slice(&inner);
    H::hash(&outer_input)
}

/// Derive an unbiased field element from HMAC output.
///
/// Used by ECDSA for deterministic nonce generation.
/// Matches C++ `crypto::get_unbiased_field_from_hmac<Hash, Fr>`.
pub fn get_unbiased_field_from_hmac<H: Hasher, P: FieldParams>(
    message: &[u8],
    key: &[u8],
) -> Field<P> {
    let hmac_output = hmac::<H>(message, key);

    // Domain separators: BLOCK_SIZE - OUTPUT_SIZE bytes
    let separator_size = H::BLOCK_SIZE - H::OUTPUT_SIZE;

    // KLO domain separator: all zeros
    let mut klo_input = vec![0u8; separator_size + H::OUTPUT_SIZE];
    klo_input[separator_size..].copy_from_slice(&hmac_output);
    let klo = H::hash(&klo_input);

    // KHI domain separator: first byte is 0x01, rest zeros
    let mut khi_input = vec![0u8; separator_size + H::OUTPUT_SIZE];
    khi_input[0] = 0x01;
    khi_input[separator_size..].copy_from_slice(&hmac_output);
    let khi = H::hash(&khi_input);

    // Read klo and khi as big-endian u256 limbs
    let read_be_limbs = |bytes: &[u8]| -> [u64; 4] {
        let d3 = u64::from_be_bytes(bytes[0..8].try_into().unwrap());
        let d2 = u64::from_be_bytes(bytes[8..16].try_into().unwrap());
        let d1 = u64::from_be_bytes(bytes[16..24].try_into().unwrap());
        let d0 = u64::from_be_bytes(bytes[24..32].try_into().unwrap());
        [d0, d1, d2, d3]
    };

    let lo_limbs = read_be_limbs(&klo);
    let hi_limbs = read_be_limbs(&khi);

    // Reduce 512-bit (hi || lo) mod p
    Field::from_u512(lo_limbs, hi_limbs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crypto::hashers::Sha256Hasher;
    use crate::ecc::curves::bn254::Bn254FrParams;

    /// RFC 4231 Test Case 1: HMAC-SHA256
    #[test]
    fn hmac_sha256_rfc4231_test_case_1() {
        let key = vec![0x0bu8; 20];
        let data = b"Hi There";
        let result = hmac::<Sha256Hasher>(data, &key);
        let expected: [u8; 32] = [
            0xb0, 0x34, 0x4c, 0x61, 0xd8, 0xdb, 0x38, 0x53,
            0x5c, 0xa8, 0xaf, 0xce, 0xaf, 0x0b, 0xf1, 0x2b,
            0x88, 0x1d, 0xc2, 0x00, 0xc9, 0x83, 0x3d, 0xa7,
            0x26, 0xe9, 0x37, 0x6c, 0x2e, 0x32, 0xcf, 0xf7,
        ];
        assert_eq!(result.as_slice(), &expected);
    }

    /// RFC 4231 Test Case 2: HMAC-SHA256 with "Jefe" key
    #[test]
    fn hmac_sha256_rfc4231_test_case_2() {
        let key = b"Jefe";
        let data = b"what do ya want for nothing?";
        let result = hmac::<Sha256Hasher>(data, key);
        let expected: [u8; 32] = [
            0x5b, 0xdc, 0xc1, 0x46, 0xbf, 0x60, 0x75, 0x4e,
            0x6a, 0x04, 0x24, 0x26, 0x08, 0x95, 0x75, 0xc7,
            0x5a, 0x00, 0x3f, 0x08, 0x9d, 0x27, 0x39, 0x83,
            0x9d, 0xec, 0x58, 0xb9, 0x64, 0xec, 0x38, 0x43,
        ];
        assert_eq!(result.as_slice(), &expected);
    }

    /// RFC 4231 Test Case 3: HMAC-SHA256 with 0xaa key
    #[test]
    fn hmac_sha256_rfc4231_test_case_3() {
        let key = vec![0xaau8; 20];
        let data = vec![0xddu8; 50];
        let result = hmac::<Sha256Hasher>(&data, &key);
        let expected: [u8; 32] = [
            0x77, 0x3e, 0xa9, 0x1e, 0x36, 0x80, 0x0e, 0x46,
            0x85, 0x4d, 0xb8, 0xeb, 0xd0, 0x91, 0x81, 0xa7,
            0x29, 0x59, 0x09, 0x8b, 0x3e, 0xf8, 0xc1, 0x22,
            0xd9, 0x63, 0x55, 0x14, 0xce, 0xd5, 0x65, 0xfe,
        ];
        assert_eq!(result.as_slice(), &expected);
    }

    #[test]
    fn get_unbiased_field_from_hmac_deterministic() {
        let msg = b"test message";
        let key = b"test key";
        let r1 = get_unbiased_field_from_hmac::<Sha256Hasher, Bn254FrParams>(msg, key);
        let r2 = get_unbiased_field_from_hmac::<Sha256Hasher, Bn254FrParams>(msg, key);
        assert_eq!(r1, r2, "same inputs should produce same field element");
        assert!(!r1.is_zero(), "result should not be zero");
    }
}
