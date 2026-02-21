use sha2::{compress256, Digest, Sha256 as Sha256Hasher};

/// Full SHA-256 hash.
pub fn sha256(input: &[u8]) -> [u8; 32] {
    let result = Sha256Hasher::digest(input);
    result.into()
}

/// SHA-256 compression function (FIPS 180-4 Section 6.2.2).
/// Required by ACIR Sha256Compression opcode.
pub fn sha256_block(h_init: &[u32; 8], input: &[u32; 16]) -> [u32; 8] {
    // Convert [u32; 16] to [u8; 64] in big-endian byte order
    let mut block = [0u8; 64];
    for (i, word) in input.iter().enumerate() {
        block[i * 4..i * 4 + 4].copy_from_slice(&word.to_be_bytes());
    }

    let mut state = *h_init;
    compress256(&mut state, &[block.into()]);
    state
}

#[cfg(test)]
mod tests {
    use super::*;

    // NIST test vector 1: "abc"
    #[test]
    fn test_nist_vector_one() {
        let input = b"abc";
        let result = sha256(input);
        let expected: [u8; 32] = [
            0xBA, 0x78, 0x16, 0xBF, 0x8F, 0x01, 0xCF, 0xEA, 0x41, 0x41, 0x40, 0xDE, 0x5D, 0xAE,
            0x22, 0x23, 0xB0, 0x03, 0x61, 0xA3, 0x96, 0x17, 0x7A, 0x9C, 0xB4, 0x10, 0xFF, 0x61,
            0xF2, 0x00, 0x15, 0xAD,
        ];
        assert_eq!(result, expected);
    }

    // NIST test vector 2: two-block message
    #[test]
    fn test_nist_vector_two() {
        let input = b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
        let result = sha256(input);
        let expected: [u8; 32] = [
            0x24, 0x8D, 0x6A, 0x61, 0xD2, 0x06, 0x38, 0xB8, 0xE5, 0xC0, 0x26, 0x93, 0x0C, 0x3E,
            0x60, 0x39, 0xA3, 0x3C, 0xE4, 0x59, 0x64, 0xFF, 0x21, 0x67, 0xF6, 0xEC, 0xED, 0xD4,
            0x19, 0xDB, 0x06, 0xC1,
        ];
        assert_eq!(result, expected);
    }

    // NIST test vector 3: single byte 0xbd
    #[test]
    fn test_nist_vector_three() {
        let input = [0xbd];
        let result = sha256(&input);
        let expected: [u8; 32] = [
            0x68, 0x32, 0x57, 0x20, 0xaa, 0xbd, 0x7c, 0x82, 0xf3, 0x0f, 0x55, 0x4b, 0x31, 0x3d,
            0x05, 0x70, 0xc9, 0x5a, 0xcc, 0xbb, 0x7d, 0xc4, 0xb5, 0xaa, 0xe1, 0x12, 0x04, 0xc0,
            0x8f, 0xfe, 0x73, 0x2b,
        ];
        assert_eq!(result, expected);
    }

    // NIST test vector 4: 4-byte input
    #[test]
    fn test_nist_vector_four() {
        let input: [u8; 4] = [0xc9, 0x8c, 0x8e, 0x55];
        let result = sha256(&input);
        let expected: [u8; 32] = [
            0x7a, 0xbc, 0x22, 0xc0, 0xae, 0x5a, 0xf2, 0x6c, 0xe9, 0x3d, 0xbb, 0x94, 0x43, 0x3a,
            0x0e, 0x0b, 0x2e, 0x11, 0x9d, 0x01, 0x4f, 0x8e, 0x7f, 0x65, 0xbd, 0x56, 0xc6, 0x1c,
            0xcc, 0xcd, 0x95, 0x04,
        ];
        assert_eq!(result, expected);
    }

    // NIST test vector 5: 1000 bytes of 'A'
    #[test]
    fn test_nist_vector_five() {
        let input = vec![b'A'; 1000];
        let result = sha256(&input);
        let expected: [u8; 32] = [
            0xc2, 0xe6, 0x86, 0x82, 0x34, 0x89, 0xce, 0xd2, 0x01, 0x7f, 0x60, 0x59, 0xb8, 0xb2,
            0x39, 0x31, 0x8b, 0x63, 0x64, 0xf6, 0xdc, 0xd8, 0x35, 0xd0, 0xa5, 0x19, 0x10, 0x5a,
            0x1e, 0xad, 0xd6, 0xe4,
        ];
        assert_eq!(result, expected);
    }

    // Compression function test: verify sha256_block matches the first block compression
    // of sha256("abc") by manually constructing the padded block.
    #[test]
    fn test_sha256_compression() {
        // SHA-256 initial hash values (FIPS 180-4 Section 5.3.3)
        let h_init: [u32; 8] = [
            0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
            0x5be0cd19,
        ];

        // "abc" padded to 512 bits (big-endian u32 words):
        // 61626380 00000000 00000000 00000000
        // 00000000 00000000 00000000 00000000
        // 00000000 00000000 00000000 00000000
        // 00000000 00000000 00000000 00000018
        let input: [u32; 16] = [
            0x61626380, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
            0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
            0x00000000, 0x00000018,
        ];

        let result = sha256_block(&h_init, &input);

        // The result should match SHA-256("abc") interpreted as big-endian u32 words
        let expected: [u32; 8] = [
            0xBA7816BF, 0x8F01CFEA, 0x414140DE, 0x5DAE2223, 0xB00361A3, 0x96177A9C, 0xB410FF61,
            0xF20015AD,
        ];
        assert_eq!(result, expected);
    }
}
