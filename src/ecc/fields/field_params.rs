/// Trait defining the parameters for a prime field in Montgomery form.
///
/// All constants use 4 x u64 limbs in little-endian order.
/// Montgomery form: elements are stored as `a * R mod p` where `R = 2^256`.
pub trait FieldParams: 'static + Send + Sync + Sized {
    /// The prime modulus p, split into 4 little-endian 64-bit limbs.
    const MODULUS: [u64; 4];

    /// R^2 mod p, used to convert into Montgomery form.
    const R_SQUARED: [u64; 4];

    /// -(p^{-1}) mod 2^64, used in Montgomery reduction.
    const R_INV: u64;

    /// Cube root of unity in Montgomery form (for endomorphism support).
    /// Zero if not applicable.
    const CUBE_ROOT: [u64; 4];

    /// Primitive root of unity in Montgomery form (for FFT).
    /// Zero if not applicable.
    const PRIMITIVE_ROOT: [u64; 4];

    /// 8 coset generators, each as 4 limbs, in Montgomery form.
    const COSET_GENERATORS_0: [u64; 8];
    const COSET_GENERATORS_1: [u64; 8];
    const COSET_GENERATORS_2: [u64; 8];
    const COSET_GENERATORS_3: [u64; 8];

    /// Whether the modulus >= 2^254. Controls which add/sub/mul path is used.
    const MODULUS_IS_BIG: bool;

    /// Whether this field supports large FFT-friendly subgroups (high 2-adicity).
    /// True for BN254 Fr (2-adicity = 28). Default false.
    const HAS_HIGH_2ADICITY: bool = false;

    // --- GLV endomorphism constants (defaults to 0, override for curves with endomorphism) ---
    const ENDO_G1_LO: u64 = 0;
    const ENDO_G1_MID: u64 = 0;
    const ENDO_G1_HI: u64 = 0;
    const ENDO_G1_HIHI: u64 = 0;
    const ENDO_G2_LO: u64 = 0;
    const ENDO_G2_MID: u64 = 0;
    const ENDO_G2_HI: u64 = 0;
    const ENDO_G2_HIHI: u64 = 0;
    const ENDO_MINUS_B1_LO: u64 = 0;
    const ENDO_MINUS_B1_MID: u64 = 0;
    const ENDO_B2_LO: u64 = 0;
    const ENDO_B2_MID: u64 = 0;
}
