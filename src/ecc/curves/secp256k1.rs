use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;
use crate::ecc::groups::curve_params::CurveParams;

// ---------------------------------------------------------------------------
// secp256k1 Base Field (Fq)
// ---------------------------------------------------------------------------

pub struct Secp256k1FqParams;

impl FieldParams for Secp256k1FqParams {
    const MODULUS: [u64; 4] = [
        0xFFFFFFFEFFFFFC2F,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
    ];
    const R_SQUARED: [u64; 4] = [8392367050913, 1, 0, 0];
    const R_INV: u64 = 15580212934572586289;
    const CUBE_ROOT: [u64; 4] = [
        0x58a4361c8e81894e,
        0x03fde1631c4b80af,
        0xf8e98978d02e3905,
        0x7a4a36aebcbb3d53,
    ];
    const PRIMITIVE_ROOT: [u64; 4] = [0, 0, 0, 0];
    const COSET_GENERATORS_0: [u64; 8] = [
        0x300000b73, 0x400000f44, 0x500001315, 0x6000016e6,
        0x700001ab7, 0x800001e88, 0x900002259, 0xa0000262a,
    ];
    const COSET_GENERATORS_1: [u64; 8] = [0; 8];
    const COSET_GENERATORS_2: [u64; 8] = [0; 8];
    const COSET_GENERATORS_3: [u64; 8] = [0; 8];
    const MODULUS_IS_BIG: bool = true; // 0xFFFF... >= 0x4000...
}

pub type Secp256k1Fq = Field<Secp256k1FqParams>;

// ---------------------------------------------------------------------------
// secp256k1 Scalar Field (Fr)
// ---------------------------------------------------------------------------

pub struct Secp256k1FrParams;

impl FieldParams for Secp256k1FrParams {
    const MODULUS: [u64; 4] = [
        0xBFD25E8CD0364141,
        0xBAAEDCE6AF48A03B,
        0xFFFFFFFFFFFFFFFE,
        0xFFFFFFFFFFFFFFFF,
    ];
    const R_SQUARED: [u64; 4] = [
        9902555850136342848,
        8364476168144746616,
        16616019711348246470,
        11342065889886772165,
    ];
    const R_INV: u64 = 5408259542528602431;
    const CUBE_ROOT: [u64; 4] = [
        0xf07deb3dc9926c9e,
        0x2c93e7ad83c6944c,
        0x73a9660652697d91,
        0x532840178558d639,
    ];
    const PRIMITIVE_ROOT: [u64; 4] = [0, 0, 0, 0];
    const COSET_GENERATORS_0: [u64; 8] = [
        0x40e4273feef0b9bb, 0x8111c8b31eba787a, 0xc13f6a264e843739, 0x16d0b997e4df5f8,
        0x419aad0cae17b4b7, 0x81c84e7fdde17376, 0xc1f5eff30dab3235, 0x22391663d74f0f4,
    ];
    const COSET_GENERATORS_1: [u64; 8] = [
        0x5a95af7e9394ded5, 0x9fe6d297e44c3e99, 0xe537f5b135039e5d, 0x2a8918ca85bafe22,
        0x6fda3be3d6725de6, 0xb52b5efd2729bdaa, 0xfa7c821677e11d6e, 0x3fcda52fc8987d33,
    ];
    const COSET_GENERATORS_2: [u64; 8] = [
        0x6, 0x7, 0x8, 0xa, 0xb, 0xc, 0xd, 0xf,
    ];
    const COSET_GENERATORS_3: [u64; 8] = [0; 8];
    const MODULUS_IS_BIG: bool = true; // 0xFFFF... >= 0x4000...
}

pub type Secp256k1Fr = Field<Secp256k1FrParams>;

// ---------------------------------------------------------------------------
// secp256k1 G1 Curve Parameters
// ---------------------------------------------------------------------------

/// secp256k1 G1: y^2 = x^3 + 7.
/// Generator coordinates stored in standard (non-Montgomery) form,
/// matching C++ which calls `.to_montgomery_form()` at runtime.
pub struct Secp256k1G1Params;

impl CurveParams for Secp256k1G1Params {
    type BaseFieldParams = Secp256k1FqParams;
    type ScalarFieldParams = Secp256k1FrParams;

    const HAS_A: bool = false;
    const A: [u64; 4] = [0, 0, 0, 0];

    /// b = 7 (placeholder, overridden by coeff_b).
    const B: [u64; 4] = [7, 0, 0, 0];

    /// Generator x in standard form.
    const GENERATOR_X: [u64; 4] = [
        0x59F2815B16F81798,
        0x029BFCDB2DCE28D9,
        0x55A06295CE870B07,
        0x79BE667EF9DCBBAC,
    ];

    /// Generator y in standard form.
    const GENERATOR_Y: [u64; 4] = [
        0x9C47D08FFB10D4B8,
        0xFD17B448A6855419,
        0x5DA4FBFC0E1108A8,
        0x483ADA7726A3C465,
    ];

    const USE_ENDOMORPHISM: bool = true;

    /// Convert from standard form at runtime, matching C++.
    fn generator_x() -> Secp256k1Fq {
        Field::from_limbs(Self::GENERATOR_X)
    }
    fn generator_y() -> Secp256k1Fq {
        Field::from_limbs(Self::GENERATOR_Y)
    }
    fn coeff_b() -> Secp256k1Fq {
        Field::from(7u64)
    }
}

pub type G1Affine = crate::ecc::groups::affine_element::AffineElement<Secp256k1G1Params>;
pub type G1Element = crate::ecc::groups::element::Element<Secp256k1G1Params>;
