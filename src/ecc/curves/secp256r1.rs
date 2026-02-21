use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;
use crate::ecc::groups::curve_params::CurveParams;

// ---------------------------------------------------------------------------
// secp256r1 (P-256) Base Field (Fq)
// ---------------------------------------------------------------------------

pub struct Secp256r1FqParams;

impl FieldParams for Secp256r1FqParams {
    const MODULUS: [u64; 4] = [
        0xFFFFFFFFFFFFFFFF,
        0x00000000FFFFFFFF,
        0x0000000000000000,
        0xFFFFFFFF00000001,
    ];
    const R_SQUARED: [u64; 4] = [3, 18446744056529682431, 18446744073709551614, 21474836477];
    const R_INV: u64 = 1;
    const CUBE_ROOT: [u64; 4] = [0, 0, 0, 0]; // Not applicable
    const PRIMITIVE_ROOT: [u64; 4] = [0, 0, 0, 0];
    const COSET_GENERATORS_0: [u64; 8] = [
        0x3, 0x4, 0x5, 0x6, 0x7, 0x8, 0x9, 0xa,
    ];
    const COSET_GENERATORS_1: [u64; 8] = [
        0xfffffffd00000000, 0xfffffffc00000000, 0xfffffffb00000000, 0xfffffffa00000000,
        0xfffffff900000000, 0xfffffff800000000, 0xfffffff700000000, 0xfffffff600000000,
    ];
    const COSET_GENERATORS_2: [u64; 8] = [
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
        0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
    ];
    const COSET_GENERATORS_3: [u64; 8] = [
        0x2fffffffc, 0x3fffffffb, 0x4fffffffa, 0x5fffffff9,
        0x6fffffff8, 0x7fffffff7, 0x8fffffff6, 0x9fffffff5,
    ];
    const MODULUS_IS_BIG: bool = true; // 0xFFFFFFFF00000001 >= 0x4000...
}

pub type Secp256r1Fq = Field<Secp256r1FqParams>;

// ---------------------------------------------------------------------------
// secp256r1 (P-256) Scalar Field (Fr)
// ---------------------------------------------------------------------------

pub struct Secp256r1FrParams;

impl FieldParams for Secp256r1FrParams {
    const MODULUS: [u64; 4] = [
        0xF3B9CAC2FC632551,
        0xBCE6FAADA7179E84,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFF00000000,
    ];
    const R_SQUARED: [u64; 4] = [
        9449762124159643298,
        5087230966250696614,
        2901921493521525849,
        7413256579398063648,
    ];
    const R_INV: u64 = 14758798090332847183;
    const CUBE_ROOT: [u64; 4] = [0, 0, 0, 0];
    const PRIMITIVE_ROOT: [u64; 4] = [0, 0, 0, 0];
    const COSET_GENERATORS_0: [u64; 8] = [
        0x55eb74ab1949fac9, 0x6231a9e81ce6d578, 0x6e77df252083b027, 0x7abe146224208ad6,
        0x8704499f27bd6585, 0x934a7edc2b5a4034, 0x9f90b4192ef71ae3, 0xabd6e9563293f592,
    ];
    const COSET_GENERATORS_1: [u64; 8] = [
        0xd5af25406e5aaa5d, 0x18c82a92c7430bd8, 0x5be12fe5202b6d53, 0x9efa35377913cece,
        0xe2133a89d1fc3049, 0x252c3fdc2ae491c4, 0x6845452e83ccf33f, 0xab5e4a80dcb554ba,
    ];
    const COSET_GENERATORS_2: [u64; 8] = [
        0x1, 0x2, 0x2, 0x2, 0x2, 0x3, 0x3, 0x3,
    ];
    const COSET_GENERATORS_3: [u64; 8] = [
        0x6fffffff9, 0x7fffffff8, 0x8fffffff7, 0x9fffffff6,
        0xafffffff5, 0xbfffffff4, 0xcfffffff3, 0xdfffffff2,
    ];
    const MODULUS_IS_BIG: bool = true; // 0xFFFFFFFF00000000 >= 0x4000...
}

pub type Secp256r1Fr = Field<Secp256r1FrParams>;

// ---------------------------------------------------------------------------
// secp256r1 G1 Curve Parameters
// ---------------------------------------------------------------------------

/// secp256r1 (P-256) G1: y^2 = x^3 - 3x + b.
/// All coordinates stored in standard (non-Montgomery) form,
/// matching C++ which calls `.to_montgomery_form()` at runtime.
pub struct Secp256r1G1Params;

impl CurveParams for Secp256r1G1Params {
    type BaseFieldParams = Secp256r1FqParams;
    type ScalarFieldParams = Secp256r1FrParams;

    const HAS_A: bool = true;

    /// a = p - 3 (standard form).
    const A: [u64; 4] = [
        0xFFFFFFFFFFFFFFFC,
        0x00000000FFFFFFFF,
        0x0000000000000000,
        0xFFFFFFFF00000001,
    ];

    /// b (standard form).
    const B: [u64; 4] = [
        0x3BCE3C3E27D2604B,
        0x651D06B0CC53B0F6,
        0xB3EBBD55769886BC,
        0x5AC635D8AA3A93E7,
    ];

    /// Generator x (standard form).
    const GENERATOR_X: [u64; 4] = [
        0xF4A13945D898C296,
        0x77037D812DEB33A0,
        0xF8BCE6E563A440F2,
        0x6B17D1F2E12C4247,
    ];

    /// Generator y (standard form).
    const GENERATOR_Y: [u64; 4] = [
        0xCBB6406837BF51F5,
        0x2BCE33576B315ECE,
        0x8EE7EB4A7C0F9E16,
        0x4FE342E2FE1A7F9B,
    ];

    const USE_ENDOMORPHISM: bool = false;

    /// Convert from standard form at runtime, matching C++.
    fn generator_x() -> Secp256r1Fq {
        Field::from_limbs(Self::GENERATOR_X)
    }
    fn generator_y() -> Secp256r1Fq {
        Field::from_limbs(Self::GENERATOR_Y)
    }
    fn coeff_a() -> Secp256r1Fq {
        Field::from_limbs(Self::A)
    }
    fn coeff_b() -> Secp256r1Fq {
        Field::from_limbs(Self::B)
    }
}

pub type G1Affine = crate::ecc::groups::affine_element::AffineElement<Secp256r1G1Params>;
pub type G1Element = crate::ecc::groups::element::Element<Secp256r1G1Params>;
