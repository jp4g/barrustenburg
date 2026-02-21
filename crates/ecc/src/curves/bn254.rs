use crate::fields::field::Field;
use crate::fields::field_params::FieldParams;
use crate::groups::curve_params::CurveParams;

// ---------------------------------------------------------------------------
// BN254 Base Field (Fq)
// ---------------------------------------------------------------------------

pub struct Bn254FqParams;

impl FieldParams for Bn254FqParams {
    const MODULUS: [u64; 4] = [
        0x3C208C16D87CFD47,
        0x97816a916871ca8d,
        0xb85045b68181585d,
        0x30644e72e131a029,
    ];
    const R_SQUARED: [u64; 4] = [
        0xF32CFC5B538AFA89,
        0xB5E71911D44501FB,
        0x47AB1EFF0A417FF6,
        0x06D89F71CAB8351F,
    ];
    const R_INV: u64 = 0x87d20782e4866389;
    const CUBE_ROOT: [u64; 4] = [
        0x71930c11d782e155,
        0xa6bb947cffbe3323,
        0xaa303344d4741444,
        0x2c3b3f0d26594943,
    ];
    const PRIMITIVE_ROOT: [u64; 4] = [0, 0, 0, 0]; // Not used for Fq
    const COSET_GENERATORS_0: [u64; 8] = [
        0x7a17caa950ad28d7, 0x4d750e37163c3674, 0x20d251c4dbcb4411, 0xf42f9552a15a51ae,
        0x4f4bc0b2b5ef64bd, 0x22a904407b7e725a, 0xf60647ce410d7ff7, 0xc9638b5c069c8d94,
    ];
    const COSET_GENERATORS_1: [u64; 8] = [
        0x1f6ac17ae15521b9, 0x29e3aca3d71c2cf7, 0x345c97cccce33835, 0x3ed582f5c2aa4372,
        0x1a4b98fbe78db996, 0x24c48424dd54c4d4, 0x2f3d6f4dd31bd011, 0x39b65a76c8e2db4f,
    ];
    const COSET_GENERATORS_2: [u64; 8] = [
        0x334bea4e696bd284, 0x99ba8dbde1e518b0, 0x29312d5a5e5edc,   0x6697d49cd2d7a508,
        0x5c65ec9f484e3a79, 0xc2d4900ec0c780a5, 0x2943337e3940c6d1, 0x8fb1d6edb1ba0cfd,
    ];
    const COSET_GENERATORS_3: [u64; 8] = [
        0x2a1f6744ce179d8e, 0x3829df06681f7cbd, 0x463456c802275bed, 0x543ece899c2f3b1c,
        0x180a96573d3d9f8,  0xf8b21270ddbb927,  0x1d9598e8a7e39857, 0x2ba010aa41eb7786,
    ];
    const MODULUS_IS_BIG: bool = false; // 0x3064... < 0x4000...

    // GLV endomorphism constants for BN254 Fq (Grumpkin's scalar field)
    // From C++ fq.hpp:82-90 — different from Fr's constants
    const ENDO_G1_LO: u64 = 0x7a7bd9d4391eb18d;
    const ENDO_G1_MID: u64 = 0x4ccef014a773d2cf;
    const ENDO_G1_HI: u64 = 0x0000000000000002;
    const ENDO_G2_LO: u64 = 0xd91d232ec7e0b3d2;
    const ENDO_G2_MID: u64 = 0x0000000000000002;
    const ENDO_MINUS_B1_LO: u64 = 0x8211bbeb7d4f1129;
    const ENDO_MINUS_B1_MID: u64 = 0x6f4d8248eeb859fc;
    const ENDO_B2_LO: u64 = 0x89d3256894d213e2;
}

pub type Fq = Field<Bn254FqParams>;

// ---------------------------------------------------------------------------
// BN254 Scalar Field (Fr)
// ---------------------------------------------------------------------------

pub struct Bn254FrParams;

impl FieldParams for Bn254FrParams {
    const MODULUS: [u64; 4] = [
        0x43E1F593F0000001,
        0x2833E84879B97091,
        0xB85045B68181585D,
        0x30644E72E131A029,
    ];
    const R_SQUARED: [u64; 4] = [
        0x1BB8E645AE216DA7,
        0x53FE3AB1E35C59E3,
        0x8C49833D53BB8085,
        0x0216D0B17F4E44A5,
    ];
    const R_INV: u64 = 0xc2e1f593efffffff;
    const CUBE_ROOT: [u64; 4] = [
        0x93e7cede4a0329b3,
        0x7d4fdca77a96c167,
        0x8be4ba08b19a750a,
        0x1cbd5653a5661c25,
    ];
    const PRIMITIVE_ROOT: [u64; 4] = [
        0x636e735580d13d9c,
        0xa22bf3742445ffd6,
        0x56452ac01eb203d8,
        0x1860ef942963f9e7,
    ];
    const COSET_GENERATORS_0: [u64; 8] = [
        0x5eef048d8fffffe7, 0xb8538a9dfffffe2,  0x3057819e4fffffdb, 0xdcedb5ba9fffffd6,
        0x8983e9d6efffffd1, 0x361a1df33fffffcc, 0xe2b0520f8fffffc7, 0x8f46862bdfffffc2,
    ];
    const COSET_GENERATORS_1: [u64; 8] = [
        0x12ee50ec1ce401d0, 0x49eac781bc44cefa, 0x307f6d866832bb01, 0x677be41c0793882a,
        0x9e785ab1a6f45454, 0xd574d1474655227e, 0xc7147dce5b5efa7,  0x436dbe728516bcd1,
    ];
    const COSET_GENERATORS_2: [u64; 8] = [
        0x29312d5a5e5ee7,   0x6697d49cd2d7a515, 0x5c65ec9f484e3a89, 0xc2d4900ec0c780b7,
        0x2943337e3940c6e5, 0x8fb1d6edb1ba0d13, 0xf6207a5d2a335342, 0x5c8f1dcca2ac9970,
    ];
    const COSET_GENERATORS_3: [u64; 8] = [
        0x463456c802275bed, 0x543ece899c2f3b1c, 0x180a96573d3d9f8,  0xf8b21270ddbb927,
        0x1d9598e8a7e39857, 0x2ba010aa41eb7786, 0x39aa886bdbf356b5, 0x47b5002d75fb35e5,
    ];
    const MODULUS_IS_BIG: bool = false; // 0x3064... < 0x4000...
    const HAS_HIGH_2ADICITY: bool = true; // BN254 Fr has 2-adicity 28

    // GLV endomorphism constants for BN254 Fr
    const ENDO_G1_LO: u64 = 0x7a7bd9d4391eb18d;
    const ENDO_G1_MID: u64 = 0x4ccef014a773d2cf;
    const ENDO_G1_HI: u64 = 0x0000000000000002;
    const ENDO_G2_LO: u64 = 0xd91d232ec7e0b3d7;
    const ENDO_G2_MID: u64 = 0x0000000000000002;
    const ENDO_MINUS_B1_LO: u64 = 0x8211bbeb7d4f1128;
    const ENDO_MINUS_B1_MID: u64 = 0x6f4d8248eeb859fc;
    const ENDO_B2_LO: u64 = 0x89d3256894d213e3;
}

pub type Fr = Field<Bn254FrParams>;

// ---------------------------------------------------------------------------
// BN254 G1 Curve Parameters
// ---------------------------------------------------------------------------

/// BN254 G1: y^2 = x^3 + 3, generator at (1, y).
/// Values stored in Montgomery form (matching C++ brace-init).
pub struct Bn254G1Params;

impl CurveParams for Bn254G1Params {
    type BaseFieldParams = Bn254FqParams;
    type ScalarFieldParams = Bn254FrParams;

    const HAS_A: bool = false;
    const A: [u64; 4] = [0, 0, 0, 0];

    /// b = 3 in Montgomery form.
    const B: [u64; 4] = [
        0x7a17caa950ad28d7,
        0x1f6ac17ae15521b9,
        0x334bea4e696bd284,
        0x2a1f6744ce179d8e,
    ];

    /// Generator x = 1, stored as zero placeholder; overridden by generator_x().
    const GENERATOR_X: [u64; 4] = [0, 0, 0, 0];

    /// Generator y in Montgomery form (from C++ g1.hpp).
    const GENERATOR_Y: [u64; 4] = [
        0xa6ba871b8b1e1b3a,
        0x14f1d651eb8e167b,
        0xccdd46def0f28c58,
        0x1c14ef83340fbe5e,
    ];

    const USE_ENDOMORPHISM: bool = true;

    /// Generator x = fq::one() = Montgomery(1), computed at runtime.
    fn generator_x() -> Fq {
        Field::one()
    }
}

pub type G1Affine = crate::groups::affine_element::AffineElement<Bn254G1Params>;
pub type G1Element = crate::groups::element::Element<Bn254G1Params>;
