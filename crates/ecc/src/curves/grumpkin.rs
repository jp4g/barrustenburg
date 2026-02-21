use crate::curves::bn254::{Bn254FqParams, Bn254FrParams};
use crate::fields::field::Field;
use crate::groups::curve_params::CurveParams;

/// Grumpkin is the "field-swapped twin" of BN254:
/// - Grumpkin Fq = BN254 Fr (same modulus)
/// - Grumpkin Fr = BN254 Fq
pub type GrumpkinFq = Field<Bn254FrParams>;
pub type GrumpkinFr = Field<Bn254FqParams>;

/// Grumpkin G1: y^2 = x^3 - 17, generator at (1, y).
/// Values stored in Montgomery form (matching C++ brace-init).
pub struct GrumpkinG1Params;

impl CurveParams for GrumpkinG1Params {
    type BaseFieldParams = Bn254FrParams; // Grumpkin Fq = BN254 Fr
    type ScalarFieldParams = Bn254FqParams; // Grumpkin Fr = BN254 Fq

    const HAS_A: bool = false;
    const A: [u64; 4] = [0, 0, 0, 0];

    /// b = -17 in Montgomery form.
    const B: [u64; 4] = [
        0xdd7056026000005a,
        0x223fa97acb319311,
        0xcc388229877910c0,
        0x034394632b724eaa,
    ];

    /// Generator x = 1, stored as zero placeholder; overridden by generator_x().
    const GENERATOR_X: [u64; 4] = [0, 0, 0, 0];

    /// Generator y in Montgomery form (from C++ grumpkin.hpp).
    const GENERATOR_Y: [u64; 4] = [
        0x11b2dff1448c41d8,
        0x23d3446f21c77dc3,
        0xaa7b8cf435dfafbb,
        0x14b34cf69dc25d68,
    ];

    const USE_ENDOMORPHISM: bool = true;

    /// Generator x = fr::one() = Montgomery(1), computed at runtime.
    fn generator_x() -> GrumpkinFq {
        Field::one()
    }
}

pub type G1Affine = crate::groups::affine_element::AffineElement<GrumpkinG1Params>;
pub type G1Element = crate::groups::element::Element<GrumpkinG1Params>;
