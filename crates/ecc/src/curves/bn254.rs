use crate::fields::field::Field;
use crate::fields::field2::Field2;
use crate::fields::field6::Field6Params;
use crate::fields::field12::Field12Params;
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

// ---------------------------------------------------------------------------
// Extension field type aliases
// ---------------------------------------------------------------------------

pub type Fq2 = Field2<Bn254FqParams>;
pub type Fq6 = crate::fields::field6::Field6<Bn254FqParams>;
pub type Fq12 = crate::fields::field12::Field12<Bn254FqParams>;

// ---------------------------------------------------------------------------
// Fq2 BN254-specific constants (from C++ fq2.hpp, __SIZEOF_INT128__ path)
// ---------------------------------------------------------------------------

impl Field2<Bn254FqParams> {
    /// Twist coefficient b' for BN254 G2: b' = b / xi = 3 / (9 + u)
    pub fn twist_coeff_b() -> Self {
        Self::new(
            Fq::from_raw([0x3bf938e377b802a8, 0x020b1b273633535d, 0x26b7edf049755260, 0x2514c6324384a86d]),
            Fq::from_raw([0x38e7ecccd1dcff67, 0x65f0b37d93ce0d3e, 0xd749d0dd22ac00aa, 0x0141b9ce4a688d4d]),
        )
    }

    /// Frobenius twist constant for G2 x-coordinate (mul_by_q).
    pub fn twist_mul_by_q_x() -> Self {
        Self::new(
            Fq::from_raw([0xb5773b104563ab30, 0x347f91c8a9aa6454, 0x7a007127242e0991, 0x1956bcd8118214ec]),
            Fq::from_raw([0x6e849f1ea0aa4757, 0xaa1c7b6d89f89141, 0xb6e713cdfae0ca3a, 0x26694fbb4e82ebc3]),
        )
    }

    /// Frobenius twist constant for G2 y-coordinate (mul_by_q).
    pub fn twist_mul_by_q_y() -> Self {
        Self::new(
            Fq::from_raw([0xe4bbdd0c2936b629, 0xbb30f162e133bacb, 0x31a9d1b6f9645366, 0x253570bea500f8dd]),
            Fq::from_raw([0xa1d77ce45ffe77c7, 0x07affd117826d1db, 0x6d16bd27bb7edc6b, 0x2c87200285defecc]),
        )
    }
}

// ---------------------------------------------------------------------------
// Fq6 parameters for BN254 (from C++ fq6.hpp, __SIZEOF_INT128__ path)
// ---------------------------------------------------------------------------

impl Field6Params for Bn254FqParams {
    /// Multiply Fq2 element by non-residue xi = 9 + u.
    /// (a0 + a1*u)(9 + u) = (9*a0 - a1) + (9*a1 + a0)*u
    fn mul_by_non_residue(a: &Field2<Self>) -> Field2<Self> {
        let mut t0 = a.c0 + a.c0; // 2*a0
        t0 = t0 + t0;              // 4*a0
        t0 = t0 + t0;              // 8*a0
        t0 = t0 + a.c0;            // 9*a0
        let mut t1 = a.c1 + a.c1; // 2*a1
        t1 = t1 + t1;              // 4*a1
        t1 = t1 + t1;              // 8*a1
        t1 = t1 + a.c1;            // 9*a1
        Field2::new(t0 - a.c1, t1 + a.c0)
    }

    fn frobenius_coeffs_c1_1() -> Field2<Self> {
        Field2::new(
            Fq::from_raw([0xb5773b104563ab30, 0x347f91c8a9aa6454, 0x7a007127242e0991, 0x1956bcd8118214ec]),
            Fq::from_raw([0x6e849f1ea0aa4757, 0xaa1c7b6d89f89141, 0xb6e713cdfae0ca3a, 0x26694fbb4e82ebc3]),
        )
    }

    fn frobenius_coeffs_c1_2() -> Field2<Self> {
        Field2::new(
            Fq::from_raw([0x3350c88e13e80b9c, 0x7dce557cdb5e56b9, 0x6001b4b8b615564a, 0x2682e617020217e0]),
            Fq::zero(),
        )
    }

    fn frobenius_coeffs_c1_3() -> Field2<Self> {
        Field2::new(
            Fq::from_raw([0xc9af22f716ad6bad, 0xb311782a4aa662b2, 0x19eeaf64e248c7f4, 0x20273e77e3439f82]),
            Fq::from_raw([0xacc02860f7ce93ac, 0x3933d5817ba76b4c, 0x69e6188b446c8467, 0x0a46036d4417cc55]),
        )
    }

    fn frobenius_coeffs_c2_1() -> Field2<Self> {
        Field2::new(
            Fq::from_raw([0x7361d77f843abe92, 0xa5bb2bd3273411fb, 0x9c941f314b3e2399, 0x15df9cddbb9fd3ec]),
            Fq::from_raw([0x5dddfd154bd8c949, 0x62cb29a5a4445b60, 0x37bc870a0c7dd2b9, 0x24830a9d3171f0fd]),
        )
    }

    fn frobenius_coeffs_c2_2() -> Field2<Self> {
        Field2::new(
            Fq::from_raw([0x71930c11d782e155, 0xa6bb947cffbe3323, 0xaa303344d4741444, 0x2c3b3f0d26594943]),
            Fq::zero(),
        )
    }

    fn frobenius_coeffs_c2_3() -> Field2<Self> {
        Field2::new(
            Fq::from_raw([0x448a93a57b6762df, 0xbfd62df528fdeadf, 0xd858f5d00e9bd47a, 0x06b03d4d3476ec58]),
            Fq::from_raw([0x2b19daf4bcc936d1, 0xa1a54e7a56f4299f, 0xb533eee05adeaef1, 0x170c812b84dda0b2]),
        )
    }
}

// ---------------------------------------------------------------------------
// Fq12 parameters for BN254 (from C++ fq12.hpp, __SIZEOF_INT128__ path)
// ---------------------------------------------------------------------------

impl Field12Params for Bn254FqParams {
    fn frobenius_coefficients_1() -> Field2<Self> {
        Field2::new(
            Fq::from_raw([0xaf9ba69633144907, 0xca6b1d7387afb78a, 0x11bded5ef08a2087, 0x02f34d751a1f3a7c]),
            Fq::from_raw([0xa222ae234c492d72, 0xd00f02a4565de15b, 0xdc2ff3a253dfc926, 0x10a75716b3899551]),
        )
    }

    fn frobenius_coefficients_2() -> Field2<Self> {
        Field2::new(
            Fq::from_raw([0xca8d800500fa1bf2, 0xf0c5d61468b39769, 0x0e201271ad0d4418, 0x04290f65bad856e6]),
            Fq::zero(),
        )
    }

    fn frobenius_coefficients_3() -> Field2<Self> {
        Field2::new(
            Fq::from_raw([0x365316184e46d97d, 0x0af7129ed4c96d9f, 0x659da72fca1009b5, 0x08116d8983a20d23]),
            Fq::from_raw([0xb1df4af7c39c1939, 0x3d9f02878a73bf7f, 0x9b2220928caf0ae0, 0x26684515eff054a6]),
        )
    }
}

// ---------------------------------------------------------------------------
// BN254 G2 types (over Fq2)
// ---------------------------------------------------------------------------

/// BN254 G2 point in Jacobian projective coordinates over Fq2.
pub struct G2Element {
    pub x: Fq2,
    pub y: Fq2,
    pub z: Fq2,
}

impl Clone for G2Element {
    #[inline]
    fn clone(&self) -> Self {
        Self { x: self.x, y: self.y, z: self.z }
    }
}

impl Copy for G2Element {}

impl std::fmt::Debug for G2Element {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "G2Element({:?}, {:?}, {:?})", self.x, self.y, self.z)
    }
}

/// BN254 G2 point in affine coordinates over Fq2.
pub struct G2AffineElement {
    pub x: Fq2,
    pub y: Fq2,
}

impl Clone for G2AffineElement {
    #[inline]
    fn clone(&self) -> Self {
        Self { x: self.x, y: self.y }
    }
}

impl Copy for G2AffineElement {}

impl std::fmt::Debug for G2AffineElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "G2AffineElement({:?}, {:?})", self.x, self.y)
    }
}

impl G2Element {
    #[inline]
    pub fn new(x: Fq2, y: Fq2, z: Fq2) -> Self {
        Self { x, y, z }
    }

    pub fn from_affine(aff: &G2AffineElement) -> Self {
        if aff.is_point_at_infinity() {
            return Self::infinity();
        }
        Self::new(aff.x, aff.y, Fq2::one())
    }

    pub fn infinity() -> Self {
        let mut r = Self::new(Fq2::zero(), Fq2::zero(), Fq2::zero());
        r.x.c0.self_set_msb();
        r
    }

    #[inline]
    pub fn is_point_at_infinity(&self) -> bool {
        self.x.c0.is_msb_set()
    }
}

impl std::ops::Neg for G2Element {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(self.x, -self.y, self.z)
    }
}

impl G2AffineElement {
    #[inline]
    pub fn new(x: Fq2, y: Fq2) -> Self {
        Self { x, y }
    }

    pub fn infinity() -> Self {
        let mut r = Self::new(Fq2::zero(), Fq2::zero());
        r.x.c0.self_set_msb();
        r
    }

    #[inline]
    pub fn is_point_at_infinity(&self) -> bool {
        self.x.c0.is_msb_set()
    }

    /// BN254 G2 generator (from C++ g2.hpp, __SIZEOF_INT128__ path).
    pub fn generator() -> Self {
        Self::new(
            Fq2::new(
                Fq::from_raw([0x8e83b5d102bc2026, 0xdceb1935497b0172, 0xfbb8264797811adf, 0x19573841af96503b]),
                Fq::from_raw([0xafb4737da84c6140, 0x6043dd5a5802d8c4, 0x09e950fc52a02f86, 0x14fef0833aea7b6b]),
            ),
            Fq2::new(
                Fq::from_raw([0x619dfa9d886be9f6, 0xfe7fd297f59e9b78, 0xff9e1a62231b7dfe, 0x28fd7eebae9e4206]),
                Fq::from_raw([0x64095b56c71856ee, 0xdc57f922327d3cbb, 0x55f935be33351076, 0x0da4a0e693fd6482]),
            ),
        )
    }
}
