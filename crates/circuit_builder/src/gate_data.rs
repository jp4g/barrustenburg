//! Gate data structures for circuit constraint specification.
//!
//! Port of `barretenberg/honk/execution_trace/gate_data.hpp`.
//! These structs describe the wire indices and scaling factors for various gate types
//! used by the Ultra circuit builder.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

/// 3-wire addition gate: a*a_scaling + b*b_scaling + c*c_scaling + const_scaling = 0
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AddTriple<P: FieldParams> {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub a_scaling: Field<P>,
    pub b_scaling: Field<P>,
    pub c_scaling: Field<P>,
    pub const_scaling: Field<P>,
}

/// 4-wire addition gate: a*a_scaling + b*b_scaling + c*c_scaling + d*d_scaling + const_scaling = 0
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AddQuad<P: FieldParams> {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub d: u32,
    pub a_scaling: Field<P>,
    pub b_scaling: Field<P>,
    pub c_scaling: Field<P>,
    pub d_scaling: Field<P>,
    pub const_scaling: Field<P>,
}

/// 4-wire mul-add gate: a*b*mul_scaling + a*a_scaling + b*b_scaling + c*c_scaling + d*d_scaling + const_scaling = 0
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MulQuad<P: FieldParams> {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub d: u32,
    pub mul_scaling: Field<P>,
    pub a_scaling: Field<P>,
    pub b_scaling: Field<P>,
    pub c_scaling: Field<P>,
    pub d_scaling: Field<P>,
    pub const_scaling: Field<P>,
}

/// Arithmetic gate with standard selector naming: q_m*a*b + q_l*a + q_r*b + q_o*c + q_c = 0
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArithmeticTriple<P: FieldParams> {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub q_m: Field<P>,
    pub q_l: Field<P>,
    pub q_r: Field<P>,
    pub q_o: Field<P>,
    pub q_c: Field<P>,
}

/// Goblin ECCVM operation tuple: stores op type, point coordinates (split into limbs), and scalar.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EccOpTuple {
    pub op: u32,
    pub x_lo: u32,
    pub x_hi: u32,
    pub y_lo: u32,
    pub y_hi: u32,
    pub z_1: u32,
    pub z_2: u32,
    pub return_is_infinity: bool,
}

/// Embedded curve point addition/subtraction: (x1, y1) +/- (x2, y2) = (x3, y3)
///
/// `sign_coefficient` is +1 for addition, -1 for subtraction (stored as field element
/// in q_1 / q_sign of the elliptic block).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EccAddGate<P: FieldParams> {
    pub x1: u32,
    pub y1: u32,
    pub x2: u32,
    pub y2: u32,
    pub x3: u32,
    pub y3: u32,
    pub sign_coefficient: Field<P>,
}

/// Embedded curve point doubling: 2 * (x1, y1) = (x3, y3)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EccDblGate {
    pub x1: u32,
    pub y1: u32,
    pub x3: u32,
    pub y3: u32,
}

/// Databus lookup gate: reads value at index from calldata/returndata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatabusLookupGate {
    pub index: u32,
    pub value: u32,
}

/// Poseidon2 external round gate data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poseidon2ExternalGate {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub d: u32,
    pub round_idx: usize,
}

/// Poseidon2 internal round gate data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poseidon2InternalGate {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub d: u32,
    pub round_idx: usize,
}
