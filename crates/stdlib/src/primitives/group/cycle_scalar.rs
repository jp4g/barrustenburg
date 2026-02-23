//! Circuit scalar type for the embedded curve's scalar field.
//!
//! Port of `barretenberg/stdlib/primitives/group/cycle_scalar.hpp` and `.cpp`.
//!
//! A `CycleScalarT<C>` represents a Grumpkin scalar (BN254 base field element)
//! as two limbs: `lo` (128 bits) and `hi` (126 bits).

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_ecc::groups::curve_params::CurveParams;

use crate::primitives::field::{validate_split_in_field_unsafe, FieldT};
use crate::primitives::witness::BuilderRef;

/// Number of bits in the scalar field (254 for both BN254 and Grumpkin).
pub const NUM_BITS: usize = 254;
/// Number of bits in the low limb.
pub const LO_BITS: usize = 128;
/// Number of bits in the high limb.
pub const HI_BITS: usize = 126;

/// Circuit representation of an embedded curve scalar field element.
///
/// The scalar is stored as `lo + hi * 2^LO_BITS` where `lo` is 128 bits
/// and `hi` is 126 bits. Range constraints on the limbs are deferred to
/// the `batch_mul` algorithm.
pub struct CycleScalarT<C: CurveParams> {
    lo: FieldT<C::BaseFieldParams>,
    hi: FieldT<C::BaseFieldParams>,
}

impl<C: CurveParams> Clone for CycleScalarT<C> {
    fn clone(&self) -> Self {
        Self {
            lo: self.lo.clone(),
            hi: self.hi.clone(),
        }
    }
}

impl<C: CurveParams> CycleScalarT<C> {
    /// Construct a constant cycle scalar from a native scalar field value.
    pub fn from_native(value: Field<C::ScalarFieldParams>) -> Self {
        let value_mont = value.from_montgomery_form();
        let (lo_limbs, hi_limbs) = decompose_lo_hi(&value_mont.data);
        Self {
            lo: FieldT::from_field(Field::from_limbs(lo_limbs)),
            hi: FieldT::from_field(Field::from_limbs(hi_limbs)),
        }
    }

    /// Construct from lo and hi field elements with validation.
    pub fn from_lo_hi(
        lo: FieldT<C::BaseFieldParams>,
        hi: FieldT<C::BaseFieldParams>,
    ) -> Self {
        let result = Self { lo, hi };
        result.validate_scalar_is_in_field();
        result
    }

    /// Construct from lo and hi field elements without validation (internal use).
    pub(crate) fn from_lo_hi_unchecked(
        lo: FieldT<C::BaseFieldParams>,
        hi: FieldT<C::BaseFieldParams>,
    ) -> Self {
        Self { lo, hi }
    }

    /// Create a cycle scalar witness from a native scalar field value.
    pub fn from_witness(
        ctx: BuilderRef<C::BaseFieldParams>,
        value: Field<C::ScalarFieldParams>,
    ) -> Self {
        let value_mont = value.from_montgomery_form();
        let (lo_limbs, hi_limbs) = decompose_lo_hi(&value_mont.data);
        let lo_field = Field::<C::BaseFieldParams>::from_limbs(lo_limbs);
        let hi_field = Field::<C::BaseFieldParams>::from_limbs(hi_limbs);

        let lo = FieldT::from_witness(ctx.clone(), lo_field);
        let hi = FieldT::from_witness(ctx, hi_field);

        let result = Self { lo, hi };
        result.validate_scalar_is_in_field();
        result
    }

    pub fn is_constant(&self) -> bool {
        self.lo.is_constant() && self.hi.is_constant()
    }

    /// Get the native scalar field value.
    pub fn get_value(&self) -> Field<C::ScalarFieldParams> {
        let lo_val = self.lo.get_value().from_montgomery_form();
        let hi_val = self.hi.get_value().from_montgomery_form();
        // Reconstruct: lo + hi * 2^LO_BITS
        let mut result = [0u64; 4];
        result[0] = lo_val.data[0];
        result[1] = lo_val.data[1];
        // hi is at most 126 bits, shifted left by 128
        result[2] = hi_val.data[0];
        result[3] = hi_val.data[1];
        Field::<C::ScalarFieldParams>::from_limbs(result)
    }

    pub fn get_context(&self) -> &Option<BuilderRef<C::BaseFieldParams>> {
        let ctx = self.lo.get_context();
        if ctx.is_some() {
            return ctx;
        }
        self.hi.get_context()
    }

    pub fn lo(&self) -> &FieldT<C::BaseFieldParams> {
        &self.lo
    }

    pub fn hi(&self) -> &FieldT<C::BaseFieldParams> {
        &self.hi
    }

    /// Validate that the scalar value is less than the scalar field modulus.
    fn validate_scalar_is_in_field(&self) {
        validate_split_in_field_unsafe(
            &self.lo,
            &self.hi,
            LO_BITS,
            &C::ScalarFieldParams::MODULUS,
        );
    }
}

/// Decompose a 256-bit value into lo (128 bits) and hi (126 bits) parts.
fn decompose_lo_hi(data: &[u64; 4]) -> ([u64; 4], [u64; 4]) {
    let lo = [data[0], data[1], 0, 0];
    let hi = [data[2], data[3], 0, 0];
    (lo, hi)
}
