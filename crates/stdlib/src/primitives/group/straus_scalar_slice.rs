//! Scalar decomposition for Straus MSM algorithm.
//!
//! Port of `barretenberg/stdlib/primitives/group/straus_scalar_slice.hpp` and `.cpp`.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::groups::curve_params::CurveParams;

use crate::primitives::field::FieldT;
use crate::primitives::witness::BuilderRef;

use super::cycle_scalar::{CycleScalarT, HI_BITS, LO_BITS};

/// Decomposed scalar slices for the Straus MSM algorithm.
pub struct StrausScalarSlices<C: CurveParams> {
    pub table_bits: usize,
    pub slices: Vec<FieldT<C::BaseFieldParams>>,
    pub slices_native: Vec<u64>,
}

impl<C: CurveParams> StrausScalarSlices<C> {
    /// Decompose a cycle scalar into bit-slices of size `table_bits`.
    pub fn new(
        ctx: BuilderRef<C::BaseFieldParams>,
        scalar: &CycleScalarT<C>,
        table_bits: usize,
    ) -> Self {
        let (lo_slices, lo_native) =
            Self::compute_scalar_slices(&ctx, scalar.lo(), LO_BITS, table_bits);
        let (hi_slices, hi_native) =
            Self::compute_scalar_slices(&ctx, scalar.hi(), HI_BITS, table_bits);

        let mut slices = lo_slices;
        slices.extend(hi_slices);
        let mut slices_native = lo_native;
        slices_native.extend(hi_native);

        Self {
            table_bits,
            slices,
            slices_native,
        }
    }

    /// Access a slice by round index.
    pub fn get(&self, index: usize) -> &FieldT<C::BaseFieldParams> {
        &self.slices[index]
    }

    fn compute_scalar_slices(
        ctx: &BuilderRef<C::BaseFieldParams>,
        scalar: &FieldT<C::BaseFieldParams>,
        num_bits: usize,
        table_bits: usize,
    ) -> (Vec<FieldT<C::BaseFieldParams>>, Vec<u64>) {
        let num_slices = (num_bits + table_bits - 1) / table_bits;
        let mut stdlib_slices = Vec::with_capacity(num_slices);
        let mut native_slices = Vec::with_capacity(num_slices);

        if num_bits == 0 {
            return (stdlib_slices, native_slices);
        }

        if scalar.is_constant() {
            let table_mask = (1u64 << table_bits) - 1u64;
            let raw_value = scalar.get_value().from_montgomery_form();
            let mut data = raw_value.data;
            for _ in 0..num_slices {
                let slice_value = data[0] & table_mask;
                stdlib_slices
                    .push(FieldT::from_field(Field::from(slice_value)));
                native_slices.push(slice_value);
                // Shift right by table_bits
                shift_right_u256(&mut data, table_bits);
            }
            return (stdlib_slices, native_slices);
        }

        // Non-constant: decompose in-circuit via the builder
        let normalized = scalar.normalize();
        let slice_indices = ctx.borrow_mut().decompose_into_default_range(
            normalized.witness_index,
            num_bits as u64,
            table_bits as u64,
            "straus_scalar_slice decompose_into_default_range",
        );

        for idx in &slice_indices {
            let slice = FieldT::from_witness_index(ctx.clone(), *idx);
            let value = slice.get_value().from_montgomery_form().data[0];
            stdlib_slices.push(slice);
            native_slices.push(value);
        }

        (stdlib_slices, native_slices)
    }
}

/// Shift a [u64; 4] right by `bits` in place.
fn shift_right_u256(data: &mut [u64; 4], bits: usize) {
    let limb_shift = bits / 64;
    let bit_shift = bits % 64;
    for i in 0..4 {
        let src = i + limb_shift;
        if src < 4 {
            data[i] = data[src] >> bit_shift;
            if bit_shift > 0 && src + 1 < 4 {
                data[i] |= data[src + 1] << (64 - bit_shift);
            }
        } else {
            data[i] = 0;
        }
    }
}
