//! Straus MSM lookup table.
//!
//! Port of `barretenberg/stdlib/primitives/group/straus_lookup_table.hpp` and `.cpp`.
//!
//! Computes a ROM-based lookup table of size `1 << table_bits` containing
//! precomputed points: `[G] + i*[P]` for i = 0..N-1, where `[G]` is an
//! offset generator to avoid point-at-infinity edge cases.

use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;
use bbrs_ecc::groups::element::Element;

use crate::primitives::field::FieldT;
use crate::primitives::witness::BuilderRef;

use super::cycle_group::CycleGroupT;

/// ROM-backed lookup table for the Straus MSM algorithm.
pub struct StrausLookupTable<C: CurveParams> {
    context: BuilderRef<C::BaseFieldParams>,
    rom_id: usize,
}

impl<C: CurveParams> StrausLookupTable<C> {
    /// Compute the native (witness) table entries as projective Elements.
    pub fn compute_native_table(
        base_point: &Element<C>,
        offset_generator: &Element<C>,
        table_bits: usize,
    ) -> Vec<Element<C>> {
        let table_size = 1usize << table_bits;
        let mut hints = Vec::with_capacity(table_size);
        hints.push(offset_generator.clone());
        for i in 1..table_size {
            let next = hints[i - 1].clone() + base_point.clone();
            hints.push(next);
        }
        hints
    }

    /// Construct a new Straus lookup table.
    pub fn new(
        context: BuilderRef<C::BaseFieldParams>,
        base_point: &CycleGroupT<C>,
        offset_generator: &CycleGroupT<C>,
        table_bits: usize,
        hints: Option<&[AffineElement<C>]>,
    ) -> Self {
        let table_size = 1usize << table_bits;
        let mut point_table: Vec<CycleGroupT<C>> = Vec::with_capacity(table_size);

        // Handle point at infinity: substitute with generator "one" for computation,
        // then conditionally correct entries afterward.
        let fallback = CycleGroupT::<C>::from_affine(AffineElement::<C>::one());
        let modded_x = FieldT::conditional_assign(
            &base_point.is_point_at_infinity(),
            &FieldT::from_field(fallback.x_val()),
            base_point.x(),
        );
        let modded_y = FieldT::conditional_assign(
            &base_point.is_point_at_infinity(),
            &FieldT::from_field(fallback.y_val()),
            base_point.y(),
        );
        let modded_base_point = CycleGroupT::<C>::from_xy_bool(
            modded_x,
            modded_y,
            super::super::bool::BoolT::from_constant(false),
            false,
        );

        let hints_available = hints.is_some()
            && !base_point.is_point_at_infinity().get_value();
        let get_hint = |i: usize| -> Option<AffineElement<C>> {
            if hints_available {
                Some(hints.unwrap()[i].clone())
            } else {
                None
            }
        };

        if base_point.is_constant()
            && !base_point.is_point_at_infinity().get_value()
        {
            // Constant base point: fix as witness for efficiency
            let modded_bp = CycleGroupT::<C>::from_constant_witness(
                context.clone(),
                &modded_base_point.get_value(),
            );
            let offset_pt = CycleGroupT::<C>::from_constant_witness(
                context.clone(),
                &offset_generator.get_value(),
            );
            point_table.push(offset_pt);
            for i in 1..table_size {
                point_table.push(
                    point_table[i - 1].unconditional_add(&modded_bp, get_hint(i - 1)),
                );
            }
        } else {
            // Non-constant: build table with x-coordinate collision checks
            let mut coordinate_check_product =
                FieldT::<C::BaseFieldParams>::from_u64(1);
            point_table.push(offset_generator.clone());
            for i in 1..table_size {
                let x_diff =
                    point_table[i - 1].x().clone() - modded_base_point.x().clone();
                coordinate_check_product =
                    coordinate_check_product * x_diff;
                point_table.push(
                    point_table[i - 1]
                        .unconditional_add(&modded_base_point, get_hint(i - 1)),
                );
            }
            coordinate_check_product
                .assert_is_not_zero("straus_lookup_table x-coordinate collision");

            // Correct for point-at-infinity input
            for i in 1..table_size {
                point_table[i] = CycleGroupT::<C>::conditional_assign(
                    &base_point.is_point_at_infinity(),
                    offset_generator,
                    &point_table[i],
                );
            }
        }

        // Construct ROM array
        let rom_id = context.borrow_mut().create_rom_array(table_size);
        for i in 0..table_size {
            // Ensure all points are witnesses (not constants)
            let pt = if point_table[i].is_constant() {
                CycleGroupT::<C>::from_constant_witness(
                    context.clone(),
                    &point_table[i].get_value(),
                )
            } else {
                point_table[i].clone()
            };
            let x_idx = pt.x().normalize().get_witness_index();
            let y_idx = pt.y().normalize().get_witness_index();
            context
                .borrow_mut()
                .set_rom_element_pair(rom_id, i, [x_idx, y_idx]);
        }

        Self { context, rom_id }
    }

    /// Read a point from the table at the given index.
    pub fn read(&self, index: &FieldT<C::BaseFieldParams>) -> CycleGroupT<C> {
        let mut idx = index.clone();
        // ROM index must be a witness
        if idx.is_constant() {
            idx = FieldT::from_witness(self.context.clone(), idx.get_value());
            idx.assert_equal(
                &FieldT::from_field(index.get_value()),
                "straus_lookup_table::read constant index",
            );
        }
        let [x_idx, y_idx] = self
            .context
            .borrow_mut()
            .read_rom_array_pair(self.rom_id, idx.get_witness_index());
        let x = FieldT::from_witness_index(self.context.clone(), x_idx);
        let y = FieldT::from_witness_index(self.context.clone(), y_idx);

        CycleGroupT::<C>::from_xy_bool(
            x,
            y,
            super::super::bool::BoolT::from_constant(false),
            false,
        )
    }
}
