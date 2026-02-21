//! Port of `InputElements` from `ultra_relation_consistency.test.cpp`.
//!
//! Provides a struct with named field references into a flat array, matching
//! the C++ `AllEntities` pattern used by relation `accumulate()` methods.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

/// Number of wire/selector/entity values. Matches C++ `InputElements::NUM_ELEMENTS = 45`.
/// We use 42 to match the actual named fields (the C++ test uses up to index 41).
pub const NUM_ELEMENTS: usize = 42;

/// Input elements for Ultra relation accumulation.
///
/// Each field is a named reference into a flat array, matching the C++ `AllEntities` pattern.
/// For verifier-side testing, these are plain field values.
pub struct InputElements<P: FieldParams> {
    pub data: [Field<P>; NUM_ELEMENTS],
}

impl<P: FieldParams> InputElements<P> {
    /// Generate random input elements (for testing).
    pub fn get_random() -> Self {
        let mut data = [Field::zero(); NUM_ELEMENTS];
        for e in data.iter_mut() {
            *e = Field::random_element();
        }
        Self { data }
    }

    /// Generate deterministic "special" input elements: 1, 2, 3, ... (for testing).
    pub fn get_special() -> Self {
        let mut data = [Field::zero(); NUM_ELEMENTS];
        for (i, e) in data.iter_mut().enumerate() {
            *e = Field::from((i + 1) as u64);
        }
        Self { data }
    }

    // Selectors
    pub fn q_c(&self) -> Field<P> { self.data[0] }
    pub fn q_l(&self) -> Field<P> { self.data[1] }
    pub fn q_r(&self) -> Field<P> { self.data[2] }
    pub fn q_o(&self) -> Field<P> { self.data[3] }
    pub fn q_4(&self) -> Field<P> { self.data[4] }
    pub fn q_m(&self) -> Field<P> { self.data[5] }
    pub fn q_arith(&self) -> Field<P> { self.data[6] }
    pub fn q_delta_range(&self) -> Field<P> { self.data[7] }
    pub fn q_elliptic(&self) -> Field<P> { self.data[8] }
    pub fn q_memory(&self) -> Field<P> { self.data[9] }
    pub fn q_nnf(&self) -> Field<P> { self.data[10] }
    pub fn q_lookup(&self) -> Field<P> { self.data[11] }
    pub fn q_poseidon2_external(&self) -> Field<P> { self.data[12] }
    pub fn q_poseidon2_internal(&self) -> Field<P> { self.data[13] }

    // Sigmas
    pub fn sigma_1(&self) -> Field<P> { self.data[14] }
    pub fn sigma_2(&self) -> Field<P> { self.data[15] }
    pub fn sigma_3(&self) -> Field<P> { self.data[16] }
    pub fn sigma_4(&self) -> Field<P> { self.data[17] }

    // IDs
    pub fn id_1(&self) -> Field<P> { self.data[18] }
    pub fn id_2(&self) -> Field<P> { self.data[19] }
    pub fn id_3(&self) -> Field<P> { self.data[20] }
    pub fn id_4(&self) -> Field<P> { self.data[21] }

    // Tables
    pub fn table_1(&self) -> Field<P> { self.data[22] }
    pub fn table_2(&self) -> Field<P> { self.data[23] }
    pub fn table_3(&self) -> Field<P> { self.data[24] }
    pub fn table_4(&self) -> Field<P> { self.data[25] }

    // Lagrange
    pub fn lagrange_first(&self) -> Field<P> { self.data[26] }
    pub fn lagrange_last(&self) -> Field<P> { self.data[27] }

    // Wires
    pub fn w_l(&self) -> Field<P> { self.data[28] }
    pub fn w_r(&self) -> Field<P> { self.data[29] }
    pub fn w_o(&self) -> Field<P> { self.data[30] }
    pub fn w_4(&self) -> Field<P> { self.data[31] }

    // Sorted accumulator
    pub fn sorted_accum(&self) -> Field<P> { self.data[32] }

    // Grand products
    pub fn z_perm(&self) -> Field<P> { self.data[33] }
    pub fn z_lookup(&self) -> Field<P> { self.data[34] }

    // Shifted wires
    pub fn w_l_shift(&self) -> Field<P> { self.data[35] }
    pub fn w_r_shift(&self) -> Field<P> { self.data[36] }
    pub fn w_o_shift(&self) -> Field<P> { self.data[37] }
    pub fn w_4_shift(&self) -> Field<P> { self.data[38] }

    // Shifted sorted accumulator
    pub fn sorted_accum_shift(&self) -> Field<P> { self.data[39] }

    // Shifted grand products
    pub fn z_perm_shift(&self) -> Field<P> { self.data[40] }
    pub fn z_lookup_shift(&self) -> Field<P> { self.data[41] }

    // Mutable setters for specific overrides in tests
    pub fn set_q_arith(&mut self, val: Field<P>) { self.data[6] = val; }
    pub fn set_q_l(&mut self, val: Field<P>) { self.data[1] = val; }
}
