//! UltraCircuitBuilder — gate methods for Ultra circuits.
//!
//! Port of `barretenberg/stdlib_circuit_builders/ultra_circuit_builder.hpp` and
//! `ultra_circuit_builder.cpp`.
//!
//! Provides the constructor, constant-variable caching, arithmetic/elliptic/range/sort
//! gate creation methods, witness fixing, finalization, tag management, and selector
//! consistency checks.

use std::collections::{BTreeMap, HashMap};

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use crate::builder_base::{CircuitBuilderBase, DEFAULT_TAG};
use crate::execution_trace::UltraExecutionTraceBlocks;
use crate::gate_data::{AddQuad, AddTriple, ArithmeticTriple, EccAddGate, EccDblGate, MulQuad};

/// The plookup range proof requires work linear in range size, thus cannot be used
/// directly for large ranges. Elements are decomposed into smaller chunks.
pub const DEFAULT_PLOOKUP_RANGE_BITNUM: u64 = 14;
/// Step size for the sorted-set delta-range checks.
pub const DEFAULT_PLOOKUP_RANGE_STEP_SIZE: u64 = 3;
/// Maximum value representable in the default plookup range.
pub const DEFAULT_PLOOKUP_RANGE_SIZE: u64 = (1 << DEFAULT_PLOOKUP_RANGE_BITNUM) - 1;

/// Number of wires per gate row.
const NUM_WIRES: usize = 4;

/// A range list: stores the set of variable indices that have been constrained
/// to lie in `[0, target_range]`, along with the tag/tau pair used for the
/// generalized permutation argument.
///
/// Port of C++ `UltraCircuitBuilder_::RangeList`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RangeList {
    pub target_range: u64,
    pub range_tag: u32,
    pub tau_tag: u32,
    pub variable_indices: Vec<u32>,
}

/// Ultra circuit builder: extends `CircuitBuilderBase` with execution trace blocks,
/// constant-variable caching, and gate methods for arithmetic, elliptic-curve,
/// range, and sort constraints.
///
/// Port of C++ `UltraCircuitBuilder_<UltraExecutionTraceBlocks>`.
#[derive(Debug, Clone)]
pub struct UltraCircuitBuilder<P: FieldParams> {
    /// Base builder with variable storage, copy constraints, tags, and error state.
    pub base: CircuitBuilderBase<P>,
    /// Storage for wires and selectors for all gate types.
    pub blocks: UltraExecutionTraceBlocks<P>,
    /// Cache: maps constant field values to their variable index, so the same
    /// constant is only added once.
    pub constant_variable_indices: HashMap<[u64; 4], u32>,
    /// Whether `finalize_circuit` has been called.
    pub circuit_finalized: bool,
    /// Range lists keyed by target range value. Each list collects all variables
    /// constrained to `[0, target_range]`. Processed during finalization into
    /// sorted-set delta-range gates.
    pub range_lists: BTreeMap<u64, RangeList>,
}

impl<P: FieldParams> UltraCircuitBuilder<P> {
    // ════════════════════════════════════════════════════════════════════
    //  Construction
    // ════════════════════════════════════════════════════════════════════

    /// Create a new UltraCircuitBuilder.
    ///
    /// Mirrors the C++ default constructor: adds a zero constant variable and
    /// initializes the tau permutation with `DEFAULT_TAG -> DEFAULT_TAG`.
    pub fn new() -> Self {
        let mut builder = Self {
            base: CircuitBuilderBase::new(),
            blocks: UltraExecutionTraceBlocks::new(),
            constant_variable_indices: HashMap::new(),
            circuit_finalized: false,
            range_lists: BTreeMap::new(),
        };
        let zero_idx = builder.put_constant_variable(Field::zero());
        builder.base.set_zero_idx(zero_idx);
        builder.base.tau_mut().insert(DEFAULT_TAG, DEFAULT_TAG);
        builder
    }

    /// Create a new UltraCircuitBuilder in write-vk mode.
    pub fn new_with_vk_mode(is_write_vk_mode: bool) -> Self {
        let mut builder = Self {
            base: CircuitBuilderBase::new_with_vk_mode(is_write_vk_mode),
            blocks: UltraExecutionTraceBlocks::new(),
            constant_variable_indices: HashMap::new(),
            circuit_finalized: false,
            range_lists: BTreeMap::new(),
        };
        let zero_idx = builder.put_constant_variable(Field::zero());
        builder.base.set_zero_idx(zero_idx);
        builder.base.tau_mut().insert(DEFAULT_TAG, DEFAULT_TAG);
        builder
    }

    /// Create from pre-existing witness values and public inputs (ACIR constructor).
    pub fn new_from_witnesses(
        witness_values: &[Field<P>],
        public_inputs: Vec<u32>,
        is_write_vk_mode: bool,
    ) -> Self {
        let mut builder = Self {
            base: CircuitBuilderBase::new_with_vk_mode(is_write_vk_mode),
            blocks: UltraExecutionTraceBlocks::new(),
            constant_variable_indices: HashMap::new(),
            circuit_finalized: false,
            range_lists: BTreeMap::new(),
        };
        for &value in witness_values {
            builder.base.add_variable(value);
        }
        builder.base.initialize_public_inputs(public_inputs);
        let zero_idx = builder.put_constant_variable(Field::zero());
        builder.base.set_zero_idx(zero_idx);
        builder.base.tau_mut().insert(DEFAULT_TAG, DEFAULT_TAG);
        builder
    }

    // ════════════════════════════════════════════════════════════════════
    //  Constant variable caching
    // ════════════════════════════════════════════════════════════════════

    /// Get or create a variable constrained to `variable`.
    ///
    /// If a variable with this value already exists in the cache, returns its
    /// index. Otherwise, creates a new variable, constrains it via `fix_witness`,
    /// and caches the mapping.
    pub fn put_constant_variable(&mut self, variable: Field<P>) -> u32 {
        let key = variable.data;
        if let Some(&idx) = self.constant_variable_indices.get(&key) {
            return idx;
        }
        let variable_index = self.base.add_variable(variable);
        self.fix_witness(variable_index, variable);
        self.constant_variable_indices.insert(key, variable_index);
        variable_index
    }

    // ════════════════════════════════════════════════════════════════════
    //  Gate creation: arithmetic
    // ════════════════════════════════════════════════════════════════════

    /// Create an addition gate: `a*a_scaling + b*b_scaling + c*c_scaling + const_scaling = 0`.
    ///
    /// Delegates to `create_big_add_gate` with the 4th wire set to zero.
    pub fn create_add_gate(&mut self, gate: &AddTriple<P>) {
        self.create_big_add_gate(
            &AddQuad {
                a: gate.a,
                b: gate.b,
                c: gate.c,
                d: self.base.zero_idx(),
                a_scaling: gate.a_scaling,
                b_scaling: gate.b_scaling,
                c_scaling: gate.c_scaling,
                d_scaling: Field::zero(),
                const_scaling: gate.const_scaling,
            },
            false,
        );
    }

    /// Create a 4-wire addition gate:
    /// `a*a_scaling + b*b_scaling + c*c_scaling + d*d_scaling + const_scaling = 0`.
    ///
    /// If `use_next_gate_w_4` is true, also adds `w_4` from the next gate row
    /// (q_arith is set to 2 instead of 1).
    pub fn create_big_add_gate(&mut self, gate: &AddQuad<P>, use_next_gate_w_4: bool) {
        self.base
            .assert_valid_variables(&[gate.a, gate.b, gate.c, gate.d]);
        self.blocks
            .arithmetic
            .block
            .populate_wires(gate.a, gate.b, gate.c, gate.d);
        self.blocks.arithmetic.block.q_m_mut().push(Field::zero());
        self.blocks
            .arithmetic
            .block
            .q_1_mut()
            .push(gate.a_scaling);
        self.blocks
            .arithmetic
            .block
            .q_2_mut()
            .push(gate.b_scaling);
        self.blocks
            .arithmetic
            .block
            .q_3_mut()
            .push(gate.c_scaling);
        self.blocks
            .arithmetic
            .block
            .q_c_mut()
            .push(gate.const_scaling);
        self.blocks
            .arithmetic
            .block
            .q_4_mut()
            .push(gate.d_scaling);
        let arith_val = if use_next_gate_w_4 {
            Field::from(2u64)
        } else {
            Field::from(1u64)
        };
        self.blocks.arithmetic.set_gate_selector(arith_val);
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
    }

    /// Create a 4-wire mul-add gate:
    /// `a*b*mul_scaling + a*a_scaling + b*b_scaling + c*c_scaling + d*d_scaling + const_scaling = 0`.
    ///
    /// If `use_next_gate_w_4` is true, also adds `w_4` from the next gate row.
    /// In that case, q_arith is set to 2 and `mul_scaling` is doubled to compensate
    /// for the factor-of-2 scaling applied by the ArithmeticRelation.
    pub fn create_big_mul_add_gate(&mut self, gate: &MulQuad<P>, use_next_gate_w_4: bool) {
        self.base
            .assert_valid_variables(&[gate.a, gate.b, gate.c, gate.d]);
        self.blocks
            .arithmetic
            .block
            .populate_wires(gate.a, gate.b, gate.c, gate.d);
        let mul_scaling = if use_next_gate_w_4 {
            gate.mul_scaling * Field::from(2u64)
        } else {
            gate.mul_scaling
        };
        self.blocks.arithmetic.block.q_m_mut().push(mul_scaling);
        self.blocks
            .arithmetic
            .block
            .q_1_mut()
            .push(gate.a_scaling);
        self.blocks
            .arithmetic
            .block
            .q_2_mut()
            .push(gate.b_scaling);
        self.blocks
            .arithmetic
            .block
            .q_3_mut()
            .push(gate.c_scaling);
        self.blocks
            .arithmetic
            .block
            .q_c_mut()
            .push(gate.const_scaling);
        self.blocks
            .arithmetic
            .block
            .q_4_mut()
            .push(gate.d_scaling);
        let arith_val = if use_next_gate_w_4 {
            Field::from(2u64)
        } else {
            Field::from(1u64)
        };
        self.blocks.arithmetic.set_gate_selector(arith_val);
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
    }

    /// Create a boolean gate: constrains `variable_index` to be 0 or 1.
    ///
    /// Generates the relation `x^2 - x = 0` via `q_m=1, q_1=-1` with both
    /// left and right wires set to the same variable.
    pub fn create_bool_gate(&mut self, variable_index: u32) {
        self.base.assert_valid_variables(&[variable_index]);
        let zero = self.base.zero_idx();
        self.blocks
            .arithmetic
            .block
            .populate_wires(variable_index, variable_index, zero, zero);
        self.blocks
            .arithmetic
            .block
            .q_m_mut()
            .push(Field::from(1u64));
        self.blocks
            .arithmetic
            .block
            .q_1_mut()
            .push(-Field::from(1u64));
        self.blocks
            .arithmetic
            .block
            .q_2_mut()
            .push(Field::zero());
        self.blocks
            .arithmetic
            .block
            .q_3_mut()
            .push(Field::zero());
        self.blocks
            .arithmetic
            .block
            .q_c_mut()
            .push(Field::zero());
        self.blocks
            .arithmetic
            .block
            .q_4_mut()
            .push(Field::zero());
        self.blocks
            .arithmetic
            .set_gate_selector(Field::from(1u64));
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
    }

    /// Create a standard arithmetic gate (3-wire, disabled 4th wire):
    /// `q_m*a*b + q_l*a + q_r*b + q_o*c + q_c = 0`.
    pub fn create_arithmetic_gate(&mut self, gate: &ArithmeticTriple<P>) {
        self.base.assert_valid_variables(&[gate.a, gate.b, gate.c]);
        let zero = self.base.zero_idx();
        self.blocks
            .arithmetic
            .block
            .populate_wires(gate.a, gate.b, gate.c, zero);
        self.blocks.arithmetic.block.q_m_mut().push(gate.q_m);
        self.blocks.arithmetic.block.q_1_mut().push(gate.q_l);
        self.blocks.arithmetic.block.q_2_mut().push(gate.q_r);
        self.blocks.arithmetic.block.q_3_mut().push(gate.q_o);
        self.blocks.arithmetic.block.q_c_mut().push(gate.q_c);
        self.blocks
            .arithmetic
            .block
            .q_4_mut()
            .push(Field::zero());
        self.blocks
            .arithmetic
            .set_gate_selector(Field::from(1u64));
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
    }

    /// Create a gate with no constraints (all selectors zero) but with
    /// potentially non-trivial wire values.
    ///
    /// Used as a "dummy" gate, e.g. to provide wire values accessible via
    /// shifts by the preceding gate.
    pub fn create_unconstrained_gate(
        &mut self,
        block_idx: UltraBlockIndex,
        idx_1: u32,
        idx_2: u32,
        idx_3: u32,
        idx_4: u32,
    ) {
        let block = self.get_block_mut(block_idx);
        block.block.populate_wires(idx_1, idx_2, idx_3, idx_4);
        block.block.q_m_mut().push(Field::zero());
        block.block.q_1_mut().push(Field::zero());
        block.block.q_2_mut().push(Field::zero());
        block.block.q_3_mut().push(Field::zero());
        block.block.q_c_mut().push(Field::zero());
        block.block.q_4_mut().push(Field::zero());
        block.set_gate_selector(Field::zero());
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Witness fixing / constant assertion
    // ════════════════════════════════════════════════════════════════════

    /// Fix a witness to a particular value by adding a gate
    /// `1 * witness + (-witness_value) = 0`.
    pub fn fix_witness(&mut self, witness_index: u32, witness_value: Field<P>) {
        self.base.assert_valid_variables(&[witness_index]);
        let zero = self.base.zero_idx();
        self.blocks
            .arithmetic
            .block
            .populate_wires(witness_index, zero, zero, zero);
        self.blocks
            .arithmetic
            .block
            .q_m_mut()
            .push(Field::zero());
        self.blocks
            .arithmetic
            .block
            .q_1_mut()
            .push(Field::from(1u64));
        self.blocks
            .arithmetic
            .block
            .q_2_mut()
            .push(Field::zero());
        self.blocks
            .arithmetic
            .block
            .q_3_mut()
            .push(Field::zero());
        self.blocks
            .arithmetic
            .block
            .q_c_mut()
            .push(-witness_value);
        self.blocks
            .arithmetic
            .block
            .q_4_mut()
            .push(Field::zero());
        self.blocks
            .arithmetic
            .set_gate_selector(Field::from(1u64));
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
    }

    /// Assert that a variable equals a constant value.
    ///
    /// Creates (or retrieves) a constant variable with value `b`, then
    /// merges the equivalence classes of `a_idx` and the constant variable.
    pub fn assert_equal_constant(&mut self, a_idx: u32, b: Field<P>, msg: &str) {
        if self.base.get_variable(a_idx) != b && !self.base.failed() {
            self.base.failure(msg.to_string());
        }
        let b_idx = self.put_constant_variable(b);
        self.base.assert_equal(a_idx, b_idx, msg);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Gate creation: elliptic curve
    // ════════════════════════════════════════════════════════════════════

    /// Create an elliptic curve addition gate.
    ///
    /// Adds either one or two gates. The gate pair has the following structure:
    ///
    /// ```text
    ///     | q_ecc | w1  | w2  | w3  | w4  |
    ///     |-------|-----|-----|-----|-----|
    ///     |    1  |  -  | x1  | y1  |  -  | --> constrained
    ///     |    0  | x2  | x3  | y3  | y2  | --> "unconstrained" (read via shifts)
    /// ```
    ///
    /// If the output of the previous gate matches the input of the current gate
    /// (`(x3,y3)_{i-1} == (x1,y1)_i`), the two are fused by setting selector values
    /// on the previous gate.
    pub fn create_ecc_add_gate(&mut self, gate: &EccAddGate<P>) {
        self.base
            .assert_valid_variables(&[gate.x1, gate.x2, gate.x3, gate.y1, gate.y2, gate.y3]);

        let block = &mut self.blocks.elliptic;

        // Can we fuse into the previous gate in the block?
        let can_fuse = block.size() > 0
            && block.block.w_r()[block.size() - 1] == gate.x1
            && block.block.w_o()[block.size() - 1] == gate.y1;

        if can_fuse {
            // Set q_sign (q_1) and q_elliptic of previous gate
            let idx = block.size() - 1;
            block.block.q_1_mut().set(idx, gate.sign_coefficient);
            block.q_elliptic_mut().set(idx, Field::from(1u64));
        } else {
            let zero = self.base.zero_idx();
            block.block.populate_wires(zero, gate.x1, gate.y1, zero);
            block.block.q_3_mut().push(Field::zero());
            block.block.q_4_mut().push(Field::zero());
            block.block.q_1_mut().push(gate.sign_coefficient);
            block.block.q_2_mut().push(Field::zero());
            block.block.q_m_mut().push(Field::zero());
            block.block.q_c_mut().push(Field::zero());
            block.set_gate_selector(Field::from(1u64));
            self.check_selector_length_consistency();
            self.base.increment_num_gates(1);
        }
        // Create the unconstrained gate with the output to be read via shifts
        self.create_unconstrained_gate(
            UltraBlockIndex::Elliptic,
            gate.x2,
            gate.x3,
            gate.y3,
            gate.y2,
        );
    }

    /// Create an elliptic curve doubling gate.
    ///
    /// Adds either one or two gates. The gate pair has the following structure:
    ///
    /// ```text
    ///     | q_ecc | w1  | w2  | w3  | w4  |
    ///     |-------|-----|-----|-----|-----|
    ///     |    1  |  -  | x1  | y1  |  -  | --> constrained (q_m=1 for doubling)
    ///     |    0  |  -  | x3  | y3  |  -  | --> "unconstrained" (read via shifts)
    /// ```
    ///
    /// If the output of the previous gate matches the input of the current gate,
    /// the two are fused.
    pub fn create_ecc_dbl_gate(&mut self, gate: &EccDblGate) {
        self.base
            .assert_valid_variables(&[gate.x1, gate.x3, gate.y1, gate.y3]);

        let block = &mut self.blocks.elliptic;

        // Can we fuse into the previous gate?
        let can_fuse = block.size() > 0
            && block.block.w_r()[block.size() - 1] == gate.x1
            && block.block.w_o()[block.size() - 1] == gate.y1;

        if can_fuse {
            let idx = block.size() - 1;
            block.q_elliptic_mut().set(idx, Field::from(1u64));
            block.block.q_m_mut().set(idx, Field::from(1u64));
        } else {
            let zero = self.base.zero_idx();
            block.block.populate_wires(zero, gate.x1, gate.y1, zero);
            block.block.q_m_mut().push(Field::from(1u64));
            block.block.q_1_mut().push(Field::zero());
            block.block.q_2_mut().push(Field::zero());
            block.block.q_3_mut().push(Field::zero());
            block.block.q_c_mut().push(Field::zero());
            block.block.q_4_mut().push(Field::zero());
            block.set_gate_selector(Field::from(1u64));
            self.check_selector_length_consistency();
            self.base.increment_num_gates(1);
        }
        let zero = self.base.zero_idx();
        self.create_unconstrained_gate(UltraBlockIndex::Elliptic, zero, gate.x3, gate.y3, zero);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Gate creation: range constraints
    // ════════════════════════════════════════════════════════════════════

    /// Create a range constraint for a variable, ensuring its value lies in `[0, 2^num_bits - 1]`.
    ///
    /// For 1-bit ranges, creates a boolean gate. For ranges ≤ DEFAULT_PLOOKUP_RANGE_BITNUM bits,
    /// creates a direct range constraint. For larger ranges, decomposes into limbs.
    pub fn create_range_constraint(&mut self, variable_index: u32, num_bits: usize, msg: &str) {
        if num_bits == 1 {
            self.create_bool_gate(variable_index);
        } else if num_bits <= DEFAULT_PLOOKUP_RANGE_BITNUM as usize {
            // Temporary fix: add arithmetic gate so the variable appears in a wire
            // (otherwise the sorted set / non-sorted set size imbalance causes issues).
            self.create_arithmetic_gate(&ArithmeticTriple {
                a: variable_index,
                b: variable_index,
                c: variable_index,
                q_m: Field::zero(),
                q_l: Field::from(1u64),
                q_r: -Field::from(1u64),
                q_o: Field::zero(),
                q_c: Field::zero(),
            });
            self.create_new_range_constraint(
                variable_index,
                (1u64 << num_bits) - 1,
                msg,
            );
        } else {
            self.decompose_into_default_range(
                variable_index,
                num_bits as u64,
                DEFAULT_PLOOKUP_RANGE_BITNUM,
                msg,
            );
        }
    }

    /// Lower-level range constraint: constrain a variable to `[0, target_range]`.
    ///
    /// Creates or reuses a `RangeList` for `target_range`, then tags the variable
    /// as belonging to that range set.
    pub fn create_new_range_constraint(
        &mut self,
        variable_index: u32,
        target_range: u64,
        msg: &str,
    ) {
        let val = self.base.get_variable(variable_index);
        let val_standard = val.from_montgomery_form();
        let is_out_of_range = val_standard.data[0] > target_range
            || val_standard.data[1] != 0
            || val_standard.data[2] != 0
            || val_standard.data[3] != 0;
        if is_out_of_range && !self.base.failed() {
            self.base.failure(msg.to_string());
        }

        if !self.range_lists.contains_key(&target_range) {
            let list = self.create_range_list(target_range);
            self.range_lists.insert(target_range, list);
        }

        let existing_tag = self.base.real_variable_tags
            [self.base.real_variable_index[variable_index as usize] as usize];
        let list_range_tag = self.range_lists[&target_range].range_tag;

        if existing_tag != list_range_tag {
            if existing_tag != DEFAULT_TAG {
                // Variable already has a tag from a different range. Find that range.
                let mut found_tag = false;
                let mut existing_range = 0u64;
                for (range, rl) in &self.range_lists {
                    if rl.range_tag == existing_tag {
                        found_tag = true;
                        existing_range = *range;
                        break;
                    }
                }
                debug_assert!(found_tag, "Tag not found in any range list");
                if existing_range < target_range {
                    // Already has a more restrictive range check — nothing to do.
                    return;
                }
                // The new constraint is more restrictive: deep-copy the variable
                // and apply the range check to the copy.
                let copied_witness = self.base.add_variable(self.base.get_variable(variable_index));
                self.create_add_gate(&AddTriple {
                    a: variable_index,
                    b: copied_witness,
                    c: self.base.zero_idx(),
                    a_scaling: Field::from(1u64),
                    b_scaling: -Field::from(1u64),
                    c_scaling: Field::zero(),
                    const_scaling: Field::zero(),
                });
                self.create_new_range_constraint(copied_witness, target_range, msg);
                return;
            }
            self.assign_tag(variable_index, list_range_tag);
            self.range_lists
                .get_mut(&target_range)
                .unwrap()
                .variable_indices
                .push(variable_index);
        }
    }

    /// Decompose a variable into limbs of `target_range_bitnum` bits each,
    /// range-check each limb, and constrain the limbs to reconstruct the original.
    ///
    /// Returns the variable indices of the limbs.
    pub fn decompose_into_default_range(
        &mut self,
        variable_index: u32,
        num_bits: u64,
        target_range_bitnum: u64,
        msg: &str,
    ) -> Vec<u32> {
        self.base.assert_valid_variables(&[variable_index]);
        assert!(num_bits > 0, "num_bits must be > 0");

        let val = self.base.get_variable(variable_index);
        let val_standard = val.from_montgomery_form();

        // Check the value fits in num_bits
        let msb = Self::get_msb_u256(&val_standard.data);
        if msb >= num_bits as u32 && !self.base.failed() {
            self.base.failure(msg.to_string());
        }

        let sublimb_mask = (1u64 << target_range_bitnum) - 1;

        let has_remainder_bits = !num_bits.is_multiple_of(target_range_bitnum);
        let num_limbs = num_bits.div_ceil(target_range_bitnum);
        let last_limb_size = num_bits - ((num_bits / target_range_bitnum) * target_range_bitnum);
        let last_limb_range = if last_limb_size > 0 {
            (1u64 << last_limb_size) - 1
        } else {
            sublimb_mask
        };

        // Extract sublimbs from the value
        let mut sublimbs = Vec::with_capacity(num_limbs as usize);
        let mut accumulator_data = val_standard.data;
        for _ in 0..num_limbs {
            sublimbs.push(accumulator_data[0] & sublimb_mask);
            // Right shift by target_range_bitnum
            Self::shift_right_u256(&mut accumulator_data, target_range_bitnum as u32);
        }

        // Create variables and range-check each limb
        let mut sublimb_indices = Vec::with_capacity(num_limbs as usize);
        for (i, &sublimb_val) in sublimbs.iter().enumerate() {
            let limb_idx = self.base.add_variable(Field::from(sublimb_val));
            sublimb_indices.push(limb_idx);
            if (i == sublimbs.len() - 1) && has_remainder_bits {
                self.create_new_range_constraint(limb_idx, last_limb_range, msg);
            } else {
                self.create_new_range_constraint(limb_idx, sublimb_mask, msg);
            }
        }

        // Create arithmetic gates to reconstruct the original value from limbs
        let num_limb_triples = num_limbs.div_ceil(3);
        let leftovers = if num_limbs.is_multiple_of(3) {
            3u64
        } else {
            num_limbs % 3
        };

        let mut acc_data = val_standard.data;
        let mut accumulator_idx = variable_index;

        for i in 0..num_limb_triples {
            let real_limbs = [
                !(i == (num_limb_triples - 1) && leftovers < 1),
                !(i == (num_limb_triples - 1) && leftovers < 2),
                !(i == (num_limb_triples - 1) && leftovers < 3),
            ];

            let round_sublimbs = [
                if real_limbs[0] { sublimbs[(3 * i) as usize] } else { 0 },
                if real_limbs[1] { sublimbs[(3 * i + 1) as usize] } else { 0 },
                if real_limbs[2] { sublimbs[(3 * i + 2) as usize] } else { 0 },
            ];
            let new_limbs = [
                if real_limbs[0] {
                    sublimb_indices[(3 * i) as usize]
                } else {
                    self.base.zero_idx()
                },
                if real_limbs[1] {
                    sublimb_indices[(3 * i + 1) as usize]
                } else {
                    self.base.zero_idx()
                },
                if real_limbs[2] {
                    sublimb_indices[(3 * i + 2) as usize]
                } else {
                    self.base.zero_idx()
                },
            ];
            let shifts = [
                target_range_bitnum * (3 * i),
                target_range_bitnum * (3 * i + 1),
                target_range_bitnum * (3 * i + 2),
            ];

            // Compute new accumulator: acc - sum of (sublimb << shift)
            let mut new_acc_data = acc_data;
            for j in 0..3 {
                let mut shifted = [0u64; 4];
                shifted[0] = round_sublimbs[j];
                Self::shift_left_u256(&mut shifted, shifts[j] as u32);
                Self::sub_u256(&mut new_acc_data, &shifted);
            }

            // Compute shift scalings as field elements: 1 << shifts[j]
            let shift_scalings: [Field<P>; 3] = std::array::from_fn(|j| {
                let mut s = [0u64; 4];
                s[0] = 1;
                Self::shift_left_u256(&mut s, shifts[j] as u32);
                Field::from_limbs(s)
            });

            self.create_big_add_gate(
                &AddQuad {
                    a: new_limbs[0],
                    b: new_limbs[1],
                    c: new_limbs[2],
                    d: accumulator_idx,
                    a_scaling: shift_scalings[0],
                    b_scaling: shift_scalings[1],
                    c_scaling: shift_scalings[2],
                    d_scaling: -Field::from(1u64),
                    const_scaling: Field::zero(),
                },
                i != num_limb_triples - 1,
            );

            accumulator_idx = self.base.add_variable(Field::from_limbs(new_acc_data));
            acc_data = new_acc_data;
        }

        sublimb_indices
    }

    // ════════════════════════════════════════════════════════════════════
    //  Gate creation: sort constraints
    // ════════════════════════════════════════════════════════════════════

    /// Create a sort constraint for a sequence of variables.
    ///
    /// Checks that neighboring differences are at most 3 (delta-range check).
    /// The variable list length must be a multiple of 4 (NUM_WIRES).
    pub fn create_sort_constraint(&mut self, variable_index: &[u32]) {
        assert_eq!(
            variable_index.len() % NUM_WIRES,
            0,
            "sort constraint variable count must be a multiple of NUM_WIRES"
        );
        self.base.assert_valid_variables(variable_index);

        for i in (0..variable_index.len()).step_by(NUM_WIRES) {
            self.blocks.delta_range.block.populate_wires(
                variable_index[i],
                variable_index[i + 1],
                variable_index[i + 2],
                variable_index[i + 3],
            );
            self.base.increment_num_gates(1);
            self.blocks.delta_range.block.q_m_mut().push(Field::zero());
            self.blocks.delta_range.block.q_1_mut().push(Field::zero());
            self.blocks.delta_range.block.q_2_mut().push(Field::zero());
            self.blocks.delta_range.block.q_3_mut().push(Field::zero());
            self.blocks.delta_range.block.q_c_mut().push(Field::zero());
            self.blocks.delta_range.block.q_4_mut().push(Field::zero());
            self.blocks
                .delta_range
                .set_gate_selector(Field::from(1u64));
            self.check_selector_length_consistency();
        }
        // Dummy gate needed because of sort widget's check of next row
        let last = variable_index[variable_index.len() - 1];
        let zero = self.base.zero_idx();
        self.create_unconstrained_gate(UltraBlockIndex::DeltaRange, last, zero, zero, zero);
    }

    /// Create a sort constraint with edge values.
    ///
    /// Like `create_sort_constraint` but additionally constrains the first element to
    /// equal `start` and the last element to equal `end`.
    pub fn create_sort_constraint_with_edges(
        &mut self,
        variable_index: &[u32],
        start: Field<P>,
        end: Field<P>,
    ) {
        assert_eq!(
            variable_index.len() % NUM_WIRES,
            0,
            "sort constraint variable count must be a multiple of NUM_WIRES"
        );
        assert!(
            variable_index.len() > NUM_WIRES,
            "sort constraint with edges requires more than NUM_WIRES variables"
        );
        self.base.assert_valid_variables(variable_index);

        // Constrain first input == start
        self.create_add_gate(&AddTriple {
            a: variable_index[0],
            b: self.base.zero_idx(),
            c: self.base.zero_idx(),
            a_scaling: Field::from(1u64),
            b_scaling: Field::zero(),
            c_scaling: Field::zero(),
            const_scaling: -start,
        });

        // Enforce range check for all but the final row
        let len = variable_index.len();
        for i in (0..len - NUM_WIRES).step_by(NUM_WIRES) {
            self.blocks.delta_range.block.populate_wires(
                variable_index[i],
                variable_index[i + 1],
                variable_index[i + 2],
                variable_index[i + 3],
            );
            self.base.increment_num_gates(1);
            self.blocks.delta_range.block.q_m_mut().push(Field::zero());
            self.blocks.delta_range.block.q_1_mut().push(Field::zero());
            self.blocks.delta_range.block.q_2_mut().push(Field::zero());
            self.blocks.delta_range.block.q_3_mut().push(Field::zero());
            self.blocks.delta_range.block.q_c_mut().push(Field::zero());
            self.blocks.delta_range.block.q_4_mut().push(Field::zero());
            self.blocks
                .delta_range
                .set_gate_selector(Field::from(1u64));
            self.check_selector_length_consistency();
        }

        // Enforce range checks of last row
        if len > NUM_WIRES {
            self.blocks.delta_range.block.populate_wires(
                variable_index[len - 4],
                variable_index[len - 3],
                variable_index[len - 2],
                variable_index[len - 1],
            );
            self.base.increment_num_gates(1);
            self.blocks.delta_range.block.q_m_mut().push(Field::zero());
            self.blocks.delta_range.block.q_1_mut().push(Field::zero());
            self.blocks.delta_range.block.q_2_mut().push(Field::zero());
            self.blocks.delta_range.block.q_3_mut().push(Field::zero());
            self.blocks.delta_range.block.q_c_mut().push(Field::zero());
            self.blocks.delta_range.block.q_4_mut().push(Field::zero());
            self.blocks
                .delta_range
                .set_gate_selector(Field::from(1u64));
            self.check_selector_length_consistency();
        }

        // Dummy gate + constrain last == end
        let zero = self.base.zero_idx();
        self.create_unconstrained_gate(
            UltraBlockIndex::DeltaRange,
            variable_index[len - 1],
            zero,
            zero,
            zero,
        );
        self.create_add_gate(&AddTriple {
            a: variable_index[len - 1],
            b: self.base.zero_idx(),
            c: self.base.zero_idx(),
            a_scaling: Field::from(1u64),
            b_scaling: Field::zero(),
            c_scaling: Field::zero(),
            const_scaling: -end,
        });
    }

    /// Put variables in the witness that aren't already used elsewhere.
    ///
    /// Pads the list to a multiple of NUM_WIRES and creates unconstrained arithmetic gates.
    pub fn create_unconstrained_gates(&mut self, variable_index: &[u32]) {
        let mut padded_list: Vec<u32> = variable_index.to_vec();
        let padding = (NUM_WIRES - (padded_list.len() % NUM_WIRES)) % NUM_WIRES;
        let zero = self.base.zero_idx();
        for _ in 0..padding {
            padded_list.push(zero);
        }
        self.base.assert_valid_variables(variable_index);
        self.base.assert_valid_variables(&padded_list);

        for i in (0..padded_list.len()).step_by(NUM_WIRES) {
            self.create_unconstrained_gate(
                UltraBlockIndex::Arithmetic,
                padded_list[i],
                padded_list[i + 1],
                padded_list[i + 2],
                padded_list[i + 3],
            );
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Tag management / generalized permutation
    // ════════════════════════════════════════════════════════════════════

    /// Get a new unique tag value.
    pub fn get_new_tag(&mut self) -> u32 {
        self.base.current_tag += 1;
        self.base.current_tag
    }

    /// Assign a tag to a variable for the generalized permutation argument.
    pub fn assign_tag(&mut self, variable_index: u32, tag: u32) {
        debug_assert!(
            tag <= self.base.current_tag,
            "tag {} exceeds current_tag {}",
            tag,
            self.base.current_tag
        );
        let real_idx = self.base.real_variable_index[variable_index as usize];
        // If already tagged with this tag, return (can happen due to copy constraints)
        if self.base.real_variable_tags[real_idx as usize] == tag {
            return;
        }
        debug_assert_eq!(
            self.base.real_variable_tags[real_idx as usize],
            DEFAULT_TAG,
            "Variable already has a non-default tag"
        );
        self.base.real_variable_tags[real_idx as usize] = tag;
    }

    /// Set tau(tag_index) = tau_index in the tag permutation.
    pub fn set_tau_at_index(&mut self, tag_index: u32, tau_index: u32) {
        self.base.tau_mut().insert(tag_index, tau_index);
    }

    /// Add a transposition to the tag permutation: tau(a) = b, tau(b) = a.
    pub fn set_tau_transposition(&mut self, tag_index_1: u32, tag_index_2: u32) {
        self.set_tau_at_index(tag_index_1, tag_index_2);
        self.set_tau_at_index(tag_index_2, tag_index_1);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Range list management
    // ════════════════════════════════════════════════════════════════════

    /// Create a new range list for `[0, target_range]`.
    ///
    /// Allocates tag/tau pairs, creates variables for each multiple of
    /// `DEFAULT_PLOOKUP_RANGE_STEP_SIZE` up to `target_range`, and places them
    /// in the witness via unconstrained gates.
    pub fn create_range_list(&mut self, target_range: u64) -> RangeList {
        let range_tag = self.get_new_tag();
        let tau_tag = self.get_new_tag();
        self.set_tau_transposition(range_tag, tau_tag);

        let mut result = RangeList {
            target_range,
            range_tag,
            tau_tag,
            variable_indices: Vec::new(),
        };

        let num_multiples = target_range / DEFAULT_PLOOKUP_RANGE_STEP_SIZE;
        result
            .variable_indices
            .reserve(num_multiples as usize + 2);

        for i in 0..=num_multiples {
            let index = self
                .base
                .add_variable(Field::from(i * DEFAULT_PLOOKUP_RANGE_STEP_SIZE));
            result.variable_indices.push(index);
            self.assign_tag(index, result.range_tag);
        }
        // Add the target_range itself
        let index = self.base.add_variable(Field::from(target_range));
        result.variable_indices.push(index);
        self.assign_tag(index, result.range_tag);

        // Need unconstrained gates so these variables appear in the witness
        self.create_unconstrained_gates(&result.variable_indices);

        result
    }

    /// Process a range list into delta-range sort constraint gates.
    ///
    /// Replaces variable indices with real indices, deduplicates, sorts by value,
    /// pads to a multiple of NUM_WIRES, tags with tau_tag, and creates sort
    /// constraint gates with edge values `[0, target_range]`.
    pub fn process_range_list(&mut self, list: &mut RangeList) {
        self.base.assert_valid_variables(&list.variable_indices);
        assert!(!list.variable_indices.is_empty());

        // Replace with real variable indices
        for x in &mut list.variable_indices {
            *x = self.base.real_variable_index[*x as usize];
        }

        // Remove duplicates
        list.variable_indices.sort();
        list.variable_indices.dedup();

        // Extract values and sort
        let mut sorted_list: Vec<u32> = list
            .variable_indices
            .iter()
            .map(|&vi| {
                let field_elem = self.base.get_variable(vi);
                field_elem.from_montgomery_form().data[0] as u32
            })
            .collect();
        sorted_list.sort();

        // Pad to a multiple of NUM_WIRES and ensure size > NUM_WIRES
        let mut padding = (NUM_WIRES - (list.variable_indices.len() % NUM_WIRES)) % NUM_WIRES;
        if list.variable_indices.len() <= NUM_WIRES {
            padding += NUM_WIRES;
        }

        let mut indices = Vec::with_capacity(padding + sorted_list.len());
        let zero = self.base.zero_idx();
        for _ in 0..padding {
            indices.push(zero);
        }
        for &sorted_value in &sorted_list {
            let index = self.base.add_variable(Field::from(sorted_value as u64));
            self.assign_tag(index, list.tau_tag);
            indices.push(index);
        }

        self.create_sort_constraint_with_edges(
            &indices,
            Field::zero(),
            Field::from(list.target_range),
        );
    }

    /// Process all range lists. Called during finalization.
    pub fn process_range_lists(&mut self) {
        // Collect keys to avoid borrow issues
        let keys: Vec<u64> = self.range_lists.keys().cloned().collect();
        for key in keys {
            let mut list = self.range_lists.remove(&key).unwrap();
            self.process_range_list(&mut list);
            self.range_lists.insert(key, list);
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  u256 bit manipulation helpers (operating on [u64; 4])
    // ════════════════════════════════════════════════════════════════════

    /// Get the MSB position of a u256 represented as [u64; 4].
    /// Returns 0 if the value is 0.
    fn get_msb_u256(data: &[u64; 4]) -> u32 {
        for i in (0..4).rev() {
            if data[i] != 0 {
                return (i as u32) * 64 + (63 - data[i].leading_zeros());
            }
        }
        0
    }

    /// Shift a u256 (as [u64; 4]) right by `shift` bits in place.
    fn shift_right_u256(data: &mut [u64; 4], shift: u32) {
        if shift == 0 {
            return;
        }
        if shift >= 256 {
            *data = [0; 4];
            return;
        }
        let word_shift = (shift / 64) as usize;
        let bit_shift = shift % 64;

        if bit_shift == 0 {
            for i in 0..4 {
                data[i] = if i + word_shift < 4 {
                    data[i + word_shift]
                } else {
                    0
                };
            }
        } else {
            for i in 0..4 {
                let lo = if i + word_shift < 4 {
                    data[i + word_shift] >> bit_shift
                } else {
                    0
                };
                let hi = if i + word_shift + 1 < 4 {
                    data[i + word_shift + 1] << (64 - bit_shift)
                } else {
                    0
                };
                data[i] = lo | hi;
            }
        }
    }

    /// Shift a u256 (as [u64; 4]) left by `shift` bits in place.
    fn shift_left_u256(data: &mut [u64; 4], shift: u32) {
        if shift == 0 {
            return;
        }
        if shift >= 256 {
            *data = [0; 4];
            return;
        }
        let word_shift = (shift / 64) as usize;
        let bit_shift = shift % 64;

        if bit_shift == 0 {
            for i in (0..4).rev() {
                data[i] = if i >= word_shift {
                    data[i - word_shift]
                } else {
                    0
                };
            }
        } else {
            for i in (0..4).rev() {
                let hi = if i >= word_shift {
                    data[i - word_shift] << bit_shift
                } else {
                    0
                };
                let lo = if i >= word_shift + 1 {
                    data[i - word_shift - 1] >> (64 - bit_shift)
                } else {
                    0
                };
                data[i] = hi | lo;
            }
        }
    }

    /// Subtract `b` from `a` (wrapping) in place for u256 as [u64; 4].
    fn sub_u256(a: &mut [u64; 4], b: &[u64; 4]) {
        let mut borrow: u64 = 0;
        for i in 0..4 {
            let (diff1, borrow1) = a[i].overflowing_sub(b[i]);
            let (diff2, borrow2) = diff1.overflowing_sub(borrow);
            a[i] = diff2;
            borrow = (borrow1 as u64) + (borrow2 as u64);
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Finalization
    // ════════════════════════════════════════════════════════════════════

    /// Finalize the circuit.
    ///
    /// If `ensure_nonzero` is true, adds gates to ensure all polynomial
    /// columns have at least one non-zero coefficient (preventing commitment
    /// to the zero polynomial). Populates the public inputs block.
    ///
    /// This is a simplified version that omits ROM/RAM processing, range list
    /// processing, and non-native field multiplication processing (those will
    /// be added in later beads).
    pub fn finalize_circuit(&mut self, ensure_nonzero: bool) {
        if !self.circuit_finalized {
            if ensure_nonzero {
                self.add_gates_to_ensure_all_polys_are_non_zero();
            }
            // NOTE: process_non_native_field_multiplications() and ROM/RAM processing
            // are omitted — they'll be added in later beads.
            self.process_range_lists();
            self.populate_public_inputs_block();
            self.circuit_finalized = true;
        }
    }

    /// Copy public input indices into the public inputs trace block.
    fn populate_public_inputs_block(&mut self) {
        let zero = self.base.zero_idx();
        let public_inputs: Vec<u32> = self.base.public_inputs().to_vec();
        for &idx in &public_inputs {
            self.blocks
                .pub_inputs
                .block
                .populate_wires(idx, idx, zero, zero);
            // All selectors are zero for public input gates.
            for sel in self.blocks.pub_inputs.block.get_selectors_mut() {
                sel.push(Field::zero());
            }
            for sel in self.blocks.pub_inputs.gate_selectors.iter_mut() {
                sel.push(Field::zero());
            }
        }
    }

    /// Ensure all polynomials have at least one non-zero coefficient.
    ///
    /// Adds dummy gates to each block type so that no selector polynomial
    /// is identically zero, which would cause issues when committing.
    pub fn add_gates_to_ensure_all_polys_are_non_zero(&mut self) {
        let zero = self.base.zero_idx();

        // ── arithmetic block: q_m, q_1, q_2, q_3, q_4 all set to 1 ──
        self.blocks
            .arithmetic
            .block
            .populate_wires(zero, zero, zero, zero);
        self.blocks
            .arithmetic
            .block
            .q_m_mut()
            .push(Field::from(1u64));
        self.blocks
            .arithmetic
            .block
            .q_1_mut()
            .push(Field::from(1u64));
        self.blocks
            .arithmetic
            .block
            .q_2_mut()
            .push(Field::from(1u64));
        self.blocks
            .arithmetic
            .block
            .q_3_mut()
            .push(Field::from(1u64));
        self.blocks
            .arithmetic
            .block
            .q_4_mut()
            .push(Field::from(1u64));
        self.blocks
            .arithmetic
            .block
            .q_c_mut()
            .push(Field::zero());
        self.blocks
            .arithmetic
            .set_gate_selector(Field::zero());
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);

        // ── delta_range block ──
        self.blocks
            .delta_range
            .block
            .populate_wires(zero, zero, zero, zero);
        self.blocks
            .delta_range
            .block
            .q_m_mut()
            .push(Field::zero());
        self.blocks
            .delta_range
            .block
            .q_1_mut()
            .push(Field::zero());
        self.blocks
            .delta_range
            .block
            .q_2_mut()
            .push(Field::zero());
        self.blocks
            .delta_range
            .block
            .q_3_mut()
            .push(Field::zero());
        self.blocks
            .delta_range
            .block
            .q_4_mut()
            .push(Field::zero());
        self.blocks
            .delta_range
            .block
            .q_c_mut()
            .push(Field::zero());
        self.blocks
            .delta_range
            .set_gate_selector(Field::from(1u64));
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
        self.create_unconstrained_gate(UltraBlockIndex::DeltaRange, zero, zero, zero, zero);

        // ── elliptic block ──
        self.blocks
            .elliptic
            .block
            .populate_wires(zero, zero, zero, zero);
        self.blocks
            .elliptic
            .block
            .q_m_mut()
            .push(Field::zero());
        self.blocks
            .elliptic
            .block
            .q_1_mut()
            .push(Field::zero());
        self.blocks
            .elliptic
            .block
            .q_2_mut()
            .push(Field::zero());
        self.blocks
            .elliptic
            .block
            .q_3_mut()
            .push(Field::zero());
        self.blocks
            .elliptic
            .block
            .q_4_mut()
            .push(Field::zero());
        self.blocks
            .elliptic
            .block
            .q_c_mut()
            .push(Field::zero());
        self.blocks
            .elliptic
            .set_gate_selector(Field::from(1u64));
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
        self.create_unconstrained_gate(UltraBlockIndex::Elliptic, zero, zero, zero, zero);

        // ── memory block ──
        self.blocks
            .memory
            .block
            .populate_wires(zero, zero, zero, zero);
        self.blocks
            .memory
            .block
            .q_m_mut()
            .push(Field::zero());
        self.blocks
            .memory
            .block
            .q_1_mut()
            .push(Field::zero());
        self.blocks
            .memory
            .block
            .q_2_mut()
            .push(Field::zero());
        self.blocks
            .memory
            .block
            .q_3_mut()
            .push(Field::zero());
        self.blocks
            .memory
            .block
            .q_4_mut()
            .push(Field::zero());
        self.blocks
            .memory
            .block
            .q_c_mut()
            .push(Field::zero());
        self.blocks
            .memory
            .set_gate_selector(Field::from(1u64));
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
        self.create_unconstrained_gate(UltraBlockIndex::Memory, zero, zero, zero, zero);

        // ── non-native field block ──
        self.blocks.nnf.block.populate_wires(zero, zero, zero, zero);
        self.blocks.nnf.block.q_m_mut().push(Field::zero());
        self.blocks.nnf.block.q_1_mut().push(Field::zero());
        self.blocks.nnf.block.q_2_mut().push(Field::zero());
        self.blocks.nnf.block.q_3_mut().push(Field::zero());
        self.blocks.nnf.block.q_4_mut().push(Field::zero());
        self.blocks.nnf.block.q_c_mut().push(Field::zero());
        self.blocks.nnf.set_gate_selector(Field::from(1u64));
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
        self.create_unconstrained_gate(UltraBlockIndex::Nnf, zero, zero, zero, zero);

        // ── q_4 and q_c non-zero: big_add_gate with d=1, q_4=1, q_c=-1 ──
        let one_idx = self.put_constant_variable(Field::one());
        self.create_big_add_gate(
            &AddQuad {
                a: zero,
                b: zero,
                c: zero,
                d: one_idx,
                a_scaling: Field::zero(),
                b_scaling: Field::zero(),
                c_scaling: Field::zero(),
                d_scaling: Field::from(1u64),
                const_scaling: -Field::from(1u64),
            },
            false,
        );

        // ── poseidon2_external block ──
        self.blocks
            .poseidon2_external
            .block
            .populate_wires(zero, zero, zero, zero);
        self.blocks
            .poseidon2_external
            .block
            .q_m_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_external
            .block
            .q_1_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_external
            .block
            .q_2_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_external
            .block
            .q_3_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_external
            .block
            .q_c_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_external
            .block
            .q_4_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_external
            .set_gate_selector(Field::from(1u64));
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
        self.create_unconstrained_gate(UltraBlockIndex::Poseidon2External, zero, zero, zero, zero);

        // ── poseidon2_internal block ──
        self.blocks
            .poseidon2_internal
            .block
            .populate_wires(zero, zero, zero, zero);
        self.blocks
            .poseidon2_internal
            .block
            .q_m_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_internal
            .block
            .q_1_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_internal
            .block
            .q_2_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_internal
            .block
            .q_3_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_internal
            .block
            .q_c_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_internal
            .block
            .q_4_mut()
            .push(Field::zero());
        self.blocks
            .poseidon2_internal
            .set_gate_selector(Field::from(1u64));
        self.check_selector_length_consistency();
        self.base.increment_num_gates(1);
        self.create_unconstrained_gate(UltraBlockIndex::Poseidon2Internal, zero, zero, zero, zero);

        // NOTE: Lookup table dummy gates are omitted — plookup will be added in a later bead.
    }

    // ════════════════════════════════════════════════════════════════════
    //  Selector consistency check
    // ════════════════════════════════════════════════════════════════════

    /// Debug helper: asserts all selectors in every block have the same length.
    ///
    /// Only active in debug builds. Each gate construction method appends values
    /// to all selectors; missing an update would cause an unsatisfiable circuit.
    pub fn check_selector_length_consistency(&self) {
        if cfg!(debug_assertions) {
            for block in self.blocks.get() {
                let all_selectors = block.get_all_selectors();
                if all_selectors.is_empty() {
                    continue;
                }
                let nominal_size = all_selectors[0].len();
                for (idx, sel) in all_selectors.iter().enumerate().skip(1) {
                    debug_assert_eq!(
                        sel.len(),
                        nominal_size,
                        "Selector {} has length {} but expected {} in block {:?}",
                        idx,
                        sel.len(),
                        nominal_size,
                        block.kind,
                    );
                }
            }
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Block access helper
    // ════════════════════════════════════════════════════════════════════

    /// Get a mutable reference to a trace block by index.
    fn get_block_mut(
        &mut self,
        idx: UltraBlockIndex,
    ) -> &mut crate::execution_trace::UltraTraceBlock<P> {
        match idx {
            UltraBlockIndex::PubInputs => &mut self.blocks.pub_inputs,
            UltraBlockIndex::Lookup => &mut self.blocks.lookup,
            UltraBlockIndex::Arithmetic => &mut self.blocks.arithmetic,
            UltraBlockIndex::DeltaRange => &mut self.blocks.delta_range,
            UltraBlockIndex::Elliptic => &mut self.blocks.elliptic,
            UltraBlockIndex::Memory => &mut self.blocks.memory,
            UltraBlockIndex::Nnf => &mut self.blocks.nnf,
            UltraBlockIndex::Poseidon2External => &mut self.blocks.poseidon2_external,
            UltraBlockIndex::Poseidon2Internal => &mut self.blocks.poseidon2_internal,
        }
    }
}

impl<P: FieldParams> Default for UltraCircuitBuilder<P> {
    fn default() -> Self {
        Self::new()
    }
}

/// Index enum for selecting which Ultra trace block to target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UltraBlockIndex {
    PubInputs,
    Lookup,
    Arithmetic,
    DeltaRange,
    Elliptic,
    Memory,
    Nnf,
    Poseidon2External,
    Poseidon2Internal,
}

// ════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::Bn254FrParams;

    type Fr = Field<Bn254FrParams>;
    type Builder = UltraCircuitBuilder<Bn254FrParams>;

    // ── Constructor tests ──────────────────────────────────────────

    #[test]
    fn test_constructor_creates_zero_variable() {
        let builder = Builder::new();
        // The constructor adds a zero constant variable and fixes it.
        assert!(builder.base.get_num_variables() > 0);
        let zero_val = builder.base.get_variable(builder.base.zero_idx());
        assert_eq!(zero_val, Fr::zero());
    }

    #[test]
    fn test_constructor_sets_tau_default() {
        let builder = Builder::new();
        assert_eq!(builder.base.tau().get(&DEFAULT_TAG), Some(&DEFAULT_TAG));
    }

    #[test]
    fn test_constructor_from_witnesses() {
        let witnesses = vec![Fr::from(10u64), Fr::from(20u64), Fr::from(30u64)];
        let pub_inputs = vec![0, 2];
        let builder = Builder::new_from_witnesses(&witnesses, pub_inputs, false);
        // The 3 witnesses + zero constant variable
        assert!(builder.base.get_num_variables() >= 4);
        assert_eq!(builder.base.get_variable(0), Fr::from(10u64));
        assert_eq!(builder.base.get_variable(1), Fr::from(20u64));
        assert_eq!(builder.base.get_variable(2), Fr::from(30u64));
        assert_eq!(builder.base.num_public_inputs(), 2);
    }

    // ── put_constant_variable tests ────────────────────────────────

    #[test]
    fn test_put_constant_variable_caches() {
        let mut builder = Builder::new();
        let idx1 = builder.put_constant_variable(Fr::from(42u64));
        let idx2 = builder.put_constant_variable(Fr::from(42u64));
        assert_eq!(idx1, idx2, "Same constant should return same index");
    }

    #[test]
    fn test_put_constant_variable_different_values() {
        let mut builder = Builder::new();
        let idx1 = builder.put_constant_variable(Fr::from(42u64));
        let idx2 = builder.put_constant_variable(Fr::from(99u64));
        assert_ne!(idx1, idx2, "Different constants should have different indices");
    }

    // ── create_add_gate tests ──────────────────────────────────────

    #[test]
    fn test_create_add_gate_basic() {
        let mut builder = Builder::new();
        let a_val = Fr::from(3u64);
        let b_val = Fr::from(7u64);
        let c_val = Fr::from(10u64);
        let a = builder.base.add_variable(a_val);
        let b = builder.base.add_variable(b_val);
        let c = builder.base.add_variable(c_val);
        // a*1 + b*1 + c*(-1) + 0 = 0  →  a + b - c = 0
        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::from(1u64),
            b_scaling: Fr::from(1u64),
            c_scaling: -Fr::from(1u64),
            const_scaling: Fr::zero(),
        });
        assert!(!builder.base.failed());
        assert!(builder.blocks.arithmetic.size() > 0);
    }

    #[test]
    fn test_create_add_gate_with_constant() {
        let mut builder = Builder::new();
        let a_val = Fr::from(5u64);
        let b_val = Fr::from(3u64);
        let c_val = Fr::from(10u64);
        let a = builder.base.add_variable(a_val);
        let b = builder.base.add_variable(b_val);
        let c = builder.base.add_variable(c_val);
        // a*1 + b*1 + c*(-1) + 2 = 0  →  a + b + 2 - c = 0  →  5+3+2-10=0
        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::from(1u64),
            b_scaling: Fr::from(1u64),
            c_scaling: -Fr::from(1u64),
            const_scaling: Fr::from(2u64),
        });
        assert!(!builder.base.failed());
    }

    // ── create_big_add_gate tests ──────────────────────────────────

    #[test]
    fn test_create_big_add_gate() {
        let mut builder = Builder::new();
        let a = builder.base.add_variable(Fr::from(1u64));
        let b = builder.base.add_variable(Fr::from(2u64));
        let c = builder.base.add_variable(Fr::from(3u64));
        let d = builder.base.add_variable(Fr::from(4u64));
        // 1 + 2 + 3 + 4 - 10 = 0
        builder.create_big_add_gate(
            &AddQuad {
                a,
                b,
                c,
                d,
                a_scaling: Fr::from(1u64),
                b_scaling: Fr::from(1u64),
                c_scaling: Fr::from(1u64),
                d_scaling: Fr::from(1u64),
                const_scaling: -Fr::from(10u64),
            },
            false,
        );
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_create_big_add_gate_use_next_gate_w4() {
        let mut builder = Builder::new();
        let a = builder.base.add_variable(Fr::from(1u64));
        let b = builder.base.add_variable(Fr::from(2u64));
        let c = builder.base.add_variable(Fr::from(3u64));
        let d = builder.base.add_variable(Fr::from(4u64));
        let initial_gates = builder.base.num_gates();
        builder.create_big_add_gate(
            &AddQuad {
                a,
                b,
                c,
                d,
                a_scaling: Fr::from(1u64),
                b_scaling: Fr::from(1u64),
                c_scaling: Fr::from(1u64),
                d_scaling: Fr::from(1u64),
                const_scaling: Fr::zero(),
            },
            true,
        );
        assert_eq!(builder.base.num_gates(), initial_gates + 1);
        // q_arith should be 2 when use_next_gate_w_4 is true
        let arith_block = &builder.blocks.arithmetic;
        let last_idx = arith_block.size() - 1;
        assert_eq!(arith_block.q_arith().get(last_idx), Fr::from(2u64));
    }

    // ── create_big_mul_add_gate tests ──────────────────────────────

    #[test]
    fn test_create_big_mul_add_gate() {
        let mut builder = Builder::new();
        let a = builder.base.add_variable(Fr::from(3u64));
        let b = builder.base.add_variable(Fr::from(4u64));
        let c = builder.base.add_variable(Fr::from(5u64));
        let d = builder.base.add_variable(Fr::from(6u64));
        // a*b*1 + a*0 + b*0 + c*0 + d*0 + (-12) = 3*4 - 12 = 0
        builder.create_big_mul_add_gate(
            &MulQuad {
                a,
                b,
                c,
                d,
                mul_scaling: Fr::from(1u64),
                a_scaling: Fr::zero(),
                b_scaling: Fr::zero(),
                c_scaling: Fr::zero(),
                d_scaling: Fr::zero(),
                const_scaling: -Fr::from(12u64),
            },
            false,
        );
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_create_big_mul_add_gate_with_next_w4() {
        let mut builder = Builder::new();
        let a = builder.base.add_variable(Fr::from(2u64));
        let b = builder.base.add_variable(Fr::from(3u64));
        let c = builder.base.add_variable(Fr::zero());
        let d = builder.base.add_variable(Fr::zero());
        builder.create_big_mul_add_gate(
            &MulQuad {
                a,
                b,
                c,
                d,
                mul_scaling: Fr::from(1u64),
                a_scaling: Fr::zero(),
                b_scaling: Fr::zero(),
                c_scaling: Fr::zero(),
                d_scaling: Fr::zero(),
                const_scaling: Fr::zero(),
            },
            true,
        );
        // q_m should be doubled (mul_scaling * 2 = 2)
        let arith_block = &builder.blocks.arithmetic;
        let last_idx = arith_block.size() - 1;
        assert_eq!(arith_block.block.q_m().get(last_idx), Fr::from(2u64));
    }

    // ── create_bool_gate tests ─────────────────────────────────────

    #[test]
    fn test_create_bool_gate_with_zero() {
        let mut builder = Builder::new();
        let a = builder.base.add_variable(Fr::zero());
        builder.create_bool_gate(a);
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_create_bool_gate_with_one() {
        let mut builder = Builder::new();
        let a = builder.base.add_variable(Fr::one());
        builder.create_bool_gate(a);
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_create_bool_gate_selectors() {
        let mut builder = Builder::new();
        let a = builder.base.add_variable(Fr::one());
        let initial_size = builder.blocks.arithmetic.size();
        builder.create_bool_gate(a);
        let idx = initial_size;
        assert_eq!(
            builder.blocks.arithmetic.block.q_m().get(idx),
            Fr::from(1u64)
        );
        assert_eq!(
            builder.blocks.arithmetic.block.q_1().get(idx),
            -Fr::from(1u64)
        );
    }

    // ── create_arithmetic_gate tests ───────────────────────────────

    #[test]
    fn test_create_arithmetic_gate() {
        let mut builder = Builder::new();
        let a = builder.base.add_variable(Fr::from(3u64));
        let b = builder.base.add_variable(Fr::from(4u64));
        let c = builder.base.add_variable(Fr::from(15u64));
        // q_m*a*b + q_l*a + q_r*b + q_o*c + q_c = 0
        // 1*3*4 + 1*3 + 0*4 + (-1)*15 + 0 = 12 + 3 - 15 = 0
        builder.create_arithmetic_gate(&ArithmeticTriple {
            a,
            b,
            c,
            q_m: Fr::from(1u64),
            q_l: Fr::from(1u64),
            q_r: Fr::zero(),
            q_o: -Fr::from(1u64),
            q_c: Fr::zero(),
        });
        assert!(!builder.base.failed());
    }

    // ── create_unconstrained_gate tests ────────────────────────────

    #[test]
    fn test_create_unconstrained_gate() {
        let mut builder = Builder::new();
        let a = builder.base.add_variable(Fr::from(10u64));
        let b = builder.base.add_variable(Fr::from(20u64));
        let initial_gates = builder.base.num_gates();
        builder.create_unconstrained_gate(UltraBlockIndex::Arithmetic, a, b, 0, 0);
        assert_eq!(builder.base.num_gates(), initial_gates + 1);
    }

    // ── fix_witness tests ──────────────────────────────────────────

    #[test]
    fn test_fix_witness() {
        let mut builder = Builder::new();
        let val = Fr::from(42u64);
        let idx = builder.base.add_variable(val);
        let initial_gates = builder.base.num_gates();
        builder.fix_witness(idx, val);
        assert_eq!(builder.base.num_gates(), initial_gates + 1);
        assert!(!builder.base.failed());
    }

    // ── assert_equal_constant tests ────────────────────────────────

    #[test]
    fn test_assert_equal_constant_matching() {
        let mut builder = Builder::new();
        let val = Fr::from(7u64);
        let idx = builder.base.add_variable(val);
        builder.assert_equal_constant(idx, val, "should match");
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_assert_equal_constant_mismatch() {
        let mut builder = Builder::new();
        let idx = builder.base.add_variable(Fr::from(7u64));
        builder.assert_equal_constant(idx, Fr::from(99u64), "mismatch");
        assert!(builder.base.failed());
        assert_eq!(builder.base.err(), "mismatch");
    }

    // ── finalize_circuit tests ─────────────────────────────────────

    #[test]
    fn test_finalize_circuit_sets_flag() {
        let mut builder = Builder::new();
        assert!(!builder.circuit_finalized);
        builder.finalize_circuit(false);
        assert!(builder.circuit_finalized);
    }

    #[test]
    fn test_finalize_circuit_idempotent() {
        let mut builder = Builder::new();
        builder.finalize_circuit(true);
        let gates_after_first = builder.base.num_gates();
        builder.finalize_circuit(true);
        assert_eq!(
            builder.base.num_gates(),
            gates_after_first,
            "Second finalize should not add more gates"
        );
    }

    #[test]
    fn test_finalize_with_ensure_nonzero() {
        let mut builder = Builder::new();
        let initial_gates = builder.base.num_gates();
        builder.finalize_circuit(true);
        assert!(
            builder.base.num_gates() > initial_gates,
            "ensure_nonzero should add gates"
        );
    }

    // ── add_gates_to_ensure_all_polys_are_non_zero tests ───────────

    #[test]
    fn test_add_gates_to_ensure_all_polys_nonzero_adds_to_all_blocks() {
        let mut builder = Builder::new();
        builder.add_gates_to_ensure_all_polys_are_non_zero();
        // Each of these blocks should have at least one gate
        assert!(builder.blocks.arithmetic.size() > 0);
        assert!(builder.blocks.delta_range.size() > 0);
        assert!(builder.blocks.elliptic.size() > 0);
        assert!(builder.blocks.memory.size() > 0);
        assert!(builder.blocks.nnf.size() > 0);
        assert!(builder.blocks.poseidon2_external.size() > 0);
        assert!(builder.blocks.poseidon2_internal.size() > 0);
    }

    // ── check_selector_length_consistency tests ────────────────────

    #[test]
    fn test_check_selector_length_consistency_passes_on_valid() {
        let mut builder = Builder::new();
        let a = builder.base.add_variable(Fr::from(1u64));
        let b = builder.base.add_variable(Fr::from(2u64));
        let c = builder.base.add_variable(Fr::from(3u64));
        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::from(1u64),
            b_scaling: Fr::from(1u64),
            c_scaling: -Fr::from(1u64),
            const_scaling: Fr::zero(),
        });
        // Should not panic
        builder.check_selector_length_consistency();
    }

    // ── populate_public_inputs_block tests ──────────────────────────

    #[test]
    fn test_populate_public_inputs_block() {
        let mut builder = Builder::new();
        let val = Fr::from(42u64);
        builder.base.add_public_variable(val);
        builder.finalize_circuit(false);
        // The pub_inputs block should have one entry
        assert_eq!(builder.blocks.pub_inputs.size(), 1);
    }

    // ═══════════════════════════════════════════════════════════════
    //  ECC gate tests
    // ═══════════════════════════════════════════════════════════════

    #[test]
    fn test_create_ecc_add_gate_basic() {
        let mut builder = Builder::new();
        // Create dummy point variables (structural test, not relation check)
        let x1 = builder.base.add_variable(Fr::from(1u64));
        let y1 = builder.base.add_variable(Fr::from(2u64));
        let x2 = builder.base.add_variable(Fr::from(3u64));
        let y2 = builder.base.add_variable(Fr::from(4u64));
        let x3 = builder.base.add_variable(Fr::from(5u64));
        let y3 = builder.base.add_variable(Fr::from(6u64));

        let initial_gates = builder.base.num_gates();
        builder.create_ecc_add_gate(&EccAddGate {
            x1,
            y1,
            x2,
            y2,
            x3,
            y3,
            sign_coefficient: Fr::from(1u64),
        });
        // Should add 2 gates (constrained + unconstrained)
        assert_eq!(builder.base.num_gates(), initial_gates + 2);
        assert_eq!(builder.blocks.elliptic.size(), 2);
    }

    #[test]
    fn test_create_ecc_add_gate_subtraction() {
        let mut builder = Builder::new();
        let x1 = builder.base.add_variable(Fr::from(10u64));
        let y1 = builder.base.add_variable(Fr::from(20u64));
        let x2 = builder.base.add_variable(Fr::from(30u64));
        let y2 = builder.base.add_variable(Fr::from(40u64));
        let x3 = builder.base.add_variable(Fr::from(50u64));
        let y3 = builder.base.add_variable(Fr::from(60u64));

        builder.create_ecc_add_gate(&EccAddGate {
            x1,
            y1,
            x2,
            y2,
            x3,
            y3,
            sign_coefficient: -Fr::from(1u64),
        });
        // q_1 (q_sign) of the first gate should be -1
        let q1_val = builder.blocks.elliptic.block.q_1().get(0);
        assert_eq!(q1_val, -Fr::from(1u64));
    }

    #[test]
    fn test_create_ecc_add_gate_chaining_fuses() {
        let mut builder = Builder::new();
        // First addition: (x1,y1) + (x2,y2) = (x3,y3)
        let x1 = builder.base.add_variable(Fr::from(1u64));
        let y1 = builder.base.add_variable(Fr::from(2u64));
        let x2 = builder.base.add_variable(Fr::from(3u64));
        let y2 = builder.base.add_variable(Fr::from(4u64));
        let x3 = builder.base.add_variable(Fr::from(5u64));
        let y3 = builder.base.add_variable(Fr::from(6u64));

        builder.create_ecc_add_gate(&EccAddGate {
            x1,
            y1,
            x2,
            y2,
            x3,
            y3,
            sign_coefficient: Fr::from(1u64),
        });
        assert_eq!(builder.blocks.elliptic.size(), 2);

        // Second addition: input (x3,y3) matches output of previous gate
        let x4 = builder.base.add_variable(Fr::from(7u64));
        let y4 = builder.base.add_variable(Fr::from(8u64));
        let x5 = builder.base.add_variable(Fr::from(9u64));
        let y5 = builder.base.add_variable(Fr::from(10u64));

        builder.create_ecc_add_gate(&EccAddGate {
            x1: x3, // Same as previous output
            y1: y3, // Same as previous output
            x2: x4,
            y2: y4,
            x3: x5,
            y3: y5,
            sign_coefficient: Fr::from(1u64),
        });
        // Should fuse: only 1 new unconstrained gate (3 total, not 4)
        assert_eq!(builder.blocks.elliptic.size(), 3);
    }

    #[test]
    fn test_create_ecc_add_gate_unchained() {
        let mut builder = Builder::new();
        // Two independent additions (can't fuse)
        for _ in 0..2 {
            let x1 = builder.base.add_variable(Fr::from(1u64));
            let y1 = builder.base.add_variable(Fr::from(2u64));
            let x2 = builder.base.add_variable(Fr::from(3u64));
            let y2 = builder.base.add_variable(Fr::from(4u64));
            let x3 = builder.base.add_variable(Fr::from(5u64));
            let y3 = builder.base.add_variable(Fr::from(6u64));
            builder.create_ecc_add_gate(&EccAddGate {
                x1,
                y1,
                x2,
                y2,
                x3,
                y3,
                sign_coefficient: Fr::from(1u64),
            });
        }
        // 2 unchained ops = 2 * 2 gates = 4
        assert_eq!(builder.blocks.elliptic.size(), 4);
    }

    #[test]
    fn test_create_ecc_dbl_gate_basic() {
        let mut builder = Builder::new();
        let x1 = builder.base.add_variable(Fr::from(1u64));
        let y1 = builder.base.add_variable(Fr::from(2u64));
        let x3 = builder.base.add_variable(Fr::from(3u64));
        let y3 = builder.base.add_variable(Fr::from(4u64));

        let initial_gates = builder.base.num_gates();
        builder.create_ecc_dbl_gate(&EccDblGate { x1, y1, x3, y3 });
        assert_eq!(builder.base.num_gates(), initial_gates + 2);
        assert_eq!(builder.blocks.elliptic.size(), 2);
        // q_m (q_is_double) should be 1
        assert_eq!(
            builder.blocks.elliptic.block.q_m().get(0),
            Fr::from(1u64)
        );
    }

    #[test]
    fn test_create_ecc_dbl_gate_chaining() {
        let mut builder = Builder::new();
        // First: add gate
        let x1 = builder.base.add_variable(Fr::from(1u64));
        let y1 = builder.base.add_variable(Fr::from(2u64));
        let x2 = builder.base.add_variable(Fr::from(3u64));
        let y2 = builder.base.add_variable(Fr::from(4u64));
        let x3 = builder.base.add_variable(Fr::from(5u64));
        let y3 = builder.base.add_variable(Fr::from(6u64));
        builder.create_ecc_add_gate(&EccAddGate {
            x1,
            y1,
            x2,
            y2,
            x3,
            y3,
            sign_coefficient: Fr::from(1u64),
        });
        assert_eq!(builder.blocks.elliptic.size(), 2);

        // Second: double the output of the first gate (should fuse)
        let x5 = builder.base.add_variable(Fr::from(7u64));
        let y5 = builder.base.add_variable(Fr::from(8u64));
        builder.create_ecc_dbl_gate(&EccDblGate {
            x1: x3,
            y1: y3,
            x3: x5,
            y3: y5,
        });
        // Fused: 3 gates total (add: 2, dbl fused: +1)
        assert_eq!(builder.blocks.elliptic.size(), 3);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Sort constraint tests
    // ═══════════════════════════════════════════════════════════════

    #[test]
    fn test_create_sort_constraint_basic() {
        let mut builder = Builder::new();
        let vals: Vec<u32> = (0..4)
            .map(|i| builder.base.add_variable(Fr::from((i + 1) as u64)))
            .collect();
        let initial_gates = builder.base.num_gates();
        builder.create_sort_constraint(&vals);
        // 1 delta_range gate + 1 dummy gate
        assert_eq!(builder.base.num_gates(), initial_gates + 2);
    }

    #[test]
    fn test_create_sort_constraint_8_elements() {
        let mut builder = Builder::new();
        let vals: Vec<u32> = (0..8)
            .map(|i| builder.base.add_variable(Fr::from((i + 1) as u64)))
            .collect();
        let initial_gates = builder.base.num_gates();
        builder.create_sort_constraint(&vals);
        // 2 delta_range gates + 1 dummy gate
        assert_eq!(builder.base.num_gates(), initial_gates + 3);
    }

    #[test]
    fn test_create_sort_constraint_with_edges_basic() {
        let mut builder = Builder::new();
        // 8 sorted values from 1 to 8
        let vals: Vec<u32> = (0..8)
            .map(|i| builder.base.add_variable(Fr::from((i + 1) as u64)))
            .collect();
        builder.create_sort_constraint_with_edges(&vals, Fr::from(1u64), Fr::from(8u64));
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_create_sort_constraint_with_edges_large() {
        let mut builder = Builder::new();
        // 28 sorted values with duplicates (from C++ test)
        let values: Vec<u64> = vec![
            1, 2, 5, 6, 7, 10, 11, 13, 16, 17, 20, 22, 22, 25, 26, 29, 29, 32, 32, 33, 35, 38,
            39, 39, 42, 42, 43, 45,
        ];
        let indices: Vec<u32> = values
            .iter()
            .map(|&v| builder.base.add_variable(Fr::from(v)))
            .collect();
        builder.create_sort_constraint_with_edges(&indices, Fr::from(1u64), Fr::from(45u64));
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_create_unconstrained_gates() {
        let mut builder = Builder::new();
        let vars: Vec<u32> = (0..7)
            .map(|i| builder.base.add_variable(Fr::from(i as u64)))
            .collect();
        let initial_gates = builder.base.num_gates();
        builder.create_unconstrained_gates(&vars);
        // 7 vars → padded to 8 → 2 unconstrained gates
        assert_eq!(builder.base.num_gates(), initial_gates + 2);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Tag management tests
    // ═══════════════════════════════════════════════════════════════

    #[test]
    fn test_get_new_tag() {
        let mut builder = Builder::new();
        let t1 = builder.get_new_tag();
        let t2 = builder.get_new_tag();
        assert_eq!(t1, 1);
        assert_eq!(t2, 2);
    }

    #[test]
    fn test_assign_tag() {
        let mut builder = Builder::new();
        let idx = builder.base.add_variable(Fr::from(42u64));
        let tag = builder.get_new_tag();
        builder.assign_tag(idx, tag);
        let real_idx = builder.base.real_variable_index[idx as usize];
        assert_eq!(builder.base.real_variable_tags[real_idx as usize], tag);
    }

    #[test]
    fn test_assign_tag_idempotent() {
        let mut builder = Builder::new();
        let idx = builder.base.add_variable(Fr::from(42u64));
        let tag = builder.get_new_tag();
        builder.assign_tag(idx, tag);
        builder.assign_tag(idx, tag); // Should not panic
        let real_idx = builder.base.real_variable_index[idx as usize];
        assert_eq!(builder.base.real_variable_tags[real_idx as usize], tag);
    }

    #[test]
    fn test_set_tau_transposition() {
        let mut builder = Builder::new();
        let t1 = builder.get_new_tag();
        let t2 = builder.get_new_tag();
        builder.set_tau_transposition(t1, t2);
        assert_eq!(builder.base.tau().get(&t1), Some(&t2));
        assert_eq!(builder.base.tau().get(&t2), Some(&t1));
    }

    // ═══════════════════════════════════════════════════════════════
    //  Range constraint tests
    // ═══════════════════════════════════════════════════════════════

    #[test]
    fn test_create_range_constraint_1_bit() {
        let mut builder = Builder::new();
        let idx = builder.base.add_variable(Fr::from(1u64));
        let initial_gates = builder.base.num_gates();
        builder.create_range_constraint(idx, 1, "1-bit range");
        // 1-bit range creates a bool gate (1 gate)
        assert_eq!(builder.base.num_gates(), initial_gates + 1);
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_create_range_constraint_small() {
        let mut builder = Builder::new();
        let idx = builder.base.add_variable(Fr::from(100u64));
        builder.create_range_constraint(idx, 8, "8-bit range");
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_create_new_range_constraint_basic() {
        let mut builder = Builder::new();
        let idx = builder.base.add_variable(Fr::from(5u64));
        builder.create_new_range_constraint(idx, 10, "test");
        assert!(!builder.base.failed());
        assert!(builder.range_lists.contains_key(&10));
    }

    #[test]
    fn test_create_new_range_constraint_out_of_range() {
        let mut builder = Builder::new();
        let idx = builder.base.add_variable(Fr::from(15u64));
        builder.create_new_range_constraint(idx, 10, "out of range");
        assert!(builder.base.failed());
        assert_eq!(builder.base.err(), "out of range");
    }

    #[test]
    fn test_create_new_range_constraint_reuses_list() {
        let mut builder = Builder::new();
        let idx1 = builder.base.add_variable(Fr::from(3u64));
        let idx2 = builder.base.add_variable(Fr::from(7u64));
        builder.create_new_range_constraint(idx1, 10, "test");
        builder.create_new_range_constraint(idx2, 10, "test");
        // Same range should reuse the same list
        assert_eq!(builder.range_lists.len(), 1);
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_create_range_list() {
        let mut builder = Builder::new();
        let list = builder.create_range_list(9);
        // Range 9, step 3: values 0, 3, 6, 9 (loop) + 9 (explicit) = 5 variables
        assert_eq!(list.target_range, 9);
        assert_eq!(list.variable_indices.len(), 5);
        assert!(list.range_tag > 0);
        assert!(list.tau_tag > 0);
        assert_ne!(list.range_tag, list.tau_tag);
    }

    #[test]
    fn test_decompose_into_default_range() {
        let mut builder = Builder::new();
        let val = Fr::from(12345u64);
        let idx = builder.base.add_variable(val);
        // Add an arithmetic gate so the variable appears in a wire
        builder.create_arithmetic_gate(&ArithmeticTriple {
            a: idx,
            b: idx,
            c: idx,
            q_m: Fr::zero(),
            q_l: Fr::from(1u64),
            q_r: -Fr::from(1u64),
            q_o: Fr::zero(),
            q_c: Fr::zero(),
        });
        let sublimb_indices =
            builder.decompose_into_default_range(idx, 16, DEFAULT_PLOOKUP_RANGE_BITNUM, "test");
        assert!(!builder.base.failed());
        // 16 bits with 14-bit range => 2 limbs
        assert_eq!(sublimb_indices.len(), 2);
    }

    #[test]
    fn test_decompose_into_default_range_large() {
        let mut builder = Builder::new();
        let val = Fr::from(1000000u64);
        let idx = builder.base.add_variable(val);
        builder.create_arithmetic_gate(&ArithmeticTriple {
            a: idx,
            b: idx,
            c: idx,
            q_m: Fr::zero(),
            q_l: Fr::from(1u64),
            q_r: -Fr::from(1u64),
            q_o: Fr::zero(),
            q_c: Fr::zero(),
        });
        let sublimb_indices =
            builder.decompose_into_default_range(idx, 28, DEFAULT_PLOOKUP_RANGE_BITNUM, "test");
        assert!(!builder.base.failed());
        // 28 bits / 14 = 2 limbs
        assert_eq!(sublimb_indices.len(), 2);
    }

    #[test]
    fn test_process_range_lists() {
        let mut builder = Builder::new();
        // Create a few range-constrained variables
        let idx1 = builder.base.add_variable(Fr::from(3u64));
        let idx2 = builder.base.add_variable(Fr::from(7u64));
        let idx3 = builder.base.add_variable(Fr::from(5u64));
        // Add them to the witness in arithmetic gates
        for &idx in &[idx1, idx2, idx3] {
            builder.create_arithmetic_gate(&ArithmeticTriple {
                a: idx,
                b: idx,
                c: idx,
                q_m: Fr::zero(),
                q_l: Fr::from(1u64),
                q_r: -Fr::from(1u64),
                q_o: Fr::zero(),
                q_c: Fr::zero(),
            });
        }
        builder.create_new_range_constraint(idx1, 10, "test");
        builder.create_new_range_constraint(idx2, 10, "test");
        builder.create_new_range_constraint(idx3, 10, "test");
        // Processing should create sort constraint gates
        let gates_before = builder.base.num_gates();
        builder.process_range_lists();
        assert!(builder.base.num_gates() > gates_before);
        assert!(!builder.base.failed());
    }

    #[test]
    fn test_finalize_processes_range_lists() {
        let mut builder = Builder::new();
        let idx = builder.base.add_variable(Fr::from(5u64));
        builder.create_arithmetic_gate(&ArithmeticTriple {
            a: idx,
            b: idx,
            c: idx,
            q_m: Fr::zero(),
            q_l: Fr::from(1u64),
            q_r: -Fr::from(1u64),
            q_o: Fr::zero(),
            q_c: Fr::zero(),
        });
        builder.create_new_range_constraint(idx, 10, "test");
        builder.finalize_circuit(false);
        assert!(builder.circuit_finalized);
        assert!(!builder.base.failed());
    }

    // ═══════════════════════════════════════════════════════════════
    //  u256 helper tests
    // ═══════════════════════════════════════════════════════════════

    #[test]
    fn test_shift_right_u256() {
        let mut data = [0xFF00u64, 0, 0, 0];
        Builder::shift_right_u256(&mut data, 8);
        assert_eq!(data[0], 0xFF);
    }

    #[test]
    fn test_shift_left_u256() {
        let mut data = [0xFFu64, 0, 0, 0];
        Builder::shift_left_u256(&mut data, 8);
        assert_eq!(data[0], 0xFF00);
    }

    #[test]
    fn test_shift_left_u256_cross_limb() {
        let mut data = [1u64, 0, 0, 0];
        Builder::shift_left_u256(&mut data, 64);
        assert_eq!(data[0], 0);
        assert_eq!(data[1], 1);
    }

    #[test]
    fn test_sub_u256() {
        let mut a = [10u64, 0, 0, 0];
        let b = [3u64, 0, 0, 0];
        Builder::sub_u256(&mut a, &b);
        assert_eq!(a[0], 7);
    }

    #[test]
    fn test_get_msb_u256() {
        assert_eq!(Builder::get_msb_u256(&[0, 0, 0, 0]), 0);
        assert_eq!(Builder::get_msb_u256(&[1, 0, 0, 0]), 0);
        assert_eq!(Builder::get_msb_u256(&[2, 0, 0, 0]), 1);
        assert_eq!(Builder::get_msb_u256(&[0xFF, 0, 0, 0]), 7);
        assert_eq!(Builder::get_msb_u256(&[0, 1, 0, 0]), 64);
    }
}
