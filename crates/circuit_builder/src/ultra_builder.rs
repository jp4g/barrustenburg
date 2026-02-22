//! UltraCircuitBuilder — arithmetic gate methods for Ultra circuits.
//!
//! Port of `barretenberg/stdlib_circuit_builders/ultra_circuit_builder.hpp` and
//! `ultra_circuit_builder.cpp` (arithmetic gate subset).
//!
//! Provides the constructor, constant-variable caching, arithmetic gate creation
//! methods, witness fixing, finalization, and selector consistency checks.

use std::collections::HashMap;

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use crate::builder_base::{CircuitBuilderBase, DEFAULT_TAG};
use crate::execution_trace::UltraExecutionTraceBlocks;
use crate::gate_data::{AddQuad, AddTriple, ArithmeticTriple, MulQuad};

/// Ultra circuit builder: extends `CircuitBuilderBase` with execution trace blocks,
/// constant-variable caching, and arithmetic gate methods.
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
            // NOTE: process_non_native_field_multiplications(), ROM/RAM processing,
            // and range list processing are omitted — they'll be added in later beads.
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
}
