//! CircuitBuilderBase — abstract base for all circuit builders.
//!
//! Port of `barretenberg/stdlib_circuit_builders/circuit_builder_base.hpp` and
//! `circuit_builder_base_impl.hpp`.
//!
//! Provides the variable storage, copy-constraint cycle management (via
//! `next_var_index` / `prev_var_index` / `real_variable_index`), public input
//! tracking, tag machinery for the generalized permutation argument, and
//! error-state bookkeeping.

use std::collections::HashMap;

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_numeric::bitop::get_msb64;

/// Default tag value — variables with this tag are NOT part of any multiset-equality check.
pub const DEFAULT_TAG: u32 = 0;

/// Sentinel: marks the *last* variable in its equivalence class (the "real" variable).
const REAL_VARIABLE: u32 = u32::MAX - 1;

/// Sentinel: marks the *first* variable in its equivalence class.
const FIRST_VARIABLE_IN_CLASS: u32 = u32::MAX - 2;

/// Base data and methods shared by all circuit builders.
///
/// In C++ this is a class template `CircuitBuilderBase<FF>`. Here it is a
/// concrete struct parameterized over `FieldParams`.
#[derive(Debug, Clone)]
pub struct CircuitBuilderBase<P: FieldParams> {
    // ── witness storage ─────────────────────────────────────────────────
    /// All witness values used by the circuit.
    variables: Vec<Field<P>>,

    /// Map from witness index to real variable index. Supports copy constraints:
    /// two witnesses with different indices can share a real variable index and
    /// thus the same witness value.
    pub real_variable_index: Vec<u32>,

    // ── copy-constraint cycle ───────────────────────────────────────────
    /// Index of the *next* variable in the equivalence-class cycle.
    /// `REAL_VARIABLE` if this is the last (real) element.
    next_var_index: Vec<u32>,

    /// Index of the *previous* variable in the equivalence-class cycle.
    /// `FIRST_VARIABLE_IN_CLASS` if this is the first element.
    prev_var_index: Vec<u32>,

    // ── public inputs ───────────────────────────────────────────────────
    /// Indices into `variables` that are designated as public inputs.
    public_inputs: Vec<u32>,

    /// Once set, no new public inputs may be added.
    public_inputs_finalized: bool,

    // ── tags ─────────────────────────────────────────────────────────────
    /// Per-real-variable tag for the generalized permutation argument.
    pub real_variable_tags: Vec<u32>,

    /// Current tag counter.
    pub current_tag: u32,

    /// Permutation on variable tags (key/values are real variable indices).
    tau: HashMap<u32, u32>,

    // ── bookkeeping ─────────────────────────────────────────────────────
    /// Index at which a witness constrained to equal 0 is stored.
    zero_idx: u32,

    /// Number of gates in the circuit.
    num_gates: usize,

    // ── error state ─────────────────────────────────────────────────────
    failed: bool,
    err: String,

    /// True when building a verification key (suppresses certain warnings).
    is_write_vk_mode: bool,
}

impl<P: FieldParams> CircuitBuilderBase<P> {
    // ════════════════════════════════════════════════════════════════════
    //  Construction
    // ════════════════════════════════════════════════════════════════════

    pub fn new() -> Self {
        Self {
            variables: Vec::new(),
            real_variable_index: Vec::new(),
            next_var_index: Vec::new(),
            prev_var_index: Vec::new(),
            public_inputs: Vec::new(),
            public_inputs_finalized: false,
            real_variable_tags: Vec::new(),
            current_tag: DEFAULT_TAG,
            tau: HashMap::new(),
            zero_idx: 0,
            num_gates: 0,
            failed: false,
            err: String::new(),
            is_write_vk_mode: false,
        }
    }

    pub fn new_with_vk_mode(is_write_vk_mode: bool) -> Self {
        Self {
            is_write_vk_mode,
            ..Self::new()
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Variable management
    // ════════════════════════════════════════════════════════════════════

    /// Add a variable (witness) to the circuit.
    ///
    /// Returns the index of the new variable.
    pub fn add_variable(&mut self, value: Field<P>) -> u32 {
        self.variables.push(value);
        let index = self.variables.len() as u32 - 1;
        self.real_variable_index.push(index);
        self.next_var_index.push(REAL_VARIABLE);
        self.prev_var_index.push(FIRST_VARIABLE_IN_CLASS);
        self.real_variable_tags.push(DEFAULT_TAG);
        index
    }

    /// Add a public variable: adds a variable and also registers it as a public input.
    pub fn add_public_variable(&mut self, value: Field<P>) -> u32 {
        let index = self.add_variable(value);
        assert!(
            !self.public_inputs_finalized,
            "Cannot add to public inputs after they have been finalized."
        );
        self.public_inputs.push(index);
        index
    }

    /// Make an existing witness variable public.
    ///
    /// Returns the position of this variable in the public inputs vector.
    pub fn set_public_input(&mut self, witness_index: u32) -> u32 {
        for &pi in &self.public_inputs {
            if pi == witness_index {
                if !self.failed {
                    self.failure("Attempted to set a public input that is already public!".into());
                }
                return 0;
            }
        }
        let public_input_index = self.public_inputs.len() as u32;
        assert!(
            !self.public_inputs_finalized,
            "Cannot add to public inputs after they have been finalized."
        );
        self.public_inputs.push(witness_index);
        public_input_index
    }

    /// Get the witness value for a variable, following copy constraints.
    #[inline]
    pub fn get_variable(&self, index: u32) -> Field<P> {
        debug_assert!((index as usize) < self.real_variable_index.len());
        let real_idx = self.real_variable_index[index as usize];
        debug_assert!((real_idx as usize) < self.variables.len());
        self.variables[real_idx as usize]
    }

    /// Set the witness value for a variable (only valid in write-vk mode).
    #[inline]
    pub fn set_variable(&mut self, index: u32, value: Field<P>) {
        assert!(self.is_write_vk_mode, "set_variable requires write-vk mode");
        debug_assert!((index as usize) < self.real_variable_index.len());
        let real_idx = self.real_variable_index[index as usize] as usize;
        debug_assert!(real_idx < self.variables.len());
        self.variables[real_idx] = value;
    }

    /// Set the witness value for a variable without mode checks.
    ///
    /// Used for testing scenarios where witness values need to be modified
    /// after circuit construction (e.g., verifying constraint soundness).
    #[inline]
    pub fn set_variable_unchecked(&mut self, index: u32, value: Field<P>) {
        let real_idx = self.real_variable_index[index as usize] as usize;
        self.variables[real_idx] = value;
    }

    /// Read-only access to all variables.
    pub fn get_variables(&self) -> &[Field<P>] {
        &self.variables
    }

    /// Number of variables (witnesses) in the circuit.
    pub fn get_num_variables(&self) -> usize {
        self.variables.len()
    }

    /// Alias used in C++ — same as `num_gates` but virtual in the original.
    pub fn get_num_finalized_gates(&self) -> usize {
        self.num_gates
    }

    // ════════════════════════════════════════════════════════════════════
    //  Copy-constraint / equivalence-class management
    // ════════════════════════════════════════════════════════════════════

    /// Join equivalence classes of `a_idx` and `b_idx`.
    ///
    /// After this call, `get_variable(a_idx)` and `get_variable(b_idx)` will
    /// return the same value. Sets `failed` if the current values disagree.
    pub fn assert_equal(&mut self, a_idx: u32, b_idx: u32, msg: &str) {
        self.assert_valid_variables(&[a_idx, b_idx]);

        let values_equal = self.get_variable(a_idx) == self.get_variable(b_idx);
        if !values_equal && !self.failed {
            self.failure(msg.to_string());
        }

        let a_real_idx = self.real_variable_index[a_idx as usize];
        let b_real_idx = self.real_variable_index[b_idx as usize];

        // Already in the same class — nothing to do.
        if a_real_idx == b_real_idx {
            return;
        }

        // Update real indices of the entire b-chain to point to a's real index.
        let b_start_idx = self.get_first_variable_in_class(b_idx);
        self.update_real_variable_indices(b_start_idx, a_real_idx);

        // Merge cycles: tie the last (real) element of b-chain to the first of a-chain.
        let a_start_idx = self.get_first_variable_in_class(a_idx);
        self.next_var_index[b_real_idx as usize] = a_start_idx;
        self.prev_var_index[a_start_idx as usize] = b_real_idx;

        // Check for tag clashes.
        let a_tag = self.real_variable_tags[a_real_idx as usize];
        let b_tag = self.real_variable_tags[b_real_idx as usize];
        let no_tag_clash =
            a_tag == DEFAULT_TAG || b_tag == DEFAULT_TAG || a_tag == b_tag;
        if !no_tag_clash && !self.failed {
            self.failure(msg.to_string());
        }
        if a_tag == DEFAULT_TAG {
            self.real_variable_tags[a_real_idx as usize] = b_tag;
        }
    }

    /// Walk backward through `prev_var_index` to find the first variable in the
    /// equivalence class containing `index`.
    pub fn get_first_variable_in_class(&self, mut index: u32) -> u32 {
        while self.prev_var_index[index as usize] != FIRST_VARIABLE_IN_CLASS {
            index = self.prev_var_index[index as usize];
        }
        index
    }

    /// Walk forward through `next_var_index`, setting all real-variable indices
    /// to `new_real_index`.
    fn update_real_variable_indices(&mut self, index: u32, new_real_index: u32) {
        let mut cur = index;
        loop {
            self.real_variable_index[cur as usize] = new_real_index;
            cur = self.next_var_index[cur as usize];
            if cur == REAL_VARIABLE {
                break;
            }
        }
    }

    /// Debug-mode check that all variable indices are valid.
    pub fn assert_valid_variables(&self, variable_indices: &[u32]) {
        if cfg!(debug_assertions) {
            for &idx in variable_indices {
                debug_assert!(
                    (idx as usize) < self.variables.len(),
                    "Variable index {} out of range (variables.len() = {})",
                    idx,
                    self.variables.len()
                );
            }
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Public inputs
    // ════════════════════════════════════════════════════════════════════

    /// Read-only access to the public inputs vector.
    pub fn public_inputs(&self) -> &[u32] {
        &self.public_inputs
    }

    /// Number of public inputs.
    pub fn num_public_inputs(&self) -> usize {
        self.public_inputs.len()
    }

    /// Prevent further additions to public inputs.
    pub fn finalize_public_inputs(&mut self) {
        self.public_inputs_finalized = true;
    }

    /// Directly initialize the public inputs vector.
    pub fn initialize_public_inputs(&mut self, public_inputs: Vec<u32>) {
        self.public_inputs = public_inputs;
    }

    // ════════════════════════════════════════════════════════════════════
    //  Gate count & subgroup size
    // ════════════════════════════════════════════════════════════════════

    /// Current number of gates.
    pub fn num_gates(&self) -> usize {
        self.num_gates
    }

    /// Increment the gate count.
    pub fn increment_num_gates(&mut self, count: usize) {
        self.num_gates += count;
    }

    /// Smallest power-of-two subgroup that fits `num_gates` gates.
    pub fn get_circuit_subgroup_size(&self, num_gates: usize) -> usize {
        if num_gates == 0 {
            return 1;
        }
        let mut log2_n = get_msb64(num_gates as u64) as usize;
        if (1usize << log2_n) != num_gates {
            log2_n += 1;
        }
        1usize << log2_n
    }

    // ════════════════════════════════════════════════════════════════════
    //  Tags / tau
    // ════════════════════════════════════════════════════════════════════

    /// Read-only access to the tag permutation.
    pub fn tau(&self) -> &HashMap<u32, u32> {
        &self.tau
    }

    /// Mutable access to the tag permutation.
    pub fn tau_mut(&mut self) -> &mut HashMap<u32, u32> {
        &mut self.tau
    }

    // ════════════════════════════════════════════════════════════════════
    //  Zero index
    // ════════════════════════════════════════════════════════════════════

    /// Index of the witness constrained to be zero.
    pub fn zero_idx(&self) -> u32 {
        self.zero_idx
    }

    /// Set the zero-variable index (used during circuit initialization).
    pub fn set_zero_idx(&mut self, value: u32) {
        self.zero_idx = value;
    }

    // ════════════════════════════════════════════════════════════════════
    //  Error state
    // ════════════════════════════════════════════════════════════════════

    /// Whether the circuit has encountered an error.
    pub fn failed(&self) -> bool {
        self.failed
    }

    /// The error message, if any.
    pub fn err(&self) -> &str {
        &self.err
    }

    /// Record a failure with the given message.
    pub fn failure(&mut self, msg: String) {
        self.failed = true;
        self.err = msg;
    }

    /// Whether the builder is in write-vk mode.
    pub fn is_write_vk_mode(&self) -> bool {
        self.is_write_vk_mode
    }
}

impl<P: FieldParams> Default for CircuitBuilderBase<P> {
    fn default() -> Self {
        Self::new()
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::Bn254FrParams;

    type Fr = Field<Bn254FrParams>;

    /// Port of C++ `UltraCircuitBuilder::BaseCase`: add a public variable, verify state.
    #[test]
    fn test_base_case() {
        let mut builder = CircuitBuilderBase::<Bn254FrParams>::new();
        let a = Fr::one();
        let idx = builder.add_public_variable(a);
        assert_eq!(idx, 0);
        assert_eq!(builder.get_variable(idx), Fr::one());
        assert_eq!(builder.num_public_inputs(), 1);
        assert_eq!(builder.public_inputs(), &[0]);
        assert!(!builder.failed());
    }

    /// Test add_variable and get_variable round-trip.
    #[test]
    fn test_add_and_get_variable() {
        let mut builder = CircuitBuilderBase::<Bn254FrParams>::new();
        let v0 = builder.add_variable(Fr::from(42u64));
        let v1 = builder.add_variable(Fr::from(99u64));
        assert_eq!(builder.get_variable(v0), Fr::from(42u64));
        assert_eq!(builder.get_variable(v1), Fr::from(99u64));
        assert_eq!(builder.get_num_variables(), 2);
    }

    /// Test assert_equal merges equivalence classes and updates real_variable_index.
    #[test]
    fn test_assert_equal_same_value() {
        let mut builder = CircuitBuilderBase::<Bn254FrParams>::new();
        let a = builder.add_variable(Fr::from(7u64));
        let b = builder.add_variable(Fr::from(7u64));

        // Before assert_equal, they have distinct real indices.
        assert_ne!(builder.real_variable_index[a as usize], builder.real_variable_index[b as usize]);

        builder.assert_equal(a, b, "test");

        // After assert_equal, they share the same real index.
        assert_eq!(builder.real_variable_index[a as usize], builder.real_variable_index[b as usize]);
        assert!(!builder.failed());
    }

    /// Test assert_equal with different values sets the failure flag.
    #[test]
    fn test_assert_equal_different_values() {
        let mut builder = CircuitBuilderBase::<Bn254FrParams>::new();
        let a = builder.add_variable(Fr::from(5u64));
        let b = builder.add_variable(Fr::from(10u64));

        builder.assert_equal(a, b, "values differ");

        assert!(builder.failed());
        assert_eq!(builder.err(), "values differ");
    }

    /// Test get_circuit_subgroup_size computes next power of two.
    #[test]
    fn test_get_circuit_subgroup_size() {
        let builder = CircuitBuilderBase::<Bn254FrParams>::new();
        assert_eq!(builder.get_circuit_subgroup_size(1), 1);
        assert_eq!(builder.get_circuit_subgroup_size(2), 2);
        assert_eq!(builder.get_circuit_subgroup_size(3), 4);
        assert_eq!(builder.get_circuit_subgroup_size(4), 4);
        assert_eq!(builder.get_circuit_subgroup_size(5), 8);
        assert_eq!(builder.get_circuit_subgroup_size(16), 16);
        assert_eq!(builder.get_circuit_subgroup_size(17), 32);
    }
}
