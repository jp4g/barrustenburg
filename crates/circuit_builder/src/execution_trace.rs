//! Execution trace block data structures for Ultra circuit builder.
//!
//! Port of `barretenberg/honk/execution_trace/execution_trace_block.hpp` and
//! `barretenberg/honk/execution_trace/ultra_execution_trace.hpp`.
//!
//! Provides the selector abstraction, execution trace blocks that hold wire and selector
//! data, and the Ultra-specific trace block layout with 9 specialized gate blocks.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

// ── Selector ──────────────────────────────────────────────────────────────────

/// A selector column in an execution trace block.
///
/// In C++ barretenberg, selectors use virtual dispatch (`Selector` → `ZeroSelector` |
/// `SlabVectorSelector`). In Rust we model this as an enum for zero-cost dispatch.
///
/// `Zero` selectors return `Field::zero()` for every index and only track their logical
/// size, saving memory for selector columns that are identically zero in a given block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Selector<P: FieldParams> {
    /// All entries are zero. Only tracks the logical length.
    Zero(usize),
    /// Backed by a Vec of field elements.
    Vec(Vec<Field<P>>),
}

impl<P: FieldParams> Selector<P> {
    /// Create a new zero selector with size 0.
    pub fn zero() -> Self {
        Selector::Zero(0)
    }

    /// Create a new vec-backed selector with size 0.
    pub fn vec() -> Self {
        Selector::Vec(Vec::new())
    }

    pub fn len(&self) -> usize {
        match self {
            Selector::Zero(n) => *n,
            Selector::Vec(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Push a field element. For `Zero` selectors, asserts the value is zero.
    pub fn push(&mut self, value: Field<P>) {
        match self {
            Selector::Zero(n) => {
                debug_assert!(
                    value == Field::zero(),
                    "cannot push non-zero value to ZeroSelector"
                );
                *n += 1;
            }
            Selector::Vec(v) => v.push(value),
        }
    }

    /// Push from an integer value. For `Zero`, asserts it is 0.
    pub fn push_int(&mut self, value: i32) {
        match self {
            Selector::Zero(n) => {
                debug_assert_eq!(value, 0, "cannot push non-zero value to ZeroSelector");
                *n += 1;
            }
            Selector::Vec(v) => v.push(Field::from(value as u64)),
        }
    }

    /// Set the value at a specific index.
    pub fn set(&mut self, idx: usize, value: Field<P>) {
        match self {
            Selector::Zero(_) => {
                debug_assert!(
                    value == Field::zero(),
                    "cannot set non-zero value in ZeroSelector"
                );
            }
            Selector::Vec(v) => v[idx] = value,
        }
    }

    /// Set the value at a specific index from an integer.
    pub fn set_int(&mut self, idx: usize, value: i32) {
        match self {
            Selector::Zero(_) => {
                debug_assert_eq!(value, 0, "cannot set non-zero value in ZeroSelector");
            }
            Selector::Vec(v) => v[idx] = Field::from(value as u64),
        }
    }

    /// Resize the selector to `new_size`, padding with zero.
    pub fn resize(&mut self, new_size: usize) {
        match self {
            Selector::Zero(n) => *n = new_size,
            Selector::Vec(v) => v.resize(new_size, Field::zero()),
        }
    }

    /// Get the element at `index`. Returns `Field::zero()` for `Zero` selectors.
    pub fn get(&self, index: usize) -> Field<P> {
        match self {
            Selector::Zero(n) => {
                debug_assert!(index < *n, "ZeroSelector index out of bounds");
                Field::zero()
            }
            Selector::Vec(v) => v[index],
        }
    }

    /// Get the last element.
    pub fn back(&self) -> Field<P> {
        match self {
            Selector::Zero(n) => {
                debug_assert!(*n > 0, "ZeroSelector is empty");
                Field::zero()
            }
            Selector::Vec(v) => *v.last().expect("selector is empty"),
        }
    }
}

impl<P: FieldParams> Default for Selector<P> {
    fn default() -> Self {
        Selector::Vec(Vec::new())
    }
}

// ── Wires ─────────────────────────────────────────────────────────────────────

/// Number of wires in the Ultra arithmetization.
pub const NUM_WIRES: usize = 4;

/// Wire columns: each wire is a vector of variable indices.
pub type Wires = [Vec<u32>; NUM_WIRES];

// ── ExecutionTraceBlock ───────────────────────────────────────────────────────

/// Base execution trace block holding wires and the 6 "non-gate" selectors
/// (q_m, q_c, q_1, q_2, q_3, q_4) that are always vec-backed.
///
/// Port of `ExecutionTraceBlock<FF, NUM_WIRES>`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExecutionTraceBlock<P: FieldParams> {
    /// Wire columns (variable indices).
    pub wires: Wires,
    /// Trace offset assigned during `compute_offsets`.
    pub trace_offset: u32,
    /// The 6 non-gate selectors: [q_m, q_c, q_1, q_2, q_3, q_4].
    pub non_gate_selectors: [Selector<P>; 6],
}

impl<P: FieldParams> Default for ExecutionTraceBlock<P> {
    fn default() -> Self {
        Self {
            wires: std::array::from_fn(|_| Vec::new()),
            trace_offset: u32::MAX,
            non_gate_selectors: std::array::from_fn(|_| Selector::vec()),
        }
    }
}

impl<P: FieldParams> ExecutionTraceBlock<P> {
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of gates (rows) in this block.
    pub fn size(&self) -> usize {
        if self.wires[0].is_empty() {
            0
        } else {
            self.wires[0].len()
        }
    }

    /// Reserve space in all wire and selector vectors.
    pub fn reserve(&mut self, size_hint: usize) {
        for wire in &mut self.wires {
            wire.reserve(size_hint);
        }
        for sel in &mut self.non_gate_selectors {
            if let Selector::Vec(v) = sel {
                v.reserve(size_hint);
            }
        }
    }

    /// Populate all 4 wires with a single row of variable indices.
    pub fn populate_wires(&mut self, idx_1: u32, idx_2: u32, idx_3: u32, idx_4: u32) {
        self.wires[0].push(idx_1);
        self.wires[1].push(idx_2);
        self.wires[2].push(idx_3);
        self.wires[3].push(idx_4);
    }

    // Convenience accessors matching C++ naming.

    pub fn w_l(&self) -> &Vec<u32> {
        &self.wires[0]
    }
    pub fn w_r(&self) -> &Vec<u32> {
        &self.wires[1]
    }
    pub fn w_o(&self) -> &Vec<u32> {
        &self.wires[2]
    }
    pub fn w_4(&self) -> &Vec<u32> {
        &self.wires[3]
    }

    pub fn w_l_mut(&mut self) -> &mut Vec<u32> {
        &mut self.wires[0]
    }
    pub fn w_r_mut(&mut self) -> &mut Vec<u32> {
        &mut self.wires[1]
    }
    pub fn w_o_mut(&mut self) -> &mut Vec<u32> {
        &mut self.wires[2]
    }
    pub fn w_4_mut(&mut self) -> &mut Vec<u32> {
        &mut self.wires[3]
    }

    // Non-gate selector accessors.

    pub fn q_m(&self) -> &Selector<P> {
        &self.non_gate_selectors[0]
    }
    pub fn q_c(&self) -> &Selector<P> {
        &self.non_gate_selectors[1]
    }
    pub fn q_1(&self) -> &Selector<P> {
        &self.non_gate_selectors[2]
    }
    pub fn q_2(&self) -> &Selector<P> {
        &self.non_gate_selectors[3]
    }
    pub fn q_3(&self) -> &Selector<P> {
        &self.non_gate_selectors[4]
    }
    pub fn q_4(&self) -> &Selector<P> {
        &self.non_gate_selectors[5]
    }

    pub fn q_m_mut(&mut self) -> &mut Selector<P> {
        &mut self.non_gate_selectors[0]
    }
    pub fn q_c_mut(&mut self) -> &mut Selector<P> {
        &mut self.non_gate_selectors[1]
    }
    pub fn q_1_mut(&mut self) -> &mut Selector<P> {
        &mut self.non_gate_selectors[2]
    }
    pub fn q_2_mut(&mut self) -> &mut Selector<P> {
        &mut self.non_gate_selectors[3]
    }
    pub fn q_3_mut(&mut self) -> &mut Selector<P> {
        &mut self.non_gate_selectors[4]
    }
    pub fn q_4_mut(&mut self) -> &mut Selector<P> {
        &mut self.non_gate_selectors[5]
    }

    /// Collect references to all selectors (non-gate + gate) for iteration.
    /// Substructs should override to include gate selectors.
    pub fn get_selectors(&self) -> Vec<&Selector<P>> {
        self.non_gate_selectors.iter().collect()
    }

    pub fn get_selectors_mut(&mut self) -> Vec<&mut Selector<P>> {
        self.non_gate_selectors.iter_mut().collect()
    }
}

// ── UltraTraceBlock ───────────────────────────────────────────────────────────

/// Identifies which gate-specific selector is active for an Ultra trace block.
///
/// Each specialized block in the C++ hierarchy overrides exactly one of the 8
/// gate selectors to be vec-backed while the remaining 7 are zero selectors.
/// We represent this with an enum rather than a class hierarchy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UltraBlockKind {
    PublicInput,
    Lookup,
    Arithmetic,
    DeltaRange,
    Elliptic,
    Memory,
    NonNativeField,
    Poseidon2External,
    Poseidon2Internal,
}

/// Index constants for the 8 gate-specific selectors within `gate_selectors`.
const Q_ARITH: usize = 0;
const Q_LOOKUP: usize = 1;
const Q_DELTA_RANGE: usize = 2;
const Q_ELLIPTIC: usize = 3;
const Q_MEMORY: usize = 4;
const Q_NNF: usize = 5;
const Q_POSEIDON2_EXTERNAL: usize = 6;
const Q_POSEIDON2_INTERNAL: usize = 7;

/// Number of gate-specific selectors in Ultra arithmetization.
pub const NUM_GATE_SELECTORS: usize = 8;

/// An Ultra execution trace block: combines the base `ExecutionTraceBlock` (6 non-gate
/// selectors + wires) with 8 gate-specific selectors. Exactly one gate selector is
/// vec-backed (determined by `kind`); the rest are zero selectors.
///
/// Port of `UltraTraceBlock` and its specialized subclasses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UltraTraceBlock<P: FieldParams> {
    pub block: ExecutionTraceBlock<P>,
    pub kind: UltraBlockKind,
    /// The 8 gate-specific selectors:
    /// [q_arith, q_lookup, q_delta_range, q_elliptic, q_memory, q_nnf,
    ///  q_poseidon2_external, q_poseidon2_internal]
    pub gate_selectors: [Selector<P>; NUM_GATE_SELECTORS],
}

impl<P: FieldParams> UltraTraceBlock<P> {
    /// Create a new Ultra trace block of the given kind.
    /// The gate selector corresponding to `kind` is initialized as vec-backed;
    /// all others are zero selectors.
    pub fn new(kind: UltraBlockKind) -> Self {
        let gate_selectors = std::array::from_fn(|i| {
            if i == Self::active_selector_index(kind) {
                Selector::vec()
            } else {
                Selector::zero()
            }
        });
        Self {
            block: ExecutionTraceBlock::default(),
            kind,
            gate_selectors,
        }
    }

    /// Return the gate selector index that is vec-backed for the given block kind.
    /// `PublicInput` has no active gate selector (all zero).
    fn active_selector_index(kind: UltraBlockKind) -> usize {
        match kind {
            // PublicInput has no gate selector — return a sentinel that won't match any index.
            UltraBlockKind::PublicInput => usize::MAX,
            UltraBlockKind::Arithmetic => Q_ARITH,
            UltraBlockKind::Lookup => Q_LOOKUP,
            UltraBlockKind::DeltaRange => Q_DELTA_RANGE,
            UltraBlockKind::Elliptic => Q_ELLIPTIC,
            UltraBlockKind::Memory => Q_MEMORY,
            UltraBlockKind::NonNativeField => Q_NNF,
            UltraBlockKind::Poseidon2External => Q_POSEIDON2_EXTERNAL,
            UltraBlockKind::Poseidon2Internal => Q_POSEIDON2_INTERNAL,
        }
    }

    /// Number of gates in this block.
    pub fn size(&self) -> usize {
        self.block.size()
    }

    /// Set the gate-specific selector value for a newly appended gate.
    /// This pushes to the active gate selector and pushes zero to all others.
    pub fn set_gate_selector(&mut self, value: Field<P>) {
        for (i, sel) in self.gate_selectors.iter_mut().enumerate() {
            if i == Self::active_selector_index(self.kind) {
                sel.push(value);
            } else {
                sel.push(Field::zero());
            }
        }
    }

    // Gate-specific selector accessors.

    pub fn q_arith(&self) -> &Selector<P> {
        &self.gate_selectors[Q_ARITH]
    }
    pub fn q_lookup(&self) -> &Selector<P> {
        &self.gate_selectors[Q_LOOKUP]
    }
    pub fn q_delta_range(&self) -> &Selector<P> {
        &self.gate_selectors[Q_DELTA_RANGE]
    }
    pub fn q_elliptic(&self) -> &Selector<P> {
        &self.gate_selectors[Q_ELLIPTIC]
    }
    pub fn q_memory(&self) -> &Selector<P> {
        &self.gate_selectors[Q_MEMORY]
    }
    pub fn q_nnf(&self) -> &Selector<P> {
        &self.gate_selectors[Q_NNF]
    }
    pub fn q_poseidon2_external(&self) -> &Selector<P> {
        &self.gate_selectors[Q_POSEIDON2_EXTERNAL]
    }
    pub fn q_poseidon2_internal(&self) -> &Selector<P> {
        &self.gate_selectors[Q_POSEIDON2_INTERNAL]
    }

    pub fn q_arith_mut(&mut self) -> &mut Selector<P> {
        &mut self.gate_selectors[Q_ARITH]
    }
    pub fn q_lookup_mut(&mut self) -> &mut Selector<P> {
        &mut self.gate_selectors[Q_LOOKUP]
    }
    pub fn q_delta_range_mut(&mut self) -> &mut Selector<P> {
        &mut self.gate_selectors[Q_DELTA_RANGE]
    }
    pub fn q_elliptic_mut(&mut self) -> &mut Selector<P> {
        &mut self.gate_selectors[Q_ELLIPTIC]
    }
    pub fn q_memory_mut(&mut self) -> &mut Selector<P> {
        &mut self.gate_selectors[Q_MEMORY]
    }
    pub fn q_nnf_mut(&mut self) -> &mut Selector<P> {
        &mut self.gate_selectors[Q_NNF]
    }
    pub fn q_poseidon2_external_mut(&mut self) -> &mut Selector<P> {
        &mut self.gate_selectors[Q_POSEIDON2_EXTERNAL]
    }
    pub fn q_poseidon2_internal_mut(&mut self) -> &mut Selector<P> {
        &mut self.gate_selectors[Q_POSEIDON2_INTERNAL]
    }

    /// Collect references to ALL selectors: 6 non-gate + 8 gate = 14 total.
    pub fn get_all_selectors(&self) -> Vec<&Selector<P>> {
        let mut sels: Vec<&Selector<P>> = self.block.non_gate_selectors.iter().collect();
        sels.extend(self.gate_selectors.iter());
        sels
    }
}

// ── UltraExecutionTraceBlocks ─────────────────────────────────────────────────

/// Number of block types in the Ultra arithmetization.
pub const NUM_ULTRA_BLOCKS: usize = 9;

/// Aggregation of all Ultra trace blocks: one per gate type.
///
/// Port of `UltraTraceBlockData` / `UltraExecutionTraceBlocks`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UltraExecutionTraceBlocks<P: FieldParams> {
    pub pub_inputs: UltraTraceBlock<P>,
    pub lookup: UltraTraceBlock<P>,
    pub arithmetic: UltraTraceBlock<P>,
    pub delta_range: UltraTraceBlock<P>,
    pub elliptic: UltraTraceBlock<P>,
    pub memory: UltraTraceBlock<P>,
    pub nnf: UltraTraceBlock<P>,
    pub poseidon2_external: UltraTraceBlock<P>,
    pub poseidon2_internal: UltraTraceBlock<P>,
}

impl<P: FieldParams> Default for UltraExecutionTraceBlocks<P> {
    fn default() -> Self {
        Self::new()
    }
}

impl<P: FieldParams> UltraExecutionTraceBlocks<P> {
    pub fn new() -> Self {
        Self {
            pub_inputs: UltraTraceBlock::new(UltraBlockKind::PublicInput),
            lookup: UltraTraceBlock::new(UltraBlockKind::Lookup),
            arithmetic: UltraTraceBlock::new(UltraBlockKind::Arithmetic),
            delta_range: UltraTraceBlock::new(UltraBlockKind::DeltaRange),
            elliptic: UltraTraceBlock::new(UltraBlockKind::Elliptic),
            memory: UltraTraceBlock::new(UltraBlockKind::Memory),
            nnf: UltraTraceBlock::new(UltraBlockKind::NonNativeField),
            poseidon2_external: UltraTraceBlock::new(UltraBlockKind::Poseidon2External),
            poseidon2_internal: UltraTraceBlock::new(UltraBlockKind::Poseidon2Internal),
        }
    }

    /// Return references to all blocks in order.
    pub fn get(&self) -> [&UltraTraceBlock<P>; NUM_ULTRA_BLOCKS] {
        [
            &self.pub_inputs,
            &self.lookup,
            &self.arithmetic,
            &self.delta_range,
            &self.elliptic,
            &self.memory,
            &self.nnf,
            &self.poseidon2_external,
            &self.poseidon2_internal,
        ]
    }

    /// Return mutable references to all blocks in order.
    pub fn get_mut(&mut self) -> [&mut UltraTraceBlock<P>; NUM_ULTRA_BLOCKS] {
        [
            &mut self.pub_inputs,
            &mut self.lookup,
            &mut self.arithmetic,
            &mut self.delta_range,
            &mut self.elliptic,
            &mut self.memory,
            &mut self.nnf,
            &mut self.poseidon2_external,
            &mut self.poseidon2_internal,
        ]
    }

    /// Return references to just the "gate blocks" (everything except pub_inputs).
    pub fn get_gate_blocks(&self) -> [&UltraTraceBlock<P>; NUM_ULTRA_BLOCKS - 1] {
        [
            &self.lookup,
            &self.arithmetic,
            &self.delta_range,
            &self.elliptic,
            &self.memory,
            &self.nnf,
            &self.poseidon2_external,
            &self.poseidon2_internal,
        ]
    }

    /// Compute trace offsets for each block. The first block starts at offset 1
    /// (row 0 is reserved for zero-row padding in Honk).
    pub fn compute_offsets(&mut self) {
        let mut offset: u32 = 1; // Row 0 is reserved
        for block in self.get_mut() {
            block.block.trace_offset = offset;
            offset += block.size() as u32;
        }
    }

    /// Total number of rows across all blocks.
    pub fn get_total_content_size(&self) -> usize {
        self.get().iter().map(|b| b.size()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::Bn254FrParams;

    type Fr = Field<Bn254FrParams>;

    #[test]
    fn test_zero_selector_returns_zero() {
        let sel: Selector<Bn254FrParams> = Selector::zero();
        assert!(sel.is_empty());
        assert_eq!(sel.len(), 0);
    }

    #[test]
    fn test_vec_selector_push_and_get() {
        let mut sel: Selector<Bn254FrParams> = Selector::vec();
        sel.push(Fr::from(42u64));
        sel.push(Fr::from(7u64));
        assert_eq!(sel.len(), 2);
        assert_eq!(sel.get(0), Fr::from(42u64));
        assert_eq!(sel.get(1), Fr::from(7u64));
        assert_eq!(sel.back(), Fr::from(7u64));
    }

    #[test]
    fn test_zero_selector_push_and_resize() {
        let mut sel: Selector<Bn254FrParams> = Selector::zero();
        sel.push(Fr::zero());
        sel.push(Fr::zero());
        assert_eq!(sel.len(), 2);
        assert_eq!(sel.get(0), Fr::zero());
        sel.resize(5);
        assert_eq!(sel.len(), 5);
    }

    #[test]
    fn test_execution_trace_block_populate_wires() {
        let mut block = ExecutionTraceBlock::<Bn254FrParams>::new();
        block.populate_wires(10, 20, 30, 40);
        block.populate_wires(11, 21, 31, 41);
        assert_eq!(block.size(), 2);
        assert_eq!(block.w_l(), &vec![10, 11]);
        assert_eq!(block.w_r(), &vec![20, 21]);
        assert_eq!(block.w_o(), &vec![30, 31]);
        assert_eq!(block.w_4(), &vec![40, 41]);
    }

    #[test]
    fn test_ultra_trace_block_kind_selectors() {
        let block = UltraTraceBlock::<Bn254FrParams>::new(UltraBlockKind::Arithmetic);
        // q_arith should be vec-backed
        assert!(matches!(block.q_arith(), Selector::Vec(_)));
        // all others should be zero
        assert!(matches!(block.q_lookup(), Selector::Zero(_)));
        assert!(matches!(block.q_delta_range(), Selector::Zero(_)));
        assert!(matches!(block.q_elliptic(), Selector::Zero(_)));
        assert!(matches!(block.q_memory(), Selector::Zero(_)));
        assert!(matches!(block.q_nnf(), Selector::Zero(_)));
        assert!(matches!(block.q_poseidon2_external(), Selector::Zero(_)));
        assert!(matches!(block.q_poseidon2_internal(), Selector::Zero(_)));
    }

    #[test]
    fn test_ultra_trace_block_public_input_all_zero_selectors() {
        let block = UltraTraceBlock::<Bn254FrParams>::new(UltraBlockKind::PublicInput);
        // PublicInput has no active gate selector — all 8 should be zero
        for sel in &block.gate_selectors {
            assert!(matches!(sel, Selector::Zero(_)));
        }
    }

    #[test]
    fn test_ultra_execution_trace_blocks_compute_offsets() {
        let mut blocks = UltraExecutionTraceBlocks::<Bn254FrParams>::new();
        // Add some gates to a few blocks
        blocks.arithmetic.block.populate_wires(0, 0, 0, 0);
        blocks.arithmetic.block.populate_wires(1, 1, 1, 1);
        blocks.lookup.block.populate_wires(2, 2, 2, 2);
        blocks.compute_offsets();
        // pub_inputs is empty but still gets offset 1
        assert_eq!(blocks.pub_inputs.block.trace_offset, 1);
        // lookup has 1 gate, so arithmetic starts at 1 + 0 (pub_inputs) + 1 (lookup)
        assert_eq!(blocks.lookup.block.trace_offset, 1); // after empty pub_inputs
        assert_eq!(blocks.arithmetic.block.trace_offset, 2); // after 1 lookup gate
        assert_eq!(blocks.delta_range.block.trace_offset, 4); // after 2 arithmetic gates
        assert_eq!(blocks.get_total_content_size(), 3);
    }

    #[test]
    fn test_ultra_execution_trace_blocks_get_gate_blocks() {
        let blocks = UltraExecutionTraceBlocks::<Bn254FrParams>::new();
        let gate_blocks = blocks.get_gate_blocks();
        assert_eq!(gate_blocks.len(), 8);
        assert_eq!(gate_blocks[0].kind, UltraBlockKind::Lookup);
        assert_eq!(gate_blocks[7].kind, UltraBlockKind::Poseidon2Internal);
    }

    #[test]
    fn test_set_gate_selector() {
        let mut block = UltraTraceBlock::<Bn254FrParams>::new(UltraBlockKind::Elliptic);
        block.set_gate_selector(Fr::from(1u64));
        // The elliptic selector should have the value
        assert_eq!(block.q_elliptic().get(0), Fr::from(1u64));
        // All others should be zero
        assert_eq!(block.q_arith().get(0), Fr::zero());
        assert_eq!(block.q_lookup().get(0), Fr::zero());
    }

    #[test]
    fn test_get_all_selectors_count() {
        let block = UltraTraceBlock::<Bn254FrParams>::new(UltraBlockKind::Memory);
        let all = block.get_all_selectors();
        assert_eq!(all.len(), 14); // 6 non-gate + 8 gate
    }
}
