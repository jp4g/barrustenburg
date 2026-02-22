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

// ════════════════════════════════════════════════════════════════════════
//  Memory (ROM/RAM) types
// ════════════════════════════════════════════════════════════════════════

/// Sentinel value for uninitialized memory cells.
pub const UNINITIALIZED_MEMORY_RECORD: u32 = u32::MAX;

/// Memory selector types for ROM/RAM gates.
///
/// Port of C++ `UltraCircuitBuilder_::MEMORY_SELECTORS`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemorySelector {
    MemNone,
    RamConsistencyCheck,
    RomConsistencyCheck,
    RamTimestampCheck,
    RomRead,
    RamRead,
    RamWrite,
}

/// Non-native field selector types.
///
/// Port of C++ `UltraCircuitBuilder_::NNF_SELECTORS`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NnfSelector {
    NnfNone,
    LimbAccumulate1,
    LimbAccumulate2,
    NonNativeField1,
    NonNativeField2,
    NonNativeField3,
}

/// ROM read record: stores index/value witness indices for a single ROM access.
///
/// Port of C++ `RomRecord`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RomRecord {
    pub index_witness: u32,
    pub value_column1_witness: u32,
    pub value_column2_witness: u32,
    pub index: u32,
    pub record_witness: u32,
    pub gate_index: usize,
}

impl PartialOrd for RomRecord {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RomRecord {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.index
            .cmp(&other.index)
            .then(self.gate_index.cmp(&other.gate_index))
    }
}

/// RAM access type: read or write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RamAccessType {
    Read,
    Write,
}

/// RAM read/write record: stores index/value/timestamp witness indices.
///
/// Port of C++ `RamRecord`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RamRecord {
    pub index_witness: u32,
    pub timestamp_witness: u32,
    pub value_witness: u32,
    pub index: u32,
    pub timestamp: u32,
    pub access_type: RamAccessType,
    pub record_witness: u32,
    pub gate_index: usize,
}

impl PartialOrd for RamRecord {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RamRecord {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.index
            .cmp(&other.index)
            .then(self.timestamp.cmp(&other.timestamp))
    }
}

/// ROM transcript: holds the state (values per cell) and access records.
///
/// Port of C++ `RomTranscript`.
#[derive(Debug, Clone)]
pub struct RomTranscript {
    /// State: each entry is [value_column1_witness, value_column2_witness].
    pub state: Vec<[u32; 2]>,
    /// Access records, sorted during finalization.
    pub records: Vec<RomRecord>,
}

/// RAM transcript: holds the state (current values) and access records.
///
/// Port of C++ `RamTranscript`.
#[derive(Debug, Clone)]
pub struct RamTranscript {
    /// Current state: one value witness per cell.
    pub state: Vec<u32>,
    /// Access records, sorted during finalization.
    pub records: Vec<RamRecord>,
    /// Running access counter (used as timestamp).
    pub access_count: u32,
}

// ════════════════════════════════════════════════════════════════════════
//  Lookup (plookup) types
// ════════════════════════════════════════════════════════════════════════

/// Column index for lookup table entries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnIdx {
    C1,
    C2,
    C3,
}

/// Lookup table read data: stores witness indices (or field values) for each column.
///
/// Port of C++ `plookup::ReadData<T>`.
#[derive(Debug, Clone)]
pub struct ReadData<T> {
    pub columns: [Vec<T>; 3],
    pub lookup_entries: Vec<Vec<T>>,
}

impl<T> Default for ReadData<T> {
    fn default() -> Self {
        Self {
            columns: [Vec::new(), Vec::new(), Vec::new()],
            lookup_entries: Vec::new(),
        }
    }
}

impl<T> std::ops::Index<ColumnIdx> for ReadData<T> {
    type Output = Vec<T>;
    fn index(&self, idx: ColumnIdx) -> &Self::Output {
        match idx {
            ColumnIdx::C1 => &self.columns[0],
            ColumnIdx::C2 => &self.columns[1],
            ColumnIdx::C3 => &self.columns[2],
        }
    }
}

impl<T> std::ops::IndexMut<ColumnIdx> for ReadData<T> {
    fn index_mut(&mut self, idx: ColumnIdx) -> &mut Self::Output {
        match idx {
            ColumnIdx::C1 => &mut self.columns[0],
            ColumnIdx::C2 => &mut self.columns[1],
            ColumnIdx::C3 => &mut self.columns[2],
        }
    }
}

/// Basic table identifier. Used to look up or create tables.
///
/// Port of C++ `plookup::BasicTableId`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BasicTableId(pub u64);

/// A basic lookup table: stores an ID, a table index, and recorded lookup entries.
///
/// Port of C++ `plookup::BasicTable` (simplified for circuit builder layer).
#[derive(Debug, Clone)]
pub struct BasicTable<P: FieldParams> {
    pub id: BasicTableId,
    pub table_index: u64,
    pub lookup_gates: Vec<Vec<Field<P>>>,
}

/// Multi-table identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MultiTableId(pub u64);

/// Multi-table: aggregates multiple basic tables with step sizes for accumulation.
///
/// Port of C++ `plookup::MultiTable` (simplified).
#[derive(Debug, Clone)]
pub struct MultiTable {
    pub id: MultiTableId,
    pub basic_table_ids: Vec<BasicTableId>,
    pub column_1_step_sizes: Vec<u64>,
    pub column_2_step_sizes: Vec<u64>,
    pub column_3_step_sizes: Vec<u64>,
}

// ════════════════════════════════════════════════════════════════════════
//  Non-native field types
// ════════════════════════════════════════════════════════════════════════

/// Cached partial non-native field multiplication.
///
/// Stores the 4-limb witness indices for two operands and intermediate results.
/// Used during finalization to create NNF gates.
///
/// Port of C++ `cached_partial_non_native_field_multiplication`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CachedPartialNnfMul {
    /// 4 limb witness indices for operand a.
    pub a: [u32; 4],
    /// 4 limb witness indices for operand b.
    pub b: [u32; 4],
    /// Intermediate lo_0 witness index.
    pub lo_0: u32,
    /// Intermediate hi_0 witness index.
    pub hi_0: u32,
    /// Intermediate hi_1 witness index.
    pub hi_1: u32,
}
