//! Plookup table type definitions.
//!
//! C++ source: plookup_tables/types.hpp
//!
//! Defines BasicTableId, MultiTableId, BasicTable, MultiTable, LookupHashTable, and ReadData.

use std::collections::HashMap;

use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};
use bbrs_ecc::fields::field::Field;
use bbrs_numeric::U256;

use super::fixed_base_params::FixedBaseParams;

// ---------------------------------------------------------------------------
// BasicTableId
// ---------------------------------------------------------------------------

/// Identifies a specific basic lookup table type.
///
/// Uses a newtype wrapper over `usize` because the C++ enum has computed
/// discriminants for the fixed-base table ranges (FIXED_BASE_0_0 through
/// FIXED_BASE_3_0 + NUM_TABLES_PER_HI_MULTITABLE).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BasicTableId(pub usize);

#[allow(non_upper_case_globals)]
impl BasicTableId {
    pub const XOR: Self = Self(0);
    pub const AND: Self = Self(1);
    pub const PEDERSEN: Self = Self(2);
    pub const AES_SPARSE_MAP: Self = Self(3);
    pub const AES_SBOX_MAP: Self = Self(4);
    pub const AES_SPARSE_NORMALIZE: Self = Self(5);
    pub const SHA256_WITNESS_NORMALIZE: Self = Self(6);
    pub const SHA256_WITNESS_SLICE_3: Self = Self(7);
    pub const SHA256_WITNESS_SLICE_7_ROTATE_4: Self = Self(8);
    pub const SHA256_WITNESS_SLICE_8_ROTATE_7: Self = Self(9);
    pub const SHA256_WITNESS_SLICE_14_ROTATE_1: Self = Self(10);
    pub const SHA256_CH_NORMALIZE: Self = Self(11);
    pub const SHA256_MAJ_NORMALIZE: Self = Self(12);
    pub const SHA256_BASE28: Self = Self(13);
    pub const SHA256_BASE28_ROTATE6: Self = Self(14);
    pub const SHA256_BASE28_ROTATE3: Self = Self(15);
    pub const SHA256_BASE16: Self = Self(16);
    pub const SHA256_BASE16_ROTATE2: Self = Self(17);
    pub const SHA256_BASE16_ROTATE6: Self = Self(18);
    pub const SHA256_BASE16_ROTATE7: Self = Self(19);
    pub const SHA256_BASE16_ROTATE8: Self = Self(20);
    pub const UINT_XOR_SLICE_6_ROTATE_0: Self = Self(21);
    pub const UINT_XOR_SLICE_2_ROTATE_0: Self = Self(22);
    pub const UINT_XOR_SLICE_4_ROTATE_0: Self = Self(23);
    pub const UINT_AND_SLICE_6_ROTATE_0: Self = Self(24);
    pub const UINT_AND_SLICE_2_ROTATE_0: Self = Self(25);
    pub const UINT_AND_SLICE_4_ROTATE_0: Self = Self(26);
    pub const SECP256K1_XLO_BASIC: Self = Self(27);
    pub const SECP256K1_XHI_BASIC: Self = Self(28);
    pub const SECP256K1_YLO_BASIC: Self = Self(29);
    pub const SECP256K1_YHI_BASIC: Self = Self(30);
    pub const SECP256K1_XYPRIME_BASIC: Self = Self(31);
    pub const SECP256K1_XLO_ENDO_BASIC: Self = Self(32);
    pub const SECP256K1_XHI_ENDO_BASIC: Self = Self(33);
    pub const SECP256K1_XYPRIME_ENDO_BASIC: Self = Self(34);
    pub const BLAKE_XOR_ROTATE0: Self = Self(35);
    pub const BLAKE_XOR_ROTATE0_SLICE5_MOD4: Self = Self(36);
    pub const BLAKE_XOR_ROTATE1: Self = Self(37);
    pub const BLAKE_XOR_ROTATE2: Self = Self(38);
    pub const BLAKE_XOR_ROTATE4: Self = Self(39);

    // Fixed-base table ranges: each range covers NUM_TABLES_PER_{LO,HI}_MULTITABLE entries
    pub const FIXED_BASE_0_0: Self = Self(40);
    pub const FIXED_BASE_1_0: Self =
        Self(Self::FIXED_BASE_0_0.0 + FixedBaseParams::NUM_TABLES_PER_LO_MULTITABLE as usize); // 55
    pub const FIXED_BASE_2_0: Self =
        Self(Self::FIXED_BASE_1_0.0 + FixedBaseParams::NUM_TABLES_PER_HI_MULTITABLE as usize); // 69
    pub const FIXED_BASE_3_0: Self =
        Self(Self::FIXED_BASE_2_0.0 + FixedBaseParams::NUM_TABLES_PER_LO_MULTITABLE as usize); // 84
    pub const HONK_DUMMY_BASIC1: Self =
        Self(Self::FIXED_BASE_3_0.0 + FixedBaseParams::NUM_TABLES_PER_HI_MULTITABLE as usize); // 98
    pub const HONK_DUMMY_BASIC2: Self = Self(Self::HONK_DUMMY_BASIC1.0 + 1);

    pub const KECCAK_INPUT: Self = Self(Self::HONK_DUMMY_BASIC2.0 + 1);
    pub const KECCAK_THETA: Self = Self(Self::KECCAK_INPUT.0 + 1);
    pub const KECCAK_RHO: Self = Self(Self::KECCAK_THETA.0 + 1);
    pub const KECCAK_CHI: Self = Self(Self::KECCAK_RHO.0 + 1);
    pub const KECCAK_OUTPUT: Self = Self(Self::KECCAK_CHI.0 + 1);
    pub const KECCAK_RHO_1: Self = Self(Self::KECCAK_OUTPUT.0 + 1);
    pub const KECCAK_RHO_2: Self = Self(Self::KECCAK_RHO_1.0 + 1);
    pub const KECCAK_RHO_3: Self = Self(Self::KECCAK_RHO_2.0 + 1);
    pub const KECCAK_RHO_4: Self = Self(Self::KECCAK_RHO_3.0 + 1);
    pub const KECCAK_RHO_5: Self = Self(Self::KECCAK_RHO_4.0 + 1);
    pub const KECCAK_RHO_6: Self = Self(Self::KECCAK_RHO_5.0 + 1);
    pub const KECCAK_RHO_7: Self = Self(Self::KECCAK_RHO_6.0 + 1);
    pub const KECCAK_RHO_8: Self = Self(Self::KECCAK_RHO_7.0 + 1);
    pub const KECCAK_RHO_9: Self = Self(Self::KECCAK_RHO_8.0 + 1);
}

// ---------------------------------------------------------------------------
// MultiTableId
// ---------------------------------------------------------------------------

/// Identifies a multi-table (a composition of basic tables).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MultiTableId(pub usize);

#[allow(non_upper_case_globals)]
impl MultiTableId {
    pub const SHA256_CH_INPUT: Self = Self(0);
    pub const SHA256_CH_OUTPUT: Self = Self(1);
    pub const SHA256_MAJ_INPUT: Self = Self(2);
    pub const SHA256_MAJ_OUTPUT: Self = Self(3);
    pub const SHA256_WITNESS_INPUT: Self = Self(4);
    pub const SHA256_WITNESS_OUTPUT: Self = Self(5);
    pub const AES_NORMALIZE: Self = Self(6);
    pub const AES_INPUT: Self = Self(7);
    pub const AES_SBOX: Self = Self(8);
    pub const FIXED_BASE_LEFT_LO: Self = Self(9);
    pub const FIXED_BASE_LEFT_HI: Self = Self(10);
    pub const FIXED_BASE_RIGHT_LO: Self = Self(11);
    pub const FIXED_BASE_RIGHT_HI: Self = Self(12);
    pub const UINT8_XOR: Self = Self(13);
    pub const UINT16_XOR: Self = Self(14);
    pub const UINT32_XOR: Self = Self(15);
    pub const UINT64_XOR: Self = Self(16);
    pub const UINT8_AND: Self = Self(17);
    pub const UINT16_AND: Self = Self(18);
    pub const UINT32_AND: Self = Self(19);
    pub const UINT64_AND: Self = Self(20);
    pub const SECP256K1_XLO: Self = Self(21);
    pub const SECP256K1_XHI: Self = Self(22);
    pub const SECP256K1_YLO: Self = Self(23);
    pub const SECP256K1_YHI: Self = Self(24);
    pub const SECP256K1_XYPRIME: Self = Self(25);
    pub const SECP256K1_XLO_ENDO: Self = Self(26);
    pub const SECP256K1_XHI_ENDO: Self = Self(27);
    pub const SECP256K1_XYPRIME_ENDO: Self = Self(28);
    pub const BLAKE_XOR: Self = Self(29);
    pub const BLAKE_XOR_ROTATE_16: Self = Self(30);
    pub const BLAKE_XOR_ROTATE_8: Self = Self(31);
    pub const BLAKE_XOR_ROTATE_7: Self = Self(32);
    pub const HONK_DUMMY_MULTI: Self = Self(33);
    pub const KECCAK_THETA_OUTPUT: Self = Self(34);
    pub const KECCAK_CHI_OUTPUT: Self = Self(35);
    pub const KECCAK_FORMAT_INPUT: Self = Self(36);
    pub const KECCAK_FORMAT_OUTPUT: Self = Self(37);
    pub const KECCAK_NORMALIZE_AND_ROTATE: Self = Self(38);

    /// Total number of multi-tables: KECCAK_NORMALIZE_AND_ROTATE + 25 (one per Keccak lane)
    pub const NUM_MULTI_TABLES: usize = Self::KECCAK_NORMALIZE_AND_ROTATE.0 + 25; // 63
}

// ---------------------------------------------------------------------------
// ColumnIdx
// ---------------------------------------------------------------------------

/// Index into the 3-column layout of a lookup table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnIdx {
    C1 = 0,
    C2 = 1,
    C3 = 2,
}

// ---------------------------------------------------------------------------
// LookupEntry
// ---------------------------------------------------------------------------

/// A single entry recording the data for a lookup gate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LookupEntry {
    /// Two key values (for 1:1, 1:2, or 2:1 lookup formats).
    pub key: [U256; 2],
    /// Two result values.
    pub value: [Fr; 2],
}

impl LookupEntry {
    pub fn new() -> Self {
        Self {
            key: [U256::ZERO, U256::ZERO],
            value: [Fr::zero(), Fr::zero()],
        }
    }

    /// Express the key-value pair as entries of a 3-column table row.
    pub fn to_table_components(&self, use_two_keys: bool) -> [Fr; 3] {
        [
            fr_from_u256(self.key[0]),
            if use_two_keys {
                fr_from_u256(self.key[1])
            } else {
                self.value[0]
            },
            if use_two_keys {
                self.value[0]
            } else {
                self.value[1]
            },
        ]
    }
}

impl PartialOrd for LookupEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for LookupEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.key[0]
            .cmp(&other.key[0])
            .then_with(|| self.key[1].cmp(&other.key[1]))
    }
}

// ---------------------------------------------------------------------------
// LookupHashTable
// ---------------------------------------------------------------------------

/// A map from (col1, col2, col3) entry to row index in a BasicTable.
///
/// Used to construct read_counts for the log-derivative lookup argument.
#[derive(Debug, Clone)]
pub struct LookupHashTable {
    index_map: HashMap<[u64; 3], usize>,
}

impl LookupHashTable {
    pub fn new() -> Self {
        Self {
            index_map: HashMap::new(),
        }
    }

    /// Initialize the entry-index map from the three columns.
    /// Uses the low limb of each field element as the hash key for simplicity.
    pub fn initialize(&mut self, column_1: &[Fr], column_2: &[Fr], column_3: &[Fr]) {
        self.index_map.clear();
        self.index_map.reserve(column_1.len());
        for i in 0..column_1.len() {
            let key = [column_1[i].data[0], column_2[i].data[0], column_3[i].data[0]];
            self.index_map.insert(key, i);
        }
    }

    /// Look up a row index by the 3-column entry.
    pub fn get(&self, key: &[Fr; 3]) -> Option<usize> {
        let hash_key = [key[0].data[0], key[1].data[0], key[2].data[0]];
        self.index_map.get(&hash_key).copied()
    }
}

impl Default for LookupHashTable {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// BasicTable
// ---------------------------------------------------------------------------

/// Function pointer type for computing table output from input keys.
pub type GetValuesFromKey = fn([u64; 2]) -> [Fr; 2];

/// A basic lookup table: stores the actual column data and lookup gate records.
///
/// Corresponds to C++ `BasicTable` in types.hpp.
#[derive(Clone)]
pub struct BasicTable {
    /// Unique ID for this table type.
    pub id: BasicTableId,
    /// Index of this table instance within the global table array.
    pub table_index: usize,
    /// Whether this table uses two input keys (2-to-1 lookup) rather than one (1-to-2).
    pub use_twin_keys: bool,

    /// Step sizes for accumulator computation in each column.
    pub column_1_step_size: Fr,
    pub column_2_step_size: Fr,
    pub column_3_step_size: Fr,

    /// Raw table data: three columns of field elements.
    pub column_1: Vec<Fr>,
    pub column_2: Vec<Fr>,
    pub column_3: Vec<Fr>,

    /// Wire data for all lookup gates created for lookups on this table.
    pub lookup_gates: Vec<LookupEntry>,

    /// Map from table entry to its row index (for constructing read counts).
    pub index_map: LookupHashTable,

    /// Function to compute output values from input key(s).
    pub get_values_from_key: GetValuesFromKey,
}

impl BasicTable {
    pub fn new() -> Self {
        Self {
            id: BasicTableId::XOR,
            table_index: 0,
            use_twin_keys: false,
            column_1_step_size: Fr::zero(),
            column_2_step_size: Fr::zero(),
            column_3_step_size: Fr::zero(),
            column_1: Vec::new(),
            column_2: Vec::new(),
            column_3: Vec::new(),
            lookup_gates: Vec::new(),
            index_map: LookupHashTable::new(),
            get_values_from_key: default_get_values,
        }
    }

    /// Initialize the index map from the current column data.
    pub fn initialize_index_map(&mut self) {
        self.index_map
            .initialize(&self.column_1, &self.column_2, &self.column_3);
    }

    /// Number of rows in this table.
    pub fn size(&self) -> usize {
        debug_assert_eq!(self.column_1.len(), self.column_2.len());
        debug_assert_eq!(self.column_2.len(), self.column_3.len());
        self.column_1.len()
    }
}

/// Default no-op values function (returns zeros).
pub fn default_get_values(_key: [u64; 2]) -> [Fr; 2] {
    [Fr::zero(), Fr::zero()]
}

// ---------------------------------------------------------------------------
// MultiTable
// ---------------------------------------------------------------------------

/// Container managing multiple BasicTables plus the data needed to combine
/// basic table outputs (e.g. limbs) into accumulators.
///
/// A MultiTable does not store raw table data. It stores basic table IDs,
/// the methods used to compute basic table entries, and metadata (coefficients
/// and step sizes) for reconstructing full values from components.
///
/// Corresponds to C++ `MultiTable` in types.hpp.
#[derive(Clone)]
pub struct MultiTable {
    /// Accumulated products of corresponding step sizes.
    pub column_1_coefficients: Vec<Fr>,
    pub column_2_coefficients: Vec<Fr>,
    pub column_3_coefficients: Vec<Fr>,

    /// Identifier for this multi-table.
    pub id: MultiTableId,

    /// IDs of the basic tables that compose this multi-table.
    pub basic_table_ids: Vec<BasicTableId>,

    /// Slice sizes for decomposing input keys.
    pub slice_sizes: Vec<u64>,

    /// Step sizes computed from coefficients (used in accumulator construction).
    pub column_1_step_sizes: Vec<Fr>,
    pub column_2_step_sizes: Vec<Fr>,
    pub column_3_step_sizes: Vec<Fr>,

    /// Functions for computing value from key for each basic table.
    pub get_table_values: Vec<GetValuesFromKey>,
}

impl MultiTable {
    /// Create a multi-table with repeated coefficient values.
    ///
    /// Each column's coefficients form a geometric series: 1, c, c^2, c^3, ...
    pub fn new_repeated(
        col_1_repeated_coeff: Fr,
        col_2_repeated_coeff: Fr,
        col_3_repeated_coeff: Fr,
        num_lookups: usize,
    ) -> Self {
        let mut col1 = Vec::with_capacity(num_lookups + 1);
        let mut col2 = Vec::with_capacity(num_lookups + 1);
        let mut col3 = Vec::with_capacity(num_lookups + 1);

        col1.push(Fr::one());
        col2.push(Fr::one());
        col3.push(Fr::one());

        for _ in 0..num_lookups {
            col1.push(*col1.last().unwrap() * col_1_repeated_coeff);
            col2.push(*col2.last().unwrap() * col_2_repeated_coeff);
            col3.push(*col3.last().unwrap() * col_3_repeated_coeff);
        }

        let mut table = Self {
            column_1_coefficients: col1,
            column_2_coefficients: col2,
            column_3_coefficients: col3,
            id: MultiTableId::SHA256_CH_INPUT,
            basic_table_ids: Vec::new(),
            slice_sizes: Vec::new(),
            column_1_step_sizes: Vec::new(),
            column_2_step_sizes: Vec::new(),
            column_3_step_sizes: Vec::new(),
            get_table_values: Vec::new(),
        };
        table.init_step_sizes();
        table
    }

    /// Create a multi-table with explicit per-lookup coefficient vectors.
    pub fn new_explicit(
        col_1_coeffs: Vec<Fr>,
        col_2_coeffs: Vec<Fr>,
        col_3_coeffs: Vec<Fr>,
    ) -> Self {
        let mut table = Self {
            column_1_coefficients: col_1_coeffs,
            column_2_coefficients: col_2_coeffs,
            column_3_coefficients: col_3_coeffs,
            id: MultiTableId::SHA256_CH_INPUT,
            basic_table_ids: Vec::new(),
            slice_sizes: Vec::new(),
            column_1_step_sizes: Vec::new(),
            column_2_step_sizes: Vec::new(),
            column_3_step_sizes: Vec::new(),
            get_table_values: Vec::new(),
        };
        table.init_step_sizes();
        table
    }

    /// Create an empty multi-table.
    pub fn empty() -> Self {
        Self {
            column_1_coefficients: Vec::new(),
            column_2_coefficients: Vec::new(),
            column_3_coefficients: Vec::new(),
            id: MultiTableId::SHA256_CH_INPUT,
            basic_table_ids: Vec::new(),
            slice_sizes: Vec::new(),
            column_1_step_sizes: Vec::new(),
            column_2_step_sizes: Vec::new(),
            column_3_step_sizes: Vec::new(),
            get_table_values: Vec::new(),
        }
    }

    /// Compute step sizes from coefficients via batch inversion.
    fn init_step_sizes(&mut self) {
        let num_lookups = self.column_1_coefficients.len();
        if num_lookups == 0 {
            return;
        }

        self.column_1_step_sizes = Vec::with_capacity(num_lookups);
        self.column_2_step_sizes = Vec::with_capacity(num_lookups);
        self.column_3_step_sizes = Vec::with_capacity(num_lookups);

        // First step size is always 1
        self.column_1_step_sizes.push(Fr::one());
        self.column_2_step_sizes.push(Fr::one());
        self.column_3_step_sizes.push(Fr::one());

        // Collect all coefficients for batch inversion
        let mut coefficient_inverses: Vec<Fr> = Vec::with_capacity(num_lookups * 3);
        coefficient_inverses.extend_from_slice(&self.column_1_coefficients);
        coefficient_inverses.extend_from_slice(&self.column_2_coefficients);
        coefficient_inverses.extend_from_slice(&self.column_3_coefficients);

        Fr::batch_invert(&mut coefficient_inverses);

        // step_size[i] = coefficient[i] * inverse(coefficient[i-1])
        for i in 1..num_lookups {
            self.column_1_step_sizes.push(
                self.column_1_coefficients[i] * coefficient_inverses[i - 1],
            );
            self.column_2_step_sizes.push(
                self.column_2_coefficients[i] * coefficient_inverses[num_lookups + i - 1],
            );
            self.column_3_step_sizes.push(
                self.column_3_coefficients[i] * coefficient_inverses[2 * num_lookups + i - 1],
            );
        }
    }
}

// ---------------------------------------------------------------------------
// ReadData
// ---------------------------------------------------------------------------

/// Container for lookup table read data (accumulators).
///
/// Stores the accumulated column values for a sequence of lookups, plus
/// the lookup entries themselves.
#[derive(Debug, Clone)]
pub struct ReadData<T: Clone> {
    columns: [Vec<T>; 3],
    pub lookup_entries: Vec<LookupEntry>,
}

impl<T: Clone> ReadData<T> {
    pub fn new() -> Self {
        Self {
            columns: [Vec::new(), Vec::new(), Vec::new()],
            lookup_entries: Vec::new(),
        }
    }

    pub fn column(&self, idx: ColumnIdx) -> &Vec<T> {
        &self.columns[idx as usize]
    }

    pub fn column_mut(&mut self, idx: ColumnIdx) -> &mut Vec<T> {
        &mut self.columns[idx as usize]
    }
}

impl<T: Clone> Default for ReadData<T> {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Helper: Fr from U256
// ---------------------------------------------------------------------------

/// Convert a U256 into an Fr field element (non-Montgomery limbs → Montgomery form).
pub fn fr_from_u256(val: U256) -> Fr {
    let limbs = val.as_words();
    Field::<Bn254FrParams>::from_limbs([limbs[0], limbs[1], limbs[2], limbs[3]])
}
