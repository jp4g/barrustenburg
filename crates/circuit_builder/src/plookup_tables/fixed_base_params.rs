//! Fixed-base scalar multiplication parameters.
//!
//! C++ source: plookup_tables/fixed_base/fixed_base_params.hpp

/// Ceiling division (const-compatible).
const fn const_ceil_div(n: u64, d: u64) -> u64 {
    (n + d - 1) / d
}

pub struct FixedBaseParams;

impl FixedBaseParams {
    pub const BITS_PER_TABLE: u64 = 9;
    pub const BITS_ON_CURVE: u64 = 254;

    pub const BITS_PER_LO_SCALAR: u64 = 128;
    pub const BITS_PER_HI_SCALAR: u64 = Self::BITS_ON_CURVE - Self::BITS_PER_LO_SCALAR; // 126

    /// Max table size: 2^BITS_PER_TABLE = 512
    pub const MAX_TABLE_SIZE: u64 = 1u64 << Self::BITS_PER_TABLE;

    pub const NUM_FIXED_BASE_MULTI_TABLES: u64 = 4;

    /// ceil(128 / 9) = 15
    pub const NUM_TABLES_PER_LO_MULTITABLE: u64 =
        const_ceil_div(Self::BITS_PER_LO_SCALAR, Self::BITS_PER_TABLE);

    /// ceil(126 / 9) = 14
    pub const NUM_TABLES_PER_HI_MULTITABLE: u64 =
        const_ceil_div(Self::BITS_PER_HI_SCALAR, Self::BITS_PER_TABLE);

    pub const MAX_NUM_TABLES_IN_MULTITABLE: u64 = if Self::NUM_TABLES_PER_LO_MULTITABLE
        > Self::NUM_TABLES_PER_HI_MULTITABLE
    {
        Self::NUM_TABLES_PER_LO_MULTITABLE
    } else {
        Self::NUM_TABLES_PER_HI_MULTITABLE
    };

    /// Returns the number of scalar mul bits for the multitable at the given index.
    /// Even indices (0, 2) are LO_SCALAR tables, odd (1, 3) are HI_SCALAR.
    pub const fn get_num_bits_of_multi_table(multitable_index: u64) -> u64 {
        if multitable_index % 2 == 0 {
            Self::BITS_PER_LO_SCALAR
        } else {
            Self::BITS_PER_HI_SCALAR
        }
    }
}
