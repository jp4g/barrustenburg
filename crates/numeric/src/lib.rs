// Numeric types and operations.
//
// Mirrors barretenberg/cpp/src/barretenberg/numeric/.
//
// - uint128: Native Rust u128 (BB only has a custom impl for 32-bit/i386)
// - uint256: Backed by crypto-bigint U256
// - uintx (uint512, uint1024): Backed by crypto-bigint concat types
// - bitop: Bit manipulation utilities
// - random: RNG wrappers

pub mod bitop;
pub mod random;
pub mod uint256;
pub mod uintx;

// Re-export primary types at the numeric level, matching BB's `numeric::uint256_t` etc.
pub use uint256::{U256, U256Ext};
pub use uintx::{U512, U1024};
