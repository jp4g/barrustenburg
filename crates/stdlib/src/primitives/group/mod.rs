pub mod cycle_group;
pub mod cycle_scalar;
pub mod straus_lookup_table;
pub mod straus_scalar_slice;

pub use cycle_group::CycleGroupT;
pub use cycle_scalar::CycleScalarT;
pub use straus_lookup_table::StrausLookupTable;
pub use straus_scalar_slice::StrausScalarSlices;

#[cfg(test)]
mod tests;
