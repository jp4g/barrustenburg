//! Circuit builder crate for the barretenberg Rust port.
//!
//! Provides gate data structures, execution trace blocks, and the circuit builder
//! base for the Ultra circuit builder, corresponding to the C++ `stdlib_circuit_builders`
//! module.

pub mod builder_base;
pub mod execution_trace;
pub mod gate_data;
pub mod ultra_builder;
