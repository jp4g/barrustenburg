//! Port of `circuit_type.hpp` — CircuitType enum.

/// Type of the circuit used in the proving system.
///
/// Port of C++ `enum CircuitType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CircuitType {
    Undefined,
    Standard,
    Ultra,
}
