//! Port of barretenberg's Honk library support.
//!
//! Provides utility functions for the Honk proving system:
//! - Grand product delta computation (public input correction)
//! - Grand product polynomial computation
//! - Log-derivative inverse polynomial computation
//! - Relation checker (debugging utility)
//! - Proof and circuit types

pub mod circuit_type;
pub mod grand_product_delta;
pub mod grand_product_library;
pub mod logderivative_library;
pub mod proof;
pub mod relation_checker;

// Re-exports
pub use circuit_type::CircuitType;
pub use grand_product_delta::{compute_public_input_delta, PERMUTATION_ARGUMENT_VALUE_SEPARATOR};
pub use grand_product_library::{compute_grand_product, GrandProductRelation};
pub use logderivative_library::{compute_logderivative_inverse, LogDerivRelation};
pub use proof::HonkProof;
