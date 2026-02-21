//! Port of `proof.hpp` — HonkProof type.

use bbrs_ecc::fields::field::Field;

/// A Honk proof is a vector of field elements.
///
/// Port of C++ `using HonkProof = std::vector<fr>;`
pub type HonkProof<P> = Vec<Field<P>>;
