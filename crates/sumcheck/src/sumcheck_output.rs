//! Port of `sumcheck_output.hpp` — Output of sumcheck protocol.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

/// Output of the sumcheck protocol.
///
/// Contains the multivariate challenge vector and the claimed evaluations
/// of all prover polynomials at the challenge point.
///
/// Port of C++ `SumcheckOutput<Flavor>`.
pub struct SumcheckOutput<P: FieldParams, AllValues> {
    /// The challenge vector u = (u_0, ..., u_{d-1}).
    pub challenge: Vec<Field<P>>,
    /// Evaluations at u of the polynomials used in sumcheck.
    pub claimed_evaluations: AllValues,
    /// Whether the verifier confirmed the evaluations (verifier only).
    pub verified: bool,
}
