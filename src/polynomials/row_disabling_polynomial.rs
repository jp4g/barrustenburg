use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;

/// Polynomial for sumcheck with disabled rows.
///
/// Tracks evaluations of the Lagrange sum L = L_{n-1} + L_{n-2} + L_{n-3} + L_{n-4}
/// at 0 and 1 through sumcheck rounds. The full polynomial is `1 - L`, which
/// vanishes at the last 4 rows and equals 1 elsewhere on the hypercube.
///
/// Port of C++ `RowDisablingPolynomial<FF>`.
pub struct RowDisablingPolynomial<P: FieldParams> {
    pub eval_at_0: Field<P>,
    pub eval_at_1: Field<P>,
}

impl<P: FieldParams> Clone for RowDisablingPolynomial<P> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<P: FieldParams> Copy for RowDisablingPolynomial<P> {}

impl<P: FieldParams> Default for RowDisablingPolynomial<P> {
    fn default() -> Self {
        Self {
            eval_at_0: Field::one(),
            eval_at_1: Field::one(),
        }
    }
}

impl<P: FieldParams> RowDisablingPolynomial<P> {
    /// Update evaluations at 0 and 1 for the given round.
    ///
    /// - Round 0: no change (both stay 1).
    /// - Round 1: eval_at_0 becomes 0.
    /// - Round >= 2: eval_at_1 *= round_challenge.
    pub fn update_evaluations(&mut self, round_challenge: Field<P>, round_idx: usize) {
        if round_idx == 1 {
            self.eval_at_0 = Field::zero();
        }
        if round_idx >= 2 {
            self.eval_at_1 = self.eval_at_1 * round_challenge;
        }
    }

    /// Evaluate `1 - L` at the full sumcheck challenge vector.
    ///
    /// Returns `1 - prod(challenges[2..log_circuit_size])`.
    pub fn evaluate_at_challenge(
        challenges: &[Field<P>],
        log_circuit_size: usize,
    ) -> Field<P> {
        let mut acc = Field::one();
        for idx in 2..log_circuit_size {
            acc = acc * challenges[idx];
        }
        Field::one() - acc
    }

    /// Evaluate `1 - L` at the sumcheck challenge with padding indicators.
    ///
    /// Returns `1 - prod_{i=2..len} (1 - indicator[i] + indicator[i] * challenges[i])`.
    pub fn evaluate_at_challenge_with_padding(
        challenges: &[Field<P>],
        padding_indicators: &[Field<P>],
    ) -> Field<P> {
        let mut acc = Field::one();
        for idx in 2..padding_indicators.len() {
            let indicator = padding_indicators[idx];
            acc = acc * (Field::one() - indicator + indicator * challenges[idx]);
        }
        Field::one() - acc
    }
}
