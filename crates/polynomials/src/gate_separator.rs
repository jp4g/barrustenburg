use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use crate::eq_polynomial::ProverEqPolynomial;

/// Gate separator polynomial for sumcheck protocol.
///
/// Computes pow_beta(X) = prod_k ((1 - X_k) + X_k * beta_k) and maintains
/// incremental partial evaluation through sumcheck rounds.
///
/// Port of C++ `GateSeparatorPolynomial<FF>`.
pub struct GateSeparatorPolynomial<P: FieldParams> {
    pub betas: Vec<Field<P>>,
    pub beta_products: Vec<Field<P>>,
    pub current_element_idx: usize,
    pub periodicity: usize,
    pub partial_evaluation_result: Field<P>,
}

impl<P: FieldParams> Clone for GateSeparatorPolynomial<P> {
    fn clone(&self) -> Self {
        Self {
            betas: self.betas.clone(),
            beta_products: self.beta_products.clone(),
            current_element_idx: self.current_element_idx,
            periodicity: self.periodicity,
            partial_evaluation_result: self.partial_evaluation_result,
        }
    }
}

impl<P: FieldParams> GateSeparatorPolynomial<P> {
    /// Full prover constructor: precomputes beta_products table.
    pub fn new(betas: Vec<Field<P>>, log_num_monomials: usize) -> Self {
        let beta_products = if betas.is_empty() {
            vec![Field::zero()]
        } else {
            ProverEqPolynomial::compute_subset_products(&betas, log_num_monomials, Field::one())
        };
        Self {
            betas,
            beta_products,
            current_element_idx: 0,
            periodicity: 2,
            partial_evaluation_result: Field::one(),
        }
    }

    /// Verifier constructor: no beta_products table needed.
    pub fn new_verifier(betas: Vec<Field<P>>) -> Self {
        Self {
            betas,
            beta_products: vec![],
            current_element_idx: 0,
            periodicity: 2,
            partial_evaluation_result: Field::one(),
        }
    }

    /// Post-challenge constructor: creates verifier poly then partially evaluates
    /// at each challenge in sequence.
    pub fn new_with_challenges(betas: Vec<Field<P>>, challenges: &[Field<P>]) -> Self {
        let mut poly = Self::new_verifier(betas);
        if !poly.betas.is_empty() {
            for &u_k in challenges {
                poly.partially_evaluate(u_k);
            }
        }
        poly
    }

    /// Returns the beta_product at the appropriate index for the current round.
    #[inline]
    pub fn at(&self, idx: usize) -> Field<P> {
        self.beta_products[(idx >> 1) * self.periodicity]
    }

    /// Returns the current beta element, or 1 if betas is empty.
    #[inline]
    pub fn current_element(&self) -> Field<P> {
        if self.betas.is_empty() {
            Field::one()
        } else {
            self.betas[self.current_element_idx]
        }
    }

    /// Evaluate ((1 - X_i) + X_i * beta_i) at X_i = challenge.
    #[inline]
    pub fn univariate_eval(&self, challenge: Field<P>) -> Field<P> {
        Field::one() + challenge * (self.betas[self.current_element_idx] - Field::one())
    }

    /// Partially evaluate: multiply partial_evaluation_result by univariate_eval(challenge).
    pub fn partially_evaluate(&mut self, challenge: Field<P>) {
        if !self.betas.is_empty() {
            let eval = self.univariate_eval(challenge);
            self.partial_evaluation_result = self.partial_evaluation_result * eval;
            self.current_element_idx += 1;
            self.periodicity *= 2;
        }
    }

    /// Partially evaluate with indicator for dummy rounds.
    /// result = (1 - indicator) * result + indicator * result * eval
    pub fn partially_evaluate_with_indicator(
        &mut self,
        challenge: Field<P>,
        indicator: Field<P>,
    ) {
        if !self.betas.is_empty() {
            let eval = self.univariate_eval(challenge);
            self.partial_evaluation_result = (Field::one() - indicator)
                * self.partial_evaluation_result
                + indicator * self.partial_evaluation_result * eval;
            self.current_element_idx += 1;
            self.periodicity *= 2;
        }
    }
}
