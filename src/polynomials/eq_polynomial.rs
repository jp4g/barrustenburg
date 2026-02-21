use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;

// ---------------------------------------------------------------------------
// VerifierEqPolynomial
// ---------------------------------------------------------------------------

/// Verifier-side eq(X, r) polynomial for O(d) evaluation without divisions.
///
/// eq(r, u) = prod_i ((1 - r_i)(1 - u_i) + r_i * u_i)
///          = prod_i (b_i + u_i * a_i)
///
/// Port of C++ `VerifierEqPolynomial<FF>`.
pub struct VerifierEqPolynomial<P: FieldParams> {
    pub r: Vec<Field<P>>,
    pub a: Vec<Field<P>>, // a_i = 2*r_i - 1
    pub b: Vec<Field<P>>, // b_i = 1 - r_i
}

impl<P: FieldParams> VerifierEqPolynomial<P> {
    pub fn new(r: Vec<Field<P>>) -> Self {
        let d = r.len();
        let mut a = vec![Field::<P>::zero(); d];
        let mut b = vec![Field::<P>::zero(); d];
        let one = Field::<P>::one();
        for i in 0..d {
            a[i] = r[i] + r[i] - one; // 2*r_i - 1
            b[i] = one - r[i]; // 1 - r_i
        }
        Self { r, a, b }
    }

    /// Evaluate eq(r, u) = prod_i (b_i + u_i * a_i).
    pub fn evaluate(&self, u: &[Field<P>]) -> Field<P> {
        assert_eq!(u.len(), self.r.len());
        let mut acc = Field::<P>::one();
        for i in 0..u.len() {
            acc = acc * (self.b[i] + u[i] * self.a[i]);
        }
        acc
    }

    /// Static evaluation without constructing the struct.
    pub fn eval(r: &[Field<P>], u: &[Field<P>]) -> Field<P> {
        assert_eq!(r.len(), u.len());
        let one = Field::<P>::one();
        let mut acc = one;
        for i in 0..r.len() {
            let ai = r[i] + r[i] - one;
            let bi = one - r[i];
            acc = acc * (bi + u[i] * ai);
        }
        acc
    }
}

// ---------------------------------------------------------------------------
// ProverEqPolynomial
// ---------------------------------------------------------------------------

/// Prover-side eq(X, r) polynomial over the Boolean hypercube.
///
/// Provides `construct(challenges, log_num_monomials) -> Vec<Field>` which builds
/// the 2^d coefficient table indexed by Boolean masks.
///
/// Port of C++ `ProverEqPolynomial<FF>`.
pub struct ProverEqPolynomial;

impl ProverEqPolynomial {
    /// Construct eq(X, r) coefficient table over {0,1}^d.
    ///
    /// Auto-selects between optimal path (when no r_i == 1) and fallback.
    pub fn construct<P: FieldParams>(
        challenges: &[Field<P>],
        log_num_monomials: usize,
    ) -> Vec<Field<P>> {
        let scaling_factor = Self::compute_scaling_factor(challenges);

        if scaling_factor.is_zero() {
            return Self::construct_eq_with_edge_cases(challenges, log_num_monomials);
        }

        // Optimal path: transform to gamma_i = r_i / (1 - r_i), then use subset products
        let gammas = Self::transform_challenge(challenges);
        Self::compute_subset_products(&gammas, log_num_monomials, scaling_factor)
    }

    /// C = prod_i (1 - r_i). If zero, at least one r_i == 1.
    pub fn compute_scaling_factor<P: FieldParams>(challenges: &[Field<P>]) -> Field<P> {
        let one = Field::<P>::one();
        let mut out = one;
        for &u_i in challenges {
            out = out * (one - u_i);
        }
        out
    }

    /// Transform r_i -> gamma_i = r_i / (1 - r_i) using batch inversion.
    fn transform_challenge<P: FieldParams>(challenges: &[Field<P>]) -> Vec<Field<P>> {
        let one = Field::<P>::one();
        let mut denominators: Vec<Field<P>> = challenges.iter().map(|&c| one - c).collect();
        Field::batch_invert(&mut denominators);
        challenges
            .iter()
            .zip(denominators.iter())
            .map(|(&c, &d_inv)| d_inv * c)
            .collect()
    }

    /// Compute subset product table: result[mask] = scaling_factor * prod_{bit k set in mask} betas[k].
    ///
    /// Port of C++ `GateSeparatorPolynomial::compute_beta_products` (single-threaded).
    pub fn compute_subset_products<P: FieldParams>(
        betas: &[Field<P>],
        log_num_monomials: usize,
        scaling_factor: Field<P>,
    ) -> Vec<Field<P>> {
        let n = 1usize << log_num_monomials;
        let mut result = vec![Field::<P>::zero(); n];
        result[0] = scaling_factor;

        for (beta_idx, &beta) in betas.iter().enumerate() {
            let window = 1usize << beta_idx;
            for j in 0..window {
                result[window + j] = beta * result[j];
            }
        }

        result
    }

    /// Fallback: direct incremental construction when some r_i == 1.
    ///
    /// Builds eq(X,r) = prod_i (b_i + a_i * X_i) by expanding one variable at a time.
    fn construct_eq_with_edge_cases<P: FieldParams>(
        r: &[Field<P>],
        log_num_monomials: usize,
    ) -> Vec<Field<P>> {
        let d = r.len();
        assert_eq!(d, log_num_monomials);
        let n = 1usize << d;
        let one = Field::<P>::one();

        // Precompute per-variable coefficients
        let mut eq_linear: Vec<Field<P>> = Vec::with_capacity(d); // a_i = 2*r_i - 1
        let mut eq_constant: Vec<Field<P>> = Vec::with_capacity(d); // b_i = 1 - r_i
        for i in 0..d {
            eq_linear.push(r[i] + r[i] - one);
            eq_constant.push(one - r[i]);
        }

        let mut current = vec![Field::<P>::one()]; // start with [1]
        current.reserve(n);

        for var_idx in 0..d {
            let a_i = eq_linear[var_idx];
            let b_i = eq_constant[var_idx];

            let prev_size = current.len();
            let mut next = vec![Field::<P>::zero(); prev_size << 1];

            for j in 0..prev_size {
                let v = current[j];
                next[j] = v * b_i; // bit_i = 0
                next[j + prev_size] = v * (b_i + a_i); // bit_i = 1
            }

            current = next;
        }

        current
    }
}
