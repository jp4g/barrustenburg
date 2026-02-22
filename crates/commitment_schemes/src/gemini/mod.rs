//! Gemini: multilinear-to-univariate reduction protocol.
//!
//! C++ source: barretenberg/cpp/src/barretenberg/commitment_schemes/gemini/gemini.hpp
//!             barretenberg/cpp/src/barretenberg/commitment_schemes/gemini/gemini_impl.hpp
//!
//! Protocol for opening several multi-linear polynomials at the same point.
//!
//! m = number of variables
//! n = 2^m
//! u = (u_0,...,u_{m-1})
//!
//! We batch polynomials with challenge rho, define A_0 = F + G/X, then fold
//! to produce d-1 univariate polynomials. The verifier reconstructs positive
//! evaluations from negative evaluations using the folding relation.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_polynomials::polynomial::Polynomial;

// ── Opening pair / claim types ───────────────────────────────────────────────

/// Opening pair (r, v) for some witness polynomial p(X) such that p(r) = v.
#[derive(Clone)]
pub struct OpeningPair<P: FieldParams> {
    pub challenge: Field<P>,
    pub evaluation: Field<P>,
}

/// Polynomial p and an opening pair (r, v) such that p(r) = v.
pub struct ProverOpeningClaim<P: FieldParams> {
    pub polynomial: Polynomial<P>,
    pub opening_pair: OpeningPair<P>,
    /// Gemini folds must be opened at both `challenge` and `-challenge`.
    /// This flag signals Shplonk to process the claim accordingly.
    pub gemini_fold: bool,
}

// ── Free functions ───────────────────────────────────────────────────────────

/// Compute powers of challenge rho: [1, rho, rho^2, ..., rho^{num_powers-1}].
///
/// C++ source: gemini::powers_of_rho
pub fn powers_of_rho<P: FieldParams>(rho: Field<P>, num_powers: usize) -> Vec<Field<P>> {
    let mut rhos = Vec::with_capacity(num_powers);
    if num_powers == 0 {
        return rhos;
    }
    rhos.push(Field::one());
    if num_powers == 1 {
        return rhos;
    }
    rhos.push(rho);
    for j in 2..num_powers {
        rhos.push(rhos[j - 1] * rho);
    }
    rhos
}

/// Compute squares of folding challenge r: [r, r^2, r^4, ..., r^{2^{num_squares-1}}].
///
/// C++ source: gemini::powers_of_evaluation_challenge
pub fn powers_of_evaluation_challenge<P: FieldParams>(
    r: Field<P>,
    num_squares: usize,
) -> Vec<Field<P>> {
    let mut squares = Vec::with_capacity(num_squares);
    if num_squares == 0 {
        return squares;
    }
    squares.push(r);
    for j in 1..num_squares {
        squares.push(squares[j - 1].sqr());
    }
    squares
}

// ── PolynomialBatcher ────────────────────────────────────────────────────────

/// Batches multilinear polynomials for the Gemini protocol.
///
/// Given unshifted polynomials f_j and to-be-shifted polynomials g_j, computes:
///   F = sum_j rho^j * f_j
///   G = sum_j rho^{k+j} * g_j
///   A_0 = F + G/X
///
/// Then computes partially evaluated polynomials:
///   A_0+(X) = F(X) + G(X)/r
///   A_0-(X) = F(X) - G(X)/r
///
/// C++ source: GeminiProver_::PolynomialBatcher
pub struct PolynomialBatcher<'a, P: FieldParams> {
    full_batched_size: usize,
    batched_unshifted: Polynomial<P>,
    batched_to_be_shifted_by_one: Polynomial<P>,
    /// References to unshifted polynomials.
    pub unshifted: Vec<&'a Polynomial<P>>,
    /// References to polynomials to be left-shifted by 1.
    pub to_be_shifted_by_one: Vec<&'a Polynomial<P>>,
}

impl<'a, P: FieldParams> PolynomialBatcher<'a, P> {
    pub fn new(full_batched_size: usize) -> Self {
        Self {
            full_batched_size,
            batched_unshifted: Polynomial::new(full_batched_size, full_batched_size, 0),
            batched_to_be_shifted_by_one: Polynomial::shiftable(full_batched_size),
            unshifted: Vec::new(),
            to_be_shifted_by_one: Vec::new(),
        }
    }

    pub fn has_unshifted(&self) -> bool {
        !self.unshifted.is_empty()
    }

    pub fn has_to_be_shifted_by_one(&self) -> bool {
        !self.to_be_shifted_by_one.is_empty()
    }

    /// Compute batched polynomial A_0 = F + G/X.
    ///
    /// C++ source: PolynomialBatcher::compute_batched
    pub fn compute_batched(&mut self, challenge: Field<P>) -> Polynomial<P> {
        let mut running_scalar = Field::<P>::one();

        // Batch unshifted: F = sum rho^j * f_j
        if self.has_unshifted() {
            for poly in &self.unshifted {
                let span = poly.as_span();
                self.batched_unshifted.add_scaled(&span, running_scalar);
                running_scalar = running_scalar * challenge;
            }
        }

        // Batch to-be-shifted: G = sum rho^{k+j} * g_j
        if self.has_to_be_shifted_by_one() {
            for poly in &self.to_be_shifted_by_one {
                let span = poly.as_span();
                self.batched_to_be_shifted_by_one
                    .add_scaled(&span, running_scalar);
                running_scalar = running_scalar * challenge;
            }
        }

        // A_0 = F + G/X
        let mut full_batched =
            Polynomial::<P>::new(self.full_batched_size, self.full_batched_size, 0);

        if self.has_unshifted() {
            full_batched += self.batched_unshifted.as_span();
        }

        if self.has_to_be_shifted_by_one() {
            // Construct G/X: the left-shift of the batched to-be-shifted polynomial.
            // (G/X)[i] = G[i+1] for i = 0, ..., n-2.
            //
            // NOTE: We construct this manually rather than using Polynomial::shifted()
            // because the Rust shifted() implementation increments start_index (right-shift
            // semantics), while the C++ implementation decrements it (left-shift / divide-by-X).
            // The C++ relies on shiftable polynomials having start_=1, which the Rust
            // Polynomial::shiftable() does not match.
            let g = &self.batched_to_be_shifted_by_one;
            let n = self.full_batched_size;
            let g_over_x_data: Vec<Field<P>> = (0..n).map(|i| g.get(i + 1)).collect();
            let g_over_x = Polynomial::from_coefficients(g_over_x_data, n);
            full_batched += g_over_x.as_span();
        }

        full_batched
    }

    /// Compute partially evaluated batch polynomials:
    ///   A_0+(X) = F(X) + G(X)/r
    ///   A_0-(X) = F(X) - G(X)/r
    ///
    /// C++ source: PolynomialBatcher::compute_partially_evaluated_batch_polynomials
    pub fn compute_partially_evaluated_batch_polynomials(
        &mut self,
        r_challenge: Field<P>,
    ) -> (Polynomial<P>, Polynomial<P>) {
        let mut a_0_pos =
            Polynomial::<P>::new(self.full_batched_size, self.full_batched_size, 0);

        if self.has_unshifted() {
            a_0_pos += self.batched_unshifted.as_span();
        }

        let mut a_0_neg = a_0_pos.clone();

        if self.has_to_be_shifted_by_one() {
            let r_inv = r_challenge.invert();
            self.batched_to_be_shifted_by_one *= r_inv; // G = G/r

            a_0_pos += self.batched_to_be_shifted_by_one.as_span(); // A_0+ += G/r
            a_0_neg -= self.batched_to_be_shifted_by_one.as_span(); // A_0- -= G/r
        }

        (a_0_pos, a_0_neg)
    }
}

// ── GeminiProver ─────────────────────────────────────────────────────────────

/// Gemini prover: computes fold polynomials and opening claims.
///
/// C++ source: GeminiProver_
pub struct GeminiProver;

impl GeminiProver {
    /// Compute the d-1 fold polynomials A_1, ..., A_{d-1} from A_0.
    ///
    /// Each fold is: A_{l+1}(X) = (1-u_l)*even(A_l)(X) + u_l*odd(A_l)(X)
    ///
    /// C++ source: GeminiProver_::compute_fold_polynomials
    pub fn compute_fold_polynomials<P: FieldParams>(
        log_n: usize,
        multilinear_challenge: &[Field<P>],
        a_0: &Polynomial<P>,
    ) -> Vec<Polynomial<P>> {
        let virtual_log_n = multilinear_challenge.len();

        // Allocate fold polynomials for the real rounds
        let mut fold_polynomials: Vec<Polynomial<P>> = Vec::with_capacity(virtual_log_n - 1);
        for l in 0..log_n - 1 {
            let n_l = 1 << (log_n - l - 1);
            fold_polynomials.push(Polynomial::new(n_l, n_l, 0));
        }

        // Iteratively fold
        // a_l_data points to the current polynomial's coefficients
        let mut prev_data: Vec<Field<P>> = a_0.data().to_vec();

        for l in 0..log_n - 1 {
            let n_l = 1 << (log_n - l - 1);
            // When virtual_log_n < log_n (high-degree attack scenario), challenges beyond
            // virtual_log_n are not meaningful. Use zero to avoid out-of-bounds access.
            // These fold polynomials are never committed or opened.
            let u_l = if l < virtual_log_n {
                multilinear_challenge[l]
            } else {
                Field::<P>::zero()
            };

            // fold(A_l)[j] = (1-u_l)*A_l[2j] + u_l*A_l[2j+1]
            //              = A_l[2j] + u_l*(A_l[2j+1] - A_l[2j])
            {
                let fold_data = fold_polynomials[l].data_mut();
                for j in 0..n_l {
                    let even = prev_data[j << 1];
                    let odd = prev_data[(j << 1) + 1];
                    fold_data[j] = even + u_l * (odd - even);
                }
            }
            prev_data = fold_polynomials[l].data().to_vec();
        }

        // Virtual rounds: after log_n - 1 real folds, the polynomial stabilizes.
        // Compute the final evaluation.
        let last = &fold_polynomials[log_n - 2];
        let u_last = if log_n - 1 < virtual_log_n {
            multilinear_challenge[log_n - 1]
        } else {
            Field::<P>::zero()
        };
        let final_eval = last.get(0) + u_last * (last.get(1) - last.get(0));

        // FOLD_{log_n} is a constant polynomial
        let mut const_fold = Polynomial::<P>::new(1, 1, 0);
        *const_fold.at_mut(0) = final_eval;
        fold_polynomials.push(const_fold);

        // FOLD_{log_n+1}, ..., FOLD_{virtual_log_n-1}: constant folds for padding
        let mut tail = Field::<P>::one();
        for k in log_n..virtual_log_n - 1 {
            tail = tail * (Field::one() - multilinear_challenge[k]);
            let mut next_const = Polynomial::<P>::new(1, 1, 0);
            *next_const.at_mut(0) = final_eval * tail;
            fold_polynomials.push(next_const);
        }

        fold_polynomials
    }

    /// Construct d+1 univariate opening claims from partially evaluated polynomials
    /// and fold polynomials.
    ///
    /// Claims: {A_0+(X), (r, A_0+(r))}, {A_0-(X), (-r, A_0-(-r))},
    ///         {A_l(X), (-r^{2^l}, A_l(-r^{2^l}))}, l = 1, ..., d-1
    ///
    /// C++ source: GeminiProver_::construct_univariate_opening_claims
    pub fn construct_univariate_opening_claims<P: FieldParams>(
        log_n: usize,
        a_0_pos: Polynomial<P>,
        a_0_neg: Polynomial<P>,
        fold_polynomials: Vec<Polynomial<P>>,
        r_challenge: Field<P>,
    ) -> Vec<ProverOpeningClaim<P>> {
        let mut claims = Vec::new();

        // A_0+(r)
        let a_0_pos_eval = a_0_pos.evaluate(&r_challenge);
        claims.push(ProverOpeningClaim {
            polynomial: a_0_pos,
            opening_pair: OpeningPair {
                challenge: r_challenge,
                evaluation: a_0_pos_eval,
            },
            gemini_fold: false,
        });

        // A_0-(-r)
        let neg_r = -r_challenge;
        let a_0_neg_eval = a_0_neg.evaluate(&neg_r);
        claims.push(ProverOpeningClaim {
            polynomial: a_0_neg,
            opening_pair: OpeningPair {
                challenge: neg_r,
                evaluation: a_0_neg_eval,
            },
            gemini_fold: false,
        });

        // Compute r^{2^l} for l = 0, ..., log_n-1
        let r_squares = powers_of_evaluation_challenge(r_challenge, log_n);

        // Claims for fold polynomials: {A_l, (-r^{2^l}, A_l(-r^{2^l}))} for l = 1, ..., d-1
        for (l, fold_poly) in fold_polynomials.into_iter().enumerate().take(log_n - 1) {
            let neg_r_sq = -r_squares[l + 1];
            let evaluation = fold_poly.evaluate(&neg_r_sq);
            claims.push(ProverOpeningClaim {
                polynomial: fold_poly,
                opening_pair: OpeningPair {
                    challenge: neg_r_sq,
                    evaluation,
                },
                gemini_fold: true,
            });
        }

        claims
    }

    /// Full Gemini prove: batch, fold, partially evaluate, construct claims.
    ///
    /// C++ source: GeminiProver_::prove (simplified: no transcript/commitments)
    pub fn prove<P: FieldParams>(
        log_n: usize,
        batcher: &mut PolynomialBatcher<P>,
        multilinear_challenge: &[Field<P>],
        rho: Field<P>,
        r_challenge: Field<P>,
    ) -> Vec<ProverOpeningClaim<P>> {
        let virtual_log_n = multilinear_challenge.len();

        // Compute A_0 = F + G/X
        let a_0 = batcher.compute_batched(rho);

        // Compute fold polynomials
        let fold_polynomials =
            Self::compute_fold_polynomials(log_n, multilinear_challenge, &a_0);

        // Compute A_0+(X) = F + G/r and A_0-(X) = F - G/r
        let (a_0_pos, a_0_neg) =
            batcher.compute_partially_evaluated_batch_polynomials(r_challenge);

        // Construct opening claims
        Self::construct_univariate_opening_claims(
            virtual_log_n,
            a_0_pos,
            a_0_neg,
            fold_polynomials,
            r_challenge,
        )
    }
}

// ── GeminiVerifier ───────────────────────────────────────────────────────────

/// Gemini verifier: reconstructs positive evaluations from negative evaluations.
///
/// C++ source: GeminiVerifier_
pub struct GeminiVerifier;

impl GeminiVerifier {
    /// Compute A_0(r), A_1(r^2), ..., A_{d-1}(r^{2^{d-1}}) from negative evaluations.
    ///
    /// The verifier recovers A_{l-1}(r^{2^{l-1}}) from:
    ///   - A_{l-1}(-r^{2^{l-1}}) (received from prover)
    ///   - A_l(r^{2^l}) (computed in previous step)
    ///
    /// Using the relation:
    ///   A_{l-1}(r^{2^{l-1}}) =
    ///     (2 * r^{2^{l-1}} * A_l(r^{2^l}) - A_{l-1}(-r^{2^{l-1}}) * (r^{2^{l-1}}*(1-u_{l-1}) - u_{l-1}))
    ///     / (r^{2^{l-1}} * (1-u_{l-1}) + u_{l-1})
    ///
    /// C++ source: GeminiVerifier_::compute_fold_pos_evaluations
    pub fn compute_fold_pos_evaluations<P: FieldParams>(
        padding_indicator_array: &[Field<P>],
        batched_evaluation: Field<P>,
        evaluation_point: &[Field<P>],  // size = virtual_log_n
        challenge_powers: &[Field<P>],  // r, r^2, r^4, ..., size = virtual_log_n
        fold_neg_evals: &[Field<P>],    // A_i(-r^{2^i}), size = virtual_log_n
        p_neg: Field<P>,
    ) -> Vec<Field<P>> {
        let virtual_log_n = evaluation_point.len();

        let mut evals: Vec<Field<P>> = fold_neg_evals.to_vec();

        let mut eval_pos_prev = batched_evaluation;

        // Add contribution of P-(-r^s) to get A_0(-r)
        evals[0] = evals[0] + p_neg;

        let mut fold_pos_evaluations = Vec::with_capacity(virtual_log_n);

        // Solve the sequence of linear equations (iterate from l = virtual_log_n down to 1)
        for l in (1..=virtual_log_n).rev() {
            let challenge_power = challenge_powers[l - 1]; // r^{2^{l-1}}
            let u = evaluation_point[l - 1];               // u_{l-1}
            let eval_neg = evals[l - 1];                   // A_{l-1}(-r^{2^{l-1}})

            // Numerator: 2 * r^{2^{l-1}} * A_l(r^{2^l}) - A_{l-1}(-r^{2^{l-1}}) * (r^{2^{l-1}}*(1-u) - u)
            let one_minus_u = Field::<P>::one() - u;
            let eval_pos =
                (challenge_power * eval_pos_prev * Field::<P>::from(2u64)
                    - eval_neg * (challenge_power * one_minus_u - u))
                    * (challenge_power * one_minus_u + u).invert();

            // If padding indicator is 1 (real round), use computed eval_pos.
            // If 0 (padding round), propagate batched_evaluation.
            let indicator = padding_indicator_array[l - 1];
            let one_minus_indicator = Field::<P>::one() - indicator;
            eval_pos_prev = indicator * eval_pos + one_minus_indicator * eval_pos_prev;

            // Store: real rounds get eval_pos_prev, padding rounds get 0
            fold_pos_evaluations.push(indicator * eval_pos_prev);
        }

        fold_pos_evaluations.reverse();
        fold_pos_evaluations
    }
}

#[cfg(test)]
mod tests;
