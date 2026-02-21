use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;

// ---------------------------------------------------------------------------
// Univariate<P, N>
// ---------------------------------------------------------------------------

/// Fixed-size polynomial represented by its evaluations on the domain {0, 1, ..., N-1}.
///
/// Port of C++ `Univariate<Fr, domain_end>` from Barretenberg's `univariate.hpp`.
pub struct Univariate<P: FieldParams, const N: usize> {
    pub evaluations: [Field<P>; N],
}

// Manual Clone/Copy because derived bounds require P: Clone/Copy unnecessarily.
impl<P: FieldParams, const N: usize> Clone for Univariate<P, N> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            evaluations: self.evaluations,
        }
    }
}

impl<P: FieldParams, const N: usize> std::fmt::Debug for Univariate<P, N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Univariate")
            .field("evaluations", &self.evaluations)
            .finish()
    }
}

impl<P: FieldParams, const N: usize> Univariate<P, N> {
    /// Construct from an array of evaluations.
    #[inline]
    pub fn new(evaluations: [Field<P>; N]) -> Self {
        Self { evaluations }
    }

    /// All evaluations set to the same scalar value.
    #[inline]
    pub fn from_scalar(value: Field<P>) -> Self {
        Self {
            evaluations: [value; N],
        }
    }

    /// All evaluations are zero.
    #[inline]
    pub fn zero() -> Self {
        Self {
            evaluations: [Field::zero(); N],
        }
    }

    /// Random evaluations (for testing).
    pub fn random() -> Self {
        let mut evals = [Field::zero(); N];
        for e in evals.iter_mut() {
            *e = Field::random_element();
        }
        Self { evaluations: evals }
    }

    /// Return the evaluation at index `i`.
    #[inline]
    pub fn value_at(&self, i: usize) -> Field<P> {
        self.evaluations[i]
    }

    /// True if every evaluation is zero.
    pub fn is_zero(&self) -> bool {
        self.evaluations.iter().all(|e| e.is_zero())
    }

    /// Domain size (number of evaluation points).
    #[inline]
    pub fn size() -> usize {
        N
    }

    // -- Scalar arithmetic helpers ------------------------------------------

    /// Pointwise addition with a scalar (applied to every evaluation).
    pub fn add_scalar(&self, s: Field<P>) -> Self {
        let mut out = *self;
        for e in out.evaluations.iter_mut() {
            *e = *e + s;
        }
        out
    }

    /// Pointwise subtraction of a scalar (applied to every evaluation).
    pub fn sub_scalar(&self, s: Field<P>) -> Self {
        let mut out = *self;
        for e in out.evaluations.iter_mut() {
            *e = *e - s;
        }
        out
    }

    /// Pointwise multiplication by a scalar (applied to every evaluation).
    pub fn mul_scalar(&self, s: Field<P>) -> Self {
        let mut out = *self;
        for e in out.evaluations.iter_mut() {
            *e = *e * s;
        }
        out
    }

    /// In-place pointwise squaring.
    pub fn self_sqr(&mut self) {
        for e in self.evaluations.iter_mut() {
            *e = *e * *e;
        }
    }

    /// Returns a new Univariate with pointwise squared evaluations.
    pub fn sqr(&self) -> Self {
        let mut out = *self;
        out.self_sqr();
        out
    }

    // -- In-place extension ---------------------------------------------------

    /// Extend evaluations in-place from {0..INITIAL-1} to {0..N-1} using Newton
    /// forward difference extrapolation. The first `INITIAL` evaluations must
    /// already be populated; the rest are overwritten.
    ///
    /// Port of C++ `Univariate::self_extend_from<INITIAL>()`.
    pub fn self_extend_from<const INITIAL: usize>(&mut self) {
        assert!(INITIAL <= N, "self_extend_from: INITIAL must be <= N");
        if INITIAL >= N {
            return;
        }

        // Build forward-difference coefficients from evaluations[0..INITIAL].
        let mut diffs = [Field::<P>::zero(); N];
        for i in 0..INITIAL {
            diffs[i] = self.evaluations[i];
        }

        // Convert evaluations to forward-difference coefficients.
        for i in 1..INITIAL {
            for j in (i..INITIAL).rev() {
                diffs[j] = diffs[j] - diffs[j - 1];
            }
        }

        // Iteratively shift the origin and record new evaluations.
        for x in 1..N {
            for j in 0..INITIAL - 1 {
                diffs[j] = diffs[j] + diffs[j + 1];
            }
            if x >= INITIAL {
                self.evaluations[x] = diffs[0];
            }
        }
    }

    // -- Barycentric evaluation ---------------------------------------------

    /// Evaluate the unique degree-(N-1) polynomial through the N evaluations at
    /// an arbitrary field point `u`, using the barycentric formula on the domain
    /// {0, 1, ..., N-1}.
    ///
    /// L_k(u) = w_k / (u - k)  /  sum_j( w_j / (u - j) )
    /// where w_k = (-1)^(N-1-k) / (k! * (N-1-k)!)
    pub fn evaluate(&self, u: Field<P>) -> Field<P> {
        if N == 0 {
            return Field::zero();
        }
        if N == 1 {
            return self.evaluations[0];
        }

        // Precompute barycentric weights: w_k = (-1)^(N-1-k) / (k! * (N-1-k)!)
        let mut weights = vec![Field::<P>::zero(); N];
        for k in 0..N {
            // k! * (N-1-k)!
            let mut denom = Field::<P>::one();
            for j in 1..=k {
                denom = denom * Field::from(j as u64);
            }
            for j in 1..=(N - 1 - k) {
                denom = denom * Field::from(j as u64);
            }
            let sign = if (N - 1 - k) % 2 == 1 {
                -Field::<P>::one()
            } else {
                Field::<P>::one()
            };
            weights[k] = sign * denom.invert();
        }

        // Check if u is one of the domain points (avoid division by zero).
        for k in 0..N {
            let domain_pt = Field::<P>::from(k as u64);
            if u == domain_pt {
                return self.evaluations[k];
            }
        }

        // Barycentric sum: f(u) = (sum_k w_k * f(k) / (u - k)) / (sum_k w_k / (u - k))
        let mut numer = Field::<P>::zero();
        let mut denom = Field::<P>::zero();
        for k in 0..N {
            let diff_inv = (u - Field::from(k as u64)).invert();
            let term = weights[k] * diff_inv;
            numer = numer + term * self.evaluations[k];
            denom = denom + term;
        }
        numer * denom.invert()
    }

    // -- Extension via Newton forward differences ---------------------------

    /// Extend evaluations from domain {0..N-1} to {0..M-1} using Newton forward
    /// difference extrapolation. `M` must be >= `N`.
    pub fn extend_to<const M: usize>(&self) -> Univariate<P, M> {
        assert!(M >= N, "extend_to: target size M must be >= N");
        let mut result = Univariate::<P, M>::zero();

        // Copy known evaluations.
        for i in 0..N {
            result.evaluations[i] = self.evaluations[i];
        }

        if M == N {
            return result;
        }

        if N == 2 {
            // Linear extrapolation: f(x) = f(0) + x * delta
            let delta = self.evaluations[1] - self.evaluations[0];
            for i in 2..M {
                result.evaluations[i] = result.evaluations[i - 1] + delta;
            }
            return result;
        }

        // General case: Newton forward differences.
        //
        // 1. Build forward difference table in-place from the N known evaluations.
        //    After this, diffs[k] = Delta^k f(0).
        // 2. Shift the origin from 0 upward one step at a time:
        //       diffs[j] += diffs[j+1]  for j = 0 .. N-2
        //    After s shifts, diffs[0] = f(s).
        // 3. We need to shift (N-1) times to reach origin (N-1), then each
        //    further shift produces the next extrapolated evaluation.
        //    We combine the initial shifts and extrapolation into a single loop
        //    over x = 1 .. M-1, recording diffs[0] = f(x) for x >= N.
        let mut diffs = vec![Field::<P>::zero(); N];
        for i in 0..N {
            diffs[i] = self.evaluations[i];
        }

        // Convert evaluations to forward-difference coefficients.
        for i in 1..N {
            for j in (i..N).rev() {
                diffs[j] = diffs[j] - diffs[j - 1];
            }
        }

        // Iteratively shift the origin and record new evaluations.
        for x in 1..M {
            for j in 0..N - 1 {
                diffs[j] = diffs[j] + diffs[j + 1];
            }
            if x >= N {
                result.evaluations[x] = diffs[0];
            }
        }

        result
    }
}

// Manual Copy (Field<P> is Copy, arrays of Copy are Copy).
impl<P: FieldParams, const N: usize> Copy for Univariate<P, N> {}

// ---------------------------------------------------------------------------
// PartialEq
// ---------------------------------------------------------------------------

impl<P: FieldParams, const N: usize> PartialEq for Univariate<P, N> {
    fn eq(&self, other: &Self) -> bool {
        self.evaluations
            .iter()
            .zip(other.evaluations.iter())
            .all(|(a, b)| a == b)
    }
}

impl<P: FieldParams, const N: usize> Eq for Univariate<P, N> {}

// ---------------------------------------------------------------------------
// Pointwise arithmetic trait impls: Univariate op Univariate
// ---------------------------------------------------------------------------

impl<P: FieldParams, const N: usize> Add for Univariate<P, N> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        let mut out = self;
        for i in 0..N {
            out.evaluations[i] = out.evaluations[i] + rhs.evaluations[i];
        }
        out
    }
}

impl<P: FieldParams, const N: usize> AddAssign for Univariate<P, N> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.evaluations[i] = self.evaluations[i] + rhs.evaluations[i];
        }
    }
}

impl<P: FieldParams, const N: usize> Sub for Univariate<P, N> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let mut out = self;
        for i in 0..N {
            out.evaluations[i] = out.evaluations[i] - rhs.evaluations[i];
        }
        out
    }
}

impl<P: FieldParams, const N: usize> SubAssign for Univariate<P, N> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.evaluations[i] = self.evaluations[i] - rhs.evaluations[i];
        }
    }
}

impl<P: FieldParams, const N: usize> Mul for Univariate<P, N> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let mut out = self;
        for i in 0..N {
            out.evaluations[i] = out.evaluations[i] * rhs.evaluations[i];
        }
        out
    }
}

impl<P: FieldParams, const N: usize> MulAssign for Univariate<P, N> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.evaluations[i] = self.evaluations[i] * rhs.evaluations[i];
        }
    }
}

impl<P: FieldParams, const N: usize> Neg for Univariate<P, N> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        let mut out = self;
        for e in out.evaluations.iter_mut() {
            *e = -*e;
        }
        out
    }
}

// ---------------------------------------------------------------------------
// UnivariateView
// ---------------------------------------------------------------------------

/// Non-owning view into univariate evaluations. Port of C++ `UnivariateView`.
pub struct UnivariateView<'a, P: FieldParams, const N: usize> {
    pub evaluations: &'a [Field<P>; N],
}

impl<'a, P: FieldParams, const N: usize> Clone for UnivariateView<'a, P, N> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, P: FieldParams, const N: usize> Copy for UnivariateView<'a, P, N> {}

impl<'a, P: FieldParams, const N: usize> UnivariateView<'a, P, N> {
    #[inline]
    pub fn new(evaluations: &'a [Field<P>; N]) -> Self {
        Self { evaluations }
    }

    #[inline]
    pub fn value_at(&self, i: usize) -> Field<P> {
        self.evaluations[i]
    }
}

impl<'a, P: FieldParams, const N: usize> From<&'a Univariate<P, N>>
    for UnivariateView<'a, P, N>
{
    #[inline]
    fn from(u: &'a Univariate<P, N>) -> Self {
        Self {
            evaluations: &u.evaluations,
        }
    }
}

// UnivariateView + Univariate -> Univariate (pointwise)
impl<'a, P: FieldParams, const N: usize> Add<Univariate<P, N>>
    for UnivariateView<'a, P, N>
{
    type Output = Univariate<P, N>;
    #[inline]
    fn add(self, rhs: Univariate<P, N>) -> Univariate<P, N> {
        let mut out = rhs;
        for i in 0..N {
            out.evaluations[i] = self.evaluations[i] + out.evaluations[i];
        }
        out
    }
}

// UnivariateView - Univariate -> Univariate (pointwise)
impl<'a, P: FieldParams, const N: usize> Sub<Univariate<P, N>>
    for UnivariateView<'a, P, N>
{
    type Output = Univariate<P, N>;
    #[inline]
    fn sub(self, rhs: Univariate<P, N>) -> Univariate<P, N> {
        let mut out = Univariate::zero();
        for i in 0..N {
            out.evaluations[i] = self.evaluations[i] - rhs.evaluations[i];
        }
        out
    }
}

// UnivariateView * Univariate -> Univariate (pointwise)
impl<'a, P: FieldParams, const N: usize> Mul<Univariate<P, N>>
    for UnivariateView<'a, P, N>
{
    type Output = Univariate<P, N>;
    #[inline]
    fn mul(self, rhs: Univariate<P, N>) -> Univariate<P, N> {
        let mut out = rhs;
        for i in 0..N {
            out.evaluations[i] = self.evaluations[i] * out.evaluations[i];
        }
        out
    }
}
