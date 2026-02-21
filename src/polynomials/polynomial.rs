use std::ops::{AddAssign, MulAssign, SubAssign};

use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;
use crate::polynomials::polynomial_arithmetic;
use crate::polynomials::polynomial_span::PolynomialSpan;
use crate::polynomials::shared_shifted_virtual_zeroes_array::SharedShiftedVirtualZeroesArray;

/// Dense univariate polynomial backed by `SharedShiftedVirtualZeroesArray`.
/// Coefficients are stored as `[c_0, c_1, ..., c_{n-1}]`.
pub struct Polynomial<P: FieldParams> {
    coefficients: SharedShiftedVirtualZeroesArray<P>,
}

// Manual Clone because derive requires P: Clone unnecessarily.
impl<P: FieldParams> Clone for Polynomial<P> {
    fn clone(&self) -> Self {
        Self {
            coefficients: self.coefficients.clone(),
        }
    }
}

// ── Constructors ──────────────────────────────────────────────────────────────

impl<P: FieldParams> Polynomial<P> {
    /// Zero polynomial of the given real `size` within a `virtual_size` starting at `start_index`.
    pub fn new(size: usize, virtual_size: usize, start_index: usize) -> Self {
        Self {
            coefficients: SharedShiftedVirtualZeroesArray::new(size, virtual_size, start_index),
        }
    }

    /// Wrap an existing coefficient vector (start_index = 0).
    pub fn from_coefficients(coeffs: Vec<Field<P>>, virtual_size: usize) -> Self {
        Self {
            coefficients: SharedShiftedVirtualZeroesArray::from_vec(coeffs, virtual_size, 0),
        }
    }

    /// Construct a polynomial from its evaluations at the given points via Lagrange interpolation.
    pub fn from_interpolation(
        points: &[Field<P>],
        evaluations: &[Field<P>],
        virtual_size: usize,
    ) -> Self {
        let coeffs = polynomial_arithmetic::compute_efficient_interpolation(evaluations, points);
        Self::from_coefficients(coeffs, virtual_size)
    }

    /// Allocate a polynomial of `virtual_size` elements (all zero) that can later be `shifted()`.
    pub fn shiftable(virtual_size: usize) -> Self {
        Self::new(virtual_size, virtual_size, 0)
    }

    /// Random polynomial with `size` non-zero coefficients.
    pub fn random(size: usize, virtual_size: usize, start_index: usize) -> Self {
        let data: Vec<Field<P>> = (0..size).map(|_| Field::random_element()).collect();
        Self {
            coefficients: SharedShiftedVirtualZeroesArray::from_vec(
                data,
                virtual_size,
                start_index,
            ),
        }
    }
}

// ── Accessors ─────────────────────────────────────────────────────────────────

impl<P: FieldParams> Polynomial<P> {
    #[inline]
    pub fn get(&self, i: usize) -> Field<P> {
        self.coefficients.get(i)
    }

    /// Mutable reference to the coefficient at logical index `i`.
    /// Panics if `i` is outside the real data range.
    #[inline]
    pub fn at_mut(&mut self, i: usize) -> &mut Field<P> {
        let start = self.start_index();
        &mut self.coefficients.data_mut()[i - start]
    }

    #[inline]
    pub fn data(&self) -> &[Field<P>] {
        self.coefficients.data()
    }

    #[inline]
    pub fn data_mut(&mut self) -> &mut [Field<P>] {
        self.coefficients.data_mut()
    }

    #[inline]
    pub fn size(&self) -> usize {
        self.coefficients.size()
    }

    #[inline]
    pub fn virtual_size(&self) -> usize {
        self.coefficients.virtual_size()
    }

    #[inline]
    pub fn start_index(&self) -> usize {
        self.coefficients.start_index()
    }

    #[inline]
    pub fn end_index(&self) -> usize {
        self.coefficients.end_index()
    }

    pub fn is_zero(&self) -> bool {
        self.data().iter().all(|c| c.is_zero())
    }

    pub fn is_empty(&self) -> bool {
        self.size() == 0
    }

    pub fn as_span(&self) -> PolynomialSpan<'_, P> {
        PolynomialSpan::new(self.data(), self.start_index())
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

impl<P: FieldParams> Polynomial<P> {
    /// Expand the real data range to cover `[new_start, new_end)`, copying existing
    /// coefficients and zero-filling the new regions.
    fn ensure_range(&mut self, new_start: usize, new_end: usize) {
        let cur_start = self.start_index();
        let cur_end = self.end_index();
        if new_start >= cur_start && new_end <= cur_end {
            return;
        }
        let actual_start = new_start.min(cur_start);
        let actual_end = new_end.max(cur_end);
        let new_size = actual_end - actual_start;
        let mut new_data = vec![Field::<P>::zero(); new_size];
        let offset = cur_start - actual_start;
        let old_data = self.data();
        for i in 0..old_data.len() {
            new_data[offset + i] = old_data[i];
        }
        self.coefficients =
            SharedShiftedVirtualZeroesArray::from_vec(new_data, self.virtual_size(), actual_start);
    }
}

// ── Arithmetic ────────────────────────────────────────────────────────────────

impl<P: FieldParams> AddAssign<PolynomialSpan<'_, P>> for Polynomial<P> {
    fn add_assign(&mut self, other: PolynomialSpan<'_, P>) {
        self.ensure_range(other.start_index, other.end_index());
        let self_start = self.start_index();
        let data = self.data_mut();
        for i in other.start_index..other.end_index() {
            data[i - self_start] = data[i - self_start] + other.get(i);
        }
    }
}

impl<P: FieldParams> SubAssign<PolynomialSpan<'_, P>> for Polynomial<P> {
    fn sub_assign(&mut self, other: PolynomialSpan<'_, P>) {
        self.ensure_range(other.start_index, other.end_index());
        let self_start = self.start_index();
        let data = self.data_mut();
        for i in other.start_index..other.end_index() {
            data[i - self_start] = data[i - self_start] - other.get(i);
        }
    }
}

impl<P: FieldParams> MulAssign<Field<P>> for Polynomial<P> {
    fn mul_assign(&mut self, scalar: Field<P>) {
        for c in self.data_mut().iter_mut() {
            *c = *c * scalar;
        }
    }
}

impl<P: FieldParams> Polynomial<P> {
    /// `self += other * scalar`
    pub fn add_scaled(&mut self, other: &PolynomialSpan<P>, scalar: Field<P>) {
        self.ensure_range(other.start_index, other.end_index());
        let self_start = self.start_index();
        let data = self.data_mut();
        for i in other.start_index..other.end_index() {
            data[i - self_start] = data[i - self_start] + other.get(i) * scalar;
        }
    }
}

// ── Evaluation ────────────────────────────────────────────────────────────────

impl<P: FieldParams> Polynomial<P> {
    /// Evaluate the polynomial at `z` using Horner's method.
    pub fn evaluate(&self, z: &Field<P>) -> Field<P> {
        polynomial_arithmetic::evaluate(self.data(), z)
    }

    /// Evaluate the multilinear extension at the given evaluation points.
    ///
    /// The polynomial is treated as having `2^n` coefficients (where `n = evaluation_points.len()`).
    /// If `shift` is true, the polynomial is conceptually shifted left by one position
    /// (coefficient i becomes coefficient i-1).
    pub fn evaluate_mle(&self, evaluation_points: &[Field<P>], shift: bool) -> Field<P> {
        let n = evaluation_points.len();
        let poly_size = 1usize << n;

        let mut tmp = vec![Field::<P>::zero(); poly_size];
        let start = self.start_index();
        if shift {
            for i in 0..poly_size {
                tmp[i] = self.get(i + 1 + start);
            }
        } else {
            for i in 0..poly_size {
                tmp[i] = self.get(i + start);
            }
        }

        let mut half = poly_size / 2;
        for j in 0..n {
            let u = evaluation_points[j];
            let one_minus_u = Field::one() - u;
            for i in 0..half {
                tmp[i] = tmp[2 * i] * one_minus_u + tmp[2 * i + 1] * u;
            }
            half /= 2;
        }
        tmp[0]
    }
}

// ── Manipulation ──────────────────────────────────────────────────────────────

impl<P: FieldParams> Polynomial<P> {
    /// Return a new polynomial whose coefficients start at index `self.start_index() + 1`,
    /// effectively dividing by X (dropping the constant term).
    pub fn shifted(&self) -> Self {
        let old_start = self.start_index();
        let old_data = self.data();
        assert!(
            old_data.len() > 1,
            "shifted() requires at least 2 backing elements"
        );
        let new_data = old_data[1..].to_vec();
        Self {
            coefficients: SharedShiftedVirtualZeroesArray::from_vec(
                new_data,
                self.virtual_size(),
                old_start + 1,
            ),
        }
    }

    /// Return a copy expanded to cover the full `[0, virtual_size)` range, with zeros
    /// where no real data exists.
    pub fn full(&self) -> Self {
        let vs = self.virtual_size();
        let mut data = vec![Field::<P>::zero(); vs];
        let start = self.start_index();
        for (i, c) in self.data().iter().enumerate() {
            data[start + i] = *c;
        }
        Self {
            coefficients: SharedShiftedVirtualZeroesArray::from_vec(data, vs, 0),
        }
    }

    /// Divide this polynomial by `(X - root)` in-place via synthetic division.
    pub fn factor_roots(&mut self, root: &Field<P>) {
        polynomial_arithmetic::factor_roots(self.data_mut(), root);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecc::curves::bn254::Bn254FrParams;

    type Fr = Field<Bn254FrParams>;

    #[test]
    fn test_from_coefficients_and_evaluate() {
        // p(X) = 1 + 2X + 3X^2, p(2) = 1 + 4 + 12 = 17
        let p =
            Polynomial::from_coefficients(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)], 4);
        assert_eq!(p.evaluate(&Fr::from(2u64)), Fr::from(17u64));
    }

    #[test]
    fn test_is_zero() {
        let p: Polynomial<Bn254FrParams> = Polynomial::new(4, 8, 0);
        assert!(p.is_zero());
    }

    #[test]
    fn test_add_assign_span() {
        let mut p =
            Polynomial::from_coefficients(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)], 4);
        let other_data = vec![Fr::from(10u64), Fr::from(20u64)];
        let span = PolynomialSpan::new(&other_data, 0);
        p += span;
        assert_eq!(p.get(0), Fr::from(11u64));
        assert_eq!(p.get(1), Fr::from(22u64));
        assert_eq!(p.get(2), Fr::from(3u64));
    }

    #[test]
    fn test_sub_assign_span() {
        let mut p = Polynomial::from_coefficients(
            vec![Fr::from(10u64), Fr::from(20u64), Fr::from(30u64)],
            4,
        );
        let other_data = vec![Fr::from(1u64), Fr::from(2u64)];
        let span = PolynomialSpan::new(&other_data, 0);
        p -= span;
        assert_eq!(p.get(0), Fr::from(9u64));
        assert_eq!(p.get(1), Fr::from(18u64));
        assert_eq!(p.get(2), Fr::from(30u64));
    }

    #[test]
    fn test_mul_assign_scalar() {
        let mut p =
            Polynomial::from_coefficients(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)], 4);
        p *= Fr::from(5u64);
        assert_eq!(p.get(0), Fr::from(5u64));
        assert_eq!(p.get(1), Fr::from(10u64));
        assert_eq!(p.get(2), Fr::from(15u64));
    }

    #[test]
    fn test_add_scaled() {
        let mut p =
            Polynomial::from_coefficients(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)], 4);
        let other_data = vec![Fr::from(10u64), Fr::from(20u64)];
        let span = PolynomialSpan::new(&other_data, 0);
        p.add_scaled(&span, Fr::from(3u64));
        assert_eq!(p.get(0), Fr::from(31u64));
        assert_eq!(p.get(1), Fr::from(62u64));
        assert_eq!(p.get(2), Fr::from(3u64));
    }

    #[test]
    fn test_shifted() {
        let p = Polynomial::from_coefficients(
            vec![
                Fr::from(1u64),
                Fr::from(2u64),
                Fr::from(3u64),
                Fr::from(4u64),
            ],
            8,
        );
        let s = p.shifted();
        assert_eq!(s.start_index(), 1);
        assert_eq!(s.size(), 3);
        assert_eq!(s.get(1), Fr::from(2u64));
        assert_eq!(s.get(2), Fr::from(3u64));
        assert_eq!(s.get(3), Fr::from(4u64));
    }

    #[test]
    fn test_full() {
        let p = Polynomial::from_coefficients(vec![Fr::from(5u64), Fr::from(6u64)], 4);
        let f = p.full();
        assert_eq!(f.start_index(), 0);
        assert_eq!(f.size(), 4);
        assert_eq!(f.get(0), Fr::from(5u64));
        assert_eq!(f.get(1), Fr::from(6u64));
        assert_eq!(f.get(2), Fr::zero());
        assert_eq!(f.get(3), Fr::zero());
    }

    #[test]
    fn test_factor_roots_and_evaluate() {
        // p(X) = (X - 2)(X - 5) = X^2 - 7X + 10
        let mut p =
            Polynomial::from_coefficients(vec![Fr::from(10u64), -Fr::from(7u64), Fr::one()], 4);
        p.factor_roots(&Fr::from(2u64));
        // Quotient is (X - 5), i.e. coeffs [-5, 1, 0]
        assert_eq!(p.get(0), -Fr::from(5u64));
        assert_eq!(p.get(1), Fr::one());
        assert_eq!(p.get(2), Fr::zero());
    }

    #[test]
    fn test_interpolation_round_trip() {
        let points = vec![
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
        ];
        // p(X) = X^2: evaluations are [1, 4, 9, 16]
        let evals = vec![
            Fr::from(1u64),
            Fr::from(4u64),
            Fr::from(9u64),
            Fr::from(16u64),
        ];
        let p = Polynomial::from_interpolation(&points, &evals, 8);
        for i in 0..points.len() {
            assert_eq!(p.evaluate(&points[i]), evals[i]);
        }
    }

    #[test]
    fn test_evaluate_mle() {
        // MLE with 2 variables (4 coefficients): f(0,0)=1, f(1,0)=2, f(0,1)=3, f(1,1)=4
        // f(u0, u1) = 1*(1-u0)*(1-u1) + 2*u0*(1-u1) + 3*(1-u0)*u1 + 4*u0*u1
        let p = Polynomial::from_coefficients(
            vec![
                Fr::from(1u64),
                Fr::from(2u64),
                Fr::from(3u64),
                Fr::from(4u64),
            ],
            4,
        );
        // Evaluate at (0, 0) -> 1
        let result = p.evaluate_mle(&[Fr::zero(), Fr::zero()], false);
        assert_eq!(result, Fr::from(1u64));
        // Evaluate at (1, 0) -> 2
        let result = p.evaluate_mle(&[Fr::one(), Fr::zero()], false);
        assert_eq!(result, Fr::from(2u64));
        // Evaluate at (1, 1) -> 4
        let result = p.evaluate_mle(&[Fr::one(), Fr::one()], false);
        assert_eq!(result, Fr::from(4u64));
    }

    #[test]
    fn test_as_span() {
        let p = Polynomial::from_coefficients(vec![Fr::from(1u64), Fr::from(2u64)], 4);
        let span = p.as_span();
        assert_eq!(span.start_index, 0);
        assert_eq!(span.size(), 2);
        assert_eq!(span.get(0), Fr::from(1u64));
        assert_eq!(span.get(1), Fr::from(2u64));
    }
}
