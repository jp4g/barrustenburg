use std::ops::{Add, AddAssign, Neg, Sub, SubAssign};

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use crate::univariate::Univariate;

/// Univariate polynomial in coefficient form: P(X) = a0 + a1*X (N=2) or a0 + a1*X + a2*X^2 (N=3).
///
/// Storage layout:
/// - N=2: `coefficients = [a0, a1, a0+a1]` (third slot caches sum for Karatsuba)
/// - N=3: `coefficients = [a0, a1+a2, a2]` (forward-difference layout for efficient evaluation)
///
/// Port of C++ `UnivariateCoefficientBasis<Fr, domain_end, has_a0_plus_a1>`.
/// Simplified: always stores the cached a0+a1 for N=2 (no bool template param).
pub struct UnivariateCoefficientBasis<P: FieldParams, const N: usize> {
    pub coefficients: [Field<P>; 3],
}

impl<P: FieldParams, const N: usize> Clone for UnivariateCoefficientBasis<P, N> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<P: FieldParams, const N: usize> Copy for UnivariateCoefficientBasis<P, N> {}

impl<P: FieldParams, const N: usize> std::fmt::Debug for UnivariateCoefficientBasis<P, N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("UnivariateCoefficientBasis")
            .field("coefficients", &self.coefficients)
            .finish()
    }
}

impl<P: FieldParams, const N: usize> PartialEq for UnivariateCoefficientBasis<P, N> {
    fn eq(&self, other: &Self) -> bool {
        if N == 2 {
            self.coefficients[0] == other.coefficients[0]
                && self.coefficients[1] == other.coefficients[1]
        } else {
            self.coefficients[0] == other.coefficients[0]
                && self.coefficients[1] == other.coefficients[1]
                && self.coefficients[2] == other.coefficients[2]
        }
    }
}

impl<P: FieldParams, const N: usize> Eq for UnivariateCoefficientBasis<P, N> {}

impl<P: FieldParams, const N: usize> UnivariateCoefficientBasis<P, N> {
    /// All coefficients zero.
    pub fn zero() -> Self {
        assert!(N == 2 || N == 3);
        Self {
            coefficients: [Field::zero(); 3],
        }
    }

    /// Random coefficients (for testing).
    pub fn random() -> Self {
        assert!(N == 2 || N == 3);
        let mut coeffs = [Field::zero(); 3];
        for i in 0..N {
            coeffs[i] = Field::random_element();
        }
        if N == 2 {
            coeffs[2] = coeffs[0] + coeffs[1];
        }
        Self { coefficients: coeffs }
    }

    // -- Serialization --------------------------------------------------------

    /// Serialize coefficients to a byte buffer (3 × 32 = 96 bytes, no length prefix).
    pub fn to_buffer(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(3 * 32);
        for coeff in &self.coefficients {
            buf.extend_from_slice(&coeff.to_be_bytes());
        }
        buf
    }

    /// Deserialize coefficients from a byte buffer (3 × 32 = 96 bytes).
    pub fn serialize_from_buffer(buffer: &[u8]) -> Self {
        let mut coefficients = [Field::zero(); 3];
        for i in 0..3 {
            let bytes: [u8; 32] = buffer[i * 32..(i + 1) * 32].try_into().unwrap();
            coefficients[i] = Field::from_be_bytes(&bytes);
        }
        Self { coefficients }
    }

    /// Construct from a degree-1 Univariate (evaluation form → coefficient form).
    ///
    /// Given f(0)=v0, f(1)=v1: a0=v0, a1=v1-v0, cache=v1.
    pub fn from_univariate(u: &Univariate<P, 2>) -> UnivariateCoefficientBasis<P, 2> {
        UnivariateCoefficientBasis::<P, 2> {
            coefficients: [
                u.evaluations[0],
                u.evaluations[1] - u.evaluations[0],
                u.evaluations[1], // a0 + a1 = v0 + (v1 - v0) = v1
            ],
        }
    }
}

// ---------------------------------------------------------------------------
// Multiplication: only for N=2, outputs N=3
// ---------------------------------------------------------------------------

impl<P: FieldParams> UnivariateCoefficientBasis<P, 2> {
    /// Karatsuba multiplication: (a0 + a1*X) * (b0 + b1*X) → degree-2.
    ///
    /// Always computes (a0+a1) inline rather than relying on the cached
    /// coefficients[2], since arithmetic ops (add/sub/neg) invalidate the cache.
    pub fn mul(&self, other: &UnivariateCoefficientBasis<P, 2>) -> UnivariateCoefficientBasis<P, 3> {
        let mut result = UnivariateCoefficientBasis::<P, 3>::zero();
        // result[0] = a0 * b0
        result.coefficients[0] = self.coefficients[0] * other.coefficients[0];
        // result[2] = a1 * b1
        result.coefficients[2] = self.coefficients[1] * other.coefficients[1];
        // result[1] = (a0+a1)*(b0+b1) - a0*b0 = a0*b1 + a1*b0 + a1*b1
        result.coefficients[1] =
            (self.coefficients[0] + self.coefficients[1])
                * (other.coefficients[0] + other.coefficients[1])
                - result.coefficients[0];
        result
    }

    /// Squaring: (a0 + a1*X)^2 → degree-2.
    ///
    /// Always computes (2*a0+a1) inline rather than relying on the cached
    /// coefficients[2], since arithmetic ops (add/sub/neg) invalidate the cache.
    pub fn sqr(&self) -> UnivariateCoefficientBasis<P, 3> {
        let mut result = UnivariateCoefficientBasis::<P, 3>::zero();
        // result[0] = a0^2
        result.coefficients[0] = self.coefficients[0] * self.coefficients[0];
        // result[2] = a1^2
        result.coefficients[2] = self.coefficients[1] * self.coefficients[1];
        // result[1] = 2*a0*a1 + a1^2 = (2*a0 + a1) * a1
        result.coefficients[1] =
            (self.coefficients[0] + self.coefficients[0] + self.coefficients[1])
                * self.coefficients[1];
        result
    }
}

// ---------------------------------------------------------------------------
// Arithmetic: UnivariateCoefficientBasis op UnivariateCoefficientBasis
// ---------------------------------------------------------------------------

impl<P: FieldParams, const N: usize> Add for UnivariateCoefficientBasis<P, N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut res = self;
        res.coefficients[0] = res.coefficients[0] + rhs.coefficients[0];
        res.coefficients[1] = res.coefficients[1] + rhs.coefficients[1];
        if N == 3 {
            res.coefficients[2] = res.coefficients[2] + rhs.coefficients[2];
        }
        res
    }
}

impl<P: FieldParams, const N: usize> AddAssign for UnivariateCoefficientBasis<P, N> {
    fn add_assign(&mut self, rhs: Self) {
        self.coefficients[0] = self.coefficients[0] + rhs.coefficients[0];
        self.coefficients[1] = self.coefficients[1] + rhs.coefficients[1];
        if N == 3 {
            self.coefficients[2] = self.coefficients[2] + rhs.coefficients[2];
        }
    }
}

impl<P: FieldParams, const N: usize> Sub for UnivariateCoefficientBasis<P, N> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut res = self;
        res.coefficients[0] = res.coefficients[0] - rhs.coefficients[0];
        res.coefficients[1] = res.coefficients[1] - rhs.coefficients[1];
        if N == 3 {
            res.coefficients[2] = res.coefficients[2] - rhs.coefficients[2];
        }
        res
    }
}

impl<P: FieldParams, const N: usize> SubAssign for UnivariateCoefficientBasis<P, N> {
    fn sub_assign(&mut self, rhs: Self) {
        self.coefficients[0] = self.coefficients[0] - rhs.coefficients[0];
        self.coefficients[1] = self.coefficients[1] - rhs.coefficients[1];
        if N == 3 {
            self.coefficients[2] = self.coefficients[2] - rhs.coefficients[2];
        }
    }
}

impl<P: FieldParams, const N: usize> Neg for UnivariateCoefficientBasis<P, N> {
    type Output = Self;
    fn neg(self) -> Self {
        let mut res = self;
        res.coefficients[0] = -res.coefficients[0];
        res.coefficients[1] = -res.coefficients[1];
        if N == 3 {
            res.coefficients[2] = -res.coefficients[2];
        }
        res
    }
}

// ---------------------------------------------------------------------------
// Scalar arithmetic
// ---------------------------------------------------------------------------

impl<P: FieldParams, const N: usize> UnivariateCoefficientBasis<P, N> {
    /// Add scalar to constant coefficient only.
    pub fn add_scalar(&self, s: Field<P>) -> Self {
        let mut res = *self;
        res.coefficients[0] = res.coefficients[0] + s;
        res
    }

    /// Subtract scalar from constant coefficient only.
    pub fn sub_scalar(&self, s: Field<P>) -> Self {
        let mut res = *self;
        res.coefficients[0] = res.coefficients[0] - s;
        res
    }

    /// Multiply all coefficients by scalar.
    pub fn mul_scalar(&self, s: Field<P>) -> Self {
        let mut res = *self;
        res.coefficients[0] = res.coefficients[0] * s;
        res.coefficients[1] = res.coefficients[1] * s;
        if N == 3 {
            res.coefficients[2] = res.coefficients[2] * s;
        }
        res
    }
}

// ---------------------------------------------------------------------------
// Conversions: coefficient basis → evaluation basis (Univariate)
// ---------------------------------------------------------------------------

impl<P: FieldParams> UnivariateCoefficientBasis<P, 2> {
    /// Convert degree-1 coefficient form to evaluation form on {0, 1, ..., M-1}.
    ///
    /// f(x) = a0 + a1*x → evaluations[i] = evaluations[i-1] + a1
    pub fn to_univariate<const M: usize>(&self) -> Univariate<P, M> {
        let mut result = Univariate::<P, M>::zero();
        result.evaluations[0] = self.coefficients[0];
        let mut prev = result.evaluations[0];
        let to_add = self.coefficients[1];
        for i in 1..M {
            prev = prev + to_add;
            result.evaluations[i] = prev;
        }
        result
    }
}

impl<P: FieldParams> UnivariateCoefficientBasis<P, 3> {
    /// Convert degree-2 coefficient form to evaluation form on {0, 1, ..., M-1}.
    ///
    /// Forward-difference layout: coefficients = [a0, a1+a2, a2].
    /// Uses second-order forward differences for O(M) evaluation.
    pub fn to_univariate<const M: usize>(&self) -> Univariate<P, M> {
        assert!(M >= 1);
        let mut result = Univariate::<P, M>::zero();
        let mut to_add = self.coefficients[1]; // a1 + a2
        let derivative = self.coefficients[2] + self.coefficients[2]; // 2*a2
        result.evaluations[0] = self.coefficients[0];
        let mut prev = result.evaluations[0];

        for i in 1..M - 1 {
            prev = prev + to_add;
            result.evaluations[i] = prev;
            to_add = to_add + derivative;
        }
        if M > 1 {
            prev = prev + to_add;
            result.evaluations[M - 1] = prev;
        }
        result
    }
}
