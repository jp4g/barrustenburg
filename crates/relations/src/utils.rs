//! Port of `utils.hpp` — `RelationUtils` utility functions for sumcheck.
//!
//! The C++ `RelationUtils<Flavor>` class provides utility methods for:
//! - Zeroing tuples of Univariates
//! - Adding nested tuples
//! - Scaling Univariates by challenges
//! - Accumulating relation evaluations
//!
//! In Rust, these are implemented as free functions and trait-based operations
//! since the Flavor parametrization is handled differently.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_polynomials::univariate::Univariate;

/// Zero out all evaluations in a Univariate.
#[inline]
pub fn zero_univariate<P: FieldParams, const N: usize>(u: &mut Univariate<P, N>) {
    for e in u.evaluations.iter_mut() {
        *e = Field::zero();
    }
}

/// Componentwise addition: `a += b` for two Univariates.
#[inline]
pub fn add_univariates<P: FieldParams, const N: usize>(
    a: &mut Univariate<P, N>,
    b: &Univariate<P, N>,
) {
    for (ai, bi) in a.evaluations.iter_mut().zip(b.evaluations.iter()) {
        *ai = *ai + *bi;
    }
}

/// Scale a Univariate by a field element: `u *= scalar`.
#[inline]
pub fn scale_univariate<P: FieldParams, const N: usize>(
    u: &mut Univariate<P, N>,
    scalar: &Field<P>,
) {
    for e in u.evaluations.iter_mut() {
        *e = *e * *scalar;
    }
}

/// Zero out an array of field values (verifier-side accumulators).
#[inline]
pub fn zero_values<P: FieldParams>(values: &mut [Field<P>]) {
    for v in values.iter_mut() {
        *v = Field::zero();
    }
}

/// Scale elements representing subrelation evaluations by separate challenges,
/// then sum them. The first subrelation is not scaled (implicitly scaled by 1).
///
/// Port of C++ `RelationUtils::scale_and_batch_elements`.
///
/// - `values`: flat array of all subrelation evaluations across all relations
/// - `subrelation_separators`: array of NUM_SUBRELATIONS - 1 challenges
///   (the first subrelation doesn't need a separator)
///
/// Returns: `values[0] + sum_{i>0} values[i] * separators[i-1]`
pub fn scale_and_batch_elements<P: FieldParams>(
    values: &[Field<P>],
    subrelation_separators: &[Field<P>],
) -> Field<P> {
    if values.is_empty() {
        return Field::zero();
    }

    let mut result = values[0];
    for (i, val) in values.iter().enumerate().skip(1) {
        result = result + *val * subrelation_separators[i - 1];
    }
    result
}

/// Macro for zeroing a tuple of Univariates (2 elements).
#[macro_export]
macro_rules! zero_univariate_tuple2 {
    ($tuple:expr) => {{
        $crate::utils::zero_univariate(&mut $tuple.0);
        $crate::utils::zero_univariate(&mut $tuple.1);
    }};
}

/// Macro for zeroing a tuple of Univariates (1 element).
#[macro_export]
macro_rules! zero_univariate_tuple1 {
    ($tuple:expr) => {{
        $crate::utils::zero_univariate(&mut $tuple.0);
    }};
}

/// Macro for adding two tuples of Univariates (2 elements each).
#[macro_export]
macro_rules! add_univariate_tuples2 {
    ($a:expr, $b:expr) => {{
        $crate::utils::add_univariates(&mut $a.0, &$b.0);
        $crate::utils::add_univariates(&mut $a.1, &$b.1);
    }};
}

/// Macro for adding two tuples of Univariates (1 element each).
#[macro_export]
macro_rules! add_univariate_tuples1 {
    ($a:expr, $b:expr) => {{
        $crate::utils::add_univariates(&mut $a.0, &$b.0);
    }};
}

/// Scale Univariates within a flat array by different challenges.
/// The first element is not scaled (implicitly multiplied by 1).
///
/// Port of C++ `RelationUtils::scale_univariates`.
pub fn scale_univariates_flat<P: FieldParams, const N: usize>(
    univariates: &mut [Univariate<P, N>],
    subrelation_separators: &[Field<P>],
) {
    for (i, u) in univariates.iter_mut().enumerate() {
        if i > 0 {
            scale_univariate(u, &subrelation_separators[i - 1]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::Bn254FrParams;

    type P = Bn254FrParams;
    type FF = Field<P>;

    #[test]
    fn test_zero_univariate() {
        let mut u = Univariate::<P, 4>::random();
        assert!(!u.is_zero());
        zero_univariate(&mut u);
        assert!(u.is_zero());
    }

    #[test]
    fn test_add_univariates() {
        let a = Univariate::<P, 3>::new([FF::from(1u64), FF::from(2u64), FF::from(3u64)]);
        let b = Univariate::<P, 3>::new([FF::from(10u64), FF::from(20u64), FF::from(30u64)]);
        let mut result = a.clone();
        add_univariates(&mut result, &b);
        assert_eq!(result.evaluations[0], FF::from(11u64));
        assert_eq!(result.evaluations[1], FF::from(22u64));
        assert_eq!(result.evaluations[2], FF::from(33u64));
    }

    #[test]
    fn test_scale_and_batch_elements() {
        // 3 subrelation values: v0, v1, v2
        // separators: s0, s1 (for subrelations 1 and 2)
        // result = v0 + v1*s0 + v2*s1
        let v0 = FF::from(5u64);
        let v1 = FF::from(7u64);
        let v2 = FF::from(11u64);
        let s0 = FF::from(2u64);
        let s1 = FF::from(3u64);

        let result = scale_and_batch_elements(&[v0, v1, v2], &[s0, s1]);
        // 5 + 7*2 + 11*3 = 5 + 14 + 33 = 52
        assert_eq!(result, FF::from(52u64));
    }

    #[test]
    fn test_scale_and_batch_single_element() {
        let v0 = FF::from(42u64);
        let result = scale_and_batch_elements(&[v0], &[]);
        assert_eq!(result, FF::from(42u64));
    }
}
