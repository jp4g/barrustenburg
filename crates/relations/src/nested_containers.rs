//! Port of `nested_containers.hpp` — type constructors for tuples of Univariates.
//!
//! The C++ code uses `TupleOfContainersOverArray` to create a `flat_tuple::tuple`
//! of `Univariate<FF, N>` values where each N comes from an array element.
//!
//! In Rust, we can't have variadic tuples keyed by const arrays at the type level.
//! Instead, we provide:
//! 1. A macro to generate tuple types from literal lengths
//! 2. An `ArrayOfValues` type (uniform array for verifier-side accumulation)
//! 3. Utility functions to work with these containers

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

/// `ArrayOfValues<FF, N>` — uniform array of FF values, one per subrelation.
///
/// Port of C++ `ArrayOfValues<FF, LENGTHS>` which is
/// `HomogeneousTupleToArray<TupleOfValues<FF, LENGTHS>>` = `std::array<FF, N>`.
///
/// Used by the verifier: each subrelation accumulates into a single FF value.
pub type ArrayOfValues<P, const N: usize> = [Field<P>; N];

/// Zero-initialize an array of field values.
pub fn zero_array_of_values<P: FieldParams, const N: usize>() -> ArrayOfValues<P, N> {
    [Field::zero(); N]
}

/// Macro to create a tuple of Univariates from a list of lengths.
///
/// # Example
/// ```ignore
/// // Creates (Univariate<P, 3>,)
/// type Rel1Accum = tuple_of_univariates!(P; 3);
///
/// // Creates (Univariate<P, 2>, Univariate<P, 5>)
/// type Rel2Accum = tuple_of_univariates!(P; 2, 5);
/// ```
#[macro_export]
macro_rules! tuple_of_univariates {
    ($P:ty; $($len:expr),+ $(,)?) => {
        ( $(bbrs_polynomials::univariate::Univariate<$P, $len>,)+ )
    };
}

/// Macro to create a zero-initialized tuple of Univariates.
#[macro_export]
macro_rules! zero_tuple_of_univariates {
    ($P:ty; $($len:expr),+ $(,)?) => {
        ( $(bbrs_polynomials::univariate::Univariate::<$P, $len>::zero(),)+ )
    };
}

#[cfg(test)]
mod tests {
    use bbrs_ecc::curves::bn254::Bn254FrParams;
    use bbrs_polynomials::univariate::Univariate;

    type P = Bn254FrParams;

    /// Port of C++ `NestedContainers::Univariate` test.
    ///
    /// The C++ test creates a `TupleOfUnivariates<FF, {0, 1, 2}>` and verifies
    /// each element matches a zero-initialized Univariate of the corresponding length.
    ///
    /// Note: Univariate<FF, 0> has no practical use but the C++ test includes it.
    /// In Rust, Univariate<P, 0> is a zero-length array which is valid.
    #[test]
    fn test_nested_containers_univariate() {
        // C++ LENGTHS = {0, 1, 2}
        // We create the tuple manually since Rust doesn't support const-array-driven tuples.
        let tuple: (
            Univariate<P, 0>,
            Univariate<P, 1>,
            Univariate<P, 2>,
        ) = (
            Univariate::zero(),
            Univariate::zero(),
            Univariate::zero(),
        );

        assert_eq!(tuple.0, Univariate::<P, 0>::zero());
        assert_eq!(tuple.1, Univariate::<P, 1>::zero());
        assert_eq!(tuple.2, Univariate::<P, 2>::zero());

        // Also test with the macro
        let from_macro = zero_tuple_of_univariates!(P; 0, 1, 2);
        assert_eq!(from_macro.0, Univariate::<P, 0>::zero());
        assert_eq!(from_macro.1, Univariate::<P, 1>::zero());
        assert_eq!(from_macro.2, Univariate::<P, 2>::zero());
    }
}
