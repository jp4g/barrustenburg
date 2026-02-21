//! Port of `grand_product_delta.hpp` — compute_public_input_delta.
//!
//! Computes the correction term for the permutation argument due to public inputs.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

/// The separator constant for the permutation argument.
/// Ensures evaluations of `id_i` / `sigma_i` on the boolean hypercube don't overlap across wires.
///
/// C++ value: `constexpr uint32_t PERMUTATION_ARGUMENT_VALUE_SEPARATOR = 1 << 28;`
pub const PERMUTATION_ARGUMENT_VALUE_SEPARATOR: u64 = 1 << 28;

/// Compute the public input delta for the permutation argument.
///
/// Port of C++ `compute_public_input_delta<Flavor>()`.
///
/// The permutation argument breaks the copy-constraint cycles for public inputs by mapping
/// σ⁰(i) = -(i+1) instead of σ⁰(i) = n+i. This creates a "delta" factor that the verifier
/// must account for:
///
/// ```text
///   Δ = ∏ᵢ (γ + xᵢ + β⋅(n+i)) / ∏ᵢ (γ + xᵢ - β⋅(1+i))
/// ```
///
/// where n = PERMUTATION_ARGUMENT_VALUE_SEPARATOR and xᵢ are the public inputs.
pub fn compute_public_input_delta<P: FieldParams>(
    public_inputs: &[Field<P>],
    beta: Field<P>,
    gamma: Field<P>,
    offset: Field<P>,
) -> Field<P> {
    if public_inputs.is_empty() {
        return Field::one();
    }

    let mut numerator = Field::<P>::one();
    let mut denominator = Field::<P>::one();

    let separator = Field::<P>::from(PERMUTATION_ARGUMENT_VALUE_SEPARATOR);
    let mut numerator_acc = gamma + beta * (separator + offset);
    let mut denominator_acc = gamma - beta * (offset + Field::one());

    for i in 0..public_inputs.len() {
        numerator = numerator * (numerator_acc + public_inputs[i]); // γ + xᵢ + β(n+i)
        denominator = denominator * (denominator_acc + public_inputs[i]); // γ + xᵢ - β(1+i)

        // Skip the final iteration's acc update (optimization from C++)
        if i < public_inputs.len() - 1 {
            numerator_acc = numerator_acc + beta;
            denominator_acc = denominator_acc - beta;
        }
    }

    numerator * denominator.invert()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_compute_public_input_delta_empty() {
        let beta = Fr::random_element();
        let gamma = Fr::random_element();
        let delta = compute_public_input_delta::<P>(&[], beta, gamma, Fr::zero());
        assert_eq!(delta, Fr::one(), "Empty public inputs should give delta = 1");
    }

    #[test]
    fn test_compute_public_input_delta_single() {
        // With a single public input x, offset=0:
        // numerator = gamma + x + beta * SEPARATOR
        // denominator = gamma + x - beta
        // delta = numerator / denominator
        let x = Fr::from(42u64);
        let beta = Fr::from(3u64);
        let gamma = Fr::from(7u64);
        let separator = Fr::from(PERMUTATION_ARGUMENT_VALUE_SEPARATOR);

        let expected_num = gamma + x + beta * separator;
        let expected_den = gamma + x - beta;
        let expected = expected_num * expected_den.invert();

        let delta = compute_public_input_delta::<P>(&[x], beta, gamma, Fr::zero());
        assert_eq!(delta, expected);
    }

    #[test]
    fn test_compute_public_input_delta_multiple() {
        // With two public inputs x0, x1, offset=0:
        // numerator = (gamma + x0 + beta * n) * (gamma + x1 + beta * (n+1))
        // denominator = (gamma + x0 - beta) * (gamma + x1 - 2*beta)
        let x0 = Fr::from(10u64);
        let x1 = Fr::from(20u64);
        let beta = Fr::from(5u64);
        let gamma = Fr::from(11u64);
        let n = Fr::from(PERMUTATION_ARGUMENT_VALUE_SEPARATOR);

        let expected_num =
            (gamma + x0 + beta * n) * (gamma + x1 + beta * (n + Fr::one()));
        let expected_den =
            (gamma + x0 - beta) * (gamma + x1 - beta * Fr::from(2u64));
        let expected = expected_num * expected_den.invert();

        let delta = compute_public_input_delta::<P>(&[x0, x1], beta, gamma, Fr::zero());
        assert_eq!(delta, expected);
    }

    #[test]
    fn test_compute_public_input_delta_with_offset() {
        // With a single public input x, offset=2:
        // numerator = gamma + x + beta * (n + 2)
        // denominator = gamma + x - beta * 3
        let x = Fr::from(100u64);
        let beta = Fr::from(7u64);
        let gamma = Fr::from(13u64);
        let offset = Fr::from(2u64);
        let n = Fr::from(PERMUTATION_ARGUMENT_VALUE_SEPARATOR);

        let expected_num = gamma + x + beta * (n + offset);
        let expected_den = gamma + x - beta * (offset + Fr::one());
        let expected = expected_num * expected_den.invert();

        let delta = compute_public_input_delta::<P>(&[x], beta, gamma, offset);
        assert_eq!(delta, expected);
    }

    #[test]
    fn test_compute_public_input_delta_random() {
        // Verify with random values that the product formula is correct
        let num_inputs = 5;
        let public_inputs: Vec<Fr> = (0..num_inputs).map(|_| Fr::random_element()).collect();
        let beta = Fr::random_element();
        let gamma = Fr::random_element();
        let offset = Fr::zero();
        let n = Fr::from(PERMUTATION_ARGUMENT_VALUE_SEPARATOR);

        // Compute expected manually
        let mut expected_num = Fr::one();
        let mut expected_den = Fr::one();
        for i in 0..num_inputs {
            let idx = Fr::from(i as u64);
            expected_num =
                expected_num * (gamma + public_inputs[i] + beta * (n + offset + idx));
            expected_den =
                expected_den * (gamma + public_inputs[i] - beta * (offset + Fr::one() + idx));
        }
        let expected = expected_num * expected_den.invert();

        let delta = compute_public_input_delta::<P>(&public_inputs, beta, gamma, offset);
        assert_eq!(delta, expected);
    }
}
