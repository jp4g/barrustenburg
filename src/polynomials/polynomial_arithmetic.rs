use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;

/// Evaluate polynomial at `z` using Horner's method.
/// Coefficients are in order [c_0, c_1, ..., c_{n-1}] so p(X) = c_0 + c_1*X + ... + c_{n-1}*X^{n-1}.
pub fn evaluate<P: FieldParams>(coeffs: &[Field<P>], z: &Field<P>) -> Field<P> {
    let n = coeffs.len();
    if n == 0 {
        return Field::zero();
    }
    let mut result = coeffs[n - 1];
    for i in (0..n - 1).rev() {
        result = result * *z + coeffs[i];
    }
    result
}

/// Synthetic division: divide the polynomial (given by `coeffs`) by (X - root) in-place.
/// After the call the leading coefficient is zero and the remaining n-1 coefficients
/// hold the quotient.  Mirrors the C++ `factor_roots`.
pub fn factor_roots<P: FieldParams>(coeffs: &mut [Field<P>], root: &Field<P>) {
    let n = coeffs.len();
    if n == 0 {
        return;
    }
    let mut work = coeffs[n - 1];
    coeffs[n - 1] = Field::zero();
    // Walk from second-highest down to 0 (wrapping arithmetic handles the i == 0 case).
    for i in (0..n - 1).rev() {
        let temp = coeffs[i];
        coeffs[i] = work;
        work = temp + work * *root;
    }
}

/// Lagrange interpolation: given evaluations `evaluations[i]` at `points[i]`, return
/// the coefficient representation of the unique polynomial of degree <= n-1 passing
/// through all points.
pub fn compute_efficient_interpolation<P: FieldParams>(
    evaluations: &[Field<P>],
    points: &[Field<P>],
) -> Vec<Field<P>> {
    let n = evaluations.len();
    assert_eq!(n, points.len());
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![evaluations[0]];
    }

    // Step 1: denominators[i] = prod_{j != i} (points[i] - points[j])
    let mut denoms = vec![Field::<P>::one(); n];
    for i in 0..n {
        for j in 0..n {
            if i != j {
                denoms[i] = denoms[i] * (points[i] - points[j]);
            }
        }
    }

    // Step 2: scaled[i] = evaluations[i] / denominators[i]
    let mut scaled = vec![Field::<P>::zero(); n];
    for i in 0..n {
        scaled[i] = evaluations[i] * denoms[i].invert();
    }

    // Step 3: result = sum_i scaled[i] * prod_{j != i}(X - points[j])
    let mut result = vec![Field::<P>::zero(); n];
    for i in 0..n {
        // Build basis polynomial L_i(X) = prod_{j != i}(X - points[j])
        let mut basis = vec![Field::<P>::zero(); n];
        basis[0] = Field::one();
        let mut deg = 0usize;
        for j in 0..n {
            if j == i {
                continue;
            }
            deg += 1;
            // Multiply current basis by (X - points[j]), high-to-low to avoid overwrites
            for k in (1..=deg).rev() {
                basis[k] = basis[k - 1] - basis[k] * points[j];
            }
            basis[0] = -basis[0] * points[j];
        }
        // Accumulate scaled[i] * basis into result
        for k in 0..n {
            result[k] = result[k] + scaled[i] * basis[k];
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecc::curves::bn254::Bn254FrParams;

    type Fr = Field<Bn254FrParams>;

    #[test]
    fn test_evaluate_empty() {
        let coeffs: Vec<Fr> = vec![];
        let z = Fr::from(5u64);
        assert_eq!(evaluate(&coeffs, &z), Fr::zero());
    }

    #[test]
    fn test_evaluate_constant() {
        let coeffs = vec![Fr::from(7u64)];
        let z = Fr::from(42u64);
        assert_eq!(evaluate(&coeffs, &z), Fr::from(7u64));
    }

    #[test]
    fn test_evaluate_linear() {
        // p(X) = 3 + 2*X, p(5) = 13
        let coeffs = vec![Fr::from(3u64), Fr::from(2u64)];
        let z = Fr::from(5u64);
        assert_eq!(evaluate(&coeffs, &z), Fr::from(13u64));
    }

    #[test]
    fn test_factor_roots() {
        // p(X) = (X - 3)(X - 7) = X^2 - 10X + 21 = 21 - 10X + X^2
        let mut coeffs = vec![Fr::from(21u64), -Fr::from(10u64), Fr::one()];
        factor_roots(&mut coeffs, &Fr::from(3u64));
        // Quotient should be (X - 7) = -7 + X, with leading coeff zeroed
        assert_eq!(coeffs[0], -Fr::from(7u64));
        assert_eq!(coeffs[1], Fr::one());
        assert_eq!(coeffs[2], Fr::zero());
    }

    #[test]
    fn test_interpolation_matches_evaluation() {
        let points = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let evals = vec![Fr::from(6u64), Fr::from(11u64), Fr::from(18u64)];
        let coeffs = compute_efficient_interpolation(&evals, &points);
        // Verify round-trip: evaluating at each point recovers the evaluation
        for i in 0..points.len() {
            assert_eq!(evaluate(&coeffs, &points[i]), evals[i]);
        }
    }
}
