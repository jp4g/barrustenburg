use bbrs_ecc::curves::bn254::Bn254FrParams;
use bbrs_ecc::fields::field::Field;
use crate::barycentric::BarycentricData;
use crate::eq_polynomial::{ProverEqPolynomial, VerifierEqPolynomial};
use crate::evaluation_domain::EvaluationDomain;
use crate::gate_separator::GateSeparatorPolynomial;
use crate::polynomial::Polynomial;
use crate::polynomial_arithmetic;
use crate::polynomial_span::PolynomialSpan;
use crate::row_disabling_polynomial::RowDisablingPolynomial;
use crate::univariate::{Univariate, UnivariateView};
use crate::univariate_coefficient_basis::UnivariateCoefficientBasis;

type Fr = Field<Bn254FrParams>;

// ═══════════════════════════════════════════════════════════════════════════════
// Polynomial tests (from polynomial.test.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn polynomial_shifted() {
    let n = 8usize;
    // shiftable(n) creates start_index=1, data at [1..n)
    let mut p: Polynomial<Bn254FrParams> = Polynomial::shiftable(n);
    for i in 1..n {
        *p.at_mut(i) = Fr::from((i + 1) as u64);
    }
    assert_eq!(p.get(0), Fr::zero()); // virtual zero at index 0
    assert_eq!(p.get(1), Fr::from(2u64));

    let s = p.shifted();
    // shifted() decrements start: 1 -> 0, so s[i] = p[i+1]
    assert_eq!(s.start_index(), 0);
    assert_eq!(s.size(), n - 1);
    for i in 0..n - 1 {
        assert_eq!(s.get(i), Fr::from((i + 2) as u64)); // s[i] = p[i+1]
    }
}

#[test]
fn polynomial_clone() {
    let n = 8usize;
    let mut p: Polynomial<Bn254FrParams> = Polynomial::new(n, n, 0);
    for i in 0..n {
        *p.at_mut(i) = Fr::from(i as u64);
    }
    let q = p.clone();
    for i in 0..n {
        assert_eq!(q.get(i), Fr::from(i as u64));
    }
    // Mutating p doesn't affect q (Rust ownership)
    *p.at_mut(0) = Fr::from(999u64);
    assert_eq!(q.get(0), Fr::zero());
}

#[test]
fn polynomial_indices() {
    let size = 4usize;
    let virtual_size = 16usize;
    let start = 3usize;
    let p: Polynomial<Bn254FrParams> = Polynomial::new(size, virtual_size, start);
    assert_eq!(p.start_index(), start);
    assert_eq!(p.end_index(), start + size);
    assert_eq!(p.size(), size);
    assert_eq!(p.virtual_size(), virtual_size);
}

#[test]
fn polynomial_add_scaled_edge_conditions() {
    // Test with offset polynomials
    let mut p1: Polynomial<Bn254FrParams> = Polynomial::new(4, 8, 0);
    for i in 0..4 {
        *p1.at_mut(i) = Fr::from(1u64);
    }
    // Create other starting at index 2
    let other_data = vec![Fr::from(2u64), Fr::from(3u64)];
    let span = PolynomialSpan::new(&other_data, 2);
    p1.add_scaled(&span, Fr::from(10u64));
    assert_eq!(p1.get(0), Fr::from(1u64));
    assert_eq!(p1.get(1), Fr::from(1u64));
    assert_eq!(p1.get(2), Fr::from(21u64)); // 1 + 2*10
    assert_eq!(p1.get(3), Fr::from(31u64)); // 1 + 3*10
}

#[test]
fn polynomial_operator_add_edge_conditions() {
    let mut p: Polynomial<Bn254FrParams> = Polynomial::new(4, 8, 0);
    for i in 0..4 {
        *p.at_mut(i) = Fr::from(10u64);
    }
    let other_data = vec![Fr::from(5u64), Fr::from(7u64)];
    let span = PolynomialSpan::new(&other_data, 1);
    p += span;
    assert_eq!(p.get(0), Fr::from(10u64));
    assert_eq!(p.get(1), Fr::from(15u64)); // 10 + 5
    assert_eq!(p.get(2), Fr::from(17u64)); // 10 + 7
    assert_eq!(p.get(3), Fr::from(10u64));
}

#[test]
fn polynomial_operator_subtract_edge_conditions() {
    let mut p: Polynomial<Bn254FrParams> = Polynomial::new(4, 8, 0);
    for i in 0..4 {
        *p.at_mut(i) = Fr::from(20u64);
    }
    let other_data = vec![Fr::from(3u64), Fr::from(7u64)];
    let span = PolynomialSpan::new(&other_data, 1);
    p -= span;
    assert_eq!(p.get(0), Fr::from(20u64));
    assert_eq!(p.get(1), Fr::from(17u64)); // 20 - 3
    assert_eq!(p.get(2), Fr::from(13u64)); // 20 - 7
    assert_eq!(p.get(3), Fr::from(20u64));
}

#[test]
fn polynomial_full() {
    let p: Polynomial<Bn254FrParams> =
        Polynomial::from_coefficients(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)], 8);
    let f = p.full();
    assert_eq!(f.start_index(), 0);
    assert_eq!(f.size(), 8);
    assert_eq!(f.get(0), Fr::from(1u64));
    assert_eq!(f.get(1), Fr::from(2u64));
    assert_eq!(f.get(2), Fr::from(3u64));
    for i in 3..8 {
        assert_eq!(f.get(i), Fr::zero());
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// polynomial_arithmetic tests
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn polynomial_arithmetic_evaluate() {
    // p(X) = 1 + 2X + 3X^2 + 4X^3
    let coeffs = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)];
    let z = Fr::from(2u64);
    // p(2) = 1 + 4 + 12 + 32 = 49
    assert_eq!(polynomial_arithmetic::evaluate(&coeffs, &z), Fr::from(49u64));
}

#[test]
fn evaluation_domain_construction() {
    let domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(16);
    assert_eq!(domain.size, 16);
    assert_eq!(domain.log2_size, 4);
    // root^16 should equal 1
    let mut r = domain.root;
    for _ in 0..16 {
        r = r * domain.root;
    }
    // Actually root is a 16th root of unity, so root^16 == 1
    let root_to_16 = {
        let mut acc = Fr::one();
        for _ in 0..16 {
            acc = acc * domain.root;
        }
        acc
    };
    assert_eq!(root_to_16, Fr::one());
    // domain * domain_inverse == 1
    assert_eq!(domain.domain * domain.domain_inverse, Fr::one());
}

#[test]
fn evaluation_domain_domain_roots() {
    let domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(32);
    // root is a 32nd root of unity
    let mut r = domain.root;
    for _ in 1..32 {
        assert_ne!(r, Fr::one(), "root should have order exactly 32");
        r = r * domain.root;
    }
    assert_eq!(r, Fr::one());
}

#[test]
fn evaluation_domain_roots() {
    let mut domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(16);
    domain.compute_lookup_table();
    let roots = domain.get_round_roots();
    assert_eq!(roots.len(), 4); // log2(16) = 4 rounds

    // Round 0: table_size = 1, should be [1]
    assert_eq!(roots[0].len(), 1);
    assert_eq!(roots[0][0], Fr::one());

    // Round 1: table_size = 2, root of order 4 (square root of domain root squared)
    assert_eq!(roots[1].len(), 2);
    assert_eq!(roots[1][0], Fr::one());
    // roots[1][1] should be a 4th root of unity
    let r4 = roots[1][1];
    let r4_cubed = r4 * r4 * r4 * r4;
    assert_eq!(r4_cubed, Fr::one());

    // Round 3: table_size = 8, root of order 16 (= domain root itself)
    assert_eq!(roots[3].len(), 8);
    assert_eq!(roots[3][0], Fr::one());
    assert_eq!(roots[3][1], domain.root);
}

#[test]
fn compute_efficient_interpolation() {
    let n = 8usize;
    let points: Vec<Fr> = (1..=n as u64).map(Fr::from).collect();
    // Evaluations of p(X) = X^2 + 2X + 1
    let evals: Vec<Fr> = points.iter().map(|&x| x * x + Fr::from(2u64) * x + Fr::one()).collect();
    let coeffs = polynomial_arithmetic::compute_efficient_interpolation(&evals, &points);
    // Verify round-trip
    for i in 0..n {
        assert_eq!(polynomial_arithmetic::evaluate(&coeffs, &points[i]), evals[i]);
    }
}

#[test]
fn compute_efficient_interpolation_domain_with_zero() {
    // Include x=0 in the domain
    let points = vec![Fr::zero(), Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
    let evals = vec![Fr::from(5u64), Fr::from(7u64), Fr::from(13u64), Fr::from(23u64)];
    let coeffs = polynomial_arithmetic::compute_efficient_interpolation(&evals, &points);
    for i in 0..points.len() {
        assert_eq!(polynomial_arithmetic::evaluate(&coeffs, &points[i]), evals[i]);
    }
}

#[test]
fn interpolation_constructor_single() {
    let points = vec![Fr::from(42u64)];
    let evals = vec![Fr::from(99u64)];
    let p: Polynomial<Bn254FrParams> = Polynomial::from_interpolation(&points, &evals, 4);
    assert_eq!(p.evaluate(&Fr::from(42u64)), Fr::from(99u64));
    // Constant polynomial — same value everywhere
    assert_eq!(p.evaluate(&Fr::from(100u64)), Fr::from(99u64));
}

#[test]
fn interpolation_constructor() {
    let n = 4usize;
    let points: Vec<Fr> = (0..n as u64).map(Fr::from).collect();
    // p(X) = 2X^3 - X + 5
    let evals: Vec<Fr> = points
        .iter()
        .map(|&x| Fr::from(2u64) * x * x * x - x + Fr::from(5u64))
        .collect();
    let p: Polynomial<Bn254FrParams> = Polynomial::from_interpolation(&points, &evals, 8);
    for i in 0..n {
        assert_eq!(p.evaluate(&points[i]), evals[i]);
    }
    // Also check at a point not in the domain
    let extra = Fr::from(10u64);
    let expected = Fr::from(2u64) * extra * extra * extra - extra + Fr::from(5u64);
    assert_eq!(p.evaluate(&extra), expected);
}

#[test]
fn evaluate_mle() {
    // 3-variable MLE (8 coefficients)
    let n = 3usize;
    let poly_size = 1usize << n; // 8
    let mut coeffs = vec![Fr::zero(); poly_size];
    for i in 0..poly_size {
        coeffs[i] = Fr::from((i * i + 1) as u64); // arbitrary values
    }
    let p: Polynomial<Bn254FrParams> = Polynomial::from_coefficients(coeffs.clone(), poly_size);

    // Evaluate at a boolean point (1, 0, 1) = index 5 (binary 101)
    let result = p.evaluate_mle(&[Fr::one(), Fr::zero(), Fr::one()], false);
    assert_eq!(result, coeffs[5]); // f(1,0,1) = coeffs[1*1 + 0*2 + 1*4] = coeffs[5]

    // Evaluate at (0,0,0) = index 0
    let result = p.evaluate_mle(&[Fr::zero(), Fr::zero(), Fr::zero()], false);
    assert_eq!(result, coeffs[0]);

    // Evaluate at (1,1,1) = index 7
    let result = p.evaluate_mle(&[Fr::one(), Fr::one(), Fr::one()], false);
    assert_eq!(result, coeffs[7]);
}

#[test]
fn move_construct_and_assign() {
    // In Rust, "move" is the default. Verify ownership transfer works.
    let p: Polynomial<Bn254FrParams> = Polynomial::from_coefficients(
        vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)],
        4,
    );
    let q = p; // move
    assert_eq!(q.get(0), Fr::from(1u64));
    assert_eq!(q.size(), 3);
}

#[test]
fn default_construct_then_assign() {
    let mut p: Polynomial<Bn254FrParams> = Polynomial::new(0, 0, 0);
    assert!(p.is_empty());
    p = Polynomial::from_coefficients(vec![Fr::from(7u64), Fr::from(8u64)], 4);
    assert_eq!(p.size(), 2);
    assert_eq!(p.get(0), Fr::from(7u64));
}

// ═══════════════════════════════════════════════════════════════════════════════
// Univariate tests (from univariate.test.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn univariate_constructors() {
    let u = Univariate::<Bn254FrParams, 4>::new([
        Fr::from(1u64),
        Fr::from(2u64),
        Fr::from(3u64),
        Fr::from(4u64),
    ]);
    assert_eq!(u.value_at(0), Fr::from(1u64));
    assert_eq!(u.value_at(3), Fr::from(4u64));

    let s = Univariate::<Bn254FrParams, 3>::from_scalar(Fr::from(5u64));
    for i in 0..3 {
        assert_eq!(s.value_at(i), Fr::from(5u64));
    }

    let z = Univariate::<Bn254FrParams, 4>::zero();
    assert!(z.is_zero());
}

#[test]
fn univariate_addition() {
    let a = Univariate::<Bn254FrParams, 3>::new([Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]);
    let b = Univariate::<Bn254FrParams, 3>::new([
        Fr::from(10u64),
        Fr::from(20u64),
        Fr::from(30u64),
    ]);
    let c = a + b;
    assert_eq!(c.value_at(0), Fr::from(11u64));
    assert_eq!(c.value_at(1), Fr::from(22u64));
    assert_eq!(c.value_at(2), Fr::from(33u64));
}

#[test]
fn univariate_multiplication() {
    let a = Univariate::<Bn254FrParams, 3>::new([Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)]);
    let b = Univariate::<Bn254FrParams, 3>::new([
        Fr::from(10u64),
        Fr::from(10u64),
        Fr::from(10u64),
    ]);
    let c = a * b;
    assert_eq!(c.value_at(0), Fr::from(20u64));
    assert_eq!(c.value_at(1), Fr::from(30u64));
    assert_eq!(c.value_at(2), Fr::from(40u64));
}

#[test]
fn univariate_view_from_univariate() {
    let u = Univariate::<Bn254FrParams, 3>::new([Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]);
    let v = UnivariateView::from(&u);
    assert_eq!(v.value_at(0), Fr::from(1u64));
    assert_eq!(v.value_at(2), Fr::from(3u64));
}

#[test]
fn univariate_from_view() {
    let u = Univariate::<Bn254FrParams, 3>::new([Fr::from(5u64), Fr::from(6u64), Fr::from(7u64)]);
    let v = UnivariateView::from(&u);
    // Construct a new Univariate from the view's data
    let u2 = Univariate::<Bn254FrParams, 3>::new(*v.evaluations);
    assert_eq!(u2, u);
}

#[test]
fn univariate_view_addition() {
    let a = Univariate::<Bn254FrParams, 3>::new([Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]);
    let b = Univariate::<Bn254FrParams, 3>::new([
        Fr::from(10u64),
        Fr::from(20u64),
        Fr::from(30u64),
    ]);
    let va = UnivariateView::from(&a);
    let c = va + b;
    assert_eq!(c.value_at(0), Fr::from(11u64));
    assert_eq!(c.value_at(1), Fr::from(22u64));
    assert_eq!(c.value_at(2), Fr::from(33u64));
}

#[test]
fn univariate_view_subtraction() {
    let a = Univariate::<Bn254FrParams, 3>::new([
        Fr::from(10u64),
        Fr::from(20u64),
        Fr::from(30u64),
    ]);
    let b = Univariate::<Bn254FrParams, 3>::new([Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]);
    let va = UnivariateView::from(&a);
    let c = va - b;
    assert_eq!(c.value_at(0), Fr::from(9u64));
    assert_eq!(c.value_at(1), Fr::from(18u64));
    assert_eq!(c.value_at(2), Fr::from(27u64));
}

#[test]
fn univariate_view_multiplication() {
    let a = Univariate::<Bn254FrParams, 3>::new([Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)]);
    let b = Univariate::<Bn254FrParams, 3>::new([Fr::from(5u64), Fr::from(5u64), Fr::from(5u64)]);
    let va = UnivariateView::from(&a);
    let c = va * b;
    assert_eq!(c.value_at(0), Fr::from(10u64));
    assert_eq!(c.value_at(1), Fr::from(15u64));
    assert_eq!(c.value_at(2), Fr::from(20u64));
}

// ═══════════════════════════════════════════════════════════════════════════════
// Univariate extension / evaluation tests
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn univariate_extend_linear() {
    // f(0) = 3, f(1) = 7 → linear: f(x) = 3 + 4x
    let u = Univariate::<Bn254FrParams, 2>::new([Fr::from(3u64), Fr::from(7u64)]);
    let ext = u.extend_to::<5>();
    assert_eq!(ext.value_at(0), Fr::from(3u64));
    assert_eq!(ext.value_at(1), Fr::from(7u64));
    assert_eq!(ext.value_at(2), Fr::from(11u64));
    assert_eq!(ext.value_at(3), Fr::from(15u64));
    assert_eq!(ext.value_at(4), Fr::from(19u64));
}

#[test]
fn univariate_extend_quadratic() {
    // f(0) = 1, f(1) = 4, f(2) = 9 → f(x) = (x+1)^2
    let u = Univariate::<Bn254FrParams, 3>::new([Fr::from(1u64), Fr::from(4u64), Fr::from(9u64)]);
    let ext = u.extend_to::<6>();
    assert_eq!(ext.value_at(3), Fr::from(16u64)); // (3+1)^2
    assert_eq!(ext.value_at(4), Fr::from(25u64)); // (4+1)^2
    assert_eq!(ext.value_at(5), Fr::from(36u64)); // (5+1)^2
}

#[test]
fn univariate_evaluate_at_point() {
    // f(0) = 1, f(1) = 3, f(2) = 7 → f(x) = x^2 + x + 1
    let u = Univariate::<Bn254FrParams, 3>::new([Fr::from(1u64), Fr::from(3u64), Fr::from(7u64)]);
    // f(5) = 25 + 5 + 1 = 31
    assert_eq!(u.evaluate(Fr::from(5u64)), Fr::from(31u64));
    // f(0) should return exact domain value
    assert_eq!(u.evaluate(Fr::from(0u64)), Fr::from(1u64));
    assert_eq!(u.evaluate(Fr::from(1u64)), Fr::from(3u64));
}

// ═══════════════════════════════════════════════════════════════════════════════
// batch_invert test
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn batch_invert() {
    let mut elements: Vec<Fr> = (1..=8).map(|i| Fr::from(i as u64)).collect();
    let originals = elements.clone();
    Fr::batch_invert(&mut elements);
    // Each element should now be the inverse of the original
    for i in 0..elements.len() {
        assert_eq!(elements[i] * originals[i], Fr::one());
    }
    // Test with zeros in the mix
    let mut with_zeros = vec![Fr::from(3u64), Fr::zero(), Fr::from(7u64), Fr::zero(), Fr::from(11u64)];
    let orig_with_zeros = with_zeros.clone();
    Fr::batch_invert(&mut with_zeros);
    assert_eq!(with_zeros[0] * orig_with_zeros[0], Fr::one());
    assert!(with_zeros[1].is_zero()); // zero stays zero
    assert_eq!(with_zeros[2] * orig_with_zeros[2], Fr::one());
    assert!(with_zeros[3].is_zero()); // zero stays zero
    assert_eq!(with_zeros[4] * orig_with_zeros[4], Fr::one());
}

// ═══════════════════════════════════════════════════════════════════════════════
// FFT / IFFT / Coset tests (from polynomial_arithmetic.test.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn fft_with_small_degree() {
    // C++ test: fft_with_small_degree — FFT[i] == evaluate(poly, ω^i) for all i
    let n = 16usize;
    let poly: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let mut fft_transform = poly.clone();

    let mut domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(n);
    domain.compute_lookup_table();
    polynomial_arithmetic::fft(&mut fft_transform, &domain);

    let mut work_root = Fr::one();
    for i in 0..n {
        let expected = polynomial_arithmetic::evaluate(&poly, &work_root);
        assert_eq!(fft_transform[i], expected, "mismatch at index {}", i);
        work_root = work_root * domain.root;
    }
}

#[test]
fn split_polynomial_fft() {
    // C++ test: split_polynomial_fft — split FFT matches monolithic FFT
    let n = 256usize;
    let poly: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();

    let mut fft_mono = poly.clone();

    // Split into 4 sub-arrays: poly_split[j] = poly[j*poly_size..(j+1)*poly_size]
    let num_poly = 4usize;
    let poly_size = n / num_poly;
    let mut split_data: Vec<Vec<Fr>> = (0..num_poly)
        .map(|j| poly[j * poly_size..(j + 1) * poly_size].to_vec())
        .collect();

    let mut domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(n);
    domain.compute_lookup_table();

    polynomial_arithmetic::fft(&mut fft_mono, &domain);

    {
        let mut split_refs: Vec<&mut [Fr]> = split_data.iter_mut().map(|v| v.as_mut_slice()).collect();
        polynomial_arithmetic::fft_split(&mut split_refs, &domain);
    }

    // Verify split matches monolithic
    for i in 0..n {
        let poly_idx = i / poly_size;
        let elem_idx = i % poly_size;
        assert_eq!(
            split_data[poly_idx][elem_idx], fft_mono[i],
            "split FFT mismatch at index {}",
            i
        );
    }

    // Also verify monolithic matches direct evaluation
    let mut work_root = Fr::one();
    for i in 0..n {
        let expected = polynomial_arithmetic::evaluate(&poly, &work_root);
        assert_eq!(fft_mono[i], expected, "FFT vs evaluate mismatch at {}", i);
        work_root = work_root * domain.root;
    }
}

#[test]
fn fft_ifft_roundtrip_large() {
    // C++ test: basic_fft — FFT→IFFT recovers original (n = 2^14)
    let n = 1 << 14;
    let original: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let mut result = original.clone();

    let mut domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(n);
    domain.compute_lookup_table();

    polynomial_arithmetic::fft(&mut result, &domain);
    polynomial_arithmetic::ifft(&mut result, &domain);

    for i in 0..n {
        assert_eq!(result[i], original[i], "roundtrip mismatch at index {}", i);
    }
}

#[test]
fn fft_ifft_consistency() {
    // C++ test: fft_ifft_consistency — FFT→IFFT roundtrip (n = 256)
    let n = 256usize;
    let original: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let mut result = original.clone();

    let mut domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(n);
    domain.compute_lookup_table();

    polynomial_arithmetic::fft(&mut result, &domain);
    polynomial_arithmetic::ifft(&mut result, &domain);

    for i in 0..n {
        assert_eq!(result[i], original[i], "mismatch at index {}", i);
    }
}

#[test]
fn split_polynomial_fft_ifft_consistency() {
    // C++ test: split_polynomial_fft_ifft_consistency
    let num_poly = 4usize;
    let poly_size = 256usize;
    let n = num_poly * poly_size;

    let mut split_data: Vec<Vec<Fr>> = (0..num_poly)
        .map(|_| (0..poly_size).map(|_| Fr::random_element()).collect())
        .collect();
    let original: Vec<Vec<Fr>> = split_data.clone();

    let mut domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(n);
    domain.compute_lookup_table();

    {
        let mut split_refs: Vec<&mut [Fr]> = split_data.iter_mut().map(|v| v.as_mut_slice()).collect();
        polynomial_arithmetic::fft_split(&mut split_refs, &domain);
        polynomial_arithmetic::ifft_split(&mut split_refs, &domain);
    }

    for j in 0..num_poly {
        for i in 0..poly_size {
            assert_eq!(
                split_data[j][i], original[j][i],
                "split roundtrip mismatch at poly {}, index {}",
                j, i
            );
        }
    }
}

#[test]
fn fft_coset_ifft_consistency() {
    // C++ test: fft_coset_ifft_consistency — coset_fft→coset_ifft roundtrip
    let n = 256usize;
    let original: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let mut result = original.clone();

    let mut domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(n);
    domain.compute_lookup_table();

    // Verify generator * generator_inverse == 1
    assert_eq!(domain.generator * domain.generator_inverse, Fr::one());

    polynomial_arithmetic::coset_fft(&mut result, &domain);
    polynomial_arithmetic::coset_ifft(&mut result, &domain);

    for i in 0..n {
        assert_eq!(result[i], original[i], "coset roundtrip mismatch at index {}", i);
    }
}

#[test]
fn split_polynomial_fft_coset_ifft_consistency() {
    // C++ test: split_polynomial_fft_coset_ifft_consistency
    let num_poly = 4usize;
    let poly_size = 256usize;
    let n = num_poly * poly_size;

    let mut split_data: Vec<Vec<Fr>> = (0..num_poly)
        .map(|_| (0..poly_size).map(|_| Fr::random_element()).collect())
        .collect();
    let original: Vec<Vec<Fr>> = split_data.clone();

    let mut domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(n);
    domain.compute_lookup_table();

    {
        let mut split_refs: Vec<&mut [Fr]> = split_data.iter_mut().map(|v| v.as_mut_slice()).collect();
        polynomial_arithmetic::coset_fft_split(&mut split_refs, &domain);
        polynomial_arithmetic::coset_ifft_split(&mut split_refs, &domain);
    }

    for j in 0..num_poly {
        for i in 0..poly_size {
            assert_eq!(
                split_data[j][i], original[j][i],
                "split coset roundtrip mismatch at poly {}, index {}",
                j, i
            );
        }
    }
}

#[test]
fn fft_coset_ifft_cross_consistency() {
    // C++ test: fft_coset_ifft_cross_consistency
    // 3 polys of sizes n, 2n, 4n. Zero-pad all to 4n, coset-fft each, sum, coset-ifft.
    let n = 2usize;

    // poly_a/b/c all have same first n coefficients; higher coefficients are zero
    let base: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let expected: Vec<Fr> = (0..n)
        .map(|i| base[i] + base[i] + base[i])
        .collect();

    let mut poly_a = vec![Fr::zero(); 4 * n];
    let mut poly_b = vec![Fr::zero(); 4 * n];
    let mut poly_c = vec![Fr::zero(); 4 * n];
    for i in 0..n {
        poly_a[i] = base[i];
        poly_b[i] = base[i];
        poly_c[i] = base[i];
    }

    let mut small_domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(n);
    let mut mid_domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(2 * n);
    let mut large_domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(4 * n);
    small_domain.compute_lookup_table();
    mid_domain.compute_lookup_table();
    large_domain.compute_lookup_table();

    // C++ applies coset_fft with the small/mid/large domain respectively,
    // but only the first domain.size elements are actually used.
    polynomial_arithmetic::coset_fft(&mut poly_a[..small_domain.size], &small_domain);
    polynomial_arithmetic::coset_fft(&mut poly_b[..mid_domain.size], &mid_domain);
    polynomial_arithmetic::coset_fft(&mut poly_c[..large_domain.size], &large_domain);

    // Sum at matching evaluation points:
    // small domain has n points, mid has 2n, large has 4n
    // We sum poly_a[i] + poly_b[2*i] + poly_c[4*i] for i in 0..n
    for i in 0..n {
        poly_a[i] = poly_a[i] + poly_c[4 * i];
        poly_a[i] = poly_a[i] + poly_b[2 * i];
    }

    polynomial_arithmetic::coset_ifft(&mut poly_a[..small_domain.size], &small_domain);

    for i in 0..n {
        assert_eq!(poly_a[i], expected[i], "cross-consistency mismatch at index {}", i);
    }
}

#[test]
fn barycentric_weight_evaluations() {
    // C++ test: barycentric_weight_evaluations
    let n = 16usize;

    let mut domain: EvaluationDomain<Bn254FrParams> = EvaluationDomain::new(n);

    let mut poly: Vec<Fr> = (0..n / 2).map(|_| Fr::random_element()).collect();
    let barycentric_poly = poly.clone();

    // Pad to n with zeros
    poly.resize(n, Fr::zero());

    // evaluation_point = 2 (in Montgomery form)
    let evaluation_point = Fr::from(2u64);

    // Barycentric evaluation using only the first n/2 evaluations
    let result = polynomial_arithmetic::compute_barycentric_evaluation(
        &barycentric_poly,
        &domain,
        &evaluation_point,
    );

    domain.compute_lookup_table();

    // IFFT to get coefficient form, then evaluate with Horner
    polynomial_arithmetic::ifft(&mut poly, &domain);
    let expected = polynomial_arithmetic::evaluate(&poly, &evaluation_point);

    assert_eq!(result, expected);
}

#[test]
fn linear_poly_product() {
    // C++ test: linear_poly_product
    let n = 64usize;
    let roots: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();

    let z = Fr::random_element();
    let mut expected = Fr::one();
    for i in 0..n {
        expected = expected * (z - roots[i]);
    }

    let dest = polynomial_arithmetic::compute_linear_polynomial_product(&roots);
    let result = polynomial_arithmetic::evaluate(&dest, &z);

    assert_eq!(result, expected);
}

#[test]
fn split_polynomial_evaluate() {
    // C++ test: split_polynomial_evaluate — chunked evaluate matches full evaluate
    let n = 256usize;
    let poly: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();

    let num_poly = 4usize;
    let poly_size = n / num_poly;
    let split_data: Vec<Vec<Fr>> = (0..num_poly)
        .map(|j| poly[j * poly_size..(j + 1) * poly_size].to_vec())
        .collect();

    let z = Fr::random_element();

    // Evaluate full polynomial
    let expected = polynomial_arithmetic::evaluate(&poly, &z);

    // Evaluate split: p(z) = p0(z) + z^{poly_size} * p1(z) + z^{2*poly_size} * p2(z) + ...
    let mut z_pow = Fr::one();
    let z_n = z.pow(&[poly_size as u64, 0, 0, 0]);
    let mut result = Fr::zero();
    for j in 0..num_poly {
        result = result + polynomial_arithmetic::evaluate(&split_data[j], &z) * z_pow;
        z_pow = z_pow * z_n;
    }

    assert_eq!(result, expected);
}

// ═══════════════════════════════════════════════════════════════════════════════
// BarycentricData tests (from barycentric.test.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn barycentric_big_domain_construction() {
    // C++ test: CompileTimeComputation — big_domain[5] == 5 for domain=2, evals=10
    let bc = BarycentricData::<Bn254FrParams>::new(2, 10);
    assert_eq!(bc.big_domain[5], Fr::from(5u64));
}

#[test]
fn barycentric_extend() {
    // C++ test: Extend — Univariate {1,2} extend_to 10 → {1..10}
    let f = Univariate::<Bn254FrParams, 2>::new([Fr::from(1u64), Fr::from(2u64)]);
    let result = f.extend_to::<10>();
    let expected = Univariate::<Bn254FrParams, 10>::new([
        Fr::from(1u64),
        Fr::from(2u64),
        Fr::from(3u64),
        Fr::from(4u64),
        Fr::from(5u64),
        Fr::from(6u64),
        Fr::from(7u64),
        Fr::from(8u64),
        Fr::from(9u64),
        Fr::from(10u64),
    ]);
    assert_eq!(result, expected);
}

#[test]
fn barycentric_self_extend() {
    // C++ test: SelfExtend — Univariate<10> with first 2 set, self_extend_from<2> → {1..10}
    let mut f = Univariate::<Bn254FrParams, 10>::new([
        Fr::from(1u64),
        Fr::from(2u64),
        Fr::zero(),
        Fr::zero(),
        Fr::zero(),
        Fr::zero(),
        Fr::zero(),
        Fr::zero(),
        Fr::zero(),
        Fr::zero(),
    ]);
    f.self_extend_from::<2>();
    let expected = Univariate::<Bn254FrParams, 10>::new([
        Fr::from(1u64),
        Fr::from(2u64),
        Fr::from(3u64),
        Fr::from(4u64),
        Fr::from(5u64),
        Fr::from(6u64),
        Fr::from(7u64),
        Fr::from(8u64),
        Fr::from(9u64),
        Fr::from(10u64),
    ]);
    assert_eq!(f, expected);
}

#[test]
fn barycentric_evaluate() {
    // C++ test: Evaluate — Univariate {1,2} evaluate at u=5 → 6
    let f = Univariate::<Bn254FrParams, 2>::new([Fr::from(1u64), Fr::from(2u64)]);
    let result = f.evaluate(Fr::from(5u64));
    assert_eq!(result, Fr::from(6u64));
}

#[test]
fn barycentric_data_2_to_3() {
    // C++ test: BarycentricData2to3 — verify tables + extension + random eval
    let bc = BarycentricData::<Bn254FrParams>::new(2, 3);

    // big_domain
    assert_eq!(bc.big_domain, vec![Fr::from(0u64), Fr::from(1u64), Fr::from(2u64)]);

    // lagrange_denominators: d_0 = 0-1 = -1, d_1 = 1-0 = 1
    assert_eq!(bc.lagrange_denominators, vec![-Fr::one(), Fr::one()]);

    // full_numerator_values: M(0) = 0*(0-1) = 0, M(1) = 1*(1-1) = 0, M(2) = 2*(2-1) = 2
    assert_eq!(
        bc.full_numerator_values,
        vec![Fr::zero(), Fr::zero(), Fr::from(2u64)]
    );

    // e1(X) = 1 + X. evaluate at random u → u + 1
    let e1 = Univariate::<Bn254FrParams, 2>::new([Fr::from(1u64), Fr::from(2u64)]);
    let u = Fr::random_element();
    let calculated = e1.evaluate(u);
    assert_eq!(u + Fr::one(), calculated);

    // Extension from 2 to 3: {1, 2} → {1, 2, 3}
    let ext = e1.extend_to::<3>();
    let expected = Univariate::<Bn254FrParams, 3>::new([Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]);
    assert_eq!(ext, expected);
}

#[test]
fn barycentric_data_5_to_6() {
    // C++ test: BarycentricData5to6 — extend degree-4 poly {1,3,25,109,321} → {..,751}
    let e1 = Univariate::<Bn254FrParams, 5>::new([
        Fr::from(1u64),
        Fr::from(3u64),
        Fr::from(25u64),
        Fr::from(109u64),
        Fr::from(321u64),
    ]);
    let ext = e1.extend_to::<6>();
    let expected = Univariate::<Bn254FrParams, 6>::new([
        Fr::from(1u64),
        Fr::from(3u64),
        Fr::from(25u64),
        Fr::from(109u64),
        Fr::from(321u64),
        Fr::from(751u64),
    ]);
    assert_eq!(ext, expected);
}

// ═══════════════════════════════════════════════════════════════════════════════
// EqPolynomial tests (from eq_polynomial.test.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

// --- Test helpers ---

/// Reference eq(r,u) = prod_i ((1 - r_i)(1 - u_i) + r_i * u_i)
fn eq_manual(r: &[Fr], u: &[Fr]) -> Fr {
    assert_eq!(r.len(), u.len());
    let one = Fr::one();
    let mut acc = one;
    for i in 0..r.len() {
        let term = (one - r[i]) * (one - u[i]) + r[i] * u[i];
        acc = acc * term;
    }
    acc
}

/// Boolean vector of length d from mask (LSB → index 0)
fn bool_vec_from_mask(d: usize, mask: u64) -> Vec<Fr> {
    (0..d)
        .map(|i| Fr::from(((mask >> i) & 1) as u64))
        .collect()
}

// --- VerifierEqPolynomial tests ---

#[test]
fn verifier_eq_initialize_coeffs() {
    // C++ test: InitializeCoeffs — a_i, b_i for r = {0, 1, 2, 3}
    let r = vec![Fr::from(0u64), Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
    let eq = VerifierEqPolynomial::<Bn254FrParams>::new(r);

    assert_eq!(eq.r.len(), 4);
    assert_eq!(eq.a.len(), 4);
    assert_eq!(eq.b.len(), 4);

    // a_i = 2*r_i - 1
    assert_eq!(eq.a[0], -Fr::one());
    assert_eq!(eq.a[1], Fr::one());
    assert_eq!(eq.a[2], Fr::from(3u64));
    assert_eq!(eq.a[3], Fr::from(5u64));

    // b_i = 1 - r_i
    assert_eq!(eq.b[0], Fr::one());
    assert_eq!(eq.b[1], Fr::zero());
    assert_eq!(eq.b[2], -Fr::one());
    assert_eq!(eq.b[3], -(Fr::from(2u64)));
}

#[test]
fn verifier_eq_evaluate_matches_manual() {
    // C++ test: EvaluateMatchesManualSmall
    let r: Vec<Fr> = (0..5).map(|i| Fr::from(i as u64)).collect();
    let u: Vec<Fr> = (5..10).map(|i| Fr::from(i as u64)).collect();

    let eq = VerifierEqPolynomial::<Bn254FrParams>::new(r.clone());
    let got = eq.evaluate(&u);
    let expect = eq_manual(&r, &u);
    assert_eq!(got, expect);
}

#[test]
fn verifier_eq_static_eval_matches_member() {
    // C++ test: StaticEvalMatchesMemberEvaluate
    let r = vec![Fr::from(2u64), Fr::from(0u64), Fr::from(5u64), Fr::from(1u64)];
    let u = vec![Fr::from(3u64), Fr::from(7u64), Fr::from(4u64), Fr::from(6u64)];

    let s = VerifierEqPolynomial::<Bn254FrParams>::eval(&r, &u);
    let eq = VerifierEqPolynomial::<Bn254FrParams>::new(r);
    let m = eq.evaluate(&u);
    assert_eq!(s, m);
}

#[test]
fn verifier_eq_symmetry() {
    // C++ test: SymmetryEqRUEqualsEqUR
    let r: Vec<Fr> = vec![Fr::from(0u64), Fr::from(2u64), Fr::from(4u64), Fr::from(6u64), Fr::from(8u64)];
    let u: Vec<Fr> = vec![Fr::from(1u64), Fr::from(3u64), Fr::from(5u64), Fr::from(7u64), Fr::from(9u64)];

    let eq_r = VerifierEqPolynomial::<Bn254FrParams>::new(r.clone());
    let eq_u = VerifierEqPolynomial::<Bn254FrParams>::new(u.clone());

    let ru = eq_r.evaluate(&u);
    let ur = eq_u.evaluate(&r);
    assert_eq!(ru, ur);
}

#[test]
fn verifier_eq_boolean_delta() {
    // C++ test: BooleanDeltaBehavior — Kronecker delta on {0,1}^5
    let d = 5usize;

    for big_r in 0..(1u64 << d) {
        let r = bool_vec_from_mask(d, big_r);
        let eq = VerifierEqPolynomial::<Bn254FrParams>::new(r);
        for big_u in 0..(1u64 << d) {
            let u = bool_vec_from_mask(d, big_u);
            let val = eq.evaluate(&u);
            if big_r == big_u {
                assert_eq!(val, Fr::one(), "R={} U={}", big_r, big_u);
            } else {
                assert_eq!(val, Fr::zero(), "R={} U={}", big_r, big_u);
            }
        }
    }
}

#[test]
fn verifier_eq_edge_cases() {
    // C++ test: EdgeCases — empty, d=1, all zeros, all ones, alternating

    // d = 0: empty product = 1
    {
        let r: Vec<Fr> = vec![];
        let u: Vec<Fr> = vec![];
        let val = VerifierEqPolynomial::<Bn254FrParams>::eval(&r, &u);
        assert_eq!(val, Fr::one());
    }

    // d = 1: explicit formula check
    {
        let r = vec![Fr::from(2u64)];
        let u = vec![Fr::from(7u64)];
        let expect = (Fr::one() - r[0]) * (Fr::one() - u[0]) + r[0] * u[0];
        let eq = VerifierEqPolynomial::<Bn254FrParams>::new(r);
        let got = eq.evaluate(&u);
        assert_eq!(got, expect);
    }

    // all zeros
    {
        let r = vec![Fr::zero(); 8];
        let u = vec![Fr::zero(); 8];
        let eq = VerifierEqPolynomial::<Bn254FrParams>::new(r);
        assert_eq!(eq.evaluate(&u), Fr::one());
    }

    // all ones
    {
        let r = vec![Fr::one(); 8];
        let u = vec![Fr::one(); 8];
        let eq = VerifierEqPolynomial::<Bn254FrParams>::new(r);
        assert_eq!(eq.evaluate(&u), Fr::one());
    }

    // alternating Boolean pattern
    {
        let r = vec![
            Fr::zero(), Fr::one(), Fr::zero(), Fr::one(),
            Fr::zero(), Fr::one(), Fr::zero(), Fr::one(),
        ];
        let u = vec![
            Fr::one(), Fr::zero(), Fr::one(), Fr::zero(),
            Fr::one(), Fr::zero(), Fr::one(), Fr::zero(),
        ];
        let eq = VerifierEqPolynomial::<Bn254FrParams>::new(r);
        assert_eq!(eq.evaluate(&u), Fr::zero());
    }
}

// --- Prover/Verifier consistency tests ---

#[test]
fn prover_eq_matches_verifier_on_boolean() {
    // C++ test: ProverTableMatchesVerifierOnBooleanPoints — d=5, all 32 Boolean points match
    let d = 5usize;
    let r: Vec<Fr> = (0..d).map(|i| Fr::from((2 * i + 7) as u64)).collect();

    let v = VerifierEqPolynomial::<Bn254FrParams>::new(r.clone());
    let peq = ProverEqPolynomial::construct(&r, d);

    for ell in 0..(1u64 << d) {
        let u = bool_vec_from_mask(d, ell);
        let got_ver = v.evaluate(&u);
        let got_prov = peq[ell as usize];
        assert_eq!(got_prov, got_ver, "ell={}", ell);
    }
}

#[test]
fn verifier_vs_prover_arbitrary_u() {
    // C++ test: VerifierVsProverForArbitraryU — transform correctness
    let d = 5usize;
    let r: Vec<Fr> = (0..d).map(|i| Fr::from((13 + i) as u64)).collect();
    let u: Vec<Fr> = (0..d).map(|i| Fr::from((17 + 2 * i) as u64)).collect();

    let v = VerifierEqPolynomial::<Bn254FrParams>::new(r.clone());
    let ver_val = v.evaluate(&u);

    // Prover-side normalized evaluation
    let one = Fr::one();
    let mut c = one;
    for i in 0..d {
        c = c * (one - r[i]);
    }

    // gamma_i = r_i / (1 - r_i)
    let gamma: Vec<Fr> = r.iter().map(|&ri| ri * (one - ri).invert()).collect();

    let mut prov_val = one;
    for i in 0..d {
        prov_val = prov_val * (one + u[i] * (gamma[i] - one));
    }
    prov_val = c * prov_val;

    assert_eq!(ver_val, prov_val);
}

#[test]
fn eq_partial_evaluation_consistency() {
    // C++ test: PartialEvaluationConsistency — d=21, sumcheck-style incremental eval
    let d = 21usize;
    let r: Vec<Fr> = (0..d).map(|_| Fr::random_element()).collect();
    let u: Vec<Fr> = (0..d).map(|_| Fr::random_element()).collect();
    let mut u_part = vec![Fr::zero(); d];

    let mut current_element = VerifierEqPolynomial::<Bn254FrParams>::eval(&r, &u_part);
    let _pol = ProverEqPolynomial::construct(&r, d);

    for i in 0..d {
        u_part[i] = Fr::one();
        let new_element = VerifierEqPolynomial::<Bn254FrParams>::eval(&r, &u_part);
        current_element = current_element + u[i] * (new_element - current_element);
        u_part[i] = u[i];
        assert_eq!(
            current_element,
            VerifierEqPolynomial::<Bn254FrParams>::eval(&r, &u_part)
        );
    }
    assert_eq!(
        current_element,
        VerifierEqPolynomial::<Bn254FrParams>::eval(&r, &u)
    );
}

#[test]
fn compute_subset_products_powers() {
    // C++ test: GateSeparatorBetaProductsOnPowers — {2,4,16} → {1,2,4,8,16,32,64,128}
    let betas = vec![Fr::from(2u64), Fr::from(4u64), Fr::from(16u64)];
    let result = ProverEqPolynomial::compute_subset_products(&betas, 3, Fr::one());

    let expected: Vec<Fr> = vec![
        Fr::from(1u64),
        Fr::from(2u64),
        Fr::from(4u64),
        Fr::from(8u64),
        Fr::from(16u64),
        Fr::from(32u64),
        Fr::from(64u64),
        Fr::from(128u64),
    ];
    assert_eq!(result, expected);
}

#[test]
fn prover_eq_all_ones() {
    // C++ test: ProverEqAllChallengesAreOnes — d=6, only mask=63 is nonzero
    let d = 6usize;
    let n = 1usize << d;
    let r = vec![Fr::one(); d];

    let coeffs = ProverEqPolynomial::construct(&r, d);
    assert_eq!(coeffs.len(), n);

    let all_ones_mask = n - 1;
    for m in 0..n {
        let got = coeffs[m];
        let expect = if m == all_ones_mask { Fr::one() } else { Fr::zero() };
        assert_eq!(got, expect, "mask={}", m);
    }
}

#[test]
fn prover_eq_some_ones() {
    // C++ test: ProverEqSomeChallengesAreOnes — d=5, r={7,1,9,1,11}, forced bits property
    let d = 5usize;
    let n = 1usize << d;
    let r = vec![
        Fr::from(7u64),
        Fr::one(),
        Fr::from(9u64),
        Fr::one(),
        Fr::from(11u64),
    ];
    let forced: Vec<usize> = vec![1, 3]; // indices where r_i == 1

    let coeffs = ProverEqPolynomial::construct(&r, d);
    assert_eq!(coeffs.len(), n);

    let verifier = VerifierEqPolynomial::<Bn254FrParams>::new(r);

    for mask in 0..(n as u64) {
        let u = bool_vec_from_mask(d, mask);
        let verifier_val = verifier.evaluate(&u);

        let has_all_forced = forced.iter().all(|&bit| ((mask >> bit) & 1) != 0);
        let table_val = coeffs[mask as usize];

        if !has_all_forced {
            assert_eq!(table_val, Fr::zero(), "mask missing forced bits, mask={}", mask);
        } else {
            assert_eq!(table_val, verifier_val, "mask={}", mask);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// GateSeparatorPolynomial tests (from gate_separator.test.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn gate_separator_full_pow_consistency() {
    // C++ test: FullPowConsistency — d=5 random betas/variables; after each
    // partially_evaluate, check partial_evaluation_result == prod_k (1-u_k + u_k*beta_k)
    let d = 5;
    let betas: Vec<Fr> = (0..d).map(|_| Fr::random_element()).collect();
    let variables: Vec<Fr> = (0..d).map(|_| Fr::random_element()).collect();

    let mut poly = GateSeparatorPolynomial::<Bn254FrParams>::new_verifier(betas.clone());

    let mut expected_eval = Fr::one();
    for i in 0..d {
        poly.partially_evaluate(variables[i]);
        expected_eval = expected_eval * (Fr::one() - variables[i] + variables[i] * betas[i]);
        assert_eq!(poly.partial_evaluation_result, expected_eval);
    }
}

#[test]
fn gate_separator_on_powers() {
    // C++ test: GateSeparatorPolynomialsOnPowers
    // betas=[2,4,16], log=3 → beta_products==[1,2,4,8,16,32,64,128]
    let betas = vec![Fr::from(2u64), Fr::from(4u64), Fr::from(16u64)];
    let poly = GateSeparatorPolynomial::<Bn254FrParams>::new(betas, 3);
    let expected: Vec<Fr> = vec![
        Fr::from(1u64),
        Fr::from(2u64),
        Fr::from(4u64),
        Fr::from(8u64),
        Fr::from(16u64),
        Fr::from(32u64),
        Fr::from(64u64),
        Fr::from(128u64),
    ];
    assert_eq!(poly.beta_products, expected);
}

#[test]
fn gate_separator_random_betas() {
    // C++ test: GateSeparatorPolynomialsOnPowersWithDifferentBeta
    // d=5 random betas; verify each beta_products[i] == product of betas[j] for set bits j in i
    let d = 5;
    let betas: Vec<Fr> = (0..d).map(|_| Fr::random_element()).collect();

    let mut expected_products = vec![Fr::zero(); 1 << d];
    for i in 0..(1usize << d) {
        expected_products[i] = Fr::one();
        for j in 0..d {
            if (i & (1 << j)) != 0 {
                expected_products[i] = expected_products[i] * betas[j];
            }
        }
    }
    let poly = GateSeparatorPolynomial::<Bn254FrParams>::new(betas, d);
    assert_eq!(poly.beta_products, expected_products);
}

#[test]
fn gate_separator_empty_betas() {
    // Empty betas: current_element()==1, partially_evaluate is no-op, result stays 1
    let poly = GateSeparatorPolynomial::<Bn254FrParams>::new(vec![], 0);
    assert_eq!(poly.current_element(), Fr::one());
    assert_eq!(poly.partial_evaluation_result, Fr::one());

    let mut poly2 = GateSeparatorPolynomial::<Bn254FrParams>::new_verifier(vec![]);
    poly2.partially_evaluate(Fr::from(42u64));
    assert_eq!(poly2.partial_evaluation_result, Fr::one());
    assert_eq!(poly2.current_element(), Fr::one());
}

#[test]
fn gate_separator_post_challenge() {
    // d=3 random; new_with_challenges result == manual product
    let d = 3;
    let betas: Vec<Fr> = (0..d).map(|_| Fr::random_element()).collect();
    let challenges: Vec<Fr> = (0..d).map(|_| Fr::random_element()).collect();

    let poly = GateSeparatorPolynomial::<Bn254FrParams>::new_with_challenges(
        betas.clone(),
        &challenges,
    );

    let mut expected = Fr::one();
    for i in 0..d {
        expected = expected * (Fr::one() - challenges[i] + challenges[i] * betas[i]);
    }
    assert_eq!(poly.partial_evaluation_result, expected);
}

// ═══════════════════════════════════════════════════════════════════════════════
// UnivariateCoefficientBasis tests (from univariate_coefficient_basis.test.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn coeff_basis_conversion() {
    // C++ test: Conversion — Random eval→coeff→eval roundtrip for degree-1
    let a0 = Fr::random_element();
    let a1 = Fr::random_element();
    let expected = Univariate::<Bn254FrParams, 2>::new([a0, a1]);
    let coeff = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&expected);
    let result: Univariate<Bn254FrParams, 2> = coeff.to_univariate();
    assert_eq!(result, expected);
}

#[test]
fn coeff_basis_addition() {
    // C++ test: Addition — {1,2}+{3,4}={4,6}, matches Univariate addition after conversion
    let f1 = Univariate::<Bn254FrParams, 2>::new([Fr::from(1u64), Fr::from(2u64)]);
    let f2 = Univariate::<Bn254FrParams, 2>::new([Fr::from(3u64), Fr::from(4u64)]);
    let f1_m = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&f1);
    let f2_m = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&f2);

    let sum = f1_m + f2_m;
    let result: Univariate<Bn254FrParams, 2> = sum.to_univariate();
    let expected = f1 + f2;
    assert_eq!(result, expected);
}

#[test]
fn coeff_basis_multiplication() {
    // C++ test: Multiplication — (1+X)*(3+X) → degree-2, compare with extend_to*extend_to
    let f1 = Univariate::<Bn254FrParams, 2>::new([Fr::from(1u64), Fr::from(2u64)]);
    let f2 = Univariate::<Bn254FrParams, 2>::new([Fr::from(3u64), Fr::from(4u64)]);
    let f1_m = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&f1);
    let f2_m = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&f2);

    let product = f1_m.mul(&f2_m);
    let result: Univariate<Bn254FrParams, 3> = product.to_univariate();
    let expected = f1.extend_to::<3>() * f2.extend_to::<3>();
    assert_eq!(result, expected);
}

#[test]
fn coeff_basis_sqr() {
    // Random f; sqr() == f*f
    let f = Univariate::<Bn254FrParams, 2>::new([Fr::random_element(), Fr::random_element()]);
    let f_m = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&f);

    let sqr_result: Univariate<Bn254FrParams, 3> = f_m.sqr().to_univariate();
    let mul_result: Univariate<Bn254FrParams, 3> = f_m.mul(&f_m).to_univariate();
    assert_eq!(sqr_result, mul_result);

    // Also verify against pointwise squaring in evaluation form
    let ext = f.extend_to::<3>();
    let expected = ext * ext;
    assert_eq!(sqr_result, expected);
}

#[test]
fn coeff_basis_scalar_ops() {
    // +Fr, -Fr, *Fr on coefficients
    let f = Univariate::<Bn254FrParams, 2>::new([Fr::from(5u64), Fr::from(11u64)]);
    let f_m = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&f);
    let s = Fr::from(3u64);

    // add_scalar: only affects constant term
    let added = f_m.add_scalar(s);
    let result_add: Univariate<Bn254FrParams, 2> = added.to_univariate();
    // a0=5, a1=11-5=6 → (a0+3)=8, a1=6 → evals: [8, 14]
    assert_eq!(result_add.value_at(0), Fr::from(8u64));
    assert_eq!(result_add.value_at(1), Fr::from(14u64));

    // sub_scalar: only affects constant term
    let subbed = f_m.sub_scalar(s);
    let result_sub: Univariate<Bn254FrParams, 2> = subbed.to_univariate();
    assert_eq!(result_sub.value_at(0), Fr::from(2u64));
    assert_eq!(result_sub.value_at(1), Fr::from(8u64));

    // mul_scalar: all coefficients
    let scaled = f_m.mul_scalar(s);
    let result_mul: Univariate<Bn254FrParams, 2> = scaled.to_univariate();
    assert_eq!(result_mul.value_at(0), Fr::from(15u64));
    assert_eq!(result_mul.value_at(1), Fr::from(33u64));

    // Verify mul after add: cache invalidation correctness.
    // In C++, add changes the type to has_a0_plus_a1=false, so mul computes inline.
    // In Rust, we always compute inline in mul/sqr to avoid stale cache issues.
    let g1 = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&Univariate::<Bn254FrParams, 2>::new([Fr::from(1u64), Fr::from(2u64)]));
    let g2 = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&Univariate::<Bn254FrParams, 2>::new([Fr::from(3u64), Fr::from(4u64)]));
    let g3 = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&Univariate::<Bn254FrParams, 2>::new([Fr::from(5u64), Fr::from(6u64)]));
    let sum = g1 + g2; // cache in coefficients[2] is now stale
    let product: Univariate<Bn254FrParams, 3> = sum.mul(&g3).to_univariate();
    // sum = {1+3=4, 2+4=6} as evals → a0=4, a1=2 as coefficients
    // g3 = {5, 6} as evals → a0=5, a1=1 as coefficients
    // (4 + 2X)(5 + X) = 20 + 4X + 10X + 2X^2 = 20 + 14X + 2X^2
    // evals: f(0)=20, f(1)=36, f(2)=60
    let sum_ext = (Univariate::<Bn254FrParams, 2>::new([Fr::from(1u64), Fr::from(2u64)])
        + Univariate::<Bn254FrParams, 2>::new([Fr::from(3u64), Fr::from(4u64)])).extend_to::<3>();
    let g3_ext = Univariate::<Bn254FrParams, 2>::new([Fr::from(5u64), Fr::from(6u64)]).extend_to::<3>();
    assert_eq!(product, sum_ext * g3_ext);
}

// ═══════════════════════════════════════════════════════════════════════════════
// RowDisablingPolynomial tests (from row_disabling_polynomial.test.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn row_disabling_default_construction() {
    let rdp = RowDisablingPolynomial::<Bn254FrParams>::default();
    assert_eq!(rdp.eval_at_0, Fr::one());
    assert_eq!(rdp.eval_at_1, Fr::one());
}

#[test]
fn row_disabling_update_evaluations() {
    // C++ test: ComputeDisabledContribution (unit part)
    // Round 0: no change. Round 1: eval_at_0=0. Round 2: eval_at_1 *= challenge.
    let c0 = Fr::random_element();
    let c1 = Fr::random_element();
    let c2 = Fr::random_element();

    let mut rdp = RowDisablingPolynomial::<Bn254FrParams>::default();

    // Round 0: both stay 1
    rdp.update_evaluations(c0, 0);
    assert_eq!(rdp.eval_at_0, Fr::one());
    assert_eq!(rdp.eval_at_1, Fr::one());

    // Round 1: eval_at_0 becomes 0
    rdp.update_evaluations(c1, 1);
    assert_eq!(rdp.eval_at_0, Fr::zero());
    assert_eq!(rdp.eval_at_1, Fr::one());

    // Round 2: eval_at_1 *= c2
    rdp.update_evaluations(c2, 2);
    assert_eq!(rdp.eval_at_0, Fr::zero());
    assert_eq!(rdp.eval_at_1, c2);
}

#[test]
fn row_disabling_evaluate_at_challenge() {
    // C++ test: ComputeDisabledContribution (eval part)
    // d=4 random challenges; result == 1 - challenges[2]*challenges[3]
    let d = 4;
    let challenges: Vec<Fr> = (0..d).map(|_| Fr::random_element()).collect();

    let result =
        RowDisablingPolynomial::<Bn254FrParams>::evaluate_at_challenge(&challenges, d);
    let expected = Fr::one() - challenges[2] * challenges[3];
    assert_eq!(result, expected);

    // Edge: d=2 → product of empty range = 1, so result = 1 - 1 = 0
    let challenges2: Vec<Fr> = (0..2).map(|_| Fr::random_element()).collect();
    let result2 =
        RowDisablingPolynomial::<Bn254FrParams>::evaluate_at_challenge(&challenges2, 2);
    assert_eq!(result2, Fr::zero());
}

// ═══════════════════════════════════════════════════════════════════════════════
// Serialization tests
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn univariate_serialization() {
    // C++ test: UnivariateTest::Serialization — round-trip through buffer
    let original = Univariate::<Bn254FrParams, 4>::random();
    let buffer = original.to_buffer();
    assert_eq!(buffer.len(), 4 * 32);
    let deserialized = Univariate::<Bn254FrParams, 4>::serialize_from_buffer(&buffer);
    assert_eq!(deserialized, original);
}

#[test]
fn coeff_basis_serialization() {
    // C++ test: UnivariateCoefficientBasisTest::Serialization — round-trip through buffer
    let original = Univariate::<Bn254FrParams, 2>::random();
    let coeff = UnivariateCoefficientBasis::<Bn254FrParams, 2>::from_univariate(&original);
    let buffer = coeff.to_buffer();
    assert_eq!(buffer.len(), 3 * 32);
    let deserialized = UnivariateCoefficientBasis::<Bn254FrParams, 2>::serialize_from_buffer(&buffer);
    let result: Univariate<Bn254FrParams, 2> = deserialized.to_univariate();
    assert_eq!(result, original);
}
