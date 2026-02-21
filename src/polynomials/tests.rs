use crate::ecc::curves::bn254::Bn254FrParams;
use crate::ecc::fields::field::Field;
use crate::polynomials::evaluation_domain::EvaluationDomain;
use crate::polynomials::polynomial::Polynomial;
use crate::polynomials::polynomial_arithmetic;
use crate::polynomials::polynomial_span::PolynomialSpan;
use crate::polynomials::univariate::{Univariate, UnivariateView};

type Fr = Field<Bn254FrParams>;

// ═══════════════════════════════════════════════════════════════════════════════
// Polynomial tests (from polynomial.test.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

#[test]
fn polynomial_shifted() {
    let n = 8usize;
    let mut p: Polynomial<Bn254FrParams> = Polynomial::shiftable(n);
    for i in 0..n {
        *p.at_mut(i) = Fr::from((i + 1) as u64);
    }
    let s = p.shifted();
    assert_eq!(s.start_index(), 1);
    assert_eq!(s.size(), n - 1);
    for i in 1..n {
        assert_eq!(s.get(i), Fr::from((i + 1) as u64));
    }
    // Virtual zero outside range
    assert_eq!(s.get(0), Fr::zero());
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
