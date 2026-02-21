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
