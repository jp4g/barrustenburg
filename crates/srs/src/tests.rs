use bbrs_ecc::curves::bn254::{G1Affine as Bn254G1Affine, G1Element as Bn254G1Element};
use bbrs_ecc::curves::grumpkin::{G1Affine as GrumpkinG1Affine, G1Element as GrumpkinG1Element};

use crate::factories::{
    Bn254CrsFactory, GrumpkinCrsFactory, MemBn254CrsFactory, MemGrumpkinCrsFactory,
};

/// Generate `n` distinct BN254 G1 affine points: G, 2G, 3G, ...
fn bn254_test_points(n: usize) -> Vec<Bn254G1Affine> {
    let generator = Bn254G1Element::from_affine(&Bn254G1Affine::one());
    let mut acc = generator;
    let mut points = Vec::with_capacity(n);
    points.push(acc.to_affine());
    for _ in 1..n {
        acc = acc + generator;
        points.push(acc.to_affine());
    }
    points
}

/// Generate `n` distinct Grumpkin G1 affine points: G, 2G, 3G, ...
fn grumpkin_test_points(n: usize) -> Vec<GrumpkinG1Affine> {
    let generator = GrumpkinG1Element::from_affine(&GrumpkinG1Affine::one());
    let mut acc = generator;
    let mut points = Vec::with_capacity(n);
    points.push(acc.to_affine());
    for _ in 1..n {
        acc = acc + generator;
        points.push(acc.to_affine());
    }
    points
}

// ===========================================================================
// MemBn254Crs tests
// ===========================================================================

#[test]
fn bn254_mem_crs_from_known_points() {
    let points = bn254_test_points(8);
    let factory = MemBn254CrsFactory::new(&points);
    let crs = factory.get_crs(8);

    assert_eq!(crs.get_monomial_size(), 8);
    assert_eq!(crs.get_monomial_points().len(), 8);

    // Verify points match what we put in
    for (a, b) in crs.get_monomial_points().iter().zip(points.iter()) {
        assert_eq!(a, b);
    }
}

#[test]
fn bn254_mem_crs_g1_identity_is_first_point() {
    let points = bn254_test_points(4);
    let factory = MemBn254CrsFactory::new(&points);
    let crs = factory.get_crs(4);

    assert_eq!(crs.get_g1_identity(), points[0]);
    assert!(crs.get_g1_identity().on_curve());
}

#[test]
fn bn254_mem_crs_verifier_crs() {
    let points = bn254_test_points(4);
    let factory = MemBn254CrsFactory::new(&points);
    let vcrs = factory.get_verifier_crs();

    assert_eq!(vcrs.get_monomial_size(), 4);
    assert_eq!(vcrs.get_g1_identity(), points[0]);
}

#[test]
#[should_panic(expected = "prover trying to get too many points")]
fn bn254_mem_crs_panics_on_too_large_degree() {
    let points = bn254_test_points(4);
    let factory = MemBn254CrsFactory::new(&points);
    let _ = factory.get_crs(5);
}

#[test]
#[should_panic(expected = "invalid g1_identity")]
fn bn254_mem_crs_panics_on_empty_points() {
    let _ = MemBn254CrsFactory::new(&[]);
}

// ===========================================================================
// MemGrumpkinCrs tests
// ===========================================================================

#[test]
fn grumpkin_mem_crs_from_known_points() {
    let points = grumpkin_test_points(8);
    let factory = MemGrumpkinCrsFactory::new(&points);
    let crs = factory.get_crs(8);

    assert_eq!(crs.get_monomial_size(), 8);
    assert_eq!(crs.get_monomial_points().len(), 8);

    for (a, b) in crs.get_monomial_points().iter().zip(points.iter()) {
        assert_eq!(a, b);
    }
}

#[test]
fn grumpkin_mem_crs_g1_identity_is_first_point() {
    let points = grumpkin_test_points(4);
    let factory = MemGrumpkinCrsFactory::new(&points);
    let crs = factory.get_crs(4);

    assert_eq!(crs.get_g1_identity(), points[0]);
    assert!(crs.get_g1_identity().on_curve());
}

#[test]
fn grumpkin_mem_crs_verifier_crs() {
    let points = grumpkin_test_points(4);
    let factory = MemGrumpkinCrsFactory::new(&points);
    let vcrs = factory.get_verifier_crs();

    assert_eq!(vcrs.get_monomial_size(), 4);
    assert_eq!(vcrs.get_g1_identity(), points[0]);
}

#[test]
#[should_panic(expected = "prover trying to get too many points")]
fn grumpkin_mem_crs_panics_on_too_large_degree() {
    let points = grumpkin_test_points(4);
    let factory = MemGrumpkinCrsFactory::new(&points);
    let _ = factory.get_crs(5);
}

#[test]
#[should_panic(expected = "invalid vector")]
fn grumpkin_mem_crs_panics_on_empty_points() {
    let _ = MemGrumpkinCrsFactory::new(&[]);
}

// ===========================================================================
// Global CRS tests
// ===========================================================================

#[test]
fn global_bn254_crs_init_and_get() {
    // Note: OnceLock can only be set once per process, so this test
    // verifies the init + get path works without panicking.
    let points = bn254_test_points(4);
    crate::global_crs::init_bn254_mem_crs_factory(&points);

    let factory = crate::global_crs::get_bn254_crs_factory();
    let crs = factory.get_crs(4);
    assert_eq!(crs.get_monomial_size(), 4);
    assert_eq!(crs.get_g1_identity(), points[0]);
}

#[test]
fn global_grumpkin_crs_init_and_get() {
    let points = grumpkin_test_points(4);
    crate::global_crs::init_grumpkin_mem_crs_factory(&points);

    let factory = crate::global_crs::get_grumpkin_crs_factory();
    let crs = factory.get_crs(4);
    assert_eq!(crs.get_monomial_size(), 4);
    assert_eq!(crs.get_g1_identity(), points[0]);
}
