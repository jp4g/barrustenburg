use std::sync::Arc;

use bbrs_ecc::curves::bn254::G1Affine as Bn254G1Affine;

use super::crs_factory::{Bn254Crs, Bn254CrsFactory};

// ---------------------------------------------------------------------------
// MemBn254Crs — in-memory CRS backed by a Vec of G1 affine points
// ---------------------------------------------------------------------------

/// In-memory CRS for BN254. Mirrors C++ `MemBn254Crs`.
///
/// Stores the monomial points as a contiguous `Vec<Bn254G1Affine>`.
/// G2 / pairing precomputation is deferred until the pairing crate lands.
pub struct MemBn254Crs {
    monomials: Vec<Bn254G1Affine>,
}

impl MemBn254Crs {
    pub fn new(points: &[Bn254G1Affine]) -> Self {
        assert!(
            !points.is_empty() && points[0].on_curve(),
            "invalid g1_identity passed to MemBn254Crs"
        );
        Self {
            monomials: points.to_vec(),
        }
    }
}

impl Bn254Crs for MemBn254Crs {
    fn get_monomial_points(&self) -> &[Bn254G1Affine] {
        &self.monomials
    }

    fn get_monomial_size(&self) -> usize {
        self.monomials.len()
    }

    fn get_g1_identity(&self) -> Bn254G1Affine {
        self.monomials[0]
    }
}

// ---------------------------------------------------------------------------
// MemBn254CrsFactory
// ---------------------------------------------------------------------------

/// Factory that creates an in-memory BN254 CRS from a pre-supplied set of points.
/// Mirrors C++ `MemBn254CrsFactory`.
pub struct MemBn254CrsFactory {
    crs: Arc<MemBn254Crs>,
}

impl MemBn254CrsFactory {
    pub fn new(points: &[Bn254G1Affine]) -> Self {
        let crs = Arc::new(MemBn254Crs::new(points));
        Self { crs }
    }
}

impl Bn254CrsFactory for MemBn254CrsFactory {
    fn get_crs(&self, degree: usize) -> Arc<dyn Bn254Crs> {
        assert!(
            self.crs.get_monomial_size() >= degree,
            "prover trying to get too many points in MemBn254CrsFactory! {} vs {}",
            self.crs.get_monomial_size(),
            degree
        );
        self.crs.clone()
    }
}
