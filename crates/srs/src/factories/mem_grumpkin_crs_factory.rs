use std::sync::Arc;

use bbrs_ecc::curves::grumpkin::G1Affine as GrumpkinG1Affine;

use super::crs_factory::{GrumpkinCrs, GrumpkinCrsFactory};

// ---------------------------------------------------------------------------
// MemGrumpkinCrs — in-memory CRS backed by a Vec of Grumpkin G1 affine points
// ---------------------------------------------------------------------------

/// In-memory CRS for Grumpkin. Mirrors C++ `MemGrumpkinCrs`.
pub struct MemGrumpkinCrs {
    monomials: Vec<GrumpkinG1Affine>,
}

impl MemGrumpkinCrs {
    pub fn new(points: &[GrumpkinG1Affine]) -> Self {
        assert!(
            !points.is_empty() && points[0].on_curve(),
            "invalid vector passed to MemGrumpkinCrs"
        );
        Self {
            monomials: points.to_vec(),
        }
    }
}

impl GrumpkinCrs for MemGrumpkinCrs {
    fn get_monomial_points(&self) -> &[GrumpkinG1Affine] {
        &self.monomials
    }

    fn get_monomial_size(&self) -> usize {
        self.monomials.len()
    }

    fn get_g1_identity(&self) -> GrumpkinG1Affine {
        self.monomials[0]
    }
}

// ---------------------------------------------------------------------------
// MemGrumpkinCrsFactory
// ---------------------------------------------------------------------------

/// Factory that creates an in-memory Grumpkin CRS from a pre-supplied set of points.
/// Mirrors C++ `MemGrumpkinCrsFactory`.
pub struct MemGrumpkinCrsFactory {
    crs: Arc<MemGrumpkinCrs>,
}

impl MemGrumpkinCrsFactory {
    pub fn new(points: &[GrumpkinG1Affine]) -> Self {
        let crs = Arc::new(MemGrumpkinCrs::new(points));
        Self { crs }
    }
}

impl GrumpkinCrsFactory for MemGrumpkinCrsFactory {
    fn get_crs(&self, degree: usize) -> Arc<dyn GrumpkinCrs> {
        assert!(
            self.crs.get_monomial_size() >= degree,
            "prover trying to get too many points in MemGrumpkinCrsFactory - {} is more than {}",
            degree,
            self.crs.get_monomial_size()
        );
        self.crs.clone()
    }
}
