use std::sync::Arc;

use bbrs_ecc::curves::bn254::G1Affine as Bn254G1Affine;
use bbrs_ecc::curves::grumpkin::G1Affine as GrumpkinG1Affine;

// ---------------------------------------------------------------------------
// Crs<BN254> specialization
// ---------------------------------------------------------------------------

/// CRS for BN254 — mirrors C++ `Crs<curve::BN254>`.
///
/// NOTE: `get_precomputed_g2_lines` and `get_g2x` are omitted because
/// pairing / G2 types are not yet ported. They will be added when the
/// pairing crate lands.
pub trait Bn254Crs: Send + Sync {
    /// Returns the monomial points (G1 affine elements) for the pippenger algorithm.
    fn get_monomial_points(&self) -> &[Bn254G1Affine];

    /// Number of monomial points.
    fn get_monomial_size(&self) -> usize;

    /// The first G1 element from the CRS, used by the Shplonk verifier.
    fn get_g1_identity(&self) -> Bn254G1Affine;
}

// ---------------------------------------------------------------------------
// Crs<Grumpkin> specialization
// ---------------------------------------------------------------------------

/// CRS for Grumpkin — mirrors C++ `Crs<curve::Grumpkin>`.
pub trait GrumpkinCrs: Send + Sync {
    /// Returns the monomial points (G1 affine elements) for the pippenger algorithm.
    fn get_monomial_points(&self) -> &[GrumpkinG1Affine];

    /// Number of monomial points.
    fn get_monomial_size(&self) -> usize;

    /// The first G1 element from the CRS.
    fn get_g1_identity(&self) -> GrumpkinG1Affine;
}

// ---------------------------------------------------------------------------
// CrsFactory — separate traits for BN254 and Grumpkin
// ---------------------------------------------------------------------------

/// Factory that produces BN254 CRS instances.
/// Mirrors C++ `CrsFactory<curve::BN254>`.
pub trait Bn254CrsFactory: Send + Sync {
    fn get_crs(&self, degree: usize) -> Arc<dyn Bn254Crs>;
    fn get_verifier_crs(&self) -> Arc<dyn Bn254Crs> {
        self.get_crs(1)
    }
}

/// Factory that produces Grumpkin CRS instances.
/// Mirrors C++ `CrsFactory<curve::Grumpkin>`.
pub trait GrumpkinCrsFactory: Send + Sync {
    fn get_crs(&self, degree: usize) -> Arc<dyn GrumpkinCrs>;
    fn get_verifier_crs(&self) -> Arc<dyn GrumpkinCrs> {
        self.get_crs(1)
    }
}
