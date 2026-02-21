use std::sync::{Arc, OnceLock};

use bbrs_ecc::curves::bn254::G1Affine as Bn254G1Affine;
use bbrs_ecc::curves::grumpkin::G1Affine as GrumpkinG1Affine;

use crate::factories::{
    Bn254CrsFactory, GrumpkinCrsFactory, MemBn254CrsFactory, MemGrumpkinCrsFactory,
};

// ---------------------------------------------------------------------------
// Global singletons — mirrors C++ anonymous-namespace statics
// ---------------------------------------------------------------------------

static BN254_CRS_FACTORY: OnceLock<Arc<dyn Bn254CrsFactory>> = OnceLock::new();
static GRUMPKIN_CRS_FACTORY: OnceLock<Arc<dyn GrumpkinCrsFactory>> = OnceLock::new();

// ---------------------------------------------------------------------------
// Initialization from memory buffers
// ---------------------------------------------------------------------------

/// Initialize the global BN254 CRS factory from in-memory G1 points.
///
/// Mirrors C++ `init_bn254_mem_crs_factory`. The G2 point parameter is
/// omitted until G2/pairing types are ported.
pub fn init_bn254_mem_crs_factory(points: &[Bn254G1Affine]) {
    let factory = Arc::new(MemBn254CrsFactory::new(points));
    let _ = BN254_CRS_FACTORY.set(factory);
}

/// Initialize the global Grumpkin CRS factory from in-memory points.
///
/// Mirrors C++ `init_grumpkin_mem_crs_factory`.
pub fn init_grumpkin_mem_crs_factory(points: &[GrumpkinG1Affine]) {
    let factory = Arc::new(MemGrumpkinCrsFactory::new(points));
    let _ = GRUMPKIN_CRS_FACTORY.set(factory);
}

// ---------------------------------------------------------------------------
// Accessors
// ---------------------------------------------------------------------------

/// Get the global BN254 CRS factory.
///
/// Panics if `init_bn254_mem_crs_factory` (or a future file-backed init)
/// has not been called.
pub fn get_bn254_crs_factory() -> Arc<dyn Bn254CrsFactory> {
    BN254_CRS_FACTORY
        .get()
        .expect("You need to initialize the global CRS with a call to init_bn254_mem_crs_factory!")
        .clone()
}

/// Get the global Grumpkin CRS factory.
///
/// Panics if `init_grumpkin_mem_crs_factory` has not been called.
pub fn get_grumpkin_crs_factory() -> Arc<dyn GrumpkinCrsFactory> {
    GRUMPKIN_CRS_FACTORY
        .get()
        .expect(
            "You need to initialize the global CRS with a call to init_grumpkin_mem_crs_factory!",
        )
        .clone()
}
