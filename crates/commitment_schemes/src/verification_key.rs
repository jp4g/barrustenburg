//! Verification keys for polynomial commitment schemes.
//!
//! Port of C++ `commitment_schemes/verification_key.hpp`.
//!
//! Provides specializations for BN254 (pairing-based) and Grumpkin (MSM-based).

use std::sync::Arc;

use bbrs_ecc::curves::bn254::{self, Fq12, G1Affine as Bn254G1Affine, G2AffineElement};
use bbrs_ecc::curves::bn254_pairing::reduced_ate_pairing;
use bbrs_ecc::curves::grumpkin::G1Affine as GrumpkinG1Affine;
use bbrs_ecc::groups::element::Element;
use bbrs_srs::factories::{Bn254Crs, GrumpkinCrs, GrumpkinCrsFactory};

/// Verifier commitment key specialized for BN254.
///
/// Provides pairing-based verification: e(P0, [1]_2) * e(P1, [x]_2) == 1.
///
/// Port of C++ `VerifierCommitmentKey<curve::BN254>`.
pub struct Bn254VerifierCommitmentKey {
    /// SRS from the global CRS factory.
    srs: Option<Arc<dyn Bn254Crs>>,
    /// The G2 SRS point [x]_2 for pairing verification.
    g2_x: G2AffineElement,
}

impl Bn254VerifierCommitmentKey {
    /// Construct with a specific G2 SRS point.
    pub fn with_g2x(g2_x: G2AffineElement) -> Self {
        Self { srs: None, g2_x }
    }

    /// Lazy-initialize the SRS from the global CRS factory.
    pub fn initialize(&mut self) {
        if self.srs.is_none() {
            self.srs = Some(
                bbrs_srs::global_crs::get_bn254_crs_factory().get_verifier_crs(),
            );
        }
    }

    /// Whether the SRS has been initialized.
    pub fn initialized(&self) -> bool {
        self.srs.is_some()
    }

    /// Get the G1 identity (first SRS point).
    pub fn get_g1_identity(&mut self) -> Bn254G1Affine {
        self.initialize();
        self.srs.as_ref().unwrap().get_g1_identity()
    }

    /// Verify a pairing equation over 2 points using the verifier SRS.
    ///
    /// Returns true iff e(P0, [1]_2) * e(P1, [x]_2) == 1_T.
    pub fn pairing_check(
        &mut self,
        p0: &Element<bn254::Bn254G1Params>,
        p1: &Element<bn254::Bn254G1Params>,
    ) -> bool {
        self.initialize();
        let p0_affine = p0.to_affine();
        let p1_affine = p1.to_affine();

        // e(P0, [1]_2)
        let e0 = reduced_ate_pairing(&p0_affine, &G2AffineElement::generator());
        // e(P1, [x]_2)
        let e1 = reduced_ate_pairing(&p1_affine, &self.g2_x);
        // Check product equals 1_T
        (e0 * e1) == Fq12::one()
    }
}

/// Verifier commitment key specialized for Grumpkin (IPA-based).
///
/// Provides access to the SRS monomial points for multi-scalar multiplication.
///
/// Port of C++ `VerifierCommitmentKey<curve::Grumpkin>`.
pub struct GrumpkinVerifierCommitmentKey {
    srs: Arc<dyn GrumpkinCrs>,
}

impl GrumpkinVerifierCommitmentKey {
    /// Construct from the global CRS factory with the given number of points.
    pub fn new(num_points: usize) -> Self {
        let factory = bbrs_srs::global_crs::get_grumpkin_crs_factory();
        Self {
            srs: factory.get_crs(num_points),
        }
    }

    /// Construct from a specific CRS factory.
    pub fn from_factory(
        num_points: usize,
        factory: &dyn GrumpkinCrsFactory,
    ) -> Self {
        Self {
            srs: factory.get_crs(num_points),
        }
    }

    /// Whether the SRS has been initialized.
    pub fn initialized(&self) -> bool {
        true // Always initialized after construction
    }

    /// Get the G1 identity point.
    pub fn get_g1_identity(&self) -> GrumpkinG1Affine {
        self.srs.get_g1_identity()
    }

    /// Get the SRS monomial points.
    pub fn get_monomial_points(&self) -> &[GrumpkinG1Affine] {
        self.srs.get_monomial_points()
    }
}
