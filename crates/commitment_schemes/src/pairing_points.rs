//! Pairing points for BN254 verification.
//!
//! Port of C++ `commitment_schemes/pairing_points.hpp`.

use bbrs_ecc::curves::bn254::{Bn254G1Params, Fr, G1Affine as Bn254G1Affine};
use bbrs_ecc::groups::element::Element;

use crate::verification_key::Bn254VerifierCommitmentKey;

/// Two EC points that represent the inputs to a pairing check.
///
/// The points may represent the output of a single partial verification or
/// the linear combination of multiple sets of pairing points (an accumulator).
///
/// Port of C++ `PairingPoints<Curve_>`.
pub struct PairingPoints {
    pub p0: Bn254G1Affine,
    pub p1: Bn254G1Affine,
}

impl PairingPoints {
    /// Create new pairing points initialized to the point at infinity.
    pub fn new() -> Self {
        Self {
            p0: Bn254G1Affine::infinity(),
            p1: Bn254G1Affine::infinity(),
        }
    }

    /// Create from two specific points.
    pub fn from_points(p0: Bn254G1Affine, p1: Bn254G1Affine) -> Self {
        Self { p0, p1 }
    }

    /// Create from an array of two points.
    pub fn from_array(points: [Bn254G1Affine; 2]) -> Self {
        Self {
            p0: points[0],
            p1: points[1],
        }
    }

    /// Index into the pairing points (0 or 1).
    pub fn get(&self, idx: usize) -> &Bn254G1Affine {
        assert!(idx < 2, "Index out of bounds");
        if idx == 0 {
            &self.p0
        } else {
            &self.p1
        }
    }

    /// Mutable index into the pairing points (0 or 1).
    pub fn get_mut(&mut self, idx: usize) -> &mut Bn254G1Affine {
        assert!(idx < 2, "Index out of bounds");
        if idx == 0 {
            &mut self.p0
        } else {
            &mut self.p1
        }
    }

    /// Aggregate with another set of pairing points using a random scalar.
    pub fn aggregate(&mut self, other: &PairingPoints) {
        assert!(
            !self.p0.is_point_at_infinity()
                && !self.p1.is_point_at_infinity()
                && !other.p0.is_point_at_infinity()
                && !other.p1.is_point_at_infinity(),
            "Shouldn't be aggregating with Point at infinity! The pairing points are probably uninitialized."
        );

        let aggregation_separator = Fr::random_element();

        // P0 = P0 + other.P0 * aggregation_separator
        let p0_proj = Element::<Bn254G1Params>::from_affine(&self.p0);
        let other_p0_proj = Element::<Bn254G1Params>::from_affine(&other.p0);
        let scaled_other_p0 = other_p0_proj.mul(&aggregation_separator);
        self.p0 = (p0_proj + scaled_other_p0).to_affine();

        // P1 = P1 + other.P1 * aggregation_separator
        let p1_proj = Element::<Bn254G1Params>::from_affine(&self.p1);
        let other_p1_proj = Element::<Bn254G1Params>::from_affine(&other.p1);
        let scaled_other_p1 = other_p1_proj.mul(&aggregation_separator);
        self.p1 = (p1_proj + scaled_other_p1).to_affine();
    }

    /// Perform the pairing check using a verifier commitment key.
    pub fn check(&self, vkey: &mut Bn254VerifierCommitmentKey) -> bool {
        let p0_proj = Element::<Bn254G1Params>::from_affine(&self.p0);
        let p1_proj = Element::<Bn254G1Params>::from_affine(&self.p1);
        vkey.pairing_check(&p0_proj, &p1_proj)
    }
}

impl Default for PairingPoints {
    fn default() -> Self {
        Self::new()
    }
}

impl PartialEq for PairingPoints {
    fn eq(&self, other: &Self) -> bool {
        self.p0 == other.p0 && self.p1 == other.p1
    }
}

impl Eq for PairingPoints {}
