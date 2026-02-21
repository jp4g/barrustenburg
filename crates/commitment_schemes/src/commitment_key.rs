//! Commitment key for polynomial commitment schemes.
//!
//! Port of C++ `commitment_schemes/commitment_key.hpp`.

use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;
use bbrs_ecc::scalar_multiplication::pippenger_msm;
use bbrs_numeric::bitop::round_up_power_2;
use bbrs_polynomials::polynomial::Polynomial;
use bbrs_polynomials::polynomial_span::PolynomialSpan;

use bbrs_ecc::curves::bn254::Bn254G1Params;
use bbrs_ecc::curves::grumpkin::GrumpkinG1Params;

/// CommitmentKey object over a pairing group G1.
///
/// Commitments are computed as C = [p(x)] = sum_i a_i * G_i where G_i is the
/// i-th element of the SRS. For BN254, the SRS is { [x^j]_1 } where x is
/// unknown. For Grumpkin, they are random points.
///
/// Port of C++ `CommitmentKey<Curve>`.
pub struct CommitmentKey<C: CurveParams> {
    /// SRS monomial points.
    srs_points: Vec<AffineElement<C>>,
    /// Dyadic (power-of-2) size of the SRS.
    pub dyadic_size: usize,
}

impl<C: CurveParams> CommitmentKey<C> {
    /// Construct from raw SRS points.
    pub fn from_points(points: Vec<AffineElement<C>>) -> Self {
        let dyadic_size = round_up_power_2(points.len() as u64) as usize;
        Self {
            srs_points: points,
            dyadic_size,
        }
    }

    /// Whether the commitment key is properly initialized.
    pub fn initialized(&self) -> bool {
        !self.srs_points.is_empty()
    }

    /// Get the SRS monomial points.
    pub fn srs_points(&self) -> &[AffineElement<C>] {
        &self.srs_points
    }

    /// Commit to a polynomial p(X) = sum_i a_i * X^i.
    ///
    /// Returns C = [p(x)] = sum_i a_i * G_i.
    pub fn commit(&self, polynomial: &Polynomial<C::ScalarFieldParams>) -> AffineElement<C> {
        let span = PolynomialSpan::new(polynomial.data(), polynomial.start_index());
        self.commit_span(&span)
    }

    /// Commit to a polynomial given as a `PolynomialSpan`.
    pub fn commit_span(
        &self,
        polynomial: &PolynomialSpan<C::ScalarFieldParams>,
    ) -> AffineElement<C> {
        let consumed_srs = polynomial.start_index + polynomial.size();
        assert!(
            consumed_srs <= self.srs_points.len(),
            "Attempting to commit to a polynomial that needs {} points with an SRS of size {}",
            consumed_srs,
            self.srs_points.len()
        );

        let scalars = polynomial.span;
        let points =
            &self.srs_points[polynomial.start_index..polynomial.start_index + scalars.len()];
        pippenger_msm(scalars, points).to_affine()
    }

    /// Batch-commit to multiple polynomials.
    ///
    /// Uses separate MSM calls for each polynomial. Returns a commitment
    /// for each input polynomial.
    pub fn batch_commit(&self, polynomials: &[&Polynomial<C::ScalarFieldParams>]) -> Vec<AffineElement<C>> {
        polynomials
            .iter()
            .map(|poly| self.commit(poly))
            .collect()
    }
}

// ── Curve-specific constructors from global CRS ───────────────────────────────

impl CommitmentKey<Bn254G1Params> {
    /// Construct a BN254 commitment key from the global CRS.
    pub fn new(num_points: usize) -> Self {
        let dyadic = round_up_power_2(num_points as u64) as usize;
        let factory = bbrs_srs::global_crs::get_bn254_crs_factory();
        let crs = factory.get_crs(dyadic);
        let points = crs.get_monomial_points().to_vec();
        Self {
            srs_points: points,
            dyadic_size: dyadic,
        }
    }
}

impl CommitmentKey<GrumpkinG1Params> {
    /// Construct a Grumpkin commitment key from the global CRS.
    pub fn new(num_points: usize) -> Self {
        let dyadic = round_up_power_2(num_points as u64) as usize;
        let factory = bbrs_srs::global_crs::get_grumpkin_crs_factory();
        let crs = factory.get_crs(dyadic);
        let points = crs.get_monomial_points().to_vec();
        Self {
            srs_points: points,
            dyadic_size: dyadic,
        }
    }
}
