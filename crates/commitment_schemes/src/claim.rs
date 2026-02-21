//! Opening claims for polynomial commitment schemes.
//!
//! Port of C++ `commitment_schemes/claim.hpp`.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;
use bbrs_polynomials::polynomial::Polynomial;

use crate::commitment_key::CommitmentKey;

/// Opening pair (r, v) for some witness polynomial p(X) such that p(r) = v.
///
/// Port of C++ `OpeningPair<Curve>`.
#[derive(Clone, Debug)]
pub struct OpeningPair<P: FieldParams> {
    /// Challenge point r.
    pub challenge: Field<P>,
    /// Evaluation v = p(r).
    pub evaluation: Field<P>,
}

impl<P: FieldParams> PartialEq for OpeningPair<P> {
    fn eq(&self, other: &Self) -> bool {
        self.challenge == other.challenge && self.evaluation == other.evaluation
    }
}

impl<P: FieldParams> Eq for OpeningPair<P> {}

/// Polynomial p and an opening pair (r, v) such that p(r) = v.
///
/// Port of C++ `ProverOpeningClaim<Curve>`.
pub struct ProverOpeningClaim<C: CurveParams> {
    /// The witness polynomial p.
    pub polynomial: Polynomial<C::ScalarFieldParams>,
    /// (challenge r, evaluation v = p(r)).
    pub opening_pair: OpeningPair<C::ScalarFieldParams>,
    /// Gemini folds must be opened at both `challenge` and `-challenge`.
    /// Instead of copying a polynomial into 2 claims, we raise this flag that
    /// turns on relevant claim processing logic in Shplonk.
    pub gemini_fold: bool,
}

/// Unverified claim (C, r, v) for some witness polynomial p(X) such that
/// C = Commit(p(X)) and p(r) = v.
///
/// Port of C++ `OpeningClaim<Curve>`.
#[derive(Clone, Debug)]
pub struct OpeningClaim<C: CurveParams> {
    /// (challenge r, evaluation v = p(r)).
    pub opening_pair: OpeningPair<C::ScalarFieldParams>,
    /// Commitment to univariate polynomial p(X).
    pub commitment: AffineElement<C>,
}

impl<C: CurveParams> PartialEq for OpeningClaim<C> {
    fn eq(&self, other: &Self) -> bool {
        self.opening_pair == other.opening_pair && self.commitment == other.commitment
    }
}

impl<C: CurveParams> Eq for OpeningClaim<C> {}

impl<C: CurveParams> OpeningClaim<C> {
    /// Inefficiently verify the claim by recomputing the commitment and
    /// evaluating the polynomial at r.
    ///
    /// Returns true iff C = Commit(p(X)) and p(r) = v.
    pub fn verify(
        &self,
        ck: &CommitmentKey<C>,
        polynomial: &Polynomial<C::ScalarFieldParams>,
    ) -> bool {
        let real_eval = polynomial.evaluate(&self.opening_pair.challenge);
        if real_eval != self.opening_pair.evaluation {
            return false;
        }
        let real_commitment = ck.commit(polynomial);
        real_commitment == self.commitment
    }
}

/// An accumulator consisting of the Shplonk evaluation challenge and vectors
/// of commitments and scalars.
///
/// This structure is used in `reduce_verify_batch_opening_claim` of KZG or IPA.
/// It always represents a zero evaluation claim.
///
/// Port of C++ `BatchOpeningClaim<Curve>`.
pub struct BatchOpeningClaim<C: CurveParams> {
    pub commitments: Vec<AffineElement<C>>,
    pub scalars: Vec<Field<C::ScalarFieldParams>>,
    pub evaluation_point: Field<C::ScalarFieldParams>,
}
