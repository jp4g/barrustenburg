//! KZG polynomial commitment scheme for BN254.
//!
//! Port of C++ `commitment_schemes/kzg/kzg.hpp`.
//! Implements the native (non-circuit) path only.

use bbrs_ecc::curves::bn254::{Bn254G1Params, Fr, G1Affine as Bn254G1Affine};
use bbrs_ecc::groups::element::Element;
use bbrs_transcript::NativeTranscript;

use crate::batch_mul::batch_mul_native;
use crate::claim::{BatchOpeningClaim, OpeningClaim, ProverOpeningClaim};
use crate::commitment_key::CommitmentKey;
use crate::pairing_points::PairingPoints;

/// KZG polynomial commitment scheme for BN254.
///
/// Port of C++ `KZG<Curve_>` specialized for BN254 (native path only).
pub struct KZG;

impl KZG {
    /// Computes the KZG commitment to an opening proof polynomial at a single
    /// evaluation point.
    ///
    /// Given a witness polynomial p(X) and an opening pair (r, v = p(r)),
    /// computes the quotient polynomial q(X) = (p(X) - v) / (X - r),
    /// commits to it, and sends the commitment to the transcript.
    ///
    /// Port of C++ `KZG::compute_opening_proof`.
    pub fn compute_opening_proof(
        ck: &CommitmentKey<Bn254G1Params>,
        opening_claim: ProverOpeningClaim<Bn254G1Params>,
        transcript: &mut NativeTranscript,
    ) {
        let mut quotient = opening_claim.polynomial;
        let pair = opening_claim.opening_pair;

        let quotient_commitment: Bn254G1Affine;
        if quotient.is_empty() {
            // Treat the empty polynomial as the zero polynomial
            quotient_commitment = Bn254G1Affine::infinity();
        } else {
            *quotient.at_mut(0) = quotient.get(0) - pair.evaluation;
            // Compute coefficients for q(X) = (p(X) - v) / (X - r)
            quotient.factor_roots(&pair.challenge);
            quotient_commitment = ck.commit(&quotient);
        }

        transcript.send_to_verifier("KZG:W", &quotient_commitment);

        // Masking challenge is used in the recursive setting to perform batch_mul.
        // Not used by the prover directly; needed for consistent transcript state
        // with the verifier in the batch path.
        let _masking_challenge: Fr = transcript.get_challenge("KZG:masking_challenge");
    }

    /// Computes the input points for the pairing check needed to verify a KZG
    /// opening claim of a single polynomial commitment.
    ///
    /// Returns `PairingPoints { P₀, P₁ }` where:
    ///   - P₀ = C − v⋅[1]₁ + r⋅[W(x)]₁
    ///   - P₁ = −[W(x)]₁
    ///
    /// Port of C++ `KZG::reduce_verify` (native path only).
    pub fn reduce_verify(
        claim: &OpeningClaim<Bn254G1Params>,
        transcript: &mut NativeTranscript,
    ) -> PairingPoints {
        let quotient_commitment: Bn254G1Affine =
            transcript.receive_from_prover("KZG:W");

        // The pairing check: e(C + r*W - v*G, [1]_2) * e(-W, [x]_2) = 1
        let mut p0 = Element::<Bn254G1Params>::from_affine(&claim.commitment);
        p0 += Element::<Bn254G1Params>::from_affine(&quotient_commitment)
            .mul(&claim.opening_pair.challenge);
        p0 -= Element::<Bn254G1Params>::one().mul(&claim.opening_pair.evaluation);

        let p1 = -quotient_commitment;

        PairingPoints::from_points(p0.to_affine(), p1)
    }

    /// Computes the input points for the pairing check needed to verify a KZG
    /// opening claim obtained from a Shplemini accumulator.
    ///
    /// The commitment C is encoded in the `batch_opening_claim` vectors:
    /// C = Σ commitments_i ⋅ scalars_i. The quotient commitment [W]₁ is
    /// appended to the vectors and the evaluation point z is added to the
    /// scalars, then a single batch_mul computes C + W⋅z.
    ///
    /// Returns `PairingPoints { P₀, P₁ }` where:
    ///   - P₀ = C + [W(x)]₁ ⋅ z
    ///   - P₁ = −[W(x)]₁
    ///
    /// Port of C++ `KZG::reduce_verify_batch_opening_claim` (native path only).
    pub fn reduce_verify_batch_opening_claim(
        mut batch_opening_claim: BatchOpeningClaim<Bn254G1Params>,
        transcript: &mut NativeTranscript,
    ) -> PairingPoints {
        let quotient_commitment: Bn254G1Affine =
            transcript.receive_from_prover("KZG:W");

        // This challenge is used to compute offset generators in batch_mul
        // in the recursive setting (not used in native path).
        let _masking_challenge: Fr =
            transcript.get_challenge("KZG:masking_challenge");

        // Place the commitment to W into commitments
        batch_opening_claim.commitments.push(quotient_commitment);
        // Update scalars by adding the Shplonk evaluation challenge z
        batch_opening_claim.scalars.push(batch_opening_claim.evaluation_point);

        // Compute C + [W]₁ ⋅ z
        let p0 = batch_mul_native::<Bn254G1Params>(
            &batch_opening_claim.commitments,
            &batch_opening_claim.scalars,
        );
        let p1 = -quotient_commitment;

        PairingPoints::from_points(p0, p1)
    }
}
