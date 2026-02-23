//! Port of `ultra_prover.hpp`/`.cpp` — Ultra Honk Prover.
//!
//! The Ultra prover orchestrates the full proof generation:
//! 1. Oink phase: commit to witness polynomials, compute challenges
//! 2. Sumcheck phase: run sumcheck protocol
//! 3. PCS phase: run Shplemini polynomial commitment opening + KZG

use bbrs_commitment_schemes::gemini::PolynomialBatcher;
use bbrs_commitment_schemes::claim::{OpeningPair, ProverOpeningClaim};
use bbrs_commitment_schemes::kzg::KZG;
use bbrs_commitment_schemes::shplonk::shplemini::ShpleminiProver;
use bbrs_ecc::curves::bn254::{Bn254FrParams, Bn254G1Params, Fr};
use bbrs_polynomials::polynomial::Polynomial;
use bbrs_transcript::NativeTranscript;

use crate::oink_prover::OinkProver;
use crate::proving_key::ProvingKey;
use crate::ultra_sumcheck::UltraSumcheckProver;

/// Ultra Honk prover.
///
/// Port of C++ `UltraProver<UltraFlavor>`.
pub struct UltraProver {
    pub proving_key: ProvingKey,
}

impl UltraProver {
    pub fn new(proving_key: ProvingKey) -> Self {
        Self { proving_key }
    }

    /// Generate a proof.
    ///
    /// Returns the serialized proof as a vector of field elements.
    pub fn prove(mut self) -> Vec<Fr> {
        let mut transcript = NativeTranscript::new();

        // Phase 1: Oink — commit to witnesses, compute challenges
        let oink_output = {
            let oink = OinkProver::new(&mut self.proving_key);
            oink.prove(&mut transcript)
        };

        let circuit_size = self.proving_key.circuit_size;
        let log_circuit_size = self.proving_key.log_circuit_size;

        // Get gate challenges for pow-polynomial
        let gate_challenges: Vec<Fr> = (0..log_circuit_size)
            .map(|i| transcript.get_challenge(&format!("Sumcheck:gate_challenge_{}", i)))
            .collect();

        // Phase 2: Sumcheck
        let sumcheck_output = {
            let sumcheck_prover = UltraSumcheckProver::new(
                circuit_size,
                &self.proving_key.polynomials,
                oink_output.alpha,
                gate_challenges,
                oink_output.relation_parameters,
            );
            sumcheck_prover.prove(&mut transcript)
        };

        // Phase 3: PCS (Shplemini + KZG)
        // Set up the polynomial batcher with all unshifted and to-be-shifted polynomials
        let mut batcher = PolynomialBatcher::<Bn254FrParams>::new(circuit_size);

        // Add unshifted polynomials (precomputed + witness excluding shifted)
        let polys = &self.proving_key.polynomials;
        let unshifted_refs: Vec<&Polynomial<Bn254FrParams>> = polys
            .get_precomputed()
            .iter()
            .chain(polys.get_witness().iter())
            .copied()
            .collect();
        for poly in &unshifted_refs {
            batcher.unshifted.push(poly);
        }

        // Add to-be-shifted polynomials
        let to_be_shifted_refs: Vec<&Polynomial<Bn254FrParams>> =
            polys.get_to_be_shifted().to_vec();
        for poly in &to_be_shifted_refs {
            batcher.to_be_shifted_by_one.push(poly);
        }

        // Run Shplemini prover
        let empty_libra: [Polynomial<Bn254FrParams>; 4] = std::array::from_fn(|_| {
            Polynomial::new(0, 0, 0)
        });
        let empty_sumcheck_round: Vec<Polynomial<Bn254FrParams>> = Vec::new();
        let empty_sumcheck_evals: Vec<[Fr; 3]> = Vec::new();

        let opening_claim = ShpleminiProver::prove(
            circuit_size,
            &mut batcher,
            &sumcheck_output.challenge,
            &self.proving_key.commitment_key,
            &mut transcript,
            &empty_libra,
            &empty_sumcheck_round,
            &empty_sumcheck_evals,
        );

        // Convert shplonk ProverOpeningClaim<FieldParams> to claim ProverOpeningClaim<CurveParams>
        let kzg_claim = ProverOpeningClaim::<Bn254G1Params> {
            polynomial: opening_claim.polynomial,
            opening_pair: OpeningPair {
                challenge: opening_claim.opening_pair.challenge,
                evaluation: opening_claim.opening_pair.evaluation,
            },
            gemini_fold: opening_claim.gemini_fold,
        };

        // Run KZG opening proof
        KZG::compute_opening_proof(
            &self.proving_key.commitment_key,
            kzg_claim,
            &mut transcript,
        );

        transcript.export_proof()
    }
}
