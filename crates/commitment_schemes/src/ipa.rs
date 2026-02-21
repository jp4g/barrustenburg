//! IPA (Inner Product Argument) commitment scheme.
//!
//! Port of C++ `commitment_schemes/ipa/ipa.hpp`.
//!
//! Uses the optimized version that only multiplies half of the elements
//! of each vector in each prover round. Works with Grumpkin curve (no pairings),
//! using SRS monomial points.

use bbrs_ecc::curves::bn254::{Bn254FqParams, Bn254FrParams};
use bbrs_ecc::curves::grumpkin::{GrumpkinFr, GrumpkinG1Params};
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::element::Element;
use bbrs_ecc::scalar_multiplication::pippenger_msm;
use bbrs_transcript::NativeTranscript;

use crate::batch_mul::batch_mul_native;
use crate::claim::{BatchOpeningClaim, OpeningClaim, OpeningPair, ProverOpeningClaim};
use crate::commitment_key::CommitmentKey;
use crate::verification_key::GrumpkinVerifierCommitmentKey;

type Fr = GrumpkinFr;
type Commitment = AffineElement<GrumpkinG1Params>;
type GroupElement = Element<GrumpkinG1Params>;
type CK = CommitmentKey<GrumpkinG1Params>;
type VK = GrumpkinVerifierCommitmentKey;

/// Convert a BN254 Fr transcript challenge to GrumpkinFr (BN254 Fq).
///
/// Transcript challenges are 127-bit values derived from Poseidon2 (which hashes
/// over BN254 Fr). Since 127 bits fits within both BN254 Fr (~254 bits) and
/// BN254 Fq (~254 bits), we reinterpret the standard-form limbs directly.
#[inline]
fn challenge_to_grumpkin_fr(challenge: Field<Bn254FrParams>) -> Fr {
    let standard = challenge.from_montgomery_form();
    Field::<Bn254FqParams>::from_limbs(standard.data)
}

/// IPA commitment scheme for a given polynomial length.
///
/// `LOG_N` is the log2 of the polynomial length (poly_length = 1 << LOG_N).
pub struct IPA<const LOG_N: usize>;

impl<const LOG_N: usize> IPA<LOG_N> {
    const POLY_LENGTH: usize = 1 << LOG_N;

    /// Compute an IPA opening proof for a single polynomial at a single evaluation point.
    ///
    /// Port of C++ `IPA::compute_opening_proof`.
    pub fn compute_opening_proof(
        ck: &CK,
        opening_claim: &ProverOpeningClaim<GrumpkinG1Params>,
        transcript: &mut NativeTranscript,
    ) {
        Self::add_prover_claim_to_hash_buffer(ck, opening_claim, transcript);
        Self::compute_opening_proof_internal(ck, opening_claim, transcript);
    }

    /// Natively verify the correctness of an IPA proof.
    ///
    /// Port of C++ `IPA::reduce_verify`.
    pub fn reduce_verify(
        vk: &VK,
        opening_claim: &OpeningClaim<GrumpkinG1Params>,
        transcript: &mut NativeTranscript,
    ) -> bool {
        Self::add_verifier_claim_to_hash_buffer(opening_claim, transcript);
        Self::reduce_verify_internal_native(vk, opening_claim, transcript)
    }

    /// Reduce a batch opening claim to a single opening claim.
    ///
    /// Computes C = sum(commitments_i * scalars_i) via MSM, returning an
    /// opening claim at (evaluation_point, 0).
    ///
    /// Port of C++ `IPA::reduce_batch_opening_claim`.
    pub fn reduce_batch_opening_claim(
        batch_opening_claim: &BatchOpeningClaim<GrumpkinG1Params>,
    ) -> OpeningClaim<GrumpkinG1Params> {
        let commitment = batch_mul_native::<GrumpkinG1Params>(
            &batch_opening_claim.commitments,
            &batch_opening_claim.scalars,
        );
        OpeningClaim {
            opening_pair: OpeningPair {
                challenge: batch_opening_claim.evaluation_point,
                evaluation: Fr::zero(),
            },
            commitment,
        }
    }

    /// Natively verify the IPA opening claim obtained from a batch accumulator.
    ///
    /// Port of C++ `IPA::reduce_verify_batch_opening_claim`.
    pub fn reduce_verify_batch_opening_claim(
        batch_opening_claim: &BatchOpeningClaim<GrumpkinG1Params>,
        vk: &VK,
        transcript: &mut NativeTranscript,
    ) -> bool {
        let opening_claim = Self::reduce_batch_opening_claim(batch_opening_claim);
        Self::add_verifier_claim_to_hash_buffer(&opening_claim, transcript);
        Self::reduce_verify_internal_native(vk, &opening_claim, transcript)
    }

    /// Evaluate the challenge polynomial at r.
    ///
    /// challenge_poly(X) = prod_{i in [k]} (1 + u_{len-i}^{-1} * X^{2^{i-1}})
    ///
    /// Port of C++ `IPA::evaluate_challenge_poly`.
    pub fn evaluate_challenge_poly(u_challenges_inv: &[Fr], r: Fr) -> Fr {
        let mut challenge_poly_eval = Fr::one();
        let mut r_pow = r;

        for i in 0..LOG_N - 1 {
            let monomial = u_challenges_inv[LOG_N - 1 - i] * r_pow;
            challenge_poly_eval = challenge_poly_eval * (Fr::one() + monomial);
            r_pow = r_pow.sqr();
        }
        // Last iteration without squaring r_pow
        let monomial = u_challenges_inv[0] * r_pow;
        challenge_poly_eval = challenge_poly_eval * (Fr::one() + monomial);
        challenge_poly_eval
    }

    /// Construct the s-vector from u_challenges_inv.
    ///
    /// s[r] = prod_i (u_i)^{-1 * r_i} where (r_i) is the binary expansion of r.
    ///
    /// Port of C++ `IPA::construct_poly_from_u_challenges_inv`.
    pub fn construct_s_vec_from_u_challenges_inv(u_challenges_inv: &[Fr]) -> Vec<Fr> {
        let poly_length = Self::POLY_LENGTH;

        let mut s_vec = vec![Fr::zero(); poly_length];
        let mut s_vec_temps = vec![Fr::zero(); poly_length / 2];

        // Double-buffering: arrange pointers so s_vec holds the final result.
        // If LOG_N is even, we swap initial assignments.
        let even_rounds = (LOG_N & 1) == 0;

        // Use raw pointers for the double-buffer swap pattern.
        // SAFETY: s_vec and s_vec_temps are disjoint heap allocations.
        let s_vec_ptr = s_vec.as_mut_ptr();
        let temps_ptr = s_vec_temps.as_mut_ptr();

        let (mut prev_ptr, mut curr_ptr) = if even_rounds {
            (s_vec_ptr, temps_ptr)
        } else {
            (temps_ptr, s_vec_ptr)
        };

        unsafe { *prev_ptr = Fr::one() };

        for i in 0..LOG_N {
            let round_size = 1 << (i + 1);
            let round_challenge = u_challenges_inv[i];
            for j in 0..round_size / 2 {
                unsafe {
                    *curr_ptr.add(j * 2) = *prev_ptr.add(j);
                    *curr_ptr.add(j * 2 + 1) = *prev_ptr.add(j) * round_challenge;
                }
            }
            std::mem::swap(&mut curr_ptr, &mut prev_ptr);
        }

        s_vec
    }

    // --- Internal methods ---

    /// Add the opening claim to the hash buffer (prover side).
    fn add_prover_claim_to_hash_buffer(
        ck: &CK,
        opening_claim: &ProverOpeningClaim<GrumpkinG1Params>,
        transcript: &mut NativeTranscript,
    ) {
        let commitment = ck.commit(&opening_claim.polynomial);
        transcript.add_to_hash_buffer("IPA:commitment", &commitment);
        transcript.add_to_hash_buffer("IPA:challenge", &opening_claim.opening_pair.challenge);
        transcript.add_to_hash_buffer("IPA:evaluation", &opening_claim.opening_pair.evaluation);
    }

    /// Add the opening claim to the hash buffer (verifier side).
    fn add_verifier_claim_to_hash_buffer(
        opening_claim: &OpeningClaim<GrumpkinG1Params>,
        transcript: &mut NativeTranscript,
    ) {
        transcript.add_to_hash_buffer("IPA:commitment", &opening_claim.commitment);
        transcript.add_to_hash_buffer("IPA:challenge", &opening_claim.opening_pair.challenge);
        transcript.add_to_hash_buffer("IPA:evaluation", &opening_claim.opening_pair.evaluation);
    }

    /// Get a challenge from the transcript and convert to GrumpkinFr.
    #[inline]
    fn get_challenge(transcript: &mut NativeTranscript, label: &str) -> Fr {
        challenge_to_grumpkin_fr(transcript.get_challenge(label))
    }

    /// Internal IPA prover.
    ///
    /// Port of C++ `IPA::compute_opening_proof_internal`.
    fn compute_opening_proof_internal(
        ck: &CK,
        opening_claim: &ProverOpeningClaim<GrumpkinG1Params>,
        transcript: &mut NativeTranscript,
    ) {
        let poly_length = Self::POLY_LENGTH;
        let polynomial = &opening_claim.polynomial;

        // Step 2: Receive challenge for the auxiliary generator
        let generator_challenge = Self::get_challenge(transcript, "IPA:generator_challenge");
        assert!(
            !generator_challenge.is_zero(),
            "The generator challenge can't be zero"
        );

        // Step 3: Compute auxiliary generator U = generator_challenge * G
        let aux_generator = GroupElement::one().mul(&generator_challenge).to_affine();

        assert!(
            poly_length > 0 && poly_length.is_power_of_two(),
            "The polynomial length should be positive and a power of two"
        );

        // Step 4: Set initial vector a to polynomial coefficients, load vector G
        let full_poly = polynomial.full();
        let mut a_vec: Vec<Fr> = full_poly.data().to_vec();
        a_vec.resize(poly_length, Fr::zero());

        let srs_elements = ck.srs_points();
        assert!(
            poly_length <= srs_elements.len(),
            "Not enough SRS points for IPA!"
        );

        let mut g_vec: Vec<Commitment> = srs_elements[..poly_length].to_vec();

        // Step 5: Compute vector b (powers of the challenge)
        let challenge = opening_claim.opening_pair.challenge;
        let mut b_vec = vec![Fr::zero(); poly_length];
        let mut b_power = Fr::one();
        for b in b_vec.iter_mut() {
            *b = b_power;
            b_power = b_power * challenge;
        }

        // Step 6: Perform IPA reduction rounds
        let mut round_size = poly_length;

        for i in 0..LOG_N {
            round_size /= 2;

            // Compute inner products
            let mut inner_prod_l = Fr::zero();
            let mut inner_prod_r = Fr::zero();
            for j in 0..round_size {
                inner_prod_l = inner_prod_l + a_vec[j] * b_vec[round_size + j];
                inner_prod_r = inner_prod_r + a_vec[round_size + j] * b_vec[j];
            }

            // Step 6.a: L_i = <a_vec_lo, G_vec_hi> + inner_prod_L * aux_generator
            let l_i = pippenger_msm::<GrumpkinG1Params>(
                &a_vec[..round_size],
                &g_vec[round_size..round_size * 2],
            ) + Element::from_affine(&aux_generator).mul(&inner_prod_l);

            // Step 6.b: R_i = <a_vec_hi, G_vec_lo> + inner_prod_R * aux_generator
            let r_i = pippenger_msm::<GrumpkinG1Params>(
                &a_vec[round_size..round_size * 2],
                &g_vec[..round_size],
            ) + Element::from_affine(&aux_generator).mul(&inner_prod_r);

            // Step 6.c: Send L_i and R_i to verifier
            let index = (LOG_N - i - 1).to_string();
            transcript.send_to_verifier(&format!("IPA:L_{index}"), &l_i.to_affine());
            transcript.send_to_verifier(&format!("IPA:R_{index}"), &r_i.to_affine());

            // Step 6.d: Receive round challenge
            let round_challenge =
                Self::get_challenge(transcript, &format!("IPA:round_challenge_{index}"));
            assert!(!round_challenge.is_zero(), "IPA round challenge is zero");
            let round_challenge_inv = round_challenge.invert();

            // Step 6.e: G_vec_new = G_vec_lo + G_vec_hi * round_challenge_inv
            for j in 0..round_size {
                let g_hi = Element::<GrumpkinG1Params>::from_affine(&g_vec[round_size + j]);
                let g_lo = Element::<GrumpkinG1Params>::from_affine(&g_vec[j]);
                g_vec[j] = (g_lo + g_hi.mul(&round_challenge_inv)).to_affine();
            }

            // Steps 6.f and 6.g: Update a_vec and b_vec
            for j in 0..round_size {
                a_vec[j] = a_vec[j] + round_challenge * a_vec[round_size + j];
                b_vec[j] = b_vec[j] + round_challenge_inv * b_vec[round_size + j];
            }
        }

        // Step 7: Send G_0 to verifier
        transcript.send_to_verifier("IPA:G_0", &g_vec[0]);

        // Step 8: Send a_0 to verifier
        transcript.send_to_verifier("IPA:a_0", &a_vec[0]);
    }

    /// Internal native IPA verifier.
    ///
    /// Port of C++ `IPA::reduce_verify_internal_native`.
    fn reduce_verify_internal_native(
        vk: &VK,
        opening_claim: &OpeningClaim<GrumpkinG1Params>,
        transcript: &mut NativeTranscript,
    ) -> bool {
        let poly_length = Self::POLY_LENGTH;

        // Step 2: Receive generator challenge and compute auxiliary generator
        let generator_challenge = Self::get_challenge(transcript, "IPA:generator_challenge");
        assert!(
            !generator_challenge.is_zero(),
            "The generator challenge can't be zero"
        );
        let aux_generator = GroupElement::one().mul(&generator_challenge).to_affine();

        // Step 3: Compute C' = C + f(beta) * U
        let c_prime = Element::from_affine(&opening_claim.commitment)
            + Element::from_affine(&aux_generator).mul(&opening_claim.opening_pair.evaluation);

        let pippenger_size = 2 * LOG_N;
        let mut round_challenges = vec![Fr::zero(); LOG_N];
        let mut msm_elements = vec![Commitment::infinity(); pippenger_size];
        let mut msm_scalars = vec![Fr::zero(); pippenger_size];

        // Step 4: Receive all L_i and R_i and populate MSM elements
        for i in 0..LOG_N {
            let index = (LOG_N - i - 1).to_string();
            let element_l: Commitment =
                transcript.receive_from_prover(&format!("IPA:L_{index}"));
            let element_r: Commitment =
                transcript.receive_from_prover(&format!("IPA:R_{index}"));
            round_challenges[i] =
                Self::get_challenge(transcript, &format!("IPA:round_challenge_{index}"));
            assert!(
                !round_challenges[i].is_zero(),
                "Round challenges can't be zero"
            );
            msm_elements[2 * i] = element_l;
            msm_elements[2 * i + 1] = element_r;
        }

        // Batch invert round challenges
        let mut round_challenges_inv = round_challenges.clone();
        Field::batch_invert(&mut round_challenges_inv);

        // Populate MSM scalars
        for i in 0..LOG_N {
            msm_scalars[2 * i] = round_challenges_inv[i];
            msm_scalars[2 * i + 1] = round_challenges[i];
        }

        // Step 5: Compute C_zero = C' + sum(u_j^{-1} * L_j + u_j * R_j)
        let lr_sums = pippenger_msm::<GrumpkinG1Params>(&msm_scalars, &msm_elements);
        let c_zero = c_prime + lr_sums;

        // Step 6: Compute b_zero succinctly
        let b_zero = Self::evaluate_challenge_poly(
            &round_challenges_inv,
            opening_claim.opening_pair.challenge,
        );

        // Step 7: Construct s-vector
        let s_vec =
            Self::construct_s_vec_from_u_challenges_inv(&round_challenges_inv[..LOG_N]);

        let srs_elements = vk.get_monomial_points();
        assert!(
            poly_length <= srs_elements.len(),
            "Not enough SRS points for IPA!"
        );

        // Step 8: Compute G_zero = <s_vec, G_vec>
        let g_zero =
            pippenger_msm::<GrumpkinG1Params>(&s_vec, &srs_elements[..poly_length]).to_affine();
        let g_zero_sent: Commitment = transcript.receive_from_prover("IPA:G_0");
        assert_eq!(
            g_zero, g_zero_sent,
            "G_0 should equal G_0 sent in transcript. IPA verification fails."
        );

        // Step 9: Receive a_zero from the prover
        let a_zero: Fr = transcript.receive_from_prover("IPA:a_0");

        // Step 10: Compute C_right = a_0 * G_s + a_0 * b_0 * U
        let right_hand_side = Element::from_affine(&g_zero).mul(&a_zero)
            + Element::from_affine(&aux_generator).mul(&(a_zero * b_zero));

        // Step 11: Check C_right == C_zero
        c_zero.to_affine() == right_hand_side.to_affine()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_polynomials::polynomial::Polynomial;

    const LOG_N: usize = 7;
    const N: usize = 1 << LOG_N;

    type TestIPA = IPA<LOG_N>;

    fn setup_crs() {
        let generator = GroupElement::from_affine(&Commitment::one());
        let mut acc = generator;
        let mut points = Vec::with_capacity(N);
        points.push(acc.to_affine());
        for _ in 1..N {
            acc = acc + generator;
            points.push(acc.to_affine());
        }
        bbrs_srs::global_crs::init_grumpkin_mem_crs_factory(&points);
    }

    fn create_ck() -> CK {
        CommitmentKey::<GrumpkinG1Params>::new(N)
    }

    fn create_vk() -> VK {
        GrumpkinVerifierCommitmentKey::new(N)
    }

    fn run_native_prove_verify(poly: &Polynomial<Bn254FqParams>, x: Fr) -> bool {
        let ck = create_ck();
        let vk = create_vk();
        let commitment = ck.commit(poly);
        let eval = poly.evaluate(&x);
        let opening_claim = OpeningClaim {
            opening_pair: OpeningPair {
                challenge: x,
                evaluation: eval,
            },
            commitment,
        };

        // Prover
        let mut prover_transcript = NativeTranscript::new();
        let prover_claim = ProverOpeningClaim {
            polynomial: poly.clone(),
            opening_pair: OpeningPair {
                challenge: x,
                evaluation: eval,
            },
            gemini_fold: false,
        };
        TestIPA::compute_opening_proof(&ck, &prover_claim, &mut prover_transcript);

        // Verifier
        let proof = prover_transcript.export_proof();
        let mut verifier_transcript = NativeTranscript::from_proof(&proof);
        TestIPA::reduce_verify(&vk, &opening_claim, &mut verifier_transcript)
    }

    #[test]
    fn commit_on_many_zero_coeff_poly_works() {
        setup_crs();
        let ck = create_ck();

        let n = 16;
        let mut coeffs = vec![Fr::zero(); n];
        coeffs[3] = Fr::random_element();
        let poly = Polynomial::from_coefficients(coeffs.clone(), n);

        let commitment = ck.commit(&poly);
        let srs = ck.srs_points();

        let mut expected = Element::<GrumpkinG1Params>::from_affine(&srs[0]).mul(&coeffs[0]);
        for i in 1..n {
            expected =
                expected + Element::<GrumpkinG1Params>::from_affine(&srs[i]).mul(&coeffs[i]);
        }
        assert_eq!(expected.to_affine(), commitment);
    }

    #[test]
    fn commit_to_zero_poly() {
        setup_crs();
        let ck = create_ck();

        let poly = Polynomial::new(N, N, 0);
        let commitment = ck.commit(&poly);
        assert!(commitment.is_point_at_infinity());

        let x = Fr::random_element();
        let eval = poly.evaluate(&x);
        assert_eq!(eval, Fr::zero());
    }

    #[test]
    fn commit_to_random_poly() {
        setup_crs();
        let ck = create_ck();

        let poly = Polynomial::random(N, N, 0);
        let commitment = ck.commit(&poly);
        let srs = ck.srs_points();
        let coeffs = poly.data();

        let mut expected = Element::<GrumpkinG1Params>::from_affine(&srs[0]).mul(&coeffs[0]);
        for i in 1..N {
            expected =
                expected + Element::<GrumpkinG1Params>::from_affine(&srs[i]).mul(&coeffs[i]);
        }
        assert_eq!(expected.to_affine(), commitment);
    }

    #[test]
    fn open_zero_polynomial() {
        setup_crs();
        let poly = Polynomial::new(N, N, 0);
        let x = Fr::random_element();
        assert!(run_native_prove_verify(&poly, x));
    }

    #[test]
    fn open_many_zeros_polynomial() {
        setup_crs();

        let mut coeffs_even = vec![Fr::zero(); N];
        let mut coeffs_odd = vec![Fr::zero(); N];
        for i in 0..N / 2 {
            coeffs_even[2 * i] = Fr::random_element();
            coeffs_odd[2 * i + 1] = Fr::random_element();
        }
        let poly_even = Polynomial::from_coefficients(coeffs_even, N);
        let poly_odd = Polynomial::from_coefficients(coeffs_odd, N);

        let x = Fr::random_element();
        assert!(run_native_prove_verify(&poly_even, x));
        assert!(run_native_prove_verify(&poly_odd, x));
    }

    #[test]
    fn open_at_zero() {
        setup_crs();
        let poly = Polynomial::random(N, N, 0);
        let x = Fr::zero();
        assert!(run_native_prove_verify(&poly, x));
    }

    #[test]
    fn open_random() {
        setup_crs();
        let poly = Polynomial::random(N, N, 0);
        let x = Fr::random_element();
        assert!(run_native_prove_verify(&poly, x));
    }

    #[test]
    fn opening_value_zero() {
        setup_crs();
        let mut poly = Polynomial::random(N, N, 0);
        let x = Fr::random_element();
        let initial_evaluation = poly.evaluate(&x);
        let change_in_linear_coefficient = initial_evaluation * x.invert();
        let data = poly.data_mut();
        data[1] = data[1] - change_in_linear_coefficient;

        assert_eq!(poly.evaluate(&x), Fr::zero());
        assert!(run_native_prove_verify(&poly, x));
    }
}
