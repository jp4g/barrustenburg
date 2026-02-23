//! Ultra verification key — commitments to precomputed polynomials.
//!
//! Port of `ultra_flavor.hpp` VerificationKey concept.

use bbrs_commitment_schemes::verification_key::Bn254VerifierCommitmentKey;
use bbrs_ecc::curves::bn254::{G1Affine, G2AffineElement};

use crate::proving_key::ProvingKey;

/// Ultra verification key containing commitments to all precomputed polynomials.
pub struct VerificationKey {
    pub circuit_size: usize,
    pub log_circuit_size: usize,
    pub num_public_inputs: usize,
    pub pub_inputs_offset: usize,
    // Selector commitments
    pub q_m: G1Affine,
    pub q_c: G1Affine,
    pub q_l: G1Affine,
    pub q_r: G1Affine,
    pub q_o: G1Affine,
    pub q_4: G1Affine,
    pub q_lookup: G1Affine,
    pub q_arith: G1Affine,
    pub q_delta_range: G1Affine,
    pub q_elliptic: G1Affine,
    pub q_memory: G1Affine,
    pub q_nnf: G1Affine,
    pub q_poseidon2_external: G1Affine,
    pub q_poseidon2_internal: G1Affine,
    // Permutation commitments
    pub sigma_1: G1Affine,
    pub sigma_2: G1Affine,
    pub sigma_3: G1Affine,
    pub sigma_4: G1Affine,
    // Identity commitments
    pub id_1: G1Affine,
    pub id_2: G1Affine,
    pub id_3: G1Affine,
    pub id_4: G1Affine,
    // Table commitments
    pub table_1: G1Affine,
    pub table_2: G1Affine,
    pub table_3: G1Affine,
    pub table_4: G1Affine,
    // Lagrange commitments
    pub lagrange_first: G1Affine,
    pub lagrange_last: G1Affine,
    // Verifier commitment key for pairing check
    pub pcs_verification_key: Bn254VerifierCommitmentKey,
}

impl VerificationKey {
    /// Compute verification key from a proving key by committing to all precomputed polynomials.
    ///
    /// Uses `G2AffineElement::infinity()` as the G2 SRS point. For a correct pairing check,
    /// use `create_with_g2x` instead with the actual G2 SRS point.
    pub fn create(proving_key: &ProvingKey) -> Self {
        let ck = &proving_key.commitment_key;

        Self {
            circuit_size: proving_key.circuit_size,
            log_circuit_size: proving_key.log_circuit_size,
            num_public_inputs: proving_key.num_public_inputs,
            pub_inputs_offset: proving_key.pub_inputs_offset,
            q_m: ck.commit(&proving_key.polynomials.q_m),
            q_c: ck.commit(&proving_key.polynomials.q_c),
            q_l: ck.commit(&proving_key.polynomials.q_l),
            q_r: ck.commit(&proving_key.polynomials.q_r),
            q_o: ck.commit(&proving_key.polynomials.q_o),
            q_4: ck.commit(&proving_key.polynomials.q_4),
            q_lookup: ck.commit(&proving_key.polynomials.q_lookup),
            q_arith: ck.commit(&proving_key.polynomials.q_arith),
            q_delta_range: ck.commit(&proving_key.polynomials.q_delta_range),
            q_elliptic: ck.commit(&proving_key.polynomials.q_elliptic),
            q_memory: ck.commit(&proving_key.polynomials.q_memory),
            q_nnf: ck.commit(&proving_key.polynomials.q_nnf),
            q_poseidon2_external: ck.commit(&proving_key.polynomials.q_poseidon2_external),
            q_poseidon2_internal: ck.commit(&proving_key.polynomials.q_poseidon2_internal),
            sigma_1: ck.commit(&proving_key.polynomials.sigma_1),
            sigma_2: ck.commit(&proving_key.polynomials.sigma_2),
            sigma_3: ck.commit(&proving_key.polynomials.sigma_3),
            sigma_4: ck.commit(&proving_key.polynomials.sigma_4),
            id_1: ck.commit(&proving_key.polynomials.id_1),
            id_2: ck.commit(&proving_key.polynomials.id_2),
            id_3: ck.commit(&proving_key.polynomials.id_3),
            id_4: ck.commit(&proving_key.polynomials.id_4),
            table_1: ck.commit(&proving_key.polynomials.table_1),
            table_2: ck.commit(&proving_key.polynomials.table_2),
            table_3: ck.commit(&proving_key.polynomials.table_3),
            table_4: ck.commit(&proving_key.polynomials.table_4),
            lagrange_first: ck.commit(&proving_key.polynomials.lagrange_first),
            lagrange_last: ck.commit(&proving_key.polynomials.lagrange_last),
            pcs_verification_key: Bn254VerifierCommitmentKey::with_g2x(
                G2AffineElement::infinity(),
            ),
        }
    }

    /// Compute verification key from a proving key with a specific G2 SRS point.
    ///
    /// Use this in tests where a custom SRS is generated.
    pub fn create_with_g2x(proving_key: &ProvingKey, g2_x: G2AffineElement) -> Self {
        let mut vk = Self::create(proving_key);
        vk.pcs_verification_key = Bn254VerifierCommitmentKey::with_g2x(g2_x);
        vk
    }

    /// Get all precomputed commitments as an ordered array (for claim batcher).
    pub fn get_all_commitments(&self) -> [G1Affine; 28] {
        [
            self.q_m, self.q_c, self.q_l, self.q_r, self.q_o, self.q_4,
            self.q_lookup, self.q_arith, self.q_delta_range, self.q_elliptic,
            self.q_memory, self.q_nnf, self.q_poseidon2_external, self.q_poseidon2_internal,
            self.sigma_1, self.sigma_2, self.sigma_3, self.sigma_4,
            self.id_1, self.id_2, self.id_3, self.id_4,
            self.table_1, self.table_2, self.table_3, self.table_4,
            self.lagrange_first, self.lagrange_last,
        ]
    }
}
