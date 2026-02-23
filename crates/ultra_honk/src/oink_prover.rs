//! Port of `oink_prover.hpp`/`.cpp` — Oink Prover for Ultra Honk.
//!
//! The Oink prover executes the pre-sumcheck commitment rounds:
//! 1. Preamble: send circuit size, public inputs
//! 2. Wire commitments: commit to w_l, w_r, w_o
//! 3. Sorted list accumulator: compute and commit to w_4, lookup read counts/tags
//! 4. Log-derivative inverse: compute and commit to lookup_inverses
//! 5. Grand product: compute and commit to z_perm
//! 6. Generate alpha challenge for sumcheck


use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};
use bbrs_flavor::ultra_flavor::ProverPolynomials;
use bbrs_honk::compute_public_input_delta;
use bbrs_relations::relation_parameters::RelationParameters;
use bbrs_transcript::NativeTranscript;

use crate::proving_key::ProvingKey;

/// Output data from the Oink prover rounds.
pub struct OinkOutput {
    /// Relation parameters with beta, gamma, eta challenges filled in.
    pub relation_parameters: RelationParameters<Fr>,
    /// The alpha challenge used to batch subrelations in sumcheck.
    pub alpha: Fr,
}

/// Oink prover for Ultra Honk.
///
/// Port of C++ `OinkProver<UltraFlavor>`.
pub struct OinkProver<'a> {
    proving_key: &'a mut ProvingKey,
}

impl<'a> OinkProver<'a> {
    pub fn new(proving_key: &'a mut ProvingKey) -> Self {
        Self { proving_key }
    }

    /// Execute the Oink rounds and return the output.
    pub fn prove(self, transcript: &mut NativeTranscript) -> OinkOutput {
        let pk = self.proving_key;

        // Round 1: Preamble — send circuit size, num public inputs, public inputs
        transcript.send_to_verifier("circuit_size", &(pk.circuit_size as u64));
        transcript.send_to_verifier("public_input_size", &(pk.num_public_inputs as u64));
        transcript.send_to_verifier("pub_inputs_offset", &(pk.pub_inputs_offset as u64));
        for (i, pi) in pk.public_inputs.iter().enumerate() {
            transcript.send_to_verifier(&format!("public_input_{}", i), pi);
        }

        // Round 2: Wire commitments — commit to w_l, w_r, w_o
        let comm_w_l = pk.commitment_key.commit(&pk.polynomials.w_l);
        let comm_w_r = pk.commitment_key.commit(&pk.polynomials.w_r);
        let comm_w_o = pk.commitment_key.commit(&pk.polynomials.w_o);
        transcript.send_to_verifier("W_L", &comm_w_l);
        transcript.send_to_verifier("W_R", &comm_w_r);
        transcript.send_to_verifier("W_O", &comm_w_o);

        // Get eta challenges for sorted list accumulator + lookup
        let eta = transcript.get_challenge("eta");
        let eta_two = transcript.get_challenge("eta_two");
        let eta_three = transcript.get_challenge("eta_three");

        // Round 3: Compute lookup read counts and tags
        compute_lookup_read_counts(pk, eta, eta_two, eta_three);

        // Commit to w_4, lookup_read_counts, lookup_read_tags
        let comm_w_4 = pk.commitment_key.commit(&pk.polynomials.w_4);
        let comm_lookup_read_counts = pk.commitment_key.commit(&pk.polynomials.lookup_read_counts);
        let comm_lookup_read_tags = pk.commitment_key.commit(&pk.polynomials.lookup_read_tags);
        transcript.send_to_verifier("W_4", &comm_w_4);
        transcript.send_to_verifier("LOOKUP_READ_COUNTS", &comm_lookup_read_counts);
        transcript.send_to_verifier("LOOKUP_READ_TAGS", &comm_lookup_read_tags);

        // Get beta, gamma challenges
        let beta = transcript.get_challenge("beta");
        let gamma = transcript.get_challenge("gamma");

        // Compute public input delta
        let pub_inputs_offset_field = Fr::from(pk.pub_inputs_offset as u64);
        let public_input_delta = compute_public_input_delta::<Bn254FrParams>(
            &pk.public_inputs,
            beta,
            gamma,
            pub_inputs_offset_field,
        );

        let beta_sqr = beta * beta;
        let beta_cube = beta_sqr * beta;

        let relation_parameters = RelationParameters {
            eta,
            eta_two,
            eta_three,
            beta,
            gamma,
            public_input_delta,
            beta_sqr,
            beta_cube,
            ..RelationParameters::default()
        };

        // Round 4: Compute lookup inverses
        compute_lookup_inverses(&mut pk.polynomials, &relation_parameters, pk.circuit_size);

        let comm_lookup_inverses = pk.commitment_key.commit(&pk.polynomials.lookup_inverses);
        transcript.send_to_verifier("LOOKUP_INVERSES", &comm_lookup_inverses);

        // Round 5: Compute grand product z_perm
        compute_grand_product_perm(&mut pk.polynomials, &relation_parameters, pk.circuit_size);

        // Set shifted polynomials
        pk.polynomials.set_shifted();

        let comm_z_perm = pk.commitment_key.commit(&pk.polynomials.z_perm);
        transcript.send_to_verifier("Z_PERM", &comm_z_perm);

        // Round 6: Get alpha challenge for sumcheck
        let alpha = transcript.get_challenge("Sumcheck:alpha");

        OinkOutput {
            relation_parameters,
            alpha,
        }
    }
}

/// Compute lookup read counts and tags from the proving key's lookup block.
fn compute_lookup_read_counts(
    pk: &mut ProvingKey,
    eta: Fr,
    eta_two: Fr,
    eta_three: Fr,
) {
    let circuit_size = pk.circuit_size;

    // Phase 1: Collect table entries (combined value, row index)
    let table_entries: Vec<(Fr, usize)> = {
        let polys = &pk.polynomials;
        let mut entries = Vec::new();
        for i in 0..circuit_size {
            let t1 = polys.table_1.get(i);
            let t2 = polys.table_2.get(i);
            let t3 = polys.table_3.get(i);
            let t4 = polys.table_4.get(i);

            if t1.is_zero() && t2.is_zero() && t3.is_zero() && t4.is_zero() {
                continue;
            }

            let combined = t1 + t2 * eta + t3 * eta_two + t4 * eta_three;
            entries.push((combined, i));
        }
        entries
    };

    // Phase 2: Collect lookup gate combined values (combined value, gate index)
    let lookup_entries: Vec<Fr> = {
        let polys = &pk.polynomials;
        let mut entries = Vec::with_capacity(circuit_size);
        for i in 0..circuit_size {
            let q_lookup = polys.q_lookup.get(i);
            if q_lookup.is_zero() {
                entries.push(Fr::zero());
                continue;
            }

            let w_1 = polys.w_l.get(i);
            let w_2 = polys.w_r.get(i);
            let w_3 = polys.w_o.get(i);
            let q_r = polys.q_r.get(i);
            let q_c = polys.q_c.get(i);

            let combined = w_1 + q_r * w_2 * eta + w_3 * eta_two + q_c * eta_three;
            entries.push(combined);
        }
        entries
    };

    // Phase 3: Match lookups to table entries (no immutable borrow alive)
    for i in 0..circuit_size {
        if pk.polynomials.q_lookup.get(i).is_zero() {
            continue;
        }

        let combined = lookup_entries[i];

        for &(table_combined, trow) in &table_entries {
            if combined == table_combined {
                *pk.polynomials.lookup_read_counts.at_mut(trow) =
                    pk.polynomials.lookup_read_counts.get(trow) + Fr::one();
                *pk.polynomials.lookup_read_tags.at_mut(trow) = Fr::one();
                break;
            }
        }
    }
}

/// Compute the grand product z_perm polynomial.
fn compute_grand_product_perm(
    polys: &mut ProverPolynomials<Bn254FrParams>,
    params: &RelationParameters<Fr>,
    circuit_size: usize,
) {
    if circuit_size <= 1 {
        return;
    }

    let beta = params.beta;
    let gamma = params.gamma;
    let iteration_size = circuit_size - 1;

    let mut numerators = vec![Fr::one(); iteration_size];
    let mut denominators = vec![Fr::one(); iteration_size];

    for i in 0..iteration_size {
        let num = (polys.w_l.get(i) + polys.id_1.get(i) * beta + gamma)
            * (polys.w_r.get(i) + polys.id_2.get(i) * beta + gamma)
            * (polys.w_o.get(i) + polys.id_3.get(i) * beta + gamma)
            * (polys.w_4.get(i) + polys.id_4.get(i) * beta + gamma);

        let den = (polys.w_l.get(i) + polys.sigma_1.get(i) * beta + gamma)
            * (polys.w_r.get(i) + polys.sigma_2.get(i) * beta + gamma)
            * (polys.w_o.get(i) + polys.sigma_3.get(i) * beta + gamma)
            * (polys.w_4.get(i) + polys.sigma_4.get(i) * beta + gamma);

        numerators[i] = num;
        denominators[i] = den;
    }

    // Running products
    for i in 0..iteration_size - 1 {
        let prev_num = numerators[i];
        numerators[i + 1] = numerators[i + 1] * prev_num;
        let prev_den = denominators[i];
        denominators[i + 1] = denominators[i + 1] * prev_den;
    }

    // Batch invert denominators
    batch_invert_in_place(&mut denominators);

    // Z[0] = 0, Z[i+1] = running_numerator[i] / running_denominator[i]
    for i in 0..iteration_size {
        *polys.z_perm.at_mut(i + 1) = numerators[i] * denominators[i];
    }
}

/// Compute lookup inverse polynomial.
fn compute_lookup_inverses(
    polys: &mut ProverPolynomials<Bn254FrParams>,
    params: &RelationParameters<Fr>,
    circuit_size: usize,
) {
    let mut denominators = vec![Fr::zero(); circuit_size];
    let mut has_inverse = vec![false; circuit_size];

    for i in 0..circuit_size {
        let q_lookup = polys.q_lookup.get(i);
        let lookup_read_counts = polys.lookup_read_counts.get(i);

        if q_lookup.is_zero() && lookup_read_counts.is_zero() {
            continue;
        }

        has_inverse[i] = true;

        let w_1 = polys.w_l.get(i);
        let w_2 = polys.w_r.get(i);
        let w_3 = polys.w_o.get(i);
        let q_r = polys.q_r.get(i);
        let q_c = polys.q_c.get(i);
        let read_denom = w_1 + params.gamma + q_r * w_2 * params.eta
            + w_3 * params.eta_two
            + q_c * params.eta_three;

        let t1 = polys.table_1.get(i);
        let t2 = polys.table_2.get(i);
        let t3 = polys.table_3.get(i);
        let t4 = polys.table_4.get(i);
        let write_denom =
            t1 + params.gamma + t2 * params.eta + t3 * params.eta_two + t4 * params.eta_three;

        denominators[i] = read_denom * write_denom;
    }

    batch_invert_nonzero(&mut denominators);

    for i in 0..circuit_size {
        if has_inverse[i] {
            *polys.lookup_inverses.at_mut(i) = denominators[i];
        }
    }
}

fn batch_invert_in_place(values: &mut [Fr]) {
    if values.is_empty() {
        return;
    }
    let n = values.len();
    let mut scratch = vec![Fr::one(); n];
    scratch[0] = values[0];
    for i in 1..n {
        scratch[i] = scratch[i - 1] * values[i];
    }
    let mut inv_acc = scratch[n - 1].invert();
    for i in (1..n).rev() {
        let temp = values[i];
        values[i] = scratch[i - 1] * inv_acc;
        inv_acc = inv_acc * temp;
    }
    values[0] = inv_acc;
}

fn batch_invert_nonzero(values: &mut [Fr]) {
    if values.is_empty() {
        return;
    }
    let nonzero_indices: Vec<usize> = values
        .iter()
        .enumerate()
        .filter(|(_, v)| !v.is_zero())
        .map(|(i, _)| i)
        .collect();
    if nonzero_indices.is_empty() {
        return;
    }
    let n = nonzero_indices.len();
    let mut scratch = vec![Fr::one(); n];
    scratch[0] = values[nonzero_indices[0]];
    for i in 1..n {
        scratch[i] = scratch[i - 1] * values[nonzero_indices[i]];
    }
    let mut inv_acc = scratch[n - 1].invert();
    for i in (1..n).rev() {
        let idx = nonzero_indices[i];
        let temp = values[idx];
        values[idx] = scratch[i - 1] * inv_acc;
        inv_acc = inv_acc * temp;
    }
    values[nonzero_indices[0]] = inv_acc;
}
