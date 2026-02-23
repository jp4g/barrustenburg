//! Proving key construction from UltraCircuitBuilder.
//!
//! Port of `decider_proving_key.hpp/.cpp` — constructs prover polynomials from
//! the circuit builder's execution trace.

use std::collections::BTreeMap;

use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
use bbrs_commitment_schemes::commitment_key::CommitmentKey;
use bbrs_ecc::curves::bn254::{Bn254FrParams, Bn254G1Params, Fr};
use bbrs_flavor::ultra_flavor::ProverPolynomials;
use bbrs_honk::PERMUTATION_ARGUMENT_VALUE_SEPARATOR;

/// The Ultra proving key bundles everything the prover needs.
pub struct ProvingKey {
    /// Circuit size (dyadic, power of 2).
    pub circuit_size: usize,
    /// Log2 of circuit size.
    pub log_circuit_size: usize,
    /// Number of public inputs.
    pub num_public_inputs: usize,
    /// Public input values.
    pub public_inputs: Vec<Fr>,
    /// Offset of the public inputs block in the trace.
    pub pub_inputs_offset: usize,
    /// The prover polynomials (all 41 columns).
    pub polynomials: ProverPolynomials<Bn254FrParams>,
    /// Commitment key for polynomial commitments.
    pub commitment_key: CommitmentKey<Bn254G1Params>,
}

impl ProvingKey {
    /// Construct a proving key from a finalized circuit builder.
    ///
    /// This extracts the execution trace into ProverPolynomials and creates
    /// the commitment key from the SRS.
    pub fn create(builder: &mut UltraCircuitBuilder<Bn254FrParams>) -> Self {
        // Finalize the circuit if not already done
        if !builder.circuit_finalized {
            builder.finalize_circuit(true);
        }

        // Compute trace offsets
        builder.blocks.compute_offsets();

        // Compute circuit size (round up to next power of 2)
        let total_content = builder.blocks.get_total_content_size();
        let circuit_size = (total_content + 1).next_power_of_two(); // +1 for zero row
        let log_circuit_size = circuit_size.trailing_zeros() as usize;

        // Collect public inputs
        let num_public_inputs = builder.base.num_public_inputs();
        let public_inputs: Vec<Fr> = builder
            .base
            .public_inputs()
            .iter()
            .map(|&idx| builder.base.get_variable(idx))
            .collect();
        let pub_inputs_offset = builder.blocks.pub_inputs.block.trace_offset as usize;

        // Create prover polynomials
        let mut polynomials = ProverPolynomials::new(circuit_size);

        // Step 1: Populate wire values and selectors from execution trace
        populate_wires_and_selectors(builder, &mut polynomials);

        // Step 2: Populate identity and sigma polynomials
        populate_id_and_sigma(builder, &mut polynomials, circuit_size);

        // Step 3: Set lagrange polynomials
        *polynomials.lagrange_first.at_mut(0) = Fr::one();
        *polynomials.lagrange_last.at_mut(circuit_size - 1) = Fr::one();

        // Step 4: Populate table polynomials from lookup tables
        populate_tables(builder, &mut polynomials);

        // Step 5: Create commitment key
        let commitment_key = CommitmentKey::<Bn254G1Params>::new(circuit_size);

        Self {
            circuit_size,
            log_circuit_size,
            num_public_inputs,
            public_inputs,
            pub_inputs_offset,
            polynomials,
            commitment_key,
        }
    }

    /// Construct a proving key from a finalized circuit builder with a custom
    /// commitment key (useful for testing without global CRS).
    pub fn create_with_ck(
        builder: &mut UltraCircuitBuilder<Bn254FrParams>,
        commitment_key: CommitmentKey<Bn254G1Params>,
    ) -> Self {
        if !builder.circuit_finalized {
            builder.finalize_circuit(true);
        }

        builder.blocks.compute_offsets();

        let total_content = builder.blocks.get_total_content_size();
        let circuit_size = (total_content + 1).next_power_of_two();
        let log_circuit_size = circuit_size.trailing_zeros() as usize;

        let num_public_inputs = builder.base.num_public_inputs();
        let public_inputs: Vec<Fr> = builder
            .base
            .public_inputs()
            .iter()
            .map(|&idx| builder.base.get_variable(idx))
            .collect();
        let pub_inputs_offset = builder.blocks.pub_inputs.block.trace_offset as usize;

        let mut polynomials = ProverPolynomials::new(circuit_size);

        populate_wires_and_selectors(builder, &mut polynomials);
        populate_id_and_sigma(builder, &mut polynomials, circuit_size);

        *polynomials.lagrange_first.at_mut(0) = Fr::one();
        *polynomials.lagrange_last.at_mut(circuit_size - 1) = Fr::one();

        populate_tables(builder, &mut polynomials);

        Self {
            circuit_size,
            log_circuit_size,
            num_public_inputs,
            public_inputs,
            pub_inputs_offset,
            polynomials,
            commitment_key,
        }
    }
}

/// Populate wire values and selector values from the execution trace blocks.
fn populate_wires_and_selectors(
    builder: &UltraCircuitBuilder<Bn254FrParams>,
    polys: &mut ProverPolynomials<Bn254FrParams>,
) {
    for block in builder.blocks.get() {
        let offset = block.block.trace_offset as usize;
        let block_size = block.size();

        for row in 0..block_size {
            let trace_row = offset + row;

            // Wires: dereference variable indices to get field values
            *polys.w_l.at_mut(trace_row) =
                builder.base.get_variable(block.block.wires[0][row]);
            *polys.w_r.at_mut(trace_row) =
                builder.base.get_variable(block.block.wires[1][row]);
            *polys.w_o.at_mut(trace_row) =
                builder.base.get_variable(block.block.wires[2][row]);
            *polys.w_4.at_mut(trace_row) =
                builder.base.get_variable(block.block.wires[3][row]);

            // Non-gate selectors: q_m, q_c, q_l, q_r, q_o, q_4
            *polys.q_m.at_mut(trace_row) = block.block.q_m().get(row);
            *polys.q_c.at_mut(trace_row) = block.block.q_c().get(row);
            *polys.q_l.at_mut(trace_row) = block.block.q_1().get(row);
            *polys.q_r.at_mut(trace_row) = block.block.q_2().get(row);
            *polys.q_o.at_mut(trace_row) = block.block.q_3().get(row);
            *polys.q_4.at_mut(trace_row) = block.block.q_4().get(row);

            // Gate selectors
            *polys.q_arith.at_mut(trace_row) = block.q_arith().get(row);
            *polys.q_lookup.at_mut(trace_row) = block.q_lookup().get(row);
            *polys.q_delta_range.at_mut(trace_row) = block.q_delta_range().get(row);
            *polys.q_elliptic.at_mut(trace_row) = block.q_elliptic().get(row);
            *polys.q_memory.at_mut(trace_row) = block.q_memory().get(row);
            *polys.q_nnf.at_mut(trace_row) = block.q_nnf().get(row);
            *polys.q_poseidon2_external.at_mut(trace_row) =
                block.q_poseidon2_external().get(row);
            *polys.q_poseidon2_internal.at_mut(trace_row) =
                block.q_poseidon2_internal().get(row);
        }
    }
}

/// Populate identity and sigma polynomials from the circuit builder's copy constraints.
fn populate_id_and_sigma(
    builder: &UltraCircuitBuilder<Bn254FrParams>,
    polys: &mut ProverPolynomials<Bn254FrParams>,
    circuit_size: usize,
) {
    let separator = Fr::from(PERMUTATION_ARGUMENT_VALUE_SEPARATOR);

    // Collect public input row positions
    let pub_inputs_offset = builder.blocks.pub_inputs.block.trace_offset as usize;
    let pub_inputs_size = builder.blocks.pub_inputs.size();
    let pub_input_rows: Vec<usize> = (0..pub_inputs_size)
        .map(|i| pub_inputs_offset + i)
        .collect();

    // Initialize identity polynomials: id_k[i] = i + k * separator
    let id_polys = [&mut polys.id_1, &mut polys.id_2, &mut polys.id_3, &mut polys.id_4];
    for (wire_idx, id_poly) in id_polys.into_iter().enumerate() {
        let wire_offset = separator * Fr::from(wire_idx as u64);
        for i in 0..circuit_size {
            *id_poly.at_mut(i) = Fr::from(i as u64) + wire_offset;
        }
    }

    // Override id_1 for public input rows
    for &pub_row in &pub_input_rows {
        *polys.id_1.at_mut(pub_row) = separator + Fr::from(pub_row as u64);
    }

    // Initialize sigma to identity (no permutation)
    let sigma_polys = [
        &mut polys.sigma_1, &mut polys.sigma_2,
        &mut polys.sigma_3, &mut polys.sigma_4,
    ];
    for (wire_idx, sigma_poly) in sigma_polys.into_iter().enumerate() {
        let wire_offset = separator * Fr::from(wire_idx as u64);
        for i in 0..circuit_size {
            *sigma_poly.at_mut(i) = Fr::from(i as u64) + wire_offset;
        }
    }

    // Override sigma_1 for public input rows
    for &pub_row in &pub_input_rows {
        *polys.sigma_1.at_mut(pub_row) = Fr::zero() - Fr::from(pub_row as u64 + 1);
    }

    // Apply copy constraints
    let mut var_occurrences: BTreeMap<u32, Vec<(usize, usize)>> = BTreeMap::new();

    for block in builder.blocks.get() {
        let offset = block.block.trace_offset as usize;
        let block_size = block.size();

        for row in 0..block_size {
            let trace_row = offset + row;
            for wire_idx in 0..4 {
                if wire_idx == 0 && pub_input_rows.contains(&trace_row) {
                    continue;
                }
                let var_idx = block.block.wires[wire_idx][row];
                let real_idx = builder.base.real_variable_index[var_idx as usize];
                var_occurrences
                    .entry(real_idx)
                    .or_insert_with(Vec::new)
                    .push((wire_idx, trace_row));
            }
        }
    }

    // Create permutation cycles in sigma
    for (_real_idx, positions) in &var_occurrences {
        if positions.len() <= 1 {
            continue;
        }
        for i in 0..positions.len() {
            let (src_wire, src_row) = positions[i];
            let (dst_wire, dst_row) = positions[(i + 1) % positions.len()];
            let dst_value =
                Fr::from(dst_row as u64) + separator * Fr::from(dst_wire as u64);
            match src_wire {
                0 => *polys.sigma_1.at_mut(src_row) = dst_value,
                1 => *polys.sigma_2.at_mut(src_row) = dst_value,
                2 => *polys.sigma_3.at_mut(src_row) = dst_value,
                3 => *polys.sigma_4.at_mut(src_row) = dst_value,
                _ => unreachable!(),
            }
        }
    }
}

/// Populate table polynomials from the builder's lookup tables.
fn populate_tables(
    builder: &UltraCircuitBuilder<Bn254FrParams>,
    polys: &mut ProverPolynomials<Bn254FrParams>,
) {
    if builder.lookup_tables.is_empty() {
        return;
    }

    let total_content = builder.blocks.get_total_content_size();
    let table_offset = total_content + 1; // +1 for zero row

    let mut row = table_offset;
    for table in &builder.lookup_tables {
        for gate in &table.lookup_gates {
            if row < polys.table_1.size() {
                if gate.len() > 0 { *polys.table_1.at_mut(row) = gate[0]; }
                if gate.len() > 1 { *polys.table_2.at_mut(row) = gate[1]; }
                if gate.len() > 2 { *polys.table_3.at_mut(row) = gate[2]; }
                if gate.len() > 3 { *polys.table_4.at_mut(row) = gate[3]; }
            }
            row += 1;
        }
    }
}
