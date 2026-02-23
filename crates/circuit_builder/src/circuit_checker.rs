//! Port of `ultra_circuit_checker.hpp/.cpp` — validates circuit satisfaction.
//!
//! Provides `UltraCircuitChecker::check()` which evaluates all Ultra relations against
//! the execution trace to validate that a circuit's witness satisfies all constraints.
//!
//! This is the validation harness used by tests to verify that circuits built with
//! the `UltraCircuitBuilder` are correctly satisfied.

use std::collections::BTreeMap;

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use bbrs_honk::relation_checker::{self, AllSubrelationFailures};
use bbrs_honk::{compute_public_input_delta, PERMUTATION_ARGUMENT_VALUE_SEPARATOR};
use bbrs_relations::relation_parameters::RelationParameters;
use bbrs_relations::ultra::input_elements::{InputElements, NUM_ELEMENTS};
use bbrs_relations::ultra::{
    arithmetic, delta_range, elliptic, logderiv_lookup, memory, non_native_field, permutation,
    poseidon2_external, poseidon2_internal,
};

use crate::ultra_builder::UltraCircuitBuilder;

/// Number of polynomial columns tracked by the circuit checker.
/// Matches `InputElements::NUM_ELEMENTS` (45).
const NUM_POLYS: usize = NUM_ELEMENTS;

/// Internal polynomial storage for the circuit checker.
///
/// Each column is a `Vec<Field<P>>` of length `circuit_size`. The column indices
/// match `InputElements` accessor indices for direct mapping.
struct CheckerPolynomials<P: FieldParams> {
    data: Vec<Vec<Field<P>>>,
}

impl<P: FieldParams> CheckerPolynomials<P> {
    fn new(circuit_size: usize) -> Self {
        Self {
            data: vec![vec![Field::zero(); circuit_size]; NUM_POLYS],
        }
    }

    /// Extract an `InputElements` for the given row index.
    fn get_row(&self, row_idx: usize) -> InputElements<P> {
        let mut ie = InputElements {
            data: [Field::zero(); NUM_ELEMENTS],
        };
        for (col, poly) in self.data.iter().enumerate() {
            ie.data[col] = poly[row_idx];
        }
        ie
    }
}

/// Column indices matching `InputElements` accessor layout.
mod col {
    // Selectors
    pub const Q_C: usize = 0;
    pub const Q_L: usize = 1;
    pub const Q_R: usize = 2;
    pub const Q_O: usize = 3;
    pub const Q_4: usize = 4;
    pub const Q_M: usize = 5;
    pub const Q_ARITH: usize = 6;
    pub const Q_DELTA_RANGE: usize = 7;
    pub const Q_ELLIPTIC: usize = 8;
    pub const Q_MEMORY: usize = 9;
    pub const Q_NNF: usize = 10;
    pub const Q_LOOKUP: usize = 11;
    pub const Q_POSEIDON2_EXTERNAL: usize = 12;
    pub const Q_POSEIDON2_INTERNAL: usize = 13;
    // Sigmas
    pub const SIGMA_1: usize = 14;
    pub const SIGMA_2: usize = 15;
    pub const SIGMA_3: usize = 16;
    pub const SIGMA_4: usize = 17;
    // IDs
    pub const ID_1: usize = 18;
    pub const ID_2: usize = 19;
    pub const ID_3: usize = 20;
    pub const ID_4: usize = 21;
    // Tables
    pub const TABLE_1: usize = 22;
    pub const TABLE_2: usize = 23;
    pub const TABLE_3: usize = 24;
    pub const TABLE_4: usize = 25;
    // Lagrange
    pub const LAGRANGE_FIRST: usize = 26;
    pub const LAGRANGE_LAST: usize = 27;
    // Wires
    pub const W_L: usize = 28;
    pub const W_R: usize = 29;
    pub const W_O: usize = 30;
    pub const W_4: usize = 31;
    // Sorted accumulator (not used by Ultra, but present in InputElements)
    pub const SORTED_ACCUM: usize = 32;
    // Grand products
    pub const Z_PERM: usize = 33;
    // z_lookup (34) not used by Ultra log-deriv lookup
    // Shifted wires
    pub const W_L_SHIFT: usize = 35;
    pub const W_R_SHIFT: usize = 36;
    pub const W_O_SHIFT: usize = 37;
    pub const W_4_SHIFT: usize = 38;
    // Shifted sorted accumulator
    pub const SORTED_ACCUM_SHIFT: usize = 39;
    // Shifted grand products
    pub const Z_PERM_SHIFT: usize = 40;
    // z_lookup_shift (41) not used
    // Lookup-specific
    pub const LOOKUP_READ_COUNTS: usize = 42;
    pub const LOOKUP_READ_TAGS: usize = 43;
    pub const LOOKUP_INVERSES: usize = 44;
}

/// UltraCircuitChecker validates that a circuit's witness satisfies all constraints.
///
/// Port of C++ `UltraCircuitChecker`.
pub struct UltraCircuitChecker;

impl UltraCircuitChecker {
    /// Check that the given circuit builder satisfies all Ultra relations.
    ///
    /// Returns `Ok(())` if all relations are satisfied, or `Err(message)` with
    /// details about which relations failed at which rows.
    ///
    /// Port of C++ `UltraCircuitChecker::check()`.
    pub fn check<P: FieldParams>(builder: &mut UltraCircuitBuilder<P>) -> Result<(), String> {
        // Finalize the circuit if not already done
        if !builder.circuit_finalized {
            builder.finalize_circuit(true);
        }

        // Compute trace offsets
        builder.blocks.compute_offsets();

        // Compute circuit size (round up to next power of 2)
        let total_content = builder.blocks.get_total_content_size();
        let circuit_size = (total_content + 1).next_power_of_two(); // +1 for zero row

        // Populate polynomial columns from the execution trace
        let mut polys = CheckerPolynomials::<P>::new(circuit_size);

        // Step 1: Populate wire values and selectors from execution trace blocks
        Self::populate_wires_and_selectors(&builder, &mut polys);

        // Step 2: Populate identity and sigma polynomials for the permutation argument
        Self::populate_id_and_sigma(&builder, &mut polys, circuit_size);

        // Step 3: Set lagrange polynomials
        polys.data[col::LAGRANGE_FIRST][0] = Field::one();
        polys.data[col::LAGRANGE_LAST][circuit_size - 1] = Field::one();

        // Step 4: Populate table polynomials from lookup tables
        Self::populate_tables(&builder, &mut polys);

        // Step 4b: Set shifted polynomials early so lookup computations can use wire shifts.
        // Z_PERM_SHIFT will be updated again after compute_grand_product_perm.
        Self::set_shifted(&mut polys, circuit_size);

        // Step 5: Generate random challenges
        let beta = Field::<P>::random_element();
        let gamma = Field::<P>::random_element();
        let eta = Field::<P>::random_element();
        let eta_two = Field::<P>::random_element();
        let eta_three = Field::<P>::random_element();

        // Step 6: Compute public input delta
        let public_input_values: Vec<Field<P>> = builder
            .base
            .public_inputs()
            .iter()
            .map(|&idx| builder.base.get_variable(idx))
            .collect();
        let pub_inputs_offset = Field::<P>::from(builder.blocks.pub_inputs.block.trace_offset as u64);
        let public_input_delta =
            compute_public_input_delta::<P>(&public_input_values, beta, gamma, pub_inputs_offset);

        let beta_sqr = beta * beta;
        let beta_cube = beta_sqr * beta;

        let params = RelationParameters {
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

        // Step 7: Compute lookup read counts and tags
        Self::compute_lookup_read_counts(&builder, &mut polys, &params);

        // Step 8: Compute grand product z_perm
        Self::compute_grand_product_perm(&mut polys, &params, circuit_size);

        // Step 9: Compute lookup inverses
        Self::compute_lookup_inverses(&mut polys, &params, circuit_size);

        // Step 10: Re-set shifted polynomials now that Z_PERM is computed
        Self::set_shifted(&mut polys, circuit_size);

        // Step 11: Check all relations
        let failures = Self::check_all_relations(&polys, &params, circuit_size);

        if failures.is_empty() {
            Ok(())
        } else {
            let mut msg = String::from("Circuit check failed:\n");
            for (name, sub_failures) in &failures {
                for (sub_idx, row) in sub_failures {
                    msg.push_str(&format!(
                        "  {} subrelation {} failed at row {}\n",
                        name, sub_idx, row
                    ));
                }
            }
            Err(msg)
        }
    }

    /// Populate wire values and selector values from the execution trace blocks.
    fn populate_wires_and_selectors<P: FieldParams>(
        builder: &UltraCircuitBuilder<P>,
        polys: &mut CheckerPolynomials<P>,
    ) {
        for block in builder.blocks.get() {
            let offset = block.block.trace_offset as usize;
            let block_size = block.size();

            for row in 0..block_size {
                let trace_row = offset + row;

                // Wires: dereference variable indices to get field values
                polys.data[col::W_L][trace_row] =
                    builder.base.get_variable(block.block.wires[0][row]);
                polys.data[col::W_R][trace_row] =
                    builder.base.get_variable(block.block.wires[1][row]);
                polys.data[col::W_O][trace_row] =
                    builder.base.get_variable(block.block.wires[2][row]);
                polys.data[col::W_4][trace_row] =
                    builder.base.get_variable(block.block.wires[3][row]);

                // Non-gate selectors: q_m, q_c, q_1(q_l), q_2(q_r), q_3(q_o), q_4
                polys.data[col::Q_M][trace_row] = block.block.q_m().get(row);
                polys.data[col::Q_C][trace_row] = block.block.q_c().get(row);
                polys.data[col::Q_L][trace_row] = block.block.q_1().get(row);
                polys.data[col::Q_R][trace_row] = block.block.q_2().get(row);
                polys.data[col::Q_O][trace_row] = block.block.q_3().get(row);
                polys.data[col::Q_4][trace_row] = block.block.q_4().get(row);

                // Gate selectors
                polys.data[col::Q_ARITH][trace_row] = block.q_arith().get(row);
                polys.data[col::Q_LOOKUP][trace_row] = block.q_lookup().get(row);
                polys.data[col::Q_DELTA_RANGE][trace_row] = block.q_delta_range().get(row);
                polys.data[col::Q_ELLIPTIC][trace_row] = block.q_elliptic().get(row);
                polys.data[col::Q_MEMORY][trace_row] = block.q_memory().get(row);
                polys.data[col::Q_NNF][trace_row] = block.q_nnf().get(row);
                polys.data[col::Q_POSEIDON2_EXTERNAL][trace_row] =
                    block.q_poseidon2_external().get(row);
                polys.data[col::Q_POSEIDON2_INTERNAL][trace_row] =
                    block.q_poseidon2_internal().get(row);
            }
        }
    }

    /// Populate identity and sigma polynomials.
    ///
    /// Identity polynomials encode the "expected" wire positions:
    ///   id_k[i] = i + k * n   (where n = separator, k = wire index)
    ///
    /// Sigma polynomials encode the copy constraints (permutation):
    ///   sigma_k[i] maps wire (k, i) to its copy target.
    ///
    /// Public input rows get special treatment:
    ///   - id_1[pub_row] = separator + pub_row  (not just pub_row)
    ///   - sigma_1[pub_row] = -(pub_row + 1)    (special virtual position)
    fn populate_id_and_sigma<P: FieldParams>(
        builder: &UltraCircuitBuilder<P>,
        polys: &mut CheckerPolynomials<P>,
        circuit_size: usize,
    ) {
        let separator = Field::<P>::from(PERMUTATION_ARGUMENT_VALUE_SEPARATOR);
        let id_cols = [col::ID_1, col::ID_2, col::ID_3, col::ID_4];
        let sigma_cols = [col::SIGMA_1, col::SIGMA_2, col::SIGMA_3, col::SIGMA_4];

        // Collect public input row positions for special handling
        let pub_inputs_offset = builder.blocks.pub_inputs.block.trace_offset as usize;
        let pub_inputs_size = builder.blocks.pub_inputs.size();
        let pub_input_rows: Vec<usize> = (0..pub_inputs_size)
            .map(|i| pub_inputs_offset + i)
            .collect();

        // Initialize identity polynomials: id_k[i] = i + k * separator
        for (wire_idx, &id_col) in id_cols.iter().enumerate() {
            let wire_offset = separator * Field::from(wire_idx as u64);
            for i in 0..circuit_size {
                polys.data[id_col][i] = Field::from(i as u64) + wire_offset;
            }
        }

        // Override id_1 for public input rows: id_1[pub_row] = separator + pub_row
        // This matches the numerator in compute_public_input_delta
        for &pub_row in &pub_input_rows {
            polys.data[col::ID_1][pub_row] = separator + Field::from(pub_row as u64);
        }

        // Initialize sigma to identity (no permutation)
        for (wire_idx, &sigma_col) in sigma_cols.iter().enumerate() {
            let wire_offset = separator * Field::from(wire_idx as u64);
            for i in 0..circuit_size {
                polys.data[sigma_col][i] = Field::from(i as u64) + wire_offset;
            }
        }

        // Override sigma_1 for public input rows: sigma_1[pub_row] = -(pub_row + 1)
        // This matches the denominator in compute_public_input_delta
        for &pub_row in &pub_input_rows {
            polys.data[col::SIGMA_1][pub_row] =
                Field::<P>::zero() - Field::from(pub_row as u64 + 1);
        }

        // Apply copy constraints: for each equivalence class, create a cycle in sigma.
        // Exclude public input wire-1 positions from the cycle since they have special sigma.
        let mut var_occurrences: BTreeMap<u32, Vec<(usize, usize)>> = BTreeMap::new();

        for block in builder.blocks.get() {
            let offset = block.block.trace_offset as usize;
            let block_size = block.size();

            for row in 0..block_size {
                let trace_row = offset + row;
                for wire_idx in 0..4 {
                    // Skip wire 0 at public input rows — handled by special sigma above
                    if wire_idx == 0 && pub_input_rows.contains(&trace_row) {
                        continue;
                    }
                    let var_idx = block.block.wires[wire_idx][row];
                    let real_idx = builder.base.real_variable_index[var_idx as usize];
                    var_occurrences
                        .entry(real_idx)
                        .or_default()
                        .push((wire_idx, trace_row));
                }
            }
        }

        // Create permutation cycles in sigma
        for (_real_idx, positions) in &var_occurrences {
            if positions.len() <= 1 {
                continue; // No copy constraint needed for singletons
            }

            // Create a cycle: position[0] → position[1] → ... → position[n-1] → position[0]
            for i in 0..positions.len() {
                let (src_wire, src_row) = positions[i];
                let (dst_wire, dst_row) = positions[(i + 1) % positions.len()];

                let dst_value =
                    Field::<P>::from(dst_row as u64) + separator * Field::from(dst_wire as u64);
                polys.data[sigma_cols[src_wire]][src_row] = dst_value;
            }
        }
    }

    /// Populate table polynomials from the builder's lookup tables.
    ///
    /// The lookup gate entries from each BasicTable are placed into the table_1..4
    /// polynomial columns at positions after the gate blocks.
    fn populate_tables<P: FieldParams>(
        builder: &UltraCircuitBuilder<P>,
        polys: &mut CheckerPolynomials<P>,
    ) {
        if builder.lookup_tables.is_empty() {
            return;
        }

        // Collect all table entries from all lookup tables.
        // Each lookup_gate entry is a Vec<Field<P>> with 3 column values,
        // and the 4th column (TABLE_4) is the table_index from the BasicTable.
        let mut table_entries: Vec<[Field<P>; 4]> = Vec::new();
        for table in &builder.lookup_tables {
            let table_idx = Field::<P>::from(table.table_index);
            for gate in &table.lookup_gates {
                let mut entry = [Field::<P>::zero(); 4];
                for (col, val) in gate.iter().enumerate() {
                    if col < 3 {
                        entry[col] = *val;
                    }
                }
                entry[3] = table_idx;
                table_entries.push(entry);
            }
        }

        // Place table entries in the table polynomials starting after the last block
        let total_content = builder.blocks.get_total_content_size();
        let table_offset = total_content + 1; // +1 for zero row

        for (i, entry) in table_entries.iter().enumerate() {
            let row = table_offset + i;
            if row < polys.data[col::TABLE_1].len() {
                polys.data[col::TABLE_1][row] = entry[0];
                polys.data[col::TABLE_2][row] = entry[1];
                polys.data[col::TABLE_3][row] = entry[2];
                polys.data[col::TABLE_4][row] = entry[3];
            }
        }
    }

    /// Compute lookup read counts and tags.
    ///
    /// For each lookup gate in the trace, find its matching table entry and
    /// increment the read count at that table row.
    fn compute_lookup_read_counts<P: FieldParams>(
        builder: &UltraCircuitBuilder<P>,
        polys: &mut CheckerPolynomials<P>,
        params: &RelationParameters<Field<P>>,
    ) {
        let lookup_block = &builder.blocks.lookup;
        let block_offset = lookup_block.block.trace_offset as usize;
        let block_size = lookup_block.size();

        if block_size == 0 {
            return;
        }

        // Build a map from table entry hash to table row position
        let total_content = builder.blocks.get_total_content_size();
        let table_offset = total_content + 1;

        // Count total table entries
        let total_table_entries: usize = builder
            .lookup_tables
            .iter()
            .map(|t| t.lookup_gates.len())
            .sum();

        // For each table row, compute the combined value and map to position
        let mut table_value_to_row: BTreeMap<[u64; 4], usize> = BTreeMap::new();
        for i in 0..total_table_entries {
            let row = table_offset + i;
            if row < polys.data[col::TABLE_1].len() {
                let key = polys.data[col::TABLE_1][row].data;
                table_value_to_row.entry(key).or_insert(row);
            }
        }

        // For each lookup gate, find its matching table entry by computing the
        // combined lookup value and matching against table entries
        for row in 0..block_size {
            let trace_row = block_offset + row;
            let q_lookup = polys.data[col::Q_LOOKUP][trace_row];

            if q_lookup.is_zero() {
                continue;
            }

            // Derive table entries from accumulator wires and step-size selectors:
            //   derived_entry_i = w_i + negative_step_size_i * w_i_shift
            let w_1 = polys.data[col::W_L][trace_row];
            let w_2 = polys.data[col::W_R][trace_row];
            let w_3 = polys.data[col::W_O][trace_row];
            let w_1_shift = polys.data[col::W_L_SHIFT][trace_row];
            let w_2_shift = polys.data[col::W_R_SHIFT][trace_row];
            let w_3_shift = polys.data[col::W_O_SHIFT][trace_row];

            let neg_step_1 = polys.data[col::Q_R][trace_row]; // q_2 = -column_1_step_size
            let neg_step_2 = polys.data[col::Q_M][trace_row]; // q_m = -column_2_step_size
            let neg_step_3 = polys.data[col::Q_C][trace_row]; // q_c = -column_3_step_size
            let table_index = polys.data[col::Q_O][trace_row]; // q_3 = table_index

            let derived_1 = w_1 + neg_step_1 * w_1_shift + params.gamma;
            let derived_2 = w_2 + neg_step_2 * w_2_shift;
            let derived_3 = w_3 + neg_step_3 * w_3_shift;

            let combined = derived_1
                + derived_2 * params.eta
                + derived_3 * params.eta_two
                + table_index * params.eta_three;

            // Search through table entries for a match
            for i in 0..total_table_entries {
                let trow = table_offset + i;
                if trow >= polys.data[col::TABLE_1].len() {
                    break;
                }
                let t1 = polys.data[col::TABLE_1][trow];
                let t2 = polys.data[col::TABLE_2][trow];
                let t3 = polys.data[col::TABLE_3][trow];
                let t4 = polys.data[col::TABLE_4][trow];

                let table_combined =
                    t1 + params.gamma + t2 * params.eta + t3 * params.eta_two + t4 * params.eta_three;

                if combined == table_combined {
                    // Increment read count at this table row
                    polys.data[col::LOOKUP_READ_COUNTS][trow] =
                        polys.data[col::LOOKUP_READ_COUNTS][trow] + Field::one();
                    polys.data[col::LOOKUP_READ_TAGS][trow] = Field::one();
                    break;
                }
            }
        }
    }

    /// Compute the grand product polynomial z_perm for the permutation argument.
    fn compute_grand_product_perm<P: FieldParams>(
        polys: &mut CheckerPolynomials<P>,
        params: &RelationParameters<Field<P>>,
        circuit_size: usize,
    ) {
        if circuit_size <= 1 {
            return;
        }

        let beta = params.beta;
        let gamma = params.gamma;

        let iteration_size = circuit_size - 1;

        // Compute numerator and denominator at each row
        let mut numerators = vec![Field::<P>::one(); iteration_size];
        let mut denominators = vec![Field::<P>::one(); iteration_size];

        for i in 0..iteration_size {
            // Numerator: ∏_k (w_k + id_k * beta + gamma)
            let num = (polys.data[col::W_L][i] + polys.data[col::ID_1][i] * beta + gamma)
                * (polys.data[col::W_R][i] + polys.data[col::ID_2][i] * beta + gamma)
                * (polys.data[col::W_O][i] + polys.data[col::ID_3][i] * beta + gamma)
                * (polys.data[col::W_4][i] + polys.data[col::ID_4][i] * beta + gamma);

            // Denominator: ∏_k (w_k + sigma_k * beta + gamma)
            let den = (polys.data[col::W_L][i] + polys.data[col::SIGMA_1][i] * beta + gamma)
                * (polys.data[col::W_R][i] + polys.data[col::SIGMA_2][i] * beta + gamma)
                * (polys.data[col::W_O][i] + polys.data[col::SIGMA_3][i] * beta + gamma)
                * (polys.data[col::W_4][i] + polys.data[col::SIGMA_4][i] * beta + gamma);

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

        // Z[0] = 0 (reserved zero row for shiftable polynomial)
        // Z[i+1] = running_numerator[i] / running_denominator[i]
        for i in 0..iteration_size {
            polys.data[col::Z_PERM][i + 1] = numerators[i] * denominators[i];
        }
    }

    /// Compute lookup inverse polynomial.
    fn compute_lookup_inverses<P: FieldParams>(
        polys: &mut CheckerPolynomials<P>,
        params: &RelationParameters<Field<P>>,
        circuit_size: usize,
    ) {
        let mut denominators = vec![Field::<P>::zero(); circuit_size];
        let mut has_inverse = vec![false; circuit_size];

        for i in 0..circuit_size {
            let q_lookup = polys.data[col::Q_LOOKUP][i];
            let lookup_read_counts = polys.data[col::LOOKUP_READ_COUNTS][i];

            // An operation exists if q_lookup != 0 or lookup_read_counts != 0
            if q_lookup.is_zero() && lookup_read_counts.is_zero() {
                continue;
            }

            has_inverse[i] = true;

            // read_term denominator: derived entries from accumulator wires + step sizes
            let w_1 = polys.data[col::W_L][i];
            let w_2 = polys.data[col::W_R][i];
            let w_3 = polys.data[col::W_O][i];
            let w_1_shift = polys.data[col::W_L_SHIFT][i];
            let w_2_shift = polys.data[col::W_R_SHIFT][i];
            let w_3_shift = polys.data[col::W_O_SHIFT][i];

            let neg_step_1 = polys.data[col::Q_R][i]; // q_2
            let neg_step_2 = polys.data[col::Q_M][i]; // q_m
            let neg_step_3 = polys.data[col::Q_C][i]; // q_c
            let table_index = polys.data[col::Q_O][i]; // q_3

            let derived_1 = w_1 + neg_step_1 * w_1_shift + params.gamma;
            let derived_2 = w_2 + neg_step_2 * w_2_shift;
            let derived_3 = w_3 + neg_step_3 * w_3_shift;
            let read_denom = derived_1
                + derived_2 * params.eta
                + derived_3 * params.eta_two
                + table_index * params.eta_three;

            // write_term denominator: table_1 + gamma + table_2*eta + table_3*eta^2 + table_4*eta^3
            let t1 = polys.data[col::TABLE_1][i];
            let t2 = polys.data[col::TABLE_2][i];
            let t3 = polys.data[col::TABLE_3][i];
            let t4 = polys.data[col::TABLE_4][i];
            let write_denom =
                t1 + params.gamma + t2 * params.eta + t3 * params.eta_two + t4 * params.eta_three;

            denominators[i] = read_denom * write_denom;
        }

        // Batch invert non-zero denominators
        batch_invert_nonzero(&mut denominators);

        // Write results
        for i in 0..circuit_size {
            if has_inverse[i] {
                polys.data[col::LOOKUP_INVERSES][i] = denominators[i];
            }
        }
    }

    /// Set shifted polynomial columns from their unshifted counterparts.
    fn set_shifted<P: FieldParams>(polys: &mut CheckerPolynomials<P>, circuit_size: usize) {
        // Shifted values: shifted[i] = unshifted[i+1]
        let shift_pairs = [
            (col::W_L, col::W_L_SHIFT),
            (col::W_R, col::W_R_SHIFT),
            (col::W_O, col::W_O_SHIFT),
            (col::W_4, col::W_4_SHIFT),
            (col::Z_PERM, col::Z_PERM_SHIFT),
            (col::SORTED_ACCUM, col::SORTED_ACCUM_SHIFT),
        ];

        for &(src, dst) in &shift_pairs {
            for i in 0..circuit_size - 1 {
                polys.data[dst][i] = polys.data[src][i + 1];
            }
            // Last row's shift wraps to zero (or stays zero)
        }
    }

    /// Check all 9 Ultra relations against the populated polynomials.
    fn check_all_relations<P: FieldParams>(
        polys: &CheckerPolynomials<P>,
        params: &RelationParameters<Field<P>>,
        circuit_size: usize,
    ) -> AllSubrelationFailures {
        let mut all_failures = AllSubrelationFailures::new();

        // Helper: check a single relation
        macro_rules! check_relation {
            ($name:expr, $mod:ident, $has_lin_dep:expr) => {{
                let linearly_independent: Vec<bool> = $mod::SUBRELATION_PARTIAL_LENGTHS
                    .iter()
                    .enumerate()
                    .map(|(i, _)| {
                        // All are linearly independent except LogDerivLookup subrelation 1
                        if $name == "LogDerivLookup" && i == 1 {
                            false
                        } else {
                            true
                        }
                    })
                    .collect();

                let failures = relation_checker::check_relation(
                    |evals: &mut Vec<Field<P>>,
                     row: &InputElements<P>,
                     params: &RelationParameters<Field<P>>,
                     scaling_factor: &Field<P>| {
                        // Convert Vec to fixed-size array for the relation accumulate function
                        let mut fixed_evals = [Field::<P>::zero(); $mod::NUM_SUBRELATIONS];
                        for (i, e) in evals.iter().enumerate().take($mod::NUM_SUBRELATIONS) {
                            fixed_evals[i] = *e;
                        }
                        $mod::accumulate(&mut fixed_evals, row, params, scaling_factor);
                        for (i, e) in fixed_evals.iter().enumerate() {
                            evals[i] = *e;
                        }
                    },
                    |i: usize| polys.get_row(i),
                    $mod::NUM_SUBRELATIONS,
                    &linearly_independent,
                    circuit_size,
                    params,
                    $has_lin_dep,
                );

                if !failures.is_empty() {
                    all_failures.insert($name.to_string(), failures);
                }
            }};
        }

        check_relation!("Arithmetic", arithmetic, false);
        check_relation!("Permutation", permutation, false);
        check_relation!("LogDerivLookup", logderiv_lookup, true);
        check_relation!("DeltaRange", delta_range, false);
        check_relation!("Elliptic", elliptic, false);
        check_relation!("Memory", memory, false);
        check_relation!("NonNativeField", non_native_field, false);
        check_relation!("Poseidon2External", poseidon2_external, false);
        check_relation!("Poseidon2Internal", poseidon2_internal, false);

        all_failures
    }
}

/// Batch-invert a slice of field elements in place.
fn batch_invert_in_place<P: FieldParams>(values: &mut [Field<P>]) {
    if values.is_empty() {
        return;
    }

    let n = values.len();
    let mut scratch = vec![Field::<P>::one(); n];
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

/// Batch-invert a slice of field elements in place, skipping zeros.
fn batch_invert_nonzero<P: FieldParams>(values: &mut [Field<P>]) {
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
    let mut scratch = vec![Field::<P>::one(); n];
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

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::Bn254FrParams;
    use bbrs_ecc::fields::field::Field;

    type P = Bn254FrParams;
    type Fr = Field<P>;

    /// Test that an empty circuit (with only the zero constant) passes the checker.
    #[test]
    fn test_empty_circuit_passes() {
        let mut builder = UltraCircuitBuilder::<P>::new();
        let result = UltraCircuitChecker::check(&mut builder);
        assert!(result.is_ok(), "Empty circuit should pass: {:?}", result);
    }

    /// Test a simple addition gate: 1 + 1 = 2.
    #[test]
    fn test_simple_addition_gate() {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let a_val = Fr::one();
        let b_val = Fr::one();
        let c_val = Fr::from(2u64);

        let a = builder.base.add_variable(a_val);
        let b = builder.base.add_variable(b_val);
        let c = builder.base.add_variable(c_val);

        use crate::gate_data::AddTriple;
        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::one(),
            b_scaling: Fr::one(),
            c_scaling: Fr::zero() - Fr::one(),
            const_scaling: Fr::zero(),
        });

        let result = UltraCircuitChecker::check(&mut builder);
        assert!(
            result.is_ok(),
            "Simple addition gate should pass: {:?}",
            result
        );
    }

    /// Test a multiplication gate: 3 * 7 = 21.
    #[test]
    fn test_simple_multiplication_gate() {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let a_val = Fr::from(3u64);
        let b_val = Fr::from(7u64);
        let c_val = Fr::from(21u64);

        let a = builder.base.add_variable(a_val);
        let b = builder.base.add_variable(b_val);
        let c = builder.base.add_variable(c_val);

        use crate::gate_data::MulQuad;
        builder.create_big_mul_add_gate(
            &MulQuad {
                a,
                b,
                c,
                d: builder.base.zero_idx(),
                mul_scaling: Fr::one(),
                a_scaling: Fr::zero(),
                b_scaling: Fr::zero(),
                c_scaling: Fr::zero() - Fr::one(),
                d_scaling: Fr::zero(),
                const_scaling: Fr::zero(),
            },
            false,
        );

        let result = UltraCircuitChecker::check(&mut builder);
        assert!(
            result.is_ok(),
            "Simple multiplication gate should pass: {:?}",
            result
        );
    }

    /// Test that an incorrect witness fails the checker.
    #[test]
    fn test_bad_witness_fails() {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let a = builder.base.add_variable(Fr::from(3u64));
        let b = builder.base.add_variable(Fr::from(7u64));
        let c = builder.base.add_variable(Fr::from(22u64)); // Wrong! Should be 21.

        use crate::gate_data::AddTriple;
        builder.create_add_gate(&AddTriple {
            a,
            b,
            c,
            a_scaling: Fr::from(7u64),  // a * 7 = 21
            b_scaling: Fr::zero(),
            c_scaling: Fr::zero() - Fr::one(), // -c
            const_scaling: Fr::zero(),
        });

        let result = UltraCircuitChecker::check(&mut builder);
        assert!(result.is_err(), "Bad witness should fail the checker");
    }

    /// Test a circuit with copy constraints (assert_equal).
    #[test]
    fn test_copy_constraints() {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let val = Fr::from(42u64);
        let a = builder.base.add_variable(val);
        let b = builder.base.add_variable(val);

        // Assert they are equal (creates copy constraint)
        builder.base.assert_equal(a, b, "test copy");

        // Use both in arithmetic gates
        use crate::gate_data::AddTriple;
        builder.create_add_gate(&AddTriple {
            a,
            b: builder.base.zero_idx(),
            c: builder.base.zero_idx(),
            a_scaling: Fr::one(),
            b_scaling: Fr::zero(),
            c_scaling: Fr::zero(),
            const_scaling: Fr::zero() - val,
        });

        builder.create_add_gate(&AddTriple {
            a: b,
            b: builder.base.zero_idx(),
            c: builder.base.zero_idx(),
            a_scaling: Fr::one(),
            b_scaling: Fr::zero(),
            c_scaling: Fr::zero(),
            const_scaling: Fr::zero() - val,
        });

        let result = UltraCircuitChecker::check(&mut builder);
        assert!(
            result.is_ok(),
            "Circuit with copy constraints should pass: {:?}",
            result
        );
    }

    /// Test a circuit with public inputs.
    #[test]
    fn test_public_inputs() {
        let mut builder = UltraCircuitBuilder::<P>::new();

        let val = Fr::from(99u64);
        let pi_idx = builder.base.add_public_variable(val);

        use crate::gate_data::AddTriple;
        builder.create_add_gate(&AddTriple {
            a: pi_idx,
            b: builder.base.zero_idx(),
            c: builder.base.zero_idx(),
            a_scaling: Fr::one(),
            b_scaling: Fr::zero(),
            c_scaling: Fr::zero(),
            const_scaling: Fr::zero() - val,
        });

        let result = UltraCircuitChecker::check(&mut builder);
        assert!(
            result.is_ok(),
            "Circuit with public inputs should pass: {:?}",
            result
        );
    }
}
