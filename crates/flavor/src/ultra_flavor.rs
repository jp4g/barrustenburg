//! Port of `ultra_flavor.hpp` — Ultra Honk Flavor.
//!
//! Defines the UltraFlavor entity types, prover polynomials, and computed constants
//! for the Ultra Honk proving system over BN254.
//!
//! Entity counts:
//! - Precomputed: 28 (selectors, sigmas, ids, tables, lagrange)
//! - Witness: 8 (wires, z_perm, lookup_inverses, lookup_read_counts, lookup_read_tags)
//! - Shifted: 5 (w_l_shift, w_r_shift, w_o_shift, w_4_shift, z_perm_shift)
//! - Total: 41

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_polynomials::polynomial::Polynomial;
use bbrs_polynomials::univariate::Univariate;

// ============================================================================
// Constants
// ============================================================================

pub const NUM_WIRES: usize = 4;
pub const NUM_PRECOMPUTED_ENTITIES: usize = 28;
pub const NUM_WITNESS_ENTITIES: usize = 8;
pub const NUM_SHIFTED_ENTITIES: usize = 5;
pub const NUM_ALL_ENTITIES: usize =
    NUM_PRECOMPUTED_ENTITIES + NUM_WITNESS_ENTITIES + NUM_SHIFTED_ENTITIES; // = 41

/// Number of relations in the UltraFlavor sumcheck.
///
/// 9 relations: Arithmetic, Permutation, LogDerivLookup, DeltaRange,
/// Elliptic, Memory, NonNativeField, Poseidon2External, Poseidon2Internal
pub const NUM_RELATIONS: usize = 9;

/// Total number of subrelations across all relations.
///
/// Arithmetic(2) + Permutation(2) + LogDerivLookup(3) + DeltaRange(4) +
/// Elliptic(2) + Memory(6) + NonNativeField(1) + Poseidon2External(4) + Poseidon2Internal(4) = 28
pub const NUM_SUBRELATIONS: usize = 28;

/// Maximum partial relation length across all subrelations.
/// From Poseidon2External/Internal which have length 7.
pub const MAX_PARTIAL_RELATION_LENGTH: usize = 7;

/// Batched relation partial length = MAX_PARTIAL_RELATION_LENGTH + 1 (non-ZK).
pub const BATCHED_RELATION_PARTIAL_LENGTH: usize = MAX_PARTIAL_RELATION_LENGTH + 1; // = 8

/// Whether this flavor uses short monomials.
pub const USE_SHORT_MONOMIALS: bool = true;

/// Whether this flavor has zero-knowledge.
pub const HAS_ZK: bool = false;

/// Whether the first row is reserved for zeros to enable shifts.
pub const HAS_ZERO_ROW: bool = true;

/// All subrelation partial lengths in flat order across all 9 relations.
pub const ALL_SUBRELATION_PARTIAL_LENGTHS: [usize; NUM_SUBRELATIONS] = [
    // Arithmetic [6, 5]
    6, 5, // Permutation [6, 3]
    6, 3, // LogDerivLookup [5, 5, 3]
    5, 5, 3, // DeltaRange [6, 6, 6, 6]
    6, 6, 6, 6, // Elliptic [6, 6]
    6, 6, // Memory [6, 6, 6, 6, 6, 6]
    6, 6, 6, 6, 6, 6, // NonNativeField [6]
    6, // Poseidon2External [7, 7, 7, 7]
    7, 7, 7, 7, // Poseidon2Internal [7, 7, 7, 7]
    7, 7, 7, 7,
];

/// Subrelation linear independence flags in flat order.
/// LogDerivLookup subrelation [1] (the lookup sum) is linearly dependent.
pub const ALL_SUBRELATION_LINEARLY_INDEPENDENT: [bool; NUM_SUBRELATIONS] = [
    // Arithmetic
    true, true, // Permutation
    true, true, // LogDerivLookup: inverse check (true), lookup sum (false), boolean check (true)
    true, false, true, // DeltaRange
    true, true, true, true, // Elliptic
    true, true, // Memory
    true, true, true, true, true, true, // NonNativeField
    true, // Poseidon2External
    true, true, true, true, // Poseidon2Internal
    true, true, true, true,
];

/// The separator value used in the permutation argument.
pub const PERMUTATION_ARGUMENT_VALUE_SEPARATOR: u64 = 1 << 28;

// ============================================================================
// AllEntities — combined flat view of all 41 entity fields
// ============================================================================

/// All entities for the UltraFlavor (41 total).
///
/// Layout: precomputed (28), witness (8), shifted (5).
pub struct AllEntities<T> {
    // Precomputed (28)
    pub q_m: T,
    pub q_c: T,
    pub q_l: T,
    pub q_r: T,
    pub q_o: T,
    pub q_4: T,
    pub q_lookup: T,
    pub q_arith: T,
    pub q_delta_range: T,
    pub q_elliptic: T,
    pub q_memory: T,
    pub q_nnf: T,
    pub q_poseidon2_external: T,
    pub q_poseidon2_internal: T,
    pub sigma_1: T,
    pub sigma_2: T,
    pub sigma_3: T,
    pub sigma_4: T,
    pub id_1: T,
    pub id_2: T,
    pub id_3: T,
    pub id_4: T,
    pub table_1: T,
    pub table_2: T,
    pub table_3: T,
    pub table_4: T,
    pub lagrange_first: T,
    pub lagrange_last: T,
    // Witness (8)
    pub w_l: T,
    pub w_r: T,
    pub w_o: T,
    pub w_4: T,
    pub z_perm: T,
    pub lookup_inverses: T,
    pub lookup_read_counts: T,
    pub lookup_read_tags: T,
    // Shifted (5)
    pub w_l_shift: T,
    pub w_r_shift: T,
    pub w_o_shift: T,
    pub w_4_shift: T,
    pub z_perm_shift: T,
}

impl<T> AllEntities<T> {
    pub fn get_all(&self) -> [&T; NUM_ALL_ENTITIES] {
        [
            &self.q_m,
            &self.q_c,
            &self.q_l,
            &self.q_r,
            &self.q_o,
            &self.q_4,
            &self.q_lookup,
            &self.q_arith,
            &self.q_delta_range,
            &self.q_elliptic,
            &self.q_memory,
            &self.q_nnf,
            &self.q_poseidon2_external,
            &self.q_poseidon2_internal,
            &self.sigma_1,
            &self.sigma_2,
            &self.sigma_3,
            &self.sigma_4,
            &self.id_1,
            &self.id_2,
            &self.id_3,
            &self.id_4,
            &self.table_1,
            &self.table_2,
            &self.table_3,
            &self.table_4,
            &self.lagrange_first,
            &self.lagrange_last,
            &self.w_l,
            &self.w_r,
            &self.w_o,
            &self.w_4,
            &self.z_perm,
            &self.lookup_inverses,
            &self.lookup_read_counts,
            &self.lookup_read_tags,
            &self.w_l_shift,
            &self.w_r_shift,
            &self.w_o_shift,
            &self.w_4_shift,
            &self.z_perm_shift,
        ]
    }

    pub fn get_all_mut(&mut self) -> [&mut T; NUM_ALL_ENTITIES] {
        [
            &mut self.q_m,
            &mut self.q_c,
            &mut self.q_l,
            &mut self.q_r,
            &mut self.q_o,
            &mut self.q_4,
            &mut self.q_lookup,
            &mut self.q_arith,
            &mut self.q_delta_range,
            &mut self.q_elliptic,
            &mut self.q_memory,
            &mut self.q_nnf,
            &mut self.q_poseidon2_external,
            &mut self.q_poseidon2_internal,
            &mut self.sigma_1,
            &mut self.sigma_2,
            &mut self.sigma_3,
            &mut self.sigma_4,
            &mut self.id_1,
            &mut self.id_2,
            &mut self.id_3,
            &mut self.id_4,
            &mut self.table_1,
            &mut self.table_2,
            &mut self.table_3,
            &mut self.table_4,
            &mut self.lagrange_first,
            &mut self.lagrange_last,
            &mut self.w_l,
            &mut self.w_r,
            &mut self.w_o,
            &mut self.w_4,
            &mut self.z_perm,
            &mut self.lookup_inverses,
            &mut self.lookup_read_counts,
            &mut self.lookup_read_tags,
            &mut self.w_l_shift,
            &mut self.w_r_shift,
            &mut self.w_o_shift,
            &mut self.w_4_shift,
            &mut self.z_perm_shift,
        ]
    }

    pub fn get_precomputed(&self) -> [&T; NUM_PRECOMPUTED_ENTITIES] {
        [
            &self.q_m, &self.q_c, &self.q_l, &self.q_r, &self.q_o, &self.q_4,
            &self.q_lookup, &self.q_arith, &self.q_delta_range, &self.q_elliptic,
            &self.q_memory, &self.q_nnf, &self.q_poseidon2_external, &self.q_poseidon2_internal,
            &self.sigma_1, &self.sigma_2, &self.sigma_3, &self.sigma_4,
            &self.id_1, &self.id_2, &self.id_3, &self.id_4,
            &self.table_1, &self.table_2, &self.table_3, &self.table_4,
            &self.lagrange_first, &self.lagrange_last,
        ]
    }

    pub fn get_witness(&self) -> [&T; NUM_WITNESS_ENTITIES] {
        [
            &self.w_l, &self.w_r, &self.w_o, &self.w_4,
            &self.z_perm, &self.lookup_inverses, &self.lookup_read_counts, &self.lookup_read_tags,
        ]
    }

    pub fn get_shifted(&self) -> [&T; NUM_SHIFTED_ENTITIES] {
        [
            &self.w_l_shift, &self.w_r_shift, &self.w_o_shift,
            &self.w_4_shift, &self.z_perm_shift,
        ]
    }

    pub fn get_wires(&self) -> [&T; NUM_WIRES] {
        [&self.w_l, &self.w_r, &self.w_o, &self.w_4]
    }

    pub fn get_sigmas(&self) -> [&T; 4] {
        [&self.sigma_1, &self.sigma_2, &self.sigma_3, &self.sigma_4]
    }

    pub fn get_ids(&self) -> [&T; 4] {
        [&self.id_1, &self.id_2, &self.id_3, &self.id_4]
    }

    pub fn get_tables(&self) -> [&T; 4] {
        [&self.table_1, &self.table_2, &self.table_3, &self.table_4]
    }

    /// Get references to polynomials that will be shifted.
    /// These are the 5 witness polynomials: w_l, w_r, w_o, w_4, z_perm.
    pub fn get_to_be_shifted(&self) -> [&T; NUM_SHIFTED_ENTITIES] {
        [&self.w_l, &self.w_r, &self.w_o, &self.w_4, &self.z_perm]
    }

    pub fn get_labels() -> &'static [&'static str; NUM_ALL_ENTITIES] {
        &[
            "q_m", "q_c", "q_l", "q_r", "q_o", "q_4",
            "q_lookup", "q_arith", "q_delta_range", "q_elliptic",
            "q_memory", "q_nnf", "q_poseidon2_external", "q_poseidon2_internal",
            "sigma_1", "sigma_2", "sigma_3", "sigma_4",
            "id_1", "id_2", "id_3", "id_4",
            "table_1", "table_2", "table_3", "table_4",
            "lagrange_first", "lagrange_last",
            "w_l", "w_r", "w_o", "w_4",
            "z_perm", "lookup_inverses", "lookup_read_counts", "lookup_read_tags",
            "w_l_shift", "w_r_shift", "w_o_shift", "w_4_shift", "z_perm_shift",
        ]
    }
}

// ============================================================================
// AllValues — AllEntities<Field<P>>
// ============================================================================

pub type AllValues<P> = AllEntities<Field<P>>;

impl<P: FieldParams> AllValues<P> {
    pub fn zero() -> Self {
        let z = Field::zero();
        AllEntities {
            q_m: z, q_c: z, q_l: z, q_r: z, q_o: z, q_4: z,
            q_lookup: z, q_arith: z, q_delta_range: z, q_elliptic: z,
            q_memory: z, q_nnf: z, q_poseidon2_external: z, q_poseidon2_internal: z,
            sigma_1: z, sigma_2: z, sigma_3: z, sigma_4: z,
            id_1: z, id_2: z, id_3: z, id_4: z,
            table_1: z, table_2: z, table_3: z, table_4: z,
            lagrange_first: z, lagrange_last: z,
            w_l: z, w_r: z, w_o: z, w_4: z,
            z_perm: z, lookup_inverses: z, lookup_read_counts: z, lookup_read_tags: z,
            w_l_shift: z, w_r_shift: z, w_o_shift: z, w_4_shift: z, z_perm_shift: z,
        }
    }
}

// ============================================================================
// ProverPolynomials — AllEntities<Polynomial<P>>
// ============================================================================

pub type ProverPolynomials<P> = AllEntities<Polynomial<P>>;

impl<P: FieldParams> ProverPolynomials<P> {
    /// Create prover polynomials of the given circuit size, all initialized to zero.
    ///
    /// Precomputed polynomials are regular (start_index = 0).
    /// Witness polynomials that need shifting are shiftable (start_index = 1).
    /// Other witness polynomials (lookup_inverses, lookup_read_counts, lookup_read_tags) are regular.
    pub fn new(circuit_size: usize) -> Self {
        let reg = || Polynomial::new(circuit_size, circuit_size, 0);
        let shift = || Polynomial::shiftable(circuit_size);
        let empty = || Polynomial::new(circuit_size, circuit_size, 0);

        AllEntities {
            // Precomputed — regular
            q_m: reg(), q_c: reg(), q_l: reg(), q_r: reg(), q_o: reg(), q_4: reg(),
            q_lookup: reg(), q_arith: reg(), q_delta_range: reg(), q_elliptic: reg(),
            q_memory: reg(), q_nnf: reg(), q_poseidon2_external: reg(), q_poseidon2_internal: reg(),
            sigma_1: reg(), sigma_2: reg(), sigma_3: reg(), sigma_4: reg(),
            id_1: reg(), id_2: reg(), id_3: reg(), id_4: reg(),
            table_1: reg(), table_2: reg(), table_3: reg(), table_4: reg(),
            lagrange_first: reg(), lagrange_last: reg(),
            // Witness — to_be_shifted are shiftable
            w_l: shift(), w_r: shift(), w_o: shift(), w_4: shift(),
            z_perm: shift(),
            // Witness — not shifted, regular
            lookup_inverses: reg(), lookup_read_counts: reg(), lookup_read_tags: reg(),
            // Shifted — set via set_shifted()
            w_l_shift: empty(), w_r_shift: empty(), w_o_shift: empty(),
            w_4_shift: empty(), z_perm_shift: empty(),
        }
    }

    pub fn get_polynomial_size(&self) -> usize {
        self.q_m.virtual_size()
    }

    /// Set shifted polynomials from their to-be-shifted counterparts.
    pub fn set_shifted(&mut self) {
        self.w_l_shift = self.w_l.shifted();
        self.w_r_shift = self.w_r.shifted();
        self.w_o_shift = self.w_o.shifted();
        self.w_4_shift = self.w_4.shifted();
        self.z_perm_shift = self.z_perm.shifted();
    }

    /// Get row values at the given index.
    pub fn get_row(&self, row_idx: usize) -> AllValues<P> {
        AllEntities {
            q_m: self.q_m.get(row_idx),
            q_c: self.q_c.get(row_idx),
            q_l: self.q_l.get(row_idx),
            q_r: self.q_r.get(row_idx),
            q_o: self.q_o.get(row_idx),
            q_4: self.q_4.get(row_idx),
            q_lookup: self.q_lookup.get(row_idx),
            q_arith: self.q_arith.get(row_idx),
            q_delta_range: self.q_delta_range.get(row_idx),
            q_elliptic: self.q_elliptic.get(row_idx),
            q_memory: self.q_memory.get(row_idx),
            q_nnf: self.q_nnf.get(row_idx),
            q_poseidon2_external: self.q_poseidon2_external.get(row_idx),
            q_poseidon2_internal: self.q_poseidon2_internal.get(row_idx),
            sigma_1: self.sigma_1.get(row_idx),
            sigma_2: self.sigma_2.get(row_idx),
            sigma_3: self.sigma_3.get(row_idx),
            sigma_4: self.sigma_4.get(row_idx),
            id_1: self.id_1.get(row_idx),
            id_2: self.id_2.get(row_idx),
            id_3: self.id_3.get(row_idx),
            id_4: self.id_4.get(row_idx),
            table_1: self.table_1.get(row_idx),
            table_2: self.table_2.get(row_idx),
            table_3: self.table_3.get(row_idx),
            table_4: self.table_4.get(row_idx),
            lagrange_first: self.lagrange_first.get(row_idx),
            lagrange_last: self.lagrange_last.get(row_idx),
            w_l: self.w_l.get(row_idx),
            w_r: self.w_r.get(row_idx),
            w_o: self.w_o.get(row_idx),
            w_4: self.w_4.get(row_idx),
            z_perm: self.z_perm.get(row_idx),
            lookup_inverses: self.lookup_inverses.get(row_idx),
            lookup_read_counts: self.lookup_read_counts.get(row_idx),
            lookup_read_tags: self.lookup_read_tags.get(row_idx),
            w_l_shift: self.w_l_shift.get(row_idx),
            w_r_shift: self.w_r_shift.get(row_idx),
            w_o_shift: self.w_o_shift.get(row_idx),
            w_4_shift: self.w_4_shift.get(row_idx),
            z_perm_shift: self.z_perm_shift.get(row_idx),
        }
    }

    /// Get row values relevant to the permutation argument (wires + sigmas + ids only).
    ///
    /// Port of C++ `get_row_for_permutation_arg()`. Used as an optimization to avoid
    /// reading all polynomials when only permutation-relevant ones are needed.
    pub fn get_row_for_permutation_arg(&self, row_idx: usize) -> AllValues<P> {
        let mut result = AllValues::<P>::zero();
        // Wires
        result.w_l = self.w_l.get(row_idx);
        result.w_r = self.w_r.get(row_idx);
        result.w_o = self.w_o.get(row_idx);
        result.w_4 = self.w_4.get(row_idx);
        // Sigmas
        result.sigma_1 = self.sigma_1.get(row_idx);
        result.sigma_2 = self.sigma_2.get(row_idx);
        result.sigma_3 = self.sigma_3.get(row_idx);
        result.sigma_4 = self.sigma_4.get(row_idx);
        // IDs
        result.id_1 = self.id_1.get(row_idx);
        result.id_2 = self.id_2.get(row_idx);
        result.id_3 = self.id_3.get(row_idx);
        result.id_4 = self.id_4.get(row_idx);
        result
    }
}

// ============================================================================
// ExtendedEdges — AllEntities<Univariate<P, N>>
// ============================================================================

pub type ExtendedEdges<P, const N: usize> = AllEntities<Univariate<P, N>>;

impl<P: FieldParams, const N: usize> ExtendedEdges<P, N> {
    pub fn zero() -> Self {
        let z = || Univariate::zero();
        AllEntities {
            q_m: z(), q_c: z(), q_l: z(), q_r: z(), q_o: z(), q_4: z(),
            q_lookup: z(), q_arith: z(), q_delta_range: z(), q_elliptic: z(),
            q_memory: z(), q_nnf: z(), q_poseidon2_external: z(), q_poseidon2_internal: z(),
            sigma_1: z(), sigma_2: z(), sigma_3: z(), sigma_4: z(),
            id_1: z(), id_2: z(), id_3: z(), id_4: z(),
            table_1: z(), table_2: z(), table_3: z(), table_4: z(),
            lagrange_first: z(), lagrange_last: z(),
            w_l: z(), w_r: z(), w_o: z(), w_4: z(),
            z_perm: z(), lookup_inverses: z(), lookup_read_counts: z(), lookup_read_tags: z(),
            w_l_shift: z(), w_r_shift: z(), w_o_shift: z(), w_4_shift: z(), z_perm_shift: z(),
        }
    }
}

// ============================================================================
// PartiallyEvaluatedMultivariates — AllEntities<Polynomial<P>>
// ============================================================================

pub type PartiallyEvaluatedMultivariates<P> = AllEntities<Polynomial<P>>;

impl<P: FieldParams> PartiallyEvaluatedMultivariates<P> {
    pub fn from_prover_polynomials(
        full_polynomials: &ProverPolynomials<P>,
        circuit_size: usize,
    ) -> Self {
        let half_size = circuit_size / 2;
        let all = full_polynomials.get_all();
        let mut result = Self::new_uniform(half_size);
        for (dest, src) in result.get_all_mut().iter_mut().zip(all.iter()) {
            let desired_size = src.end_index() / 2 + src.end_index() % 2;
            **dest = Polynomial::new(desired_size, half_size, 0);
        }
        result
    }

    fn new_uniform(size: usize) -> Self {
        let p = || Polynomial::new(size, size, 0);
        AllEntities {
            q_m: p(), q_c: p(), q_l: p(), q_r: p(), q_o: p(), q_4: p(),
            q_lookup: p(), q_arith: p(), q_delta_range: p(), q_elliptic: p(),
            q_memory: p(), q_nnf: p(), q_poseidon2_external: p(), q_poseidon2_internal: p(),
            sigma_1: p(), sigma_2: p(), sigma_3: p(), sigma_4: p(),
            id_1: p(), id_2: p(), id_3: p(), id_4: p(),
            table_1: p(), table_2: p(), table_3: p(), table_4: p(),
            lagrange_first: p(), lagrange_last: p(),
            w_l: p(), w_r: p(), w_o: p(), w_4: p(),
            z_perm: p(), lookup_inverses: p(), lookup_read_counts: p(), lookup_read_tags: p(),
            w_l_shift: p(), w_r_shift: p(), w_o_shift: p(), w_4_shift: p(), z_perm_shift: p(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_ultra_flavor_entity_counts() {
        assert_eq!(NUM_PRECOMPUTED_ENTITIES, 28);
        assert_eq!(NUM_WITNESS_ENTITIES, 8);
        assert_eq!(NUM_SHIFTED_ENTITIES, 5);
        assert_eq!(NUM_ALL_ENTITIES, 41);
    }

    #[test]
    fn test_ultra_flavor_constants() {
        assert_eq!(MAX_PARTIAL_RELATION_LENGTH, 7);
        assert_eq!(BATCHED_RELATION_PARTIAL_LENGTH, 8);
        assert_eq!(NUM_SUBRELATIONS, 28);
        assert_eq!(NUM_RELATIONS, 9);
        assert_eq!(ALL_SUBRELATION_PARTIAL_LENGTHS.len(), NUM_SUBRELATIONS);
        assert_eq!(ALL_SUBRELATION_LINEARLY_INDEPENDENT.len(), NUM_SUBRELATIONS);
    }

    #[test]
    fn test_ultra_all_entities_get_all_size() {
        let values = AllValues::<P>::zero();
        assert_eq!(values.get_all().len(), NUM_ALL_ENTITIES);
        assert_eq!(values.get_precomputed().len(), NUM_PRECOMPUTED_ENTITIES);
        assert_eq!(values.get_witness().len(), NUM_WITNESS_ENTITIES);
        assert_eq!(values.get_shifted().len(), NUM_SHIFTED_ENTITIES);
        assert_eq!(values.get_wires().len(), NUM_WIRES);
        assert_eq!(values.get_sigmas().len(), 4);
        assert_eq!(values.get_ids().len(), 4);
        assert_eq!(values.get_tables().len(), 4);
        assert_eq!(values.get_to_be_shifted().len(), NUM_SHIFTED_ENTITIES);
    }

    #[test]
    fn test_ultra_all_entities_labels() {
        let labels = AllEntities::<Fr>::get_labels();
        assert_eq!(labels.len(), NUM_ALL_ENTITIES);
        assert_eq!(labels[0], "q_m");
        assert_eq!(labels[1], "q_c");
        assert_eq!(labels[27], "lagrange_last");
        assert_eq!(labels[28], "w_l");
        assert_eq!(labels[35], "lookup_read_tags");
        assert_eq!(labels[36], "w_l_shift");
        assert_eq!(labels[40], "z_perm_shift");
    }

    #[test]
    fn test_ultra_prover_polynomials_construction() {
        let circuit_size = 8;
        let polys = ProverPolynomials::<P>::new(circuit_size);

        // Precomputed polys: regular (start_index = 0)
        assert_eq!(polys.q_m.start_index(), 0);
        assert_eq!(polys.q_m.virtual_size(), circuit_size);

        // Shiftable witness polys: start_index = 1
        assert_eq!(polys.w_l.start_index(), 1);
        assert_eq!(polys.w_r.start_index(), 1);
        assert_eq!(polys.z_perm.start_index(), 1);

        // Non-shiftable witness polys: regular (start_index = 0)
        assert_eq!(polys.lookup_inverses.start_index(), 0);
        assert_eq!(polys.lookup_read_counts.start_index(), 0);
    }

    #[test]
    fn test_ultra_prover_polynomials_set_shifted() {
        let circuit_size = 8;
        let mut polys = ProverPolynomials::<P>::new(circuit_size);

        *polys.w_l.at_mut(1) = Fr::from(10u64);
        *polys.w_l.at_mut(2) = Fr::from(20u64);
        *polys.z_perm.at_mut(1) = Fr::from(100u64);
        *polys.z_perm.at_mut(2) = Fr::from(200u64);

        polys.set_shifted();

        assert_eq!(polys.w_l_shift.get(0), Fr::from(10u64));
        assert_eq!(polys.w_l_shift.get(1), Fr::from(20u64));
        assert_eq!(polys.z_perm_shift.get(0), Fr::from(100u64));
        assert_eq!(polys.z_perm_shift.get(1), Fr::from(200u64));
    }

    #[test]
    fn test_ultra_prover_polynomials_get_row() {
        let circuit_size = 4;
        let mut polys = ProverPolynomials::<P>::new(circuit_size);

        *polys.q_arith.at_mut(0) = Fr::from(1u64);
        *polys.q_arith.at_mut(1) = Fr::from(2u64);
        *polys.w_l.at_mut(1) = Fr::from(10u64);
        *polys.sigma_1.at_mut(1) = Fr::from(50u64);

        polys.set_shifted();

        let row1 = polys.get_row(1);
        assert_eq!(row1.q_arith, Fr::from(2u64));
        assert_eq!(row1.w_l, Fr::from(10u64));
        assert_eq!(row1.sigma_1, Fr::from(50u64));
    }

    #[test]
    fn test_ultra_extended_edges_zero() {
        let edges = ExtendedEdges::<P, 8>::zero();
        assert!(edges.q_m.is_zero());
        assert!(edges.w_l.is_zero());
        assert!(edges.z_perm_shift.is_zero());
        assert_eq!(edges.get_all().len(), NUM_ALL_ENTITIES);
    }
}
