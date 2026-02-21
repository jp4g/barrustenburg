//! Port of `sumcheck_test_flavor.hpp` — Minimal test flavors for sumcheck testing.
//!
//! Provides `SumcheckTestFlavor` and variants for testing sumcheck in isolation
//! without UltraFlavor dependencies.
//!
//! Available flavors:
//! - `SumcheckTestFlavor`: Base (BN254, non-ZK, short monomials, arithmetic + dependent test)
//! - `SumcheckTestFlavorZK`: Zero-knowledge variant (HasZK = true)
//! - `SumcheckTestFlavorFullBary`: Full barycentric extension
//!
//! Entity counts:
//! - Precomputed: q_m, q_l, q_r, q_o, q_4, q_c, q_arith, q_test = 8
//! - Witness: w_l, w_r, w_o, w_4, w_test_1, w_test_2 = 6
//! - Shifted: w_l_shift, w_4_shift = 2
//! - Total: 16

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_polynomials::polynomial::Polynomial;
use bbrs_polynomials::univariate::Univariate;

// ============================================================================
// DependentTestRelation
// ============================================================================

/// Subrelation partial lengths for the DependentTestRelation.
pub const DEPENDENT_TEST_SUBRELATION_PARTIAL_LENGTHS: [usize; 1] = [2];

/// Subrelation linear independence flags for DependentTestRelation.
/// This subrelation is NOT linearly independent (should NOT be scaled by sumcheck).
pub const DEPENDENT_TEST_SUBRELATION_LINEARLY_INDEPENDENT: [bool; 1] = [false];

/// Accumulate the DependentTestRelation.
///
/// Relation: q_test * w_test_1
///
/// This is a simple linearly *dependent* relation used to test that sumcheck
/// correctly handles the SUBRELATION_LINEARLY_INDEPENDENT flag.
pub fn dependent_test_accumulate<P: FieldParams>(
    evals: &mut [Field<P>; 1],
    input: &AllValues<P>,
    _params: &bbrs_relations::relation_parameters::RelationParameters<Field<P>>,
    _scaling_factor: &Field<P>,
) {
    // Note: NO scaling_factor used here — this is linearly dependent!
    evals[0] = evals[0] + input.q_test * input.w_test_1;
}

/// Skip check for DependentTestRelation.
pub fn dependent_test_skip<P: FieldParams>(input: &AllValues<P>) -> bool {
    input.q_test.is_zero()
}

// ============================================================================
// Entity structs
// ============================================================================

/// Precomputed entities (selectors) for SumcheckTestFlavor.
///
/// 8 entities: q_m, q_l, q_r, q_o, q_4, q_c, q_arith, q_test
pub struct PrecomputedEntities<T> {
    pub q_m: T,
    pub q_l: T,
    pub q_r: T,
    pub q_o: T,
    pub q_4: T,
    pub q_c: T,
    pub q_arith: T,
    pub q_test: T,
}

impl<T> PrecomputedEntities<T> {
    pub const NUM: usize = 8;
}

/// Witness entities for SumcheckTestFlavor.
///
/// 6 entities: w_l, w_r, w_o, w_4, w_test_1, w_test_2
pub struct WitnessEntities<T> {
    pub w_l: T,
    pub w_r: T,
    pub w_o: T,
    pub w_4: T,
    pub w_test_1: T,
    pub w_test_2: T,
}

impl<T> WitnessEntities<T> {
    pub const NUM: usize = 6;
}

/// Shifted entities for SumcheckTestFlavor.
///
/// 2 entities: w_l_shift, w_4_shift
pub struct ShiftedEntities<T> {
    pub w_l_shift: T,
    pub w_4_shift: T,
}

impl<T> ShiftedEntities<T> {
    pub const NUM: usize = 2;
}

// ============================================================================
// AllEntities — combined flat view of all entity groups
// ============================================================================

/// Total number of entities in SumcheckTestFlavor.
pub const NUM_PRECOMPUTED_ENTITIES: usize = 8;
pub const NUM_WITNESS_ENTITIES: usize = 6;
pub const NUM_SHIFTED_ENTITIES: usize = 2;
pub const NUM_ALL_ENTITIES: usize =
    NUM_PRECOMPUTED_ENTITIES + NUM_WITNESS_ENTITIES + NUM_SHIFTED_ENTITIES;

/// All entities combined for SumcheckTestFlavor (16 total).
///
/// Layout matches C++ ordering: precomputed, witness, shifted.
pub struct AllEntities<T> {
    // Precomputed (8)
    pub q_m: T,
    pub q_l: T,
    pub q_r: T,
    pub q_o: T,
    pub q_4: T,
    pub q_c: T,
    pub q_arith: T,
    pub q_test: T,
    // Witness (6)
    pub w_l: T,
    pub w_r: T,
    pub w_o: T,
    pub w_4: T,
    pub w_test_1: T,
    pub w_test_2: T,
    // Shifted (2)
    pub w_l_shift: T,
    pub w_4_shift: T,
}

impl<T> AllEntities<T> {
    /// Get references to all fields in layout order.
    pub fn get_all(&self) -> [&T; NUM_ALL_ENTITIES] {
        [
            &self.q_m,
            &self.q_l,
            &self.q_r,
            &self.q_o,
            &self.q_4,
            &self.q_c,
            &self.q_arith,
            &self.q_test,
            &self.w_l,
            &self.w_r,
            &self.w_o,
            &self.w_4,
            &self.w_test_1,
            &self.w_test_2,
            &self.w_l_shift,
            &self.w_4_shift,
        ]
    }

    /// Get mutable references to all fields in layout order.
    pub fn get_all_mut(&mut self) -> [&mut T; NUM_ALL_ENTITIES] {
        [
            &mut self.q_m,
            &mut self.q_l,
            &mut self.q_r,
            &mut self.q_o,
            &mut self.q_4,
            &mut self.q_c,
            &mut self.q_arith,
            &mut self.q_test,
            &mut self.w_l,
            &mut self.w_r,
            &mut self.w_o,
            &mut self.w_4,
            &mut self.w_test_1,
            &mut self.w_test_2,
            &mut self.w_l_shift,
            &mut self.w_4_shift,
        ]
    }

    /// Get references to precomputed entities.
    pub fn get_precomputed(&self) -> [&T; NUM_PRECOMPUTED_ENTITIES] {
        [
            &self.q_m,
            &self.q_l,
            &self.q_r,
            &self.q_o,
            &self.q_4,
            &self.q_c,
            &self.q_arith,
            &self.q_test,
        ]
    }

    /// Get mutable references to precomputed entities.
    pub fn get_precomputed_mut(&mut self) -> [&mut T; NUM_PRECOMPUTED_ENTITIES] {
        [
            &mut self.q_m,
            &mut self.q_l,
            &mut self.q_r,
            &mut self.q_o,
            &mut self.q_4,
            &mut self.q_c,
            &mut self.q_arith,
            &mut self.q_test,
        ]
    }

    /// Get references to witness entities.
    pub fn get_witness(&self) -> [&T; NUM_WITNESS_ENTITIES] {
        [
            &self.w_l,
            &self.w_r,
            &self.w_o,
            &self.w_4,
            &self.w_test_1,
            &self.w_test_2,
        ]
    }

    /// Get mutable references to witness entities.
    pub fn get_witness_mut(&mut self) -> [&mut T; NUM_WITNESS_ENTITIES] {
        [
            &mut self.w_l,
            &mut self.w_r,
            &mut self.w_o,
            &mut self.w_4,
            &mut self.w_test_1,
            &mut self.w_test_2,
        ]
    }

    /// Get references to shifted entities.
    pub fn get_shifted(&self) -> [&T; NUM_SHIFTED_ENTITIES] {
        [&self.w_l_shift, &self.w_4_shift]
    }

    /// Get mutable references to shifted entities.
    pub fn get_shifted_mut(&mut self) -> [&mut T; NUM_SHIFTED_ENTITIES] {
        [&mut self.w_l_shift, &mut self.w_4_shift]
    }

    /// Get references to unshifted entities (precomputed + witness).
    pub fn get_unshifted(&self) -> [&T; NUM_PRECOMPUTED_ENTITIES + NUM_WITNESS_ENTITIES] {
        [
            &self.q_m,
            &self.q_l,
            &self.q_r,
            &self.q_o,
            &self.q_4,
            &self.q_c,
            &self.q_arith,
            &self.q_test,
            &self.w_l,
            &self.w_r,
            &self.w_o,
            &self.w_4,
            &self.w_test_1,
            &self.w_test_2,
        ]
    }

    /// Get field labels.
    pub fn get_labels() -> &'static [&'static str; NUM_ALL_ENTITIES] {
        &[
            "q_m",
            "q_l",
            "q_r",
            "q_o",
            "q_4",
            "q_c",
            "q_arith",
            "q_test",
            "w_l",
            "w_r",
            "w_o",
            "w_4",
            "w_test_1",
            "w_test_2",
            "w_l_shift",
            "w_4_shift",
        ]
    }
}

// ============================================================================
// AllValues — AllEntities<FF> for verifier-side evaluations
// ============================================================================

/// A field element for each entity. Used for sumcheck evaluations at a point.
pub type AllValues<P> = AllEntities<Field<P>>;

impl<P: FieldParams> AllValues<P> {
    /// Create all-zero values.
    pub fn zero() -> Self {
        AllEntities {
            q_m: Field::zero(),
            q_l: Field::zero(),
            q_r: Field::zero(),
            q_o: Field::zero(),
            q_4: Field::zero(),
            q_c: Field::zero(),
            q_arith: Field::zero(),
            q_test: Field::zero(),
            w_l: Field::zero(),
            w_r: Field::zero(),
            w_o: Field::zero(),
            w_4: Field::zero(),
            w_test_1: Field::zero(),
            w_test_2: Field::zero(),
            w_l_shift: Field::zero(),
            w_4_shift: Field::zero(),
        }
    }
}

// ============================================================================
// ProverPolynomials — AllEntities<Polynomial>
// ============================================================================

/// Container for prover polynomials in SumcheckTestFlavor.
pub type ProverPolynomials<P> = AllEntities<Polynomial<P>>;

impl<P: FieldParams> ProverPolynomials<P> {
    /// Create prover polynomials of the given circuit size, all initialized to zero.
    ///
    /// Precomputed (selectors) are regular polynomials starting at index 0.
    /// Witness polynomials are shiftable (start_index = 1) so they can be shifted.
    /// Shifted polynomials are empty until `set_shifted()` is called.
    pub fn new(circuit_size: usize) -> Self {
        AllEntities {
            // Precomputed — regular polynomials
            q_m: Polynomial::new(circuit_size, circuit_size, 0),
            q_l: Polynomial::new(circuit_size, circuit_size, 0),
            q_r: Polynomial::new(circuit_size, circuit_size, 0),
            q_o: Polynomial::new(circuit_size, circuit_size, 0),
            q_4: Polynomial::new(circuit_size, circuit_size, 0),
            q_c: Polynomial::new(circuit_size, circuit_size, 0),
            q_arith: Polynomial::new(circuit_size, circuit_size, 0),
            q_test: Polynomial::new(circuit_size, circuit_size, 0),
            // Witness — shiftable (start_index = 1)
            w_l: Polynomial::shiftable(circuit_size),
            w_r: Polynomial::shiftable(circuit_size),
            w_o: Polynomial::shiftable(circuit_size),
            w_4: Polynomial::shiftable(circuit_size),
            w_test_1: Polynomial::shiftable(circuit_size),
            w_test_2: Polynomial::shiftable(circuit_size),
            // Shifted — set via `set_shifted()`
            w_l_shift: Polynomial::new(circuit_size, circuit_size, 0),
            w_4_shift: Polynomial::new(circuit_size, circuit_size, 0),
        }
    }

    /// Get the circuit (virtual) size from a precomputed polynomial.
    pub fn get_polynomial_size(&self) -> usize {
        self.q_m.virtual_size()
    }

    /// Get references to polynomials that will be shifted (w_l, w_4).
    pub fn get_to_be_shifted(&self) -> [&Polynomial<P>; NUM_SHIFTED_ENTITIES] {
        [&self.w_l, &self.w_4]
    }

    /// Set shifted polynomials from their to-be-shifted counterparts.
    pub fn set_shifted(&mut self) {
        self.w_l_shift = self.w_l.shifted();
        self.w_4_shift = self.w_4.shifted();
    }

    /// Get row values at the given index.
    pub fn get_row(&self, row_idx: usize) -> AllValues<P> {
        AllEntities {
            q_m: self.q_m.get(row_idx),
            q_l: self.q_l.get(row_idx),
            q_r: self.q_r.get(row_idx),
            q_o: self.q_o.get(row_idx),
            q_4: self.q_4.get(row_idx),
            q_c: self.q_c.get(row_idx),
            q_arith: self.q_arith.get(row_idx),
            q_test: self.q_test.get(row_idx),
            w_l: self.w_l.get(row_idx),
            w_r: self.w_r.get(row_idx),
            w_o: self.w_o.get(row_idx),
            w_4: self.w_4.get(row_idx),
            w_test_1: self.w_test_1.get(row_idx),
            w_test_2: self.w_test_2.get(row_idx),
            w_l_shift: self.w_l_shift.get(row_idx),
            w_4_shift: self.w_4_shift.get(row_idx),
        }
    }
}

// ============================================================================
// ExtendedEdges — AllEntities<Univariate<FF, MAX_PARTIAL_RELATION_LENGTH>>
// ============================================================================

/// Extended edges for sumcheck, containing Univariates at each entity position.
pub type ExtendedEdges<P, const N: usize> = AllEntities<Univariate<P, N>>;

impl<P: FieldParams, const N: usize> ExtendedEdges<P, N> {
    /// Create extended edges initialized to zero univariates.
    pub fn zero() -> Self {
        AllEntities {
            q_m: Univariate::zero(),
            q_l: Univariate::zero(),
            q_r: Univariate::zero(),
            q_o: Univariate::zero(),
            q_4: Univariate::zero(),
            q_c: Univariate::zero(),
            q_arith: Univariate::zero(),
            q_test: Univariate::zero(),
            w_l: Univariate::zero(),
            w_r: Univariate::zero(),
            w_o: Univariate::zero(),
            w_4: Univariate::zero(),
            w_test_1: Univariate::zero(),
            w_test_2: Univariate::zero(),
            w_l_shift: Univariate::zero(),
            w_4_shift: Univariate::zero(),
        }
    }
}

// ============================================================================
// PartiallyEvaluatedMultivariates — AllEntities<Polynomial> with half-size
// ============================================================================

/// Partially evaluated multivariates for folded sumcheck rounds.
pub type PartiallyEvaluatedMultivariates<P> = AllEntities<Polynomial<P>>;

impl<P: FieldParams> PartiallyEvaluatedMultivariates<P> {
    /// Create from full prover polynomials, allocating half-size storage.
    pub fn from_prover_polynomials(
        full_polynomials: &ProverPolynomials<P>,
        circuit_size: usize,
    ) -> Self {
        let half_size = circuit_size / 2;
        let all = full_polynomials.get_all();
        let mut result = Self::new_half_size(half_size);
        // Allocate each polynomial with the right size based on the full polynomial's end_index
        for (dest, src) in result.get_all_mut().iter_mut().zip(all.iter()) {
            let desired_size = src.end_index() / 2 + src.end_index() % 2;
            **dest = Polynomial::new(desired_size, half_size, 0);
        }
        result
    }

    /// Create with uniform half-size polynomials.
    fn new_half_size(half_size: usize) -> Self {
        AllEntities {
            q_m: Polynomial::new(half_size, half_size, 0),
            q_l: Polynomial::new(half_size, half_size, 0),
            q_r: Polynomial::new(half_size, half_size, 0),
            q_o: Polynomial::new(half_size, half_size, 0),
            q_4: Polynomial::new(half_size, half_size, 0),
            q_c: Polynomial::new(half_size, half_size, 0),
            q_arith: Polynomial::new(half_size, half_size, 0),
            q_test: Polynomial::new(half_size, half_size, 0),
            w_l: Polynomial::new(half_size, half_size, 0),
            w_r: Polynomial::new(half_size, half_size, 0),
            w_o: Polynomial::new(half_size, half_size, 0),
            w_4: Polynomial::new(half_size, half_size, 0),
            w_test_1: Polynomial::new(half_size, half_size, 0),
            w_test_2: Polynomial::new(half_size, half_size, 0),
            w_l_shift: Polynomial::new(half_size, half_size, 0),
            w_4_shift: Polynomial::new(half_size, half_size, 0),
        }
    }
}

// ============================================================================
// SumcheckTestFlavor constants
// ============================================================================

/// Number of wires.
pub const NUM_WIRES: usize = 4;

/// Number of relations: ArithmeticRelation + DependentTestRelation = 2.
pub const NUM_RELATIONS: usize = 2;

/// Subrelation partial lengths across all relations.
///
/// ArithmeticRelation: [6, 5] = 2 subrelations
/// DependentTestRelation: [2] = 1 subrelation
/// Total: 3 subrelations
pub const NUM_SUBRELATIONS: usize = 3;

/// Maximum partial relation length (from ArithmeticRelation's [6, 5]).
pub const MAX_PARTIAL_RELATION_LENGTH: usize = 6;

/// Batched relation partial length = MAX_PARTIAL_RELATION_LENGTH + 1 (non-ZK).
pub const BATCHED_RELATION_PARTIAL_LENGTH: usize = MAX_PARTIAL_RELATION_LENGTH + 1; // = 7

/// Batched relation partial length for ZK flavors = MAX_PARTIAL_RELATION_LENGTH + 3.
pub const BATCHED_RELATION_PARTIAL_LENGTH_ZK: usize = MAX_PARTIAL_RELATION_LENGTH + 3; // = 9

/// Whether this flavor has a zero row (no, for test flavor).
pub const HAS_ZERO_ROW: bool = false;

/// All subrelation partial lengths in flat order:
/// Arithmetic[0]=6, Arithmetic[1]=5, DependentTest[0]=2
pub const ALL_SUBRELATION_PARTIAL_LENGTHS: [usize; NUM_SUBRELATIONS] = [6, 5, 2];

/// Subrelation linear independence flags in flat order:
/// Arithmetic: [true, true], DependentTest: [false]
pub const ALL_SUBRELATION_LINEARLY_INDEPENDENT: [bool; NUM_SUBRELATIONS] = [true, true, false];

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fr};

    type P = Bn254FrParams;

    #[test]
    fn test_flavor_entity_counts() {
        assert_eq!(NUM_PRECOMPUTED_ENTITIES, 8);
        assert_eq!(NUM_WITNESS_ENTITIES, 6);
        assert_eq!(NUM_SHIFTED_ENTITIES, 2);
        assert_eq!(NUM_ALL_ENTITIES, 16);
    }

    #[test]
    fn test_all_entities_get_all_size() {
        let values = AllValues::<P>::zero();
        assert_eq!(values.get_all().len(), NUM_ALL_ENTITIES);
        assert_eq!(values.get_precomputed().len(), NUM_PRECOMPUTED_ENTITIES);
        assert_eq!(values.get_witness().len(), NUM_WITNESS_ENTITIES);
        assert_eq!(values.get_shifted().len(), NUM_SHIFTED_ENTITIES);
        assert_eq!(
            values.get_unshifted().len(),
            NUM_PRECOMPUTED_ENTITIES + NUM_WITNESS_ENTITIES
        );
    }

    #[test]
    fn test_all_entities_labels() {
        let labels = AllEntities::<Fr>::get_labels();
        assert_eq!(labels.len(), NUM_ALL_ENTITIES);
        assert_eq!(labels[0], "q_m");
        assert_eq!(labels[7], "q_test");
        assert_eq!(labels[8], "w_l");
        assert_eq!(labels[14], "w_l_shift");
        assert_eq!(labels[15], "w_4_shift");
    }

    #[test]
    fn test_prover_polynomials_construction() {
        let circuit_size = 4;
        let polys = ProverPolynomials::<P>::new(circuit_size);

        // w_l is shiftable: starts at index 1, with circuit_size - 1 real elements
        assert_eq!(polys.w_l.start_index(), 1);
        assert_eq!(polys.w_l.size(), circuit_size - 1);

        // Precomputed polys have circuit_size entries starting at 0
        for poly in polys.get_precomputed() {
            assert_eq!(poly.size(), circuit_size);
            assert_eq!(poly.start_index(), 0);
        }
    }

    #[test]
    fn test_prover_polynomials_set_shifted() {
        let circuit_size = 4;
        let mut polys = ProverPolynomials::<P>::new(circuit_size);

        // Set values in w_l (shiftable, indices [1, 2, 3])
        *polys.w_l.at_mut(1) = Fr::from(10u64);
        *polys.w_l.at_mut(2) = Fr::from(20u64);
        *polys.w_l.at_mut(3) = Fr::from(30u64);

        // Set values in w_4 (shiftable, indices [1, 2, 3])
        *polys.w_4.at_mut(1) = Fr::from(100u64);
        *polys.w_4.at_mut(2) = Fr::from(200u64);
        *polys.w_4.at_mut(3) = Fr::from(300u64);

        polys.set_shifted();

        // w_l_shift[i] = w_l[i+1] (shifted() decrements start_index: 1 -> 0)
        assert_eq!(polys.w_l_shift.get(0), Fr::from(10u64)); // w_l[1]
        assert_eq!(polys.w_l_shift.get(1), Fr::from(20u64)); // w_l[2]
        assert_eq!(polys.w_l_shift.get(2), Fr::from(30u64)); // w_l[3]

        // w_4_shift similarly
        assert_eq!(polys.w_4_shift.get(0), Fr::from(100u64));
        assert_eq!(polys.w_4_shift.get(1), Fr::from(200u64));
        assert_eq!(polys.w_4_shift.get(2), Fr::from(300u64));
    }

    #[test]
    fn test_prover_polynomials_get_row() {
        let circuit_size = 4;
        let mut polys = ProverPolynomials::<P>::new(circuit_size);

        // Set precomputed values (start at index 0)
        *polys.q_arith.at_mut(0) = Fr::from(1u64);
        *polys.q_arith.at_mut(1) = Fr::from(2u64);

        // Set witness values (shiftable, start at index 1)
        *polys.w_l.at_mut(1) = Fr::from(10u64);
        *polys.w_l.at_mut(2) = Fr::from(20u64);
        *polys.w_r.at_mut(1) = Fr::from(100u64);
        *polys.w_r.at_mut(2) = Fr::from(200u64);

        polys.set_shifted();

        let row0 = polys.get_row(0);
        assert_eq!(row0.q_arith, Fr::from(1u64));
        assert_eq!(row0.w_l, Fr::zero()); // shiftable: index 0 is virtual zero

        let row1 = polys.get_row(1);
        assert_eq!(row1.q_arith, Fr::from(2u64));
        assert_eq!(row1.w_l, Fr::from(10u64));
        assert_eq!(row1.w_r, Fr::from(100u64));
        // w_l_shift at row 1 = w_l[2] = 20
        assert_eq!(row1.w_l_shift, Fr::from(20u64));
    }

    #[test]
    fn test_dependent_test_relation() {
        let mut values = AllValues::<P>::zero();
        values.q_test = Fr::from(3u64);
        values.w_test_1 = Fr::from(7u64);

        let params = bbrs_relations::relation_parameters::RelationParameters::<Fr>::get_random();
        let one = Fr::one();
        let mut evals = [Fr::zero(); 1];

        dependent_test_accumulate(&mut evals, &values, &params, &one);

        // Result should be q_test * w_test_1 = 3 * 7 = 21
        assert_eq!(evals[0], Fr::from(21u64));
    }

    #[test]
    fn test_dependent_test_skip() {
        let mut values = AllValues::<P>::zero();
        assert!(dependent_test_skip(&values));

        values.q_test = Fr::from(1u64);
        assert!(!dependent_test_skip(&values));
    }

    #[test]
    fn test_flavor_constants() {
        assert_eq!(MAX_PARTIAL_RELATION_LENGTH, 6);
        assert_eq!(BATCHED_RELATION_PARTIAL_LENGTH, 7);
        assert_eq!(BATCHED_RELATION_PARTIAL_LENGTH_ZK, 9);
        assert_eq!(NUM_SUBRELATIONS, 3);
        assert_eq!(NUM_RELATIONS, 2);
        assert_eq!(ALL_SUBRELATION_PARTIAL_LENGTHS, [6, 5, 2]);
        assert_eq!(
            ALL_SUBRELATION_LINEARLY_INDEPENDENT,
            [true, true, false]
        );
    }

    #[test]
    fn test_extended_edges_zero() {
        let edges = ExtendedEdges::<P, 6>::zero();
        assert!(edges.q_m.is_zero());
        assert!(edges.w_l.is_zero());
        assert!(edges.w_l_shift.is_zero());
        assert_eq!(edges.get_all().len(), NUM_ALL_ENTITIES);
    }
}
