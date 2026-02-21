//! Port of `logderivative_library.hpp` — compute log-derivative inverse polynomials.
//!
//! Provides traits and functions for log-derivative lookup and permutation arguments.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_polynomials::polynomial::Polynomial;

use bbrs_relations::relation_parameters::RelationParameters;

/// Trait for log-derivative lookup/permutation relations.
///
/// Port of the C++ pattern where relations like `LogDerivLookupRelation` define
/// methods for computing read/write terms and managing the inverse polynomial.
pub trait LogDerivRelation<P: FieldParams> {
    /// The type providing row values from the prover polynomials.
    type AllValues;

    /// The type containing all prover polynomials.
    type ProverPolynomials;

    /// Number of read terms in the relation.
    const READ_TERMS: usize;

    /// Number of write terms in the relation.
    const WRITE_TERMS: usize;

    /// Check whether a lookup/permutation operation exists at this row.
    fn operation_exists_at_row(row: &Self::AllValues) -> bool;

    /// Compute a read term at a given row for read_index.
    fn compute_read_term(
        row: &Self::AllValues,
        params: &RelationParameters<Field<P>>,
        read_index: usize,
    ) -> Field<P>;

    /// Compute a write term at a given row for write_index.
    fn compute_write_term(
        row: &Self::AllValues,
        params: &RelationParameters<Field<P>>,
        write_index: usize,
    ) -> Field<P>;

    /// Get a mutable reference to the inverse polynomial.
    fn get_inverse_polynomial_mut(polys: &mut Self::ProverPolynomials) -> &mut Polynomial<P>;

    /// Get the inverse polynomial's start index (offset).
    fn get_inverse_start_index(polys: &Self::ProverPolynomials) -> usize;

    /// Get the inverse polynomial size.
    fn get_inverse_size(polys: &Self::ProverPolynomials) -> usize;

    /// Get row values at a given index.
    fn get_row(polys: &Self::ProverPolynomials, idx: usize) -> Self::AllValues;
}

/// Compute the inverse polynomial I(X) required for log-derivative lookups.
///
/// Port of C++ `compute_logderivative_inverse<FF, Relation, Polynomials>()`.
///
/// For each row i where an operation exists:
///   I[i] = 1 / (∏ read_term[j] * ∏ write_term[k])
///
/// If no operation exists at row i, I[i] = 0.
pub fn compute_logderivative_inverse<P, R>(
    polynomials: &mut R::ProverPolynomials,
    relation_parameters: &RelationParameters<Field<P>>,
    circuit_size: usize,
) where
    P: FieldParams,
    R: LogDerivRelation<P>,
{
    let offset = R::get_inverse_start_index(polynomials);

    // First pass: compute the product of all read/write terms at each row
    // We need to collect row data first since we need immutable access to polys for get_row
    let mut denominator_products = vec![Field::<P>::zero(); circuit_size];
    let mut has_inverse = vec![false; circuit_size];

    for i in 0..circuit_size {
        let row = R::get_row(polynomials, i + offset);
        if !R::operation_exists_at_row(&row) {
            continue;
        }
        has_inverse[i] = true;

        let mut denominator = Field::<P>::one();
        for read_idx in 0..R::READ_TERMS {
            denominator = denominator * R::compute_read_term(&row, relation_parameters, read_idx);
        }
        for write_idx in 0..R::WRITE_TERMS {
            denominator =
                denominator * R::compute_write_term(&row, relation_parameters, write_idx);
        }
        denominator_products[i] = denominator;
    }

    // Batch invert all non-zero denominators
    batch_invert_nonzero(&mut denominator_products);

    // Write results to the inverse polynomial
    let inv_poly = R::get_inverse_polynomial_mut(polynomials);
    for i in 0..circuit_size {
        if has_inverse[i] {
            *inv_poly.at_mut(i) = denominator_products[i];
        }
    }
}

/// Batch-invert a slice of field elements in place, skipping zeros.
fn batch_invert_nonzero<P: FieldParams>(values: &mut [Field<P>]) {
    if values.is_empty() {
        return;
    }

    // Collect non-zero indices
    let nonzero_indices: Vec<usize> = values
        .iter()
        .enumerate()
        .filter(|(_, v)| !v.is_zero())
        .map(|(i, _)| i)
        .collect();

    if nonzero_indices.is_empty() {
        return;
    }

    // Forward pass: compute running product of non-zero values
    let n = nonzero_indices.len();
    let mut scratch = vec![Field::<P>::one(); n];
    scratch[0] = values[nonzero_indices[0]];
    for i in 1..n {
        scratch[i] = scratch[i - 1] * values[nonzero_indices[i]];
    }

    // Invert the final product
    let mut inv_acc = scratch[n - 1].invert();

    // Backward pass: propagate inverses
    for i in (1..n).rev() {
        let idx = nonzero_indices[i];
        let temp = values[idx];
        values[idx] = scratch[i - 1] * inv_acc;
        inv_acc = inv_acc * temp;
    }
    values[nonzero_indices[0]] = inv_acc;
}
