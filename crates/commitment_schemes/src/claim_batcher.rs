//! Claim batching logic for Shplemini verification.
//!
//! Port of C++ `commitment_schemes/claim_batcher.hpp`.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;

/// A batch of commitments and evaluations with a batching scalar.
pub struct Batch<C: CurveParams> {
    pub commitments: Vec<AffineElement<C>>,
    pub evaluations: Vec<Field<C::ScalarFieldParams>>,
    /// Scalar used for batching the claims, excluding the power of batching challenge rho.
    pub scalar: Field<C::ScalarFieldParams>,
}

impl<C: CurveParams> Default for Batch<C> {
    fn default() -> Self {
        Self {
            commitments: Vec::new(),
            evaluations: Vec::new(),
            scalar: Field::zero(),
        }
    }
}

/// An interleaved batch of commitments grouped for interleaving.
pub struct InterleavedBatch<C: CurveParams> {
    pub commitments_groups: Vec<Vec<AffineElement<C>>>,
    pub evaluations: Vec<Field<C::ScalarFieldParams>>,
    pub scalars_pos: Vec<Field<C::ScalarFieldParams>>,
    pub scalars_neg: Vec<Field<C::ScalarFieldParams>>,
    pub shplonk_denominator: Field<C::ScalarFieldParams>,
}

impl<C: CurveParams> Default for InterleavedBatch<C> {
    fn default() -> Self {
        Self {
            commitments_groups: Vec::new(),
            evaluations: Vec::new(),
            scalars_pos: Vec::new(),
            scalars_neg: Vec::new(),
            shplonk_denominator: Field::zero(),
        }
    }
}

/// Logic to support batching opening claims for unshifted and shifted
/// polynomials in Shplemini.
///
/// Stores references to the commitments/evaluations of unshifted and shifted
/// polynomials to be batched opened via Shplemini. Aggregates the commitments
/// and batching scalars for each batch into the corresponding containers for
/// Shplemini. Computes the batched evaluation. Contains logic for computing
/// the per-batch scalars used to batch each set of claims.
///
/// Port of C++ `ClaimBatcher_<Curve>`.
pub struct ClaimBatcher<C: CurveParams> {
    /// Commitments and evaluations of unshifted polynomials.
    pub unshifted: Option<Batch<C>>,
    /// Commitments of to-be-shifted-by-1 polys, evals of their shifts.
    pub shifted: Option<Batch<C>>,
    /// Commitments to groups of polynomials to be combined by interleaving
    /// and evaluations of the resulting interleaved polynomials.
    pub interleaved: Option<InterleavedBatch<C>>,
}

impl<C: CurveParams> Default for ClaimBatcher<C> {
    fn default() -> Self {
        Self {
            unshifted: None,
            shifted: None,
            interleaved: None,
        }
    }
}

impl<C: CurveParams> ClaimBatcher<C> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn get_unshifted(&self) -> Option<&Batch<C>> {
        self.unshifted.as_ref()
    }

    pub fn get_shifted(&self) -> Option<&Batch<C>> {
        self.shifted.as_ref()
    }

    pub fn get_interleaved(&self) -> Option<&InterleavedBatch<C>> {
        self.interleaved.as_ref()
    }

    pub fn get_groups_to_be_interleaved_size(&self) -> u32 {
        match &self.interleaved {
            Some(il) => {
                if il.commitments_groups.is_empty() {
                    0
                } else {
                    il.commitments_groups[0].len() as u32
                }
            }
            None => 0,
        }
    }

    pub fn get_unshifted_batch_scalar(&self) -> Field<C::ScalarFieldParams> {
        match &self.unshifted {
            Some(batch) => batch.scalar,
            None => Field::zero(),
        }
    }

    /// Compute scalars used to batch each set of claims, excluding the
    /// contribution from batching challenge rho.
    ///
    /// Computes scalars s_0, s_1 given by:
    ///   s_0 = (1/(z-r) + nu * 1/(z+r))
    ///   s_1 = r^{-1} * (1/(z-r) - nu * 1/(z+r))
    ///
    /// # Arguments
    /// * `inverted_vanishing_evals` - Precomputed inverted vanishing evaluations
    /// * `nu_challenge` - Shplonk batching challenge nu
    /// * `r_challenge` - Gemini evaluation challenge r
    pub fn compute_scalars_for_each_batch(
        &mut self,
        inverted_vanishing_evals: &[Field<C::ScalarFieldParams>],
        nu_challenge: &Field<C::ScalarFieldParams>,
        r_challenge: &Field<C::ScalarFieldParams>,
    ) {
        let inverse_vanishing_eval_pos = inverted_vanishing_evals[0];
        let inverse_vanishing_eval_neg = inverted_vanishing_evals[1];

        if let Some(ref mut batch) = self.unshifted {
            // (1/(z-r) + nu/(z+r))
            batch.scalar = inverse_vanishing_eval_pos + *nu_challenge * inverse_vanishing_eval_neg;
        }

        if let Some(ref mut batch) = self.shifted {
            // r^{-1} * (1/(z-r) - nu/(z+r))
            batch.scalar = r_challenge.invert()
                * (inverse_vanishing_eval_pos - *nu_challenge * inverse_vanishing_eval_neg);
        }

        let group_size = self.get_groups_to_be_interleaved_size();
        if let Some(ref mut il) = self.interleaved {
            // The index uses 2 * get_msb(group_size) to find the right denominator
            let interleaving_denominator_index = if group_size > 0 {
                2 * (32 - (group_size.leading_zeros() as usize) - 1)
            } else {
                0
            };

            assert!(
                group_size % 2 == 0,
                "Interleaved groups size must be even"
            );

            let mut r_shift_pos = Field::one();
            let mut r_shift_neg = Field::one();
            il.shplonk_denominator = inverted_vanishing_evals[interleaving_denominator_index];

            for i in 0..group_size as usize {
                il.scalars_pos.push(r_shift_pos);
                il.scalars_neg.push(r_shift_neg);
                if i < (group_size as usize) - 1 {
                    r_shift_pos *= *r_challenge;
                    r_shift_neg *= -*r_challenge;
                }
            }
        }
    }

    /// Append the commitments and scalars from each batch of claims to
    /// Shplemini vectors that will subsequently be inputs to batch mul.
    /// Update the batched evaluation and the running batching challenge
    /// (power of rho) in place.
    pub fn update_batch_mul_inputs_and_batched_evaluation(
        &self,
        commitments: &mut Vec<AffineElement<C>>,
        scalars: &mut Vec<Field<C::ScalarFieldParams>>,
        batched_evaluation: &mut Field<C::ScalarFieldParams>,
        rho: &Field<C::ScalarFieldParams>,
        shplonk_batching_pos: Field<C::ScalarFieldParams>,
        shplonk_batching_neg: Field<C::ScalarFieldParams>,
    ) {
        let mut rho_power = Field::one();

        // Helper closure: aggregate claim data from a batch
        let mut aggregate = |batch: &Batch<C>, rho_power: &mut Field<C::ScalarFieldParams>| {
            for (commitment, evaluation) in
                batch.commitments.iter().zip(batch.evaluations.iter())
            {
                commitments.push(*commitment);
                scalars.push(-batch.scalar * *rho_power);
                *batched_evaluation += *evaluation * *rho_power;
                *rho_power *= *rho;
            }
        };

        if let Some(ref batch) = self.unshifted {
            aggregate(batch, &mut rho_power);
        }

        if let Some(ref batch) = self.shifted {
            aggregate(batch, &mut rho_power);
        }

        if let Some(ref il) = self.interleaved {
            let group_size = self.get_groups_to_be_interleaved_size() as usize;
            assert!(group_size % 2 == 0, "Interleaved groups size must be even");

            let mut group_idx = 0;
            for j in 0..il.commitments_groups.len() {
                for i in 0..group_size {
                    commitments.push(il.commitments_groups[j][i]);
                    scalars.push(
                        -rho_power
                            * il.shplonk_denominator
                            * (shplonk_batching_pos * il.scalars_pos[i]
                                + shplonk_batching_neg * il.scalars_neg[i]),
                    );
                }
                *batched_evaluation += il.evaluations[group_idx] * rho_power;
                if j != il.commitments_groups.len() - 1 {
                    rho_power *= *rho;
                }
                group_idx += 1;
            }
        }
    }
}
