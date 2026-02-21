// Batched Affine Addition — reduce sequences of affine points via repeated pairwise addition.
//
// Ported from C++ `ecc/batched_affine_addition/batched_affine_addition.hpp/cpp`.
//
// Given a flat array of points and a list of sequence lengths, reduces each sequence
// to a single point using repeated rounds of pairwise affine addition. All inversions
// within a round are batched via Montgomery's trick (1 inversion + 3N muls for N pairs).
//
// Single-threaded implementation (C++ uses parallel_for for thread distribution).

use crate::groups::affine_element::AffineElement;
use crate::groups::curve_params::{BaseField, CurveParams};

/// Add two affine points given a pre-computed denominator = 1/(x2 - x1).
#[inline]
fn affine_add_with_denominator<C: CurveParams>(
    p1: &AffineElement<C>,
    p2: &AffineElement<C>,
    denominator: &BaseField<C>,
) -> AffineElement<C> {
    let lambda = *denominator * (p2.y - p1.y);
    let x3 = lambda.sqr() - p2.x - p1.x;
    let y3 = lambda * (p1.x - x3) - p1.y;
    AffineElement::new(x3, y3)
}

/// Batch compute 1/(x2 - x1) for all pairs in the given sequences.
///
/// Returns a Vec of denominators, one per pair.
fn batch_compute_slope_inverses<C: CurveParams>(
    points: &[AffineElement<C>],
    sequence_counts: &[usize],
) -> Vec<BaseField<C>> {
    // Count total pairs
    let total_pairs: usize = sequence_counts.iter().map(|&c| c >> 1).sum();
    if total_pairs == 0 {
        return Vec::new();
    }

    let mut denominators: Vec<BaseField<C>> = Vec::with_capacity(total_pairs);
    let mut differences: Vec<BaseField<C>> = Vec::with_capacity(total_pairs);

    // Forward pass: accumulate products of (x2 - x1)
    let mut accumulator = BaseField::<C>::one();
    let mut point_idx = 0;

    for &count in sequence_counts {
        let num_pairs = count >> 1;
        for _ in 0..num_pairs {
            let x1 = points[point_idx].x;
            let x2 = points[point_idx + 1].x;
            point_idx += 2;

            let diff = x2 - x1;
            differences.push(diff);
            denominators.push(accumulator);
            accumulator = accumulator * diff;
        }
        // Skip unpaired point (odd count)
        if count & 1 == 1 {
            point_idx += 1;
        }
    }

    // Invert the full product
    let mut inverse = accumulator.invert();

    // Backward pass: extract individual denominators
    for i in (0..total_pairs).rev() {
        denominators[i] = denominators[i] * inverse;
        inverse = inverse * differences[i];
    }

    denominators
}

/// Reduce each sequence of affine points to a single point via repeated pairwise addition.
///
/// `points` is a flat array of points grouped by sequence. `sequence_counts[i]` gives the
/// length of the i-th sequence. Points are modified in place; after return, the first
/// `sequence_counts.len()` positions contain the reduced points.
///
/// Returns a Vec of reduced points (one per sequence).
pub fn batched_affine_add_in_place<C: CurveParams>(
    points: &mut Vec<AffineElement<C>>,
    sequence_counts: &mut Vec<usize>,
) -> Vec<AffineElement<C>> {
    batched_affine_add_recursive::<C>(points, sequence_counts);

    // After full reduction, each sequence is length 1. Extract results.
    let mut results = Vec::with_capacity(sequence_counts.len());
    let mut idx = 0;
    for &count in sequence_counts.iter() {
        debug_assert_eq!(count, 1);
        results.push(points[idx]);
        idx += count;
    }
    results
}

/// Internal recursive pairwise reduction.
fn batched_affine_add_recursive<C: CurveParams>(
    points: &mut [AffineElement<C>],
    sequence_counts: &mut [usize],
) {
    let total_points: usize = sequence_counts.iter().sum();
    if total_points <= 1 {
        return;
    }

    // Check if any sequence needs reduction
    let needs_reduction = sequence_counts.iter().any(|&c| c > 1);
    if !needs_reduction {
        return;
    }

    // Batch compute all 1/(x2 - x1) denominators
    let denominators = batch_compute_slope_inverses::<C>(points, sequence_counts);

    // Pairwise addition: write results compactly into the points array
    let mut point_idx = 0;
    let mut result_idx = 0;
    let mut pair_idx = 0;
    let mut more_additions = false;

    for count in sequence_counts.iter_mut() {
        let num_pairs = *count >> 1;
        let overflow = *count & 1 == 1;

        for _ in 0..num_pairs {
            let p1 = points[point_idx];
            let p2 = points[point_idx + 1];
            let denom = &denominators[pair_idx];
            points[result_idx] = affine_add_with_denominator::<C>(&p1, &p2, denom);
            point_idx += 2;
            result_idx += 1;
            pair_idx += 1;
        }

        // Carry forward unpaired point
        if overflow {
            points[result_idx] = points[point_idx];
            point_idx += 1;
            result_idx += 1;
        }

        let updated_count = num_pairs + overflow as usize;
        *count = updated_count;
        more_additions = more_additions || updated_count > 1;
    }

    // Recurse on the compacted array
    if more_additions {
        batched_affine_add_recursive::<C>(&mut points[..result_idx], sequence_counts);
    }
}
