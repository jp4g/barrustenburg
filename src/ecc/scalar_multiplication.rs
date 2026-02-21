// Pippenger Multi-Scalar Multiplication
//
// Ported from C++ `ecc/scalar_multiplication/scalar_multiplication.hpp/cpp`.
//
// Single-threaded implementation of Pippenger's bucket method for computing
// sum(scalar[i] * point[i]). Two paths:
// - Small path (< threshold): Jacobian buckets, simple add_assign_affine
// - Large path (>= threshold): Affine buckets with batch inversion trick

use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;
use crate::ecc::groups::affine_element::AffineElement;
use crate::ecc::groups::curve_params::{BaseField, CurveParams};
use crate::ecc::groups::element::Element;

// ---------------------------------------------------------------------------
// BitVector — packed bit vector for tracking bucket occupancy
// ---------------------------------------------------------------------------

struct BitVector {
    data: Vec<u64>,
}

impl BitVector {
    fn new(num_bits: usize) -> Self {
        let num_words = (num_bits + 63) / 64;
        Self {
            data: vec![0u64; num_words],
        }
    }

    #[inline]
    fn get(&self, idx: usize) -> bool {
        let word = idx >> 6;
        let bit = idx & 63;
        (self.data[word] >> bit) & 1 == 1
    }

    #[inline]
    fn set(&mut self, idx: usize, value: bool) {
        let word = idx >> 6;
        let bit = idx & 63;
        if value {
            self.data[word] |= 1u64 << bit;
        } else {
            self.data[word] &= !(1u64 << bit);
        }
    }

    fn clear(&mut self) {
        for w in self.data.iter_mut() {
            *w = 0;
        }
    }
}

// ---------------------------------------------------------------------------
// get_scalar_slice — extract bits from a 4-limb scalar
// ---------------------------------------------------------------------------

/// Extract `slice_size` bits from a 4-limb scalar at round `round`.
/// Scalar must NOT be in Montgomery form. Bits are extracted MSB-first:
/// at round 0 we take the topmost slice, at round (num_rounds-1) the lowest.
#[inline]
fn get_scalar_slice(scalar: &[u64; 4], round: usize, slice_size: usize, num_bits: usize) -> u32 {
    let hi_bit = num_bits - round * slice_size;
    let last_slice = hi_bit < slice_size;
    let target_slice_size = if last_slice { hi_bit } else { slice_size };
    let lo_bit = if last_slice { 0 } else { hi_bit - slice_size };

    let start_limb = lo_bit / 64;
    let end_limb = hi_bit / 64;
    let lo_slice_offset = lo_bit & 63;
    let lo_slice_bits = target_slice_size.min(64 - lo_slice_offset);
    let hi_slice_bits = target_slice_size - lo_slice_bits;

    let lo_mask = if lo_slice_bits == 64 {
        u64::MAX
    } else {
        (1u64 << lo_slice_bits) - 1
    };
    let lo_slice = (scalar[start_limb] >> lo_slice_offset) & lo_mask;

    let hi_mask = if hi_slice_bits == 0 {
        0u64
    } else if hi_slice_bits == 64 {
        u64::MAX
    } else {
        (1u64 << hi_slice_bits) - 1
    };
    let hi_slice = if end_limb < 4 {
        scalar[end_limb] & hi_mask
    } else {
        0
    };

    (lo_slice as u32) + ((hi_slice as u32) << lo_slice_bits)
}

/// Compute number of significant bits in the scalar field modulus.
fn num_bits_in_scalar_field<C: CurveParams>() -> usize {
    let modulus = C::ScalarFieldParams::MODULUS;
    for i in (0..4).rev() {
        if modulus[i] != 0 {
            return (i * 64) + (64 - modulus[i].leading_zeros() as usize);
        }
    }
    0
}

// ---------------------------------------------------------------------------
// Radix sort — sorts u64 schedule entries by lower 32 bits
// ---------------------------------------------------------------------------

/// 4-pass 8-bit radix sort on the lower 32 bits of u64 entries.
/// Returns the number of entries whose lower 32 bits are zero.
fn radix_sort(schedule: &mut [u64]) -> usize {
    if schedule.is_empty() {
        return 0;
    }
    let n = schedule.len();
    let mut scratch = vec![0u64; n];

    // 4 passes of 8 bits each, covering the lower 32 bits
    for pass in 0..4u32 {
        let shift = pass * 8;
        let mut counts = [0usize; 256];

        // Count
        for &entry in schedule.iter() {
            let byte = ((entry >> shift) & 0xFF) as usize;
            counts[byte] += 1;
        }

        // Prefix sum (exclusive)
        let mut offsets = [0usize; 256];
        for i in 1..256 {
            offsets[i] = offsets[i - 1] + counts[i - 1];
        }

        // Scatter into scratch
        for &entry in schedule.iter() {
            let byte = ((entry >> shift) & 0xFF) as usize;
            scratch[offsets[byte]] = entry;
            offsets[byte] += 1;
        }

        // Copy back
        schedule.copy_from_slice(&scratch);
    }

    // Count zero entries (lower 32 bits == 0)
    let mut num_zero = 0;
    for &entry in schedule.iter() {
        if (entry & 0xFFFFFFFF) == 0 {
            num_zero += 1;
        } else {
            break; // sorted, so zeros are at the front
        }
    }
    num_zero
}

// ---------------------------------------------------------------------------
// add_affine_points — batch pairwise affine addition via Montgomery inversion
// ---------------------------------------------------------------------------

/// Given `points[0..num_points]` arranged as pairs (points[0]+points[1], points[2]+points[3], ...),
/// compute all pairwise sums using a single field inversion.
/// Results are stored in `points[num_points/2..]`.
///
/// Matches C++ `MSM::add_affine_points`.
fn add_affine_points<C: CurveParams>(
    points: &mut [AffineElement<C>],
    num_points: usize,
    scratch_space: &mut [BaseField<C>],
) {
    let mut accumulator = BaseField::<C>::one();

    // Forward pass: accumulate (x2 - x1) products, storing intermediate state
    for i in (0..num_points).step_by(2) {
        scratch_space[i >> 1] = points[i].x + points[i + 1].x; // x2 + x1 (saved for later)
        points[i + 1].x = points[i + 1].x - points[i].x; // x2 - x1 (dx)
        points[i + 1].y = points[i + 1].y - points[i].y; // y2 - y1 (dy)
        points[i + 1].y = points[i + 1].y * accumulator; // dy * accumulator_old
        accumulator = accumulator * points[i + 1].x; // accumulator *= dx
    }

    debug_assert!(
        !accumulator.is_zero(),
        "attempted to invert zero in add_affine_points"
    );
    accumulator = accumulator.invert();

    // Backward pass: extract individual inverses and compute additions
    let half = num_points / 2;
    let mut i = num_points - 2;
    loop {
        // lambda = dy * (1/dx) = (points[i+1].y * accumulator_remaining)
        points[i + 1].y = points[i + 1].y * accumulator;
        accumulator = accumulator * points[i + 1].x; // restore: accumulator *= dx

        // x3 = lambda^2 - (x1 + x2)
        points[i + 1].x = points[i + 1].y.sqr();
        points[half + (i >> 1)].x = points[i + 1].x - scratch_space[i >> 1];

        // y3 = lambda * (x1 - x3) - y1
        points[i].x = points[i].x - points[half + (i >> 1)].x;
        points[i].x = points[i].x * points[i + 1].y;
        points[half + (i >> 1)].y = points[i].x - points[i].y;

        if i < 2 {
            break;
        }
        i -= 2;
    }
}

// ---------------------------------------------------------------------------
// consume_point_schedule — fill affine buckets with batch addition
// ---------------------------------------------------------------------------

const BATCH_SIZE: usize = 2048;
const BATCH_OVERFLOW: usize = 2;

/// Iteratively add points into affine buckets using the batch inversion trick.
fn consume_point_schedule<C: CurveParams>(
    schedule: &[u64],
    points: &[AffineElement<C>],
    buckets: &mut [AffineElement<C>],
    bitvector: &mut BitVector,
) {
    let num_points = schedule.len();
    if num_points == 0 {
        return;
    }

    let mut scratch_points: Vec<AffineElement<C>> = vec![AffineElement::infinity(); BATCH_SIZE + BATCH_OVERFLOW];
    let mut scratch_fields: Vec<BaseField<C>> =
        vec![BaseField::<C>::zero(); (BATCH_SIZE + BATCH_OVERFLOW) / 2];
    let mut bucket_destinations: Vec<usize> = vec![0; (BATCH_SIZE + BATCH_OVERFLOW) / 2];

    let mut point_it: usize = 0;
    let mut affine_input_it: usize = 0;

    loop {
        // Step 1: Fill scratch space with pairs to add
        while (affine_input_it + 1) < BATCH_SIZE && point_it < num_points.saturating_sub(1) {
            let lhs_schedule = schedule[point_it];
            let rhs_schedule = schedule[point_it + 1];
            let lhs_bucket = (lhs_schedule & 0xFFFFFFFF) as usize;
            let rhs_bucket = (rhs_schedule & 0xFFFFFFFF) as usize;
            let lhs_point = (lhs_schedule >> 32) as usize;
            let rhs_point = (rhs_schedule >> 32) as usize;

            let has_bucket = bitvector.get(lhs_bucket);
            let buckets_match = lhs_bucket == rhs_bucket;
            let do_affine_add = buckets_match || has_bucket;

            if do_affine_add {
                scratch_points[affine_input_it] = points[lhs_point];
                scratch_points[affine_input_it + 1] = if buckets_match {
                    points[rhs_point]
                } else {
                    buckets[lhs_bucket]
                };
                bucket_destinations[affine_input_it >> 1] = lhs_bucket;

                // Update bitvector: bucket stays occupied only if it had an accumulator AND both
                // points went to the same bucket (so the accumulator wasn't consumed)
                bitvector.set(lhs_bucket, has_bucket && buckets_match);

                affine_input_it += 2;
                point_it += if buckets_match { 2 } else { 1 };
            } else {
                // No affine add: cache point into bucket
                buckets[lhs_bucket] = points[lhs_point];
                bitvector.set(lhs_bucket, true);
                point_it += 1;
            }
        }

        // Handle last remaining point
        if point_it == num_points - 1 {
            let lhs_schedule = schedule[point_it];
            let lhs_bucket = (lhs_schedule & 0xFFFFFFFF) as usize;
            let lhs_point = (lhs_schedule >> 32) as usize;
            let has_bucket = bitvector.get(lhs_bucket);

            if has_bucket {
                scratch_points[affine_input_it] = points[lhs_point];
                scratch_points[affine_input_it + 1] = buckets[lhs_bucket];
                bitvector.set(lhs_bucket, false);
                bucket_destinations[affine_input_it >> 1] = lhs_bucket;
                affine_input_it += 2;
            } else {
                buckets[lhs_bucket] = points[lhs_point];
                bitvector.set(lhs_bucket, true);
            }
            point_it += 1;
        }

        // Step 2: Batch add all queued pairs
        let num_to_add = affine_input_it;
        if num_to_add >= 2 {
            add_affine_points::<C>(&mut scratch_points, num_to_add, &mut scratch_fields);
        }

        // Step 3: Process addition results — feed back into buckets or re-queue
        let half = num_to_add / 2;
        let mut new_scratch_it: usize = 0;
        let mut output_it: usize = 0;

        while output_it < half.saturating_sub(1) {
            let lhs_bucket = bucket_destinations[output_it];
            let rhs_bucket = bucket_destinations[output_it + 1];
            let has_bucket = bitvector.get(lhs_bucket);
            let buckets_match = lhs_bucket == rhs_bucket;
            let do_affine_add = buckets_match || has_bucket;

            if do_affine_add {
                scratch_points[new_scratch_it] = scratch_points[half + output_it];
                scratch_points[new_scratch_it + 1] = if buckets_match {
                    scratch_points[half + output_it + 1]
                } else {
                    buckets[lhs_bucket]
                };
                bucket_destinations[new_scratch_it >> 1] = lhs_bucket;
                bitvector.set(lhs_bucket, has_bucket && buckets_match);
                new_scratch_it += 2;
                output_it += if buckets_match { 2 } else { 1 };
            } else {
                buckets[lhs_bucket] = scratch_points[half + output_it];
                bitvector.set(lhs_bucket, true);
                output_it += 1;
            }
        }

        // Handle last output point
        if output_it == half.saturating_sub(1) && half > 0 {
            let lhs_bucket = bucket_destinations[output_it];
            let has_bucket = bitvector.get(lhs_bucket);
            if has_bucket {
                scratch_points[new_scratch_it] = scratch_points[half + output_it];
                scratch_points[new_scratch_it + 1] = buckets[lhs_bucket];
                bitvector.set(lhs_bucket, false);
                bucket_destinations[new_scratch_it >> 1] = lhs_bucket;
                new_scratch_it += 2;
            } else {
                buckets[lhs_bucket] = scratch_points[half + output_it];
                bitvector.set(lhs_bucket, true);
            }
        }

        affine_input_it = new_scratch_it;

        // If we've processed all input points and have no queued additions, we're done
        if point_it >= num_points && new_scratch_it == 0 {
            break;
        }
    }
}

// ---------------------------------------------------------------------------
// accumulate_buckets — prefix-sum reduction of bucket array
// ---------------------------------------------------------------------------

/// Compute sum(i * bucket[i]) for i=1..num_buckets using prefix-sum approach.
/// Works for both affine and Jacobian bucket arrays.
fn accumulate_affine_buckets<C: CurveParams>(
    buckets: &[AffineElement<C>],
    bitvector: &BitVector,
    num_buckets: usize,
) -> Element<C> {
    // Find highest non-empty bucket
    let mut starting_index = num_buckets - 1;
    while starting_index > 0 && !bitvector.get(starting_index) {
        starting_index -= 1;
    }
    if starting_index == 0 {
        return Element::infinity();
    }

    let mut prefix_sum = Element::from_affine(&buckets[starting_index]);
    let mut sum = prefix_sum;

    for i in (1..starting_index).rev() {
        if bitvector.get(i) {
            prefix_sum.add_assign_affine(&buckets[i]);
        }
        sum = sum + prefix_sum;
    }
    sum
}

fn accumulate_jacobian_buckets<C: CurveParams>(
    buckets: &[Element<C>],
    bitvector: &BitVector,
    num_buckets: usize,
) -> Element<C> {
    let mut starting_index = num_buckets - 1;
    while starting_index > 0 && !bitvector.get(starting_index) {
        starting_index -= 1;
    }
    if starting_index == 0 {
        return Element::infinity();
    }

    let mut prefix_sum = buckets[starting_index];
    let mut sum = prefix_sum;

    for i in (1..starting_index).rev() {
        if bitvector.get(i) {
            prefix_sum = prefix_sum + buckets[i];
        }
        sum = sum + prefix_sum;
    }
    sum
}

// ---------------------------------------------------------------------------
// evaluate_small_pippenger_round — Jacobian path
// ---------------------------------------------------------------------------

fn evaluate_small_pippenger_round<C: CurveParams>(
    scalars: &[[u64; 4]],
    points: &[AffineElement<C>],
    nonzero_indices: &[usize],
    round_index: usize,
    jacobian_buckets: &mut [Element<C>],
    bucket_bitvector: &mut BitVector,
    previous_result: Element<C>,
    bits_per_slice: usize,
    num_bits: usize,
) -> Element<C> {
    let num_buckets = jacobian_buckets.len();

    // Add points to Jacobian buckets
    for &idx in nonzero_indices {
        let bucket_index =
            get_scalar_slice(&scalars[idx], round_index, bits_per_slice, num_bits) as usize;
        debug_assert!(bucket_index < num_buckets);
        if bucket_index > 0 {
            if bucket_bitvector.get(bucket_index) {
                jacobian_buckets[bucket_index].add_assign_affine(&points[idx]);
            } else {
                jacobian_buckets[bucket_index] = Element::from_affine(&points[idx]);
                bucket_bitvector.set(bucket_index, true);
            }
        }
    }

    let round_output = accumulate_jacobian_buckets::<C>(jacobian_buckets, bucket_bitvector, num_buckets);
    bucket_bitvector.clear();

    // Double previous result by appropriate number of bits
    let num_rounds = (num_bits + bits_per_slice - 1) / bits_per_slice;
    let num_doublings = if round_index == num_rounds - 1 && num_bits % bits_per_slice != 0 {
        num_bits % bits_per_slice
    } else {
        bits_per_slice
    };

    let mut result = previous_result;
    for _ in 0..num_doublings {
        result.self_dbl();
    }
    result = result + round_output;
    result
}

// ---------------------------------------------------------------------------
// evaluate_pippenger_round — affine trick path
// ---------------------------------------------------------------------------

fn evaluate_pippenger_round<C: CurveParams>(
    scalars: &[[u64; 4]],
    points: &[AffineElement<C>],
    nonzero_indices: &[usize],
    round_index: usize,
    affine_buckets: &mut [AffineElement<C>],
    bucket_bitvector: &mut BitVector,
    schedule_buf: &mut [u64],
    previous_result: Element<C>,
    bits_per_slice: usize,
    num_bits: usize,
) -> Element<C> {
    let size = nonzero_indices.len();

    // Build round schedule: (point_index << 32) | bucket_index
    for i in 0..size {
        let bucket_index =
            get_scalar_slice(&scalars[nonzero_indices[i]], round_index, bits_per_slice, num_bits)
                as u64;
        schedule_buf[i] = bucket_index | ((nonzero_indices[i] as u64) << 32);
    }

    // Sort by bucket index (lower 32 bits) and count zeros
    let work_schedule = &mut schedule_buf[..size];
    let num_zero = radix_sort(work_schedule);
    let round_size = size - num_zero;

    let round_output = if round_size > 0 {
        let point_schedule = &work_schedule[num_zero..];
        consume_point_schedule::<C>(point_schedule, points, affine_buckets, bucket_bitvector);
        let result =
            accumulate_affine_buckets::<C>(affine_buckets, bucket_bitvector, affine_buckets.len());
        bucket_bitvector.clear();
        result
    } else {
        Element::infinity()
    };

    // Double previous result
    let num_rounds = (num_bits + bits_per_slice - 1) / bits_per_slice;
    let num_doublings = if round_index == num_rounds - 1 && num_bits % bits_per_slice != 0 {
        num_bits % bits_per_slice
    } else {
        bits_per_slice
    };

    let mut result = previous_result;
    for _ in 0..num_doublings {
        result.self_dbl();
    }
    result = result + round_output;
    result
}

// ---------------------------------------------------------------------------
// use_affine_trick — decide whether batch affine is beneficial
// ---------------------------------------------------------------------------

fn use_affine_trick(num_points: usize, num_buckets: usize, num_bits: usize) -> bool {
    if num_points < 128 {
        return false;
    }
    let cost_of_inversion = num_bits + (num_bits + 3) / 4 + 14;
    let cost_saving_per_op: usize = 5;
    let extra_jacobian_cost: usize = 5;

    let group_op_saving =
        num_points * cost_saving_per_op + num_buckets * extra_jacobian_cost;
    let inversion_cost = ((num_points as f64).log2() * cost_of_inversion as f64) as usize;

    group_op_saving > inversion_cost
}

// ---------------------------------------------------------------------------
// pippenger_msm — main entry point
// ---------------------------------------------------------------------------

/// Compute sum(scalars[i] * points[i]) using Pippenger's bucket method.
pub fn pippenger_msm<C: CurveParams>(
    scalars: &[Field<C::ScalarFieldParams>],
    points: &[AffineElement<C>],
) -> Element<C> {
    assert_eq!(scalars.len(), points.len());
    let n = scalars.len();
    if n == 0 {
        return Element::infinity();
    }

    // Convert scalars out of Montgomery form and collect nonzero indices
    let raw_scalars: Vec<[u64; 4]> = scalars
        .iter()
        .map(|s| s.from_montgomery_form().data)
        .collect();

    let nonzero_indices: Vec<usize> = raw_scalars
        .iter()
        .enumerate()
        .filter(|(_, s)| s[0] != 0 || s[1] != 0 || s[2] != 0 || s[3] != 0)
        .map(|(i, _)| i)
        .collect();

    let num_nonzero = nonzero_indices.len();

    // Trivial cases
    if num_nonzero == 0 {
        return Element::infinity();
    }
    if num_nonzero == 1 {
        let idx = nonzero_indices[0];
        return Element::from_affine(&points[idx]).mul_without_endomorphism(&scalars[idx]);
    }

    let num_bits = num_bits_in_scalar_field::<C>();
    let bits_per_slice = get_optimal_log_num_buckets(num_nonzero, num_bits);
    let num_buckets = 1usize << bits_per_slice;
    let num_rounds = (num_bits + bits_per_slice - 1) / bits_per_slice;

    if use_affine_trick(num_nonzero, num_buckets, num_bits) {
        // Affine trick path
        let mut affine_buckets = vec![AffineElement::<C>::infinity(); num_buckets];
        let mut bitvector = BitVector::new(num_buckets);
        let mut schedule_buf = vec![0u64; num_nonzero];
        let mut result = Element::infinity();

        for round in 0..num_rounds {
            result = evaluate_pippenger_round::<C>(
                &raw_scalars,
                points,
                &nonzero_indices,
                round,
                &mut affine_buckets,
                &mut bitvector,
                &mut schedule_buf,
                result,
                bits_per_slice,
                num_bits,
            );
        }
        result
    } else {
        // Jacobian path (small point count)
        let mut jacobian_buckets = vec![Element::<C>::infinity(); num_buckets];
        let mut bitvector = BitVector::new(num_buckets);
        let mut result = Element::infinity();

        for round in 0..num_rounds {
            result = evaluate_small_pippenger_round::<C>(
                &raw_scalars,
                points,
                &nonzero_indices,
                round,
                &mut jacobian_buckets,
                &mut bitvector,
                result,
                bits_per_slice,
                num_bits,
            );
        }
        result
    }
}

/// Heuristic cost-model for optimal bucket width (matches C++ `get_optimal_log_num_buckets`).
fn get_optimal_log_num_buckets(num_points: usize, num_bits: usize) -> usize {
    let mut cached_cost = usize::MAX;
    let mut target = 1;
    for bit_slice in 1..20 {
        let num_rounds = (num_bits + bit_slice - 1) / bit_slice;
        let num_buckets = 1usize << bit_slice;
        let addition_cost = num_rounds * num_points;
        let bucket_cost = num_rounds * num_buckets * 5;
        let total_cost = addition_cost + bucket_cost;
        if total_cost < cached_cost {
            cached_cost = total_cost;
            target = bit_slice;
        }
    }
    target
}

/// Naive MSM fallback: sum of individual scalar multiplications.
pub fn naive_msm<C: CurveParams>(
    scalars: &[Field<C::ScalarFieldParams>],
    points: &[AffineElement<C>],
) -> Element<C> {
    assert_eq!(scalars.len(), points.len());
    let mut acc = Element::<C>::infinity();
    for (s, p) in scalars.iter().zip(points.iter()) {
        let proj = Element::from_affine(p);
        acc += proj.mul_without_endomorphism(s);
    }
    acc
}

// ---------------------------------------------------------------------------
// Internal test helpers
// ---------------------------------------------------------------------------

#[cfg(test)]
pub(crate) fn get_scalar_slice_test(
    scalar: &[u64; 4],
    round: usize,
    slice_size: usize,
    num_bits: usize,
) -> u32 {
    get_scalar_slice(scalar, round, slice_size, num_bits)
}

#[cfg(test)]
pub(crate) fn radix_sort_test(schedule: &mut [u64]) -> usize {
    radix_sort(schedule)
}
