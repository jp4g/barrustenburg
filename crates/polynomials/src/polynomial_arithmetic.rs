use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use crate::evaluation_domain::EvaluationDomain;

// ═══════════════════════════════════════════════════════════════════════════════
// Bit-reversal
// ═══════════════════════════════════════════════════════════════════════════════

/// Reverse the lower `bit_length` bits of `x`.
/// Matches C++ `reverse_bits` (polynomial_arithmetic.cpp:34-41).
#[inline]
fn reverse_bits(x: u32, bit_length: u32) -> u32 {
    let mut x = x;
    x = ((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1);
    x = ((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2);
    x = ((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4);
    x = ((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8);
    ((x >> 16) | (x << 16)) >> (32 - bit_length)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Core FFT engines
// ═══════════════════════════════════════════════════════════════════════════════

/// Single-array Cooley-Tukey FFT. Reads from `coeffs`, writes to `target`.
/// `round_roots` is the twiddle-factor lookup table (one vec per butterfly round).
fn fft_inner<P: FieldParams>(
    coeffs: &[Field<P>],
    target: &mut [Field<P>],
    round_roots: &[Vec<Field<P>>],
    n: usize,
) {
    let log2_n = n.trailing_zeros();

    // Stage 1: bit-reversal + first butterfly
    for i in 0..n / 2 {
        let swap_1 = reverse_bits((2 * i) as u32, log2_n) as usize;
        let swap_2 = reverse_bits((2 * i + 1) as u32, log2_n) as usize;
        target[2 * i] = coeffs[swap_1] + coeffs[swap_2];
        target[2 * i + 1] = coeffs[swap_1] - coeffs[swap_2];
    }

    // Butterfly rounds
    let mut m = 2usize;
    while m < n {
        let round_idx = m.trailing_zeros() as usize;
        let block_mask = m - 1;
        let index_mask = !block_mask;

        for i in 0..n / 2 {
            let k1 = (i & index_mask) << 1;
            let j1 = i & block_mask;
            let temp = round_roots[round_idx][j1] * target[k1 + j1 + m];
            target[k1 + j1 + m] = target[k1 + j1] - temp;
            target[k1 + j1] = target[k1 + j1] + temp;
        }

        m <<= 1;
    }
}

/// Multi-polynomial interleaved Cooley-Tukey FFT.
/// `polys` contains sub-arrays whose concatenation forms the full polynomial.
/// Element at logical index `i` lives at `polys[i / poly_size][i % poly_size]`.
fn fft_inner_split<P: FieldParams>(
    polys: &mut [&mut [Field<P>]],
    round_roots: &[Vec<Field<P>>],
    n: usize,
) {
    let num_polys = polys.len();
    let poly_size = n / num_polys;
    let log2_poly_size = poly_size.trailing_zeros() as usize;
    let poly_mask = poly_size - 1;
    let log2_n = n.trailing_zeros();

    let mut scratch = vec![Field::<P>::zero(); n];

    // Stage 1: bit-reversal + first butterfly
    for i in 0..n / 2 {
        let swap_1 = reverse_bits((2 * i) as u32, log2_n) as usize;
        let swap_2 = reverse_bits((2 * i + 1) as u32, log2_n) as usize;

        let val_1 = polys[swap_1 >> log2_poly_size][swap_1 & poly_mask];
        let val_2 = polys[swap_2 >> log2_poly_size][swap_2 & poly_mask];

        scratch[2 * i] = val_1 + val_2;
        scratch[2 * i + 1] = val_1 - val_2;
    }

    // Handle tiny domain (n <= 2)
    if n <= 2 {
        for i in 0..n {
            polys[i >> log2_poly_size][i & poly_mask] = scratch[i];
        }
        return;
    }

    // Butterfly rounds
    let mut m = 2usize;
    while m < n {
        let round_idx = m.trailing_zeros() as usize;
        let block_mask = m - 1;
        let index_mask = !block_mask;

        if m != n >> 1 {
            // Non-final round: work in scratch
            for i in 0..n / 2 {
                let k1 = (i & index_mask) << 1;
                let j1 = i & block_mask;
                let temp = round_roots[round_idx][j1] * scratch[k1 + j1 + m];
                scratch[k1 + j1 + m] = scratch[k1 + j1] - temp;
                scratch[k1 + j1] = scratch[k1 + j1] + temp;
            }
        } else {
            // Final round: write back to polys
            for i in 0..n / 2 {
                let k1 = (i & index_mask) << 1;
                let j1 = i & block_mask;

                let idx_1 = k1 + j1;
                let idx_2 = k1 + j1 + m;

                let temp = round_roots[round_idx][j1] * scratch[idx_2];
                polys[idx_2 >> log2_poly_size][idx_2 & poly_mask] = scratch[idx_1] - temp;
                polys[idx_1 >> log2_poly_size][idx_1 & poly_mask] = scratch[idx_1] + temp;
            }
        }

        m <<= 1;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Public FFT / IFFT functions
// ═══════════════════════════════════════════════════════════════════════════════

/// In-place forward FFT.
pub fn fft<P: FieldParams>(coeffs: &mut [Field<P>], domain: &EvaluationDomain<P>) {
    let n = domain.size;
    let mut target = vec![Field::<P>::zero(); n];
    fft_inner(coeffs, &mut target, domain.get_round_roots(), n);
    coeffs[..n].copy_from_slice(&target);
}

/// Forward FFT with separate output buffer.
pub fn fft_with_target<P: FieldParams>(
    coeffs: &[Field<P>],
    target: &mut [Field<P>],
    domain: &EvaluationDomain<P>,
) {
    fft_inner(coeffs, target, domain.get_round_roots(), domain.size);
}

/// Split-array forward FFT.
pub fn fft_split<P: FieldParams>(
    polys: &mut [&mut [Field<P>]],
    domain: &EvaluationDomain<P>,
) {
    fft_inner_split(polys, domain.get_round_roots(), domain.size);
}

/// In-place inverse FFT.
pub fn ifft<P: FieldParams>(coeffs: &mut [Field<P>], domain: &EvaluationDomain<P>) {
    let n = domain.size;
    let mut target = vec![Field::<P>::zero(); n];
    fft_inner(coeffs, &mut target, domain.get_inverse_round_roots(), n);
    for i in 0..n {
        coeffs[i] = target[i] * domain.domain_inverse;
    }
}

/// Inverse FFT with separate output buffer.
pub fn ifft_with_target<P: FieldParams>(
    coeffs: &[Field<P>],
    target: &mut [Field<P>],
    domain: &EvaluationDomain<P>,
) {
    let n = domain.size;
    fft_inner(coeffs, target, domain.get_inverse_round_roots(), n);
    for i in 0..n {
        target[i] = target[i] * domain.domain_inverse;
    }
}

/// Split-array inverse FFT.
pub fn ifft_split<P: FieldParams>(
    polys: &mut [&mut [Field<P>]],
    domain: &EvaluationDomain<P>,
) {
    let n = domain.size;
    let num_polys = polys.len();
    let poly_size = n / num_polys;
    let log2_poly_size = poly_size.trailing_zeros() as usize;
    let poly_mask = poly_size - 1;

    fft_inner_split(polys, domain.get_inverse_round_roots(), n);
    for i in 0..n {
        polys[i >> log2_poly_size][i & poly_mask] =
            polys[i >> log2_poly_size][i & poly_mask] * domain.domain_inverse;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Coset FFT / IFFT operations
// ═══════════════════════════════════════════════════════════════════════════════

/// Scale coefficients by successive powers of a generator:
/// `target[i] = coeffs[i] * generator_start * generator_shift^i`.
fn scale_by_generator<P: FieldParams>(
    coeffs: &[Field<P>],
    target: &mut [Field<P>],
    generator_start: Field<P>,
    generator_shift: Field<P>,
    size: usize,
) {
    let mut work_generator = generator_start;
    for i in 0..size {
        target[i] = coeffs[i] * work_generator;
        work_generator = work_generator * generator_shift;
    }
}

/// In-place coset FFT: scale by generator powers, then FFT.
pub fn coset_fft<P: FieldParams>(coeffs: &mut [Field<P>], domain: &EvaluationDomain<P>) {
    let n = domain.size;
    let mut scaled = vec![Field::<P>::zero(); n];
    scale_by_generator(coeffs, &mut scaled, Field::one(), domain.generator, n);
    let mut target = vec![Field::<P>::zero(); n];
    fft_inner(&scaled, &mut target, domain.get_round_roots(), n);
    coeffs[..n].copy_from_slice(&target);
}

/// Coset FFT with separate output buffer.
pub fn coset_fft_with_target<P: FieldParams>(
    coeffs: &[Field<P>],
    target: &mut [Field<P>],
    domain: &EvaluationDomain<P>,
) {
    let n = domain.size;
    let mut scaled = vec![Field::<P>::zero(); n];
    scale_by_generator(coeffs, &mut scaled, Field::one(), domain.generator, n);
    fft_inner(&scaled, target, domain.get_round_roots(), n);
}

/// Split-array coset FFT.
pub fn coset_fft_split<P: FieldParams>(
    polys: &mut [&mut [Field<P>]],
    domain: &EvaluationDomain<P>,
) {
    let num_polys = polys.len();
    let poly_size = domain.size / num_polys;
    let generator_pow_n = domain.generator.pow(&[poly_size as u64, 0, 0, 0]);
    let mut generator_start = Field::<P>::one();

    for i in 0..num_polys {
        let mut scaled = vec![Field::<P>::zero(); poly_size];
        scale_by_generator(polys[i], &mut scaled, generator_start, domain.generator, poly_size);
        polys[i][..poly_size].copy_from_slice(&scaled);
        generator_start = generator_start * generator_pow_n;
    }
    fft_inner_split(polys, domain.get_round_roots(), domain.size);
}

/// In-place coset inverse FFT: IFFT, then scale by generator_inverse powers.
pub fn coset_ifft<P: FieldParams>(coeffs: &mut [Field<P>], domain: &EvaluationDomain<P>) {
    ifft(coeffs, domain);
    let n = domain.size;
    let mut scaled = vec![Field::<P>::zero(); n];
    scale_by_generator(coeffs, &mut scaled, Field::one(), domain.generator_inverse, n);
    coeffs[..n].copy_from_slice(&scaled);
}

/// Split-array coset inverse FFT.
pub fn coset_ifft_split<P: FieldParams>(
    polys: &mut [&mut [Field<P>]],
    domain: &EvaluationDomain<P>,
) {
    let n = domain.size;
    let num_polys = polys.len();
    let poly_size = n / num_polys;

    ifft_split(polys, domain);

    let generator_inv_pow_n = domain.generator_inverse.pow(&[poly_size as u64, 0, 0, 0]);
    let mut generator_start = Field::<P>::one();

    for i in 0..num_polys {
        let mut scaled = vec![Field::<P>::zero(); poly_size];
        scale_by_generator(
            polys[i],
            &mut scaled,
            generator_start,
            domain.generator_inverse,
            poly_size,
        );
        polys[i][..poly_size].copy_from_slice(&scaled);
        generator_start = generator_start * generator_inv_pow_n;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Helper functions
// ═══════════════════════════════════════════════════════════════════════════════

/// Evaluate a polynomial given in Lagrange basis (evaluations at domain roots of unity)
/// at an arbitrary point `z` using the barycentric formula.
///
/// `coeffs` are the evaluations f(ω^0), f(ω^1), ..., f(ω^{num_coeffs-1}).
pub fn compute_barycentric_evaluation<P: FieldParams>(
    coeffs: &[Field<P>],
    domain: &EvaluationDomain<P>,
    z: &Field<P>,
) -> Field<P> {
    let num_coeffs = coeffs.len();

    // numerator = (z^n - 1) / n
    let mut numerator = *z;
    for _ in 0..domain.log2_size {
        numerator = numerator.sqr();
    }
    numerator = numerator - Field::<P>::one();
    numerator = numerator * domain.domain_inverse;

    // denominators[i] = z * ω^{-i} - 1
    let mut denoms = vec![Field::<P>::zero(); num_coeffs];
    denoms[0] = *z - Field::<P>::one();
    let mut work_root = domain.root_inverse;
    for i in 1..num_coeffs {
        denoms[i] = work_root * *z - Field::<P>::one();
        work_root = work_root * domain.root_inverse;
    }

    Field::<P>::batch_invert(&mut denoms);

    let mut result = Field::<P>::zero();
    for i in 0..num_coeffs {
        result = result + coeffs[i] * denoms[i];
    }

    result * numerator
}

/// Compute the coefficient form of ∏_{i=0}^{n-1} (X - roots[i]).
/// Returns a Vec of length n+1 in ascending coefficient order:
/// `result[0]` = constant term, `result[n]` = 1 (leading coefficient).
pub fn compute_linear_polynomial_product<P: FieldParams>(
    roots: &[Field<P>],
) -> Vec<Field<P>> {
    let n = roots.len();
    let mut dest = vec![Field::<P>::zero(); n + 1];

    let mut scratch = roots.to_vec();

    // Compute sum of roots
    let mut sum = Field::<P>::zero();
    for i in 0..n {
        sum = sum + scratch[i];
    }

    dest[n] = Field::<P>::one();
    dest[n - 1] = -sum;

    let mut constant = Field::<P>::one();
    for i in 0..n.saturating_sub(1) {
        let mut temp = Field::<P>::zero();
        for j in 0..n - 1 - i {
            // scratch[j] = roots[j] * sum(scratch[j+1..n-i])
            let mut partial_sum = Field::<P>::zero();
            for k in (j + 1)..(n - i) {
                partial_sum = partial_sum + scratch[k];
            }
            scratch[j] = roots[j] * partial_sum;
            temp = temp + scratch[j];
        }
        dest[n - 2 - i] = temp * constant;
        constant = -constant;
    }

    dest
}

// ═══════════════════════════════════════════════════════════════════════════════
// Polynomial evaluation / interpolation
// ═══════════════════════════════════════════════════════════════════════════════

/// Evaluate polynomial at `z` using Horner's method.
/// Coefficients are in order [c_0, c_1, ..., c_{n-1}] so p(X) = c_0 + c_1*X + ... + c_{n-1}*X^{n-1}.
pub fn evaluate<P: FieldParams>(coeffs: &[Field<P>], z: &Field<P>) -> Field<P> {
    let n = coeffs.len();
    if n == 0 {
        return Field::zero();
    }
    let mut result = coeffs[n - 1];
    for i in (0..n - 1).rev() {
        result = result * *z + coeffs[i];
    }
    result
}

/// Synthetic division: divide the polynomial (given by `coeffs`) by (X - root) in-place.
/// After the call the leading coefficient is zero and the remaining n-1 coefficients
/// hold the quotient.  Mirrors the C++ `factor_roots`.
pub fn factor_roots<P: FieldParams>(coeffs: &mut [Field<P>], root: &Field<P>) {
    let n = coeffs.len();
    if n == 0 {
        return;
    }
    let mut work = coeffs[n - 1];
    coeffs[n - 1] = Field::zero();
    // Walk from second-highest down to 0 (wrapping arithmetic handles the i == 0 case).
    for i in (0..n - 1).rev() {
        let temp = coeffs[i];
        coeffs[i] = work;
        work = temp + work * *root;
    }
}

/// Lagrange interpolation: given evaluations `evaluations[i]` at `points[i]`, return
/// the coefficient representation of the unique polynomial of degree <= n-1 passing
/// through all points.
pub fn compute_efficient_interpolation<P: FieldParams>(
    evaluations: &[Field<P>],
    points: &[Field<P>],
) -> Vec<Field<P>> {
    let n = evaluations.len();
    assert_eq!(n, points.len());
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![evaluations[0]];
    }

    // Step 1: denominators[i] = prod_{j != i} (points[i] - points[j])
    let mut denoms = vec![Field::<P>::one(); n];
    for i in 0..n {
        for j in 0..n {
            if i != j {
                denoms[i] = denoms[i] * (points[i] - points[j]);
            }
        }
    }

    // Step 2: scaled[i] = evaluations[i] / denominators[i]
    let mut scaled = vec![Field::<P>::zero(); n];
    for i in 0..n {
        scaled[i] = evaluations[i] * denoms[i].invert();
    }

    // Step 3: result = sum_i scaled[i] * prod_{j != i}(X - points[j])
    let mut result = vec![Field::<P>::zero(); n];
    for i in 0..n {
        // Build basis polynomial L_i(X) = prod_{j != i}(X - points[j])
        let mut basis = vec![Field::<P>::zero(); n];
        basis[0] = Field::one();
        let mut deg = 0usize;
        for j in 0..n {
            if j == i {
                continue;
            }
            deg += 1;
            // Multiply current basis by (X - points[j]), high-to-low to avoid overwrites
            for k in (1..=deg).rev() {
                basis[k] = basis[k - 1] - basis[k] * points[j];
            }
            basis[0] = -basis[0] * points[j];
        }
        // Accumulate scaled[i] * basis into result
        for k in 0..n {
            result[k] = result[k] + scaled[i] * basis[k];
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::Bn254FrParams;

    type Fr = Field<Bn254FrParams>;

    #[test]
    fn test_evaluate_empty() {
        let coeffs: Vec<Fr> = vec![];
        let z = Fr::from(5u64);
        assert_eq!(evaluate(&coeffs, &z), Fr::zero());
    }

    #[test]
    fn test_evaluate_constant() {
        let coeffs = vec![Fr::from(7u64)];
        let z = Fr::from(42u64);
        assert_eq!(evaluate(&coeffs, &z), Fr::from(7u64));
    }

    #[test]
    fn test_evaluate_linear() {
        // p(X) = 3 + 2*X, p(5) = 13
        let coeffs = vec![Fr::from(3u64), Fr::from(2u64)];
        let z = Fr::from(5u64);
        assert_eq!(evaluate(&coeffs, &z), Fr::from(13u64));
    }

    #[test]
    fn test_factor_roots() {
        // p(X) = (X - 3)(X - 7) = X^2 - 10X + 21 = 21 - 10X + X^2
        let mut coeffs = vec![Fr::from(21u64), -Fr::from(10u64), Fr::one()];
        factor_roots(&mut coeffs, &Fr::from(3u64));
        // Quotient should be (X - 7) = -7 + X, with leading coeff zeroed
        assert_eq!(coeffs[0], -Fr::from(7u64));
        assert_eq!(coeffs[1], Fr::one());
        assert_eq!(coeffs[2], Fr::zero());
    }

    #[test]
    fn test_interpolation_matches_evaluation() {
        let points = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let evals = vec![Fr::from(6u64), Fr::from(11u64), Fr::from(18u64)];
        let coeffs = compute_efficient_interpolation(&evals, &points);
        // Verify round-trip: evaluating at each point recovers the evaluation
        for i in 0..points.len() {
            assert_eq!(evaluate(&coeffs, &points[i]), evals[i]);
        }
    }
}
