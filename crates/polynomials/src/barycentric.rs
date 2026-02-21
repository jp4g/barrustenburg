use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

/// Precomputed tables for barycentric interpolation/extension on domain {0, 1, ..., domain_size-1}.
///
/// Port of C++ `BarycentricDataRunTime<Fr, domain_end, num_evals>`.
pub struct BarycentricData<P: FieldParams> {
    pub domain_size: usize,
    pub num_evals: usize,
    pub big_domain_size: usize,
    pub big_domain: Vec<Field<P>>,
    pub lagrange_denominators: Vec<Field<P>>,
    /// Row-major: index `k * domain_size + j`. For k < domain_size, entries come from
    /// batch-inverting the lagrange_denominators (which are zero-padded for those rows).
    pub precomputed_denominator_inverses: Vec<Field<P>>,
    /// M(k) = prod_{j=0}^{domain_size-1} (k - j) for k in [0, num_evals).
    /// Zero for k < domain_size.
    pub full_numerator_values: Vec<Field<P>>,
}

impl<P: FieldParams> BarycentricData<P> {
    pub fn new(domain_size: usize, num_evals: usize) -> Self {
        let big_domain_size = std::cmp::max(domain_size, num_evals);

        // big_domain = [0, 1, 2, ..., big_domain_size - 1]
        let big_domain: Vec<Field<P>> = (0..big_domain_size)
            .map(|i| Field::from(i as u64))
            .collect();

        // lagrange_denominators: d_i = prod_{j != i} (i - j) for i in [0, domain_size)
        let mut lagrange_denominators = vec![Field::<P>::one(); domain_size];
        for i in 0..domain_size {
            for j in 0..domain_size {
                if j != i {
                    lagrange_denominators[i] =
                        lagrange_denominators[i] * (big_domain[i] - big_domain[j]);
                }
            }
        }

        // precomputed_denominator_inverses: for k in [domain_size, num_evals), j in [0, domain_size):
        //   result[k * domain_size + j] = lagrange_denominators[j] * (big_domain[k] - big_domain[j])
        // then batch invert the whole array.
        let total = domain_size * num_evals;
        let mut denom_inverses = vec![Field::<P>::zero(); total];

        if num_evals > 1 {
            for k in domain_size..num_evals {
                for j in 0..domain_size {
                    let mut inv = lagrange_denominators[j];
                    inv = inv * (big_domain[k] - big_domain[j]);
                    denom_inverses[k * domain_size + j] = inv;
                }
            }
        }

        Field::batch_invert(&mut denom_inverses);

        // full_numerator_values: M(k) = prod_{j=0}^{domain_size-1} (k - j)
        let mut full_numerator_values = vec![Field::<P>::one(); num_evals];
        for i in 0..num_evals {
            let v_i = Field::<P>::from(i as u64);
            for j in 0..domain_size {
                full_numerator_values[i] = full_numerator_values[i] * (v_i - big_domain[j]);
            }
        }

        Self {
            domain_size,
            num_evals,
            big_domain_size,
            big_domain,
            lagrange_denominators,
            precomputed_denominator_inverses: denom_inverses,
            full_numerator_values,
        }
    }
}
