use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;

/// Precomputed root-of-unity tables for NTT/FFT over a prime field.
///
/// Ported from Barretenberg's `EvaluationDomain<FF>` (evaluation_domain.hpp/cpp).
/// Stores the primitive root of unity for a power-of-2 domain, along with
/// twiddle-factor lookup tables for each butterfly round.
pub struct EvaluationDomain<P: FieldParams> {
    /// Domain size (must be a power of 2).
    pub size: usize,
    /// log2(size).
    pub log2_size: usize,
    /// Primitive root of unity of order `size`.
    pub root: Field<P>,
    /// Inverse of `root`.
    pub root_inverse: Field<P>,
    /// `size` as a field element.
    pub domain: Field<P>,
    /// Inverse of `domain`.
    pub domain_inverse: Field<P>,
    /// Multiplicative generator (coset offset), i.e. `coset_generators[0]`.
    pub generator: Field<P>,
    /// Inverse of `generator`.
    pub generator_inverse: Field<P>,
    /// Inverse of 4 as a field element.
    pub four_inverse: Field<P>,
    /// Twiddle-factor tables for forward NTT, one per butterfly round.
    /// `round_roots[i]` has `2^i` entries.
    round_roots: Vec<Vec<Field<P>>>,
    /// Twiddle-factor tables for inverse NTT, one per butterfly round.
    /// `inverse_round_roots[i]` has `2^i` entries.
    inverse_round_roots: Vec<Vec<Field<P>>>,
}

impl<P: FieldParams> EvaluationDomain<P> {
    /// Create a new evaluation domain for the given power-of-2 size.
    ///
    /// Computes the primitive root of unity, its inverse, domain element and
    /// inverse, and the coset generator. Lookup tables are **not** computed
    /// until [`compute_lookup_table`] is called.
    ///
    /// # Panics
    /// Panics if `domain_size` is not a power of 2.
    pub fn new(domain_size: usize) -> Self {
        assert!(
            domain_size.is_power_of_two(),
            "domain_size must be a power of 2, got {}",
            domain_size
        );

        let log2_size = domain_size.trailing_zeros() as usize;

        // For fields without high 2-adicity (e.g. Grumpkin Fq), there is no
        // meaningful root of unity; C++ falls back to Field::one().
        let root = if P::HAS_HIGH_2ADICITY {
            Field::<P>::get_root_of_unity(log2_size)
        } else {
            Field::<P>::one()
        };

        let root_inverse = root.invert();

        let domain = Field::<P>::from(domain_size as u64);
        let domain_inverse = domain.invert();

        // multiplicative_generator() == coset_generators[0], stored in Montgomery form.
        let generator = Field::<P>::from_raw([
            P::COSET_GENERATORS_0[0],
            P::COSET_GENERATORS_1[0],
            P::COSET_GENERATORS_2[0],
            P::COSET_GENERATORS_3[0],
        ]);
        let generator_inverse = generator.invert();

        let four_inverse = Field::<P>::from(4u64).invert();

        Self {
            size: domain_size,
            log2_size,
            root,
            root_inverse,
            domain,
            domain_inverse,
            generator,
            generator_inverse,
            four_inverse,
            round_roots: Vec::new(),
            inverse_round_roots: Vec::new(),
        }
    }

    /// Fill the twiddle-factor lookup tables for every butterfly round.
    ///
    /// For round `i` (0-indexed), the butterfly group size is `m = 2^(i+1)` and
    /// the table holds `m / 2 = 2^i` successive powers of the `m`-th root of
    /// unity (forward and inverse).
    pub fn compute_lookup_table(&mut self) {
        self.round_roots.resize(self.log2_size, Vec::new());
        self.inverse_round_roots.resize(self.log2_size, Vec::new());

        for i in 0..self.log2_size {
            let table_size = 1usize << i; // m/2 = 2^i

            // Derive the root of order m = 2^(i+1) by squaring the domain root
            // (log2_size - i - 1) times.
            let mut round_root = self.root;
            for _ in 0..(self.log2_size - i - 1) {
                round_root = round_root.sqr();
            }

            let mut inv_round_root = self.root_inverse;
            for _ in 0..(self.log2_size - i - 1) {
                inv_round_root = inv_round_root.sqr();
            }

            self.round_roots[i] = compute_lookup_table_single(round_root, table_size);
            self.inverse_round_roots[i] =
                compute_lookup_table_single(inv_round_root, table_size);
        }
    }

    /// Access the forward-NTT twiddle tables.
    pub fn get_round_roots(&self) -> &Vec<Vec<Field<P>>> {
        &self.round_roots
    }

    /// Access the inverse-NTT twiddle tables.
    pub fn get_inverse_round_roots(&self) -> &Vec<Vec<Field<P>>> {
        &self.inverse_round_roots
    }
}

/// Build a single lookup table: `[1, r, r^2, r^3, ..., r^(size-1)]`.
fn compute_lookup_table_single<P: FieldParams>(
    input_root: Field<P>,
    size: usize,
) -> Vec<Field<P>> {
    let mut table = Vec::with_capacity(size);
    if size == 0 {
        return table;
    }
    table.push(Field::<P>::one());
    for i in 1..size {
        table.push(table[i - 1] * input_root);
    }
    table
}
