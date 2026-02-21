//! Port of `relation_parameters.hpp` — parameters used by grand product relations.
//!
//! Contains eta, beta, gamma and other challenge values needed by
//! permutation, lookup, and memory relations.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

/// Number of binary limbs in Goblin Translator.
pub const NUM_BINARY_LIMBS_IN_GOBLIN_TRANSLATOR: usize = 4;

/// Number of native limbs in Goblin Translator.
pub const NUM_NATIVE_LIMBS_IN_GOBLIN_TRANSLATOR: usize = 1;

/// Number of challenge powers in Goblin Translator.
pub const NUM_CHALLENGE_POWERS_IN_GOBLIN_TRANSLATOR: usize = 4;

/// Total limbs = binary + native.
pub const NUM_TOTAL_LIMBS: usize =
    NUM_BINARY_LIMBS_IN_GOBLIN_TRANSLATOR + NUM_NATIVE_LIMBS_IN_GOBLIN_TRANSLATOR;

/// Number of fields to fold.
pub const NUM_TO_FOLD: usize = 6;

/// Container for parameters used by the grand product (permutation, lookup) Honk relations.
///
/// Port of C++ `RelationParameters<T>`.
#[derive(Clone, Debug)]
pub struct RelationParameters<T: Clone> {
    /// Lookup + Aux Memory
    pub eta: T,
    /// Lookup + Aux Memory
    pub eta_two: T,
    /// Lookup + Aux Memory
    pub eta_three: T,
    /// Permutation + Lookup
    pub beta: T,
    /// Permutation + Lookup
    pub gamma: T,
    /// Permutation
    pub public_input_delta: T,
    /// beta^2
    pub beta_sqr: T,
    /// beta^3
    pub beta_cube: T,
    /// Used in ECCVM set membership gadget.
    pub eccvm_set_permutation_delta: T,
    /// Translator accumulated result (4 limbs).
    pub accumulated_result: [T; NUM_BINARY_LIMBS_IN_GOBLIN_TRANSLATOR],
    /// Translator evaluation input x (5 limbs).
    pub evaluation_input_x: [T; NUM_TOTAL_LIMBS],
    /// Translator batching challenge v (4 × 5 matrix).
    pub batching_challenge_v: [[T; NUM_TOTAL_LIMBS]; NUM_CHALLENGE_POWERS_IN_GOBLIN_TRANSLATOR],
}

impl<P: FieldParams> RelationParameters<Field<P>> {
    /// All parameters initialized to zero.
    pub fn new() -> Self {
        Self {
            eta: Field::zero(),
            eta_two: Field::zero(),
            eta_three: Field::zero(),
            beta: Field::zero(),
            gamma: Field::zero(),
            public_input_delta: Field::zero(),
            beta_sqr: Field::zero(),
            beta_cube: Field::zero(),
            eccvm_set_permutation_delta: Field::zero(),
            accumulated_result: [Field::zero(); NUM_BINARY_LIMBS_IN_GOBLIN_TRANSLATOR],
            evaluation_input_x: [Field::zero(); NUM_TOTAL_LIMBS],
            batching_challenge_v: [[Field::zero(); NUM_TOTAL_LIMBS];
                NUM_CHALLENGE_POWERS_IN_GOBLIN_TRANSLATOR],
        }
    }

    /// Returns a reference array of the 6 fields to fold:
    /// [eta, eta_two, eta_three, beta, gamma, public_input_delta]
    pub fn get_to_fold(&self) -> [&Field<P>; NUM_TO_FOLD] {
        [
            &self.eta,
            &self.eta_two,
            &self.eta_three,
            &self.beta,
            &self.gamma,
            &self.public_input_delta,
        ]
    }

    /// Generate random parameters (for testing).
    ///
    /// Port of C++ `RelationParameters::get_random()`.
    pub fn get_random() -> Self {
        let beta = Field::<P>::random_element();
        let beta_sqr = beta * beta;
        let beta_cube = beta_sqr * beta;
        let gamma = Field::<P>::random_element();

        // eccvm_set_permutation_delta = γ·(γ + β²)·(γ + 2β²)·(γ + 3β²)
        let eccvm_set_permutation_delta = gamma
            * (gamma + beta_sqr)
            * (gamma + beta_sqr + beta_sqr)
            * (gamma + beta_sqr + beta_sqr + beta_sqr);

        let accumulated_result = [
            Field::random_element(),
            Field::random_element(),
            Field::random_element(),
            Field::random_element(),
        ];

        let evaluation_input_x = [
            Field::random_element(),
            Field::random_element(),
            Field::random_element(),
            Field::random_element(),
            Field::random_element(),
        ];

        let batching_challenge_v = [
            [
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
            ],
            [
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
            ],
            [
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
            ],
            [
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
                Field::random_element(),
            ],
        ];

        Self {
            eta: Field::random_element(),
            eta_two: Field::random_element(),
            eta_three: Field::random_element(),
            beta,
            gamma,
            public_input_delta: Field::random_element(),
            beta_sqr,
            beta_cube,
            eccvm_set_permutation_delta,
            accumulated_result,
            evaluation_input_x,
            batching_challenge_v,
        }
    }
}

impl<P: FieldParams> Default for RelationParameters<Field<P>> {
    fn default() -> Self {
        Self::new()
    }
}
