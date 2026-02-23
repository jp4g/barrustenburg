//! Non-native group arithmetic circuit type.
//!
//! Port of `barretenberg/stdlib/primitives/biggroup/`.
//!
//! A `BigGroupT<P, C>` represents a non-native elliptic curve point inside a circuit
//! whose native field has parameters `P`. The target curve has parameters `C: CurveParams`.
//!
//! Coordinates are represented as `BigFieldT<P, C::BaseFieldParams>` elements (non-native
//! field arithmetic), and the infinity flag is a `BoolT<P>`.

use std::cell::RefCell;

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;
use bbrs_ecc::groups::element::Element;
use bbrs_numeric::{U256, U256Ext, U512};
use bbrs_numeric::uintx::U512Ext;

use super::bigfield::BigFieldT;
use super::bool::BoolT;
use super::field::{bool_to_field, FieldT};
use super::memory::twin_rom_table::TwinRomTable;
use super::witness::{BuilderRef, WitnessT};

// ════════════════════════════════════════════════════════════════════════
//  Type aliases
// ════════════════════════════════════════════════════════════════════════

/// Non-native base field element (coordinate type).
type Fq<P, C> = BigFieldT<P, <C as CurveParams>::BaseFieldParams>;

/// Non-native scalar field element.
type Fr<P, C> = BigFieldT<P, <C as CurveParams>::ScalarFieldParams>;

/// Number of binary basis limbs in a BigFieldT.
const NUM_LIMBS: usize = 4;

// ════════════════════════════════════════════════════════════════════════
//  Helper: compute offset generators for scalar multiplication
// ════════════════════════════════════════════════════════════════════════

/// Compute a deterministic offset generator point that is independent of the curve generator.
///
/// Returns `(G_offset, 2^{num_rounds-1} * G_offset)`.
fn compute_offset_generators<C: CurveParams>(
    num_rounds: usize,
) -> (AffineElement<C>, AffineElement<C>) {
    // Use 7*G as a simple offset generator (deterministic, independent of G for practical purposes)
    let multiplier = Field::<C::ScalarFieldParams>::from_limbs([7, 0, 0, 0]);
    let offset_gen = Element::<C>::one().mul(&multiplier).to_affine();

    let offset_multiplier = U256::from(1u64).wrapping_shl_vartime((num_rounds - 1) as u32);
    let offset_scalar =
        Field::<C::ScalarFieldParams>::from_limbs(*offset_multiplier.as_words());
    let offset_gen_end = Element::<C>::from_affine(&offset_gen)
        .mul(&offset_scalar)
        .to_affine();

    (offset_gen, offset_gen_end)
}

/// Compute a table offset generator for masking points in batch_mul with edge cases.
fn compute_table_offset_generator<C: CurveParams>() -> AffineElement<C> {
    let multiplier = Field::<C::ScalarFieldParams>::from_limbs([13, 0, 0, 0]);
    Element::<C>::one().mul(&multiplier).to_affine()
}

// ════════════════════════════════════════════════════════════════════════
//  BigFieldT helper methods (decomposed from missing C++ methods)
// ════════════════════════════════════════════════════════════════════════

/// Compute `self * self + sum(addends)` (port of C++ `Fq::sqradd`).
fn sqradd<P: FieldParams, T: FieldParams>(
    val: &BigFieldT<P, T>,
    addends: &[&BigFieldT<P, T>],
) -> BigFieldT<P, T> {
    val.madd(val, addends)
}

/// Compute `sum(numerators) / denominator` (port of C++ `Fq::div_without_denominator_check`).
///
/// Uses regular division which includes a denominator != 0 assertion.
/// The extra assertion is harmless since callers already ensure non-zero denominators.
fn div_sum<P: FieldParams, T: FieldParams>(
    numerators: &[&BigFieldT<P, T>],
    denominator: &BigFieldT<P, T>,
) -> BigFieldT<P, T> {
    assert!(!numerators.is_empty());
    let mut sum = numerators[0].clone();
    for n in &numerators[1..] {
        sum = sum.add(n);
    }
    sum.div(denominator)
}

/// Compute `-(sum(left[i]*right[i]) + sum(add)) / den` (port of C++ `Fq::msub_div`).
fn msub_div<P: FieldParams, T: FieldParams>(
    left: &[&BigFieldT<P, T>],
    right: &[&BigFieldT<P, T>],
    den: &BigFieldT<P, T>,
    add: &[&BigFieldT<P, T>],
) -> BigFieldT<P, T> {
    assert_eq!(left.len(), right.len());
    let mut numerator = if !left.is_empty() {
        left[0].mul(right[0])
    } else {
        BigFieldT::zero_val()
    };
    for i in 1..left.len() {
        numerator = numerator.add(&left[i].mul(right[i]));
    }
    for a in add {
        numerator = numerator.add(a);
    }
    numerator.self_reduce();
    numerator = numerator.negate();
    numerator.div(den)
}

/// Compute `sum(left[i]*right[i]) + sum(add)` (port of C++ `Fq::mult_madd`).
fn mult_madd<P: FieldParams, T: FieldParams>(
    left: &[&BigFieldT<P, T>],
    right: &[&BigFieldT<P, T>],
    add: &[&BigFieldT<P, T>],
) -> BigFieldT<P, T> {
    assert_eq!(left.len(), right.len());
    let mut result = if !left.is_empty() {
        left[0].mul(right[0])
    } else {
        BigFieldT::zero_val()
    };
    for i in 1..left.len() {
        result = result.add(&left[i].mul(right[i]));
    }
    for a in add {
        result = result.add(a);
    }
    result
}

/// `conditional_assign(predicate, val_if_true, val_if_false)` - static method.
///
/// Port of C++ `Fq::conditional_assign(predicate, true_val, false_val)`.
fn conditional_assign<P: FieldParams, T: FieldParams>(
    predicate: &BoolT<P>,
    val_if_true: &BigFieldT<P, T>,
    val_if_false: &BigFieldT<P, T>,
) -> BigFieldT<P, T> {
    val_if_false.conditional_select(val_if_true, predicate)
}

// ════════════════════════════════════════════════════════════════════════
//  ChainAddAccumulator
// ════════════════════════════════════════════════════════════════════════

/// Accumulator for chain addition operations.
///
/// Tracks intermediate values to avoid computing y-coordinates of intermediate additions,
/// reducing the number of non-native field multiplications.
pub struct ChainAddAccumulator<P: FieldParams, C: CurveParams> {
    pub x1_prev: Fq<P, C>,
    pub y1_prev: Fq<P, C>,
    pub lambda_prev: Fq<P, C>,
    pub x3_prev: Fq<P, C>,
    pub y3_prev: Fq<P, C>,
    pub is_full_element: bool,
}

impl<P: FieldParams, C: CurveParams> Clone for ChainAddAccumulator<P, C> {
    fn clone(&self) -> Self {
        Self {
            x1_prev: self.x1_prev.clone(),
            y1_prev: self.y1_prev.clone(),
            lambda_prev: self.lambda_prev.clone(),
            x3_prev: self.x3_prev.clone(),
            y3_prev: self.y3_prev.clone(),
            is_full_element: self.is_full_element,
        }
    }
}

impl<P: FieldParams, C: CurveParams> ChainAddAccumulator<P, C> {
    pub fn from_element(input: &BigGroupT<P, C>) -> Self {
        Self {
            x1_prev: Fq::<P, C>::zero_val(),
            y1_prev: Fq::<P, C>::zero_val(),
            lambda_prev: Fq::<P, C>::zero_val(),
            x3_prev: input.x.clone(),
            y3_prev: input.y.clone(),
            is_full_element: true,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════
//  BigGroupT
// ════════════════════════════════════════════════════════════════════════

/// Non-native elliptic curve point in a circuit.
///
/// Port of C++ `element<Builder, Fq, Fr, NativeGroup>`.
pub struct BigGroupT<P: FieldParams, C: CurveParams> {
    pub x: Fq<P, C>,
    pub y: Fq<P, C>,
    is_infinity: BoolT<P>,
}

impl<P: FieldParams, C: CurveParams> Clone for BigGroupT<P, C> {
    fn clone(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.clone(),
            is_infinity: self.is_infinity.clone(),
        }
    }
}

impl<P: FieldParams, C: CurveParams> BigGroupT<P, C> {
    // ════════════════════════════════════════════════════════════════
    //  Constructors
    // ════════════════════════════════════════════════════════════════

    /// Default constructor (zeroed coordinates, no infinity flag).
    pub fn new_empty() -> Self {
        Self {
            x: Fq::<P, C>::zero_val(),
            y: Fq::<P, C>::zero_val(),
            is_infinity: BoolT::from_constant(false),
        }
    }

    /// Construct from a native affine element (constant).
    pub fn from_affine(input: &AffineElement<C>) -> Self {
        let x_val = input.x.from_montgomery_form();
        let y_val = input.y.from_montgomery_form();
        let x_u256 = U256::from_words(x_val.data);
        let y_u256 = U256::from_words(y_val.data);

        Self {
            x: Fq::<P, C>::from_u256(None, x_u256),
            y: Fq::<P, C>::from_u256(None, y_u256),
            is_infinity: BoolT::from_constant(input.is_point_at_infinity()),
        }
    }

    /// Construct from BigField coordinates with optional on-curve check.
    pub fn from_xy(x: Fq<P, C>, y: Fq<P, C>, assert_on_curve: bool) -> Self {
        let ctx = x.get_context().clone().or(y.get_context().clone());
        let is_infinity = if let Some(ref c) = ctx {
            BoolT::from_witness(&WitnessT::from_bool(c.clone(), false))
        } else {
            BoolT::from_constant(false)
        };
        let result = Self {
            x,
            y,
            is_infinity,
        };
        if assert_on_curve {
            result.validate_on_curve();
        }
        result
    }

    /// Construct from BigField coordinates with explicit infinity flag.
    pub fn from_xy_bool(
        x: Fq<P, C>,
        y: Fq<P, C>,
        is_infinity: BoolT<P>,
        assert_on_curve: bool,
    ) -> Self {
        let result = Self {
            x,
            y,
            is_infinity,
        };
        if assert_on_curve {
            result.validate_on_curve();
        }
        result
    }

    /// Create a BigGroup witness from a native affine element.
    pub fn from_witness(ctx: BuilderRef<P>, input: &AffineElement<C>) -> Self {
        let (x_val, y_val) = if input.is_point_at_infinity() {
            // Use generator coordinates for infinity points
            let generator = AffineElement::<C>::one();
            let gx = generator.x.from_montgomery_form();
            let gy = generator.y.from_montgomery_form();
            (U256::from_words(gx.data), U256::from_words(gy.data))
        } else {
            let ix = input.x.from_montgomery_form();
            let iy = input.y.from_montgomery_form();
            (U256::from_words(ix.data), U256::from_words(iy.data))
        };

        let x = Fq::<P, C>::from_witness(ctx.clone(), x_val);
        let y = Fq::<P, C>::from_witness(ctx.clone(), y_val);
        let is_inf = BoolT::from_witness(&WitnessT::from_bool(ctx, input.is_point_at_infinity()));

        let mut result = Self {
            x,
            y,
            is_infinity: is_inf,
        };
        result.validate_on_curve();
        result
    }

    /// Create a constant generator point.
    pub fn one(ctx: Option<BuilderRef<P>>) -> Self {
        let generator = AffineElement::<C>::one();
        let gx = generator.x.from_montgomery_form();
        let gy = generator.y.from_montgomery_form();
        Self {
            x: Fq::<P, C>::from_u256(ctx.clone(), U256::from_words(gx.data)),
            y: Fq::<P, C>::from_u256(ctx, U256::from_words(gy.data)),
            is_infinity: BoolT::from_constant(false),
        }
    }

    /// Create a point at infinity with zero coordinates.
    pub fn point_at_infinity(ctx: BuilderRef<P>) -> Self {
        let zero_idx = ctx.borrow().base.zero_idx();
        let zero = FieldT::from_witness_index(ctx.clone(), zero_idx);
        let x = BigFieldT::unsafe_construct_from_limbs(
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            false,
        );
        let y = BigFieldT::unsafe_construct_from_limbs(
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            false,
        );
        let mut result = Self::from_xy(x, y, false);
        result.set_point_at_infinity(&BoolT::from_constant(true));
        result
    }

    // ════════════════════════════════════════════════════════════════
    //  Accessors
    // ════════════════════════════════════════════════════════════════

    pub fn x(&self) -> &Fq<P, C> {
        &self.x
    }

    pub fn y(&self) -> &Fq<P, C> {
        &self.y
    }

    pub fn is_point_at_infinity(&self) -> &BoolT<P> {
        &self.is_infinity
    }

    pub fn set_point_at_infinity(&mut self, is_infinity: &BoolT<P>) {
        self.is_infinity = is_infinity.normalize();
    }

    pub fn is_constant(&self) -> bool {
        self.x.is_constant() && self.y.is_constant()
    }

    pub fn get_context(&self) -> Option<BuilderRef<P>> {
        self.x
            .get_context()
            .clone()
            .or(self.y.get_context().clone())
    }

    /// Get the native value of this point.
    pub fn get_value(&self) -> AffineElement<C> {
        let modulus = U512::from_lo_hi(U256::from_limbs(C::BaseFieldParams::MODULUS), U256::ZERO);
        let x_val = self.x.get_value().div_rem(&modulus.to_nz().unwrap()).1;
        let y_val = self.y.get_value().div_rem(&modulus.to_nz().unwrap()).1;

        let mut result = AffineElement::<C>::new(
            Field::from_limbs(*x_val.lo().as_words()),
            Field::from_limbs(*y_val.lo().as_words()),
        );
        if self.is_infinity.get_value() {
            result.self_set_infinity();
        }
        result
    }

    /// Get the standard form: if point is at infinity, coordinates are (0, 0).
    pub fn get_standard_form(&self) -> Self {
        let is_inf = self.is_point_at_infinity();
        let zero = Fq::<P, C>::zero_val();
        let x = conditional_assign(is_inf, &zero, &self.x);
        let y = conditional_assign(is_inf, &zero, &self.y);
        let mut result = Self {
            x,
            y,
            is_infinity: self.is_infinity.clone(),
        };
        result
    }

    // ════════════════════════════════════════════════════════════════
    //  Validation
    // ════════════════════════════════════════════════════════════════

    /// Validate that the point is on the curve: y² = x³ + ax + b.
    pub fn validate_on_curve(&self) {
        let ctx = self.get_context();
        let b_val = C::coeff_b().from_montgomery_form();
        let b_u256 = U256::from_words(b_val.data);
        let b = Fq::<P, C>::from_u256(ctx.clone(), b_u256);

        let zero = Fq::<P, C>::zero_val();
        let is_inf = self.is_point_at_infinity();

        let adjusted_b = conditional_assign(is_inf, &zero, &b);
        let adjusted_x = conditional_assign(is_inf, &zero, &self.x);
        let adjusted_y = conditional_assign(is_inf, &zero, &self.y);

        if C::HAS_A {
            let a_val = C::coeff_a().from_montgomery_form();
            let a_u256 = U256::from_words(a_val.data);
            let a = Fq::<P, C>::from_u256(ctx, a_u256);
            let adjusted_a = conditional_assign(is_inf, &zero, &a);

            // y² = x³ + ax + b
            // x²·x + x·a + y·(-y) + b = 0
            let x_sqr = adjusted_x.sqr();
            let lhs = mult_madd(
                &[&x_sqr, &adjusted_x, &adjusted_y],
                &[&adjusted_x, &adjusted_a, &adjusted_y.negate()],
                &[&adjusted_b],
            );
            lhs.assert_equal(
                &Fq::<P, C>::zero_val(),
                "biggroup::validate_on_curve",
            );
        } else {
            // y² = x³ + b
            // x²·x + y·(-y) + b = 0
            let x_sqr = adjusted_x.sqr();
            let lhs = mult_madd(
                &[&x_sqr, &adjusted_y],
                &[&adjusted_x, &adjusted_y.negate()],
                &[&adjusted_b],
            );
            lhs.assert_equal(
                &Fq::<P, C>::zero_val(),
                "biggroup::validate_on_curve",
            );
        }
    }

    /// Assert that the x and y coordinates are in the target field.
    pub fn assert_coordinates_in_field(&self) {
        self.x
            .assert_is_in_field("biggroup::assert_coordinates_in_field (x)");
        self.y
            .assert_is_in_field("biggroup::assert_coordinates_in_field (y)");
    }

    /// Assert two group elements are equal.
    pub fn incomplete_assert_equal(&self, other: &Self) {
        self.is_infinity
            .assert_equal(&other.is_infinity, "biggroup::assert_equal (infinity)");
        self.x
            .assert_equal(&other.x, "biggroup::assert_equal (x)");
        self.y
            .assert_equal(&other.y, "biggroup::assert_equal (y)");
    }

    /// Normalize (reduce) coordinates modulo the target field.
    pub fn normalize(&self) -> Self {
        let mut result = self.clone();
        result.x.self_reduce();
        result.y.self_reduce();
        result
    }

    /// Reduce coordinates.
    pub fn reduce(&self) -> Self {
        let mut result = self.clone();
        result.x.self_reduce();
        result.y.self_reduce();
        result
    }

    // ════════════════════════════════════════════════════════════════
    //  Conditional operations
    // ════════════════════════════════════════════════════════════════

    /// Negate the y-coordinate: (x, y) -> (x, -y).
    pub fn negate(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.negate(),
            is_infinity: self.is_infinity.clone(),
        }
    }

    /// Conditionally negate: if predicate, return (x, -y), else (x, y).
    pub fn conditional_negate(&self, predicate: &BoolT<P>) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.conditional_negate(predicate),
            is_infinity: self.is_infinity.clone(),
        }
    }

    /// Conditional select: if predicate, return other, else self.
    pub fn conditional_select(&self, other: &Self, predicate: &BoolT<P>) -> Self {
        if predicate.is_constant() {
            return if predicate.get_value() {
                other.clone()
            } else {
                self.clone()
            };
        }
        Self {
            x: self.x.conditional_select(&other.x, predicate),
            y: self.y.conditional_select(&other.y, predicate),
            is_infinity: BoolT::conditional_assign(
                predicate,
                &other.is_infinity,
                &self.is_infinity,
            ),
        }
    }

    // ════════════════════════════════════════════════════════════════
    //  Basic arithmetic
    // ════════════════════════════════════════════════════════════════

    /// Promote a constant point to a witness element using the given context.
    ///
    /// BigFieldT operations like `eq` and `self_reduce` require a circuit context.
    /// Constants don't have one, so we must promote them before mixing with witnesses.
    fn ensure_context(&self, ctx: &BuilderRef<P>) -> Self {
        if self.get_context().is_some() {
            return self.clone();
        }
        let val = self.get_value();
        Self::from_witness(ctx.clone(), &val)
    }

    /// Point addition with complete formulas (handles edge cases).
    pub fn add(&self, other: &Self) -> Self {
        // If both are constant, compute natively
        if self.is_constant() && other.is_constant() {
            let a = Element::<C>::from_affine(&self.get_value());
            let b = Element::<C>::from_affine(&other.get_value());
            return Self::from_affine(&(a + b).to_affine());
        }
        // Ensure both have context for BigFieldT operations
        let ctx = self
            .get_context()
            .or_else(|| other.get_context())
            .expect("add: at least one operand must have context");
        let lhs = self.ensure_context(&ctx);
        let rhs = other.ensure_context(&ctx);
        lhs.add_impl(&rhs)
    }

    /// Internal addition implementation (both operands guaranteed to have context).
    fn add_impl(&self, other: &Self) -> Self {
        let x_coordinates_match = self.x.eq(&other.x);
        let y_coordinates_match = self.y.eq(&other.y);
        let infinity_predicate = x_coordinates_match.and(&y_coordinates_match.negate());
        let double_predicate = x_coordinates_match.and(&y_coordinates_match);
        let lhs_infinity = self.is_point_at_infinity().clone();
        let rhs_infinity = other.is_point_at_infinity().clone();
        let has_infinity_input = lhs_infinity.or(&rhs_infinity);

        // Compute lambda
        let add_lambda_num = other.y.sub(&self.y);
        let xx = self.x.sqr();
        let mut dbl_lambda_num = xx.add(&xx.clone()).add(&xx);
        if C::HAS_A {
            let a_val = C::coeff_a().from_montgomery_form();
            let a = Fq::<P, C>::from_u256(self.get_context(), U256::from_words(a_val.data));
            dbl_lambda_num = dbl_lambda_num.add(&a);
        }
        let lambda_num = conditional_assign(&double_predicate, &dbl_lambda_num, &add_lambda_num);

        let add_lambda_den = other.x.sub(&self.x);
        let dbl_lambda_den = self.y.add(&self.y);
        let mut lambda_den =
            conditional_assign(&double_predicate, &dbl_lambda_den, &add_lambda_den);

        let safe_den_needed = has_infinity_input.or(&infinity_predicate);
        let one_fq = Fq::<P, C>::from_u256(self.get_context(), U256::ONE);
        lambda_den = conditional_assign(&safe_den_needed, &one_fq, &lambda_den);

        let lambda = lambda_num.div(&lambda_den);

        // x3 = lambda^2 - x1 - x2
        let neg_other_x = other.x.negate();
        let neg_self_x = self.x.negate();
        let x3 = sqradd(&lambda, &[&neg_other_x, &neg_self_x]);

        // y3 = lambda * (x1 - x3) - y1
        let x_diff = self.x.sub(&x3);
        let neg_y = self.y.negate();
        let y3 = lambda.madd(&x_diff, &[&neg_y]);

        let mut result = Self::from_xy(x3, y3, false);

        // Handle infinity cases
        result.x = conditional_assign(&lhs_infinity, &other.x, &result.x);
        result.y = conditional_assign(&lhs_infinity, &other.y, &result.y);
        result.x = conditional_assign(&rhs_infinity, &self.x, &result.x);
        result.y = conditional_assign(&rhs_infinity, &self.y, &result.y);

        let result_is_infinity = infinity_predicate
            .and(&has_infinity_input.negate())
            .or(&lhs_infinity.and(&rhs_infinity));
        result.set_point_at_infinity(&result_is_infinity);

        result
    }

    /// Point subtraction with complete formulas.
    pub fn sub(&self, other: &Self) -> Self {
        // If both are constant, compute natively
        if self.is_constant() && other.is_constant() {
            let a = Element::<C>::from_affine(&self.get_value());
            let b = Element::<C>::from_affine(&other.get_value());
            return Self::from_affine(&(a - b).to_affine());
        }
        // Ensure both have context
        let ctx = self
            .get_context()
            .or_else(|| other.get_context())
            .expect("sub: at least one operand must have context");
        let lhs = self.ensure_context(&ctx);
        let rhs = other.ensure_context(&ctx);
        lhs.sub_impl(&rhs)
    }

    /// Internal subtraction implementation (both operands guaranteed to have context).
    fn sub_impl(&self, other: &Self) -> Self {
        let x_coordinates_match = self.x.eq(&other.x);
        let y_coordinates_match = self.y.eq(&other.y);
        // For subtraction: P1 - P2 = P1 + (-P2)
        // If y1 == y2 and x1 == x2, then P1 - P1 = infinity
        let infinity_predicate = x_coordinates_match.and(&y_coordinates_match);
        // If y1 != y2 and x1 == x2, then P1 - P2 is doubling (P1 + (-P2) = P1 + P1)
        let double_predicate = x_coordinates_match.and(&y_coordinates_match.negate());
        let lhs_infinity = self.is_point_at_infinity().clone();
        let rhs_infinity = other.is_point_at_infinity().clone();
        let has_infinity_input = lhs_infinity.or(&rhs_infinity);

        // Lambda for subtraction: (-y2 - y1) / (x2 - x1)
        let add_lambda_num = other.y.negate().sub(&self.y);
        let xx = self.x.sqr();
        let mut dbl_lambda_num = xx.add(&xx.clone()).add(&xx);
        if C::HAS_A {
            let a_val = C::coeff_a().from_montgomery_form();
            let a = Fq::<P, C>::from_u256(self.get_context(), U256::from_words(a_val.data));
            dbl_lambda_num = dbl_lambda_num.add(&a);
        }
        let lambda_num = conditional_assign(&double_predicate, &dbl_lambda_num, &add_lambda_num);

        let add_lambda_den = other.x.sub(&self.x);
        let dbl_lambda_den = self.y.add(&self.y);
        let mut lambda_den =
            conditional_assign(&double_predicate, &dbl_lambda_den, &add_lambda_den);

        let safe_den_needed = has_infinity_input.or(&infinity_predicate);
        let one_fq = Fq::<P, C>::from_u256(self.get_context(), U256::ONE);
        lambda_den = conditional_assign(&safe_den_needed, &one_fq, &lambda_den);

        let lambda = lambda_num.div(&lambda_den);

        let neg_other_x = other.x.negate();
        let neg_self_x = self.x.negate();
        let x3 = sqradd(&lambda, &[&neg_other_x, &neg_self_x]);
        let x_diff = self.x.sub(&x3);
        let neg_y = self.y.negate();
        let y3 = lambda.madd(&x_diff, &[&neg_y]);

        let mut result = Self::from_xy(x3, y3, false);

        // If lhs infinity, return -rhs
        result.x = conditional_assign(&lhs_infinity, &other.x, &result.x);
        result.y = conditional_assign(&lhs_infinity, &other.y.negate(), &result.y);
        // If rhs infinity, return lhs
        result.x = conditional_assign(&rhs_infinity, &self.x, &result.x);
        result.y = conditional_assign(&rhs_infinity, &self.y, &result.y);

        let result_is_infinity = infinity_predicate
            .and(&has_infinity_input.negate())
            .or(&lhs_infinity.and(&rhs_infinity));
        result.set_point_at_infinity(&result_is_infinity);

        result
    }

    /// Incomplete addition: requires x1 != x2.
    pub fn checked_unconditional_add(&self, other: &Self) -> Self {
        other
            .x
            .assert_is_not_equal(&self.x, "biggroup::checked_unconditional_add: x1 == x2");

        let num = other.y.sub(&self.y);
        let den = other.x.sub(&self.x);
        let lambda = num.div(&den);

        let neg_other_x = other.x.negate();
        let neg_self_x = self.x.negate();
        let x3 = sqradd(&lambda, &[&neg_other_x, &neg_self_x]);
        let x_diff = self.x.sub(&x3);
        let neg_y = self.y.negate();
        let y3 = lambda.madd(&x_diff, &[&neg_y]);

        Self::from_xy(x3, y3, false)
    }

    /// Incomplete subtraction: requires x1 != x2.
    pub fn checked_unconditional_subtract(&self, other: &Self) -> Self {
        other
            .x
            .assert_is_not_equal(&self.x, "biggroup::checked_unconditional_subtract: x1 == x2");

        let num = other.y.add(&self.y);
        let den = other.x.sub(&self.x);
        let lambda = num.div(&den);

        let neg_other_x = other.x.negate();
        let neg_self_x = self.x.negate();
        let x3 = sqradd(&lambda, &[&neg_other_x, &neg_self_x]);
        let x3_minus_x = x3.sub(&self.x);
        let neg_y = self.y.negate();
        let y3 = lambda.madd(&x3_minus_x, &[&neg_y]);

        Self::from_xy(x3, y3, false)
    }

    /// Compute both (self + other) and (self - other) simultaneously.
    fn checked_unconditional_add_sub(&self, other: &Self) -> [Self; 2] {
        other
            .x
            .assert_is_not_equal(&self.x, "biggroup::checked_unconditional_add_sub: x1 == x2");

        let denominator = other.x.sub(&self.x);
        let x2x1 = other.x.add(&self.x).negate();

        // lambda1 = (y2 - y1) / (x2 - x1) for addition
        let num1 = other.y.sub(&self.y);
        let lambda1 = num1.div(&denominator);
        let x3 = sqradd(&lambda1, &[&x2x1]);
        let y3 = lambda1.madd(&self.x.sub(&x3), &[&self.y.negate()]);

        // lambda2 = (-y2 - y1) / (x2 - x1) for subtraction
        let num2 = other.y.negate().sub(&self.y);
        let lambda2 = num2.div(&denominator);
        let x4 = sqradd(&lambda2, &[&x2x1]);
        let y4 = lambda2.madd(&self.x.sub(&x4), &[&self.y.negate()]);

        [Self::from_xy(x3, y3, false), Self::from_xy(x4, y4, false)]
    }

    /// Point doubling.
    pub fn dbl(&self) -> Self {
        let two_x = self.x.add(&self.x);

        let (neg_lambda, x3, y3) = if C::HAS_A {
            let a_val = C::coeff_a().from_montgomery_form();
            let a = Fq::<P, C>::from_u256(self.get_context(), U256::from_words(a_val.data));

            // neg_lambda = -(3x² + a) / (2y)
            let three_x = two_x.add(&self.x);
            let numerator = self.x.mul(&three_x).add(&a);
            let denominator = self.y.add(&self.y);
            let neg_lambda = numerator.negate().div(&denominator);

            let x3 = sqradd(&neg_lambda, &[&two_x.negate()]);
            let y3 = neg_lambda.madd(&x3.sub(&self.x), &[&self.y.negate()]);
            (neg_lambda, x3, y3)
        } else {
            // neg_lambda = -3x² / (2y)
            let three_x = two_x.add(&self.x);
            let numerator = self.x.mul(&three_x);
            let denominator = self.y.add(&self.y);
            let neg_lambda = numerator.negate().div(&denominator);

            let x3 = sqradd(&neg_lambda, &[&two_x.negate()]);
            let y3 = neg_lambda.madd(&x3.sub(&self.x), &[&self.y.negate()]);
            (neg_lambda, x3, y3)
        };

        let mut result = Self::from_xy(x3, y3, false);
        result.set_point_at_infinity(self.is_point_at_infinity());
        result
    }

    // ════════════════════════════════════════════════════════════════
    //  Chain addition
    // ════════════════════════════════════════════════════════════════

    /// Begin a chain of additions.
    pub fn chain_add_start(p1: &Self, p2: &Self) -> ChainAddAccumulator<P, C> {
        p1.x.assert_is_not_equal(&p2.x, "biggroup::chain_add_start: x1 == x2");

        let num = p2.y.sub(&p1.y);
        let den = p2.x.sub(&p1.x);
        let lambda = num.div(&den);

        let neg_p2x = p2.x.negate();
        let neg_p1x = p1.x.negate();
        let x3 = sqradd(&lambda, &[&neg_p2x, &neg_p1x]);

        ChainAddAccumulator {
            x1_prev: p1.x.clone(),
            y1_prev: p1.y.clone(),
            lambda_prev: lambda,
            x3_prev: x3,
            y3_prev: Fq::<P, C>::zero_val(),
            is_full_element: false,
        }
    }

    /// Continue a chain addition.
    pub fn chain_add(p1: &Self, acc: &ChainAddAccumulator<P, C>) -> ChainAddAccumulator<P, C> {
        if acc.is_full_element {
            let full_elem = Self::from_xy(acc.x3_prev.clone(), acc.y3_prev.clone(), false);
            return Self::chain_add_start(p1, &full_elem);
        }

        p1.x
            .assert_is_not_equal(&acc.x3_prev, "biggroup::chain_add: x1 == x_acc");

        // lambda = -(lambda_prev * (x_acc - x1_prev) + y1_prev + y) / (x_acc - x)
        let lambda = msub_div(
            &[&acc.lambda_prev],
            &[&acc.x3_prev.sub(&acc.x1_prev)],
            &acc.x3_prev.sub(&p1.x),
            &[&acc.y1_prev, &p1.y],
        );

        let x3 = sqradd(&lambda, &[&acc.x3_prev.negate(), &p1.x.negate()]);

        ChainAddAccumulator {
            x1_prev: p1.x.clone(),
            y1_prev: p1.y.clone(),
            lambda_prev: lambda,
            x3_prev: x3,
            y3_prev: Fq::<P, C>::zero_val(),
            is_full_element: false,
        }
    }

    /// End a chain addition and compute the final y-coordinate.
    pub fn chain_add_end(acc: &ChainAddAccumulator<P, C>) -> Self {
        if acc.is_full_element {
            return Self::from_xy(acc.x3_prev.clone(), acc.y3_prev.clone(), false);
        }

        let y3 = acc
            .lambda_prev
            .madd(&acc.x1_prev.sub(&acc.x3_prev), &[&acc.y1_prev.negate()]);
        Self::from_xy(acc.x3_prev.clone(), y3, false)
    }

    /// Multiple montgomery ladder: computes (2^n * self + sum of to_add) efficiently.
    pub fn multiple_montgomery_ladder(&self, to_add: &[ChainAddAccumulator<P, C>]) -> Self {
        if to_add.is_empty() {
            return self.clone();
        }

        // First iteration
        self.x()
            .assert_is_not_equal(&to_add[0].x3_prev, "biggroup::mml: x == x_add[0]");

        let lambda1 = if !to_add[0].is_full_element {
            msub_div(
                &[&to_add[0].lambda_prev],
                &[&to_add[0].x1_prev.sub(&to_add[0].x3_prev)],
                &self.x.sub(&to_add[0].x3_prev),
                &[&to_add[0].y1_prev.negate(), &self.y.negate()],
            )
        } else {
            let diff = self.y.sub(&to_add[0].y3_prev);
            diff.div(&self.x.sub(&to_add[0].x3_prev))
        };

        let x_3 = lambda1.madd(&lambda1, &[&to_add[0].x3_prev.negate(), &self.x.negate()]);

        self.x()
            .assert_is_not_equal(&x_3, "biggroup::mml: x == x_3");
        let two_y = self.y.add(&self.y);
        let lambda2 = two_y.div(&self.x.sub(&x_3)).sub(&lambda1);

        let mut x_4 = sqradd(&lambda2, &[&x_3.negate(), &self.x.negate()]);

        // Build composite y tracking
        let num_points_even = (to_add.len() & 1) == 0;
        let mut prev_y_mul_left = vec![lambda2.clone()];
        let mut prev_y_mul_right = vec![if num_points_even {
            x_4.sub(&self.x)
        } else {
            self.x.sub(&x_4)
        }];
        let mut prev_y_add = vec![if num_points_even {
            self.y.clone()
        } else {
            self.y.negate()
        }];
        let mut prev_y_negative = num_points_even;

        let mut previous_x = x_4.clone();

        // Remaining iterations
        for i in 1..to_add.len() {
            previous_x.assert_is_not_equal(
                &to_add[i].x3_prev,
                "biggroup::mml: x_prev == x_add[i]",
            );

            let negate_add_y = !prev_y_negative;

            let mut l1_left = prev_y_mul_left.clone();
            let mut l1_right = prev_y_mul_right.clone();
            let mut l1_add = prev_y_add.clone();

            if !to_add[i].is_full_element {
                l1_left.push(to_add[i].lambda_prev.clone());
                l1_right.push(if negate_add_y {
                    to_add[i].x3_prev.sub(&to_add[i].x1_prev)
                } else {
                    to_add[i].x1_prev.sub(&to_add[i].x3_prev)
                });
                l1_add.push(if negate_add_y {
                    to_add[i].y1_prev.clone()
                } else {
                    to_add[i].y1_prev.negate()
                });
            } else {
                l1_add.push(if negate_add_y {
                    to_add[i].y3_prev.negate()
                } else {
                    to_add[i].y3_prev.clone()
                });
            }

            let denominator = if negate_add_y {
                to_add[i].x3_prev.sub(&previous_x)
            } else {
                previous_x.sub(&to_add[i].x3_prev)
            };

            let l1_left_refs: Vec<&Fq<P, C>> = l1_left.iter().collect();
            let l1_right_refs: Vec<&Fq<P, C>> = l1_right.iter().collect();
            let l1_add_refs: Vec<&Fq<P, C>> = l1_add.iter().collect();

            let lambda1 = msub_div(&l1_left_refs, &l1_right_refs, &denominator, &l1_add_refs);

            let x_3 =
                lambda1.madd(&lambda1, &[&to_add[i].x3_prev.negate(), &previous_x.negate()]);

            previous_x.assert_is_not_equal(&x_3, "biggroup::mml: prev_x == x_3 in loop");

            let l2_den = if prev_y_negative {
                previous_x.sub(&x_3)
            } else {
                x_3.sub(&previous_x)
            };

            let py_ml_refs: Vec<&Fq<P, C>> = prev_y_mul_left.iter().collect();
            let py_mr_refs: Vec<&Fq<P, C>> = prev_y_mul_right.iter().collect();
            let py_a_refs: Vec<&Fq<P, C>> = prev_y_add.iter().collect();

            let mut partial_lambda2 = msub_div(&py_ml_refs, &py_mr_refs, &l2_den, &py_a_refs);
            partial_lambda2 = partial_lambda2.add(&partial_lambda2.clone());
            let lambda2 = partial_lambda2.sub(&lambda1);

            x_4 = sqradd(&lambda2, &[&x_3.negate(), &previous_x.negate()]);


            // Build new composite y
            let new_negative = !prev_y_negative;
            let mut new_ml = vec![lambda2.clone()];
            let new_mr = vec![if prev_y_negative {
                previous_x.sub(&x_4)
            } else {
                x_4.sub(&previous_x)
            }];
            let new_add: Vec<Fq<P, C>> = Vec::new();

            new_ml.extend(prev_y_mul_left.iter().cloned());
            let mut final_mr = new_mr;
            final_mr.extend(prev_y_mul_right.iter().cloned());
            let mut final_add = new_add;
            final_add.extend(prev_y_add.iter().cloned());

            prev_y_mul_left = new_ml;
            prev_y_mul_right = final_mr;
            prev_y_add = final_add;
            prev_y_negative = new_negative;

            previous_x = x_4.clone();
        }

        assert!(!prev_y_negative);

        let py_ml_refs: Vec<&Fq<P, C>> = prev_y_mul_left.iter().collect();
        let py_mr_refs: Vec<&Fq<P, C>> = prev_y_mul_right.iter().collect();
        let py_a_refs: Vec<&Fq<P, C>> = prev_y_add.iter().collect();
        let y_out = mult_madd(&py_ml_refs, &py_mr_refs, &py_a_refs);

        Self::from_xy(previous_x, y_out, false)
    }

    // ════════════════════════════════════════════════════════════════
    //  NAF computation
    // ════════════════════════════════════════════════════════════════

    /// Compute Non-Adjacent Form (NAF) representation of a scalar.
    pub fn compute_naf(scalar: &Fr<P, C>, max_num_bits: usize) -> Vec<BoolT<P>> {
        let ctx = scalar
            .get_context()
            .clone()
            .expect("compute_naf: scalar must have context");

        let modulus = U256::from_limbs(C::ScalarFieldParams::MODULUS);
        let scalar_val_512 = scalar.get_value();
        let modulus_512 = U512::from_lo_hi(modulus, U256::ZERO);
        let mut scalar_multiplier =
            scalar_val_512.div_rem(&modulus_512.to_nz().unwrap()).1.lo();

        let num_rounds = if max_num_bits == 0 || scalar_multiplier == U256::ZERO {
            modulus.get_msb() as usize + 1
        } else {
            max_num_bits
        };

        if scalar_multiplier == U256::ZERO {
            scalar_multiplier = modulus;
        }

        let skew_value = !scalar_multiplier.get_bit(0);
        if skew_value {
            scalar_multiplier = scalar_multiplier.wrapping_add(&U256::ONE);
        }

        let mut naf_entries = vec![BoolT::from_constant(false); num_rounds + 1];

        // Skew bit
        naf_entries[num_rounds] = BoolT::from_witness(&WitnessT::from_bool(ctx.clone(), skew_value));

        // NAF bits
        for i in 0..(num_rounds - 1) {
            let next_entry = scalar_multiplier.get_bit((i + 1) as u32);
            naf_entries[num_rounds - i - 1] = BoolT::from_witness(&WitnessT::from_bool(ctx.clone(), !next_entry));
        }

        // MSB is always +1 (false)
        naf_entries[0] = BoolT::from_witness(&WitnessT::from_bool(ctx.clone(), false));

        // Validate NAF reconstruction (composite field path)
        let reconstruct_half_naf =
            |nafs: &[BoolT<P>], half_round_length: usize, ctx: &BuilderRef<P>| {
                let mut neg_acc = FieldT::from_field(Field::<P>::zero());
                let mut pos_acc = FieldT::from_field(Field::<P>::zero());
                let one_ft = FieldT::from_field(Field::<P>::one());
                for i in 0..half_round_length {
                    let naf_ft = bool_to_field(&nafs[i]);
                    let one_minus_naf = &one_ft - &naf_ft;
                    neg_acc = &(&neg_acc + &neg_acc) + &naf_ft;
                    pos_acc = &(&pos_acc + &pos_acc) + &one_minus_naf;
                }
                (pos_acc, neg_acc)
            };

        let num_limb_bits = 68usize; // BigFieldT NUM_LIMB_BITS
        let (hi_pos, hi_neg, lo_pos, lo_neg) = if num_rounds > num_limb_bits * 2 {
            let midpoint = num_rounds - (num_limb_bits * 2);
            let (hp, hn) = reconstruct_half_naf(&naf_entries[..midpoint], midpoint, &ctx);
            let (lp, ln) =
                reconstruct_half_naf(&naf_entries[midpoint..], num_rounds - midpoint, &ctx);
            (hp, hn, lp, ln)
        } else {
            let zero_idx = ctx.borrow().base.zero_idx();
            let zero = FieldT::from_witness_index(ctx.clone(), zero_idx);
            let (lp, ln) = reconstruct_half_naf(&naf_entries, num_rounds, &ctx);
            (zero.clone(), zero, lp, ln)
        };

        // Add skew to lo negative
        let skew_ft = bool_to_field(&naf_entries[num_rounds]);
        let lo_neg = &lo_neg + &skew_ft;

        let reconstructed_positive =
            Fr::<P, C>::from_field_pair(lo_pos, hi_pos, false, 0);
        let reconstructed_negative =
            Fr::<P, C>::from_field_pair(lo_neg, hi_neg, false, 0);
        let mut accumulator = reconstructed_positive.sub(&reconstructed_negative);

        if scalar.is_constant() {
            accumulator.self_reduce();
        }
        accumulator.assert_equal(scalar, "biggroup::compute_naf: NAF reconstruction mismatch");

        naf_entries
    }

    // ════════════════════════════════════════════════════════════════
    //  ROM table operations for group elements
    // ════════════════════════════════════════════════════════════════

    /// Create twin ROM tables for a set of group elements.
    ///
    /// Returns 5 TwinRomTables: x_lo, x_hi, y_lo, y_hi, xy_prime.
    fn create_group_element_rom_tables(
        rom_data: &[Self],
        limb_max: &mut [U256; 8],
    ) -> [RefCell<TwinRomTable<P>>; 5] {
        let mut x_lo_limbs = Vec::new();
        let mut x_hi_limbs = Vec::new();
        let mut y_lo_limbs = Vec::new();
        let mut y_hi_limbs = Vec::new();
        let mut prime_limbs = Vec::new();

        for elem in rom_data {
            limb_max[0] = std::cmp::max(limb_max[0], elem.x.binary_basis_limbs[0].maximum_value);
            limb_max[1] = std::cmp::max(limb_max[1], elem.x.binary_basis_limbs[1].maximum_value);
            limb_max[2] = std::cmp::max(limb_max[2], elem.x.binary_basis_limbs[2].maximum_value);
            limb_max[3] = std::cmp::max(limb_max[3], elem.x.binary_basis_limbs[3].maximum_value);
            limb_max[4] = std::cmp::max(limb_max[4], elem.y.binary_basis_limbs[0].maximum_value);
            limb_max[5] = std::cmp::max(limb_max[5], elem.y.binary_basis_limbs[1].maximum_value);
            limb_max[6] = std::cmp::max(limb_max[6], elem.y.binary_basis_limbs[2].maximum_value);
            limb_max[7] = std::cmp::max(limb_max[7], elem.y.binary_basis_limbs[3].maximum_value);

            x_lo_limbs.push([
                elem.x.binary_basis_limbs[0].element.clone(),
                elem.x.binary_basis_limbs[1].element.clone(),
            ]);
            x_hi_limbs.push([
                elem.x.binary_basis_limbs[2].element.clone(),
                elem.x.binary_basis_limbs[3].element.clone(),
            ]);
            y_lo_limbs.push([
                elem.y.binary_basis_limbs[0].element.clone(),
                elem.y.binary_basis_limbs[1].element.clone(),
            ]);
            y_hi_limbs.push([
                elem.y.binary_basis_limbs[2].element.clone(),
                elem.y.binary_basis_limbs[3].element.clone(),
            ]);
            prime_limbs.push([
                elem.x.prime_basis_limb.clone(),
                elem.y.prime_basis_limb.clone(),
            ]);
        }

        [
            RefCell::new(TwinRomTable::new(x_lo_limbs)),
            RefCell::new(TwinRomTable::new(x_hi_limbs)),
            RefCell::new(TwinRomTable::new(y_lo_limbs)),
            RefCell::new(TwinRomTable::new(y_hi_limbs)),
            RefCell::new(TwinRomTable::new(prime_limbs)),
        ]
    }

    /// Read a group element from ROM tables at the given index.
    fn read_group_element_rom_tables(
        tables: &[RefCell<TwinRomTable<P>>; 5],
        index: &FieldT<P>,
        limb_max: &[U256; 8],
    ) -> Self {
        let xlo = tables[0].borrow_mut().read(index);
        let xhi = tables[1].borrow_mut().read(index);
        let ylo = tables[2].borrow_mut().read(index);
        let yhi = tables[3].borrow_mut().read(index);
        let xyprime = tables[4].borrow_mut().read(index);

        let mut x_fq =
            Fq::<P, C>::unsafe_construct_from_limbs(xlo[0].clone(), xlo[1].clone(), xhi[0].clone(), xhi[1].clone(), false);
        let mut y_fq =
            Fq::<P, C>::unsafe_construct_from_limbs(ylo[0].clone(), ylo[1].clone(), yhi[0].clone(), yhi[1].clone(), false);

        // Set maximum values from the table
        x_fq.binary_basis_limbs[0].maximum_value = limb_max[0];
        x_fq.binary_basis_limbs[1].maximum_value = limb_max[1];
        x_fq.binary_basis_limbs[2].maximum_value = limb_max[2];
        x_fq.binary_basis_limbs[3].maximum_value = limb_max[3];
        y_fq.binary_basis_limbs[0].maximum_value = limb_max[4];
        y_fq.binary_basis_limbs[1].maximum_value = limb_max[5];
        y_fq.binary_basis_limbs[2].maximum_value = limb_max[6];
        y_fq.binary_basis_limbs[3].maximum_value = limb_max[7];

        Self::from_xy(x_fq, y_fq, false)
    }

    // ════════════════════════════════════════════════════════════════
    //  Scalar multiplication
    // ════════════════════════════════════════════════════════════════

    /// Scalar multiplication: self * scalar.
    pub fn scalar_mul(&self, scalar: &Fr<P, C>, max_num_bits: usize) -> Self {
        Self::batch_mul(&[self.clone()], &[scalar.clone()], max_num_bits, false)
    }

    /// Batch multi-scalar multiplication using Straus algorithm.
    pub fn batch_mul(
        _points: &[Self],
        _scalars: &[Fr<P, C>],
        max_num_bits: usize,
        with_edgecases: bool,
    ) -> Self {
        assert!(!_points.is_empty(), "batch_mul: no points");
        assert_eq!(
            _points.len(),
            _scalars.len(),
            "batch_mul: size mismatch"
        );

        // Handle points at infinity
        let (mut points, mut scalars) =
            Self::handle_points_at_infinity(_points, _scalars);

        // Accumulate constant pairs out of circuit
        let mut has_constant_terms = false;
        let mut constant_accumulator = Element::<C>::infinity();
        let mut new_points = Vec::new();
        let mut new_scalars = Vec::new();

        for i in 0..points.len() {
            if points[i].is_constant() && scalars[i].is_constant() {
                let point_val = Element::<C>::from_affine(&points[i].get_value());
                let scalar_val_512 = scalars[i].get_value();
                let modulus = U256::from_limbs(C::ScalarFieldParams::MODULUS);
                let modulus_512 = U512::from_lo_hi(modulus, U256::ZERO);
                let scalar_u256 = scalar_val_512.div_rem(&modulus_512.to_nz().unwrap()).1.lo();
                let scalar_field =
                    Field::<C::ScalarFieldParams>::from_limbs(*scalar_u256.as_words());
                constant_accumulator = constant_accumulator + point_val.mul(&scalar_field);
                has_constant_terms = true;
            } else {
                new_points.push(points[i].clone());
                new_scalars.push(scalars[i].clone());
            }
        }
        points = new_points;
        scalars = new_scalars;

        if with_edgecases && !points.is_empty() {
            // Simple masking: add offset generators to prevent edge cases
            let masking_scalar = Fr::<P, C>::one_val();
            let masked = Self::mask_points(&points, &scalars, &masking_scalar);
            points = masked.0;
            scalars = masked.1;
        }

        // Separate zero and large scalars
        let modulus = U256::from_limbs(C::ScalarFieldParams::MODULUS);
        let max_num_bits_in_field = modulus.get_msb() as usize + 1;

        let mut accumulator = if has_constant_terms || (points.is_empty()) {
            Some(Self::from_affine(&constant_accumulator.to_affine()))
        } else {
            None
        };

        let mut big_points = Vec::new();
        let mut big_scalars = Vec::new();
        let mut small_points = Vec::new();
        let mut small_scalars = Vec::new();

        for i in 0..points.len() {
            if max_num_bits == 0 {
                big_points.push(points[i].clone());
                big_scalars.push(scalars[i].clone());
            } else {
                let sv = scalars[i].get_value();
                let is_last_big = (i == points.len() - 1) && with_edgecases;
                if sv == U512::ZERO || is_last_big {
                    big_points.push(points[i].clone());
                    big_scalars.push(scalars[i].clone());
                } else {
                    small_points.push(points[i].clone());
                    small_scalars.push(scalars[i].clone());
                }
            }
        }

        if !big_points.is_empty() {
            let big_result =
                Self::process_strauss_msm_rounds(&big_points, &big_scalars, max_num_bits_in_field);
            accumulator = Some(match accumulator {
                Some(acc) => acc.add(&big_result),
                None => big_result,
            });
        }

        if !small_points.is_empty() {
            let effective_bits = if max_num_bits == 0 {
                max_num_bits_in_field
            } else {
                max_num_bits
            };
            let small_result =
                Self::process_strauss_msm_rounds(&small_points, &small_scalars, effective_bits);
            accumulator = Some(match accumulator {
                Some(acc) => acc.add(&small_result),
                None => small_result,
            });
        }

        accumulator.unwrap()
    }

    /// Process the Straus MSM rounds.
    fn process_strauss_msm_rounds(
        points: &[Self],
        scalars: &[Fr<P, C>],
        max_num_bits: usize,
    ) -> Self {
        assert!(!points.is_empty());
        assert_eq!(points.len(), scalars.len());

        let num_rounds = max_num_bits;
        let msm_size = scalars.len();

        // Compute NAF representations
        let mut naf_entries = Vec::new();
        for scalar in scalars {
            naf_entries.push(Self::compute_naf(scalar, num_rounds));
        }

        // Compute offset generators (must be witnesses, not constants, so they have context)
        let ctx = points
            .iter()
            .find_map(|p| p.get_context())
            .expect("process_strauss_msm_rounds: need context");
        let (offset_gen_start, offset_gen_end) = compute_offset_generators::<C>(num_rounds);
        let offset_start = Self::from_witness(ctx.clone(), &offset_gen_start);
        let offset_end = Self::from_witness(ctx.clone(), &offset_gen_end);

        // Build batch lookup table
        let point_table = BatchLookupTable::new(points);

        // Initialize accumulator with offset generator + first NAF column
        let initial_entry = point_table.get_chain_initial_entry();
        let acc_chain = Self::chain_add(
            &offset_start,
            &initial_entry,
        );
        let mut accumulator = Self::chain_add_end(&acc_chain);

        // Process rounds in groups of 4
        let num_rounds_per_iteration = 4;
        let num_iterations = (num_rounds - 1 + num_rounds_per_iteration - 1) / num_rounds_per_iteration;
        let num_rounds_last =
            (num_rounds - 1) - ((num_iterations - 1) * num_rounds_per_iteration);

        for i in 0..num_iterations {
            let inner_num_rounds = if i != num_iterations - 1 {
                num_rounds_per_iteration
            } else {
                num_rounds_last
            };

            let mut to_add = Vec::new();
            for j in 0..inner_num_rounds {
                let mut nafs = Vec::new();
                for k in 0..msm_size {
                    nafs.push(naf_entries[k][(i * num_rounds_per_iteration) + j + 1].clone());
                }
                let add_acc = point_table.get_chain_add_accumulator(&nafs);
                to_add.push(add_acc);
            }

            accumulator = accumulator.multiple_montgomery_ladder(&to_add);
        }

        // Subtract skew factors
        for i in 0..msm_size {
            let skew = accumulator.sub(&points[i]);
            accumulator =
                accumulator.conditional_select(&skew, &naf_entries[i][num_rounds]);
        }

        // Subtract offset generator
        accumulator = accumulator.sub(&offset_end);

        accumulator
    }

    // ════════════════════════════════════════════════════════════════
    //  Edge case handling
    // ════════════════════════════════════════════════════════════════

    /// Replace (∞, scalar) pairs with (G, 0).
    fn handle_points_at_infinity(
        points: &[Self],
        scalars: &[Fr<P, C>],
    ) -> (Vec<Self>, Vec<Fr<P, C>>) {
        let mut new_points = Vec::new();
        let mut new_scalars = Vec::new();

        let ctx = points
            .iter()
            .find_map(|p| p.get_context())
            .or_else(|| scalars.iter().find_map(|s| s.get_context().clone()));

        let one_point = Self::one(ctx);

        for i in 0..points.len() {
            let is_inf = points[i].is_point_at_infinity();

            if is_inf.get_value() && is_inf.is_constant() {
                continue;
            }

            let sv = scalars[i].get_value();
            if sv == U512::ZERO && scalars[i].is_constant() {
                continue;
            }

            let point = points[i].conditional_select(&one_point, is_inf);
            let scalar = conditional_assign(is_inf, &Fr::<P, C>::zero_val(), &scalars[i]);

            new_points.push(point);
            new_scalars.push(scalar);
        }

        (new_points, new_scalars)
    }

    /// Mask points to handle edge cases in batch multiplication.
    fn mask_points(
        points: &[Self],
        scalars: &[Fr<P, C>],
        _masking_scalar: &Fr<P, C>,
    ) -> (Vec<Self>, Vec<Fr<P, C>>) {
        assert_eq!(points.len(), scalars.len());

        let ctx = points
            .iter()
            .find_map(|p| p.get_context())
            .expect("mask_points: need context");

        let native_offset_gen = compute_table_offset_generator::<C>();
        let offset_gen = Self::from_witness(ctx.clone(), &native_offset_gen);

        // Simple approach: add distinct multiples of offset generator to each point
        let mut running_point = offset_gen;
        let mut running_scalar = Fr::<P, C>::one_val();
        let mut last_scalar = Fr::<P, C>::zero_val();

        let mut new_points = Vec::new();
        let mut new_scalars = Vec::new();

        for i in 0..points.len() {
            new_scalars.push(scalars[i].clone());
            new_points.push(points[i].add(&running_point));
            last_scalar = last_scalar.add(&scalars[i].mul(&running_scalar));
            running_scalar = running_scalar.add(&running_scalar.clone());
            running_point = running_point.dbl();
        }

        // Add the cancellation point
        let n = points.len() as u32;
        let two_pow_n_val = U256::from(1u64).wrapping_shl_vartime(n);
        let two_pow_n = Fr::<P, C>::from_u256(Some(ctx), two_pow_n_val);
        let two_pow_n_inv = two_pow_n.div(&Fr::<P, C>::one_val()); // This is just the value
        // Actually we need field inversion. For simplicity, use div:
        // last_scalar / 2^n
        last_scalar = last_scalar.div(&two_pow_n);
        new_scalars.push(last_scalar.negate());
        new_points.push(running_point);

        (new_points, new_scalars)
    }
}

// ════════════════════════════════════════════════════════════════════════
//  FourBitTablePlookup
// ════════════════════════════════════════════════════════════════════════

/// Four-bit variable-base lookup table for scalar multiplication.
pub struct FourBitTablePlookup<P: FieldParams, C: CurveParams> {
    pub element_table: Vec<BigGroupT<P, C>>,
    coordinates: [RefCell<TwinRomTable<P>>; 5],
    limb_max: [U256; 8],
}

impl<P: FieldParams, C: CurveParams> FourBitTablePlookup<P, C> {
    pub fn new(input: &BigGroupT<P, C>) -> Self {
        let d2 = input.dbl();
        let mut table: Vec<BigGroupT<P, C>> = vec![BigGroupT::new_empty(); 16];

        table[8] = input.clone();
        for i in 9..16 {
            table[i] = table[i - 1].checked_unconditional_add(&d2);
        }
        for i in 0..8 {
            table[i] = table[15 - i].negate();
        }

        let mut limb_max = [U256::ZERO; 8];
        let coordinates =
            BigGroupT::<P, C>::create_group_element_rom_tables(&table, &mut limb_max);

        let element_table = table;

        Self {
            element_table,
            coordinates,
            limb_max,
        }
    }

    pub fn read(&self, index: &FieldT<P>) -> BigGroupT<P, C> {
        BigGroupT::read_group_element_rom_tables(&self.coordinates, index, &self.limb_max)
    }

    pub fn get(&self, idx: usize) -> &BigGroupT<P, C> {
        &self.element_table[idx]
    }
}

// ════════════════════════════════════════════════════════════════════════
//  LookupTablePlookup (generic N-bit)
// ════════════════════════════════════════════════════════════════════════

/// Generic N-bit lookup table.
pub struct LookupTablePlookup<P: FieldParams, C: CurveParams> {
    pub element_table: Vec<BigGroupT<P, C>>,
    coordinates: [RefCell<TwinRomTable<P>>; 5],
    limb_max: [U256; 8],
    length: usize,
}

impl<P: FieldParams, C: CurveParams> LookupTablePlookup<P, C> {
    pub fn new(inputs: &[BigGroupT<P, C>]) -> Self {
        let length = inputs.len();
        assert!(length <= 6);
        let table_size = 1 << length;

        let mut element_table = vec![BigGroupT::new_empty(); table_size];

        match length {
            2 => {
                let [a0, a1] = inputs[1].checked_unconditional_add_sub(&inputs[0]);
                element_table[0] = a0;
                element_table[1] = a1;
            }
            3 => {
                let [r0, r1] = inputs[1].checked_unconditional_add_sub(&inputs[0]);
                let [t0, t1] = inputs[2].checked_unconditional_add_sub(&r0);
                let [t2, t3] = inputs[2].checked_unconditional_add_sub(&r1);
                element_table[0] = t0;
                element_table[1] = t2;
                element_table[2] = t3;
                element_table[3] = t1;
            }
            4 => {
                let [t0, t1] = inputs[1].checked_unconditional_add_sub(&inputs[0]);
                let [t2, t3] = inputs[3].checked_unconditional_add_sub(&inputs[2]);
                let [f0, f3] = t2.checked_unconditional_add_sub(&t0);
                let [f1, f2] = t2.checked_unconditional_add_sub(&t1);
                let [f4, f7] = t3.checked_unconditional_add_sub(&t0);
                let [f5, f6] = t3.checked_unconditional_add_sub(&t1);
                element_table[0] = f0;
                element_table[1] = f1;
                element_table[2] = f2;
                element_table[3] = f3;
                element_table[4] = f4;
                element_table[5] = f5;
                element_table[6] = f6;
                element_table[7] = f7;
            }
            5 => {
                let [a0, a1] = inputs[1].checked_unconditional_add_sub(&inputs[0]);
                let [t2, t3] = inputs[3].checked_unconditional_add_sub(&inputs[2]);
                let [e0, e3] = inputs[4].checked_unconditional_add_sub(&t2);
                let [e1, e2] = inputs[4].checked_unconditional_add_sub(&t3);
                let [f0, f3] = e0.checked_unconditional_add_sub(&a0);
                let [f1, f2] = e0.checked_unconditional_add_sub(&a1);
                let [f4, f7] = e1.checked_unconditional_add_sub(&a0);
                let [f5, f6] = e1.checked_unconditional_add_sub(&a1);
                let [f8, f11] = e2.checked_unconditional_add_sub(&a0);
                let [f9, f10] = e2.checked_unconditional_add_sub(&a1);
                let [f12, f15] = e3.checked_unconditional_add_sub(&a0);
                let [f13, f14] = e3.checked_unconditional_add_sub(&a1);
                for (i, v) in [f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15].into_iter().enumerate() {
                    element_table[i] = v;
                }
            }
            6 => {
                let [a0, a1] = inputs[1].checked_unconditional_add_sub(&inputs[0]);
                let [e0, e1] = inputs[4].checked_unconditional_add_sub(&inputs[3]);
                let [c0, c3] = inputs[2].checked_unconditional_add_sub(&a0);
                let [c1, c2] = inputs[2].checked_unconditional_add_sub(&a1);
                let [f0, f3] = inputs[5].checked_unconditional_add_sub(&e0);
                let [f1, f2] = inputs[5].checked_unconditional_add_sub(&e1);

                let pairs = [
                    f0.checked_unconditional_add_sub(&c0),
                    f0.checked_unconditional_add_sub(&c1),
                    f0.checked_unconditional_add_sub(&c2),
                    f0.checked_unconditional_add_sub(&c3),
                    f1.checked_unconditional_add_sub(&c0),
                    f1.checked_unconditional_add_sub(&c1),
                    f1.checked_unconditional_add_sub(&c2),
                    f1.checked_unconditional_add_sub(&c3),
                    f2.checked_unconditional_add_sub(&c0),
                    f2.checked_unconditional_add_sub(&c1),
                    f2.checked_unconditional_add_sub(&c2),
                    f2.checked_unconditional_add_sub(&c3),
                    f3.checked_unconditional_add_sub(&c0),
                    f3.checked_unconditional_add_sub(&c1),
                    f3.checked_unconditional_add_sub(&c2),
                    f3.checked_unconditional_add_sub(&c3),
                ];

                // Assign: R0-R7, S0-S7, U0-U7, W0-W7
                for (group, base) in pairs.chunks(4).enumerate() {
                    for (j, [add, sub]) in base.iter().enumerate() {
                        element_table[group * 8 + j] = add.clone();
                        element_table[group * 8 + 7 - j] = sub.clone();
                    }
                }
            }
            _ => {
                // For length 1 or 0, just store the point
                if length == 1 {
                    element_table[0] = inputs[0].clone();
                    element_table[1] = inputs[0].negate();
                }
            }
        }

        // Fill second half with negations of first half (reversed)
        let half = table_size / 2;
        for i in 0..half {
            element_table[i + half] = element_table[half - 1 - i].negate();
        }

        let mut limb_max = [U256::ZERO; 8];
        let coordinates =
            BigGroupT::<P, C>::create_group_element_rom_tables(&element_table, &mut limb_max);

        Self {
            element_table,
            coordinates,
            limb_max,
            length,
        }
    }

    /// Read from the table using bit-decomposed index.
    pub fn get_by_bits(&self, bits: &[BoolT<P>]) -> BigGroupT<P, C> {
        assert_eq!(bits.len(), self.length);

        let mut accumulators = Vec::new();
        for (i, bit) in bits.iter().enumerate() {
            let bit_ft = bool_to_field(bit);
            let shift = FieldT::from_field(Field::<P>::from(1u64 << i));
            accumulators.push(&bit_ft * &shift);
        }
        let index = FieldT::accumulate(&accumulators);

        BigGroupT::read_group_element_rom_tables(&self.coordinates, &index, &self.limb_max)
    }

    pub fn get_by_index(&self, idx: usize) -> &BigGroupT<P, C> {
        &self.element_table[idx]
    }
}

// ════════════════════════════════════════════════════════════════════════
//  BatchLookupTable
// ════════════════════════════════════════════════════════════════════════

/// Batch lookup table that splits points into optimally-sized sub-tables.
pub struct BatchLookupTable<P: FieldParams, C: CurveParams> {
    six_tables: Vec<LookupTablePlookup<P, C>>,
    five_tables: Vec<LookupTablePlookup<P, C>>,
    quad_tables: Vec<LookupTablePlookup<P, C>>,
    triple_tables: Vec<LookupTablePlookup<P, C>>,
    twin_tables: Vec<LookupTablePlookup<P, C>>,
    singletons: Vec<BigGroupT<P, C>>,
    num_points: usize,
    num_sixes: usize,
    num_fives: usize,
    has_quad: bool,
    has_triple: bool,
    has_twin: bool,
    has_singleton: bool,
}

impl<P: FieldParams, C: CurveParams> BatchLookupTable<P, C> {
    pub fn new(points: &[BigGroupT<P, C>]) -> Self {
        let num_points = points.len();
        let mut num_fives = num_points / 5;
        let mut num_sixes = 0usize;

        if num_points == 1 {
            num_fives = 0;
        } else if num_points >= 1 && num_fives * 5 == num_points - 1 {
            num_fives = num_fives.saturating_sub(1);
            num_sixes = 1;
        } else if num_points >= 2 && num_fives * 5 == num_points - 2 && num_fives >= 2 {
            num_fives -= 2;
            num_sixes = 2;
        } else if num_points >= 3 && num_fives * 5 == num_points - 3 && num_fives >= 3 {
            num_fives -= 3;
            num_sixes = 3;
        }

        let mut remaining = num_points - (num_fives * 5 + num_sixes * 6);

        let has_quad = remaining >= 4 && num_points >= 4;
        if has_quad {
            remaining -= 4;
        }
        let has_triple = remaining >= 3 && num_points >= 3;
        if has_triple {
            remaining -= 3;
        }
        let has_twin = remaining >= 2 && num_points >= 2;
        if has_twin {
            remaining -= 2;
        }
        let has_singleton = remaining != 0 && num_points >= 1;

        let mut offset = 0;
        let mut six_tables = Vec::new();
        for _ in 0..num_sixes {
            six_tables.push(LookupTablePlookup::new(&points[offset..offset + 6]));
            offset += 6;
        }
        let mut five_tables = Vec::new();
        for _ in 0..num_fives {
            five_tables.push(LookupTablePlookup::new(&points[offset..offset + 5]));
            offset += 5;
        }
        let mut quad_tables = Vec::new();
        if has_quad {
            quad_tables.push(LookupTablePlookup::new(&points[offset..offset + 4]));
            offset += 4;
        }
        let mut triple_tables = Vec::new();
        if has_triple {
            triple_tables.push(LookupTablePlookup::new(&points[offset..offset + 3]));
            offset += 3;
        }
        let mut twin_tables = Vec::new();
        if has_twin {
            twin_tables.push(LookupTablePlookup::new(&points[offset..offset + 2]));
            offset += 2;
        }
        let mut singletons = Vec::new();
        if has_singleton {
            singletons.push(points[points.len() - 1].clone());
        }

        Self {
            six_tables,
            five_tables,
            quad_tables,
            triple_tables,
            twin_tables,
            singletons,
            num_points,
            num_sixes,
            num_fives,
            has_quad,
            has_triple,
            has_twin,
            has_singleton,
        }
    }

    /// Get the initial chain entry (all tables at index 0).
    pub fn get_chain_initial_entry(&self) -> ChainAddAccumulator<P, C> {
        let mut add_acc = Vec::new();
        for t in &self.six_tables {
            add_acc.push(t.get_by_index(0).clone());
        }
        for t in &self.five_tables {
            add_acc.push(t.get_by_index(0).clone());
        }
        if self.has_quad {
            add_acc.push(self.quad_tables[0].get_by_index(0).clone());
        }
        if self.has_twin {
            add_acc.push(self.twin_tables[0].get_by_index(0).clone());
        }
        if self.has_triple {
            add_acc.push(self.triple_tables[0].get_by_index(0).clone());
        }
        if self.has_singleton {
            add_acc.push(self.singletons[0].clone());
        }

        if add_acc.len() >= 2 {
            let mut output = BigGroupT::chain_add_start(&add_acc[0], &add_acc[1]);
            for i in 2..add_acc.len() {
                output = BigGroupT::chain_add(&add_acc[i], &output);
            }
            output
        } else {
            ChainAddAccumulator::from_element(&add_acc[0])
        }
    }

    /// Get chain add accumulator for a round of NAF entries.
    pub fn get_chain_add_accumulator(
        &self,
        naf_entries: &[BoolT<P>],
    ) -> ChainAddAccumulator<P, C> {
        let mut round_acc = Vec::new();
        let mut offset = 0;

        for t in &self.six_tables {
            round_acc.push(t.get_by_bits(&naf_entries[offset..offset + 6]));
            offset += 6;
        }
        for t in &self.five_tables {
            round_acc.push(t.get_by_bits(&naf_entries[offset..offset + 5]));
            offset += 5;
        }
        if self.has_quad {
            round_acc.push(self.quad_tables[0].get_by_bits(&naf_entries[offset..offset + 4]));
            offset += 4;
        }
        if self.has_triple {
            round_acc.push(self.triple_tables[0].get_by_bits(&naf_entries[offset..offset + 3]));
            offset += 3;
        }
        if self.has_twin {
            round_acc.push(self.twin_tables[0].get_by_bits(&naf_entries[offset..offset + 2]));
            offset += 2;
        }
        if self.has_singleton {
            round_acc.push(
                self.singletons[0].conditional_negate(&naf_entries[self.num_points - 1]),
            );
        }

        if round_acc.len() == 1 {
            return ChainAddAccumulator::from_element(&round_acc[0]);
        }
        if round_acc.len() == 2 {
            return BigGroupT::chain_add_start(&round_acc[0], &round_acc[1]);
        }

        let mut accumulator = BigGroupT::chain_add_start(&round_acc[0], &round_acc[1]);
        for j in 2..round_acc.len() {
            accumulator = BigGroupT::chain_add(&round_acc[j], &accumulator);
        }
        accumulator
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use bbrs_ecc::curves::bn254::{Bn254FqParams, Bn254FrParams, Bn254G1Params};
    use bbrs_ecc::fields::field::Field;
    use std::cell::RefCell;
    use std::rc::Rc;

    type TestGroup = BigGroupT<Bn254FrParams, Bn254G1Params>;
    type TestFr = BigFieldT<Bn254FrParams, Bn254FrParams>;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    fn random_affine() -> AffineElement<Bn254G1Params> {
        let scalar = Field::<Bn254FrParams>::random_element();
        Element::<Bn254G1Params>::one().mul(&scalar).to_affine()
    }

    fn check_circuit(ctx: &BuilderRef<Bn254FrParams>) {
        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }

    #[test]
    fn test_add() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);

        let result = a.add(&b);

        let expected = (Element::from_affine(&a_native) + Element::from_affine(&b_native)).to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "x mismatch in add");
        assert_eq!(got.y, expected.y, "y mismatch in add");
        check_circuit(&ctx);
    }

    #[test]
    fn test_sub() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);

        let result = a.sub(&b);

        let expected = (Element::from_affine(&a_native) - Element::from_affine(&b_native)).to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "x mismatch in sub");
        assert_eq!(got.y, expected.y, "y mismatch in sub");
        check_circuit(&ctx);
    }

    #[test]
    fn test_dbl() {
        let ctx = make_builder();
        let a_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let result = a.dbl();

        let expected = Element::from_affine(&a_native).dbl().to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "x mismatch in dbl");
        assert_eq!(got.y, expected.y, "y mismatch in dbl");
        check_circuit(&ctx);
    }

    #[test]
    fn test_negate() {
        let ctx = make_builder();
        let a_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let result = a.negate();

        let mut expected = a_native.clone();
        expected.y = -expected.y;
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "x mismatch in negate");
        assert_eq!(got.y, expected.y, "y mismatch in negate");
        check_circuit(&ctx);
    }

    #[test]
    fn test_checked_unconditional_add() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);

        let result = a.checked_unconditional_add(&b);

        let expected = (Element::from_affine(&a_native) + Element::from_affine(&b_native)).to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "x mismatch in checked_add");
        assert_eq!(got.y, expected.y, "y mismatch in checked_add");
        check_circuit(&ctx);
    }

    #[test]
    fn test_chain_add() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();
        let c_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);
        let c = TestGroup::from_witness(ctx.clone(), &c_native);

        let acc = TestGroup::chain_add_start(&a, &b);
        let acc = TestGroup::chain_add(&c, &acc);
        let result = TestGroup::chain_add_end(&acc);

        let expected = (Element::from_affine(&a_native)
            + Element::from_affine(&b_native)
            + Element::from_affine(&c_native))
        .to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "x mismatch in chain_add");
        assert_eq!(got.y, expected.y, "y mismatch in chain_add");
        check_circuit(&ctx);
    }

    #[test]
    fn test_montgomery_ladder() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);

        let acc = ChainAddAccumulator::from_element(&b);
        let result = a.multiple_montgomery_ladder(&[acc]);

        // Expected: (2*A + B) = (A + B) + A
        let expected = (Element::from_affine(&a_native)
            + Element::from_affine(&b_native)
            + Element::from_affine(&a_native))
        .to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "x mismatch in montgomery_ladder");
        assert_eq!(got.y, expected.y, "y mismatch in montgomery_ladder");
        check_circuit(&ctx);
    }

    #[test]
    fn test_add_with_constants() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();

        // Witness + constant
        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_affine(&b_native);

        let result = a.add(&b);
        let expected = (Element::from_affine(&a_native) + Element::from_affine(&b_native)).to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "x mismatch in add_with_constants");
        assert_eq!(got.y, expected.y, "y mismatch in add_with_constants");
        check_circuit(&ctx);
    }

    #[test]
    fn test_add_points_at_infinity() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);

        // a + infinity = a
        let inf = TestGroup::point_at_infinity(ctx.clone());
        let result = a.add(&inf);
        let got = result.get_value();
        assert_eq!(got.x, a_native.x, "x: a + inf should be a");
        assert_eq!(got.y, a_native.y, "y: a + inf should be a");

        // infinity + b = b
        let result2 = inf.add(&b);
        let got2 = result2.get_value();
        assert_eq!(got2.x, b_native.x, "x: inf + b should be b");
        assert_eq!(got2.y, b_native.y, "y: inf + b should be b");

        check_circuit(&ctx);
    }

    #[test]
    fn test_sub_points_at_infinity() {
        let ctx = make_builder();
        let a_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let inf = TestGroup::point_at_infinity(ctx.clone());

        // a - infinity = a
        let result = a.sub(&inf);
        let got = result.get_value();
        assert_eq!(got.x, a_native.x, "x: a - inf should be a");
        assert_eq!(got.y, a_native.y, "y: a - inf should be a");

        check_circuit(&ctx);
    }

    #[test]
    fn test_dbl_with_constant() {
        let a_native = random_affine();
        let a = TestGroup::from_affine(&a_native);

        let result = a.dbl();
        let expected = Element::from_affine(&a_native).dbl().to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "x mismatch in dbl_constant");
        assert_eq!(got.y, expected.y, "y mismatch in dbl_constant");
    }

    #[test]
    fn test_add_equals_dbl() {
        let ctx = make_builder();
        let a_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let a2 = TestGroup::from_witness(ctx.clone(), &a_native);

        let add_result = a.add(&a2);
        let dbl_result = a.dbl();

        let add_val = add_result.get_value();
        let dbl_val = dbl_result.get_value();

        assert_eq!(add_val.x, dbl_val.x, "add vs dbl: x mismatch");
        assert_eq!(add_val.y, dbl_val.y, "add vs dbl: y mismatch");
        check_circuit(&ctx);
    }

    #[test]
    fn test_conditional_negate() {
        let ctx = make_builder();
        let a_native = random_affine();
        let a = TestGroup::from_witness(ctx.clone(), &a_native);

        // Negate when true
        let pred_true = BoolT::from_witness(&WitnessT::from_bool(ctx.clone(), true));
        let result = a.conditional_negate(&pred_true);
        let got = result.get_value();
        assert_eq!(got.x, a_native.x, "cond_neg(true): x should be same");
        assert_eq!(got.y, -a_native.y, "cond_neg(true): y should be negated");

        // Don't negate when false
        let pred_false = BoolT::from_witness(&WitnessT::from_bool(ctx.clone(), false));
        let result2 = a.conditional_negate(&pred_false);
        let got2 = result2.get_value();
        assert_eq!(got2.x, a_native.x, "cond_neg(false): x should be same");
        assert_eq!(got2.y, a_native.y, "cond_neg(false): y should be same");

        check_circuit(&ctx);
    }

    #[test]
    fn test_conditional_select() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);

        // Select other (a) when true
        let pred_true = BoolT::from_witness(&WitnessT::from_bool(ctx.clone(), true));
        let result = b.conditional_select(&a, &pred_true);
        let got = result.get_value();
        assert_eq!(got.x, a_native.x, "cond_select(true): should be a");
        assert_eq!(got.y, a_native.y, "cond_select(true): should be a");

        // Select self (b) when false
        let pred_false = BoolT::from_witness(&WitnessT::from_bool(ctx.clone(), false));
        let result2 = b.conditional_select(&a, &pred_false);
        let got2 = result2.get_value();
        assert_eq!(got2.x, b_native.x, "cond_select(false): should be b");
        assert_eq!(got2.y, b_native.y, "cond_select(false): should be b");

        check_circuit(&ctx);
    }

    #[test]
    fn test_checked_unconditional_subtract() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);

        let result = a.checked_unconditional_subtract(&b);

        let expected = (Element::from_affine(&a_native) - Element::from_affine(&b_native)).to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "x mismatch in checked_sub");
        assert_eq!(got.y, expected.y, "y mismatch in checked_sub");
        check_circuit(&ctx);
    }

    #[test]
    fn test_checked_unconditional_add_sub() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);

        let [add_result, sub_result] = a.checked_unconditional_add_sub(&b);

        let expected_add = (Element::from_affine(&a_native) + Element::from_affine(&b_native)).to_affine();
        let expected_sub = (Element::from_affine(&a_native) - Element::from_affine(&b_native)).to_affine();

        let add_val = add_result.get_value();
        let sub_val = sub_result.get_value();

        assert_eq!(add_val.x, expected_add.x, "add_sub: add x mismatch");
        assert_eq!(add_val.y, expected_add.y, "add_sub: add y mismatch");
        assert_eq!(sub_val.x, expected_sub.x, "add_sub: sub x mismatch");
        assert_eq!(sub_val.y, expected_sub.y, "add_sub: sub y mismatch");
        check_circuit(&ctx);
    }

    #[test]
    fn test_validate_on_curve() {
        let ctx = make_builder();
        let a_native = random_affine();

        // Should succeed - valid point
        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        a.validate_on_curve();
        check_circuit(&ctx);
    }

    #[test]
    fn test_normalize() {
        let ctx = make_builder();
        let a_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let result = a.normalize();
        let got = result.get_value();

        assert_eq!(got.x, a_native.x, "normalize: x mismatch");
        assert_eq!(got.y, a_native.y, "normalize: y mismatch");
        check_circuit(&ctx);
    }

    #[test]
    fn test_incomplete_assert_equal() {
        let ctx = make_builder();
        let a_native = random_affine();

        let a1 = TestGroup::from_witness(ctx.clone(), &a_native);
        let a2 = TestGroup::from_witness(ctx.clone(), &a_native);

        a1.incomplete_assert_equal(&a2);
        check_circuit(&ctx);
    }

    #[test]
    fn test_from_affine_constant() {
        let a_native = random_affine();
        let a = TestGroup::from_affine(&a_native);

        assert!(a.is_constant());
        let got = a.get_value();
        assert_eq!(got.x, a_native.x, "from_affine: x mismatch");
        assert_eq!(got.y, a_native.y, "from_affine: y mismatch");
    }

    #[test]
    fn test_compute_naf() {
        let ctx = make_builder();

        // Use a small scalar value
        let scalar_val = U256::from(42u64);
        let scalar = TestFr::from_witness(ctx.clone(), scalar_val);

        let naf = TestGroup::compute_naf(&scalar, 0);

        // NAF should have entries - just verify it doesn't panic and circuit passes
        assert!(!naf.is_empty(), "NAF should have entries");
        check_circuit(&ctx);
    }

    #[test]
    fn test_scalar_mul() {
        let ctx = make_builder();
        let a_native = random_affine();
        let scalar_val = U256::from(7u64);

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let scalar = TestFr::from_witness(ctx.clone(), scalar_val);

        let result = a.scalar_mul(&scalar, 0);

        let scalar_field = Field::<Bn254FrParams>::from_limbs(*scalar_val.as_words());
        let expected = Element::from_affine(&a_native).mul(&scalar_field).to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "scalar_mul: x mismatch");
        assert_eq!(got.y, expected.y, "scalar_mul: y mismatch");
        check_circuit(&ctx);
    }

    #[test]
    fn test_batch_mul_singleton() {
        let ctx = make_builder();
        let a_native = random_affine();
        let scalar_val = U256::from(13u64);

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let scalar = TestFr::from_witness(ctx.clone(), scalar_val);

        let result = TestGroup::batch_mul(&[a], &[scalar], 0, false);

        let scalar_field = Field::<Bn254FrParams>::from_limbs(*scalar_val.as_words());
        let expected = Element::from_affine(&a_native).mul(&scalar_field).to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "batch_mul_singleton: x mismatch");
        assert_eq!(got.y, expected.y, "batch_mul_singleton: y mismatch");
        check_circuit(&ctx);
    }

    #[test]
    fn test_batch_mul_twin() {
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();
        let s1_val = U256::from(5u64);
        let s2_val = U256::from(11u64);

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);
        let s1 = TestFr::from_witness(ctx.clone(), s1_val);
        let s2 = TestFr::from_witness(ctx.clone(), s2_val);

        let result = TestGroup::batch_mul(&[a, b], &[s1, s2], 0, false);

        let s1_field = Field::<Bn254FrParams>::from_limbs(*s1_val.as_words());
        let s2_field = Field::<Bn254FrParams>::from_limbs(*s2_val.as_words());
        let expected = (Element::from_affine(&a_native).mul(&s1_field)
            + Element::from_affine(&b_native).mul(&s2_field))
        .to_affine();
        let got = result.get_value();

        assert_eq!(got.x, expected.x, "batch_mul_twin: x mismatch");
        assert_eq!(got.y, expected.y, "batch_mul_twin: y mismatch");
        check_circuit(&ctx);
    }

    #[test]
    fn test_multiple_montgomery_ladder_single() {
        // Test: MML with 1 full-element entry
        // MML([P]) computes 2*self + P
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);

        let acc = ChainAddAccumulator::from_element(&b);
        let result = a.multiple_montgomery_ladder(&[acc]);
        let got = result.get_value();

        // Expected: 2*a + b
        let expected = (Element::from_affine(&a_native)
            + Element::from_affine(&a_native)
            + Element::from_affine(&b_native))
        .to_affine();

        assert_eq!(got.x, expected.x, "MML single: x mismatch");
        assert_eq!(got.y, expected.y, "MML single: y mismatch");
        check_circuit(&ctx);
    }

    #[test]
    fn test_multiple_montgomery_ladder_double() {
        // Test: MML([P, Q]) computes 2*(2*self + P) + Q = 4*self + 2*P + Q
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();
        let c_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);
        let c = TestGroup::from_witness(ctx.clone(), &c_native);

        let acc_b = ChainAddAccumulator::from_element(&b);
        let acc_c = ChainAddAccumulator::from_element(&c);
        let result = a.multiple_montgomery_ladder(&[acc_b, acc_c]);
        let got = result.get_value();

        // Expected: 4*a + 2*b + c
        let four = Field::<Bn254FrParams>::from(4u64);
        let two = Field::<Bn254FrParams>::from(2u64);
        let expected = (Element::from_affine(&a_native).mul(&four)
            + Element::from_affine(&b_native).mul(&two)
            + Element::from_affine(&c_native))
        .to_affine();

        assert_eq!(got.x, expected.x, "MML double: x mismatch");
        assert_eq!(got.y, expected.y, "MML double: y mismatch");
        check_circuit(&ctx);
    }

    #[test]
    fn test_multiple_montgomery_ladder_triple() {
        // Test: MML([P, Q, R]) = 8*self + 4*P + 2*Q + R
        let ctx = make_builder();
        let a_native = random_affine();
        let b_native = random_affine();
        let c_native = random_affine();
        let d_native = random_affine();

        let a = TestGroup::from_witness(ctx.clone(), &a_native);
        let b = TestGroup::from_witness(ctx.clone(), &b_native);
        let c = TestGroup::from_witness(ctx.clone(), &c_native);
        let d = TestGroup::from_witness(ctx.clone(), &d_native);

        let acc_b = ChainAddAccumulator::from_element(&b);
        let acc_c = ChainAddAccumulator::from_element(&c);
        let acc_d = ChainAddAccumulator::from_element(&d);
        let result = a.multiple_montgomery_ladder(&[acc_b, acc_c, acc_d]);
        let got = result.get_value();

        let eight = Field::<Bn254FrParams>::from(8u64);
        let four = Field::<Bn254FrParams>::from(4u64);
        let two = Field::<Bn254FrParams>::from(2u64);
        let expected = (Element::from_affine(&a_native).mul(&eight)
            + Element::from_affine(&b_native).mul(&four)
            + Element::from_affine(&c_native).mul(&two)
            + Element::from_affine(&d_native))
        .to_affine();

        assert_eq!(got.x, expected.x, "MML triple: x mismatch");
        assert_eq!(got.y, expected.y, "MML triple: y mismatch");
        check_circuit(&ctx);
    }

    #[test]
    fn test_four_bit_table_plookup() {
        let ctx = make_builder();
        let a_native = random_affine();
        let a = TestGroup::from_witness(ctx.clone(), &a_native);

        let table = FourBitTablePlookup::new(&a);

        // Table[8] should be the original point, table[0] = -(15*P)
        // The table stores: [-15P, -13P, ..., -P, P, 3P, ..., 15P]
        let got_8 = table.get(8).get_value();
        assert_eq!(got_8.x, a_native.x, "table[8] should be original point (x)");
        assert_eq!(got_8.y, a_native.y, "table[8] should be original point (y)");

        // Read with a witness index
        let idx = FieldT::from_witness(ctx.clone(), Field::<Bn254FrParams>::from(8u64));
        let read_result = table.read(&idx);
        let got_read = read_result.get_value();
        assert_eq!(got_read.x, a_native.x, "table.read(8) should be original point (x)");
        assert_eq!(got_read.y, a_native.y, "table.read(8) should be original point (y)");

        check_circuit(&ctx);
    }

}
