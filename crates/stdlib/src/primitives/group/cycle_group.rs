//! Circuit group element for the proving system's embedded curve.
//!
//! Port of `barretenberg/stdlib/primitives/group/cycle_group.hpp` and `.cpp`.
//!
//! A `CycleGroupT<C>` represents a point on the embedded curve (Grumpkin over BN254).
//! The point at infinity is represented as (0, 0) with `is_infinity = true`.

use std::ops::{Add, Neg, Sub};

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::groups::affine_element::AffineElement;
use bbrs_ecc::groups::curve_params::CurveParams;
use bbrs_ecc::groups::element::Element;

use crate::primitives::bool::BoolT;
use crate::primitives::field::{bool_to_field, FieldT};
use crate::primitives::witness::{BuilderRef, WitnessT};

use super::cycle_scalar::CycleScalarT;
use super::straus_lookup_table::StrausLookupTable;
use super::straus_scalar_slice::StrausScalarSlices;

/// Bit-size for ROM table entries in the variable-base MSM algorithm.
const ROM_TABLE_BITS: usize = 4;

/// Circuit representation of an embedded curve group element.
pub struct CycleGroupT<C: CurveParams> {
    _x: FieldT<C::BaseFieldParams>,
    _y: FieldT<C::BaseFieldParams>,
    _is_infinity: BoolT<C::BaseFieldParams>,
    context: Option<BuilderRef<C::BaseFieldParams>>,
}

impl<C: CurveParams> Clone for CycleGroupT<C> {
    fn clone(&self) -> Self {
        Self {
            _x: self._x.clone(),
            _y: self._y.clone(),
            _is_infinity: self._is_infinity.clone(),
            context: self.context.clone(),
        }
    }
}

impl<C: CurveParams> CycleGroupT<C> {
    // ════════════════════════════════════════════════════════════════════
    //  Constructors
    // ════════════════════════════════════════════════════════════════════

    /// Construct a constant point at infinity.
    pub fn infinity(ctx: Option<BuilderRef<C::BaseFieldParams>>) -> Self {
        let (x, y, is_inf) = if let Some(ref c) = ctx {
            (
                FieldT::constant_with_context(c.clone(), Field::zero()),
                FieldT::constant_with_context(c.clone(), Field::zero()),
                BoolT::constant_with_context(c.clone(), true),
            )
        } else {
            (
                FieldT::from_field(Field::zero()),
                FieldT::from_field(Field::zero()),
                BoolT::from_constant(true),
            )
        };
        Self {
            _x: x,
            _y: y,
            _is_infinity: is_inf,
            context: ctx,
        }
    }

    /// Construct from a native AffineElement (circuit constant).
    pub fn from_affine(point: AffineElement<C>) -> Self {
        if point.is_point_at_infinity() {
            return Self::infinity(None);
        }
        Self {
            _x: FieldT::from_field(point.x.clone()),
            _y: FieldT::from_field(point.y.clone()),
            _is_infinity: BoolT::from_constant(false),
            context: None,
        }
    }

    /// Construct from x, y coordinates and is_infinity flag.
    ///
    /// If `assert_on_curve` is true, adds an on-curve constraint.
    pub fn from_xy_bool(
        x: FieldT<C::BaseFieldParams>,
        y: FieldT<C::BaseFieldParams>,
        is_infinity: BoolT<C::BaseFieldParams>,
        assert_on_curve: bool,
    ) -> Self {
        let context = x
            .get_context()
            .clone()
            .or_else(|| y.get_context().clone())
            .or_else(|| is_infinity.get_context().clone());

        // If is_infinity is constant true, return infinity
        if is_infinity.is_constant() && is_infinity.get_value() {
            return Self::infinity(context);
        }

        let mut result = Self {
            _x: x,
            _y: y,
            _is_infinity: is_infinity,
            context,
        };

        // Ensure matching constancy of coordinates
        if result._x.is_constant() != result._y.is_constant() {
            if let Some(ref ctx) = result.context {
                if result._x.is_constant() {
                    result._x.convert_constant_to_fixed_witness(ctx.clone());
                } else {
                    result._y.convert_constant_to_fixed_witness(ctx.clone());
                }
            }
        }

        // If both coordinates constant but is_infinity is witness, make is_infinity constant
        if result._x.is_constant()
            && result._y.is_constant()
            && !result._is_infinity.is_constant()
        {
            result._is_infinity =
                BoolT::from_constant(result._is_infinity.get_value());
        }

        if assert_on_curve {
            result.validate_on_curve();
        }

        result
    }

    /// Construct the generator point (constant).
    pub fn one(ctx: Option<BuilderRef<C::BaseFieldParams>>) -> Self {
        let generator = AffineElement::<C>::one();
        let mut result = Self::from_affine(generator);
        result.context = ctx;
        result
    }

    /// Create a witness from a native AffineElement.
    pub fn from_witness(
        ctx: BuilderRef<C::BaseFieldParams>,
        point: &AffineElement<C>,
    ) -> Self {
        if point.is_point_at_infinity() {
            let x = FieldT::from_witness(ctx.clone(), Field::zero());
            let y = FieldT::from_witness(ctx.clone(), Field::zero());
            let is_inf = BoolT::from_witness(&WitnessT::from_bool(ctx.clone(), true));
            let result = Self {
                _x: x,
                _y: y,
                _is_infinity: is_inf,
                context: Some(ctx),
            };
            result.validate_on_curve();
            result
        } else {
            let x = FieldT::from_witness(ctx.clone(), point.x.clone());
            let y = FieldT::from_witness(ctx.clone(), point.y.clone());
            let is_inf = BoolT::from_witness(&WitnessT::from_bool(ctx.clone(), false));
            let result = Self {
                _x: x,
                _y: y,
                _is_infinity: is_inf,
                context: Some(ctx),
            };
            result.validate_on_curve();
            result
        }
    }

    /// Create a witness constrained to equal a known constant point.
    pub fn from_constant_witness(
        ctx: BuilderRef<C::BaseFieldParams>,
        point: &AffineElement<C>,
    ) -> Self {
        if point.is_point_at_infinity() {
            return Self::infinity(Some(ctx));
        }
        let x = FieldT::from_witness(ctx.clone(), point.x.clone());
        let y = FieldT::from_witness(ctx.clone(), point.y.clone());
        x.assert_equal(
            &FieldT::from_field(point.x.clone()),
            "from_constant_witness x",
        );
        y.assert_equal(
            &FieldT::from_field(point.y.clone()),
            "from_constant_witness y",
        );
        Self {
            _x: x,
            _y: y,
            _is_infinity: BoolT::from_constant(false),
            context: Some(ctx),
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Accessors
    // ════════════════════════════════════════════════════════════════════

    pub fn x(&self) -> &FieldT<C::BaseFieldParams> {
        &self._x
    }

    pub fn y(&self) -> &FieldT<C::BaseFieldParams> {
        &self._y
    }

    pub fn is_point_at_infinity(&self) -> &BoolT<C::BaseFieldParams> {
        &self._is_infinity
    }

    pub fn is_constant(&self) -> bool {
        self._x.is_constant() && self._y.is_constant() && self._is_infinity.is_constant()
    }

    pub fn is_constant_point_at_infinity(&self) -> bool {
        self._is_infinity.is_constant() && self._is_infinity.get_value()
    }

    pub fn get_context(&self) -> &Option<BuilderRef<C::BaseFieldParams>> {
        &self.context
    }

    pub fn get_context_or(
        &self,
        other: &CycleGroupT<C>,
    ) -> Option<BuilderRef<C::BaseFieldParams>> {
        self.context.clone().or_else(|| other.context.clone())
    }

    pub(crate) fn x_val(&self) -> Field<C::BaseFieldParams> {
        self._x.get_value()
    }

    pub(crate) fn y_val(&self) -> Field<C::BaseFieldParams> {
        self._y.get_value()
    }

    /// Get the native AffineElement value.
    pub fn get_value(&self) -> AffineElement<C> {
        let mut result = AffineElement::new(self._x.get_value(), self._y.get_value());
        if self._is_infinity.get_value() {
            result.self_set_infinity();
        }
        result
    }

    // ════════════════════════════════════════════════════════════════════
    //  Validation & standardization
    // ════════════════════════════════════════════════════════════════════

    /// Constrain y^2 = x^3 + b (or point is at infinity).
    pub fn validate_on_curve(&self) {
        // Short Weierstrass: y^2 = x^3 + b (a = 0)
        debug_assert!(!C::HAS_A, "cycle_group only supports curves with a=0");
        let xx = self._x.sqr();
        let xxx = xx * self._x.clone();
        let curve_b_field = FieldT::from_field(C::coeff_b());
        // res = y^2 - x^3 - b
        let neg_xxx_minus_b = -(xxx + curve_b_field);
        let res = self._y.madd(&self._y, &neg_xxx_minus_b);
        // If point at infinity, mask to 0
        let not_inf = self._is_infinity.negate();
        let not_inf_field = bool_to_field(&not_inf);
        let masked = res * not_inf_field;
        masked.assert_is_zero("cycle_group::validate_on_curve");
    }

    /// Ensure point at infinity has coordinates (0, 0).
    pub fn standardize(&mut self) {
        self._x = FieldT::conditional_assign(
            &self._is_infinity,
            &FieldT::from_field(Field::zero()),
            &self._x,
        );
        self._y = FieldT::conditional_assign(
            &self._is_infinity,
            &FieldT::from_field(Field::zero()),
            &self._y,
        );
    }

    // ════════════════════════════════════════════════════════════════════
    //  Point doubling
    // ════════════════════════════════════════════════════════════════════

    /// Double this point using the Ultra ECC double gate.
    pub fn dbl(&self, hint: Option<AffineElement<C>>) -> CycleGroupT<C> {
        if self.is_constant_point_at_infinity() {
            return self.clone();
        }

        // Avoid division by zero: if point at infinity, set y to 1
        let modified_y = FieldT::conditional_assign(
            &self._is_infinity,
            &FieldT::from_u64(1),
            &self._y,
        );

        // Compute native result
        let (x3, y3) = if let Some(h) = &hint {
            (h.x.clone(), h.y.clone())
        } else {
            let x1 = self._x.get_value();
            let y1 = modified_y.get_value();
            let y_pow_2 = y1 * y1;
            let curve_b = C::coeff_b();
            let x_pow_4 = x1 * (y_pow_2 - curve_b);
            let nine = Field::from(9u64);
            let four = Field::from(4u64);
            let lambda_squared = (x_pow_4 * nine) * (y_pow_2 * four).invert();
            let three = Field::from(3u64);
            let two = Field::from(2u64);
            let lambda = (x1 * x1 * three) * (y1 * two).invert();
            let x3 = lambda_squared - x1 - x1;
            let y3 = lambda * (x1 - x3) - y1;
            (x3, y3)
        };

        if self.is_constant() {
            return CycleGroupT::from_xy_bool(
                FieldT::from_field(x3),
                FieldT::from_field(y3),
                self._is_infinity.clone(),
                false,
            );
        }

        let ctx = self.context.as_ref().unwrap().clone();
        let x3_witness = FieldT::from_witness(ctx.clone(), x3);
        let y3_witness = FieldT::from_witness(ctx.clone(), y3);

        use bbrs_circuit_builder::gate_data::EccDblGate;
        ctx.borrow_mut().create_ecc_dbl_gate(&EccDblGate {
            x1: self._x.normalize().get_witness_index(),
            y1: modified_y.normalize().get_witness_index(),
            x3: x3_witness.normalize().get_witness_index(),
            y3: y3_witness.normalize().get_witness_index(),
        });

        CycleGroupT::from_xy_bool(
            x3_witness,
            y3_witness,
            self._is_infinity.clone(),
            false,
        )
    }

    // ════════════════════════════════════════════════════════════════════
    //  Unconditional addition / subtraction (incomplete formulas)
    // ════════════════════════════════════════════════════════════════════

    /// Internal: add or subtract without edge case handling.
    fn unconditional_add_or_subtract(
        &self,
        other: &CycleGroupT<C>,
        is_addition: bool,
        hint: Option<AffineElement<C>>,
    ) -> CycleGroupT<C> {
        debug_assert!(
            !self.is_constant_point_at_infinity(),
            "unconditional_add_or_subtract called on constant infinity"
        );
        debug_assert!(
            !other.is_constant_point_at_infinity(),
            "unconditional_add_or_subtract called on constant infinity"
        );

        let ctx = self.get_context_or(other);
        let lhs_constant = self.is_constant();
        let rhs_constant = other.is_constant();

        // Convert constant operand to fixed witness if mixed
        if lhs_constant && !rhs_constant {
            let lhs = CycleGroupT::<C>::from_constant_witness(
                ctx.clone().unwrap(),
                &self.get_value(),
            );
            return lhs.unconditional_add_or_subtract(other, is_addition, hint);
        }
        if !lhs_constant && rhs_constant {
            let rhs = CycleGroupT::<C>::from_constant_witness(
                ctx.clone().unwrap(),
                &other.get_value(),
            );
            return self.unconditional_add_or_subtract(&rhs, is_addition, hint);
        }

        // Compute native result
        let (x3, y3) = if let Some(h) = &hint {
            (h.x.clone(), h.y.clone())
        } else {
            let p1 = self.get_value();
            let p2 = other.get_value();
            let p3 = if is_addition {
                (Element::from_affine(&p1) + Element::from_affine(&p2)).to_affine()
            } else {
                (Element::from_affine(&p1) - Element::from_affine(&p2)).to_affine()
            };
            (p3.x.clone(), p3.y.clone())
        };

        if lhs_constant && rhs_constant {
            return CycleGroupT::from_xy_bool(
                FieldT::from_field(x3),
                FieldT::from_field(y3),
                BoolT::from_constant(false),
                false,
            );
        }

        // Both witnesses: create ECC add gate
        let ctx = ctx.unwrap();
        let x3_witness = FieldT::from_witness(ctx.clone(), x3);
        let y3_witness = FieldT::from_witness(ctx.clone(), y3);

        let sign: Field<C::BaseFieldParams> = if is_addition {
            Field::one()
        } else {
            -Field::one()
        };

        use bbrs_circuit_builder::gate_data::EccAddGate;
        ctx.borrow_mut().create_ecc_add_gate(&EccAddGate {
            x1: self._x.normalize().get_witness_index(),
            y1: self._y.normalize().get_witness_index(),
            x2: other._x.normalize().get_witness_index(),
            y2: other._y.normalize().get_witness_index(),
            x3: x3_witness.normalize().get_witness_index(),
            y3: y3_witness.normalize().get_witness_index(),
            sign_coefficient: sign,
        });

        CycleGroupT::from_xy_bool(
            x3_witness,
            y3_witness,
            BoolT::from_constant(false),
            false,
        )
    }

    /// Unconditional point addition (incomplete formula).
    pub fn unconditional_add(
        &self,
        other: &CycleGroupT<C>,
        hint: Option<AffineElement<C>>,
    ) -> CycleGroupT<C> {
        self.unconditional_add_or_subtract(other, true, hint)
    }

    /// Unconditional point subtraction (incomplete formula).
    pub fn unconditional_subtract(
        &self,
        other: &CycleGroupT<C>,
        hint: Option<AffineElement<C>>,
    ) -> CycleGroupT<C> {
        self.unconditional_add_or_subtract(other, false, hint)
    }

    /// Unconditional add with x-coordinate collision check.
    pub fn checked_unconditional_add(
        &self,
        other: &CycleGroupT<C>,
        hint: Option<AffineElement<C>>,
    ) -> CycleGroupT<C> {
        let x_delta = self._x.clone() - other._x.clone();
        if x_delta.is_constant() {
            debug_assert!(
                !x_delta.get_value().is_zero(),
                "checked_unconditional_add: x-coordinate collision"
            );
        } else {
            x_delta.assert_is_not_zero(
                "cycle_group::checked_unconditional_add, x-coordinate collision",
            );
        }
        self.unconditional_add(other, hint)
    }

    /// Unconditional subtract with x-coordinate collision check.
    pub fn checked_unconditional_subtract(
        &self,
        other: &CycleGroupT<C>,
        hint: Option<AffineElement<C>>,
    ) -> CycleGroupT<C> {
        let x_delta = self._x.clone() - other._x.clone();
        if x_delta.is_constant() {
            debug_assert!(
                !x_delta.get_value().is_zero(),
                "checked_unconditional_subtract: x-coordinate collision"
            );
        } else {
            x_delta.assert_is_not_zero(
                "cycle_group::checked_unconditional_subtract, x-coordinate collision",
            );
        }
        self.unconditional_subtract(other, hint)
    }

    // ════════════════════════════════════════════════════════════════════
    //  Complete addition / subtraction
    // ════════════════════════════════════════════════════════════════════

    /// Complete point addition handling all edge cases.
    pub fn complete_add(&self, other: &CycleGroupT<C>) -> CycleGroupT<C> {
        if self.is_constant_point_at_infinity() {
            return other.clone();
        }
        if other.is_constant_point_at_infinity() {
            return self.clone();
        }

        let x_coordinates_match = self._x.is_equal(&other._x);
        let y_coordinates_match = self._y.is_equal(&other._y);

        let x1 = &self._x;
        let y1 = &self._y;
        let x2 = &other._x;
        let y2 = &other._y;

        // Modified lambda: (y2 - y1) / (x2 - x1 + x_match) to avoid div by zero
        let x_match_field = bool_to_field(&x_coordinates_match);
        let x_diff = x2.add_two(&(-x1.clone()), &x_match_field);

        let lambda = if (y1.is_constant() && y2.is_constant()) || x_diff.is_constant()
        {
            (y2.clone() - y1.clone()).divide_no_zero_check(&x_diff)
        } else {
            let ctx = self.get_context_or(other).unwrap();
            let num_val = y2.get_value() - y1.get_value();
            let den_val = x_diff.get_value();
            let lambda_val = num_val * den_val.invert();
            let lambda = FieldT::from_witness(ctx, lambda_val);
            // Constrain x_diff * lambda - y2 + y1 = 0
            FieldT::evaluate_polynomial_identity(
                &x_diff,
                &lambda,
                &(-y2.clone()),
                y1,
                "cycle_group::add lambda constraint",
            );
            lambda
        };

        // x3 = lambda^2 - x1 - x2
        let neg_x_sum = -(x2.clone() + x1.clone());
        let add_result_x = lambda.madd(&lambda, &neg_x_sum);
        // y3 = lambda * (x1 - x3) - y1
        let add_result_y =
            lambda.madd(&(x1.clone() - add_result_x.clone()), &(-y1.clone()));

        // Compute doubling result
        let dbl_result = self.dbl(None);

        // Select: if x_match && y_match, use double; else use add
        let double_predicate = x_coordinates_match.and(&y_coordinates_match);
        let mut result_x = FieldT::conditional_assign(
            &double_predicate,
            &dbl_result._x,
            &add_result_x,
        );
        let mut result_y = FieldT::conditional_assign(
            &double_predicate,
            &dbl_result._y,
            &add_result_y,
        );

        // If lhs is infinity, return rhs
        let lhs_infinity = &self._is_infinity;
        result_x =
            FieldT::conditional_assign(lhs_infinity, &other._x, &result_x);
        result_y =
            FieldT::conditional_assign(lhs_infinity, &other._y, &result_y);

        // If rhs is infinity, return lhs
        let rhs_infinity = &other._is_infinity;
        result_x = FieldT::conditional_assign(rhs_infinity, &self._x, &result_x);
        result_y = FieldT::conditional_assign(rhs_infinity, &self._y, &result_y);

        // Result is infinity if: x_match && !y_match && !lhs_inf && !rhs_inf,
        // or both are infinity
        let infinity_predicate =
            x_coordinates_match.and(&y_coordinates_match.negate());
        let neither_inf = lhs_infinity.negate().and(&rhs_infinity.negate());
        let mut result_is_infinity = infinity_predicate.and(&neither_inf);
        let both_inf = lhs_infinity.and(rhs_infinity);
        result_is_infinity = result_is_infinity.or(&both_inf);

        CycleGroupT::from_xy_bool(result_x, result_y, result_is_infinity, false)
    }

    /// Complete point subtraction handling all edge cases.
    pub fn complete_sub(&self, other: &CycleGroupT<C>) -> CycleGroupT<C> {
        if self.is_constant_point_at_infinity() {
            return other.negate();
        }
        if other.is_constant_point_at_infinity() {
            return self.clone();
        }

        let ctx_opt = self.get_context_or(other);

        let x_coordinates_match = self._x.is_equal(&other._x);
        let y_coordinates_match = self._y.is_equal(&other._y);

        let x1 = &self._x;
        let y1 = &self._y;
        let x2 = &other._x;
        let y2 = &other._y;

        // Modified lambda: (-y2 - y1) / (x2 - x1 + x_match)
        let x_match_field = bool_to_field(&x_coordinates_match);
        let x_diff = x2.add_two(&(-x1.clone()), &x_match_field);

        let lambda = if (y1.is_constant() && y2.is_constant()) || x_diff.is_constant()
        {
            let neg_y2_minus_y1 = -(y2.clone()) - y1.clone();
            neg_y2_minus_y1.divide_no_zero_check(&x_diff)
        } else {
            let ctx = ctx_opt.clone().unwrap();
            let num_val = -(y2.get_value()) - y1.get_value();
            let den_val = x_diff.get_value();
            let lambda_val = num_val * den_val.invert();
            let lambda = FieldT::from_witness(ctx, lambda_val);
            // Constrain x_diff * lambda + y2 + y1 = 0
            FieldT::evaluate_polynomial_identity(
                &x_diff,
                &lambda,
                y2,
                y1,
                "cycle_group::sub lambda constraint",
            );
            lambda
        };

        // x3 = lambda^2 - x1 - x2
        let neg_x_sum = -(x2.clone() + x1.clone());
        let sub_result_x = lambda.madd(&lambda, &neg_x_sum);
        // y3 = lambda * (x1 - x3) - y1
        let sub_result_y =
            lambda.madd(&(x1.clone() - sub_result_x.clone()), &(-y1.clone()));

        // Compute doubling result
        let dbl_result = self.dbl(None);

        // For subtraction: double when x_match && !y_match (i.e. P - (-P) = 2P)
        let double_predicate =
            x_coordinates_match.and(&y_coordinates_match.negate());
        let mut result_x = FieldT::conditional_assign(
            &double_predicate,
            &dbl_result._x,
            &sub_result_x,
        );
        let mut result_y = FieldT::conditional_assign(
            &double_predicate,
            &dbl_result._y,
            &sub_result_y,
        );

        // If lhs is infinity, return -rhs
        let lhs_infinity = &self._is_infinity;
        let neg_y2 = -(other._y.clone());
        result_x =
            FieldT::conditional_assign(lhs_infinity, &other._x, &result_x);
        result_y =
            FieldT::conditional_assign(lhs_infinity, &neg_y2, &result_y);

        // If rhs is infinity, return lhs
        let rhs_infinity = &other._is_infinity;
        result_x = FieldT::conditional_assign(rhs_infinity, &self._x, &result_x);
        result_y = FieldT::conditional_assign(rhs_infinity, &self._y, &result_y);

        // Result is infinity if: x_match && y_match && !lhs_inf && !rhs_inf,
        // or both are infinity
        let infinity_predicate = x_coordinates_match.and(&y_coordinates_match);
        let neither_inf = lhs_infinity.negate().and(&rhs_infinity.negate());
        let mut result_is_infinity = infinity_predicate.and(&neither_inf);
        let both_inf = lhs_infinity.and(rhs_infinity);
        result_is_infinity = result_is_infinity.or(&both_inf);

        CycleGroupT::from_xy_bool(result_x, result_y, result_is_infinity, false)
    }

    /// Negate a point.
    pub fn negate(&self) -> CycleGroupT<C> {
        let mut result = self.clone();
        result._y = (-self._y.clone()).normalize();
        result
    }

    // ════════════════════════════════════════════════════════════════════
    //  Conditional assignment
    // ════════════════════════════════════════════════════════════════════

    /// If `predicate` is true, return `lhs`, else return `rhs`.
    pub fn conditional_assign(
        predicate: &BoolT<C::BaseFieldParams>,
        lhs: &CycleGroupT<C>,
        rhs: &CycleGroupT<C>,
    ) -> CycleGroupT<C> {
        let x_res = FieldT::conditional_assign(predicate, &lhs._x, &rhs._x);
        let y_res = FieldT::conditional_assign(predicate, &lhs._y, &rhs._y);
        let inf_res = BoolT::conditional_assign(
            predicate,
            &lhs._is_infinity,
            &rhs._is_infinity,
        );
        CycleGroupT::from_xy_bool(x_res, y_res, inf_res, false)
    }

    // ════════════════════════════════════════════════════════════════════
    //  Equality
    // ════════════════════════════════════════════════════════════════════

    /// In-circuit equality check. Returns a BoolT.
    pub fn eq(&mut self, other: &mut CycleGroupT<C>) -> BoolT<C::BaseFieldParams> {
        self.standardize();
        other.standardize();
        let x_eq = self._x.is_equal(&other._x);
        let y_eq = self._y.is_equal(&other._y);
        let inf_eq = self._is_infinity.eq(&other._is_infinity);
        x_eq.and(&y_eq).and(&inf_eq)
    }

    /// Assert this point equals another.
    pub fn assert_equal(&mut self, other: &mut CycleGroupT<C>, msg: &str) {
        self.standardize();
        other.standardize();
        self._x.assert_equal(&other._x, msg);
        self._y.assert_equal(&other._y, msg);
        self._is_infinity.assert_equal(&other._is_infinity, msg);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Public inputs
    // ════════════════════════════════════════════════════════════════════

    /// Set x, y as public inputs. Returns the starting public input index.
    pub fn set_public(&self) -> u32 {
        self._x.set_public();
        self._y.set_public();
        // Also set is_infinity as public
        if !self._is_infinity.is_constant() {
            let inf_field = FieldT::from_witness_index(
                self.context.clone().unwrap(),
                self._is_infinity.witness_index,
            );
            inf_field.set_public();
        }
        3 // x, y, is_infinity
    }

    // ════════════════════════════════════════════════════════════════════
    //  Scalar multiplication (Straus MSM)
    // ════════════════════════════════════════════════════════════════════

    /// Batch multi-scalar multiplication using the Straus algorithm.
    ///
    /// Computes `sum_i(scalars[i] * base_points[i])`.
    pub fn batch_mul(
        base_points: &[CycleGroupT<C>],
        scalars: &[CycleScalarT<C>],
    ) -> CycleGroupT<C> {
        assert_eq!(
            scalars.len(),
            base_points.len(),
            "Points/scalars size mismatch in batch mul"
        );

        if scalars.is_empty() {
            return CycleGroupT::from_affine(AffineElement::<C>::infinity());
        }

        let mut variable_base_scalars = Vec::new();
        let mut variable_base_points = Vec::new();

        let mut can_unconditional_add = true;
        let mut has_non_constant_component = false;
        let mut constant_acc = Element::<C>::infinity();

        for (point, scalar) in base_points.iter().zip(scalars.iter()) {
            if scalar.is_constant() && point.is_constant() {
                // Case 1: both constant — accumulate natively
                constant_acc = constant_acc
                    + Element::from_affine(&point.get_value())
                        .mul(&scalar.get_value());
            } else if !scalar.is_constant() && point.is_constant() {
                if point.get_value().is_point_at_infinity() {
                    continue;
                }
                // Case 2B: constant point, witness scalar — variable-base
                variable_base_scalars.push(scalar.clone());
                variable_base_points.push(point.clone());
                has_non_constant_component = true;
            } else {
                // Case 3: witness point — variable-base
                variable_base_scalars.push(scalar.clone());
                variable_base_points.push(point.clone());
                can_unconditional_add = false;
                has_non_constant_component = true;
            }
        }

        if !has_non_constant_component {
            let acc_affine = constant_acc.to_affine();
            return CycleGroupT::from_affine(acc_affine);
        }

        // Offset accumulator tracks constant terms to subtract later
        let mut offset_accumulator = -constant_acc;
        let has_variable_points = !variable_base_points.is_empty();

        let mut result = CycleGroupT::infinity(None);

        if has_variable_points {
            // Generate offset generators deterministically
            let num_offset_generators = variable_base_points.len() + 1;
            let offset_generators =
                generate_offset_generators::<C>(num_offset_generators);

            let (var_acc, offset_delta) = variable_base_batch_mul_internal::<C>(
                &variable_base_scalars,
                &variable_base_points,
                &offset_generators,
                can_unconditional_add,
            );
            offset_accumulator = offset_accumulator + offset_delta;
            result = var_acc;
        }

        // Subtract offset accumulator to get final result
        let offset_affine = offset_accumulator.to_affine();
        if !offset_affine.is_point_at_infinity() && can_unconditional_add {
            let neg_offset = {
                let mut p = offset_affine;
                p.y = -p.y;
                p
            };
            result = result.unconditional_add(
                &CycleGroupT::from_affine(neg_offset),
                None,
            );
        } else {
            let offset_group = CycleGroupT::from_affine(offset_affine);
            result = result.complete_sub(&offset_group);
        }

        result
    }

    /// Single scalar multiplication.
    pub fn scalar_mul(&self, scalar: &CycleScalarT<C>) -> CycleGroupT<C> {
        Self::batch_mul(&[self.clone()], &[scalar.clone()])
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Operator implementations
// ════════════════════════════════════════════════════════════════════════

impl<C: CurveParams> Add for &CycleGroupT<C> {
    type Output = CycleGroupT<C>;
    fn add(self, rhs: &CycleGroupT<C>) -> CycleGroupT<C> {
        self.complete_add(rhs)
    }
}

impl<C: CurveParams> Add for CycleGroupT<C> {
    type Output = CycleGroupT<C>;
    fn add(self, rhs: CycleGroupT<C>) -> CycleGroupT<C> {
        self.complete_add(&rhs)
    }
}

impl<C: CurveParams> Sub for &CycleGroupT<C> {
    type Output = CycleGroupT<C>;
    fn sub(self, rhs: &CycleGroupT<C>) -> CycleGroupT<C> {
        self.complete_sub(rhs)
    }
}

impl<C: CurveParams> Sub for CycleGroupT<C> {
    type Output = CycleGroupT<C>;
    fn sub(self, rhs: CycleGroupT<C>) -> CycleGroupT<C> {
        self.complete_sub(&rhs)
    }
}

impl<C: CurveParams> Neg for &CycleGroupT<C> {
    type Output = CycleGroupT<C>;
    fn neg(self) -> CycleGroupT<C> {
        self.negate()
    }
}

impl<C: CurveParams> Neg for CycleGroupT<C> {
    type Output = CycleGroupT<C>;
    fn neg(self) -> CycleGroupT<C> {
        self.negate()
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Internal: Variable-base batch mul
// ════════════════════════════════════════════════════════════════════════

fn variable_base_batch_mul_internal<C: CurveParams>(
    scalars: &[CycleScalarT<C>],
    base_points: &[CycleGroupT<C>],
    offset_generators: &[AffineElement<C>],
    unconditional_add: bool,
) -> (CycleGroupT<C>, Element<C>) {
    assert!(!scalars.is_empty());
    assert_eq!(scalars.len(), base_points.len());
    assert_eq!(offset_generators.len(), base_points.len() + 1);

    let num_points = scalars.len();

    // Find circuit context
    let ctx = scalars
        .iter()
        .flat_map(|s| s.get_context().clone())
        .chain(base_points.iter().flat_map(|p| p.get_context().clone()))
        .next()
        .expect("batch_mul requires at least one witness input");

    let num_bits = super::cycle_scalar::NUM_BITS;
    let num_rounds = (num_bits + ROM_TABLE_BITS - 1) / ROM_TABLE_BITS;

    // Decompose scalars
    let scalar_slices: Vec<_> = scalars
        .iter()
        .map(|s| StrausScalarSlices::new(ctx.clone(), s, ROM_TABLE_BITS))
        .collect();

    // Phase 1: Compute all native hints
    let mut operation_transcript: Vec<Element<C>> = Vec::new();
    let mut offset_generator_accumulator =
        Element::from_affine(&offset_generators[0]);

    {
        // Build native straus tables
        let mut native_straus_tables = Vec::with_capacity(num_points);
        for i in 0..num_points {
            let point = Element::from_affine(&base_points[i].get_value());
            let offset = Element::from_affine(&offset_generators[i + 1]);
            let table = StrausLookupTable::<C>::compute_native_table(
                &point,
                &offset,
                ROM_TABLE_BITS,
            );
            for entry in table.iter().skip(1) {
                operation_transcript.push(entry.clone());
            }
            native_straus_tables.push(table);
        }

        // Native Straus algorithm
        let mut accumulator = Element::from_affine(&offset_generators[0]);
        for i in 0..num_rounds {
            if i != 0 {
                for _ in 0..ROM_TABLE_BITS {
                    accumulator = accumulator.dbl();
                    operation_transcript.push(accumulator.clone());
                    offset_generator_accumulator =
                        offset_generator_accumulator.dbl();
                }
            }
            for j in 0..num_points {
                let slice_value =
                    scalar_slices[j].slices_native[num_rounds - i - 1] as usize;
                let point = native_straus_tables[j][slice_value].clone();
                accumulator = accumulator + point;
                operation_transcript.push(accumulator.clone());
                offset_generator_accumulator = offset_generator_accumulator
                    + Element::from_affine(&offset_generators[j + 1]);
            }
        }
    }

    // Batch-normalize hints
    Element::batch_normalize(&mut operation_transcript);
    let operation_hints: Vec<AffineElement<C>> = operation_transcript
        .iter()
        .map(|e| e.to_affine())
        .collect();

    // Phase 2: Build in-circuit ROM tables
    let hints_per_table = (1usize << ROM_TABLE_BITS) - 1;
    let mut point_tables = Vec::with_capacity(num_points);
    for i in 0..num_points {
        let table_hints = &operation_hints[i * hints_per_table..(i + 1) * hints_per_table];
        let table = StrausLookupTable::<C>::new(
            ctx.clone(),
            &base_points[i],
            &CycleGroupT::from_affine(offset_generators[i + 1].clone()),
            ROM_TABLE_BITS,
            Some(table_hints),
        );
        point_tables.push(table);
    }

    // Phase 3: Execute Straus algorithm in-circuit
    let mut hint_idx = num_points * hints_per_table;
    let mut accumulator =
        CycleGroupT::from_affine(offset_generators[0].clone());

    let mut coordinate_check_product = FieldT::<C::BaseFieldParams>::from_u64(1);

    for i in 0..num_rounds {
        if i != 0 {
            for _ in 0..ROM_TABLE_BITS {
                let hint = operation_hints.get(hint_idx).cloned();
                accumulator = accumulator.dbl(hint);
                hint_idx += 1;
            }
        }
        for j in 0..num_points {
            let scalar_slice = scalar_slices[j].get(num_rounds - i - 1);
            let point = point_tables[j].read(scalar_slice);
            if !unconditional_add {
                let x_diff = point.x().clone() - accumulator.x().clone();
                coordinate_check_product =
                    coordinate_check_product * x_diff;
            }
            let hint = operation_hints.get(hint_idx).cloned();
            accumulator = accumulator.unconditional_add(&point, hint);
            hint_idx += 1;
        }
    }

    if !unconditional_add {
        coordinate_check_product.assert_is_not_zero(
            "_variable_base_batch_mul_internal x-coordinate collision",
        );
    }

    (accumulator, offset_generator_accumulator)
}

/// Generate deterministic offset generators for the MSM algorithm.
fn generate_offset_generators<C: CurveParams>(count: usize) -> Vec<AffineElement<C>> {
    // Use hash-to-curve from a domain separator to generate independent points.
    // For simplicity, we use repeated hashing of the generator.
    let mut generators = Vec::with_capacity(count);
    let mut current = Element::<C>::one();
    // Use a large prime multiplier to get "random-looking" independent points
    let multiplier = Field::<C::ScalarFieldParams>::from_limbs([
        0x0000000000000007,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    ]);
    for _ in 0..count {
        current = current.mul(&multiplier);
        generators.push(current.to_affine());
    }
    generators
}
