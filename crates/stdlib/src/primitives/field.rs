//! Circuit field element type.
//!
//! Port of `barretenberg/stdlib/primitives/field/field.hpp` and `field.cpp`.
//!
//! A `FieldT<P>` represents a circuit-level field element with lazy normalization.
//! Each element tracks a witness index plus multiplicative and additive scaling constants:
//!
//! ```text
//! value = witness[witness_index] * multiplicative_constant + additive_constant
//! ```
//!
//! Constants (with `witness_index == IS_CONSTANT`) store their value in `additive_constant`.
//! Operations between constants create no gates.

use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::rc::Rc;

use bbrs_circuit_builder::gate_data::{AddQuad, AddTriple, ArithmeticTriple, MulQuad};
use bbrs_circuit_builder::ultra_builder::{self, UltraCircuitBuilder};
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::witness::{BuilderRef, WitnessT, IS_CONSTANT};

/// Compute MSB position of a u256 represented as [u64; 4] limbs.
fn get_msb_u256(data: &[u64; 4]) -> u32 {
    for i in (0..4).rev() {
        if data[i] != 0 {
            return (i as u32) * 64 + (63 - data[i].leading_zeros());
        }
    }
    0
}

// ════════════════════════════════════════════════════════════════════════
//  Context validation
// ════════════════════════════════════════════════════════════════════════

/// Validate that two optional builder references point to the same builder.
/// Returns the non-None reference (or None if both are None).
fn validate_two_contexts<P: FieldParams>(
    a: &Option<BuilderRef<P>>,
    b: &Option<BuilderRef<P>>,
) -> Option<BuilderRef<P>> {
    match (a, b) {
        (None, None) => None,
        (Some(c), None) | (None, Some(c)) => Some(c.clone()),
        (Some(c1), Some(c2)) => {
            assert!(
                Rc::ptr_eq(c1, c2),
                "Field elements belong to different builders"
            );
            Some(c1.clone())
        }
    }
}

/// Validate that three optional builder references point to the same builder.
fn validate_three_contexts<P: FieldParams>(
    a: &Option<BuilderRef<P>>,
    b: &Option<BuilderRef<P>>,
    c: &Option<BuilderRef<P>>,
) -> Option<BuilderRef<P>> {
    let ab = validate_two_contexts(a, b);
    validate_two_contexts(&ab, c)
}

/// Validate four optional builder references.
fn validate_four_contexts<P: FieldParams>(
    a: &Option<BuilderRef<P>>,
    b: &Option<BuilderRef<P>>,
    c: &Option<BuilderRef<P>>,
    d: &Option<BuilderRef<P>>,
) -> Option<BuilderRef<P>> {
    let ab = validate_two_contexts(a, b);
    let cd = validate_two_contexts(c, d);
    validate_two_contexts(&ab, &cd)
}

// ════════════════════════════════════════════════════════════════════════
//  FieldT
// ════════════════════════════════════════════════════════════════════════

/// A circuit field element with lazy normalization.
///
/// Port of C++ `field_t<Builder>`.
pub struct FieldT<P: FieldParams> {
    pub(crate) context: Option<BuilderRef<P>>,
    pub(crate) additive_constant: Field<P>,
    pub(crate) multiplicative_constant: Field<P>,
    pub(crate) witness_index: u32,
}

impl<P: FieldParams> Clone for FieldT<P> {
    fn clone(&self) -> Self {
        Self {
            context: self.context.clone(),
            additive_constant: self.additive_constant,
            multiplicative_constant: self.multiplicative_constant,
            witness_index: self.witness_index,
        }
    }
}

impl<P: FieldParams> FieldT<P> {
    // ════════════════════════════════════════════════════════════════════
    //  Constructors
    // ════════════════════════════════════════════════════════════════════

    /// Create a constant zero field element (no builder context).
    pub fn default() -> Self {
        Self {
            context: None,
            additive_constant: Field::zero(),
            multiplicative_constant: Field::one(),
            witness_index: IS_CONSTANT,
        }
    }

    /// Create a constant zero field element with builder context.
    pub fn with_context(ctx: BuilderRef<P>) -> Self {
        Self {
            context: Some(ctx),
            additive_constant: Field::zero(),
            multiplicative_constant: Field::one(),
            witness_index: IS_CONSTANT,
        }
    }

    /// Create a constant field element with builder context and value.
    pub fn constant_with_context(ctx: BuilderRef<P>, value: Field<P>) -> Self {
        Self {
            context: Some(ctx),
            additive_constant: value,
            multiplicative_constant: Field::one(),
            witness_index: IS_CONSTANT,
        }
    }

    /// Create a constant from a native field element (no context).
    pub fn from_field(value: Field<P>) -> Self {
        Self {
            context: None,
            additive_constant: value,
            multiplicative_constant: Field::one(),
            witness_index: IS_CONSTANT,
        }
    }

    /// Create a constant from a u64 value.
    pub fn from_u64(value: u64) -> Self {
        Self::from_field(Field::from(value))
    }

    /// Create a field element from a witness.
    pub fn from_witness_t(witness: &WitnessT<P>) -> Self {
        Self {
            context: Some(witness.context.clone()),
            additive_constant: Field::zero(),
            multiplicative_constant: Field::one(),
            witness_index: witness.witness_index,
        }
    }

    /// Create a field element by witness index.
    pub fn from_witness_index(ctx: BuilderRef<P>, witness_index: u32) -> Self {
        Self {
            context: Some(ctx),
            additive_constant: Field::zero(),
            multiplicative_constant: Field::one(),
            witness_index,
        }
    }

    /// Create a new witness from a native field value. The witness is "free"
    /// (not yet constrained).
    pub fn from_witness(ctx: BuilderRef<P>, value: Field<P>) -> Self {
        Self::from_witness_t(&WitnessT::new(ctx, value))
    }

    /// Create a new witness that is a copy of another field element,
    /// constrained to be equal.
    pub fn copy_as_new_witness(ctx: BuilderRef<P>, other: &FieldT<P>) -> Self {
        let result = Self::from_witness(ctx, other.get_value());
        result.assert_equal(other, "field_t::copy_as_new_witness, assert_equal");
        result
    }

    // ════════════════════════════════════════════════════════════════════
    //  Value access
    // ════════════════════════════════════════════════════════════════════

    /// Get the actual value represented by this field element.
    ///
    /// For constants, returns `additive_constant`.
    /// For witnesses, returns `witness * multiplicative_constant + additive_constant`.
    pub fn get_value(&self) -> Field<P> {
        if !self.is_constant() {
            let ctx = self.context.as_ref().expect("non-constant must have context");
            let witness = ctx.borrow().base.get_variable(self.witness_index);
            self.multiplicative_constant * witness + self.additive_constant
        } else {
            debug_assert!(self.multiplicative_constant == Field::one());
            self.additive_constant
        }
    }

    /// Get the builder context (if any).
    pub fn get_context(&self) -> &Option<BuilderRef<P>> {
        &self.context
    }

    /// Check if this field element is a circuit constant.
    pub fn is_constant(&self) -> bool {
        self.witness_index == IS_CONSTANT
    }

    /// Check if this field element is normalized
    /// (multiplicative_constant == 1 and additive_constant == 0, or is constant).
    pub fn is_normalized(&self) -> bool {
        self.is_constant()
            || (self.multiplicative_constant == Field::one()
                && self.additive_constant == Field::zero())
    }

    /// Check if two field elements share the same witness index.
    pub fn witness_indices_match(a: &FieldT<P>, b: &FieldT<P>) -> bool {
        a.witness_index == b.witness_index
    }

    // ════════════════════════════════════════════════════════════════════
    //  Normalization
    // ════════════════════════════════════════════════════════════════════

    /// Return a normalized copy: creates a gate so the result has
    /// `multiplicative_constant == 1` and `additive_constant == 0`.
    pub fn normalize(&self) -> Self {
        if self.is_normalized() {
            return self.clone();
        }
        let ctx = self.context.as_ref().expect("normalize requires context").clone();
        let mut builder = ctx.borrow_mut();

        let value = builder.base.get_variable(self.witness_index);
        let normalized_value = value * self.multiplicative_constant + self.additive_constant;
        let result_index = builder.base.add_variable(normalized_value);
        let zero_idx = builder.base.zero_idx();

        builder.create_add_gate(&AddTriple {
            a: self.witness_index,
            b: zero_idx,
            c: result_index,
            a_scaling: self.multiplicative_constant,
            b_scaling: Field::zero(),
            c_scaling: -Field::one(),
            const_scaling: self.additive_constant,
        });

        Self {
            context: self.context.clone(),
            additive_constant: Field::zero(),
            multiplicative_constant: Field::one(),
            witness_index: result_index,
        }
    }

    /// Get the normalized witness index. May add a gate for normalization.
    pub fn get_witness_index(&self) -> u32 {
        self.normalize().witness_index
    }

    // ════════════════════════════════════════════════════════════════════
    //  Public input / witness fixing
    // ════════════════════════════════════════════════════════════════════

    /// Make this field element a public input. Returns the public input position.
    pub fn set_public(&self) -> u32 {
        assert!(!self.is_constant(), "cannot set a constant as public input");
        let normalized = self.normalize();
        let ctx = self.context.as_ref().unwrap().clone();
        ctx.borrow_mut().base.set_public_input(normalized.witness_index)
    }

    /// Fix the witness value: constrains it with a selector so it cannot be
    /// changed by the prover.
    pub fn fix_witness(&mut self) {
        assert!(!self.is_constant(), "cannot fix a constant");
        assert!(self.context.is_some(), "fix_witness requires context");
        *self = self.normalize();
        let ctx = self.context.as_ref().unwrap().clone();
        let value = self.get_value();
        ctx.borrow_mut().fix_witness(self.witness_index, value);
    }

    /// Convert a constant field element to a fixed witness.
    pub fn convert_constant_to_fixed_witness(&mut self, ctx: BuilderRef<P>) {
        assert!(self.is_constant(), "must be constant");
        let value = self.additive_constant;
        *self = Self::from_witness(ctx.clone(), value);
        let current_value = self.get_value();
        ctx.borrow_mut().fix_witness(self.witness_index, current_value);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Arithmetic: addition
    // ════════════════════════════════════════════════════════════════════

    fn add_impl(&self, other: &FieldT<P>) -> FieldT<P> {
        let ctx = validate_two_contexts(&self.context, &other.context);
        assert!(
            ctx.is_some() || (self.is_constant() && other.is_constant()),
            "non-constant elements require context"
        );

        if Self::witness_indices_match(self, other) && !self.is_constant() {
            // Same witness: just update scaling factors
            FieldT {
                context: ctx,
                additive_constant: self.additive_constant + other.additive_constant,
                multiplicative_constant: self.multiplicative_constant
                    + other.multiplicative_constant,
                witness_index: self.witness_index,
            }
        } else if self.is_constant() && other.is_constant() {
            FieldT {
                context: ctx,
                additive_constant: self.additive_constant + other.additive_constant,
                multiplicative_constant: Field::one(),
                witness_index: IS_CONSTANT,
            }
        } else if !self.is_constant() && other.is_constant() {
            FieldT {
                context: ctx,
                additive_constant: self.additive_constant + other.additive_constant,
                multiplicative_constant: self.multiplicative_constant,
                witness_index: self.witness_index,
            }
        } else if self.is_constant() && !other.is_constant() {
            FieldT {
                context: ctx,
                additive_constant: self.additive_constant + other.additive_constant,
                multiplicative_constant: other.multiplicative_constant,
                witness_index: other.witness_index,
            }
        } else {
            // Both are distinct witnesses: create an add gate
            let result_index = {
                let ctx_ref = ctx.as_ref().unwrap();
                let mut builder = ctx_ref.borrow_mut();
                let left = builder.base.get_variable(self.witness_index);
                let right = builder.base.get_variable(other.witness_index);
                let result_value = left * self.multiplicative_constant
                    + right * other.multiplicative_constant
                    + self.additive_constant
                    + other.additive_constant;
                let result_index = builder.base.add_variable(result_value);

                builder.create_add_gate(&AddTriple {
                    a: self.witness_index,
                    b: other.witness_index,
                    c: result_index,
                    a_scaling: self.multiplicative_constant,
                    b_scaling: other.multiplicative_constant,
                    c_scaling: -Field::one(),
                    const_scaling: self.additive_constant + other.additive_constant,
                });
                result_index
            };

            FieldT {
                context: ctx,
                additive_constant: Field::zero(),
                multiplicative_constant: Field::one(),
                witness_index: result_index,
            }
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Arithmetic: multiplication
    // ════════════════════════════════════════════════════════════════════

    fn mul_impl(&self, other: &FieldT<P>) -> FieldT<P> {
        let ctx = validate_two_contexts(&self.context, &other.context);
        assert!(
            ctx.is_some() || (self.is_constant() && other.is_constant()),
            "non-constant elements require context"
        );

        if self.is_constant() && other.is_constant() {
            FieldT {
                context: ctx,
                additive_constant: self.additive_constant * other.additive_constant,
                multiplicative_constant: Field::one(),
                witness_index: IS_CONSTANT,
            }
        } else if !self.is_constant() && other.is_constant() {
            FieldT {
                context: ctx,
                additive_constant: self.additive_constant * other.additive_constant,
                multiplicative_constant: self.multiplicative_constant * other.additive_constant,
                witness_index: self.witness_index,
            }
        } else if self.is_constant() && !other.is_constant() {
            FieldT {
                context: ctx,
                additive_constant: self.additive_constant * other.additive_constant,
                multiplicative_constant: other.multiplicative_constant * self.additive_constant,
                witness_index: other.witness_index,
            }
        } else {
            // Both are witnesses: create arithmetic gate
            let q_m = self.multiplicative_constant * other.multiplicative_constant;
            let q_l = self.multiplicative_constant * other.additive_constant;
            let q_r = self.additive_constant * other.multiplicative_constant;
            let q_c = self.additive_constant * other.additive_constant;

            let result_index = {
                let ctx_ref = ctx.as_ref().unwrap();
                let mut builder = ctx_ref.borrow_mut();
                let left = builder.base.get_variable(self.witness_index);
                let right = builder.base.get_variable(other.witness_index);

                let result_value = left * right * q_m + left * q_l + right * q_r + q_c;
                let result_index = builder.base.add_variable(result_value);

                builder.create_arithmetic_gate(&ArithmeticTriple {
                    a: self.witness_index,
                    b: other.witness_index,
                    c: result_index,
                    q_m,
                    q_l,
                    q_r,
                    q_o: -Field::one(),
                    q_c,
                });
                result_index
            };

            FieldT {
                context: ctx,
                additive_constant: Field::zero(),
                multiplicative_constant: Field::one(),
                witness_index: result_index,
            }
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Arithmetic: division
    // ════════════════════════════════════════════════════════════════════

    /// Division without checking that the divisor is non-zero.
    pub fn divide_no_zero_check(&self, other: &FieldT<P>) -> FieldT<P> {
        let ctx = validate_two_contexts(&self.context, &other.context);
        assert!(
            ctx.is_some() || (self.is_constant() && other.is_constant()),
            "non-constant elements require context"
        );

        if self.is_constant() && other.is_constant() {
            let multiplier = if other.additive_constant != Field::zero() {
                other.additive_constant.invert()
            } else {
                Field::one()
            };
            FieldT {
                context: ctx,
                additive_constant: self.additive_constant * multiplier,
                multiplicative_constant: Field::one(),
                witness_index: IS_CONSTANT,
            }
        } else if !self.is_constant() && other.is_constant() {
            let multiplier = if other.additive_constant != Field::zero() {
                other.additive_constant.invert()
            } else {
                Field::one()
            };
            FieldT {
                context: ctx,
                additive_constant: self.additive_constant * multiplier,
                multiplicative_constant: self.multiplicative_constant * multiplier,
                witness_index: self.witness_index,
            }
        } else if self.is_constant() && !other.is_constant() {
            if self.get_value() == Field::zero() {
                FieldT {
                    context: ctx,
                    additive_constant: Field::zero(),
                    multiplicative_constant: Field::one(),
                    witness_index: IS_CONSTANT,
                }
            } else {
                let result_index = {
                    let ctx_ref = ctx.as_ref().unwrap();
                    let mut builder = ctx_ref.borrow_mut();

                    let numerator = self.get_value();
                    let denom = other.get_value_with_builder(&builder);
                    let denom_inv = if denom.is_zero() {
                        Field::zero()
                    } else {
                        denom.invert()
                    };
                    let out = numerator * denom_inv;
                    let result_index = builder.base.add_variable(out);

                    let q_m = other.multiplicative_constant;
                    let q_l = other.additive_constant;
                    let q_c = -self.get_value();

                    builder.create_arithmetic_gate(&ArithmeticTriple {
                        a: result_index,
                        b: other.witness_index,
                        c: result_index,
                        q_m,
                        q_l,
                        q_r: Field::zero(),
                        q_o: Field::zero(),
                        q_c,
                    });
                    result_index
                };

                FieldT {
                    context: ctx,
                    additive_constant: Field::zero(),
                    multiplicative_constant: Field::one(),
                    witness_index: result_index,
                }
            }
        } else {
            // Both witnesses
            let result_index = {
                let ctx_ref = ctx.as_ref().unwrap();
                let mut builder = ctx_ref.borrow_mut();

                let numerator = self.get_value_with_builder(&builder);
                let denom = other.get_value_with_builder(&builder);
                let denom_inv = if denom.is_zero() {
                    Field::zero()
                } else {
                    denom.invert()
                };
                let out = numerator * denom_inv;
                let result_index = builder.base.add_variable(out);

                builder.create_arithmetic_gate(&ArithmeticTriple {
                    a: result_index,
                    b: other.witness_index,
                    c: self.witness_index,
                    q_m: other.multiplicative_constant,
                    q_l: other.additive_constant,
                    q_r: Field::zero(),
                    q_o: -self.multiplicative_constant,
                    q_c: -self.additive_constant,
                });
                result_index
            };

            FieldT {
                context: ctx,
                additive_constant: Field::zero(),
                multiplicative_constant: Field::one(),
                witness_index: result_index,
            }
        }
    }

    /// Helper to get value when we already hold a borrow on the builder.
    fn get_value_with_builder(&self, builder: &UltraCircuitBuilder<P>) -> Field<P> {
        if self.is_constant() {
            self.additive_constant
        } else {
            let witness = builder.base.get_variable(self.witness_index);
            self.multiplicative_constant * witness + self.additive_constant
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Advanced arithmetic
    // ════════════════════════════════════════════════════════════════════

    /// Square: `self * self`.
    pub fn sqr(&self) -> Self {
        self.mul_impl(self)
    }

    /// Raise to a power with a constant u32 exponent.
    pub fn pow(&self, exponent: u32) -> Self {
        if self.is_constant() {
            return Self::from_field(
                self.get_value().pow(&[exponent as u64, 0, 0, 0]),
            );
        }
        if exponent == 0 {
            return Self::from_field(Field::one());
        }

        let mut accumulator_initialized = false;
        let mut accumulator = Self::default();
        let mut running_power = self.clone();
        let mut shifted = exponent;

        while shifted != 0 {
            if shifted & 1 != 0 {
                if !accumulator_initialized {
                    accumulator = running_power.clone();
                    accumulator_initialized = true;
                } else {
                    accumulator = accumulator.mul_impl(&running_power);
                }
            }
            if shifted >= 2 {
                running_power = running_power.sqr();
            }
            shifted >>= 1;
        }
        accumulator
    }

    /// Multiplicative inverse: `1 / self`.
    pub fn invert(&self) -> Self {
        if self.is_constant() {
            assert!(
                !self.get_value().is_zero(),
                "field_t::invert denominator is constant 0"
            );
        }

        if let Some(ctx) = &self.context {
            if self.get_value().is_zero() && !ctx.borrow().base.failed() {
                ctx.borrow_mut()
                    .base
                    .failure("field_t::invert denominator is 0".to_string());
            }
        }

        Self::from_field(Field::one()).divide_no_zero_check(self)
    }

    /// Efficiently compute `self * to_mul + to_add` using a big_mul_add gate.
    pub fn madd(&self, to_mul: &FieldT<P>, to_add: &FieldT<P>) -> FieldT<P> {
        let ctx = validate_three_contexts(&self.context, &to_mul.context, &to_add.context);

        if self.is_constant() || to_mul.is_constant() {
            return self.mul_impl(to_mul).add_impl(to_add);
        }

        let result_index = {
            let ctx_ref = ctx.as_ref().unwrap();
            let mut builder = ctx_ref.borrow_mut();

            let mul_scaling = self.multiplicative_constant * to_mul.multiplicative_constant;
            let a_scaling = self.multiplicative_constant * to_mul.additive_constant;
            let b_scaling = to_mul.multiplicative_constant * self.additive_constant;
            let c_scaling = to_add.multiplicative_constant;
            let const_scaling =
                self.additive_constant * to_mul.additive_constant + to_add.additive_constant;

            let a = if self.is_constant() {
                Field::zero()
            } else {
                builder.base.get_variable(self.witness_index)
            };
            let b = if to_mul.is_constant() {
                Field::zero()
            } else {
                builder.base.get_variable(to_mul.witness_index)
            };
            let c = if to_add.is_constant() {
                Field::zero()
            } else {
                builder.base.get_variable(to_add.witness_index)
            };

            let out = a * b * mul_scaling + a * a_scaling + b * b_scaling + c * c_scaling + const_scaling;
            let result_index = builder.base.add_variable(out);

            let zero_idx = builder.base.zero_idx();
            builder.create_big_mul_add_gate(
                &MulQuad {
                    a: if self.is_constant() { zero_idx } else { self.witness_index },
                    b: if to_mul.is_constant() { zero_idx } else { to_mul.witness_index },
                    c: if to_add.is_constant() { zero_idx } else { to_add.witness_index },
                    d: result_index,
                    mul_scaling,
                    a_scaling,
                    b_scaling,
                    c_scaling,
                    d_scaling: -Field::one(),
                    const_scaling,
                },
                false,
            );

            result_index
        };

        FieldT {
            context: ctx,
            additive_constant: Field::zero(),
            multiplicative_constant: Field::one(),
            witness_index: result_index,
        }
    }

    /// Efficiently compute `self + add_b + add_c` using a big_mul_add gate.
    pub fn add_two(&self, add_b: &FieldT<P>, add_c: &FieldT<P>) -> FieldT<P> {
        if self.is_constant() || add_b.is_constant() || add_c.is_constant() {
            return self.add_impl(add_b).add_impl(add_c);
        }

        let ctx = validate_three_contexts(&self.context, &add_b.context, &add_c.context);

        let result_index = {
            let ctx_ref = ctx.as_ref().unwrap();
            let mut builder = ctx_ref.borrow_mut();

            let a_scaling = self.multiplicative_constant;
            let b_scaling = add_b.multiplicative_constant;
            let c_scaling = add_c.multiplicative_constant;
            let const_scaling =
                self.additive_constant + add_b.additive_constant + add_c.additive_constant;

            let a = builder.base.get_variable(self.witness_index);
            let b = builder.base.get_variable(add_b.witness_index);
            let c = builder.base.get_variable(add_c.witness_index);
            let out = a * a_scaling + b * b_scaling + c * c_scaling + const_scaling;
            let result_index = builder.base.add_variable(out);

            let zero_idx = builder.base.zero_idx();
            builder.create_big_mul_add_gate(
                &MulQuad {
                    a: if self.is_constant() { zero_idx } else { self.witness_index },
                    b: if add_b.is_constant() { zero_idx } else { add_b.witness_index },
                    c: if add_c.is_constant() { zero_idx } else { add_c.witness_index },
                    d: result_index,
                    mul_scaling: Field::zero(),
                    a_scaling,
                    b_scaling,
                    c_scaling,
                    d_scaling: -Field::one(),
                    const_scaling,
                },
                false,
            );

            result_index
        };

        FieldT {
            context: ctx,
            additive_constant: Field::zero(),
            multiplicative_constant: Field::one(),
            witness_index: result_index,
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Assertions
    // ════════════════════════════════════════════════════════════════════

    /// Assert that `self` equals zero. Creates an arithmetic constraint.
    pub fn assert_is_zero(&self, msg: &str) {
        if self.is_constant() {
            assert!(
                self.additive_constant == Field::zero(),
                "{}",
                msg
            );
            return;
        }

        let ctx = self.context.as_ref().unwrap().clone();
        {
            let builder = ctx.borrow();
            if self.get_value_with_builder(&builder) != Field::zero() && !builder.base.failed() {
                drop(builder);
                ctx.borrow_mut().base.failure(msg.to_string());
            }
        }

        let mut builder = ctx.borrow_mut();
        let zero_idx = builder.base.zero_idx();
        builder.create_arithmetic_gate(&ArithmeticTriple {
            a: self.witness_index,
            b: zero_idx,
            c: zero_idx,
            q_m: Field::zero(),
            q_l: self.multiplicative_constant,
            q_r: Field::zero(),
            q_o: Field::zero(),
            q_c: self.additive_constant,
        });
    }

    /// Assert that `self` is not zero by proving it has an inverse.
    pub fn assert_is_not_zero(&self, msg: &str) {
        if self.is_constant() {
            assert!(
                self.additive_constant != Field::zero(),
                "{}",
                msg
            );
            return;
        }

        let ctx = self.context.as_ref().unwrap().clone();
        let value = self.get_value();

        if value == Field::zero() && !ctx.borrow().base.failed() {
            ctx.borrow_mut().base.failure(msg.to_string());
        }

        let inverse_value = if value.is_zero() {
            Field::zero()
        } else {
            value.invert()
        };

        let inverse = WitnessT::new(ctx.clone(), inverse_value);

        // Mark inverse as used
        {
            let mut builder = ctx.borrow_mut();
            builder.base.assert_equal(
                inverse.witness_index,
                inverse.witness_index,
                "mark_witness_as_used",
            );
        }

        // Create constraint: (self.v * self.mul + self.add) * inverse.v == 1
        let mut builder = ctx.borrow_mut();
        let zero_idx = builder.base.zero_idx();
        builder.create_arithmetic_gate(&ArithmeticTriple {
            a: self.witness_index,
            b: inverse.witness_index,
            c: zero_idx,
            q_m: self.multiplicative_constant,
            q_l: Field::zero(),
            q_r: self.additive_constant,
            q_o: Field::zero(),
            q_c: -Field::one(),
        });
    }

    /// Assert that `self == rhs`.
    pub fn assert_equal(&self, rhs: &FieldT<P>, msg: &str) {
        let ctx = validate_two_contexts(&self.context, &rhs.context);

        if self.is_constant() && rhs.is_constant() {
            assert!(
                self.get_value() == rhs.get_value(),
                "field_t::assert_equal: constants are not equal"
            );
            return;
        }

        let ctx_ref = ctx.as_ref().unwrap();

        if self.is_constant() {
            ctx_ref
                .borrow_mut()
                .assert_equal_constant(rhs.get_witness_index(), self.get_value(), msg);
        } else if rhs.is_constant() {
            ctx_ref
                .borrow_mut()
                .assert_equal_constant(self.get_witness_index(), rhs.get_value(), msg);
        } else if self.is_normalized() || rhs.is_normalized() {
            let lhs_idx = self.get_witness_index();
            let rhs_idx = rhs.get_witness_index();
            ctx_ref.borrow_mut().base.assert_equal(lhs_idx, rhs_idx, msg);
        } else {
            // Both non-normalized witnesses: single add gate for a - b = 0
            let values_equal = self.get_value() == rhs.get_value();
            if !values_equal && !ctx_ref.borrow().base.failed() {
                ctx_ref.borrow_mut().base.failure(msg.to_string());
            }

            let mut builder = ctx_ref.borrow_mut();
            let zero_idx = builder.base.zero_idx();
            builder.create_add_gate(&AddTriple {
                a: self.witness_index,
                b: rhs.witness_index,
                c: zero_idx,
                a_scaling: self.multiplicative_constant,
                b_scaling: -rhs.multiplicative_constant,
                c_scaling: Field::zero(),
                const_scaling: self.additive_constant - rhs.additive_constant,
            });
        }
    }

    /// Assert that `self != rhs`.
    pub fn assert_not_equal(&self, rhs: &FieldT<P>, msg: &str) {
        let neg_rhs = -(rhs.clone());
        let diff = self.add_impl(&neg_rhs);
        diff.assert_is_not_zero(msg);
    }

    /// Assert that `self` is in the given set.
    pub fn assert_is_in_set(&self, set: &[FieldT<P>], msg: &str) {
        let neg_first = -(set[0].clone());
        let mut product = self.add_impl(&neg_first);
        for item in &set[1..] {
            let neg_item = -(item.clone());
            product = product.mul_impl(&self.add_impl(&neg_item));
        }
        product.assert_is_zero(msg);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Conditional assignment & equality
    // ════════════════════════════════════════════════════════════════════

    /// If `predicate` is true, return `lhs`; otherwise return `rhs`.
    ///
    /// result = rhs + predicate * (lhs - rhs)
    pub fn conditional_assign(
        predicate: &super::bool::BoolT<P>,
        lhs: &FieldT<P>,
        rhs: &FieldT<P>,
    ) -> FieldT<P> {
        let predicate_field = bool_to_field(predicate);
        let diff = lhs.clone() - rhs.clone();
        diff.madd(&predicate_field, rhs)
    }

    /// Returns a `BoolT` that is true iff `self == other` (in-circuit).
    pub fn is_equal(&self, other: &FieldT<P>) -> super::bool::BoolT<P> {
        use super::bool::BoolT;

        let diff = self.clone() - other.clone();
        if diff.is_constant() {
            return BoolT::from_constant(diff.get_value().is_zero());
        }

        let ctx = diff.context.as_ref().unwrap().clone();
        let value = diff.get_value();
        let is_zero_val = value.is_zero();
        let inverse_val = if value.is_zero() {
            Field::zero()
        } else {
            value.invert()
        };

        let is_zero_witness = WitnessT::new(ctx.clone(), Field::from(is_zero_val as u64));
        let inverse_witness = WitnessT::new(ctx.clone(), inverse_val);

        let is_zero_field = FieldT::from_witness_t(&is_zero_witness);
        let inverse_field = FieldT::from_witness_t(&inverse_witness);

        // Constraint 1: diff * inverse = 1 - is_zero
        // i.e. diff * inverse + is_zero - 1 = 0
        let one = FieldT::<P>::from_u64(1);
        FieldT::evaluate_polynomial_identity(
            &diff,
            &inverse_field,
            &is_zero_field,
            &(-one),
            "field_t::is_equal: diff * inverse + is_zero - 1 = 0",
        );

        // Constraint 2: diff * is_zero = 0
        FieldT::evaluate_polynomial_identity(
            &diff,
            &is_zero_field,
            &FieldT::default(),
            &FieldT::default(),
            "field_t::is_equal: diff * is_zero = 0",
        );

        // Constrain is_zero to be boolean
        ctx.borrow_mut().create_bool_gate(is_zero_field.normalize().witness_index);

        BoolT::from_witness_index_unsafe(ctx, is_zero_field.normalize().witness_index)
    }

    // ════════════════════════════════════════════════════════════════════
    //  Range constraints
    // ════════════════════════════════════════════════════════════════════

    /// Constrain the normalized value to be less than 2^num_bits.
    pub fn create_range_constraint(&self, num_bits: usize, msg: &str) {
        if num_bits == 0 {
            self.assert_is_zero("0-bit range_constraint on non-zero field_t.");
            return;
        }
        if self.is_constant() {
            // Check that the value fits in num_bits
            let val = self.get_value().from_montgomery_form();
            let msb = get_msb_u256(&val.data);
            assert!(msb < num_bits as u32, "{}", msg);
        } else {
            let ctx = self.context.as_ref().unwrap().clone();
            let normalized = self.normalize();
            ctx.borrow_mut().decompose_into_default_range(
                normalized.witness_index,
                num_bits as u64,
                ultra_builder::DEFAULT_PLOOKUP_RANGE_BITNUM,
                msg,
            );
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Polynomial identity evaluation
    // ════════════════════════════════════════════════════════════════════

    /// Constrain `a + b + c + d = 0` using a big_add gate.
    pub fn evaluate_linear_identity(
        a: &FieldT<P>,
        b: &FieldT<P>,
        c: &FieldT<P>,
        d: &FieldT<P>,
        msg: &str,
    ) {
        let ctx = validate_four_contexts(&a.context, &b.context, &c.context, &d.context);

        if a.is_constant() && b.is_constant() && c.is_constant() && d.is_constant() {
            assert!(
                (a.get_value() + b.get_value() + c.get_value() + d.get_value()).is_zero(),
                "{}",
                msg
            );
            return;
        }

        let ctx_ref = ctx.as_ref().unwrap();

        let identity_holds =
            (a.get_value() + b.get_value() + c.get_value() + d.get_value()).is_zero();
        if !identity_holds && !ctx_ref.borrow().base.failed() {
            ctx_ref.borrow_mut().base.failure(msg.to_string());
        }

        let const_scaling =
            a.additive_constant + b.additive_constant + c.additive_constant + d.additive_constant;

        let mut builder = ctx_ref.borrow_mut();
        let zero_idx = builder.base.zero_idx();
        builder.create_big_add_gate(
            &AddQuad {
                a: if a.is_constant() { zero_idx } else { a.witness_index },
                b: if b.is_constant() { zero_idx } else { b.witness_index },
                c: if c.is_constant() { zero_idx } else { c.witness_index },
                d: if d.is_constant() { zero_idx } else { d.witness_index },
                a_scaling: a.multiplicative_constant,
                b_scaling: b.multiplicative_constant,
                c_scaling: c.multiplicative_constant,
                d_scaling: d.multiplicative_constant,
                const_scaling,
            },
            false,
        );
    }

    /// Constrain `a * b + c + d = 0` using a big_mul_add gate.
    pub fn evaluate_polynomial_identity(
        a: &FieldT<P>,
        b: &FieldT<P>,
        c: &FieldT<P>,
        d: &FieldT<P>,
        msg: &str,
    ) {
        if a.is_constant() && b.is_constant() && c.is_constant() && d.is_constant() {
            assert!(
                (a.get_value() * b.get_value() + c.get_value() + d.get_value()).is_zero(),
                "{}",
                msg
            );
            return;
        }

        let ctx = validate_four_contexts(&a.context, &b.context, &c.context, &d.context);
        let ctx_ref = ctx.as_ref().unwrap();

        let identity_holds =
            (a.get_value() * b.get_value() + c.get_value() + d.get_value()).is_zero();
        if !identity_holds && !ctx_ref.borrow().base.failed() {
            ctx_ref.borrow_mut().base.failure(msg.to_string());
        }

        let mul_scaling = a.multiplicative_constant * b.multiplicative_constant;
        let a_scaling = a.multiplicative_constant * b.additive_constant;
        let b_scaling = b.multiplicative_constant * a.additive_constant;
        let c_scaling = c.multiplicative_constant;
        let d_scaling = d.multiplicative_constant;
        let const_scaling =
            a.additive_constant * b.additive_constant + c.additive_constant + d.additive_constant;

        let mut builder = ctx_ref.borrow_mut();
        let zero_idx = builder.base.zero_idx();
        builder.create_big_mul_add_gate(
            &MulQuad {
                a: if a.is_constant() { zero_idx } else { a.witness_index },
                b: if b.is_constant() { zero_idx } else { b.witness_index },
                c: if c.is_constant() { zero_idx } else { c.witness_index },
                d: if d.is_constant() { zero_idx } else { d.witness_index },
                mul_scaling,
                a_scaling,
                b_scaling,
                c_scaling,
                d_scaling,
                const_scaling,
            },
            false,
        );
    }

    /// Efficiently compute the sum of a vector of field elements.
    pub fn accumulate(input: &[FieldT<P>]) -> FieldT<P> {
        if input.is_empty() {
            return FieldT::from_field(Field::zero());
        }
        if input.len() == 1 {
            return input[0].normalize();
        }

        let mut accumulator: Vec<FieldT<P>> = Vec::new();
        let mut constant_term = FieldT::from_field(Field::zero());

        for element in input {
            if element.is_constant() {
                constant_term = constant_term.add_impl(element);
            } else {
                accumulator.push(element.clone());
            }
        }
        if accumulator.is_empty() {
            return constant_term;
        }
        accumulator[0] = accumulator[0].add_impl(&constant_term);

        // Extract context from the first witness
        let ctx = accumulator[0].context.as_ref().unwrap().clone();

        // Compute total output value
        let mut output = Field::zero();
        for acc in &accumulator {
            output = output + acc.get_value();
        }

        let num_elements = accumulator.len();
        // Pad to multiple of 3
        let padding = if num_elements % 3 == 0 {
            0
        } else {
            3 - (num_elements % 3)
        };
        for _ in 0..padding {
            let zero_idx = ctx.borrow().base.zero_idx();
            accumulator.push(FieldT::from_witness_index(ctx.clone(), zero_idx));
        }

        let padded_len = accumulator.len();
        let num_gates = padded_len / 3;
        let last_gate_idx = num_gates - 1;

        let total = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), output));
        let mut accumulating_total = total.clone();

        for i in 0..last_gate_idx {
            let mut builder = ctx.borrow_mut();
            let const_scaling = accumulator[3 * i].additive_constant
                + accumulator[3 * i + 1].additive_constant
                + accumulator[3 * i + 2].additive_constant;

            builder.create_big_add_gate(
                &AddQuad {
                    a: accumulator[3 * i].witness_index,
                    b: accumulator[3 * i + 1].witness_index,
                    c: accumulator[3 * i + 2].witness_index,
                    d: accumulating_total.witness_index,
                    a_scaling: accumulator[3 * i].multiplicative_constant,
                    b_scaling: accumulator[3 * i + 1].multiplicative_constant,
                    c_scaling: accumulator[3 * i + 2].multiplicative_constant,
                    d_scaling: -Field::one(),
                    const_scaling,
                },
                true, // use_next_gate_w_4
            );

            let new_total = accumulating_total.get_value_with_builder(&builder)
                - accumulator[3 * i].get_value_with_builder(&builder)
                - accumulator[3 * i + 1].get_value_with_builder(&builder)
                - accumulator[3 * i + 2].get_value_with_builder(&builder);
            let new_total_idx = builder.base.add_variable(new_total);
            drop(builder);
            accumulating_total = FieldT::from_witness_index(ctx.clone(), new_total_idx);
        }

        // Last gate (no use_next_gate_w_4)
        {
            let mut builder = ctx.borrow_mut();
            let i = last_gate_idx;
            let const_scaling = accumulator[3 * i].additive_constant
                + accumulator[3 * i + 1].additive_constant
                + accumulator[3 * i + 2].additive_constant;

            builder.create_big_add_gate(
                &AddQuad {
                    a: accumulator[3 * i].witness_index,
                    b: accumulator[3 * i + 1].witness_index,
                    c: accumulator[3 * i + 2].witness_index,
                    d: accumulating_total.witness_index,
                    a_scaling: accumulator[3 * i].multiplicative_constant,
                    b_scaling: accumulator[3 * i + 1].multiplicative_constant,
                    c_scaling: accumulator[3 * i + 2].multiplicative_constant,
                    d_scaling: -Field::one(),
                    const_scaling,
                },
                false,
            );
        }

        total.normalize()
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Operator implementations
// ════════════════════════════════════════════════════════════════════════

impl<P: FieldParams> Add for FieldT<P> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.add_impl(&rhs)
    }
}

impl<P: FieldParams> Add<&FieldT<P>> for FieldT<P> {
    type Output = FieldT<P>;
    fn add(self, rhs: &FieldT<P>) -> FieldT<P> {
        self.add_impl(rhs)
    }
}

impl<P: FieldParams> Add<&FieldT<P>> for &FieldT<P> {
    type Output = FieldT<P>;
    fn add(self, rhs: &FieldT<P>) -> FieldT<P> {
        self.add_impl(rhs)
    }
}

impl<P: FieldParams> Sub for FieldT<P> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self.add_impl(&(-rhs))
    }
}

impl<P: FieldParams> Sub<&FieldT<P>> for &FieldT<P> {
    type Output = FieldT<P>;
    fn sub(self, rhs: &FieldT<P>) -> FieldT<P> {
        let neg_rhs = -(rhs.clone());
        self.add_impl(&neg_rhs)
    }
}

impl<P: FieldParams> Mul for FieldT<P> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.mul_impl(&rhs)
    }
}

impl<P: FieldParams> Mul<&FieldT<P>> for &FieldT<P> {
    type Output = FieldT<P>;
    fn mul(self, rhs: &FieldT<P>) -> FieldT<P> {
        self.mul_impl(rhs)
    }
}

impl<P: FieldParams> Div for FieldT<P> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        rhs.assert_is_not_zero("field_t::operator/ divisor is 0");
        self.divide_no_zero_check(&rhs)
    }
}

impl<P: FieldParams> Div<&FieldT<P>> for &FieldT<P> {
    type Output = FieldT<P>;
    fn div(self, rhs: &FieldT<P>) -> FieldT<P> {
        rhs.assert_is_not_zero("field_t::operator/ divisor is 0");
        self.divide_no_zero_check(rhs)
    }
}

impl<P: FieldParams> Neg for FieldT<P> {
    type Output = Self;
    fn neg(mut self) -> Self {
        self.additive_constant = -self.additive_constant;
        if !self.is_constant() {
            self.multiplicative_constant = -self.multiplicative_constant;
        }
        self
    }
}

impl<P: FieldParams> AddAssign for FieldT<P> {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.add_impl(&rhs);
    }
}

impl<P: FieldParams> SubAssign for FieldT<P> {
    fn sub_assign(&mut self, rhs: Self) {
        let neg_rhs = -rhs;
        *self = self.add_impl(&neg_rhs);
    }
}

impl<P: FieldParams> MulAssign for FieldT<P> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.mul_impl(&rhs);
    }
}

impl<P: FieldParams> DivAssign for FieldT<P> {
    fn div_assign(&mut self, rhs: Self) {
        rhs.assert_is_not_zero("field_t::operator/= divisor is 0");
        *self = self.divide_no_zero_check(&rhs);
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Utility functions
// ════════════════════════════════════════════════════════════════════════

/// Convert a `BoolT` into a `FieldT` (0 or 1).
pub fn bool_to_field<P: FieldParams>(b: &super::bool::BoolT<P>) -> FieldT<P> {
    if b.is_constant() {
        FieldT::from_field(Field::from(b.get_value() as u64))
    } else {
        let normalized = b.normalize();
        FieldT::from_witness_index(
            normalized.context.as_ref().unwrap().clone(),
            normalized.witness_index,
        )
    }
}

/// Validate that `lo + hi * 2^lo_bits < field_modulus`.
///
/// Uses a borrow-subtraction technique. Range constraints on lo and hi are
/// assumed to be enforced externally (e.g., by the batch_mul algorithm).
pub fn validate_split_in_field_unsafe<P: FieldParams>(
    lo: &FieldT<P>,
    hi: &FieldT<P>,
    lo_bits: usize,
    field_modulus: &[u64; 4],
) {
    let total_bits = get_msb_u256(field_modulus) as usize + 1;
    let hi_bits = total_bits - lo_bits;

    // Split the field modulus: r_lo = modulus mod 2^lo_bits, r_hi = modulus >> lo_bits
    let r_lo_limbs = u256_mask_lo(field_modulus, lo_bits);
    let r_hi_limbs = u256_shr(field_modulus, lo_bits);

    // Check if we need to borrow: borrow = 1 if lo_value > r_lo
    let lo_val = lo.get_value().from_montgomery_form();
    let need_borrow = u256_gt(&lo_val.data, &r_lo_limbs);

    let borrow = if lo.is_constant() {
        FieldT::from_field(Field::from(need_borrow as u64))
    } else {
        let ctx = lo.context.as_ref().unwrap().clone();
        let borrow_field =
            FieldT::from_witness(ctx.clone(), Field::from(need_borrow as u64));
        ctx.borrow_mut().create_new_range_constraint(
            borrow_field.normalize().witness_index,
            1,
            "borrow",
        );
        borrow_field
    };

    let r_hi_field = FieldT::<P>::from_field(Field::from_limbs(r_hi_limbs));
    let r_lo_field = FieldT::<P>::from_field(Field::from_limbs(r_lo_limbs));
    let shift = {
        let mut s = [0u64; 4];
        let limb_idx = lo_bits / 64;
        let bit_idx = lo_bits % 64;
        if limb_idx < 4 {
            s[limb_idx] = 1u64 << bit_idx;
        }
        FieldT::<P>::from_field(Field::from_limbs(s))
    };

    // Hi range check = r_hi - hi - borrow
    // Lo range check = r_lo - lo + borrow * 2^lo_bits
    let hi_diff = (r_hi_field - hi.clone()) - borrow.clone();
    let lo_diff = (r_lo_field - lo.clone()) + borrow.mul_impl(&shift);

    hi_diff.create_range_constraint(hi_bits, "validate_split_in_field_unsafe hi");
    lo_diff.create_range_constraint(lo_bits, "validate_split_in_field_unsafe lo");
}

/// Mask a [u64; 4] to keep only the lower `bits` bits.
fn u256_mask_lo(v: &[u64; 4], bits: usize) -> [u64; 4] {
    let mut result = [0u64; 4];
    let full_limbs = bits / 64;
    let remaining = bits % 64;
    for i in 0..full_limbs.min(4) {
        result[i] = v[i];
    }
    if full_limbs < 4 && remaining > 0 {
        result[full_limbs] = v[full_limbs] & ((1u64 << remaining) - 1);
    }
    result
}

/// Shift a [u64; 4] right by `bits` bits.
fn u256_shr(v: &[u64; 4], bits: usize) -> [u64; 4] {
    let mut result = [0u64; 4];
    let limb_shift = bits / 64;
    let bit_shift = bits % 64;
    for i in 0..(4 - limb_shift) {
        result[i] = v[i + limb_shift] >> bit_shift;
        if bit_shift > 0 && i + limb_shift + 1 < 4 {
            result[i] |= v[i + limb_shift + 1] << (64 - bit_shift);
        }
    }
    result
}

/// Compare two [u64; 4] values: returns true if a > b.
fn u256_gt(a: &[u64; 4], b: &[u64; 4]) -> bool {
    for i in (0..4).rev() {
        if a[i] > b[i] {
            return true;
        }
        if a[i] < b[i] {
            return false;
        }
    }
    false
}

// ════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_ecc::curves::bn254::Bn254FrParams;

    type Fr = Field<Bn254FrParams>;
    type FrField = FieldT<Bn254FrParams>;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    fn check_circuit(builder: &BuilderRef<Bn254FrParams>) -> Result<(), String> {
        UltraCircuitChecker::check(&mut builder.borrow_mut())
    }

    fn num_gates(builder: &BuilderRef<Bn254FrParams>) -> usize {
        builder.borrow().base.num_gates()
    }

    // ── Constructor tests ─────────────────────────────────────────────

    #[test]
    fn test_constructor_from_witness() {
        let builder = make_builder();
        let value = Fr::random_element();
        let w = WitnessT::new(builder.clone(), value);
        let a = FrField::from_witness_t(&w);
        assert_eq!(a.get_value(), value);
        assert!(!a.is_constant());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_field_constant_from_u64() {
        let a = FrField::from_u64(42);
        assert_eq!(a.get_value(), Fr::from(42u64));
        assert!(a.is_constant());
    }

    #[test]
    fn test_field_constant_from_field() {
        let val = Fr::random_element();
        let a = FrField::from_field(val);
        assert_eq!(a.get_value(), val);
        assert!(a.is_constant());
    }

    #[test]
    fn test_field_from_witness_index() {
        let builder = make_builder();
        let val = Fr::from(99u64);
        let idx = builder.borrow_mut().base.add_variable(val);
        let a = FrField::from_witness_index(builder.clone(), idx);
        assert_eq!(a.get_value(), val);
        assert!(!a.is_constant());
    }

    #[test]
    fn test_field_from_witness() {
        let builder = make_builder();
        let val = Fr::random_element();
        let a = FrField::from_witness(builder.clone(), val);
        assert_eq!(a.get_value(), val);
        assert!(!a.is_constant());
    }

    #[test]
    fn test_field_is_constant() {
        let a = FrField::from_u64(5);
        assert!(a.is_constant());

        let builder = make_builder();
        let b = FrField::from_witness(builder, Fr::from(5u64));
        assert!(!b.is_constant());
    }

    #[test]
    fn test_field_is_normalized() {
        let a = FrField::from_u64(5);
        assert!(a.is_normalized()); // constants are always normalized

        let builder = make_builder();
        let b = FrField::from_witness(builder, Fr::from(5u64));
        assert!(b.is_normalized()); // fresh witness has mul=1, add=0

        // Multiply by constant changes scaling factors
        let c = b * FrField::from_u64(3);
        assert!(!c.is_normalized()); // mul=3, add=0
    }

    #[test]
    fn test_field_get_value() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(10u64));
        assert_eq!(a.get_value(), Fr::from(10u64));

        // After scaling by constant, value should still be correct
        let b = a * FrField::from_u64(3);
        assert_eq!(b.get_value(), Fr::from(30u64));
    }

    #[test]
    fn test_field_witness_indices_match() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let b = a.clone();
        assert!(FrField::witness_indices_match(&a, &b));

        let c = FrField::from_witness(builder, Fr::from(5u64));
        assert!(!FrField::witness_indices_match(&a, &c));
    }

    // ── Addition tests ────────────────────────────────────────────────

    #[test]
    fn test_add() {
        let builder = make_builder();

        // witness + witness
        let a_val = Fr::random_element();
        let b_val = Fr::random_element();
        let a = FrField::from_witness(builder.clone(), a_val);
        let b = FrField::from_witness(builder.clone(), b_val);
        let c = a.clone() + b.clone();
        assert_eq!(c.get_value(), a_val + b_val);

        // witness + constant
        let d = a.clone() + FrField::from_u64(5);
        assert_eq!(d.get_value(), a_val + Fr::from(5u64));

        // constant + constant
        let e = FrField::from_u64(3) + FrField::from_u64(7);
        assert_eq!(e.get_value(), Fr::from(10u64));
        assert!(e.is_constant());

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_same_witness() {
        let builder = make_builder();
        let val = Fr::from(5u64);
        let a = FrField::from_witness(builder.clone(), val);
        let b = a.clone() + a.clone();
        // same witness index, just double the multiplicative constant
        assert_eq!(b.get_value(), val + val);
        assert!(FrField::witness_indices_match(&a, &b));
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_constant_witness() {
        let builder = make_builder();
        let val = Fr::from(10u64);
        let constant = FrField::from_u64(7);
        let witness = FrField::from_witness(builder.clone(), val);
        let result = constant + witness;
        assert_eq!(result.get_value(), val + Fr::from(7u64));
        assert!(!result.is_constant());
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Subtraction tests ─────────────────────────────────────────────

    #[test]
    fn test_sub() {
        let builder = make_builder();
        let a_val = Fr::random_element();
        let b_val = Fr::random_element();
        let a = FrField::from_witness(builder.clone(), a_val);
        let b = FrField::from_witness(builder.clone(), b_val);
        let c = a - b;
        assert_eq!(c.get_value(), a_val - b_val);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_sub_constants() {
        let a = FrField::from_u64(10);
        let b = FrField::from_u64(3);
        let c = a - b;
        assert_eq!(c.get_value(), Fr::from(7u64));
        assert!(c.is_constant());
    }

    // ── Multiplication tests ──────────────────────────────────────────

    #[test]
    fn test_mul() {
        let builder = make_builder();
        let a_val = Fr::random_element();
        let b_val = Fr::random_element();
        let a = FrField::from_witness(builder.clone(), a_val);
        let b = FrField::from_witness(builder.clone(), b_val);
        let c = a * b;
        assert_eq!(c.get_value(), a_val * b_val);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_mul_constants() {
        let a = FrField::from_u64(6);
        let b = FrField::from_u64(7);
        let c = a * b;
        assert_eq!(c.get_value(), Fr::from(42u64));
        assert!(c.is_constant());
    }

    #[test]
    fn test_mul_constant_witness() {
        let builder = make_builder();
        let val = Fr::from(5u64);
        let a = FrField::from_witness(builder.clone(), val);
        let b = FrField::from_u64(3);
        let c = a * b;
        assert_eq!(c.get_value(), Fr::from(15u64));
        assert!(!c.is_constant());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_mul_by_zero() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::random_element());
        let zero = FrField::from_u64(0);
        let c = a * zero;
        assert_eq!(c.get_value(), Fr::zero());
    }

    // ── Division tests ────────────────────────────────────────────────

    #[test]
    fn test_div() {
        let builder = make_builder();

        // witness / witness
        let a_val = Fr::random_element();
        let b_val = Fr::random_element();
        let a = FrField::from_witness(builder.clone(), a_val);
        let b = FrField::from_witness(builder.clone(), b_val);
        let c = a.clone() / b.clone();
        assert_eq!(c.get_value(), a_val * b_val.invert());

        // witness / constant
        let d = a.clone() / FrField::from_u64(5);
        assert_eq!(d.get_value(), a_val * Fr::from(5u64).invert());

        // constant / witness
        let e = FrField::from_u64(10) / b.clone();
        assert_eq!(e.get_value(), Fr::from(10u64) * b_val.invert());

        // 0 / witness
        let f = FrField::from_u64(0) / b;
        assert_eq!(f.get_value(), Fr::zero());
        assert!(f.is_constant()); // numerator 0 produces constant 0

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_div_constant_witness() {
        let builder = make_builder();
        let val = Fr::from(15u64);
        let a = FrField::from_witness(builder.clone(), val);
        let b = FrField::from_u64(3);
        let c = a / b;
        assert_eq!(c.get_value(), Fr::from(5u64));
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_div_edge_cases() {
        // const / const
        let a = FrField::from_u64(10);
        let b = FrField::from_u64(2);
        let c = a / b;
        assert_eq!(c.get_value(), Fr::from(5u64));
    }

    // ── Negation tests ────────────────────────────────────────────────

    #[test]
    fn test_neg() {
        let builder = make_builder();
        let val = Fr::from(7u64);
        let a = FrField::from_witness(builder.clone(), val);
        let b = -a;
        assert_eq!(b.get_value(), -val);
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Sqr / pow / invert tests ──────────────────────────────────────

    #[test]
    fn test_sqr() {
        let builder = make_builder();
        let val = Fr::from(5u64);
        let a = FrField::from_witness(builder.clone(), val);
        let b = a.sqr();
        assert_eq!(b.get_value(), Fr::from(25u64));
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_pow_zero_exponent() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let b = a.pow(0);
        assert_eq!(b.get_value(), Fr::one());
    }

    #[test]
    fn test_pow_one_exponent() {
        let builder = make_builder();
        let val = Fr::from(5u64);
        let a = FrField::from_witness(builder.clone(), val);
        let b = a.pow(1);
        assert_eq!(b.get_value(), val);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_pow_constant_base() {
        let a = FrField::from_u64(3);
        let b = a.pow(4);
        assert_eq!(b.get_value(), Fr::from(81u64));
        assert!(b.is_constant());
    }

    #[test]
    fn test_pow() {
        let builder = make_builder();
        let val = Fr::from(2u64);
        let a = FrField::from_witness(builder.clone(), val);

        let b = a.pow(10);
        assert_eq!(b.get_value(), Fr::from(1024u64));
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_invert() {
        let builder = make_builder();
        let val = Fr::from(5u64);
        let a = FrField::from_witness(builder.clone(), val);
        let b = a.invert();
        assert_eq!(b.get_value(), val.invert());

        // constant invert
        let c = FrField::from_u64(7);
        let d = c.invert();
        assert_eq!(d.get_value(), Fr::from(7u64).invert());

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    #[should_panic(expected = "field_t::invert denominator is constant 0")]
    fn test_invert_zero() {
        let a = FrField::from_u64(0);
        let _ = a.invert();
    }

    // ── Increment tests ───────────────────────────────────────────────

    #[test]
    fn test_prefix_increment() {
        let builder = make_builder();
        let mut a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        a += FrField::from_u64(1);
        assert_eq!(a.get_value(), Fr::from(6u64));
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_postfix_increment() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let old = a.clone();
        let b = a + FrField::from_u64(1);
        assert_eq!(old.get_value(), Fr::from(5u64));
        assert_eq!(b.get_value(), Fr::from(6u64));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Compound assignment tests ─────────────────────────────────────

    #[test]
    fn test_compound_assign() {
        let builder = make_builder();
        let mut a = FrField::from_witness(builder.clone(), Fr::from(10u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(3u64));

        a += b.clone();
        assert_eq!(a.get_value(), Fr::from(13u64));

        a -= FrField::from_u64(3);
        assert_eq!(a.get_value(), Fr::from(10u64));

        a *= FrField::from_u64(2);
        assert_eq!(a.get_value(), Fr::from(20u64));

        assert!(check_circuit(&builder).is_ok());
    }

    // ── Madd / add_two tests ──────────────────────────────────────────

    #[test]
    fn test_madd() {
        let builder = make_builder();
        let a_val = Fr::random_element();
        let b_val = Fr::random_element();
        let c_val = Fr::random_element();

        // All witnesses
        let a = FrField::from_witness(builder.clone(), a_val);
        let b = FrField::from_witness(builder.clone(), b_val);
        let c = FrField::from_witness(builder.clone(), c_val);
        let result = a.madd(&b, &c);
        assert_eq!(result.get_value(), a_val * b_val + c_val);

        // With constant multiplicand
        let d = FrField::from_u64(5);
        let result2 = a_val_field(&builder, a_val).madd(&d, &c_val_field(&builder, c_val));
        assert_eq!(result2.get_value(), a_val * Fr::from(5u64) + c_val);

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_madd_all_constants() {
        let a = FrField::from_u64(3);
        let b = FrField::from_u64(4);
        let c = FrField::from_u64(5);
        let result = a.madd(&b, &c);
        assert_eq!(result.get_value(), Fr::from(17u64)); // 3*4 + 5
        assert!(result.is_constant());
    }

    #[test]
    fn test_add_two() {
        let builder = make_builder();
        let a_val = Fr::random_element();
        let b_val = Fr::random_element();
        let c_val = Fr::random_element();

        let a = FrField::from_witness(builder.clone(), a_val);
        let b = FrField::from_witness(builder.clone(), b_val);
        let c = FrField::from_witness(builder.clone(), c_val);
        let result = a.add_two(&b, &c);
        assert_eq!(result.get_value(), a_val + b_val + c_val);
        assert!(check_circuit(&builder).is_ok());
    }

    // helpers for tests
    fn a_val_field(
        builder: &BuilderRef<Bn254FrParams>,
        val: Fr,
    ) -> FrField {
        FrField::from_witness(builder.clone(), val)
    }
    fn c_val_field(
        builder: &BuilderRef<Bn254FrParams>,
        val: Fr,
    ) -> FrField {
        FrField::from_witness(builder.clone(), val)
    }

    #[test]
    fn test_madd_add_two_gate_count() {
        // Test gate count for madd and add_two with various const/witness combos
        // When all 3 are witnesses, madd uses 1 big_mul_add gate
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(3u64));
        let c = FrField::from_witness(builder.clone(), Fr::from(4u64));

        let before = num_gates(&builder);
        let _ = a.madd(&b, &c);
        let after = num_gates(&builder);
        assert_eq!(after - before, 1, "madd(w,w,w) should use 1 gate");

        // add_two with all witnesses: 1 gate
        let d = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let e = FrField::from_witness(builder.clone(), Fr::from(6u64));
        let f = FrField::from_witness(builder.clone(), Fr::from(7u64));
        let before = num_gates(&builder);
        let _ = d.add_two(&e, &f);
        let after = num_gates(&builder);
        assert_eq!(after - before, 1, "add_two(w,w,w) should use 1 gate");

        assert!(check_circuit(&builder).is_ok());
    }

    // ── Normalize tests ───────────────────────────────────────────────

    #[test]
    fn test_normalize() {
        let builder = make_builder();
        let val = Fr::from(5u64);
        let a = FrField::from_witness(builder.clone(), val);

        // Multiply by constant to make non-normalized
        let b = a * FrField::from_u64(3);
        assert!(!b.is_normalized());
        assert_eq!(b.get_value(), Fr::from(15u64));

        let c = b.normalize();
        assert!(c.is_normalized());
        assert_eq!(c.get_value(), Fr::from(15u64));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Assert tests ──────────────────────────────────────────────────

    #[test]
    fn test_assert_equal() {
        let builder = make_builder();
        let val = Fr::random_element();

        // const == const
        let a = FrField::from_field(val);
        let b = FrField::from_field(val);
        a.assert_equal(&b, "const equal");

        // witness == same-value witness
        let c = FrField::from_witness(builder.clone(), val);
        let d = FrField::from_witness(builder.clone(), val);
        c.assert_equal(&d, "witness equal");

        // witness == constant
        let e = FrField::from_witness(builder.clone(), val);
        let f = FrField::from_field(val);
        e.assert_equal(&f, "witness == const");

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_equal_gate_count() {
        let builder = make_builder();
        let val = Fr::from(7u64);

        // Both normalized: copy constraint only (no gate)
        let a = FrField::from_witness(builder.clone(), val);
        let b = FrField::from_witness(builder.clone(), val);
        let before = num_gates(&builder);
        a.assert_equal(&b, "test");
        let after = num_gates(&builder);
        // assert_equal with both normalized uses a copy constraint, 0 gates
        assert_eq!(after - before, 0, "normalized assert_equal = 0 gates");

        // One constant: assert_equal_constant
        let c = FrField::from_witness(builder.clone(), val);
        let d = FrField::from_field(val);
        let before = num_gates(&builder);
        c.assert_equal(&d, "test");
        let after = num_gates(&builder);
        // assert_equal_constant adds 1 gate (fix_witness for new constant variable)
        assert_eq!(after - before, 1, "const assert_equal = 1 gate (fix_witness)");

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_is_zero() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::zero());
        a.assert_is_zero("should be zero");
        assert!(check_circuit(&builder).is_ok());

        // Constant zero
        let b = FrField::from_u64(0);
        b.assert_is_zero("const zero");
    }

    #[test]
    fn test_assert_is_not_zero() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        a.assert_is_not_zero("should not be zero");
        assert!(check_circuit(&builder).is_ok());

        // Constant non-zero
        let b = FrField::from_u64(7);
        b.assert_is_not_zero("const non-zero");
    }

    #[test]
    fn test_assert_not_equal() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(7u64));
        a.assert_not_equal(&b, "should differ");
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_is_in_set() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(3u64));
        let set = vec![
            FrField::from_u64(1),
            FrField::from_u64(2),
            FrField::from_u64(3),
            FrField::from_u64(4),
        ];
        a.assert_is_in_set(&set, "should be in set");
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_is_in_set_fails() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let set = vec![
            FrField::from_u64(1),
            FrField::from_u64(2),
            FrField::from_u64(3),
        ];
        a.assert_is_in_set(&set, "should not be in set");
        // Circuit should fail
        assert!(check_circuit(&builder).is_err());
    }

    // ── Copy / fix tests ──────────────────────────────────────────────

    #[test]
    fn test_copy_as_new_witness() {
        let builder = make_builder();
        let val = Fr::random_element();
        let a = FrField::from_witness(builder.clone(), val);
        let b = FrField::copy_as_new_witness(builder.clone(), &a);
        assert_eq!(b.get_value(), val);
        assert!(!FrField::witness_indices_match(&a, &b));
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_fix_witness() {
        let builder = make_builder();
        let val = Fr::from(42u64);
        let mut a = FrField::from_witness(builder.clone(), val);
        a.fix_witness();
        assert_eq!(a.get_value(), val);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_set_public() {
        let builder = make_builder();
        let val = Fr::from(7u64);
        let a = FrField::from_witness(builder.clone(), val);
        let pi_idx = a.set_public();
        assert_eq!(pi_idx, 0);
        assert_eq!(builder.borrow().base.num_public_inputs(), 1);
    }

    #[test]
    fn test_convert_constant_to_fixed_witness() {
        let builder = make_builder();
        let val = Fr::from(99u64);
        let mut a = FrField::from_field(val);
        assert!(a.is_constant());
        a.convert_constant_to_fixed_witness(builder.clone());
        assert!(!a.is_constant());
        assert_eq!(a.get_value(), val);
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Arithmetic chain tests ────────────────────────────────────────

    #[test]
    fn test_add_mul_with_constants() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(1u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let c = FrField::from_witness(builder.clone(), Fr::from(3u64));
        let d = FrField::from_witness(builder.clone(), Fr::from(4u64));

        // Various combinations
        let result = (a.clone() + b.clone()) * (c.clone() + d.clone());
        assert_eq!(result.get_value(), Fr::from(21u64)); // (1+2)*(3+4) = 21

        // With constants
        let e = a.clone() * FrField::from_u64(5) + FrField::from_u64(3);
        assert_eq!(e.get_value(), Fr::from(8u64)); // 1*5 + 3

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_multiplicative_constant_regression() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(3u64));

        // Multiplying by constant should track multiplicative_constant
        let b = a.clone() * FrField::from_u64(5);
        assert_eq!(b.multiplicative_constant, Fr::from(5u64));
        assert_eq!(b.additive_constant, Fr::zero());

        // Adding constant should track additive_constant
        let c = a.clone() + FrField::from_u64(7);
        assert_eq!(c.additive_constant, Fr::from(7u64));
        assert_eq!(c.multiplicative_constant, Fr::one());

        // Combined: a * 5 + 7
        let d = b + FrField::from_u64(7);
        assert_eq!(d.get_value(), Fr::from(22u64)); // 3*5 + 7
    }

    #[test]
    fn test_arithmetic_chain() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(3u64));

        // Chain: (a + b) * a - b = (2+3)*2 - 3 = 7
        let result = (a.clone() + b.clone()) * a.clone() - b;
        assert_eq!(result.get_value(), Fr::from(7u64));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Fibonacci test ────────────────────────────────────────────────

    #[test]
    fn test_field_fibonacci() {
        let builder = make_builder();
        let mut a = FrField::from_witness(builder.clone(), Fr::from(1u64));
        let mut b = FrField::from_witness(builder.clone(), Fr::from(1u64));

        for _ in 0..16 {
            let c = a.clone() + b.clone();
            a = b;
            b = c;
        }
        // fib(18) = 4181 (starting from 1, 1)
        // Actually: fib sequence 1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584
        // After 16 additions: b = fib(18) = 4181 (wait let me recalculate)
        // f(1)=1,f(2)=1,f(3)=2,f(4)=3,f(5)=5,f(6)=8,f(7)=13,f(8)=21,f(9)=34,f(10)=55,f(11)=89,f(12)=144,f(13)=233,f(14)=377,f(15)=610,f(16)=987,f(17)=1597,f(18)=2584
        // Start a=f(1)=1, b=f(2)=1. After 16 iterations a=f(17)=1597, b=f(18)=2584
        assert_eq!(b.get_value(), Fr::from(2584u64));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Pythagorean test ──────────────────────────────────────────────

    #[test]
    fn test_field_pythagorean() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(3u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(4u64));
        let c = FrField::from_witness(builder.clone(), Fr::from(5u64));

        let lhs = a.sqr() + b.sqr(); // 9 + 16 = 25
        let rhs = c.sqr(); // 25
        lhs.assert_equal(&rhs, "pythagorean");
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Larger circuit test ───────────────────────────────────────────

    #[test]
    fn test_larger_circuit() {
        let builder = make_builder();

        // Create a chain of multiplications
        let mut acc = FrField::from_witness(builder.clone(), Fr::from(2u64));
        for i in 0..100 {
            let next = FrField::from_witness(builder.clone(), Fr::from((i + 3) as u64));
            acc = acc * next;
        }
        // Just verify the circuit is valid
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Range constraint test ─────────────────────────────────────────

    #[test]
    fn test_create_range_constraint() {
        let builder = make_builder();

        // 2-bit range: value must be < 4
        let a = FrField::from_witness(builder.clone(), Fr::from(3u64));
        a.create_range_constraint(2, "2-bit range");

        // 8-bit range: value must be < 256
        let b = FrField::from_witness(builder.clone(), Fr::from(255u64));
        b.create_range_constraint(8, "8-bit range");

        assert!(check_circuit(&builder).is_ok());
    }

    // ── Polynomial identity tests ─────────────────────────────────────

    #[test]
    fn test_evaluate_linear_identity() {
        let builder = make_builder();
        // a + b + c + d = 0
        let a = FrField::from_witness(builder.clone(), Fr::from(3u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let c = FrField::from_witness(builder.clone(), -Fr::from(8u64));
        let d = FrField::from_u64(0);
        FrField::evaluate_linear_identity(&a, &b, &c, &d, "linear identity");
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_evaluate_polynomial_identity() {
        let builder = make_builder();
        // a * b + c + d = 0 => 3 * 5 + (-15) + 0 = 0
        let a = FrField::from_witness(builder.clone(), Fr::from(3u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let c = FrField::from_witness(builder.clone(), -Fr::from(15u64));
        let d = FrField::from_u64(0);
        FrField::evaluate_polynomial_identity(&a, &b, &c, &d, "poly identity");
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Accumulate tests ──────────────────────────────────────────────

    #[test]
    fn test_accumulate() {
        let builder = make_builder();

        // Empty
        let result = FrField::accumulate(&[]);
        assert_eq!(result.get_value(), Fr::zero());

        // Single element
        let a = FrField::from_witness(builder.clone(), Fr::from(42u64));
        let result = FrField::accumulate(&[a]);
        assert_eq!(result.get_value(), Fr::from(42u64));

        // Multiple elements
        let mut elements = Vec::new();
        let mut expected = Fr::zero();
        for i in 1..=10 {
            let val = Fr::from(i as u64);
            elements.push(FrField::from_witness(builder.clone(), val));
            expected = expected + val;
        }
        let result = FrField::accumulate(&elements);
        assert_eq!(result.get_value(), expected);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_accumulate_with_constants() {
        let builder = make_builder();
        let elements = vec![
            FrField::from_witness(builder.clone(), Fr::from(1u64)),
            FrField::from_u64(2),
            FrField::from_witness(builder.clone(), Fr::from(3u64)),
            FrField::from_u64(4),
        ];
        let result = FrField::accumulate(&elements);
        assert_eq!(result.get_value(), Fr::from(10u64));
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_accumulate_gate_count() {
        let builder = make_builder();
        let mut elements = Vec::new();
        for i in 1..=9 {
            elements.push(FrField::from_witness(builder.clone(), Fr::from(i as u64)));
        }
        let before = num_gates(&builder);
        let _ = FrField::accumulate(&elements);
        let after = num_gates(&builder);
        // 9 elements / 3 = 3 gates (total is already normalized, no extra gate)
        assert_eq!(after - before, 3, "accumulate(9 witnesses) = 3 gates");
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Division by zero handling ─────────────────────────────────────

    #[test]
    fn test_div_zero_check() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let b = FrField::from_witness(builder.clone(), Fr::zero());

        // Division by zero should set circuit failure
        let _ = a / b;
        assert!(builder.borrow().base.failed());
    }
}
