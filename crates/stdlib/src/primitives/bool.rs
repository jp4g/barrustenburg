//! Circuit boolean type.
//!
//! Port of `barretenberg/stdlib/primitives/bool/bool.hpp` and `bool.cpp`.
//!
//! A `BoolT<P>` represents a circuit-level boolean with lazy negation.
//! Each element tracks a witness index, the native boolean value, and an
//! `witness_inverted` flag. The actual value is:
//!
//! ```text
//! value = witness_bool XOR witness_inverted
//!       = w + i - 2*i*w   (as a field element)
//! ```
//!
//! Negation flips `witness_inverted` without creating a gate, making NOT free.
//! Binary operations (AND, OR, XOR, ==) create a single arithmetic gate when
//! both operands are witnesses.

use std::ops::{BitAnd, BitOr, BitXor, Not};
use std::rc::Rc;

use bbrs_circuit_builder::gate_data::ArithmeticTriple;
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::witness::{BuilderRef, WitnessT, IS_CONSTANT};

// ════════════════════════════════════════════════════════════════════════
//  Helpers
// ════════════════════════════════════════════════════════════════════════

/// Convert a small signed integer to a field element.
/// Handles negative values by using field negation.
fn field_from_i32<P: FieldParams>(v: i32) -> Field<P> {
    if v >= 0 {
        Field::from(v as u64)
    } else {
        -Field::from((-v) as u64)
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Context validation (for BoolT)
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
                "Bool elements belong to different builders"
            );
            Some(c1.clone())
        }
    }
}

// ════════════════════════════════════════════════════════════════════════
//  BoolT
// ════════════════════════════════════════════════════════════════════════

/// A circuit boolean element with lazy negation.
///
/// Port of C++ `bool_t<Builder>`.
pub struct BoolT<P: FieldParams> {
    pub(crate) context: Option<BuilderRef<P>>,
    pub(crate) witness_bool: bool,
    pub(crate) witness_inverted: bool,
    pub(crate) witness_index: u32,
}

impl<P: FieldParams> Clone for BoolT<P> {
    fn clone(&self) -> Self {
        Self {
            context: self.context.clone(),
            witness_bool: self.witness_bool,
            witness_inverted: self.witness_inverted,
            witness_index: self.witness_index,
        }
    }
}

impl<P: FieldParams> BoolT<P> {
    // ════════════════════════════════════════════════════════════════════
    //  Constructors
    // ════════════════════════════════════════════════════════════════════

    /// Create a constant `BoolT` from a native `bool` value (no builder context).
    pub fn from_constant(value: bool) -> Self {
        Self {
            context: None,
            witness_bool: value,
            witness_inverted: false,
            witness_index: IS_CONSTANT,
        }
    }

    /// Create a constant false `BoolT` with a builder context.
    pub fn with_context(ctx: BuilderRef<P>) -> Self {
        Self {
            context: Some(ctx),
            witness_bool: false,
            witness_inverted: false,
            witness_index: IS_CONSTANT,
        }
    }

    /// Create a constant `BoolT` from a native `bool` with a builder context.
    pub fn constant_with_context(ctx: BuilderRef<P>, value: bool) -> Self {
        Self {
            context: Some(ctx),
            witness_bool: value,
            witness_inverted: false,
            witness_index: IS_CONSTANT,
        }
    }

    /// Create a `BoolT` from a witness. The witness value is constrained
    /// to be 0 or 1 via a bool gate (`x^2 = x`).
    ///
    /// Panics if the witness value is not 0 or 1.
    pub fn from_witness(value: &WitnessT<P>) -> Self {
        assert!(
            value.witness == Field::zero() || value.witness == Field::one(),
            "bool_t: witness value is not 0 or 1"
        );
        let witness_bool = value.witness == Field::one();
        value
            .context
            .borrow_mut()
            .create_bool_gate(value.witness_index);
        Self {
            context: Some(value.context.clone()),
            witness_bool,
            witness_inverted: false,
            witness_index: value.witness_index,
        }
    }

    /// Create a `BoolT` from a witness index that is **known** to already
    /// be constrained as boolean. No constraint is added.
    ///
    /// # Safety (logical)
    /// The caller must ensure the witness at `witness_index` is already
    /// constrained to be 0 or 1. Only an out-of-circuit assertion is performed.
    pub fn from_witness_index_unsafe(ctx: BuilderRef<P>, witness_index: u32) -> Self {
        assert_ne!(witness_index, IS_CONSTANT);
        let value = ctx.borrow().base.get_variable(witness_index);
        let v_sq = value * value;
        assert!(
            v_sq - value == Field::zero(),
            "bool_t: creating a witness bool from a non-boolean value"
        );
        let witness_bool = value == Field::one();
        Self {
            context: Some(ctx),
            witness_bool,
            witness_inverted: false,
            witness_index,
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Value access
    // ════════════════════════════════════════════════════════════════════

    /// Get the actual boolean value (accounting for inversion).
    pub fn get_value(&self) -> bool {
        self.witness_bool ^ self.witness_inverted
    }

    /// Check if this element is a circuit constant.
    pub fn is_constant(&self) -> bool {
        self.witness_index == IS_CONSTANT
    }

    /// Check if the witness is inverted.
    pub fn is_inverted(&self) -> bool {
        if self.is_constant() {
            debug_assert!(!self.witness_inverted);
        }
        self.witness_inverted
    }

    /// Get the builder context (if any).
    pub fn get_context(&self) -> &Option<BuilderRef<P>> {
        &self.context
    }

    /// Get the witness index of the normalized boolean.
    /// This normalizes first to ensure the witness contains the actual value.
    pub fn get_witness_index(&self) -> u32 {
        self.normalize().witness_index
    }

    // ════════════════════════════════════════════════════════════════════
    //  Normalization
    // ════════════════════════════════════════════════════════════════════

    /// Normalize: if `witness_inverted` is true, create a gate to produce
    /// a new witness containing the actual boolean value.
    ///
    /// Returns a normalized copy (with `witness_inverted == false`).
    pub fn normalize(&self) -> Self {
        if self.is_constant() {
            debug_assert!(!self.witness_inverted);
            return self.clone();
        }
        if !self.witness_inverted {
            return self.clone();
        }

        let ctx = self.context.as_ref().expect("normalize requires context").clone();
        let mut builder = ctx.borrow_mut();

        let value = self.witness_bool ^ self.witness_inverted;
        let new_witness = builder
            .base
            .add_variable(Field::from(value as u64));

        let inverted = self.witness_inverted as i32;
        let q_l: Field<P> = field_from_i32(1 - 2 * inverted);
        let q_c: Field<P> = field_from_i32(inverted);
        let q_o = -Field::<P>::one();

        let zero_idx = builder.base.zero_idx();
        builder.create_arithmetic_gate(&ArithmeticTriple {
            a: self.witness_index,
            b: zero_idx,
            c: new_witness,
            q_m: Field::zero(),
            q_l,
            q_r: Field::zero(),
            q_o,
            q_c,
        });

        Self {
            context: self.context.clone(),
            witness_bool: value,
            witness_inverted: false,
            witness_index: new_witness,
        }
    }

    /// Fix witness value so prover cannot change it.
    pub fn fix_witness(&self) {
        assert!(!self.is_constant());
        let ctx = self.context.as_ref().expect("fix_witness requires context");
        ctx.borrow_mut()
            .fix_witness(self.witness_index, Field::from(self.get_value() as u64));
    }

    // ════════════════════════════════════════════════════════════════════
    //  Boolean operations
    // ════════════════════════════════════════════════════════════════════

    /// Bitwise AND: `self & other`.
    ///
    /// For two witnesses, creates one arithmetic gate.
    /// With constant operands, optimizes to no gates.
    pub fn and(&self, other: &BoolT<P>) -> BoolT<P> {
        let ctx = validate_two_contexts(&self.context, &other.context);
        let left = self.witness_inverted ^ self.witness_bool;
        let right = other.witness_inverted ^ other.witness_bool;
        let result_bool = left && right;

        if !self.is_constant() && !other.is_constant() {
            let ctx_ref = ctx.as_ref().expect("witnesses require context").clone();
            let mut builder = ctx_ref.borrow_mut();
            let value = Field::<P>::from(result_bool as u64);
            let result_index = builder.base.add_variable(value);

            let i_a = self.witness_inverted as i32;
            let i_b = other.witness_inverted as i32;

            // q_m*w_a*w_b + q_l*w_a + q_r*w_b + q_o*result + q_c = 0
            let q_m: Field<P> = field_from_i32(1 - 2 * i_b - 2 * i_a + 4 * i_a * i_b);
            let q_l: Field<P> = field_from_i32(i_b * (1 - 2 * i_a));
            let q_r: Field<P> = field_from_i32(i_a * (1 - 2 * i_b));
            let q_o = -Field::<P>::one();
            let q_c: Field<P> = field_from_i32(i_a * i_b);

            builder.create_arithmetic_gate(&ArithmeticTriple {
                a: self.witness_index,
                b: other.witness_index,
                c: result_index,
                q_m,
                q_l,
                q_r,
                q_o,
                q_c,
            });

            BoolT {
                context: self.context.clone(),
                witness_bool: result_bool,
                witness_inverted: false,
                witness_index: result_index,
            }
        } else if !self.is_constant() && other.is_constant() {
            debug_assert!(!other.witness_inverted);
            if other.witness_bool {
                self.clone()
            } else {
                other.clone()
            }
        } else if self.is_constant() && !other.is_constant() {
            debug_assert!(!self.witness_inverted);
            if self.witness_bool {
                other.clone()
            } else {
                self.clone()
            }
        } else {
            // both constant
            BoolT {
                context: ctx,
                witness_bool: result_bool,
                witness_inverted: false,
                witness_index: IS_CONSTANT,
            }
        }
    }

    /// Bitwise OR: `self | other`.
    pub fn or(&self, other: &BoolT<P>) -> BoolT<P> {
        let ctx = validate_two_contexts(&self.context, &other.context);
        let result_bool =
            (self.witness_bool ^ self.witness_inverted) | (other.witness_bool ^ other.witness_inverted);

        if !self.is_constant() && !other.is_constant() {
            let ctx_ref = ctx.as_ref().expect("witnesses require context").clone();
            let mut builder = ctx_ref.borrow_mut();
            let value = Field::<P>::from(result_bool as u64);
            let result_index = builder.base.add_variable(value);

            let i_a = self.witness_inverted as i32;
            let i_b = other.witness_inverted as i32;

            let q_m: Field<P> = field_from_i32(-(1 - 2 * i_b) * (1 - 2 * i_a));
            let q_l: Field<P> = field_from_i32((1 - 2 * i_a) * (1 - i_b));
            let q_r: Field<P> = field_from_i32((1 - i_a) * (1 - 2 * i_b));
            let q_o = -Field::<P>::one();
            let q_c: Field<P> = field_from_i32(i_a + i_b - i_a * i_b);

            builder.create_arithmetic_gate(&ArithmeticTriple {
                a: self.witness_index,
                b: other.witness_index,
                c: result_index,
                q_m,
                q_l,
                q_r,
                q_o,
                q_c,
            });

            BoolT {
                context: self.context.clone(),
                witness_bool: result_bool,
                witness_inverted: false,
                witness_index: result_index,
            }
        } else if !self.is_constant() && other.is_constant() {
            debug_assert!(!other.witness_inverted);
            if other.witness_bool {
                other.clone()
            } else {
                self.clone()
            }
        } else if self.is_constant() && !other.is_constant() {
            debug_assert!(!self.witness_inverted);
            if self.witness_bool {
                self.clone()
            } else {
                other.clone()
            }
        } else {
            BoolT {
                context: ctx,
                witness_bool: result_bool,
                witness_inverted: false,
                witness_index: IS_CONSTANT,
            }
        }
    }

    /// Bitwise XOR: `self ^ other`.
    pub fn xor(&self, other: &BoolT<P>) -> BoolT<P> {
        let ctx = validate_two_contexts(&self.context, &other.context);
        let result_bool =
            (self.witness_bool ^ self.witness_inverted) ^ (other.witness_bool ^ other.witness_inverted);

        if !self.is_constant() && !other.is_constant() {
            let ctx_ref = ctx.as_ref().expect("witnesses require context").clone();
            let mut builder = ctx_ref.borrow_mut();
            let value = Field::<P>::from(result_bool as u64);
            let result_index = builder.base.add_variable(value);

            let i_a = self.witness_inverted as i32;
            let i_b = other.witness_inverted as i32;
            let aux_prod = (1 - 2 * i_b) * (1 - 2 * i_a);

            let q_m: Field<P> = field_from_i32(-2 * aux_prod);
            let q_l: Field<P> = field_from_i32(aux_prod);
            let q_r: Field<P> = field_from_i32(aux_prod);
            let q_o = -Field::<P>::one();
            let q_c: Field<P> = field_from_i32(i_a + i_b - 2 * i_a * i_b);

            builder.create_arithmetic_gate(&ArithmeticTriple {
                a: self.witness_index,
                b: other.witness_index,
                c: result_index,
                q_m,
                q_l,
                q_r,
                q_o,
                q_c,
            });

            BoolT {
                context: self.context.clone(),
                witness_bool: result_bool,
                witness_inverted: false,
                witness_index: result_index,
            }
        } else if !self.is_constant() && other.is_constant() {
            debug_assert!(!other.witness_inverted);
            if other.witness_bool {
                self.negate()
            } else {
                self.clone()
            }
        } else if self.is_constant() && !other.is_constant() {
            debug_assert!(!self.witness_inverted);
            if self.witness_bool {
                other.negate()
            } else {
                other.clone()
            }
        } else {
            BoolT {
                context: ctx,
                witness_bool: result_bool,
                witness_inverted: false,
                witness_index: IS_CONSTANT,
            }
        }
    }

    /// Logical NOT: `!self`.
    ///
    /// For constants, flips the value. For witnesses, flips `witness_inverted`.
    /// No gate is created.
    pub fn negate(&self) -> BoolT<P> {
        let mut result = self.clone();
        if result.is_constant() {
            debug_assert!(!self.witness_inverted);
            result.witness_bool = !result.witness_bool;
        } else {
            result.witness_inverted = !result.witness_inverted;
        }
        result
    }

    /// Equality check: `self == other` (returns a `BoolT`).
    pub fn eq(&self, other: &BoolT<P>) -> BoolT<P> {
        let ctx = validate_two_contexts(&self.context, &other.context);
        let result_bool =
            (self.witness_bool ^ self.witness_inverted) == (other.witness_bool ^ other.witness_inverted);

        if !self.is_constant() && !other.is_constant() {
            let ctx_ref = ctx.as_ref().expect("witnesses require context").clone();
            let mut builder = ctx_ref.borrow_mut();
            let value = Field::<P>::from(result_bool as u64);
            let result_index = builder.base.add_variable(value);

            let i_a = self.witness_inverted as i32;
            let i_b = other.witness_inverted as i32;
            let aux_prod = (1 - 2 * i_b) * (1 - 2 * i_a);

            let q_m: Field<P> = field_from_i32(2 * aux_prod);
            let q_l: Field<P> = field_from_i32(-aux_prod);
            let q_r: Field<P> = field_from_i32(-aux_prod);
            let q_o = -Field::<P>::one();
            let q_c: Field<P> = field_from_i32(1 - i_a - i_b + 2 * i_a * i_b);

            builder.create_arithmetic_gate(&ArithmeticTriple {
                a: self.witness_index,
                b: other.witness_index,
                c: result_index,
                q_m,
                q_l,
                q_r,
                q_o,
                q_c,
            });

            BoolT {
                context: self.context.clone(),
                witness_bool: result_bool,
                witness_inverted: false,
                witness_index: result_index,
            }
        } else if !self.is_constant() && other.is_constant() {
            debug_assert!(!other.witness_inverted);
            if other.witness_bool {
                self.clone()
            } else {
                self.negate()
            }
        } else if self.is_constant() && !other.is_constant() {
            debug_assert!(!self.witness_inverted);
            if self.witness_bool {
                other.clone()
            } else {
                other.negate()
            }
        } else {
            BoolT {
                context: ctx,
                witness_bool: result_bool,
                witness_inverted: false,
                witness_index: IS_CONSTANT,
            }
        }
    }

    /// Inequality check: `self != other` (equivalent to XOR for booleans).
    pub fn neq(&self, other: &BoolT<P>) -> BoolT<P> {
        self.xor(other)
    }

    // ════════════════════════════════════════════════════════════════════
    //  Implication and assertions
    // ════════════════════════════════════════════════════════════════════

    /// Logical implication: `self => other` (equivalent to `!self || other`).
    pub fn implies(&self, other: &BoolT<P>) -> BoolT<P> {
        self.negate().or(other)
    }

    /// Bidirectional implication (iff): `self <=> other` (equivalent to `!(self ^ other)`).
    pub fn implies_both_ways(&self, other: &BoolT<P>) -> BoolT<P> {
        self.xor(other).negate()
    }

    /// Assert that `self => other` must hold, i.e. `(!self || other) == true`.
    pub fn must_imply(&self, other: &BoolT<P>, msg: &str) {
        self.implies(other).assert_equal_bool(true, msg);
    }

    /// Assert that this bool equals the given constant value.
    pub fn assert_equal_bool(&self, rhs_value: bool, msg: &str) {
        let rhs = BoolT::from_constant(rhs_value);
        self.assert_equal(&rhs, msg);
    }

    /// Assert that two booleans are equal.
    pub fn assert_equal(&self, rhs: &BoolT<P>, msg: &str) {
        let ctx = validate_two_contexts(&self.context, &rhs.context);

        if self.is_constant() && rhs.is_constant() {
            assert_eq!(
                self.get_value(),
                rhs.get_value(),
                "bool_t::assert_equal: constants are not equal"
            );
        } else if self.is_constant() {
            debug_assert!(!self.witness_inverted);
            let lhs_value = if rhs.witness_inverted {
                !self.witness_bool
            } else {
                self.witness_bool
            };
            let ctx_ref = ctx.as_ref().expect("assert_equal requires context");
            ctx_ref.borrow_mut().assert_equal_constant(
                rhs.witness_index,
                Field::from(lhs_value as u64),
                msg,
            );
        } else if rhs.is_constant() {
            debug_assert!(!rhs.witness_inverted);
            let rhs_value = if self.witness_inverted {
                !rhs.witness_bool
            } else {
                rhs.witness_bool
            };
            let ctx_ref = ctx.as_ref().expect("assert_equal requires context");
            ctx_ref.borrow_mut().assert_equal_constant(
                self.witness_index,
                Field::from(rhs_value as u64),
                msg,
            );
        } else {
            // Both are witnesses. Need to normalize if exactly one is inverted.
            let mut left = self.clone();
            let mut right = rhs.clone();
            if self.witness_inverted ^ rhs.witness_inverted {
                left = left.normalize();
                right = right.normalize();
            }
            let ctx_ref = ctx.as_ref().expect("assert_equal requires context");
            ctx_ref
                .borrow_mut()
                .base
                .assert_equal(left.witness_index, right.witness_index, msg);
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Conditional assignment
    // ════════════════════════════════════════════════════════════════════

    /// Conditional assignment: if `predicate` is true, return `lhs`, else `rhs`.
    /// Always returns a normalized result.
    pub fn conditional_assign(predicate: &BoolT<P>, lhs: &BoolT<P>, rhs: &BoolT<P>) -> BoolT<P> {
        if predicate.is_constant() {
            let result = if predicate.get_value() {
                lhs.clone()
            } else {
                rhs.clone()
            };
            return result.normalize();
        }

        let same = lhs.witness_index == rhs.witness_index;
        let witness_same = same && !lhs.is_constant() && (lhs.witness_inverted == rhs.witness_inverted);
        let const_same = same && lhs.is_constant() && (lhs.witness_bool == rhs.witness_bool);
        if witness_same || const_same {
            return lhs.normalize();
        }

        // (predicate && lhs) || (!predicate && rhs)
        predicate.and(lhs).or(&predicate.negate().and(rhs)).normalize()
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Operator implementations
// ════════════════════════════════════════════════════════════════════════

impl<P: FieldParams> BitAnd for &BoolT<P> {
    type Output = BoolT<P>;
    fn bitand(self, rhs: Self) -> BoolT<P> {
        self.and(rhs)
    }
}

impl<P: FieldParams> BitAnd for BoolT<P> {
    type Output = BoolT<P>;
    fn bitand(self, rhs: Self) -> BoolT<P> {
        self.and(&rhs)
    }
}

impl<P: FieldParams> BitOr for &BoolT<P> {
    type Output = BoolT<P>;
    fn bitor(self, rhs: Self) -> BoolT<P> {
        self.or(rhs)
    }
}

impl<P: FieldParams> BitOr for BoolT<P> {
    type Output = BoolT<P>;
    fn bitor(self, rhs: Self) -> BoolT<P> {
        self.or(&rhs)
    }
}

impl<P: FieldParams> BitXor for &BoolT<P> {
    type Output = BoolT<P>;
    fn bitxor(self, rhs: Self) -> BoolT<P> {
        self.xor(rhs)
    }
}

impl<P: FieldParams> BitXor for BoolT<P> {
    type Output = BoolT<P>;
    fn bitxor(self, rhs: Self) -> BoolT<P> {
        self.xor(&rhs)
    }
}

impl<P: FieldParams> Not for &BoolT<P> {
    type Output = BoolT<P>;
    fn not(self) -> BoolT<P> {
        BoolT::negate(self)
    }
}

impl<P: FieldParams> Not for BoolT<P> {
    type Output = BoolT<P>;
    fn not(self) -> BoolT<P> {
        BoolT::negate(&self)
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
    use bbrs_ecc::curves::bn254::Bn254FrParams;
    use std::cell::RefCell;

    type B = BoolT<Bn254FrParams>;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    fn check_circuit(builder: &BuilderRef<Bn254FrParams>) -> bool {
        UltraCircuitChecker::check(&mut builder.borrow_mut()).is_ok()
    }

    fn num_gates(builder: &BuilderRef<Bn254FrParams>) -> usize {
        builder.borrow().base.num_gates()
    }

    /// All 8 input combinations: (is_const, value, is_inverted).
    struct BoolInput {
        is_const: bool,
        value: bool,
        is_inverted: bool,
    }

    fn all_inputs() -> Vec<BoolInput> {
        let mut inputs = Vec::new();
        for is_const in [false, true] {
            for value in [false, true] {
                for is_inverted in [false, true] {
                    inputs.push(BoolInput {
                        is_const,
                        value,
                        is_inverted,
                    });
                }
            }
        }
        inputs
    }

    fn create_bool_ct(input: &BoolInput, builder: &BuilderRef<Bn254FrParams>) -> B {
        let b = if input.is_const {
            B::from_constant(input.value)
        } else {
            B::from_witness(&WitnessT::from_bool(builder.clone(), input.value))
        };
        if input.is_inverted { b.negate() } else { b }
    }

    /// Generic binary operation test: exhaustively tests all 64 input combinations.
    fn test_binary_op(
        op_name: &str,
        op: impl Fn(&B, &B) -> B,
        expected_op: impl Fn(bool, bool) -> bool,
    ) {
        let builder = make_builder();
        for lhs in all_inputs().iter() {
            for rhs in all_inputs().iter() {
                let a = create_bool_ct(lhs, &builder);
                let b = create_bool_ct(rhs, &builder);

                let gates_before = num_gates(&builder);
                let c = op(&a, &b);

                let expected = expected_op(
                    lhs.value ^ lhs.is_inverted,
                    rhs.value ^ rhs.is_inverted,
                );
                assert_eq!(
                    c.get_value(),
                    expected,
                    "Failed on {} with lhs={{const={}, val={}, inv={}}}, rhs={{const={}, val={}, inv={}}}",
                    op_name,
                    lhs.is_const, lhs.value, lhs.is_inverted,
                    rhs.is_const, rhs.value, rhs.is_inverted,
                );

                if a.is_constant() && b.is_constant() {
                    assert!(c.is_constant());
                }
                if !a.is_constant() && !b.is_constant() {
                    assert!(!c.is_constant());
                }

                let gate_diff = num_gates(&builder) - gates_before;
                let expected_gates = if !a.is_constant() && !b.is_constant() {
                    1
                } else {
                    0
                };
                assert_eq!(
                    gate_diff, expected_gates,
                    "{}: expected {} gates, got {} gates",
                    op_name, expected_gates, gate_diff,
                );
            }
        }
        assert!(check_circuit(&builder));
    }

    // ─── Test 1: Construct from constant bool ───────────────────────────

    #[test]
    fn test_construct_from_const_bool() {
        let builder = make_builder();
        let gates_before = num_gates(&builder);
        let a_true = B::from_constant(true);
        let a_false = B::from_constant(false);
        assert!(a_true.get_value());
        assert!(!a_false.get_value());
        assert!(a_true.is_constant() && a_false.is_constant());
        assert!(!a_true.is_inverted() && !a_false.is_inverted());
        assert_eq!(num_gates(&builder) - gates_before, 0);
    }

    // ─── Test 2: Construct from witness index (unsafe) ──────────────────

    #[test]
    fn test_construct_from_witness_index() {
        let builder = make_builder();
        let gates_before = num_gates(&builder);
        let idx_zero = builder.borrow_mut().base.add_variable(Field::zero());
        let idx_one = builder.borrow_mut().base.add_variable(Field::one());

        let b_zero = B::from_witness_index_unsafe(builder.clone(), idx_zero);
        assert!(!b_zero.get_value());

        let b_one = B::from_witness_index_unsafe(builder.clone(), idx_one);
        assert!(b_one.get_value());

        assert_eq!(num_gates(&builder) - gates_before, 0);
    }

    #[test]
    #[should_panic(expected = "bool_t: creating a witness bool from a non-boolean value")]
    fn test_construct_from_witness_index_non_bool() {
        let builder = make_builder();
        let idx_bad = builder.borrow_mut().base.add_variable(Field::from(15u64));
        let _b = B::from_witness_index_unsafe(builder.clone(), idx_bad);
    }

    // ─── Test 3: Construct from witness ─────────────────────────────────

    #[test]
    fn test_construct_from_witness() {
        let builder = make_builder();
        let gates_before = num_gates(&builder);

        let a_true = B::from_witness(&WitnessT::from_bool(builder.clone(), true));
        let a_false = B::from_witness(&WitnessT::from_bool(builder.clone(), false));
        assert!(a_true.get_value());
        assert!(!a_false.get_value());
        assert!(!a_true.is_constant() && !a_false.is_constant());
        assert!(!a_true.is_inverted() && !a_false.is_inverted());
        // Each witness bool creates 1 bool gate
        assert_eq!(num_gates(&builder) - gates_before, 2);
        assert!(check_circuit(&builder));
    }

    #[test]
    #[should_panic(expected = "bool_t: witness value is not 0 or 1")]
    fn test_construct_from_witness_non_bool() {
        let builder = make_builder();
        let w = WitnessT::new(builder.clone(), Field::from(5u64));
        let _b = B::from_witness(&w);
    }

    // ─── Test 4: AND ────────────────────────────────────────────────────

    #[test]
    fn test_and() {
        test_binary_op("AND", |a, b| a.and(b), |a, b| a && b);
    }

    // ─── Test 5: OR ─────────────────────────────────────────────────────

    #[test]
    fn test_or() {
        test_binary_op("OR", |a, b| a.or(b), |a, b| a || b);
    }

    // ─── Test 6: XOR ────────────────────────────────────────────────────

    #[test]
    fn test_xor() {
        test_binary_op("XOR", |a, b| a.xor(b), |a, b| a ^ b);
    }

    // ─── Test 7: EQ ─────────────────────────────────────────────────────

    #[test]
    fn test_eq() {
        test_binary_op("EQ", |a, b| a.eq(b), |a, b| a == b);
    }

    // ─── Test 8: NEQ ────────────────────────────────────────────────────

    #[test]
    fn test_neq() {
        test_binary_op("NEQ", |a, b| a.neq(b), |a, b| a != b);
    }

    // ─── Test 9: Implies ────────────────────────────────────────────────

    #[test]
    fn test_implies() {
        test_binary_op("IMPLIES", |a, b| a.implies(b), |a, b| !a || b);
    }

    // ─── Test 10: Implies both ways ─────────────────────────────────────

    #[test]
    fn test_implies_both_ways() {
        test_binary_op(
            "IFF",
            |a, b| a.implies_both_ways(b),
            |a, b| !(a ^ b),
        );
    }

    // ─── Test 11: NOT ───────────────────────────────────────────────────

    #[test]
    fn test_not() {
        let builder = make_builder();
        for input in all_inputs().iter() {
            let a = create_bool_ct(input, &builder);
            let gates_before = num_gates(&builder);
            let b = a.negate();
            // NOT should create zero gates
            assert_eq!(num_gates(&builder) - gates_before, 0);
            assert_eq!(b.get_value(), !(input.value ^ input.is_inverted));
        }
        assert!(check_circuit(&builder));
    }

    // ─── Test 12: Normalize ─────────────────────────────────────────────

    #[test]
    fn test_normalize() {
        for input in all_inputs().iter() {
            let builder = make_builder();
            let a = create_bool_ct(input, &builder);
            let gates_before = num_gates(&builder);

            let c = a.normalize();
            assert_eq!(c.get_value(), a.get_value());
            assert!(!c.is_inverted());

            let gate_diff = num_gates(&builder) - gates_before;
            let expected_gates = if !a.is_constant() && input.is_inverted {
                1
            } else {
                0
            };
            assert_eq!(gate_diff, expected_gates);
            assert!(check_circuit(&builder));
        }
    }

    // ─── Test 13: Assert equal ──────────────────────────────────────────

    #[test]
    fn test_assert_equal() {
        for lhs in all_inputs().iter() {
            for rhs in all_inputs().iter() {
                let builder = make_builder();
                let a = create_bool_ct(lhs, &builder);
                let b = create_bool_ct(rhs, &builder);

                let values_match = a.get_value() == b.get_value();

                if a.is_constant() && b.is_constant() {
                    if !values_match {
                        // Should panic
                        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            a.assert_equal(&b, "test");
                        }));
                        assert!(result.is_err());
                    }
                    // Skip checking circuit for const-const mismatch
                } else if !a.is_constant() && !b.is_constant() {
                    a.assert_equal(&b, "test");
                    // CircuitChecker doesn't verify permutation for copy constraints,
                    // so builder.failed() indicates the failure
                    assert_eq!(builder.borrow().base.failed(), !values_match);
                } else {
                    a.assert_equal(&b, "test");
                    assert_eq!(check_circuit(&builder), values_match);
                }
            }
        }
    }

    // ─── Test 14: Must imply ────────────────────────────────────────────

    #[test]
    fn test_must_imply() {
        for lhs in all_inputs().iter() {
            for rhs in all_inputs().iter() {
                let builder = make_builder();
                let a = create_bool_ct(lhs, &builder);
                let b = create_bool_ct(rhs, &builder);

                let a_val = lhs.value ^ lhs.is_inverted;
                let b_val = rhs.value ^ rhs.is_inverted;
                let implication_holds = !a_val || b_val;

                if a.is_constant() && b.is_constant() && !implication_holds {
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        a.must_imply(&b, "test");
                    }));
                    assert!(result.is_err());
                } else {
                    a.must_imply(&b, "test must_imply");
                    // Check circuit validity
                    let circuit_valid = check_circuit(&builder);
                    let builder_failed = builder.borrow().base.failed();
                    // Circuit should pass iff implication holds (or builder marks failure)
                    assert!(
                        (circuit_valid && implication_holds)
                            || (!circuit_valid && !implication_holds)
                            || builder_failed,
                        "must_imply: a={}, b={}, expected={}, circuit_valid={}, builder_failed={}",
                        a_val, b_val, implication_holds, circuit_valid, builder_failed,
                    );
                }
            }
        }
    }

    // ─── Test 15: Conditional assign ────────────────────────────────────

    #[test]
    fn test_conditional_assign() {
        for pred_input in all_inputs().iter() {
            for lhs in all_inputs().iter() {
                for rhs in all_inputs().iter() {
                    let builder = make_builder();
                    let condition = create_bool_ct(pred_input, &builder);
                    let a = create_bool_ct(lhs, &builder);
                    let b = create_bool_ct(rhs, &builder);

                    let result = B::conditional_assign(&condition, &a, &b);

                    let expected = if condition.get_value() {
                        a.get_value()
                    } else {
                        b.get_value()
                    };
                    assert_eq!(result.get_value(), expected);
                    assert!(!result.is_inverted());
                    assert!(check_circuit(&builder));
                }
            }
        }
    }

    // ─── Test 16: Operator overloads ────────────────────────────────────

    #[test]
    fn test_operator_overloads() {
        let builder = make_builder();
        let a = B::from_witness(&WitnessT::from_bool(builder.clone(), true));
        let b = B::from_witness(&WitnessT::from_bool(builder.clone(), false));

        // BitAnd
        let c = &a & &b;
        assert!(!c.get_value());

        // BitOr
        let c = &a | &b;
        assert!(c.get_value());

        // BitXor
        let c = &a ^ &b;
        assert!(c.get_value());

        // Not
        let c = !&a;
        assert!(!c.get_value());

        // Owned operators
        let a2 = a.clone();
        let b2 = b.clone();
        let c = a2 & b2;
        assert!(!c.get_value());

        assert!(check_circuit(&builder));
    }

    // ─── Test 17: Simple proof (compound operations) ────────────────────

    #[test]
    fn test_simple_proof() {
        let builder = make_builder();
        let a = B::from_witness(&WitnessT::from_bool(builder.clone(), true));
        let b = B::from_witness(&WitnessT::from_bool(builder.clone(), false));

        // a ^ b should be true
        let c = a.xor(&b);
        assert!(c.get_value());

        // a & b should be false
        let d = a.and(&b);
        assert!(!d.get_value());

        // a | b should be true
        let e = a.or(&b);
        assert!(e.get_value());

        // a implies a should be true
        let f = a.implies(&a);
        assert!(f.get_value());

        // a implies_both_ways a should be true
        let g = a.implies_both_ways(&a);
        assert!(g.get_value());

        // Chain: (a & (b | c)) ^ (d.implies(f)) should equal ((a&b) | (a&c)) ^ (!d | f)
        // This is a tautology
        let lhs = a.and(&(b.or(&c))).xor(&d.implies(&f));
        let rhs = (a.and(&b).or(&a.and(&c))).xor(&d.negate().or(&f));
        let equivalent = lhs.implies_both_ways(&rhs);
        assert!(equivalent.get_value());

        assert!(check_circuit(&builder));
    }

    // ─── Test 18: Boolean algebra identity ──────────────────────────────

    #[test]
    fn test_boolean_algebra_identities() {
        // Test De Morgan's laws and distributive properties
        for a_val in [false, true] {
            for b_val in [false, true] {
                let builder = make_builder();
                let a = B::from_witness(&WitnessT::from_bool(builder.clone(), a_val));
                let b = B::from_witness(&WitnessT::from_bool(builder.clone(), b_val));

                // De Morgan: !(a & b) == (!a | !b)
                let lhs = a.and(&b).negate();
                let rhs = a.negate().or(&b.negate());
                assert_eq!(lhs.get_value(), rhs.get_value());

                // De Morgan: !(a | b) == (!a & !b)
                let lhs = a.or(&b).negate();
                let rhs = a.negate().and(&b.negate());
                assert_eq!(lhs.get_value(), rhs.get_value());

                // Distributive: a & (a | b) == a
                let lhs = a.and(&a.or(&b));
                assert_eq!(lhs.get_value(), a.get_value());

                // Complement: a ^ a == false
                let lhs = a.xor(&a);
                assert!(!lhs.get_value());

                // Identity: a | false == a
                let f = B::from_constant(false);
                let lhs = a.or(&f);
                assert_eq!(lhs.get_value(), a.get_value());

                // Identity: a & true == a
                let t = B::from_constant(true);
                let lhs = a.and(&t);
                assert_eq!(lhs.get_value(), a.get_value());

                assert!(check_circuit(&builder));
            }
        }
    }
}
