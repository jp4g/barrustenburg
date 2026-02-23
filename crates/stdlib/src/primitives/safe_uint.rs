//! Circuit safe unsigned integer type.
//!
//! Port of `barretenberg/stdlib/primitives/safe_uint/safe_uint.hpp` and `safe_uint.cpp`.
//!
//! A `SafeUintT<P>` wraps a `FieldT<P>` and tracks a `current_max` (U256) to
//! ensure positive integer operations never overflow the field modulus.
//!
//! Unlike the regular uint circuit types, safe_uint operations model positive
//! integer arithmetic (not modular). Every operation checks that the tracked
//! maximum value stays within `[0, modulus - 1]`.

use std::ops::{Add, AddAssign, Mul, MulAssign};

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_numeric::uint256::U256Ext;
use bbrs_numeric::uintx::U512Ext;
use bbrs_numeric::{U256, U512};
use super::bool::BoolT;
use super::field::FieldT;
use super::witness::{BuilderRef, WitnessT};

// ════════════════════════════════════════════════════════════════════════
//  Helpers
// ════════════════════════════════════════════════════════════════════════

/// Convert a `Field<P>` (Montgomery form) to a `U256` integer value.
fn field_to_u256<P: FieldParams>(f: &Field<P>) -> U256 {
    U256::from_limbs(f.from_montgomery_form().data)
}

/// Convert a `U256` integer value to a `Field<P>` (Montgomery form).
fn u256_to_field<P: FieldParams>(v: &U256) -> Field<P> {
    Field::from_limbs(v.limbs())
}

/// Convert a `BoolT` to a `FieldT` (0 or 1 as a field element).
fn bool_to_field<P: FieldParams>(b: &BoolT<P>) -> FieldT<P> {
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

/// Compute `field_t::is_zero()` — returns a `BoolT` constrained to be
/// true iff the field element is zero. Creates 2 arithmetic gates.
fn field_is_zero<P: FieldParams>(f: &FieldT<P>) -> BoolT<P> {
    use bbrs_circuit_builder::gate_data::ArithmeticTriple;

    if f.is_constant() {
        return BoolT::from_constant(f.get_value() == Field::zero());
    }

    let ctx = f.context.as_ref().unwrap().clone();
    let value = f.get_value();
    let is_zero_val = value == Field::zero();
    let inverse_val = if value.is_zero() {
        Field::zero()
    } else {
        value.invert()
    };

    // Create witnesses
    let is_zero_witness = WitnessT::from_bool(ctx.clone(), is_zero_val);
    let is_zero_bool = BoolT::from_witness(&is_zero_witness);
    let inverse_witness = WitnessT::new(ctx.clone(), inverse_val);

    // Normalize field element for constraint
    let normalized = f.normalize();

    let mut builder = ctx.borrow_mut();
    let zero_idx = builder.base.zero_idx();

    // Constraint 1: value * inverse + is_zero - 1 = 0
    // i.e., value * inverse = 1 - is_zero
    builder.create_big_mul_add_gate(
        &bbrs_circuit_builder::gate_data::MulQuad {
            a: normalized.witness_index,
            b: inverse_witness.witness_index,
            c: is_zero_bool.witness_index,
            d: zero_idx,
            mul_scaling: Field::one(),
            a_scaling: Field::zero(),
            b_scaling: Field::zero(),
            c_scaling: Field::one(),
            d_scaling: Field::zero(),
            const_scaling: -Field::one(),
        },
        false,
    );

    // Constraint 2: value * is_zero = 0
    builder.create_arithmetic_gate(&ArithmeticTriple {
        a: normalized.witness_index,
        b: is_zero_bool.witness_index,
        c: zero_idx,
        q_m: Field::one(),
        q_l: Field::zero(),
        q_r: Field::zero(),
        q_o: Field::zero(),
        q_c: Field::zero(),
    });

    drop(builder);
    is_zero_bool
}

// ════════════════════════════════════════════════════════════════════════
//  SafeUintT
// ════════════════════════════════════════════════════════════════════════

/// A circuit positive integer element with overflow tracking.
///
/// Port of C++ `safe_uint_t<Builder>`.
pub struct SafeUintT<P: FieldParams> {
    pub value: FieldT<P>,
    pub current_max: U256,
}

impl<P: FieldParams> Clone for SafeUintT<P> {
    fn clone(&self) -> Self {
        Self {
            value: self.value.clone(),
            current_max: self.current_max,
        }
    }
}

impl<P: FieldParams> SafeUintT<P> {
    /// Maximum bit number for range constraints (modulus MSB position).
    pub fn max_bit_num() -> usize {
        U256::from_limbs(P::MODULUS).get_msb() as usize
    }

    /// Maximum representable value: `modulus - 1`.
    pub fn max_value() -> U256 {
        U256::from_limbs(P::MODULUS).wrapping_sub(&U256::ONE)
    }

    // ════════════════════════════════════════════════════════════════════
    //  Private constructor
    // ════════════════════════════════════════════════════════════════════

    /// Create without range check (internal use only).
    /// Verifies `current_max <= MAX_VALUE` at build time.
    fn new_unsafe(value: FieldT<P>, current_max: U256) -> Self {
        assert!(
            current_max <= Self::max_value(),
            "exceeded modulus in safe_uint class"
        );
        Self { value, current_max }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Constructors
    // ════════════════════════════════════════════════════════════════════

    /// Default: value = 0, max = 0.
    pub fn default() -> Self {
        Self {
            value: FieldT::from_u64(0),
            current_max: U256::ZERO,
        }
    }

    /// Create from a field element with a range constraint of `bit_num` bits.
    ///
    /// The value is constrained to `[0, 2^bit_num - 1]`.
    pub fn from_field_with_range(value: FieldT<P>, bit_num: usize, description: &str) -> Self {
        assert!(
            bit_num <= Self::max_bit_num(),
            "bit_num exceeds max for safe_uint"
        );
        value.create_range_constraint(
            bit_num,
            &format!("safe_uint_t range constraint failure: {}", description),
        );
        let current_max = if bit_num == 0 {
            U256::ZERO
        } else {
            U256::ONE
                .wrapping_shl_vartime(bit_num as u32)
                .wrapping_sub(&U256::ONE)
        };
        Self { value, current_max }
    }

    /// Create from a constant field value.
    /// The max is set to the value itself (tight bound).
    pub fn from_constant_field(const_value: Field<P>) -> Self {
        let max = field_to_u256(&const_value);
        Self {
            value: FieldT::from_field(const_value),
            current_max: max,
        }
    }

    /// Create from a constant U256 value.
    pub fn from_constant_u256(const_value: U256) -> Self {
        let field_val = u256_to_field::<P>(&const_value);
        Self {
            value: FieldT::from_field(field_val),
            current_max: const_value,
        }
    }

    /// Create from a constant u32 value.
    pub fn from_constant_u32(const_value: u32) -> Self {
        Self::from_constant_u256(U256::from_limbs([const_value as u64, 0, 0, 0]))
    }

    /// Create from a `BoolT`. Max is 1 (no redundant range check needed).
    pub fn from_bool(other: &BoolT<P>) -> Self {
        Self {
            value: bool_to_field(other),
            current_max: U256::ONE,
        }
    }

    /// Create a constant witness constrained to equal the given value.
    pub fn create_constant_witness(ctx: BuilderRef<P>, value: Field<P>) -> Self {
        let _ = WitnessT::create_constant_witness(ctx, value);
        let max = field_to_u256(&value);
        Self::new_unsafe(FieldT::from_field(value), max)
    }

    /// Create from a witness index. The value is range-constrained to `bit_num` bits.
    pub fn from_witness_index(ctx: BuilderRef<P>, witness_index: u32, bit_num: usize) -> Self {
        let field = FieldT::from_witness_index(ctx, witness_index);
        Self::from_field_with_range(field, bit_num, "from_witness_index")
    }

    // ════════════════════════════════════════════════════════════════════
    //  Value access
    // ════════════════════════════════════════════════════════════════════

    /// Get the actual value as a native field element.
    pub fn get_value(&self) -> Field<P> {
        self.value.get_value()
    }

    /// Get the builder context (if any).
    pub fn get_context(&self) -> &Option<BuilderRef<P>> {
        self.value.get_context()
    }

    /// Check if this element is a circuit constant.
    pub fn is_constant(&self) -> bool {
        self.value.is_constant()
    }

    /// Normalize: return a copy where `multiplicative_constant == 1`
    /// and `additive_constant == 0`.
    pub fn normalize(&self) -> Self {
        Self {
            value: self.value.normalize(),
            current_max: self.current_max,
        }
    }

    /// Make this value a public input.
    pub fn set_public(&self) -> u32 {
        self.value.set_public()
    }

    /// Get the normalized witness index.
    pub fn get_witness_index(&self) -> u32 {
        self.value.get_witness_index()
    }

    // ════════════════════════════════════════════════════════════════════
    //  Arithmetic: addition
    // ════════════════════════════════════════════════════════════════════

    fn add_impl(&self, other: &SafeUintT<P>) -> SafeUintT<P> {
        let new_max = self.current_max.wrapping_add(&other.current_max);
        // If wrapping occurred, new_max < either operand
        if new_max < self.current_max || new_max < other.current_max {
            panic!("exceeded modulus in safe_uint class");
        }
        Self::new_unsafe(
            self.value.clone() + other.value.clone(),
            new_max,
        )
    }

    /// Efficiently compute `self + add_a + add_b`.
    pub fn add_two(&self, add_a: &SafeUintT<P>, add_b: &SafeUintT<P>) -> SafeUintT<P> {
        let sum = self
            .current_max
            .wrapping_add(&add_a.current_max)
            .wrapping_add(&add_b.current_max);
        assert!(
            sum <= Self::max_value(),
            "Exceeded modulus in add_two"
        );
        let new_val = self.value.add_two(&add_a.value, &add_b.value);
        Self::new_unsafe(new_val, sum)
    }

    /// Efficiently compute `self * to_mul + to_add`.
    pub fn madd(&self, to_mul: &SafeUintT<P>, to_add: &SafeUintT<P>) -> SafeUintT<P> {
        let wide_product = U512::from_lo_hi(self.current_max, U256::ZERO)
            .wrapping_mul(&U512::from_lo_hi(to_mul.current_max, U256::ZERO));
        let wide_sum = wide_product
            .wrapping_add(&U512::from_lo_hi(to_add.current_max, U256::ZERO));
        assert!(
            wide_sum <= U512::from_lo_hi(Self::max_value(), U256::ZERO),
            "Exceeded modulus in madd"
        );
        let new_val = self.value.madd(&to_mul.value, &to_add.value);
        Self::new_unsafe(new_val, wide_sum.lo())
    }

    // ════════════════════════════════════════════════════════════════════
    //  Arithmetic: subtraction
    // ════════════════════════════════════════════════════════════════════

    /// Subtraction with a pre-determined bound on the difference size.
    pub fn subtract(
        &self,
        other: &SafeUintT<P>,
        difference_bit_size: usize,
        description: &str,
    ) -> SafeUintT<P> {
        assert!(
            difference_bit_size <= Self::max_bit_num(),
            "difference bit size exceeds max"
        );
        assert!(
            !(self.value.is_constant() && other.value.is_constant()),
            "subtract requires at least one witness"
        );

        let difference_val = self.value.clone() - other.value.clone();
        let difference = SafeUintT::from_field_with_range(
            difference_val,
            difference_bit_size,
            &format!("subtract: {}", description),
        );

        if difference.current_max.wrapping_add(&other.current_max) > Self::max_value() {
            panic!("maximum value exceeded in safe_uint subtract");
        }
        difference
    }

    fn sub_impl(&self, other: &SafeUintT<P>) -> SafeUintT<P> {
        // If both are constants and underflow, panic
        if self.value.is_constant() && other.value.is_constant() {
            let self_u256 = field_to_u256(&self.value.get_value());
            let other_u256 = field_to_u256(&other.value.get_value());
            assert!(
                self_u256 >= other_u256,
                "underflow in safe_uint constant subtraction"
            );
        }

        let difference_val = self.value.clone() - other.value.clone();
        let bit_size = (self.current_max.get_msb() + 1) as usize;
        let difference = SafeUintT::from_field_with_range(
            difference_val,
            bit_size,
            "- operator",
        );

        if difference.current_max.wrapping_add(&other.current_max) > Self::max_value() {
            panic!("maximum value exceeded in safe_uint minus operator");
        }
        difference
    }

    // ════════════════════════════════════════════════════════════════════
    //  Arithmetic: multiplication
    // ════════════════════════════════════════════════════════════════════

    fn mul_impl(&self, other: &SafeUintT<P>) -> SafeUintT<P> {
        let wide = self.current_max.widening_mul(&other.current_max);
        let (lo, hi) = wide.split();
        assert!(
            hi == U256::ZERO,
            "exceeded modulus in safe_uint class"
        );
        Self::new_unsafe(
            self.value.clone() * other.value.clone(),
            lo,
        )
    }

    // ════════════════════════════════════════════════════════════════════
    //  Arithmetic: division
    // ════════════════════════════════════════════════════════════════════

    /// Division with pre-determined bounds on quotient and remainder bit sizes.
    pub fn divide(
        &self,
        other: &SafeUintT<P>,
        quotient_bit_size: usize,
        remainder_bit_size: usize,
        description: &str,
    ) -> SafeUintT<P> {
        self.divide_with_fn(
            other,
            quotient_bit_size,
            remainder_bit_size,
            description,
            |val, divisor| {
                let divisor_nz = divisor.to_nz().expect("division by zero");
                val.div_rem(&divisor_nz)
            },
        )
    }

    /// Division with custom quotient/remainder function.
    pub fn divide_with_fn(
        &self,
        other: &SafeUintT<P>,
        quotient_bit_size: usize,
        remainder_bit_size: usize,
        description: &str,
        get_quotient: impl Fn(U256, U256) -> (U256, U256),
    ) -> SafeUintT<P> {
        assert!(
            !self.value.is_constant(),
            "divide requires witness dividend"
        );
        assert!(
            quotient_bit_size <= Self::max_bit_num(),
            "quotient bit size exceeds max"
        );
        assert!(
            remainder_bit_size <= Self::max_bit_num(),
            "remainder bit size exceeds max"
        );

        let val = field_to_u256(&self.value.get_value());
        let divisor_val = field_to_u256(&other.value.get_value());
        let (quotient_val, remainder_val) = get_quotient(val, divisor_val);

        let ctx = self.value.get_context().as_ref().unwrap().clone();
        let quotient_field = FieldT::from_witness(ctx.clone(), u256_to_field(&quotient_val));
        let remainder_field = FieldT::from_witness(ctx, u256_to_field(&remainder_val));

        let quotient = SafeUintT::from_field_with_range(
            quotient_field,
            quotient_bit_size,
            &format!("divide method quotient: {}", description),
        );
        let remainder = SafeUintT::from_field_with_range(
            remainder_field,
            remainder_bit_size,
            &format!("divide method remainder: {}", description),
        );

        // quotient * other + remainder = self (implicitly checks no overflow)
        let int_val = quotient.clone() * other.clone() + remainder.clone();

        // Constrain remainder < divisor via: other - remainder - 1 >= 0
        let remainder_plus_one = remainder + SafeUintT::from_constant_u32(1);
        let _ = other.sub_impl(&remainder_plus_one);

        self.assert_equal(&int_val, "divide method quotient and/or remainder incorrect");

        quotient
    }

    fn div_impl(&self, other: &SafeUintT<P>) -> SafeUintT<P> {
        assert!(
            !self.value.is_constant(),
            "divide requires witness dividend"
        );

        let val = field_to_u256(&self.value.get_value());
        let divisor_val = field_to_u256(&other.value.get_value());
        let divisor_nz = divisor_val.to_nz().expect("division by zero");
        let (quotient_val, remainder_val) = val.div_rem(&divisor_nz);

        let ctx = self.value.get_context().as_ref().unwrap().clone();
        let quotient_field = FieldT::from_witness(ctx.clone(), u256_to_field(&quotient_val));
        let remainder_field = FieldT::from_witness(ctx, u256_to_field(&remainder_val));

        let q_bit_size = (self.current_max.get_msb() + 1) as usize;
        let r_bit_size = (other.current_max.get_msb() + 1) as usize;

        let quotient = SafeUintT::from_field_with_range(
            quotient_field,
            q_bit_size,
            "/ operator quotient",
        );
        let remainder = SafeUintT::from_field_with_range(
            remainder_field,
            r_bit_size,
            "/ operator remainder",
        );

        // quotient * other + remainder = self (implicitly checks no overflow)
        let int_val = quotient.clone() * other.clone() + remainder.clone();

        // Constrain remainder < divisor via: other - remainder - 1 >= 0
        let remainder_plus_one = remainder + SafeUintT::from_constant_u32(1);
        let _ = other.sub_impl(&remainder_plus_one);

        self.assert_equal(&int_val, "/ operator quotient and/or remainder incorrect");

        quotient
    }

    // ════════════════════════════════════════════════════════════════════
    //  Comparison and assertions
    // ════════════════════════════════════════════════════════════════════

    /// Circuit equality check: returns a `BoolT` constrained to represent
    /// whether `self == other`.
    pub fn eq(&self, other: &SafeUintT<P>) -> BoolT<P> {
        let diff = self.value.clone() - other.value.clone();
        field_is_zero(&diff)
    }

    /// Circuit inequality check.
    pub fn neq(&self, other: &SafeUintT<P>) -> BoolT<P> {
        self.eq(other).negate()
    }

    /// Check if the value is zero (returns a constrained `BoolT`).
    pub fn is_zero(&self) -> BoolT<P> {
        field_is_zero(&self.value)
    }

    /// Assert that `self == rhs`.
    pub fn assert_equal(&self, rhs: &SafeUintT<P>, msg: &str) {
        self.value.assert_equal(&rhs.value, msg);
    }

    /// Assert that `self` is zero.
    pub fn assert_is_zero(&self, msg: &str) {
        self.value.assert_is_zero(msg);
    }

    /// Assert that `self` is not zero.
    pub fn assert_is_not_zero(&self, msg: &str) {
        self.value.assert_is_not_zero(msg);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Bit operations
    // ════════════════════════════════════════════════════════════════════

    /// Decompose the value into three slices: `[lo, slice, hi]` where:
    /// - `lo` = bits `[0, lsb)`
    /// - `slice` = bits `[lsb, msb]`
    /// - `hi` = bits `[msb+1, ...)`
    ///
    /// Constrains: `lo + (slice << lsb) + (hi << (msb+1)) == self`.
    pub fn slice(&self, msb: u8, lsb: u8) -> [SafeUintT<P>; 3] {
        assert!(msb >= lsb, "msb must be >= lsb");
        let max_bit = Self::max_bit_num();
        assert!((msb as usize) < max_bit, "msb exceeds max bit length");

        let value_u256 = field_to_u256(&self.get_value());
        assert!(
            value_u256 < U256::ONE.wrapping_shl_vartime(max_bit as u32),
            "value exceeds max no-wrap bit length"
        );

        let msb_plus_one = (msb as u32) + 1;

        let lo_mask = U256::ONE.wrapping_shl_vartime(lsb as u32).wrapping_sub(&U256::ONE);
        let lo = value_u256.bitand(&lo_mask);

        let slice_bits = msb_plus_one - (lsb as u32);
        let slice_mask = U256::ONE.wrapping_shl_vartime(slice_bits).wrapping_sub(&U256::ONE);
        let slice_val = value_u256.wrapping_shr_vartime(lsb as u32).bitand(&slice_mask);

        let hi_bits = max_bit as u32 - msb_plus_one;
        let hi = value_u256.wrapping_shr_vartime(msb_plus_one);

        let (lo_wit, slice_wit, hi_wit);
        if self.value.is_constant() {
            lo_wit = SafeUintT::from_constant_u256(lo);
            slice_wit = SafeUintT::from_constant_u256(slice_val);
            hi_wit = SafeUintT::from_constant_u256(hi);
        } else {
            let ctx = self.value.get_context().as_ref().unwrap().clone();
            lo_wit = SafeUintT::from_field_with_range(
                FieldT::from_witness(ctx.clone(), u256_to_field(&lo)),
                lsb as usize,
                "slice lo_wit",
            );
            slice_wit = SafeUintT::from_field_with_range(
                FieldT::from_witness(ctx.clone(), u256_to_field(&slice_val)),
                slice_bits as usize,
                "slice slice_wit",
            );
            hi_wit = SafeUintT::from_field_with_range(
                FieldT::from_witness(ctx, u256_to_field(&hi)),
                hi_bits as usize,
                "slice hi_wit",
            );
        }

        // Constrain: hi * (1 << msb_plus_one) + lo + slice * (1 << lsb) == self
        let shift_hi = SafeUintT::from_constant_u256(
            U256::ONE.wrapping_shl_vartime(msb_plus_one),
        );
        let shift_lsb = SafeUintT::from_constant_u256(
            U256::ONE.wrapping_shl_vartime(lsb as u32),
        );
        let reconstructed = hi_wit.clone() * shift_hi + lo_wit.clone() + slice_wit.clone() * shift_lsb;
        self.assert_equal(&reconstructed, "slice decomposition constraint");

        [lo_wit, slice_wit, hi_wit]
    }

    // ════════════════════════════════════════════════════════════════════
    //  Conditional assignment
    // ════════════════════════════════════════════════════════════════════

    /// If `predicate` is true, return `lhs`; otherwise return `rhs`.
    pub fn conditional_assign(
        predicate: &BoolT<P>,
        lhs: &SafeUintT<P>,
        rhs: &SafeUintT<P>,
    ) -> SafeUintT<P> {
        let predicate_field = bool_to_field(predicate);
        let diff = lhs.value.clone() - rhs.value.clone();
        let new_val = diff.madd(&predicate_field, &rhs.value);
        let new_max = if lhs.current_max > rhs.current_max {
            lhs.current_max
        } else {
            rhs.current_max
        };
        Self::new_unsafe(new_val, new_max)
    }

    // ════════════════════════════════════════════════════════════════════
    //  Conversions
    // ════════════════════════════════════════════════════════════════════

    /// Convert to a `BoolT` (explicit cast).
    pub fn to_bool(&self) -> BoolT<P> {
        if self.value.is_constant() {
            BoolT::from_constant(self.value.get_value() == Field::one())
        } else {
            let normalized = self.value.normalize();
            BoolT::from_witness_index_unsafe(
                normalized.context.as_ref().unwrap().clone(),
                normalized.witness_index,
            )
        }
    }

    /// Convert to a `FieldT`.
    pub fn to_field(&self) -> FieldT<P> {
        self.value.clone()
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Operator implementations
// ════════════════════════════════════════════════════════════════════════

impl<P: FieldParams> Add for SafeUintT<P> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.add_impl(&rhs)
    }
}

impl<P: FieldParams> Add<&SafeUintT<P>> for &SafeUintT<P> {
    type Output = SafeUintT<P>;
    fn add(self, rhs: &SafeUintT<P>) -> SafeUintT<P> {
        self.add_impl(rhs)
    }
}

impl<P: FieldParams> Mul for SafeUintT<P> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.mul_impl(&rhs)
    }
}

impl<P: FieldParams> Mul<&SafeUintT<P>> for &SafeUintT<P> {
    type Output = SafeUintT<P>;
    fn mul(self, rhs: &SafeUintT<P>) -> SafeUintT<P> {
        self.mul_impl(rhs)
    }
}

impl<P: FieldParams> std::ops::Sub for SafeUintT<P> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self.sub_impl(&rhs)
    }
}

impl<P: FieldParams> std::ops::Sub<&SafeUintT<P>> for &SafeUintT<P> {
    type Output = SafeUintT<P>;
    fn sub(self, rhs: &SafeUintT<P>) -> SafeUintT<P> {
        self.sub_impl(rhs)
    }
}

impl<P: FieldParams> std::ops::Div for SafeUintT<P> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self.div_impl(&rhs)
    }
}

impl<P: FieldParams> std::ops::Div<&SafeUintT<P>> for &SafeUintT<P> {
    type Output = SafeUintT<P>;
    fn div(self, rhs: &SafeUintT<P>) -> SafeUintT<P> {
        self.div_impl(rhs)
    }
}

impl<P: FieldParams> AddAssign for SafeUintT<P> {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.add_impl(&rhs);
    }
}

impl<P: FieldParams> MulAssign for SafeUintT<P> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.mul_impl(&rhs);
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use std::rc::Rc;
    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use bbrs_ecc::curves::bn254::Bn254FrParams;

    type Fr = Field<Bn254FrParams>;
    type SuintCt = SafeUintT<Bn254FrParams>;
    type FrField = FieldT<Bn254FrParams>;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    fn check_circuit(builder: &BuilderRef<Bn254FrParams>) -> Result<(), String> {
        UltraCircuitChecker::check(&mut builder.borrow_mut())
    }

    // ── Test 1: Constructor with value out of range ──────────────────

    #[test]
    fn test_constructor_value_out_of_range_fails() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(100u64));
        // 100 does not fit in 2 bits (max 3)
        let _b = SuintCt::from_field_with_range(a, 2, "b");
        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 2: Constructor with value in range ──────────────────────

    #[test]
    fn test_constructor_value_in_range() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(100u64));
        let _b = SuintCt::from_field_with_range(a, 7, "b");
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 3: Multiply ─────────────────────────────────────────────

    #[test]
    fn test_multiply() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(9u64));
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, 4, "d");
        let e = c * d;
        assert_eq!(field_to_u256(&e.get_value()), U256::from_limbs([18, 0, 0, 0]));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 4: Multiply overflow detection ──────────────────────────

    #[test]
    #[should_panic(expected = "exceeded modulus in safe_uint class")]
    fn test_multiply_overflow() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let mut c = SuintCt::from_field_with_range(a.clone(), 2, "c");
        let d = SuintCt::from_field_with_range(a, 2, "d");
        // 3^160 < modulus < 3^161, so 159 iterations should succeed
        for _ in 0..159 {
            c = c * d.clone();
        }
        // 160th should overflow
        let _ = c * d;
    }

    // ── Test 5: Multiply constant overflow detection ─────────────────

    #[test]
    #[should_panic(expected = "exceeded modulus in safe_uint class")]
    fn test_multiply_constant_overflow() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let mut c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_constant_field(Fr::from(2u64));
        // Constants have tighter max (2 instead of 3), so more iterations before overflow
        // 2^254 > modulus, so 252 iterations succeed but 253rd fails
        for _ in 0..252 {
            c = c * d.clone();
        }
        let _ = c * d;
    }

    // ── Test 6: Add ──────────────────────────────────────────────────

    #[test]
    fn test_add() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(9u64));
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, 4, "d");
        let e = c + d;
        assert_eq!(field_to_u256(&e.get_value()), U256::from_limbs([11, 0, 0, 0]));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 7: Add overflow detection ───────────────────────────────

    #[test]
    #[should_panic(expected = "exceeded modulus in safe_uint class")]
    fn test_add_overflow() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let mut c = SuintCt::from_field_with_range(a.clone(), 2, "c");
        let d = SuintCt::from_field_with_range(a, 2, "d");
        // Build c up near the modulus via multiplication
        for _ in 0..159 {
            c = c * d.clone();
        }
        // c + c should overflow
        let _ = c.clone() + c;
    }

    // ── Test 8: Subtract method ──────────────────────────────────────

    #[test]
    fn test_subtract() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(9u64));
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, 4, "d");
        // 9 - 2 = 7, fits in 3 bits
        let e = d.subtract(&c, 3, "d - c");
        assert_eq!(field_to_u256(&e.get_value()), U256::from_limbs([7, 0, 0, 0]));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 9: Subtract result out of range ─────────────────────────

    #[test]
    fn test_subtract_result_out_of_range() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(9u64));
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, 4, "d");
        // 9 - 2 = 7, does NOT fit in 2 bits
        let _e = d.subtract(&c, 2, "d - c");
        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 10: Subtract underflow general case ─────────────────────

    #[test]
    fn test_subtract_underflow_general() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(0u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(1u64));
        let c = SuintCt::from_field_with_range(a, 0, "c");
        let d = SuintCt::from_field_with_range(b, 1, "d");
        // 0 - 1 underflows, range constraint should catch it
        let _e = c.subtract(&d, SuintCt::max_bit_num(), "");
        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 11: Subtract underflow special case ─────────────────────

    #[test]
    #[should_panic(expected = "maximum value exceeded in safe_uint subtract")]
    fn test_subtract_underflow_special() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let half_mod = {
            let m = U256::from_limbs(Bn254FrParams::MODULUS);
            let two = U256::from_limbs([2, 0, 0, 0]);
            let nz = two.to_nz().unwrap();
            let (q, _) = m.div_rem(&nz);
            u256_to_field::<Bn254FrParams>(&q)
        };
        let b = FrField::from_witness(builder.clone(), half_mod);
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, SuintCt::max_bit_num(), "d");
        let _ = c.subtract(&d, SuintCt::max_bit_num(), "");
    }

    // ── Test 12: Minus operator ──────────────────────────────────────

    #[test]
    fn test_minus_operator() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(9u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let c = SuintCt::from_field_with_range(a, 4, "c");
        let d = SuintCt::from_field_with_range(b, 2, "d");
        let e = c - d;
        assert_eq!(field_to_u256(&e.get_value()), U256::from_limbs([7, 0, 0, 0]));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 13: Minus operator valid on zero ────────────────────────

    #[test]
    fn test_minus_operator_valid_on_zero() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, 3, "d");
        // 2 - 2 = 0, should not underflow
        let e = c - d;
        assert_eq!(field_to_u256(&e.get_value()), U256::ZERO);
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 14: Minus underflow general case 1 ──────────────────────

    #[test]
    fn test_minus_underflow_general1() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let half_mod = {
            let m = U256::from_limbs(Bn254FrParams::MODULUS);
            let two = U256::from_limbs([2, 0, 0, 0]);
            let nz = two.to_nz().unwrap();
            let (q, _) = m.div_rem(&nz);
            u256_to_field::<Bn254FrParams>(&q)
        };
        let b = FrField::from_witness(builder.clone(), half_mod);
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, SuintCt::max_bit_num(), "d");
        // c - d underflows, range constraint catches it (c.current_max = 3, small bit range)
        let _e = c - d;
        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 15: Minus underflow general case 2 ──────────────────────

    #[test]
    fn test_minus_underflow_general2() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(3u64));
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, 3, "d");
        // 2 - 3 underflows
        let _e = c - d;
        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 16: Minus underflow special case 1 ──────────────────────

    #[test]
    #[should_panic(expected = "maximum value exceeded in safe_uint minus operator")]
    fn test_minus_underflow_special1() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(1u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(0u64));
        let c = SuintCt::from_field_with_range(a, SuintCt::max_bit_num(), "c");
        let d = SuintCt::from_field_with_range(b, SuintCt::max_bit_num(), "d");
        // Even though 1 - 0 is not underflow, max ranges overlap too much
        let _ = c - d;
    }

    // ── Test 17: Minus underflow special case 2 ──────────────────────

    #[test]
    #[should_panic(expected = "maximum value exceeded in safe_uint minus operator")]
    fn test_minus_underflow_special2() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(0u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(1u64));
        let c = SuintCt::from_field_with_range(a, SuintCt::max_bit_num(), "c");
        let d = SuintCt::from_field_with_range(b, SuintCt::max_bit_num(), "d");
        let _ = c - d;
    }

    // ── Test 18: Divide method ───────────────────────────────────────

    #[test]
    fn test_divide_method() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(9u64));
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, 4, "d");
        // 9 / 2 = 4 remainder 1; quotient fits in 3 bits, remainder in 1 bit
        let q = d.divide(&c, 3, 1, "d/c");
        assert_eq!(field_to_u256(&q.get_value()), U256::from_limbs([4, 0, 0, 0]));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 19: Divide method quotient range too small ──────────────

    #[test]
    fn test_divide_quotient_range_too_small() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(32u64));
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, 6, "d");
        // 32 / 2 = 16, needs 5 bits but we only allow 4
        let _q = d.divide(&c, 4, 1, "d/c");
        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 20: Divide remainder too large ──────────────────────────

    #[test]
    #[should_panic]
    fn test_divide_remainder_too_large() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let c = SuintCt::from_field_with_range(a, 3, "c");
        let modulus_minus_one_div_3 = {
            let m = SuintCt::max_value();
            let three = U256::from_limbs([3, 0, 0, 0]);
            let nz = three.to_nz().unwrap();
            let (q, _) = m.div_rem(&nz);
            q
        };
        let d = SuintCt::from_constant_u256(modulus_minus_one_div_3);
        let _ = c / d;
    }

    // ── Test 21: Divide quotient/remainder incorrect ─────────────────

    #[test]
    fn test_divide_quotient_remainder_incorrect() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(19u64));
        let c = SuintCt::from_field_with_range(a, 3, "c");
        let d = SuintCt::from_field_with_range(b, 5, "d");
        // Supply wrong quotient/remainder: 2, 3 (correct is 3, 4)
        let _q = d.divide_with_fn(&c, 3, 2, "d/c", |_val, _div| {
            (U256::from_limbs([2, 0, 0, 0]), U256::from_limbs([3, 0, 0, 0]))
        });
        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 22: Divide quotient/remainder mod r ─────────────────────

    #[test]
    fn test_divide_quotient_remainder_mod_r() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(5u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(19u64));
        let c = SuintCt::from_field_with_range(a, 3, "c");
        let d = SuintCt::from_field_with_range(b, 5, "d");
        // Supply quotient = (19/5 mod r), remainder = 0 — correct mod r but not integer division
        let _q = d.divide_with_fn(&c, 3, 1, "d/c", |val, div| {
            let fr_quot = u256_to_field::<Bn254FrParams>(&val)
                * u256_to_field::<Bn254FrParams>(&div).invert();
            (field_to_u256(&fr_quot), U256::ZERO)
        });
        // The field quotient is huge, so it should fail the 3-bit range check
        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 23: Div operator ────────────────────────────────────────

    #[test]
    fn test_div_operator() {
        let builder = make_builder();
        let a_field = FrField::from_witness(builder.clone(), Fr::from(1000u64));
        let mut a = SuintCt::from_field_with_range(a_field, 10, "a");
        let b = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), Fr::from(2u64)),
            2,
            "b",
        );
        a = a / b;
        assert_eq!(field_to_u256(&a.get_value()), U256::from_limbs([500, 0, 0, 0]));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 24: Divide operator comprehensive ───────────────────────

    #[test]
    fn test_divide_operator_success() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(9u64));
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, 4, "d");
        let _q = d / c;
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_divide_operator_quotient_too_small() {
        let builder = make_builder();
        let a = FrField::from_witness(builder.clone(), Fr::from(2u64));
        let b = FrField::from_witness(builder.clone(), Fr::from(32u64));
        let c = SuintCt::from_field_with_range(a, 2, "c");
        let d = SuintCt::from_field_with_range(b, 5, "d");
        // 32 / 2 = 16, quotient needs more bits than d allows
        let _q = d / c;
        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 25: Slice ───────────────────────────────────────────────

    #[test]
    fn test_slice() {
        let builder = make_builder();
        // 0b11110110101001011 = 126283
        let a = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), Fr::from(126283u64)),
            17,
            "a",
        );
        let slices = a.slice(10, 3);
        // lo = bits[0,3) = 0b011 = 3
        assert_eq!(field_to_u256(&slices[0].get_value()), U256::from_limbs([3, 0, 0, 0]));
        // slice = bits[3,10] = 0b10101001 = 169
        assert_eq!(field_to_u256(&slices[1].get_value()), U256::from_limbs([169, 0, 0, 0]));
        // hi = bits[11,...) = 0b111101 = 61
        assert_eq!(field_to_u256(&slices[2].get_value()), U256::from_limbs([61, 0, 0, 0]));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 26: Slice equal msb/lsb ─────────────────────────────────

    #[test]
    fn test_slice_equal_msb_lsb() {
        let builder = make_builder();
        // 0b11110110101001011 = 126283
        let a = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), Fr::from(126283u64)),
            17,
            "a",
        );
        let slices = a.slice(6, 6);
        // lo = bits[0,6) = 0b001011 = 11
        assert_eq!(field_to_u256(&slices[0].get_value()), U256::from_limbs([11, 0, 0, 0]));
        // slice = bits[6,6] = bit 6 = 1
        assert_eq!(field_to_u256(&slices[1].get_value()), U256::from_limbs([1, 0, 0, 0]));
        // hi = bits[7,...) = 0b1111011010 = 986
        assert_eq!(field_to_u256(&slices[2].get_value()), U256::from_limbs([986, 0, 0, 0]));
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 27: Slice random ────────────────────────────────────────

    #[test]
    fn test_slice_random() {
        let builder = make_builder();
        let lsb: u8 = 106;
        let msb: u8 = 189;
        // Generate a random value that fits in 252 bits
        let a_fr = Fr::random_element();
        let a_u256 = field_to_u256(&a_fr);
        let mask_252 = U256::ONE.wrapping_shl_vartime(252).wrapping_sub(&U256::ONE);
        let a_masked = a_u256.bitand(&mask_252);
        let a_field_val = u256_to_field::<Bn254FrParams>(&a_masked);
        let a = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), a_field_val),
            252,
            "a",
        );
        let slices = a.slice(msb, lsb);

        let lo_mask = U256::ONE.wrapping_shl_vartime(lsb as u32).wrapping_sub(&U256::ONE);
        let expected0 = a_masked.bitand(&lo_mask);
        let slice_mask = U256::ONE
            .wrapping_shl_vartime((msb - lsb) as u32 + 1)
            .wrapping_sub(&U256::ONE);
        let expected1 = a_masked.wrapping_shr_vartime(lsb as u32).bitand(&slice_mask);
        let expected2 = a_masked.wrapping_shr_vartime((msb as u32) + 1);

        assert_eq!(field_to_u256(&slices[0].get_value()), expected0);
        assert_eq!(field_to_u256(&slices[1].get_value()), expected1);
        assert_eq!(field_to_u256(&slices[2].get_value()), expected2);
        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test 28: Div remainder constraint (operator) ─────────────────

    #[test]
    fn test_operator_div_remainder_constraint() {
        let builder = make_builder();
        let val = U256::from_limbs([5, 0, 0, 0]);
        let val_field = u256_to_field::<Bn254FrParams>(&val);

        let a = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), val_field),
            32,
            "a",
        );
        let b = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), val_field),
            32,
            "b",
        );

        // Manually supply wrong witnesses: quotient = 0, remainder = val
        // This should fail because remainder must be < divisor
        let quotient_val = U256::ZERO;
        let remainder_val = val;

        let quotient_field = FrField::from_witness(builder.clone(), u256_to_field(&quotient_val));
        let remainder_field = FrField::from_witness(builder.clone(), u256_to_field(&remainder_val));
        let q_bits = (a.current_max.get_msb() + 1) as usize;
        let quotient = SuintCt::from_field_with_range(quotient_field, q_bits, "q");
        let remainder = SuintCt::from_field_with_range(remainder_field, q_bits, "r");

        let int_val = quotient * b.clone() + remainder.clone();

        // Constrain remainder < divisor
        let delta = b - remainder - SuintCt::from_constant_u32(1);
        let delta_field = FieldT::from_witness_index(
            delta.value.get_context().as_ref().unwrap().clone(),
            delta.value.get_witness_index(),
        );
        delta_field.create_range_constraint(
            (a.current_max.get_msb() + 1) as usize,
            "delta range",
        );

        a.assert_equal(&int_val, "div constraint");

        // Should fail: remainder = val = divisor, not < divisor
        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 29: Div remainder constraint (divide method) ────────────

    #[test]
    fn test_div_remainder_constraint_method() {
        let builder = make_builder();
        let val = U256::from_limbs([5, 0, 0, 0]);
        let val_field = u256_to_field::<Bn254FrParams>(&val);

        let a = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), val_field),
            32,
            "a",
        );
        let b = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), val_field),
            32,
            "b",
        );

        // Supply bad witnesses: quotient = 0, remainder = val
        let _q = a.divide_with_fn(&b, 32, 32, "", |val, _div| {
            (U256::ZERO, val)
        });

        assert!(check_circuit(&builder).is_err());
    }

    // ── Test 30: Is zero ─────────────────────────────────────────────

    #[test]
    fn test_is_zero() {
        let builder = make_builder();
        let a = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), Fr::from(0u64)),
            8,
            "a",
        );
        let b = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), Fr::from(5u64)),
            8,
            "b",
        );
        assert!(a.is_zero().get_value());
        assert!(!b.is_zero().get_value());
        assert!(check_circuit(&builder).is_ok());

        // Constant cases
        let c = SuintCt::from_constant_u32(0);
        assert!(c.is_zero().get_value());
        let d = SuintCt::from_constant_u32(42);
        assert!(!d.is_zero().get_value());
    }

    // ── Test 31: Assert equal ────────────────────────────────────────

    #[test]
    fn test_assert_equal() {
        let builder = make_builder();
        let a = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), Fr::from(42u64)),
            8,
            "a",
        );
        let b = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), Fr::from(42u64)),
            8,
            "b",
        );
        a.assert_equal(&b, "should be equal");
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_equal_fails() {
        let builder = make_builder();
        let a = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), Fr::from(42u64)),
            8,
            "a",
        );
        let b = SuintCt::from_field_with_range(
            FrField::from_witness(builder.clone(), Fr::from(43u64)),
            8,
            "b",
        );
        a.assert_equal(&b, "should fail");
        assert!(builder.borrow().base.failed());
    }
}
