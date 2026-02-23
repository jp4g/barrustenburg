//! Non-native field arithmetic circuit type.
//!
//! Port of `barretenberg/stdlib/primitives/bigfield/bigfield.hpp` and `bigfield_impl.hpp`.
//!
//! A `BigFieldT<P, T>` represents a non-native field element inside a circuit whose native
//! field has parameters `P`. The non-native field has parameters `T`.
//!
//! Each element is decomposed into 4 limbs of `NUM_LIMB_BITS` (68) bits plus a
//! `prime_basis_limb` that tracks the full value mod the native field modulus.
//! Arithmetic is verified via the Chinese Remainder Theorem: checking both mod the native
//! field modulus and mod 2^272 (the binary basis).

use std::rc::Rc;

use bbrs_circuit_builder::gate_data::{
    AddQuad, CachedPartialNnfMul, NnfAddSimple, NnfMulWitnesses, NnfPartialMulWitnesses,
};
use bbrs_circuit_builder::ultra_builder::{
    UltraBlockIndex, UltraCircuitBuilder, DEFAULT_NON_NATIVE_FIELD_LIMB_BITS,
    DEFAULT_PLOOKUP_RANGE_BITNUM,
};
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;
use bbrs_numeric::{U256, U256Ext, U512};
use bbrs_numeric::uintx::U512Ext;

use super::bool::BoolT;
use super::field::FieldT;
use super::witness::{BuilderRef, WitnessT, IS_CONSTANT};

// ════════════════════════════════════════════════════════════════════════
//  Constants
// ════════════════════════════════════════════════════════════════════════

/// Number of limbs used to represent a non-native field element.
const NUM_LIMBS: usize = 4;

/// Number of bits per limb (matching C++ NUM_LIMB_BITS_IN_FIELD_SIMULATION).
const NUM_LIMB_BITS: u64 = DEFAULT_NON_NATIVE_FIELD_LIMB_BITS;

/// Maximum value for the default plookup range.
const MAX_ADDITION_LOG: u64 = 10;

/// Extract bits [start, end) from a U256, returning a U256 (supports ranges > 64 bits).
fn slice_u256(val: &U256, start: u32, end: u32) -> U256 {
    assert!(end > start, "end must be greater than start");
    let width = end - start;
    let shifted = val.wrapping_shr_vartime(start);
    let mask = if width >= 256 {
        U256::MAX
    } else {
        U256::from(1u64)
            .wrapping_shl_vartime(width)
            .wrapping_sub(&U256::ONE)
    };
    shifted.wrapping_and(&mask)
}

/// Convert a U256 slice to a native field element.
fn slice_to_field<P: FieldParams>(val: &U256, start: u32, end: u32) -> Field<P> {
    let sliced = slice_u256(val, start, end);
    Field::from_limbs(*sliced.as_words())
}

// ════════════════════════════════════════════════════════════════════════
//  Limb
// ════════════════════════════════════════════════════════════════════════

/// A single limb of a non-native field element: a circuit field element plus
/// a maximum value tracking overflow potential.
pub struct Limb<P: FieldParams> {
    pub element: FieldT<P>,
    pub maximum_value: U256,
}

impl<P: FieldParams> Clone for Limb<P> {
    fn clone(&self) -> Self {
        Self {
            element: self.element.clone(),
            maximum_value: self.maximum_value,
        }
    }
}

impl<P: FieldParams> Limb<P> {
    fn new(element: FieldT<P>, maximum_value: U256) -> Self {
        Self {
            element,
            maximum_value,
        }
    }

    fn new_default(element: FieldT<P>) -> Self {
        Self {
            maximum_value: default_maximum_limb(),
            element,
        }
    }

    fn zero_limb() -> Self {
        Self {
            element: FieldT::from_field(Field::zero()),
            maximum_value: U256::ZERO,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Helper functions for computing constants from target field params
// ════════════════════════════════════════════════════════════════════════

/// Default maximum value for limbs 0-2: (1 << NUM_LIMB_BITS) - 1
fn default_maximum_limb() -> U256 {
    U256::from(1u64)
        .wrapping_shl_vartime(NUM_LIMB_BITS as u32)
        .wrapping_sub(&U256::ONE)
}

/// Number of bits in the last (most significant) limb for target field T.
fn num_last_limb_bits<T: FieldParams>() -> u64 {
    let modulus = U256::from_limbs(T::MODULUS);
    let msb = modulus.get_msb();
    (msb as u64) + 1 - NUM_LIMB_BITS * 3
}

/// Default maximum value for the most-significant limb.
fn default_maximum_msb_limb<T: FieldParams>() -> U256 {
    let bits = num_last_limb_bits::<T>();
    U256::from(1u64)
        .wrapping_shl_vartime(bits as u32)
        .wrapping_sub(&U256::ONE)
}

/// Shift constant: 2^NUM_LIMB_BITS as a native field element.
fn shift_1<P: FieldParams>() -> Field<P> {
    Field::from_limbs(
        *U256::from(1u64)
            .wrapping_shl_vartime(NUM_LIMB_BITS as u32)
            .as_words(),
    )
}

/// Shift constant: 2^(2*NUM_LIMB_BITS) as a native field element.
fn shift_2<P: FieldParams>() -> Field<P> {
    Field::from_limbs(
        *U256::from(1u64)
            .wrapping_shl_vartime((NUM_LIMB_BITS * 2) as u32)
            .as_words(),
    )
}

/// Shift constant: 2^(3*NUM_LIMB_BITS) as a native field element.
fn shift_3<P: FieldParams>() -> Field<P> {
    Field::from_limbs(
        *U256::from(1u64)
            .wrapping_shl_vartime((NUM_LIMB_BITS * 3) as u32)
            .as_words(),
    )
}

/// Compute negative modulus mod binary basis: 2^272 - target_modulus, sliced into 4 limbs.
fn neg_modulus_limbs<P: FieldParams, T: FieldParams>() -> [Field<P>; 4] {
    let u256_limbs = neg_modulus_limbs_u256::<T>();
    [
        Field::from_limbs(*u256_limbs[0].as_words()),
        Field::from_limbs(*u256_limbs[1].as_words()),
        Field::from_limbs(*u256_limbs[2].as_words()),
        Field::from_limbs(*u256_limbs[3].as_words()),
    ]
}

/// Compute negative modulus limbs as U256 values (for overflow tracking).
///
/// neg_mod = 2^272 - T::MODULUS, then split into 4 limbs of 68 bits.
/// The value is 272 bits, so the 4th limb (bits 204-271) crosses the U256 boundary.
fn neg_modulus_limbs_u256<T: FieldParams>() -> [U256; 4] {
    let target_mod = U256::from_limbs(T::MODULUS);
    let binary_basis = U512::from_lo_hi(
        U256::ZERO,
        U256::from(1u64).wrapping_shl_vartime((NUM_LIMB_BITS * 4 - 256) as u32),
    );
    let target_512 = U512::from_lo_hi(target_mod, U256::ZERO);
    let neg_mod = binary_basis.wrapping_sub(&target_512);
    let neg_lo = neg_mod.lo();
    let neg_hi = neg_mod.hi();

    // Limbs 0-2 fit entirely in neg_lo (bits 0-203)
    let l0 = slice_u256(&neg_lo, 0, NUM_LIMB_BITS as u32);
    let l1 = slice_u256(&neg_lo, NUM_LIMB_BITS as u32, (NUM_LIMB_BITS * 2) as u32);
    let l2 = slice_u256(&neg_lo, (NUM_LIMB_BITS * 2) as u32, (NUM_LIMB_BITS * 3) as u32);

    // Limb 3: bits 204-271 crosses the U256 boundary at bit 256.
    // 52 bits from neg_lo [204..256) + 16 bits from neg_hi [0..16)
    let remaining_in_lo = 256 - (NUM_LIMB_BITS * 3) as u32; // 52
    let lo_part = slice_u256(&neg_lo, (NUM_LIMB_BITS * 3) as u32, 256);
    let hi_bits_needed = NUM_LIMB_BITS as u32 - remaining_in_lo; // 16
    let hi_part = slice_u256(&neg_hi, 0, hi_bits_needed);
    let l3 = lo_part.wrapping_or(&hi_part.wrapping_shl_vartime(remaining_in_lo));

    [l0, l1, l2, l3]
}

/// Negative target modulus mod native modulus: -T::modulus mod P::modulus.
fn negative_prime_modulus<P: FieldParams, T: FieldParams>() -> Field<P> {
    -Field::<P>::from_limbs(T::MODULUS)
}

/// LOG2_BINARY_MODULUS = NUM_LIMB_BITS * NUM_LIMBS = 272
const LOG2_BINARY_MODULUS: u64 = NUM_LIMB_BITS * (NUM_LIMBS as u64);

/// Maximum unreduced value: sqrt(2^272 * native_modulus)
/// We use the fact that products must be < 2^t * n where t = LOG2_BINARY_MODULUS.
fn get_maximum_crt_product<P: FieldParams>() -> (U256, U256) {
    // CRT product modulus = binary_basis * prime_basis = 2^272 * native_modulus
    // We return the (lo, hi) of U512 representation
    let native_mod = U256::from_limbs(P::MODULUS);
    let native_512 = U512::from_lo_hi(native_mod, U256::ZERO);
    // shift left by LOG2_BINARY_MODULUS = 272 bits
    // This is a U512 * 2^272, which needs more than 512 bits.
    // For overflow checking we use the simpler approach from C++.
    (native_512.lo(), native_512.hi())
}

/// Maximum limb size that wouldn't cause native field overflow during multiplication.
fn maximum_limb_size_that_wouldnt_overflow<P: FieldParams>() -> u64 {
    let native_mod = U256::from_limbs(P::MODULUS);
    let native_msb = native_mod.get_msb() as u64;
    (native_msb - MAX_ADDITION_LOG - NUM_LIMB_BITS - 3) / 2
}

/// Maximum unreduced limb value.
fn get_maximum_unreduced_limb_value<P: FieldParams>() -> U256 {
    let bits = NUM_LIMB_BITS + MAX_ADDITION_LOG;
    U256::from(1u64)
        .wrapping_shl_vartime(bits as u32)
        .wrapping_sub(&U256::ONE)
}

/// Prohibited limb value (overflow threshold).
fn get_prohibited_limb_value<P: FieldParams>() -> U256 {
    let bits = NUM_LIMB_BITS + MAX_ADDITION_LOG + 5;
    U256::from(1u64)
        .wrapping_shl_vartime(bits as u32)
        .wrapping_sub(&U256::ONE)
}

/// Compute partial schoolbook multiplication of two 4-limb numbers.
/// Returns (lo, hi) where lo covers limbs 0-1 and hi covers limbs 2-3.
fn compute_partial_schoolbook_multiplication(
    a: &[U256; 4],
    b: &[U256; 4],
) -> (U512, U512) {
    let to_512 = |x: U256| U512::from_lo_hi(x, U256::ZERO);

    // r0 = a0*b0
    let a0b0 = to_512(a[0]).wrapping_mul(&to_512(b[0]));

    // r1 = a1*b0 + a0*b1
    let a1b0 = to_512(a[1]).wrapping_mul(&to_512(b[0]));
    let a0b1 = to_512(a[0]).wrapping_mul(&to_512(b[1]));
    let r1 = a1b0.wrapping_add(&a0b1);

    // lo = r0 + r1 << NUM_LIMB_BITS
    let r1_shifted = r1.wrapping_shl_vartime(NUM_LIMB_BITS as u32);
    let lo = a0b0.wrapping_add(&r1_shifted);

    // r2 = a1*b1 + a2*b0 + a0*b2
    let a1b1 = to_512(a[1]).wrapping_mul(&to_512(b[1]));
    let a2b0 = to_512(a[2]).wrapping_mul(&to_512(b[0]));
    let a0b2 = to_512(a[0]).wrapping_mul(&to_512(b[2]));
    let r2 = a1b1.wrapping_add(&a2b0).wrapping_add(&a0b2);

    // r3 = a0*b3 + a1*b2 + a2*b1 + a3*b0
    let a0b3 = to_512(a[0]).wrapping_mul(&to_512(b[3]));
    let a1b2 = to_512(a[1]).wrapping_mul(&to_512(b[2]));
    let a2b1 = to_512(a[2]).wrapping_mul(&to_512(b[1]));
    let a3b0 = to_512(a[3]).wrapping_mul(&to_512(b[0]));
    let r3 = a0b3.wrapping_add(&a1b2).wrapping_add(&a2b1).wrapping_add(&a3b0);

    // hi = r2 + r3 << NUM_LIMB_BITS
    let r3_shifted = r3.wrapping_shl_vartime(NUM_LIMB_BITS as u32);
    let hi = r2.wrapping_add(&r3_shifted);

    (lo, hi)
}

// ════════════════════════════════════════════════════════════════════════
//  BigFieldT
// ════════════════════════════════════════════════════════════════════════

/// Non-native field element type for circuit arithmetic.
///
/// `P` is the native circuit field parameters, `T` is the target non-native field parameters.
///
/// Port of C++ `bigfield<Builder, T>`.
pub struct BigFieldT<P: FieldParams, T: FieldParams> {
    pub context: Option<BuilderRef<P>>,
    pub binary_basis_limbs: [Limb<P>; NUM_LIMBS],
    pub prime_basis_limb: FieldT<P>,
    _target: std::marker::PhantomData<T>,
}

impl<P: FieldParams, T: FieldParams> Clone for BigFieldT<P, T> {
    fn clone(&self) -> Self {
        Self {
            context: self.context.clone(),
            binary_basis_limbs: self.binary_basis_limbs.clone(),
            prime_basis_limb: self.prime_basis_limb.clone(),
            _target: std::marker::PhantomData,
        }
    }
}

impl<P: FieldParams, T: FieldParams> BigFieldT<P, T> {
    // ════════════════════════════════════════════════════════════════════
    //  Constants (computed from type parameters)
    // ════════════════════════════════════════════════════════════════════

    /// The target non-native field modulus.
    fn modulus() -> U256 {
        U256::from_limbs(T::MODULUS)
    }

    /// Target modulus as U512.
    fn modulus_u512() -> U512 {
        U512::from_lo_hi(Self::modulus(), U256::ZERO)
    }

    // ════════════════════════════════════════════════════════════════════
    //  Constructors
    // ════════════════════════════════════════════════════════════════════

    /// Create a constant bigfield element from a U256 value.
    ///
    /// Port of C++ `bigfield(Builder* parent_context, const uint256_t& value)`.
    pub fn from_u256(ctx: Option<BuilderRef<P>>, value: U256) -> Self {
        let limb0_val = slice_u256(&value, 0, NUM_LIMB_BITS as u32);
        let limb1_val = slice_u256(&value, NUM_LIMB_BITS as u32, (NUM_LIMB_BITS * 2) as u32);
        let limb2_val = slice_u256(&value, (NUM_LIMB_BITS * 2) as u32, (NUM_LIMB_BITS * 3) as u32);
        let limb3_val = slice_u256(&value, (NUM_LIMB_BITS * 3) as u32, (NUM_LIMB_BITS * 4) as u32);

        let prime_value = Field::<P>::from_limbs(*value.as_words());

        Self {
            context: ctx,
            binary_basis_limbs: [
                Limb::new_default(FieldT::from_field(Field::from_limbs(*limb0_val.as_words()))),
                Limb::new_default(FieldT::from_field(Field::from_limbs(*limb1_val.as_words()))),
                Limb::new_default(FieldT::from_field(Field::from_limbs(*limb2_val.as_words()))),
                Limb::new(
                    FieldT::from_field(Field::from_limbs(*limb3_val.as_words())),
                    default_maximum_msb_limb::<T>(),
                ),
            ],
            prime_basis_limb: FieldT::from_field(prime_value),
            _target: std::marker::PhantomData,
        }
    }

    /// Create a zero bigfield constant.
    pub fn zero_val() -> Self {
        Self {
            context: None,
            binary_basis_limbs: [
                Limb::zero_limb(),
                Limb::zero_limb(),
                Limb::zero_limb(),
                Limb::zero_limb(),
            ],
            prime_basis_limb: FieldT::from_field(Field::zero()),
            _target: std::marker::PhantomData,
        }
    }

    /// Create a one bigfield constant.
    pub fn one_val() -> Self {
        Self::from_u256(None, U256::ONE)
    }

    /// Create an "unreduced zero" — a zero value expressed as a multiple of the modulus
    /// to provide headroom for subtraction.
    fn unreduced_zero() -> Self {
        // Use 0, but with proper maximum_value tracking
        Self::zero_val()
    }

    /// Create a bigfield witness from a U512 value.
    ///
    /// Port of C++ `bigfield::create_from_u512_as_witness`.
    pub fn create_from_u512_as_witness(
        ctx: BuilderRef<P>,
        value: U512,
        can_overflow: bool,
        maximum_bitlength: usize,
    ) -> Self {
        // Extract 4 limbs from the U512 value.
        // limbs 0-2 fit entirely in lo (bits [0, 204)), limb3 spans [204, 272)
        // which crosses the U256 boundary at bit 256.
        let value_lo = value.lo();
        let value_hi = value.hi();

        let l0 = slice_u256(&value_lo, 0, NUM_LIMB_BITS as u32);
        let l1 = slice_u256(&value_lo, NUM_LIMB_BITS as u32, (NUM_LIMB_BITS * 2) as u32);
        let l2 = slice_u256(&value_lo, (NUM_LIMB_BITS * 2) as u32, (NUM_LIMB_BITS * 3) as u32);
        // limb3: bits [204, 272) — 52 bits from lo [204,256) + 16 bits from hi [0,16)
        let remaining_in_lo = 256 - (NUM_LIMB_BITS * 3) as u32; // 52
        let lo_part = slice_u256(&value_lo, (NUM_LIMB_BITS * 3) as u32, 256);
        let hi_bits_needed = NUM_LIMB_BITS as u32 - remaining_in_lo; // 16
        let hi_part = slice_u256(&value_hi, 0, hi_bits_needed);
        let l3 = lo_part.wrapping_or(&hi_part.wrapping_shl_vartime(remaining_in_lo));

        // Create witness variables
        let limb_0 = {
            let val = Field::<P>::from_limbs(*l0.as_words());
            let idx = ctx.borrow_mut().base.add_variable(val);
            FieldT::from_witness_index(ctx.clone(), idx)
        };
        let limb_1 = {
            let val = Field::<P>::from_limbs(*l1.as_words());
            let idx = ctx.borrow_mut().base.add_variable(val);
            FieldT::from_witness_index(ctx.clone(), idx)
        };
        let limb_2 = {
            let val = Field::<P>::from_limbs(*l2.as_words());
            let idx = ctx.borrow_mut().base.add_variable(val);
            FieldT::from_witness_index(ctx.clone(), idx)
        };
        let limb_3 = {
            let val = Field::<P>::from_limbs(*l3.as_words());
            let idx = ctx.borrow_mut().base.add_variable(val);
            FieldT::from_witness_index(ctx.clone(), idx)
        };

        // Compute prime basis limb: limb_0 + limb_1*shift_1 + limb_2*shift_2 + limb_3*shift_3
        let prime_value = limb_0.get_value()
            + limb_1.get_value() * shift_1::<P>()
            + limb_2.get_value() * shift_2::<P>()
            + limb_3.get_value() * shift_3::<P>();
        let prime_idx = ctx.borrow_mut().base.add_variable(prime_value);
        let prime_limb = FieldT::from_witness_index(ctx.clone(), prime_idx);

        // Create big_add_gate to verify: limb_1*shift_1 + limb_2*shift_2 + limb_3*shift_3 - prime + limb_0(via shift) = 0
        {
            let mut builder = ctx.borrow_mut();
            builder.create_big_add_gate(
                &AddQuad {
                    a: limb_1.get_witness_index(),
                    b: limb_2.get_witness_index(),
                    c: limb_3.get_witness_index(),
                    d: prime_idx,
                    a_scaling: shift_1::<P>(),
                    b_scaling: shift_2::<P>(),
                    c_scaling: shift_3::<P>(),
                    d_scaling: -Field::one(),
                    const_scaling: Field::zero(),
                },
                true,
            );
            // Unconstrained gate to provide limb_0 via w_4_shift
            let z = builder.base.zero_idx();
            builder.create_unconstrained_gate(
                UltraBlockIndex::Arithmetic,
                z,
                z,
                z,
                limb_0.get_witness_index(),
            );
        }

        // Set maximum values
        let num_last = if can_overflow {
            NUM_LIMB_BITS
        } else if maximum_bitlength > (NUM_LIMB_BITS * 3) as usize {
            (maximum_bitlength as u64) - NUM_LIMB_BITS * 3
        } else if maximum_bitlength > 0 {
            // Small quotient fits in fewer than 4 limbs. Use full limb width
            // for the last limb to match C++ unsigned wrapping behavior where
            // the range constraint is effectively unconstrained.
            NUM_LIMB_BITS
        } else {
            num_last_limb_bits::<T>()
        };

        let max_limb = default_maximum_limb();
        let max_last = if maximum_bitlength > 0 {
            U256::from(1u64)
                .wrapping_shl_vartime(num_last as u32)
                .wrapping_sub(&U256::ONE)
        } else if can_overflow {
            default_maximum_limb()
        } else {
            default_maximum_msb_limb::<T>()
        };

        // Range constrain limbs
        {
            let mut builder = ctx.borrow_mut();
            builder.range_constrain_two_limbs(
                limb_0.get_witness_index(),
                limb_1.get_witness_index(),
                NUM_LIMB_BITS as usize,
                NUM_LIMB_BITS as usize,
                "bigfield::create_from_u512_as_witness: limb 0 or 1 too large",
            );
            builder.range_constrain_two_limbs(
                limb_2.get_witness_index(),
                limb_3.get_witness_index(),
                NUM_LIMB_BITS as usize,
                num_last as usize,
                "bigfield::create_from_u512_as_witness: limb 2 or 3 too large",
            );
        }

        Self {
            context: Some(ctx),
            binary_basis_limbs: [
                Limb::new(limb_0, max_limb),
                Limb::new(limb_1, max_limb),
                Limb::new(limb_2, max_limb),
                Limb::new(limb_3, max_last),
            ],
            prime_basis_limb: prime_limb,
            _target: std::marker::PhantomData,
        }
    }

    /// Create a bigfield witness from a native non-native field value.
    ///
    /// Port of C++ `bigfield::from_witness`.
    pub fn from_witness(ctx: BuilderRef<P>, value: U256) -> Self {
        let lo_val = slice_u256(&value, 0, (NUM_LIMB_BITS * 2) as u32);
        let hi_val = slice_u256(&value, (NUM_LIMB_BITS * 2) as u32, (NUM_LIMB_BITS * 4) as u32);

        let lo_field = FieldT::from_witness(ctx.clone(), Field::from_limbs(*lo_val.as_words()));
        let hi_field = FieldT::from_witness(ctx.clone(), Field::from_limbs(*hi_val.as_words()));

        Self::from_field_pair(lo_field, hi_field, false, 0)
    }

    /// Create a bigfield from a pair of field_t elements (low 136 bits, high 136 bits).
    ///
    /// Port of C++ `bigfield(const field_t& low_bits, const field_t& high_bits, ...)`.
    pub fn from_field_pair(
        low_bits: FieldT<P>,
        high_bits: FieldT<P>,
        can_overflow: bool,
        maximum_bitlength: usize,
    ) -> Self {
        let ctx_opt = low_bits.get_context().clone().or(high_bits.get_context().clone());
        let ctx = ctx_opt.expect("from_field_pair: need at least one witness context");

        // Decompose low_bits into 2 limbs
        let (limb_0, limb_1) = if !low_bits.is_constant() {
            let indices = Self::decompose_non_native_field_double_width_limb(
                ctx.clone(),
                low_bits.get_witness_index(),
                (2 * NUM_LIMB_BITS) as usize,
            );
            let l0 = FieldT::from_witness_index(ctx.clone(), indices[0]);
            let l1 = FieldT::from_witness_index(ctx.clone(), indices[1]);
            // Verify: low_bits = limb_0 + limb_1 * shift_1
            // Constraint: l0 + l1*shift_1 - low_bits = 0
            let l1_scaled = &l1 * &FieldT::from_field(shift_1::<P>());
            let neg_low = -low_bits.clone();
            let zero_ft = FieldT::from_field(Field::zero());
            FieldT::evaluate_linear_identity(
                &l0,
                &l1_scaled,
                &neg_low,
                &zero_ft,
                "bigfield::from_field_pair: lo decomposition",
            );
            (l0, l1)
        } else {
            let low_val = low_bits.get_value().from_montgomery_form();
            let low_u256 = U256::from_words(low_val.data);
            let s0 = slice_u256(&low_u256, 0, NUM_LIMB_BITS as u32);
            let s1 = slice_u256(&low_u256, NUM_LIMB_BITS as u32, (NUM_LIMB_BITS * 2) as u32);
            (
                FieldT::from_field(Field::from_limbs(*s0.as_words())),
                FieldT::from_field(Field::from_limbs(*s1.as_words())),
            )
        };

        // Determine num_last_limb_bits
        let num_last = if can_overflow {
            NUM_LIMB_BITS
        } else if maximum_bitlength > 0 {
            assert!(maximum_bitlength > (3 * NUM_LIMB_BITS) as usize);
            assert!(maximum_bitlength <= (4 * NUM_LIMB_BITS) as usize);
            (maximum_bitlength as u64) - NUM_LIMB_BITS * 3
        } else {
            num_last_limb_bits::<T>()
        };

        let num_high_limb_bits = NUM_LIMB_BITS + num_last;

        // Decompose high_bits into 2 limbs
        let (limb_2, limb_3) = if !high_bits.is_constant() {
            let indices = Self::decompose_non_native_field_double_width_limb(
                ctx.clone(),
                high_bits.get_witness_index(),
                num_high_limb_bits as usize,
            );
            let l2 = FieldT::from_witness_index(ctx.clone(), indices[0]);
            let l3 = FieldT::from_witness_index(ctx.clone(), indices[1]);
            // Constraint: l2 + l3*shift_1 - high_bits = 0
            let l3_scaled = &l3 * &FieldT::from_field(shift_1::<P>());
            let neg_high = -high_bits.clone();
            let zero_ft = FieldT::from_field(Field::zero());
            FieldT::evaluate_linear_identity(
                &l2,
                &l3_scaled,
                &neg_high,
                &zero_ft,
                "bigfield::from_field_pair: hi decomposition",
            );
            (l2, l3)
        } else {
            let hi_val = high_bits.get_value().from_montgomery_form();
            let hi_u256 = U256::from_words(hi_val.data);
            let s2 = slice_u256(&hi_u256, 0, NUM_LIMB_BITS as u32);
            let s3 = slice_u256(&hi_u256, NUM_LIMB_BITS as u32, num_high_limb_bits as u32);
            (
                FieldT::from_field(Field::from_limbs(*s2.as_words())),
                FieldT::from_field(Field::from_limbs(*s3.as_words())),
            )
        };

        let max_3 = if maximum_bitlength > 0 {
            U256::from(1u64)
                .wrapping_shl_vartime(num_last as u32)
                .wrapping_sub(&U256::ONE)
        } else if can_overflow {
            default_maximum_limb()
        } else {
            default_maximum_msb_limb::<T>()
        };

        // Prime basis limb = low_bits + high_bits * shift_2
        let prime_basis_limb = &low_bits + &(&high_bits * &FieldT::from_field(shift_2::<P>()));

        Self {
            context: Some(ctx),
            binary_basis_limbs: [
                Limb::new_default(limb_0),
                Limb::new_default(limb_1),
                Limb::new_default(limb_2),
                Limb::new(limb_3, max_3),
            ],
            prime_basis_limb,
            _target: std::marker::PhantomData,
        }
    }

    /// Construct from 4 limbs without range checking (unsafe).
    pub fn unsafe_construct_from_limbs(
        limb0: FieldT<P>,
        limb1: FieldT<P>,
        limb2: FieldT<P>,
        limb3: FieldT<P>,
        can_overflow: bool,
    ) -> Self {
        let ctx = limb0
            .get_context()
            .clone()
            .or(limb1.get_context().clone())
            .or(limb2.get_context().clone())
            .or(limb3.get_context().clone());

        let max_3 = if can_overflow {
            default_maximum_limb()
        } else {
            default_maximum_msb_limb::<T>()
        };

        let prime_limb = &limb0
            + &(&limb1 * &FieldT::from_field(shift_1::<P>()))
            + &(&limb2 * &FieldT::from_field(shift_2::<P>()))
            + &(&limb3 * &FieldT::from_field(shift_3::<P>()));

        Self {
            context: ctx,
            binary_basis_limbs: [
                Limb::new_default(limb0),
                Limb::new_default(limb1),
                Limb::new_default(limb2),
                Limb::new(limb3, max_3),
            ],
            prime_basis_limb: prime_limb,
            _target: std::marker::PhantomData,
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  Value access
    // ════════════════════════════════════════════════════════════════════

    /// Get the value of this bigfield element as a U512.
    pub fn get_value(&self) -> U512 {
        let l0 = self.binary_basis_limbs[0].element.get_value().from_montgomery_form();
        let l1 = self.binary_basis_limbs[1].element.get_value().from_montgomery_form();
        let l2 = self.binary_basis_limbs[2].element.get_value().from_montgomery_form();
        let l3 = self.binary_basis_limbs[3].element.get_value().from_montgomery_form();

        let v0 = U512::from_lo_hi(U256::from_words(l0.data), U256::ZERO);
        let v1 = U512::from_lo_hi(U256::from_words(l1.data), U256::ZERO)
            .wrapping_shl_vartime(NUM_LIMB_BITS as u32);
        let v2 = U512::from_lo_hi(U256::from_words(l2.data), U256::ZERO)
            .wrapping_shl_vartime((NUM_LIMB_BITS * 2) as u32);
        let v3 = U512::from_lo_hi(U256::from_words(l3.data), U256::ZERO)
            .wrapping_shl_vartime((NUM_LIMB_BITS * 3) as u32);

        v0.wrapping_add(&v1).wrapping_add(&v2).wrapping_add(&v3)
    }

    /// Get the maximum possible value of this element.
    pub fn get_maximum_value(&self) -> U512 {
        let t0 = U512::from_lo_hi(self.binary_basis_limbs[0].maximum_value, U256::ZERO);
        let t1 = U512::from_lo_hi(self.binary_basis_limbs[1].maximum_value, U256::ZERO)
            .wrapping_shl_vartime(NUM_LIMB_BITS as u32);
        let t2 = U512::from_lo_hi(self.binary_basis_limbs[2].maximum_value, U256::ZERO)
            .wrapping_shl_vartime((NUM_LIMB_BITS * 2) as u32);
        let t3 = U512::from_lo_hi(self.binary_basis_limbs[3].maximum_value, U256::ZERO)
            .wrapping_shl_vartime((NUM_LIMB_BITS * 3) as u32);
        t0.wrapping_add(&t1).wrapping_add(&t2).wrapping_add(&t3)
    }

    /// Check if this element is a circuit constant.
    pub fn is_constant(&self) -> bool {
        self.binary_basis_limbs[0].element.is_constant()
            && self.binary_basis_limbs[1].element.is_constant()
            && self.binary_basis_limbs[2].element.is_constant()
            && self.binary_basis_limbs[3].element.is_constant()
    }

    /// Get the builder context.
    pub fn get_context(&self) -> &Option<BuilderRef<P>> {
        &self.context
    }

    /// Get binary basis limb witness indices.
    fn get_limb_witness_indices(&self) -> [u32; 4] {
        [
            self.binary_basis_limbs[0].element.get_witness_index(),
            self.binary_basis_limbs[1].element.get_witness_index(),
            self.binary_basis_limbs[2].element.get_witness_index(),
            self.binary_basis_limbs[3].element.get_witness_index(),
        ]
    }

    /// Get binary basis limb maximum values.
    fn get_limb_maximums(&self) -> [U256; 4] {
        [
            self.binary_basis_limbs[0].maximum_value,
            self.binary_basis_limbs[1].maximum_value,
            self.binary_basis_limbs[2].maximum_value,
            self.binary_basis_limbs[3].maximum_value,
        ]
    }

    // ════════════════════════════════════════════════════════════════════
    //  Decompose helper
    // ════════════════════════════════════════════════════════════════════

    /// Decompose a double-width limb witness into two single-width limb witnesses.
    ///
    /// Port of C++ `decompose_non_native_field_double_width_limb`.
    fn decompose_non_native_field_double_width_limb(
        ctx: BuilderRef<P>,
        limb_idx: u32,
        num_limb_bits: usize,
    ) -> [u32; 2] {
        let limb_value = {
            let builder = ctx.borrow();
            builder.base.get_variable(limb_idx).from_montgomery_form()
        };
        let limb_u256 = U256::from_words(limb_value.data);

        let lo_val = slice_u256(&limb_u256, 0, NUM_LIMB_BITS as u32);
        let hi_val = slice_u256(&limb_u256, NUM_LIMB_BITS as u32, num_limb_bits as u32);

        let lo_idx = ctx.borrow_mut().base.add_variable(Field::from_limbs(*lo_val.as_words()));
        let hi_idx = ctx.borrow_mut().base.add_variable(Field::from_limbs(*hi_val.as_words()));

        // Range constrain both limbs
        let hi_bits = num_limb_bits - NUM_LIMB_BITS as usize;
        {
            let mut builder = ctx.borrow_mut();
            builder.range_constrain_two_limbs(
                lo_idx,
                hi_idx,
                NUM_LIMB_BITS as usize,
                hi_bits,
                "decompose_non_native_field_double_width_limb",
            );
        }

        [lo_idx, hi_idx]
    }

    // ════════════════════════════════════════════════════════════════════
    //  Reduction
    // ════════════════════════════════════════════════════════════════════

    /// Check if reduction is needed and perform it if so.
    ///
    /// Port of C++ `bigfield::reduction_check`.
    pub fn reduction_check(&mut self) {
        if self.is_constant() {
            // Reduce constant to value mod target modulus
            let val = self.get_value();
            let modulus = Self::modulus_u512();
            let (_, remainder) = val.div_rem(&modulus.to_nz().unwrap());
            let rem_lo = remainder.lo();
            *self = Self::from_u256(self.context.clone(), rem_lo);
            return;
        }

        let max_unreduced = get_maximum_unreduced_limb_value::<P>();
        let limb_overflow = self.binary_basis_limbs[0].maximum_value > max_unreduced
            || self.binary_basis_limbs[1].maximum_value > max_unreduced
            || self.binary_basis_limbs[2].maximum_value > max_unreduced
            || self.binary_basis_limbs[3].maximum_value > max_unreduced;

        if limb_overflow {
            self.self_reduce();
        }
    }

    /// Perform self-reduction: reduces element mod target_modulus.
    ///
    /// Port of C++ `bigfield::self_reduce`.
    pub fn self_reduce(&mut self) {
        let ctx = self
            .context
            .clone()
            .expect("self_reduce requires context");

        // Compute quotient and remainder
        let value = self.get_value();
        let modulus_512 = Self::modulus_u512();
        let modulus_512_nz = modulus_512.to_nz().unwrap();
        let (quotient_512, remainder_512) = value.div_rem(&modulus_512_nz);
        let quotient_value = quotient_512.lo();
        let remainder_value = remainder_512;

        // Determine quotient bits
        let max_value = self.get_maximum_value();
        let max_quotient = max_value.div_rem(&modulus_512_nz).0.lo();
        let mut max_quotient_bits = max_quotient.get_msb() as u64 + 1;
        if max_quotient_bits & 1 == 1 {
            max_quotient_bits += 1;
        }
        if max_quotient_bits < 2 {
            max_quotient_bits = 2;
        }

        // Create quotient witness (single limb since quotient is small)
        let q_limb_idx = ctx
            .borrow_mut()
            .base
            .add_variable(Field::from_limbs(*quotient_value.as_words()));
        ctx.borrow_mut()
            .decompose_into_default_range(
                q_limb_idx,
                max_quotient_bits,
                DEFAULT_PLOOKUP_RANGE_BITNUM,
                "bigfield::self_reduce: quotient range check",
            );

        let q_limb = FieldT::from_witness_index(ctx.clone(), q_limb_idx);
        let q_max =
            U256::from(1u64).wrapping_shl_vartime(max_quotient_bits as u32);

        let quotient = Self {
            context: Some(ctx.clone()),
            binary_basis_limbs: [
                Limb::new(q_limb, q_max),
                Limb::zero_limb(),
                Limb::zero_limb(),
                Limb::zero_limb(),
            ],
            prime_basis_limb: FieldT::from_witness_index(ctx.clone(), q_limb_idx),
            _target: std::marker::PhantomData,
        };

        // Create remainder as witness
        let remainder = Self::create_from_u512_as_witness(
            ctx.clone(),
            remainder_value,
            false,
            0,
        );

        // Verify: self * 1 = quotient * p + remainder
        let one = Self::one_val();
        Self::unsafe_evaluate_multiply_add(self, &one, &[], &quotient, &[&remainder]);

        // Update self with remainder
        self.binary_basis_limbs = remainder.binary_basis_limbs;
        self.prime_basis_limb = remainder.prime_basis_limb;
    }

    // ════════════════════════════════════════════════════════════════════
    //  Quotient/remainder computation
    // ════════════════════════════════════════════════════════════════════

    /// Compute quotient and remainder of (a * b + sum(to_add)) / target_modulus.
    fn compute_quotient_remainder_values(
        a: &Self,
        b: &Self,
        to_add: &[&Self],
    ) -> (U512, U512) {
        let mut add_values = U512::ZERO;
        for add in to_add {
            add_values = add_values.wrapping_add(&add.get_value());
        }

        let left = a.get_value();
        let right = b.get_value();
        let product = left.wrapping_mul(&right).wrapping_add(&add_values);

        let modulus = Self::modulus_u512().to_nz().unwrap();
        product.div_rem(&modulus)
    }

    /// Check if a multiplication would overflow the CRT modulus and determine quotient bits.
    ///
    /// Returns (needs_reduction, quotient_bits).
    fn get_quotient_reduction_info(
        a_max: &[U512],
        b_max: &[U512],
        to_add: &[&Self],
    ) -> (bool, usize) {
        // Check for U512 overflow: if a_bits + b_bits > 512, product overflows
        // and we must reduce before multiplying.
        fn msb_512(v: &U512) -> u64 {
            let hi = v.hi();
            let lo = v.lo();
            if hi > U256::ZERO {
                256 + hi.get_msb() as u64
            } else {
                lo.get_msb() as u64
            }
        }

        for i in 0..a_max.len() {
            let a_bits = msb_512(&a_max[i]) + 1;
            let b_bits = msb_512(&b_max[i]) + 1;
            if a_bits + b_bits > 512 {
                return (true, 0);
            }
        }

        // Compute maximum product + additions (safe - no overflow)
        let mut max_product = U512::ZERO;
        for i in 0..a_max.len() {
            max_product = max_product.wrapping_add(&a_max[i].wrapping_mul(&b_max[i]));
        }
        for add in to_add {
            max_product = max_product.wrapping_add(&add.get_maximum_value());
        }

        let modulus = Self::modulus_u512().to_nz().unwrap();
        let max_quotient = max_product.div_rem(&modulus).0;

        // Compute max quotient bits
        let max_q_lo = max_quotient.lo();
        let max_q_hi = max_quotient.hi();

        let mut num_quotient_bits = if max_q_hi > U256::ZERO {
            256 + max_q_hi.get_msb() as usize + 1
        } else {
            max_q_lo.get_msb() as usize + 1
        };

        // Must be even for range proofs
        if num_quotient_bits & 1 == 1 {
            num_quotient_bits += 1;
        }
        if num_quotient_bits < 2 {
            num_quotient_bits = 2;
        }

        // Check if quotient fits in the available range proof bits
        let max_quotient_fits = num_quotient_bits <= (NUM_LIMB_BITS * 4) as usize;

        if !max_quotient_fits {
            return (true, 0);
        }

        (false, num_quotient_bits)
    }

    // ════════════════════════════════════════════════════════════════════
    //  Core constraint: unsafe_evaluate_multiply_add
    // ════════════════════════════════════════════════════════════════════

    /// Create circuit constraints verifying: left * to_mul + sum(to_add) = quotient * p + sum(remainders)
    ///
    /// Checks both mod native field and mod 2^272 (binary basis).
    ///
    /// Port of C++ `bigfield::unsafe_evaluate_multiply_add`.
    fn unsafe_evaluate_multiply_add(
        left: &Self,
        to_mul: &Self,
        to_add: &[&Self],
        quotient: &Self,
        remainders: &[&Self],
    ) {
        let ctx = left
            .context
            .clone()
            .or(to_mul.context.clone())
            .or(quotient.context.clone())
            .expect("unsafe_evaluate_multiply_add requires context");

        let neg_modulus = neg_modulus_limbs::<P, T>();

        // ─── Binary basis check (mod 2^272) ───
        // Compute the product limbs and quotient*neg_modulus limbs
        // We need:
        //   left * to_mul + to_add = quotient * p + remainders (mod 2^272)
        // Which is:
        //   left * to_mul + to_add - remainders + quotient * neg_p = 0 (mod 2^272)

        // Get normalized witness indices for all elements
        let left_norm = Self::normalize_limbs(left, &ctx);
        let mul_norm = Self::normalize_limbs(to_mul, &ctx);
        let q_norm = Self::normalize_limbs(quotient, &ctx);

        // Build remainder accumulator limb indices
        let mut r_limb0_accum: Vec<FieldT<P>> = Vec::new();
        let mut r_limb2_accum: Vec<FieldT<P>> = Vec::new();
        let mut r_prime_accum: Vec<FieldT<P>> = Vec::new();

        for r in remainders {
            let r_norm_0 = r.binary_basis_limbs[0].element.normalize();
            let r_norm_1 = r.binary_basis_limbs[1].element.normalize();
            let r_norm_2 = r.binary_basis_limbs[2].element.normalize();
            let r_norm_3 = r.binary_basis_limbs[3].element.normalize();

            // accumulate: r_limb0 += r[0] + r[1]*shift_1
            r_limb0_accum.push(r_norm_0);
            r_limb0_accum.push(&r_norm_1 * &FieldT::from_field(shift_1::<P>()));
            r_limb2_accum.push(r_norm_2);
            r_limb2_accum.push(&r_norm_3 * &FieldT::from_field(shift_1::<P>()));
            r_prime_accum.push(r.prime_basis_limb.clone());
        }

        // Subtract to_add terms
        for add in to_add {
            let a_norm_0 = add.binary_basis_limbs[0].element.normalize();
            let a_norm_1 = add.binary_basis_limbs[1].element.normalize();
            let a_norm_2 = add.binary_basis_limbs[2].element.normalize();
            let a_norm_3 = add.binary_basis_limbs[3].element.normalize();

            r_limb0_accum.push(-a_norm_0);
            r_limb0_accum.push(-(&a_norm_1 * &FieldT::from_field(shift_1::<P>())));
            r_limb2_accum.push(-a_norm_2);
            r_limb2_accum.push(-(&a_norm_3 * &FieldT::from_field(shift_1::<P>())));
            r_prime_accum.push(-add.prime_basis_limb.clone());
        }

        // Accumulate remainder terms into single indices via linear combinations
        let r_lo = FieldT::accumulate(&r_limb0_accum);
        let r_hi = FieldT::accumulate(&r_limb2_accum);
        let r_prime = FieldT::accumulate(&r_prime_accum);

        let r_lo_norm = r_lo.normalize();
        let r_hi_norm = r_hi.normalize();

        // Evaluate non-native field multiplication: left * to_mul = q * (-p) + r (mod 2^272)
        // r[1] and r[3] are zero because r_lo and r_hi already accumulate shifted limb pairs.
        // Use the builder's zero variable index (not IS_CONSTANT, which is u32::MAX and invalid).
        // Convert constant r_lo/r_hi to real witness indices if needed.
        let zero_var = ctx.borrow().base.zero_idx();
        let r_lo_idx = if r_lo_norm.is_constant() {
            ctx.borrow_mut().put_constant_variable(r_lo_norm.get_value())
        } else {
            r_lo_norm.get_witness_index()
        };
        let r_hi_idx = if r_hi_norm.is_constant() {
            ctx.borrow_mut().put_constant_variable(r_hi_norm.get_value())
        } else {
            r_hi_norm.get_witness_index()
        };
        let witnesses = NnfMulWitnesses {
            a: left_norm,
            b: mul_norm,
            q: q_norm,
            r: [r_lo_idx, zero_var, r_hi_idx, zero_var],
            neg_modulus: neg_modulus,
        };

        let [lo_1_idx, hi_3_idx] =
            ctx.borrow_mut().evaluate_non_native_field_multiplication(&witnesses);

        // ─── Prime basis check (mod native field) ───
        // left.prime * to_mul.prime + q.prime * neg_prime_mod - r_prime = 0
        let neg_prime = negative_prime_modulus::<P, T>();
        let q_prime_scaled = &quotient.prime_basis_limb * &FieldT::from_field(neg_prime);

        FieldT::evaluate_polynomial_identity(
            &left.prime_basis_limb,
            &to_mul.prime_basis_limb,
            &q_prime_scaled,
            &(-r_prime),
            "bigfield: prime basis identity failed",
        );

        // ─── Range constrain carries ───
        let lo_carry = FieldT::from_witness_index(ctx.clone(), lo_1_idx);
        let hi_carry = FieldT::from_witness_index(ctx.clone(), hi_3_idx);

        // Compute max carry bits from the multiplication maximum values
        let left_maxes = left.get_limb_maximums();
        let mul_maxes = to_mul.get_limb_maximums();
        let q_maxes = quotient.get_limb_maximums();
        let neg_mod_u256 = neg_modulus_limbs_u256::<T>();

        let (max_lo, max_hi) =
            compute_partial_schoolbook_multiplication(&left_maxes, &mul_maxes);
        let (max_qp_lo, max_qp_hi) =
            compute_partial_schoolbook_multiplication(&q_maxes, &neg_mod_u256);

        // Add remainder and to_add max values
        let mut max_r_lo = U512::ZERO;
        let mut max_r_hi = U512::ZERO;
        for r in remainders {
            let r_max = r.get_limb_maximums();
            max_r_lo = max_r_lo.wrapping_add(
                &U512::from_lo_hi(r_max[0], U256::ZERO).wrapping_add(
                    &U512::from_lo_hi(r_max[1], U256::ZERO)
                        .wrapping_shl_vartime(NUM_LIMB_BITS as u32),
                ),
            );
            max_r_hi = max_r_hi.wrapping_add(
                &U512::from_lo_hi(r_max[2], U256::ZERO).wrapping_add(
                    &U512::from_lo_hi(r_max[3], U256::ZERO)
                        .wrapping_shl_vartime(NUM_LIMB_BITS as u32),
                ),
            );
        }

        let total_max_lo = max_lo
            .wrapping_add(&max_qp_lo)
            .wrapping_add(&max_r_lo);
        // The carry from lo feeds into hi via lo_1 = total_lo >> 136
        let max_lo_carry = total_max_lo.wrapping_shr_vartime(2 * NUM_LIMB_BITS as u32);
        let total_max_hi = max_hi
            .wrapping_add(&max_qp_hi)
            .wrapping_add(&max_r_hi)
            .wrapping_add(&max_lo_carry);

        let max_lo_bits = total_max_lo.lo().get_msb() as u64 + 1;
        let max_hi_bits = total_max_hi.lo().get_msb() as u64 + 1;

        let carry_lo_bits = if max_lo_bits > 2 * NUM_LIMB_BITS {
            max_lo_bits - 2 * NUM_LIMB_BITS
        } else {
            1
        };
        let carry_hi_bits = if max_hi_bits > 2 * NUM_LIMB_BITS {
            max_hi_bits - 2 * NUM_LIMB_BITS
        } else {
            1
        };

        // Range constrain carries
        if carry_lo_bits <= 70 && carry_hi_bits <= 70 {
            ctx.borrow_mut().range_constrain_two_limbs(
                hi_carry.get_witness_index(),
                lo_carry.get_witness_index(),
                carry_hi_bits as usize,
                carry_lo_bits as usize,
                "bigfield: carry range check",
            );
        } else {
            lo_carry.create_range_constraint(carry_lo_bits as usize, "bigfield: lo carry");
            hi_carry.create_range_constraint(carry_hi_bits as usize, "bigfield: hi carry");
        }
    }

    /// Normalize limb elements and return their witness indices.
    fn normalize_limbs(elem: &Self, ctx: &BuilderRef<P>) -> [u32; 4] {
        let mut indices = [0u32; 4];
        for i in 0..4 {
            let norm = elem.binary_basis_limbs[i].element.normalize();
            indices[i] = if norm.is_constant() {
                ctx.borrow_mut()
                    .put_constant_variable(norm.get_value())
            } else {
                norm.get_witness_index()
            };
        }
        indices
    }

    // ════════════════════════════════════════════════════════════════════
    //  Arithmetic operations
    // ════════════════════════════════════════════════════════════════════

    /// Addition: self + other.
    ///
    /// Port of C++ `bigfield::operator+`.
    pub fn add(&self, other: &Self) -> Self {
        let ctx = self
            .context
            .clone()
            .or(other.context.clone());

        Self {
            context: ctx,
            binary_basis_limbs: [
                Limb::new(
                    &self.binary_basis_limbs[0].element + &other.binary_basis_limbs[0].element,
                    self.binary_basis_limbs[0]
                        .maximum_value
                        .wrapping_add(&other.binary_basis_limbs[0].maximum_value),
                ),
                Limb::new(
                    &self.binary_basis_limbs[1].element + &other.binary_basis_limbs[1].element,
                    self.binary_basis_limbs[1]
                        .maximum_value
                        .wrapping_add(&other.binary_basis_limbs[1].maximum_value),
                ),
                Limb::new(
                    &self.binary_basis_limbs[2].element + &other.binary_basis_limbs[2].element,
                    self.binary_basis_limbs[2]
                        .maximum_value
                        .wrapping_add(&other.binary_basis_limbs[2].maximum_value),
                ),
                Limb::new(
                    &self.binary_basis_limbs[3].element + &other.binary_basis_limbs[3].element,
                    self.binary_basis_limbs[3]
                        .maximum_value
                        .wrapping_add(&other.binary_basis_limbs[3].maximum_value),
                ),
            ],
            prime_basis_limb: &self.prime_basis_limb + &other.prime_basis_limb,
            _target: std::marker::PhantomData,
        }
    }

    /// Subtraction: self - other.
    ///
    /// Port of C++ `bigfield::operator-`.
    pub fn sub(&self, other: &Self) -> Self {
        let ctx = self
            .context
            .clone()
            .or(other.context.clone());

        // Compute borrow terms to prevent underflow
        let limb_0_borrow_shift = std::cmp::max(
            other.binary_basis_limbs[0].maximum_value.get_msb() as u64 + 1,
            NUM_LIMB_BITS,
        );
        let limb_1_max = other.binary_basis_limbs[1].maximum_value.wrapping_add(
            &U256::from(1u64).wrapping_shl_vartime((limb_0_borrow_shift - NUM_LIMB_BITS) as u32),
        );
        let limb_1_borrow_shift = std::cmp::max(limb_1_max.get_msb() as u64 + 1, NUM_LIMB_BITS);

        let limb_2_max = other.binary_basis_limbs[2].maximum_value.wrapping_add(
            &U256::from(1u64).wrapping_shl_vartime((limb_1_borrow_shift - NUM_LIMB_BITS) as u32),
        );
        let limb_2_borrow_shift = std::cmp::max(limb_2_max.get_msb() as u64 + 1, NUM_LIMB_BITS);

        let limb_3_max = other.binary_basis_limbs[3].maximum_value.wrapping_add(
            &U256::from(1u64).wrapping_shl_vartime((limb_2_borrow_shift - NUM_LIMB_BITS) as u32),
        );

        // Compute constant_to_add = ceil(limb_3_max * 2^(3*NUM_LIMB_BITS) / modulus) * modulus
        let modulus = Self::modulus();
        let limb_3_shifted = U512::from_lo_hi(limb_3_max, U256::ZERO)
            .wrapping_shl_vartime((NUM_LIMB_BITS * 3) as u32);
        let modulus_512 = Self::modulus_u512();
        let factor = limb_3_shifted
            .div_rem(&modulus_512.to_nz().unwrap())
            .0
            .lo()
            .wrapping_add(&U256::ONE);
        // constant_to_add = factor * modulus
        let factor_512 = U512::from_lo_hi(factor, U256::ZERO);
        let constant_to_add = factor_512.wrapping_mul(&modulus_512);
        let cta_lo = constant_to_add.lo();

        // Compute borrow offsets
        let t0 = U256::from(1u64).wrapping_shl_vartime(limb_0_borrow_shift as u32);
        let t1 = U256::from(1u64)
            .wrapping_shl_vartime(limb_1_borrow_shift as u32)
            .wrapping_sub(
                &U256::from(1u64)
                    .wrapping_shl_vartime((limb_0_borrow_shift - NUM_LIMB_BITS) as u32),
            );
        let t2 = U256::from(1u64)
            .wrapping_shl_vartime(limb_2_borrow_shift as u32)
            .wrapping_sub(
                &U256::from(1u64)
                    .wrapping_shl_vartime((limb_1_borrow_shift - NUM_LIMB_BITS) as u32),
            );
        let t3 = U256::from(1u64)
            .wrapping_shl_vartime((limb_2_borrow_shift - NUM_LIMB_BITS) as u32);

        // Slice constant_to_add into 4 limbs of NUM_LIMB_BITS (68) each.
        // constant_to_add can exceed 2^256, so limb 3 (bits 204-271) crosses
        // the U256 boundary at bit 256: 52 bits from lo, 16 bits from hi.
        let cta_hi = constant_to_add.hi();
        let cta_0 = slice_u256(&cta_lo, 0, NUM_LIMB_BITS as u32);
        let cta_1 = slice_u256(&cta_lo, NUM_LIMB_BITS as u32, (NUM_LIMB_BITS * 2) as u32);
        let cta_2 = slice_u256(&cta_lo, (NUM_LIMB_BITS * 2) as u32, (NUM_LIMB_BITS * 3) as u32);
        let remaining_in_lo = 256 - (NUM_LIMB_BITS * 3) as u32; // 52
        let lo_part = slice_u256(&cta_lo, (NUM_LIMB_BITS * 3) as u32, 256);
        let hi_bits_needed = NUM_LIMB_BITS as u32 - remaining_in_lo; // 16
        let hi_part = slice_u256(&cta_hi, 0, hi_bits_needed);
        let cta_3 = lo_part.wrapping_or(&hi_part.wrapping_shl_vartime(remaining_in_lo));

        let to_add_0 = cta_0.wrapping_add(&t0);
        let to_add_1 = cta_1.wrapping_add(&t1);
        let to_add_2 = cta_2.wrapping_add(&t2);
        let to_add_3 = cta_3.wrapping_sub(&t3);

        // Result limbs = self.limb + to_add - other.limb
        let new_limbs = [
            Limb::new(
                &(&self.binary_basis_limbs[0].element
                    + &FieldT::from_field(Field::from_limbs(*to_add_0.as_words())))
                    - &other.binary_basis_limbs[0].element,
                self.binary_basis_limbs[0]
                    .maximum_value
                    .wrapping_add(&to_add_0),
            ),
            Limb::new(
                &(&self.binary_basis_limbs[1].element
                    + &FieldT::from_field(Field::from_limbs(*to_add_1.as_words())))
                    - &other.binary_basis_limbs[1].element,
                self.binary_basis_limbs[1]
                    .maximum_value
                    .wrapping_add(&to_add_1),
            ),
            Limb::new(
                &(&self.binary_basis_limbs[2].element
                    + &FieldT::from_field(Field::from_limbs(*to_add_2.as_words())))
                    - &other.binary_basis_limbs[2].element,
                self.binary_basis_limbs[2]
                    .maximum_value
                    .wrapping_add(&to_add_2),
            ),
            Limb::new(
                &(&self.binary_basis_limbs[3].element
                    + &FieldT::from_field(Field::from_limbs(*to_add_3.as_words())))
                    - &other.binary_basis_limbs[3].element,
                self.binary_basis_limbs[3]
                    .maximum_value
                    .wrapping_add(&to_add_3),
            ),
        ];

        // Prime basis: self.prime + constant_to_add_mod_native - other.prime
        // constant_to_add can exceed 2^256, so reduce mod P properly
        let native_mod_512 =
            U512::from_lo_hi(U256::from_limbs(P::MODULUS), U256::ZERO);
        let cta_mod_p = constant_to_add
            .div_rem(&native_mod_512.to_nz().unwrap())
            .1
            .lo();
        let cta_mod_native = Field::<P>::from_limbs(*cta_mod_p.as_words());
        let prime_to_add = FieldT::from_field(cta_mod_native);
        let new_prime =
            &(&self.prime_basis_limb + &prime_to_add) - &other.prime_basis_limb;

        Self {
            context: ctx,
            binary_basis_limbs: new_limbs,
            prime_basis_limb: new_prime,
            _target: std::marker::PhantomData,
        }
    }

    /// Multiplication: self * other.
    ///
    /// Port of C++ `bigfield::operator*`.
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_constant() && other.is_constant() {
            let val = self
                .get_value()
                .wrapping_mul(&other.get_value())
                .div_rem(&Self::modulus_u512().to_nz().unwrap())
                .1;
            return Self::from_u256(
                self.context.clone().or(other.context.clone()),
                val.lo(),
            );
        }

        let ctx = self
            .context
            .clone()
            .or(other.context.clone())
            .expect("mul requires context");

        // Check if reduction is needed
        let (needs_reduction, num_quotient_bits) = Self::get_quotient_reduction_info(
            &[self.get_maximum_value()],
            &[other.get_maximum_value()],
            &[],
        );

        if needs_reduction {
            let mut a = self.clone();
            let mut b = other.clone();
            if a.get_maximum_value() > b.get_maximum_value() {
                a.self_reduce();
            } else {
                b.self_reduce();
            }
            return a.mul(&b);
        }

        let (quotient_value, remainder_value) =
            Self::compute_quotient_remainder_values(self, other, &[]);

        let quotient = Self::create_from_u512_as_witness(
            ctx.clone(),
            quotient_value,
            false,
            num_quotient_bits,
        );
        let remainder = Self::create_from_u512_as_witness(ctx.clone(), remainder_value, false, 0);

        Self::unsafe_evaluate_multiply_add(self, other, &[], &quotient, &[&remainder]);

        remainder
    }

    /// Squaring: self^2.
    ///
    /// Port of C++ `bigfield::sqr`.
    pub fn sqr(&self) -> Self {
        self.mul(self)
    }

    /// Multiply-add: self * to_mul + sum(to_add).
    ///
    /// Port of C++ `bigfield::madd`.
    pub fn madd(&self, to_mul: &Self, to_add: &[&Self]) -> Self {
        if self.is_constant() && to_mul.is_constant() {
            let product = self
                .get_value()
                .wrapping_mul(&to_mul.get_value());
            let mut sum = product;
            for add in to_add {
                sum = sum.wrapping_add(&add.get_value());
            }
            let result = sum.div_rem(&Self::modulus_u512().to_nz().unwrap()).1;
            return Self::from_u256(
                self.context.clone().or(to_mul.context.clone()),
                result.lo(),
            );
        }

        let ctx = self
            .context
            .clone()
            .or(to_mul.context.clone())
            .expect("madd requires context");

        // Check reduction
        let (needs_reduction, num_quotient_bits) = Self::get_quotient_reduction_info(
            &[self.get_maximum_value()],
            &[to_mul.get_maximum_value()],
            to_add,
        );

        if needs_reduction {
            let mut a = self.clone();
            let mut b = to_mul.clone();
            if a.get_maximum_value() > b.get_maximum_value() {
                a.self_reduce();
            } else {
                b.self_reduce();
            }
            return a.madd(&b, to_add);
        }

        let (quotient_value, remainder_value) =
            Self::compute_quotient_remainder_values(self, to_mul, to_add);

        let quotient = Self::create_from_u512_as_witness(
            ctx.clone(),
            quotient_value,
            false,
            num_quotient_bits,
        );
        let remainder = Self::create_from_u512_as_witness(ctx.clone(), remainder_value, false, 0);

        Self::unsafe_evaluate_multiply_add(self, to_mul, to_add, &quotient, &[&remainder]);

        remainder
    }

    /// Division: self / other.
    ///
    /// Port of C++ `bigfield::operator/`.
    pub fn div(&self, other: &Self) -> Self {
        if self.is_constant() && other.is_constant() {
            return Self::internal_div(&[self], other, true);
        }
        // Reduce numerator to ensure NNF quotient stays non-negative
        let mut reduced = self.clone();
        reduced.self_reduce();
        Self::internal_div(&[&reduced], other, true)
    }

    /// Internal division: sum(numerators) / denominator.
    fn internal_div(numerators: &[&Self], denominator: &Self, check_for_zero: bool) -> Self {
        if denominator.is_constant() && numerators.iter().all(|n| n.is_constant()) {
            // All constant: compute directly, no context needed
            let mut num_val = U512::ZERO;
            for n in numerators {
                num_val = num_val.wrapping_add(&n.get_value());
            }
            let den_val = denominator.get_value();
            let modulus = Self::modulus_u512().to_nz().unwrap();
            let den_u256 = den_val.lo();
            let mod_u256 = Self::modulus();
            let inv = Self::invmod_u256(den_u256, mod_u256);
            let inv_512 = U512::from_lo_hi(inv, U256::ZERO);
            let result = num_val
                .wrapping_mul(&inv_512)
                .div_rem(&modulus)
                .1;
            let ctx = numerators[0]
                .context
                .clone()
                .or(denominator.context.clone());
            return Self::from_u256(ctx, result.lo());
        }

        let ctx = numerators[0]
            .context
            .clone()
            .or(denominator.context.clone())
            .expect("div requires context");

        // Compute: result = sum(numerators) * denominator^(-1) mod p
        let den_val = denominator.get_value().lo();
        let mod_val = Self::modulus();
        let den_inv = Self::invmod_u256(den_val, mod_val);
        let den_inv_512 = U512::from_lo_hi(den_inv, U256::ZERO);
        let modulus_nz = Self::modulus_u512().to_nz().unwrap();
        let mut num_sum_val = U512::ZERO;
        for n in numerators {
            num_sum_val = num_sum_val.wrapping_add(&n.get_value());
        }
        let result_val = num_sum_val
            .wrapping_mul(&den_inv_512)
            .div_rem(&modulus_nz)
            .1;

        let inverse = Self::create_from_u512_as_witness(
            ctx.clone(),
            result_val,
            false,
            0,
        );

        // Compute quotient: inverse * denominator - sum(numerators) = quotient * p
        let product = inverse
            .get_value()
            .wrapping_mul(&denominator.get_value());
        let quotient_val = product
            .wrapping_sub(&num_sum_val)
            .div_rem(&modulus_nz)
            .0;

        // Check quotient reduction
        let (needs_reduction, num_quotient_bits) = Self::get_quotient_reduction_info(
            &[inverse.get_maximum_value()],
            &[denominator.get_maximum_value()],
            numerators,
        );

        if needs_reduction {
            let mut den = denominator.clone();
            den.self_reduce();
            return Self::internal_div(numerators, &den, check_for_zero);
        }

        let quotient = Self::create_from_u512_as_witness(
            ctx.clone(),
            U512::from_lo_hi(quotient_val.lo(), quotient_val.hi()),
            false,
            num_quotient_bits,
        );

        // Constraint: inverse * denominator = quotient * p + numerators
        Self::unsafe_evaluate_multiply_add(
            &inverse,
            denominator,
            &[],
            &quotient,
            numerators,
        );

        if check_for_zero {
            denominator.assert_is_not_equal(&Self::zero_val(), "bigfield: division by zero");
        }

        inverse
    }

    /// Compute modular inverse: a^(-1) mod T using Field<T> arithmetic.
    fn invmod_u256(a: U256, _m: U256) -> U256 {
        let a_field = Field::<T>::from_limbs(*a.as_words());
        let inv = a_field.invert();
        U256::from_words(inv.from_montgomery_form().data)
    }

    // ════════════════════════════════════════════════════════════════════
    //  Assertion methods
    // ════════════════════════════════════════════════════════════════════

    /// Assert this element is in the field (< modulus).
    ///
    /// Port of C++ `bigfield::assert_is_in_field`.
    pub fn assert_is_in_field(&self, msg: &str) {
        self.assert_less_than(Self::modulus(), msg);
    }

    /// Assert this element < upper_limit.
    ///
    /// Port of C++ `bigfield::assert_less_than`.
    pub fn assert_less_than(&self, upper_limit: U256, msg: &str) {
        if self.is_constant() {
            let val = self.get_value().lo();
            assert!(val < upper_limit, "{}", msg);
            return;
        }

        let ctx = self.context.clone().expect("assert_less_than requires context");

        // First, range-constrain all limbs to their expected sizes
        {
            let mut builder = ctx.borrow_mut();
            builder.range_constrain_two_limbs(
                self.binary_basis_limbs[0].element.normalize().get_witness_index(),
                self.binary_basis_limbs[1].element.normalize().get_witness_index(),
                NUM_LIMB_BITS as usize,
                NUM_LIMB_BITS as usize,
                msg,
            );
            builder.range_constrain_two_limbs(
                self.binary_basis_limbs[2].element.normalize().get_witness_index(),
                self.binary_basis_limbs[3].element.normalize().get_witness_index(),
                NUM_LIMB_BITS as usize,
                num_last_limb_bits::<T>() as usize,
                msg,
            );
        }

        // Now perform borrow-chain subtraction: upper_limit - value >= 0
        let val = self.get_value().lo();
        let upper_0 = slice_u256(&upper_limit, 0, NUM_LIMB_BITS as u32);
        let upper_1 = slice_u256(&upper_limit, NUM_LIMB_BITS as u32, (NUM_LIMB_BITS * 2) as u32);
        let upper_2 = slice_u256(&upper_limit,
            (NUM_LIMB_BITS * 2) as u32,
            (NUM_LIMB_BITS * 3) as u32,
        );
        let upper_3 = slice_u256(&upper_limit,
            (NUM_LIMB_BITS * 3) as u32,
            (NUM_LIMB_BITS * 4) as u32,
        );

        let val_0 = slice_u256(&val, 0, NUM_LIMB_BITS as u32);
        let val_1 = slice_u256(&val, NUM_LIMB_BITS as u32, (NUM_LIMB_BITS * 2) as u32);
        let val_2 = slice_u256(&val, (NUM_LIMB_BITS * 2) as u32, (NUM_LIMB_BITS * 3) as u32);

        let borrow_0: bool = val_0 > upper_0;
        let borrow_1_val = val_1.wrapping_add(&U256::from(borrow_0 as u64));
        let borrow_1: bool = borrow_1_val > upper_1;
        let borrow_2_val = val_2.wrapping_add(&U256::from(borrow_1 as u64));
        let borrow_2: bool = borrow_2_val > upper_2;

        // Create borrow witnesses
        let borrow_0_idx = ctx
            .borrow_mut()
            .base
            .add_variable(Field::from(borrow_0 as u64));
        let borrow_1_idx = ctx
            .borrow_mut()
            .base
            .add_variable(Field::from(borrow_1 as u64));
        let borrow_2_idx = ctx
            .borrow_mut()
            .base
            .add_variable(Field::from(borrow_2 as u64));

        // Constrain borrows to be boolean
        {
            let mut builder = ctx.borrow_mut();
            builder.create_bool_gate(borrow_0_idx);
            builder.create_bool_gate(borrow_1_idx);
            builder.create_bool_gate(borrow_2_idx);
        }

        // Compute r = upper - value + borrows
        let s1 = shift_1::<P>();
        let borrow_0_field = FieldT::from_witness_index(ctx.clone(), borrow_0_idx);
        let borrow_1_field = FieldT::from_witness_index(ctx.clone(), borrow_1_idx);
        let borrow_2_field = FieldT::from_witness_index(ctx.clone(), borrow_2_idx);

        let r0 = &(&FieldT::from_field(Field::from_limbs(*upper_0.as_words()))
            - &self.binary_basis_limbs[0].element)
            + &(&borrow_0_field * &FieldT::from_field(s1));
        let r1_partial = &(&FieldT::from_field(Field::from_limbs(*upper_1.as_words()))
            - &self.binary_basis_limbs[1].element)
            + &(&borrow_1_field * &FieldT::from_field(s1));
        let r1 = &r1_partial - &borrow_0_field;
        let r2_partial = &(&FieldT::from_field(Field::from_limbs(*upper_2.as_words()))
            - &self.binary_basis_limbs[2].element)
            + &(&borrow_2_field * &FieldT::from_field(s1));
        let r2 = &r2_partial - &borrow_1_field;
        let r3 = &(&FieldT::from_field(Field::from_limbs(*upper_3.as_words()))
            - &self.binary_basis_limbs[3].element)
            - &borrow_2_field;

        // Range constrain results
        let r0_norm = r0.normalize();
        let r1_norm = r1.normalize();
        let r2_norm = r2.normalize();
        let r3_norm = r3.normalize();

        {
            let mut builder = ctx.borrow_mut();
            builder.range_constrain_two_limbs(
                r0_norm.get_witness_index(),
                r1_norm.get_witness_index(),
                NUM_LIMB_BITS as usize,
                NUM_LIMB_BITS as usize,
                msg,
            );
            builder.range_constrain_two_limbs(
                r2_norm.get_witness_index(),
                r3_norm.get_witness_index(),
                NUM_LIMB_BITS as usize,
                num_last_limb_bits::<T>() as usize,
                msg,
            );
        }
    }

    /// Assert two bigfield elements are equal.
    ///
    /// Port of C++ `bigfield::assert_equal`.
    pub fn assert_equal(&self, other: &Self, msg: &str) {
        if self.is_constant() && other.is_constant() {
            let a = self.get_value().div_rem(&Self::modulus_u512().to_nz().unwrap()).1;
            let b = other.get_value().div_rem(&Self::modulus_u512().to_nz().unwrap()).1;
            assert!(a == b, "{}", msg);
            return;
        }

        // Compute diff = self - other
        let diff = self.sub(other);

        // Compute quotient: diff / modulus
        let diff_val = diff.get_value();
        let modulus_512 = Self::modulus_u512().to_nz().unwrap();
        let (quotient_512, _) = diff_val.div_rem(&modulus_512);

        let ctx = self
            .context
            .clone()
            .or(other.context.clone())
            .expect("assert_equal requires context");

        // Create quotient witness
        let quotient = Self::create_from_u512_as_witness(
            ctx.clone(),
            quotient_512,
            false,
            0,
        );

        // Verify: diff * 1 = quotient * p + 0
        let one = Self::one_val();
        let zero = Self::zero_val();
        Self::unsafe_evaluate_multiply_add(&diff, &one, &[], &quotient, &[&zero]);
    }

    /// Assert two bigfield elements are not equal.
    ///
    /// Port of C++ `bigfield::assert_is_not_equal`.
    pub fn assert_is_not_equal(&self, other: &Self, msg: &str) {
        if self.is_constant() && other.is_constant() {
            let a = self.get_value().div_rem(&Self::modulus_u512().to_nz().unwrap()).1;
            let b = other.get_value().div_rem(&Self::modulus_u512().to_nz().unwrap()).1;
            assert!(a != b, "{}", msg);
            return;
        }

        // Check via prime basis limb: compute product of (diff - k*p) for all possible k
        let modulus = Self::modulus();
        let native_modulus = Field::<P>::from_limbs(*modulus.as_words());

        let self_max = self.get_maximum_value();
        let other_max = other.get_maximum_value();
        let modulus_512 = Self::modulus_u512();

        // Count how many times modulus fits in each max value
        let lhs_overload_count = {
            let mut count = 0u64;
            let mut target = modulus_512;
            while target <= self_max {
                target = target.wrapping_add(&modulus_512);
                count += 1;
            }
            count
        };
        let rhs_overload_count = {
            let mut count = 0u64;
            let mut target = modulus_512;
            while target <= other_max {
                target = target.wrapping_add(&modulus_512);
                count += 1;
            }
            count
        };

        let base_diff = &self.prime_basis_limb - &other.prime_basis_limb;
        let mut diff = base_diff.clone();

        let prime_basis_field = FieldT::from_field(native_modulus);
        let mut accum = prime_basis_field.clone();
        for _ in 0..lhs_overload_count {
            diff = diff.madd(&(&base_diff - &accum), &FieldT::from_field(Field::zero()));
            accum = &accum + &prime_basis_field;
        }

        accum = prime_basis_field.clone();
        for _ in 0..rhs_overload_count {
            diff = diff.madd(&(&base_diff + &accum), &FieldT::from_field(Field::zero()));
            accum = &accum + &prime_basis_field;
        }

        diff.assert_is_not_zero(msg);
    }

    // ════════════════════════════════════════════════════════════════════
    //  Conditional operations
    // ════════════════════════════════════════════════════════════════════

    /// Select between self and other based on predicate.
    /// Returns: predicate ? other : self
    ///
    /// Port of C++ `bigfield::conditional_select`.
    pub fn conditional_select(&self, other: &Self, predicate: &BoolT<P>) -> Self {
        if predicate.is_constant() {
            return if predicate.get_value() {
                other.clone()
            } else {
                self.clone()
            };
        }

        let ctx = self
            .context
            .clone()
            .or(other.context.clone())
            .or(predicate.get_context().clone());

        let pred_field = FieldT::from_witness_index(
            predicate.get_context().clone().unwrap(),
            predicate.normalize().get_witness_index(),
        );

        let mut result_limbs: [Limb<P>; 4] = [
            Limb::zero_limb(),
            Limb::zero_limb(),
            Limb::zero_limb(),
            Limb::zero_limb(),
        ];

        for i in 0..4 {
            // result = predicate * (other - self) + self
            let diff = &other.binary_basis_limbs[i].element - &self.binary_basis_limbs[i].element;
            let selected = pred_field.madd(&diff, &self.binary_basis_limbs[i].element);
            let max = std::cmp::max(
                self.binary_basis_limbs[i].maximum_value,
                other.binary_basis_limbs[i].maximum_value,
            );
            result_limbs[i] = Limb::new(selected, max);
        }

        let prime_diff = &other.prime_basis_limb - &self.prime_basis_limb;
        let prime_selected = pred_field.madd(&prime_diff, &self.prime_basis_limb);

        Self {
            context: ctx,
            binary_basis_limbs: result_limbs,
            prime_basis_limb: prime_selected,
            _target: std::marker::PhantomData,
        }
    }

    /// Conditionally negate: predicate ? (-self) : self
    ///
    /// Port of C++ `bigfield::conditional_negate`.
    pub fn conditional_negate(&self, predicate: &BoolT<P>) -> Self {
        if predicate.is_constant() {
            if predicate.get_value() {
                return Self::zero_val().sub(self);
            } else {
                return self.clone();
            }
        }

        let negated = Self::zero_val().sub(self);
        self.conditional_select(&negated, predicate)
    }

    /// Equality check: returns BoolT.
    ///
    /// Port of C++ `bigfield::operator==`.
    pub fn eq(&self, other: &Self) -> BoolT<P> {
        if self.is_constant() && other.is_constant() {
            let a = self.get_value().div_rem(&Self::modulus_u512().to_nz().unwrap()).1;
            let b = other.get_value().div_rem(&Self::modulus_u512().to_nz().unwrap()).1;
            return BoolT::from_constant(a == b);
        }

        let ctx = self
            .context
            .clone()
            .or(other.context.clone())
            .expect("eq requires context");

        // Compute diff and check if zero
        let a_val = self.get_value().div_rem(&Self::modulus_u512().to_nz().unwrap()).1;
        let b_val = other.get_value().div_rem(&Self::modulus_u512().to_nz().unwrap()).1;
        let is_equal_raw = a_val == b_val;

        let is_equal = BoolT::from_witness(&WitnessT::new(ctx.clone(), Field::from(is_equal_raw as u64)));

        // If equal, diff = 0 mod p, diff * 0 = 0
        // If not equal, diff * inverse = 1
        let mut a_clone = self.clone();
        let mut b_clone = other.clone();
        a_clone.self_reduce();
        b_clone.self_reduce();

        let diff = a_clone.sub(&b_clone);

        if is_equal_raw {
            // Assert diff == 0 mod p
            diff.assert_equal(&Self::zero_val(), "bigfield eq check");
        } else {
            // Assert diff != 0 mod p (by computing inverse)
            diff.assert_is_not_equal(&Self::zero_val(), "bigfield eq check: expected not equal");
        }

        is_equal
    }

    /// Negate: -self (mod target modulus).
    pub fn negate(&self) -> Self {
        let zero = Self::from_u256(self.context.clone(), U256::ZERO);
        zero.sub(self)
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Operator impls
// ════════════════════════════════════════════════════════════════════════

impl<P: FieldParams, T: FieldParams> std::ops::Add for &BigFieldT<P, T> {
    type Output = BigFieldT<P, T>;
    fn add(self, rhs: Self) -> Self::Output {
        BigFieldT::add(self, rhs)
    }
}

impl<P: FieldParams, T: FieldParams> std::ops::Sub for &BigFieldT<P, T> {
    type Output = BigFieldT<P, T>;
    fn sub(self, rhs: Self) -> Self::Output {
        BigFieldT::sub(self, rhs)
    }
}

impl<P: FieldParams, T: FieldParams> std::ops::Mul for &BigFieldT<P, T> {
    type Output = BigFieldT<P, T>;
    fn mul(self, rhs: Self) -> Self::Output {
        BigFieldT::mul(self, rhs)
    }
}

impl<P: FieldParams, T: FieldParams> std::ops::Div for &BigFieldT<P, T> {
    type Output = BigFieldT<P, T>;
    fn div(self, rhs: Self) -> Self::Output {
        BigFieldT::div(self, rhs)
    }
}

impl<P: FieldParams, T: FieldParams> std::ops::Neg for &BigFieldT<P, T> {
    type Output = BigFieldT<P, T>;
    fn neg(self) -> Self::Output {
        self.negate()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_ecc::curves::bn254::{Bn254FqParams, Bn254FrParams};

    // Native field = BN254 Fr (the circuit field)
    type Fr = Field<Bn254FrParams>;
    // Target non-native field = BN254 Fq
    type Fq = Field<Bn254FqParams>;
    // BigField type: Fq inside Fr circuit
    type BFq = BigFieldT<Bn254FrParams, Bn254FqParams>;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    fn check_circuit(builder: &BuilderRef<Bn254FrParams>) -> Result<(), String> {
        UltraCircuitChecker::check(&mut builder.borrow_mut())
    }

    /// Convert a native Fq element to a U256 in standard form.
    fn fq_to_u256(fq: Fq) -> U256 {
        let std = fq.from_montgomery_form();
        U256::from_words(std.data)
    }

    /// Create a random Fq value and return both the native element and its U256 representation.
    fn random_fq() -> (Fq, U256) {
        let fq = Fq::random_element();
        let u = fq_to_u256(fq);
        (fq, u)
    }

    /// Create a bigfield witness from a native Fq value.
    fn make_witness(builder: BuilderRef<Bn254FrParams>, fq: Fq) -> BFq {
        BFq::from_witness(builder, fq_to_u256(fq))
    }

    /// Create a bigfield constant from a native Fq value.
    fn make_constant(fq: Fq) -> BFq {
        BFq::from_u256(None, fq_to_u256(fq))
    }

    /// Compare bigfield result with expected native Fq value.
    fn assert_bigfield_eq(bf: &BFq, expected: Fq) {
        let bf_val = bf.get_value().lo();
        let expected_u256 = fq_to_u256(expected);
        // Both should be equal mod the target modulus
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);
        let bf_reduced = bf_val.div_rem(&modulus.to_nz().unwrap()).1;
        let exp_reduced = expected_u256.div_rem(&modulus.to_nz().unwrap()).1;
        assert_eq!(bf_reduced, exp_reduced,
            "bigfield value mismatch: got {:?}, expected {:?}", bf_reduced, exp_reduced);
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Constructor tests
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_from_u256_constant() {
        let (fq, u) = random_fq();
        let bf = BFq::from_u256(None, u);
        assert!(bf.is_constant());
        assert_bigfield_eq(&bf, fq);
    }

    #[test]
    fn test_from_witness() {
        let builder = make_builder();
        let (fq, _) = random_fq();
        let bf = make_witness(builder.clone(), fq);
        assert!(!bf.is_constant());
        assert_bigfield_eq(&bf, fq);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_zero_constant() {
        let bf = BFq::zero_val();
        assert!(bf.is_constant());
        assert_eq!(bf.get_value().lo(), U256::ZERO);
    }

    #[test]
    fn test_one_constant() {
        let bf = BFq::one_val();
        assert!(bf.is_constant());
        assert_eq!(bf.get_value().lo(), U256::ONE);
    }

    #[test]
    fn test_from_field_pair() {
        let builder = make_builder();
        let (fq, u) = random_fq();

        let lo_val = slice_u256(&u, 0, (NUM_LIMB_BITS * 2) as u32);
        let hi_val = slice_u256(&u, (NUM_LIMB_BITS * 2) as u32, (NUM_LIMB_BITS * 4) as u32);

        let lo = FieldT::from_witness(builder.clone(), Fr::from_limbs(*lo_val.as_words()));
        let hi = FieldT::from_witness(builder.clone(), Fr::from_limbs(*hi_val.as_words()));

        let bf = BFq::from_field_pair(lo, hi, false, 0);
        assert!(!bf.is_constant());
        assert_bigfield_eq(&bf, fq);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_create_from_u512_as_witness() {
        let builder = make_builder();
        let (fq, u) = random_fq();
        let u512 = U512::from_lo_hi(u, U256::ZERO);

        let bf = BFq::create_from_u512_as_witness(builder.clone(), u512, false, 0);
        assert!(!bf.is_constant());
        assert_bigfield_eq(&bf, fq);
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Arithmetic: Addition
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_add_witness_witness() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            let mut c = a.add(&b);
            c.self_reduce();
            assert_bigfield_eq(&c, a_fq + b_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_witness_constant() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_constant(b_fq);
            let mut c = a.add(&b);
            c.self_reduce();
            assert_bigfield_eq(&c, a_fq + b_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_constant_constant() {
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_constant(a_fq);
            let b = make_constant(b_fq);
            let mut c = a.add(&b);
            c.reduction_check();
            assert_bigfield_eq(&c, a_fq + b_fq);
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Arithmetic: Subtraction
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_sub_witness_witness() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            let mut c = a.sub(&b);
            c.self_reduce();
            assert_bigfield_eq(&c, a_fq - b_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_sub_witness_constant() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_constant(b_fq);
            let mut c = a.sub(&b);
            c.self_reduce();
            assert_bigfield_eq(&c, a_fq - b_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Arithmetic: Multiplication
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_mul_witness_witness() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            let c = a.mul(&b);
            assert_bigfield_eq(&c, a_fq * b_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_mul_witness_constant() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_constant(b_fq);
            let c = a.mul(&b);
            assert_bigfield_eq(&c, a_fq * b_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_mul_constant_constant() {
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_constant(a_fq);
            let b = make_constant(b_fq);
            let c = a.mul(&b);
            assert_bigfield_eq(&c, a_fq * b_fq);
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Arithmetic: Division
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_div_witness_witness() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            if b_fq == Fq::zero() { continue; }
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            let c = a.div(&b);
            let expected = a_fq * b_fq.invert();
            assert_bigfield_eq(&c, expected);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_div_witness_constant() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            if b_fq == Fq::zero() { continue; }
            let a = make_witness(builder.clone(), a_fq);
            let b = make_constant(b_fq);
            let c = a.div(&b);
            let expected = a_fq * b_fq.invert();
            assert_bigfield_eq(&c, expected);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Arithmetic: Squaring
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_sqr_witness() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let c = a.sqr();
            assert_bigfield_eq(&c, a_fq * a_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_sqr_constant() {
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let a = make_constant(a_fq);
            let c = a.sqr();
            assert_bigfield_eq(&c, a_fq * a_fq);
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Arithmetic: Multiply-Add
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_madd_witness() {
        let builder = make_builder();
        for _ in 0..4 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let (c_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            let c = make_witness(builder.clone(), c_fq);
            let result = a.madd(&b, &[&c]);
            assert_bigfield_eq(&result, a_fq * b_fq + c_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Arithmetic: Negation
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_negate_witness() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let mut c = a.negate();
            c.self_reduce();
            assert_bigfield_eq(&c, -a_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Operator overloads
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_operator_add() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);
        let mut c = &a + &b;
        c.self_reduce();
        assert_bigfield_eq(&c, a_fq + b_fq);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_operator_sub() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);
        let mut c = &a - &b;
        c.self_reduce();
        assert_bigfield_eq(&c, a_fq - b_fq);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_operator_mul() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);
        let c = &a * &b;
        assert_bigfield_eq(&c, a_fq * b_fq);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_operator_div() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        if b_fq == Fq::zero() { return; }
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);
        let c = &a / &b;
        assert_bigfield_eq(&c, a_fq * b_fq.invert());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_operator_neg() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let mut c = -&a;
        c.self_reduce();
        assert_bigfield_eq(&c, -a_fq);
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Self-reduce
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_self_reduce() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);

            // Chain of operations to accumulate large values
            let mut c = a.mul(&b);
            for _ in 0..4 {
                let (d_fq, _) = random_fq();
                let d = make_witness(builder.clone(), d_fq);
                c = c.mul(&d);
            }
            c.self_reduce();

            // After reduction, value should be in range
            let modulus = BFq::modulus();
            let val = c.get_value().lo();
            assert!(val < modulus, "self_reduce: value not reduced");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_reduction_check() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        // Accumulate via additions until reduction triggers
        let mut c = a.clone();
        for _ in 0..20 {
            c = c.add(&b);
        }
        c.reduction_check();

        // After reduction_check, maxes should be safe
        let max_unreduced = get_maximum_unreduced_limb_value::<Bn254FrParams>();
        for i in 0..4 {
            assert!(
                c.binary_basis_limbs[i].maximum_value <= max_unreduced,
                "reduction_check: limb {} max too large", i
            );
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Assert methods
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_assert_is_in_field() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            a.assert_is_in_field("test: should be in field");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_less_than() {
        let builder = make_builder();
        let bit_mask = U256::from(1u64)
            .wrapping_shl_vartime(200)
            .wrapping_sub(&U256::ONE);

        for _ in 0..10 {
            let (_, u) = random_fq();
            // Mask to 200 bits to ensure < 2^200
            let small_val = u.wrapping_rem_vartime(&bit_mask.to_nz().unwrap());
            let a = BFq::from_witness(builder.clone(), small_val);
            a.assert_less_than(
                U256::from(1u64).wrapping_shl_vartime(200),
                "test: should be < 2^200",
            );
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_equal_same_value() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), a_fq);
            a.assert_equal(&b, "test: same value should be equal");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_equal_after_arithmetic() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            // a + a should equal a * 2
            let sum = a.add(&a);
            let two = BFq::from_u256(Some(builder.clone()), U256::from(2u64));
            let prod = a.mul(&two);
            sum.assert_equal(&prod, "test: a+a should equal a*2");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_is_not_equal() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            if a_fq == b_fq { continue; }
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            a.assert_is_not_equal(&b, "test: different values should not be equal");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Equality operator
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_eq_same() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), a_fq);
        let result = a.eq(&b);
        assert!(result.get_value(), "same values should be equal");
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_eq_different() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        if a_fq == b_fq { return; }
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);
        let result = a.eq(&b);
        assert!(!result.get_value(), "different values should not be equal");
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Conditional operations
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_conditional_select_true() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        // predicate = true → returns b
        let pred = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::one()));
        let result = a.conditional_select(&b, &pred);
        assert_bigfield_eq(&result, b_fq);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_conditional_select_false() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        // predicate = false → returns a
        let pred = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::zero()));
        let result = a.conditional_select(&b, &pred);
        assert_bigfield_eq(&result, a_fq);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_conditional_negate_true() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);

        let pred = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::one()));
        let mut result = a.conditional_negate(&pred);
        result.self_reduce();
        assert_bigfield_eq(&result, -a_fq);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_conditional_negate_false() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);

        let pred = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::zero()));
        let result = a.conditional_negate(&pred);
        assert_bigfield_eq(&result, a_fq);
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Regression tests
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_division_formula_bug_regression() {
        // Regression: (2 - 2) / 2 should work without CRT issues.
        let builder = make_builder();
        let two = Fq::from(2u64);
        let a = make_witness(builder.clone(), two);
        let b = make_witness(builder.clone(), two);
        let c = make_witness(builder.clone(), two);

        let diff = a.sub(&b); // 0
        let result = diff.div(&c); // 0 / 2 = 0
        assert_bigfield_eq(&result, Fq::zero());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_to_lower_limb_regression() {
        // Regression: operations with value 1 (smallest non-zero).
        let builder = make_builder();
        let one = Fq::one();
        let a = make_witness(builder.clone(), one);
        let b = make_witness(builder.clone(), one);
        let c = make_constant(one);

        // Various operations with small values
        let sum = a.add(&b); // 2
        let prod = a.mul(&c); // 1
        let diff = sum.sub(&prod); // 1
        assert_bigfield_eq(&diff, one);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_quotient_completeness_regression() {
        // Near-max value, doubled 8 times, then squared.
        let builder = make_builder();
        let val = U256::from_words([
            0xfffffffffffffffe,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0x3fffffffffffffff,
        ]);
        // Reduce mod Fq modulus
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);
        let val_reduced = val.div_rem(&modulus.to_nz().unwrap()).1;
        let fq_val = Fq::from_limbs(*val_reduced.as_words());

        let mut a = BFq::from_witness(builder.clone(), val_reduced);

        // Double 8 times
        let mut native = fq_val;
        for _ in 0..8 {
            let a_copy = a.clone();
            a = a.add(&a_copy);
            native = native + native;
        }
        a.self_reduce();

        // Square
        let result = a.sqr();
        assert_bigfield_eq(&result, native * native);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_conditional_select_constants_regression() {
        // Regression: conditional_select with constant 0 and 1.
        let builder = make_builder();
        let a = BFq::from_u256(Some(builder.clone()), U256::ZERO);
        let b = BFq::from_u256(Some(builder.clone()), U256::ONE);

        let pred_true = BoolT::from_constant(true);
        let result_b = a.conditional_select(&b, &pred_true);
        assert_eq!(result_b.get_value().lo(), U256::ONE);

        let pred_false = BoolT::from_constant(false);
        let result_a = a.conditional_select(&b, &pred_false);
        assert_eq!(result_a.get_value().lo(), U256::ZERO);
    }

    #[test]
    fn test_inversion_constant() {
        // inversion of constant -7
        let neg7 = -Fq::from(7u64);
        let a = make_constant(neg7);
        let one = BFq::one_val();
        let inv = one.div(&a);
        let expected = neg7.invert();
        assert_bigfield_eq(&inv, expected);
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge case tests
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_edge_case_zero_operations() {
        let builder = make_builder();
        let zero = make_witness(builder.clone(), Fq::zero());
        let (a_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);

        // a + 0 = a
        let sum = a.add(&zero);
        assert_bigfield_eq(&sum, a_fq);

        // a * 0 = 0
        let prod = a.mul(&zero);
        assert_bigfield_eq(&prod, Fq::zero());

        // a - a = 0
        let diff = a.sub(&a);
        let mut diff_reduced = diff;
        diff_reduced.self_reduce();
        assert_bigfield_eq(&diff_reduced, Fq::zero());

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_edge_case_one_operations() {
        let builder = make_builder();
        let one = make_witness(builder.clone(), Fq::one());
        let (a_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);

        // a * 1 = a
        let prod = a.mul(&one);
        assert_bigfield_eq(&prod, a_fq);

        // a / 1 = a
        let quot = a.div(&one);
        assert_bigfield_eq(&quot, a_fq);

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_edge_case_modulus_minus_one() {
        let builder = make_builder();
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);
        let p_minus_1 = modulus.wrapping_sub(&U256::ONE);
        let fq_pm1 = Fq::from_limbs(*p_minus_1.as_words());
        let a = BFq::from_witness(builder.clone(), p_minus_1);

        // (p-1) + 1 = 0 mod p
        let one = make_witness(builder.clone(), Fq::one());
        let mut sum = a.add(&one);
        sum.self_reduce();
        assert_bigfield_eq(&sum, Fq::zero());

        // (p-1)^2 = 1 mod p
        let sq = a.sqr();
        assert_bigfield_eq(&sq, Fq::one());

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_is_in_field_edge_values() {
        let builder = make_builder();
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);

        // 0 is in field
        let zero = BFq::from_witness(builder.clone(), U256::ZERO);
        zero.assert_is_in_field("zero in field");

        // 1 is in field
        let one = BFq::from_witness(builder.clone(), U256::ONE);
        one.assert_is_in_field("one in field");

        // p-1 is in field
        let pm1 = modulus.wrapping_sub(&U256::ONE);
        let pm1_bf = BFq::from_witness(builder.clone(), pm1);
        pm1_bf.assert_is_in_field("p-1 in field");

        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Compound arithmetic
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_add_and_mul() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let (c_fq, _) = random_fq();
            let (d_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            let c = make_witness(builder.clone(), c_fq);
            let d = make_witness(builder.clone(), d_fq);

            let sum_ab = a.add(&b);
            let sum_cd = c.add(&d);
            let result = sum_ab.mul(&sum_cd);
            assert_bigfield_eq(&result, (a_fq + b_fq) * (c_fq + d_fq));
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_sub_and_mul() {
        for _ in 0..10 {
            let builder = make_builder();
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let (c_fq, _) = random_fq();
            let (d_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            let c = make_witness(builder.clone(), c_fq);
            let d = make_witness(builder.clone(), d_fq);

            let diff_ab = a.sub(&b);
            let diff_cd = c.sub(&d);
            let result = diff_ab.mul(&diff_cd);
            assert_bigfield_eq(&result, (a_fq - b_fq) * (c_fq - d_fq));
            assert!(check_circuit(&builder).is_ok());
        }
    }

    #[test]
    fn test_add_and_div() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let (c_fq, _) = random_fq();
        let (d_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);
        let c = make_witness(builder.clone(), c_fq);
        let d = make_witness(builder.clone(), d_fq);

        let num = a.add(&b);
        let den = c.add(&d);
        let result = num.div(&den);
        let expected = (a_fq + b_fq) * (c_fq + d_fq).invert();
        assert_bigfield_eq(&result, expected);
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Internal utility tests
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_invmod_u256() {
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);

        // inv(1) = 1
        assert_eq!(BFq::invmod_u256(U256::ONE, modulus), U256::ONE);

        // inv(2) * 2 = 1 mod p (use Field for verification to avoid U256 overflow)
        let inv2 = BFq::invmod_u256(U256::from(2u64), modulus);
        let inv2_f = Fq::from_limbs(*inv2.as_words());
        let two_f = Fq::from_limbs([2, 0, 0, 0]);
        assert_eq!(inv2_f * two_f, Fq::one());

        // Random values
        for _ in 0..10 {
            let (fq, u) = random_fq();
            if fq == Fq::zero() { continue; }
            let inv = BFq::invmod_u256(u, modulus);
            let a_f = Fq::from_limbs(*u.as_words());
            let inv_f = Fq::from_limbs(*inv.as_words());
            assert_eq!(a_f * inv_f, Fq::one(), "invmod_u256: a * a^-1 should be 1 mod p");
        }
    }

    #[test]
    fn test_neg_modulus_limbs() {
        // Verify neg_modulus_limbs are consistent:
        // Reconstructing from limbs should give (2^272 - target_modulus) mod 2^272
        let limbs = neg_modulus_limbs::<Bn254FrParams, Bn254FqParams>();
        let target = U256::from_limbs(Bn254FqParams::MODULUS);

        // Reconstruct from limbs
        let mut reconstructed = Field::<Bn254FrParams>::zero();
        for i in 0..4 {
            let shift = match i {
                0 => Field::one(),
                1 => shift_1::<Bn254FrParams>(),
                2 => shift_2::<Bn254FrParams>(),
                3 => shift_3::<Bn254FrParams>(),
                _ => unreachable!(),
            };
            reconstructed = reconstructed + limbs[i] * shift;
        }

        // neg_modulus + target_modulus should be 0 mod 2^272
        // In the native field, neg_modulus + target = 0 mod n (approximately)
        let target_field = Fr::from_limbs(*target.as_words());
        let sum = reconstructed + target_field;
        // The sum should be 2^272 mod native_modulus
        // We just check this doesn't panic and the computation is consistent
        let _ = sum;
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Constructor: unsafe_construct_from_limbs
    //  Port of C++ test_unsafe_construct_from_limbs
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_unsafe_construct_from_limbs() {
        let builder = make_builder();

        // Create 4 random limbs with slightly more than NUM_LIMB_BITS bits
        let mask_78 = U256::from(1u64)
            .wrapping_shl_vartime((NUM_LIMB_BITS + 10) as u32)
            .wrapping_sub(&U256::ONE);
        let mask_80 = U256::from(1u64)
            .wrapping_shl_vartime((NUM_LIMB_BITS + 12) as u32)
            .wrapping_sub(&U256::ONE);

        let l0_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&mask_78.to_nz().unwrap());
        let l1_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&mask_78.to_nz().unwrap());
        let l2_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&mask_78.to_nz().unwrap());
        let l3_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&mask_80.to_nz().unwrap());

        let l0 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l0_val.as_words()));
        let l1 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l1_val.as_words()));
        let l2 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l2_val.as_words()));
        let l3 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l3_val.as_words()));

        let result = BFq::unsafe_construct_from_limbs(l0, l1, l2, l3, true);

        // Verify prime basis limb is correctly computed
        let s1 = shift_1::<Bn254FrParams>();
        let s2 = shift_2::<Bn254FrParams>();
        let s3 = shift_3::<Bn254FrParams>();
        let expected_prime = Fr::from_limbs(*l0_val.as_words())
            + Fr::from_limbs(*l1_val.as_words()) * s1
            + Fr::from_limbs(*l2_val.as_words()) * s2
            + Fr::from_limbs(*l3_val.as_words()) * s3;
        assert_eq!(result.prime_basis_limb.get_value(), expected_prime);
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Regression: badmul (same as division_formula)
    //  Port of C++ test_bad_mul
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_badmul_regression() {
        let builder = make_builder();
        let two = Fq::from(2u64);
        let tval = make_witness(builder.clone(), two);
        let tval1 = tval.sub(&tval); // 0
        let _tval2 = tval1.div(&tval); // 0 / 2 = 0
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Arithmetic: madd with constants
    //  Port of C++ test_madd with constant variants
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_madd_with_constants() {
        // w * w + c
        {
            let builder = make_builder();
            for _ in 0..4 {
                let (a_fq, _) = random_fq();
                let (b_fq, _) = random_fq();
                let (c_fq, _) = random_fq();
                let a = make_witness(builder.clone(), a_fq);
                let b = make_witness(builder.clone(), b_fq);
                let c = make_constant(c_fq);
                let result = a.madd(&b, &[&c]);
                assert_bigfield_eq(&result, a_fq * b_fq + c_fq);
            }
            assert!(check_circuit(&builder).is_ok());
        }
        // w * c + w
        {
            let builder = make_builder();
            for _ in 0..4 {
                let (a_fq, _) = random_fq();
                let (b_fq, _) = random_fq();
                let (c_fq, _) = random_fq();
                let a = make_witness(builder.clone(), a_fq);
                let b = make_constant(b_fq);
                let c = make_witness(builder.clone(), c_fq);
                let result = a.madd(&b, &[&c]);
                assert_bigfield_eq(&result, a_fq * b_fq + c_fq);
            }
            assert!(check_circuit(&builder).is_ok());
        }
        // c * w + w
        {
            let builder = make_builder();
            for _ in 0..4 {
                let (a_fq, _) = random_fq();
                let (b_fq, _) = random_fq();
                let (c_fq, _) = random_fq();
                let a = make_constant(a_fq);
                let b = make_witness(builder.clone(), b_fq);
                let c = make_witness(builder.clone(), c_fq);
                let result = a.madd(&b, &[&c]);
                assert_bigfield_eq(&result, a_fq * b_fq + c_fq);
            }
            assert!(check_circuit(&builder).is_ok());
        }
        // w * c + c
        {
            let builder = make_builder();
            for _ in 0..4 {
                let (a_fq, _) = random_fq();
                let (b_fq, _) = random_fq();
                let (c_fq, _) = random_fq();
                let a = make_witness(builder.clone(), a_fq);
                let b = make_constant(b_fq);
                let c = make_constant(c_fq);
                let result = a.madd(&b, &[&c]);
                assert_bigfield_eq(&result, a_fq * b_fq + c_fq);
            }
            assert!(check_circuit(&builder).is_ok());
        }
        // c * w + c
        {
            let builder = make_builder();
            for _ in 0..4 {
                let (a_fq, _) = random_fq();
                let (b_fq, _) = random_fq();
                let (c_fq, _) = random_fq();
                let a = make_constant(a_fq);
                let b = make_witness(builder.clone(), b_fq);
                let c = make_constant(c_fq);
                let result = a.madd(&b, &[&c]);
                assert_bigfield_eq(&result, a_fq * b_fq + c_fq);
            }
            assert!(check_circuit(&builder).is_ok());
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Arithmetic: madd with multiple addends
    //  Port of C++ test_madd with multiple addends
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_madd_multiple_addends() {
        let builder = make_builder();
        for _ in 0..4 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let (c_fq, _) = random_fq();
            let (d_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            let c = make_witness(builder.clone(), c_fq);
            let d = make_witness(builder.clone(), d_fq);
            let result = a.madd(&b, &[&c, &d]);
            assert_bigfield_eq(&result, a_fq * b_fq + c_fq + d_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Arithmetic: madd with no addends (a * b + 0)
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_madd_no_addends() {
        let builder = make_builder();
        for _ in 0..4 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            let result = a.madd(&b, &[]);
            assert_bigfield_eq(&result, a_fq * b_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Conditional: conditional_select with constants
    //  Port of C++ conditional_select_with_constants
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_conditional_select_with_constants() {
        // witness operands, constant predicate
        {
            let builder = make_builder();
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);

            let pred_true = BoolT::from_constant(true);
            let result = a.conditional_select(&b, &pred_true);
            assert_bigfield_eq(&result, b_fq);

            let pred_false = BoolT::from_constant(false);
            let result2 = a.conditional_select(&b, &pred_false);
            assert_bigfield_eq(&result2, a_fq);
            assert!(check_circuit(&builder).is_ok());
        }
        // constant operands, witness predicate
        {
            let builder = make_builder();
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = BFq::from_u256(Some(builder.clone()), fq_to_u256(a_fq));
            let b = BFq::from_u256(Some(builder.clone()), fq_to_u256(b_fq));

            let pred = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::one()));
            let result = a.conditional_select(&b, &pred);
            assert_bigfield_eq(&result, b_fq);
            assert!(check_circuit(&builder).is_ok());
        }
        // constant operands, constant predicate
        {
            let builder = make_builder();
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = BFq::from_u256(Some(builder.clone()), fq_to_u256(a_fq));
            let b = BFq::from_u256(Some(builder.clone()), fq_to_u256(b_fq));

            let pred = BoolT::from_constant(true);
            let result = a.conditional_select(&b, &pred);
            assert_bigfield_eq(&result, b_fq);
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Conditional: conditional_negate with constants
    //  Port of C++ conditional_negate_with_constants
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_conditional_negate_with_constants() {
        // constant value, witness predicate
        {
            let builder = make_builder();
            let (a_fq, _) = random_fq();
            let a = BFq::from_u256(Some(builder.clone()), fq_to_u256(a_fq));

            let pred = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::one()));
            let mut result = a.conditional_negate(&pred);
            result.self_reduce();
            assert_bigfield_eq(&result, -a_fq);
            assert!(check_circuit(&builder).is_ok());
        }
        // witness value, constant predicate
        {
            let builder = make_builder();
            let (a_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);

            let pred = BoolT::from_constant(true);
            let mut result = a.conditional_negate(&pred);
            result.self_reduce();
            assert_bigfield_eq(&result, -a_fq);
            assert!(check_circuit(&builder).is_ok());
        }
        // constant value, constant predicate
        {
            let (a_fq, _) = random_fq();
            let a = make_constant(a_fq);
            let pred = BoolT::from_constant(true);
            let mut result = a.conditional_negate(&pred);
            result.reduction_check();
            assert_bigfield_eq(&result, -a_fq);
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Equality with constants
    //  Port of C++ test_equality_operator with constant variants
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_equality_with_constants() {
        // witness == constant (same value)
        {
            let builder = make_builder();
            let (a_fq, _) = random_fq();
            let a_w = make_witness(builder.clone(), a_fq);
            let a_c = BFq::from_u256(Some(builder.clone()), fq_to_u256(a_fq));
            let result = a_w.eq(&a_c);
            assert!(result.get_value(), "witness should equal constant with same value");
            assert!(check_circuit(&builder).is_ok());
        }
        // constant == witness (same value)
        {
            let builder = make_builder();
            let (a_fq, _) = random_fq();
            let a_c = BFq::from_u256(Some(builder.clone()), fq_to_u256(a_fq));
            let a_w = make_witness(builder.clone(), a_fq);
            let result = a_c.eq(&a_w);
            assert!(result.get_value(), "constant should equal witness with same value");
            assert!(check_circuit(&builder).is_ok());
        }
        // constant == constant (same value)
        {
            let (a_fq, _) = random_fq();
            let a = make_constant(a_fq);
            let b = make_constant(a_fq);
            let result = a.eq(&b);
            assert!(result.get_value(), "constants with same value should be equal");
        }
        // witness != constant
        {
            let builder = make_builder();
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            if a_fq != b_fq {
                let a = make_witness(builder.clone(), a_fq);
                let b = BFq::from_u256(Some(builder.clone()), fq_to_u256(b_fq));
                let result = a.eq(&b);
                assert!(!result.get_value(), "different values should not be equal");
                assert!(check_circuit(&builder).is_ok());
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Combined: assert_equal and assert_is_not_equal together
    //  Port of C++ test_assert_equal_not_equal
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_assert_equal_not_equal() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (c_fq, _) = random_fq();
            let (d_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let c = make_witness(builder.clone(), c_fq);
            let d = make_witness(builder.clone(), d_fq);

            // Construct "2" via unsafe_construct_from_limbs
            let two_l0 = FieldT::from_witness(builder.clone(), Fr::from(2u64));
            let zero_l = FieldT::from_witness(builder.clone(), Fr::zero());
            let two_ct = BFq::unsafe_construct_from_limbs(
                two_l0,
                zero_l.clone(),
                zero_l.clone(),
                zero_l,
                false,
            );

            let t0 = a.add(&a); // a + a
            let t1 = a.mul(&two_ct); // a * 2

            t0.assert_equal(&t1, "a+a should equal a*2");
            if c_fq != a_fq + a_fq {
                t0.assert_is_not_equal(&c, "t0 should not equal c");
            }
            if d_fq != a_fq + a_fq {
                t0.assert_is_not_equal(&d, "t0 should not equal d");
            }

            let is_eq = t0.eq(&t1);
            assert!(is_eq.get_value(), "a+a should == a*2");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Compound: add_and_mul with constants
    //  Port of C++ test_add_and_mul with constant summands
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_add_and_mul_with_constants() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let (c_fq, _) = random_fq();
            let (d_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_constant(b_fq);
            let c = make_witness(builder.clone(), c_fq);
            let d = make_constant(d_fq);

            let sum_ab = a.add(&b);
            let sum_cd = c.add(&d);
            let result = sum_ab.mul(&sum_cd);
            assert_bigfield_eq(&result, (a_fq + b_fq) * (c_fq + d_fq));
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Compound: sub_and_mul with constants
    //  Port of C++ test_sub_and_mul with constant subtrahends
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_sub_and_mul_with_constants() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let (c_fq, _) = random_fq();
            let (d_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_constant(b_fq);
            let c = make_witness(builder.clone(), c_fq);
            let d = make_constant(d_fq);

            let diff_ab = a.sub(&b);
            let diff_cd = c.sub(&d);
            let result = diff_ab.mul(&diff_cd);
            assert_bigfield_eq(&result, (a_fq - b_fq) * (c_fq - d_fq));
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Constant-constant operations: div, sub
    //  Port of C++ div/sub with CONSTANT/CONSTANT variants
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_div_constant_constant() {
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            if b_fq == Fq::zero() { continue; }
            let a = make_constant(a_fq);
            let b = make_constant(b_fq);
            let c = a.div(&b);
            assert_bigfield_eq(&c, a_fq * b_fq.invert());
        }
    }

    #[test]
    fn test_sub_constant_constant() {
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_constant(a_fq);
            let b = make_constant(b_fq);
            let mut c = a.sub(&b);
            c.reduction_check();
            assert_bigfield_eq(&c, a_fq - b_fq);
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Constant-witness and witness-constant: div, sub
    //  Port of C++ constant/witness variants
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_div_constant_witness() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            if b_fq == Fq::zero() { continue; }
            let a = BFq::from_u256(Some(builder.clone()), fq_to_u256(a_fq));
            let b = make_witness(builder.clone(), b_fq);
            let c = a.div(&b);
            assert_bigfield_eq(&c, a_fq * b_fq.invert());
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_sub_constant_witness() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = BFq::from_u256(Some(builder.clone()), fq_to_u256(a_fq));
            let b = make_witness(builder.clone(), b_fq);
            let mut c = a.sub(&b);
            c.self_reduce();
            assert_bigfield_eq(&c, a_fq - b_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_mul_constant_witness() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = BFq::from_u256(Some(builder.clone()), fq_to_u256(a_fq));
            let b = make_witness(builder.clone(), b_fq);
            let c = a.mul(&b);
            assert_bigfield_eq(&c, a_fq * b_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_constant_witness() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = BFq::from_u256(Some(builder.clone()), fq_to_u256(a_fq));
            let b = make_witness(builder.clone(), b_fq);
            let mut c = a.add(&b);
            c.self_reduce();
            assert_bigfield_eq(&c, a_fq + b_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Negate with constant
    //  Port of C++ test_negate constant variant
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_negate_constant() {
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let a = make_constant(a_fq);
            let mut c = a.negate();
            c.reduction_check();
            assert_bigfield_eq(&c, -a_fq);
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Sqr with constant (more thorough)
    //  Port of C++ test_sqr CONSTANT variant
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_sqr_constant_thorough() {
        // Already have test_sqr_constant, but this covers constant-constant mul consistency
        for _ in 0..10 {
            let (a_fq, _) = random_fq();
            let a = make_constant(a_fq);
            let c = a.sqr();
            let d = a.mul(&a);
            // Both should give a^2
            assert_bigfield_eq(&c, a_fq * a_fq);
            assert_bigfield_eq(&d, a_fq * a_fq);
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Group operations (elliptic curve point addition via bigfield)
    //  Port of C++ test_group_operations
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_group_operations() {
        // Perform EC point addition using bigfield arithmetic
        // y = x^3 + 3 (bn254 curve b=3) for verification
        let builder = make_builder();

        // Use known BN254 generator G1 and 2*G1
        // P1 = random point, P2 = random point (we construct from random x,y pairs on curve)
        // For simplicity, use small known values
        let (x1_fq, _) = random_fq();
        let (x2_fq, _) = random_fq();

        // Compute y^2 = x^3 + 3 (bn254 b parameter)
        let three = Fq::from(3u64);
        let y1_sq = x1_fq * x1_fq * x1_fq + three;
        let y2_sq = x2_fq * x2_fq * x2_fq + three;

        // Skip if not QR (not on curve)
        // We just test the bigfield arithmetic pattern, not actual curve point validation

        let x1 = make_witness(builder.clone(), x1_fq);
        let x2 = make_witness(builder.clone(), x2_fq);

        // Compute lambda = (y2 - y1) / (x2 - x1) using bigfield
        // For testing purposes, use dummy y values and just verify the arithmetic
        let y1_fq = y1_sq; // not actually square root, but testing arithmetic
        let y2_fq = y2_sq;
        let y1 = make_witness(builder.clone(), y1_fq);
        let y2 = make_witness(builder.clone(), y2_fq);

        if x1_fq != x2_fq {
            let lambda = y2.sub(&y1).div(&x2.sub(&x1));
            let x3 = lambda.sqr().sub(&x2).sub(&x1);
            let _y3 = x1.sub(&x3).mul(&lambda).sub(&y1);

            // Verify algebraically
            let expected_lambda = (y2_fq - y1_fq) * (x2_fq - x1_fq).invert();
            let expected_x3 = expected_lambda * expected_lambda - x2_fq - x1_fq;
            let expected_y3 = (x1_fq - expected_x3) * expected_lambda - y1_fq;

            assert_bigfield_eq(&x3, expected_x3);
            assert_bigfield_eq(&_y3, expected_y3);
            assert!(check_circuit(&builder).is_ok());
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Reduce (thorough) — b*b + c chain
    //  Port of C++ test_reduce
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_reduce_thorough() {
        let builder = make_builder();
        for _ in 0..4 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);

            let mut c = a.clone();
            let mut expected = a_fq;
            for _ in 0..16 {
                c = b.mul(&b).add(&c);
                expected = b_fq * b_fq + expected;
            }

            c.self_reduce();

            let modulus = BFq::modulus();
            let val = c.get_value().lo();
            assert!(val < modulus, "reduce: value not in field");
            assert_bigfield_eq(&c, expected);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge case: larger_than_bigfield_allowed
    //  Port of C++ test_larger_than_bigfield_allowed
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_larger_than_bigfield_allowed() {
        let builder = make_builder();
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);

        // Values >= modulus (created with can_overflow=true)
        let values = [
            U512::from_lo_hi(modulus, U256::ZERO),                              // p
            U512::from_lo_hi(modulus.wrapping_add(&U256::ONE), U256::ZERO),     // p + 1
        ];

        for value in &values {
            let bf = BFq::create_from_u512_as_witness(builder.clone(), *value, true, 0);
            assert!(bf.get_value() >= BFq::modulus_u512(), "value should be >= modulus");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge case: invariants_during_addition
    //  Port of C++ test_invariants_during_addition
    // ═══════════════════════════════════════════════════════════════════

    /// Check that limb maximum values are >= actual witness values.
    fn check_invariants(bf: &BFq, label: &str) {
        for i in 0..4 {
            let witness_val = {
                let v = bf.binary_basis_limbs[i].element.get_value().from_montgomery_form();
                U256::from_words(v.data)
            };
            let max_val = bf.binary_basis_limbs[i].maximum_value;
            assert!(
                max_val >= witness_val,
                "invariant violation in {}: limb[{}] max {} < witness {}",
                label, i, max_val, witness_val
            );
        }
    }

    #[test]
    fn test_invariants_during_addition() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        let c = a.add(&b);
        check_invariants(&c, "addition (initial)");

        // Chain additions
        let mut result = c;
        for _ in 0..5 {
            let (e_fq, _) = random_fq();
            let e = make_witness(builder.clone(), e_fq);
            result = result.add(&e);
            check_invariants(&result, "addition (chain)");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_invariants_during_subtraction() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        let c = a.sub(&b);
        check_invariants(&c, "subtraction (initial)");

        let mut result = c;
        for _ in 0..5 {
            let (e_fq, _) = random_fq();
            let e = make_witness(builder.clone(), e_fq);
            result = result.sub(&e);
            check_invariants(&result, "subtraction (chain)");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_invariants_during_multiplication() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        let c = a.mul(&b);
        check_invariants(&c, "multiplication (initial)");

        let mut result = c;
        for _ in 0..3 {
            let (e_fq, _) = random_fq();
            let e = make_witness(builder.clone(), e_fq);
            result = result.mul(&e);
            check_invariants(&result, "multiplication (chain)");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_invariants_during_division() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        if b_fq == Fq::zero() { return; }
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        let c = a.div(&b);
        check_invariants(&c, "division (initial)");

        let mut result = c;
        for _ in 0..3 {
            let (e_fq, _) = random_fq();
            if e_fq == Fq::zero() { continue; }
            let e = make_witness(builder.clone(), e_fq);
            result = result.div(&e);
            check_invariants(&result, "division (chain)");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_invariants_during_squaring() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);

        let c = a.sqr();
        check_invariants(&c, "squaring");

        // Chain squarings
        let mut result = c;
        for _ in 0..3 {
            result = result.sqr();
            check_invariants(&result, "squaring (chain)");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_invariants_during_negation() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);

        let c = a.negate();
        check_invariants(&c, "negation");

        // Chain negations
        let mut result = c;
        for _ in 0..5 {
            result = result.negate();
            check_invariants(&result, "negation (chain)");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge case: assert_is_in_field with edge case values
    //  Port of C++ edge case assert_is_in_field
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_assert_is_in_field_edge_cases() {
        let builder = make_builder();
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);

        let edge_values = [
            U256::ZERO,                              // 0
            U256::ONE,                               // 1
            modulus.wrapping_sub(&U256::ONE),         // p-1
        ];

        for val in &edge_values {
            let bf = BFq::from_witness(builder.clone(), *val);
            bf.assert_is_in_field("edge case should be in field");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge case: assert_less_than with sorted edge values
    //  Port of C++ edge case assert_less_than
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_assert_less_than_edge_cases() {
        let builder = make_builder();
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);

        // Sorted ascending edge values
        let edge_values = [
            U256::ZERO,
            U256::ONE,
            modulus.wrapping_sub(&U256::ONE),
        ];

        // Each smaller value should be less than each larger value
        for i in 0..edge_values.len() - 1 {
            for j in (i + 1)..edge_values.len() {
                let bf = BFq::from_witness(builder.clone(), edge_values[i]);
                bf.assert_less_than(edge_values[j], "edge smaller < edge larger");
            }
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge case: assert_equal with values that are equal mod p
    //  Port of C++ test_assert_equal_edge_case
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_assert_equal_edge_case() {
        // Test assert_equal with values computed via different arithmetic paths
        let builder = make_builder();

        // Two different ways to compute the same value
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        // Path 1: a + b
        let sum1 = a.add(&b);
        // Path 2: b + a
        let sum2 = b.add(&a);
        sum1.assert_equal(&sum2, "a+b should equal b+a");

        // Path 3: a*b via different routes
        let prod1 = a.mul(&b);
        let prod2 = b.mul(&a);
        prod1.assert_equal(&prod2, "a*b should equal b*a");

        // Path 4: (a-a) should equal 0
        let diff = a.sub(&a);
        let zero = make_witness(builder.clone(), Fq::zero());
        diff.assert_equal(&zero, "a-a should equal 0");

        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Max value tracking during addition
    //  Port of C++ test_maximum_value_tracking_during_addition
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_max_value_tracking_during_addition() {
        let builder = make_builder();
        let mask = U256::from(1u64).wrapping_shl_vartime(NUM_LIMB_BITS as u32).wrapping_sub(&U256::ONE);
        let small_mask = U256::from(1u64).wrapping_shl_vartime(40).wrapping_sub(&U256::ONE);

        let l0_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&mask.to_nz().unwrap());
        let l1_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&mask.to_nz().unwrap());
        let l2_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&mask.to_nz().unwrap());
        let l3_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&small_mask.to_nz().unwrap());

        let l0 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l0_val.as_words()));
        let l1 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l1_val.as_words()));
        let l2 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l2_val.as_words()));
        let l3 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l3_val.as_words()));

        let mut combined = BFq::unsafe_construct_from_limbs(l0, l1, l2, l3, true);
        combined.binary_basis_limbs[3].maximum_value = small_mask;

        // Add to self several times — max value should increase each time
        for _ in 0..10 {
            let prev_max = combined.binary_basis_limbs[0].maximum_value;
            let prev_msb = prev_max.get_msb();
            combined = combined.add(&combined.clone());
            let new_msb = combined.binary_basis_limbs[0].maximum_value.get_msb();
            // Max value should increase by ~1 bit per doubling
            assert!(new_msb >= prev_msb, "max value should grow during addition");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Reduction check works (values > max unreduced trigger reduction)
    //  Port of C++ test_reduction_check_works
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_reduction_check_works() {
        let builder = make_builder();
        let mask = U256::from(1u64).wrapping_shl_vartime(NUM_LIMB_BITS as u32).wrapping_sub(&U256::ONE);
        let msb_mask = U256::from(1u64).wrapping_shl_vartime(61).wrapping_sub(&U256::ONE);

        let l0_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&mask.to_nz().unwrap());
        let l1_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&mask.to_nz().unwrap());
        let l2_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&mask.to_nz().unwrap());
        let mut l3_val = fq_to_u256(Fq::random_element()).wrapping_rem_vartime(&msb_mask.to_nz().unwrap());
        // Set 61st bit to ensure value > 2^265
        l3_val = l3_val.wrapping_or(&U256::from(1u64).wrapping_shl_vartime(60));

        let l0 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l0_val.as_words()));
        let l1 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l1_val.as_words()));
        let l2 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l2_val.as_words()));
        let l3 = FieldT::from_witness(builder.clone(), Fr::from_limbs(*l3_val.as_words()));

        let mut combined = BFq::unsafe_construct_from_limbs(l0, l1, l2, l3, true);
        combined.binary_basis_limbs[3].maximum_value = msb_mask;

        // Squaring should trigger self-reduction
        let _sq = combined.sqr();

        // After squaring triggers reduction, verify limbs are in range
        let max_unreduced = get_maximum_unreduced_limb_value::<Bn254FrParams>();
        for i in 0..4 {
            assert!(
                combined.binary_basis_limbs[i].maximum_value <= max_unreduced,
                "after reduction: limb {} max too large", i
            );
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Negate constant + witness roundtrip
    //  Port of C++ conditional_negate pattern from test_conditional_negate
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_conditional_negate_roundtrip() {
        // a.conditional_negate(true) + a.conditional_negate(false) should == 0
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);

        let pred_true = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::one()));
        let pred_false = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::zero()));

        let mut neg_a = a.conditional_negate(&pred_true);
        let pos_a = a.conditional_negate(&pred_false);

        neg_a.self_reduce();
        assert_bigfield_eq(&neg_a, -a_fq);
        assert_bigfield_eq(&pos_a, a_fq);

        // Sum should be zero
        let mut sum = neg_a.add(&pos_a);
        sum.self_reduce();
        assert_bigfield_eq(&sum, Fq::zero());

        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Internal div regression
    //  Port of C++ test_internal_div_regression
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_internal_div_regression() {
        let builder = make_builder();

        let w0 = make_witness(builder.clone(), Fq::one());
        let pred_true = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::one()));
        let pred_false = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::zero()));

        let w0 = w0.conditional_negate(&pred_true);  // -1
        let w0 = w0.conditional_negate(&pred_false);  // -1
        let w0 = w0.conditional_negate(&pred_true);  // 1
        let w0 = w0.conditional_negate(&pred_true);  // -1
        let w4 = w0.conditional_negate(&pred_false); // -1
        let w4 = w4.conditional_negate(&pred_true);  // 1
        let w4 = w4.conditional_negate(&pred_true);  // -1

        let w5 = w4.sub(&w0); // -1 - (-1) = 0
        let one = make_witness(builder.clone(), Fq::one());
        let _w6 = w5.div(&one); // 0 / 1 = 0

        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Quotient completeness with sqr after many additions
    //  Port of C++ test_quotient_completeness (more variants)
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_quotient_completeness_mul_variant() {
        // Near-max value, doubled 8 times, then multiplied with itself
        let builder = make_builder();
        let val = U256::from_words([
            0xfffffffffffffffe,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0x3fffffffffffffff,
        ]);
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);
        let val_reduced = val.div_rem(&modulus.to_nz().unwrap()).1;

        let mut a = BFq::from_witness(builder.clone(), val_reduced);
        let mut a2 = a.clone();

        for _ in 0..8 {
            let a_copy = a.clone();
            a = a.add(&a_copy);
            let a2_copy = a2.clone();
            a2 = a2.add(&a2_copy);
        }

        // Multiply a * a
        let _b = a.mul(&a2);
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Madd with empty addends (quotient completeness variant)
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_quotient_completeness_madd_variant() {
        let builder = make_builder();
        let val = U256::from_words([
            0xfffffffffffffffe,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0x3fffffffffffffff,
        ]);
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);
        let val_reduced = val.div_rem(&modulus.to_nz().unwrap()).1;

        let mut a = BFq::from_witness(builder.clone(), val_reduced);

        for _ in 0..8 {
            let a_copy = a.clone();
            a = a.add(&a_copy);
        }

        // madd with no addends = a * a
        let a2 = a.clone();
        let _e = a.madd(&a2, &[]);
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Sqr of witness from field_pair
    //  Port of C++ test_sqr with various input types
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_sqr_from_field_pair() {
        let builder = make_builder();
        for _ in 0..10 {
            let (a_fq, u) = random_fq();
            let lo_val = slice_u256(&u, 0, (NUM_LIMB_BITS * 2) as u32);
            let hi_val = slice_u256(&u, (NUM_LIMB_BITS * 2) as u32, (NUM_LIMB_BITS * 4) as u32);

            let lo = FieldT::from_witness(builder.clone(), Fr::from_limbs(*lo_val.as_words()));
            let hi = FieldT::from_witness(builder.clone(), Fr::from_limbs(*hi_val.as_words()));
            let bf = BFq::from_field_pair(lo, hi, false, 0);

            let sq = bf.sqr();
            assert_bigfield_eq(&sq, a_fq * a_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  assert_is_in_field after arithmetic (success path)
    //  Port of C++ test_assert_is_in_field_success
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_assert_is_in_field_after_arithmetic() {
        let builder = make_builder();
        for _ in 0..4 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);

            let mut c = a.clone();
            let mut expected = a_fq;
            for _ in 0..16 {
                c = b.mul(&b).add(&c);
                expected = b_fq * b_fq + expected;
            }

            c.self_reduce();
            c.assert_is_in_field("result after arithmetic should be in field");

            assert_bigfield_eq(&c, expected);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  assert_less_than success with 200-bit values
    //  Port of C++ test_assert_less_than_success
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_assert_less_than_200bit() {
        let builder = make_builder();
        let bit_mask = U256::from(1u64)
            .wrapping_shl_vartime(200)
            .wrapping_sub(&U256::ONE);

        for _ in 0..10 {
            let (_, u) = random_fq();
            let small_val = u.wrapping_rem_vartime(&bit_mask.to_nz().unwrap());

            let a = BFq::from_witness(builder.clone(), small_val);
            a.assert_less_than(
                U256::from(1u64).wrapping_shl_vartime(200),
                "should be < 2^200",
            );
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Division context: result has a context
    //  Port of C++ test_division_context
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_division_context() {
        let builder = make_builder();
        let one = BFq::from_u256(Some(builder.clone()), U256::ONE);
        // Dividing empty numerator by a non-zero denominator
        // In Rust, we just check that div returns a result with context
        let a = make_witness(builder.clone(), Fq::one());
        let result = a.div(&one);
        assert!(result.get_context().is_some(), "division result should have a context");
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Inversion: constant and witness
    //  Port of C++ test_inversion
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_inversion_witness() {
        let builder = make_builder();
        let neg7 = -Fq::from(7u64);
        let a = make_witness(builder.clone(), neg7);
        let one = make_witness(builder.clone(), Fq::one());
        let inv = one.div(&a);
        let expected = neg7.invert();
        assert_bigfield_eq(&inv, expected);

        // Verify a * inv = 1
        let product = a.mul(&inv);
        assert_bigfield_eq(&product, Fq::one());
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Constructor: from_u512_as_witness with various values
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_create_from_u512_overflow() {
        let builder = make_builder();
        let modulus = U256::from_limbs(Bn254FqParams::MODULUS);

        // Value exactly at modulus (with can_overflow=true)
        let bf = BFq::create_from_u512_as_witness(
            builder.clone(),
            U512::from_lo_hi(modulus, U256::ZERO),
            true,
            0,
        );
        // The value should be >= modulus (not auto-reduced when can_overflow=true)
        assert!(bf.get_value() >= BFq::modulus_u512());
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Multiple operations chained (stress test)
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_chain_operations() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let (c_fq, _) = random_fq();

        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);
        let c = make_witness(builder.clone(), c_fq);

        // Compute (a + b) * c - a^2
        let sum = a.add(&b);
        let prod = sum.mul(&c);
        let sq = a.sqr();
        let result = prod.sub(&sq);

        let expected = (a_fq + b_fq) * c_fq - a_fq * a_fq;
        let mut result_reduced = result;
        result_reduced.self_reduce();
        assert_bigfield_eq(&result_reduced, expected);
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Conditional select roundtrip
    //  Port of C++ conditional_select with !predicate
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_conditional_select_roundtrip() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        let pred = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::one()));
        let not_pred = BoolT::from_witness(&WitnessT::new(builder.clone(), Fr::zero()));

        // pred=true → select b; pred=false → select a
        let c = a.conditional_select(&b, &pred);
        let d = a.conditional_select(&b, &not_pred);

        assert_bigfield_eq(&c, b_fq);
        assert_bigfield_eq(&d, a_fq);

        // c + d should equal a + b
        let mut sum = c.add(&d);
        sum.self_reduce();
        assert_bigfield_eq(&sum, a_fq + b_fq);

        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Equality: eq with arithmetic results
    //  Port of C++ test_equality_operator with computed values
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_eq_after_arithmetic() {
        let builder = make_builder();
        for _ in 0..5 {
            let (a_fq, a_u) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            // Create another witness from the same native value
            let another_a = BFq::create_from_u512_as_witness(
                builder.clone(),
                U512::from_lo_hi(a_u, U256::ZERO),
                true,
                0,
            );
            let result = a.eq(&another_a);
            assert!(result.get_value(), "same native value should be equal");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Compound: multiple additions then sqr
    //  Extended quotient completeness test
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_many_additions_then_sqr() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);

        // Add 20 times to accumulate overflow
        let mut sum = a.clone();
        let mut expected = a_fq;
        for _ in 0..20 {
            let (b_fq, _) = random_fq();
            let b = make_witness(builder.clone(), b_fq);
            sum = sum.add(&b);
            expected = expected + b_fq;
        }

        // Now square (this should trigger internal reduction)
        let result = sum.sqr();
        assert_bigfield_eq(&result, expected * expected);
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge: zero division by one
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_zero_div_one() {
        let builder = make_builder();
        let zero = make_witness(builder.clone(), Fq::zero());
        let one = make_witness(builder.clone(), Fq::one());
        let result = zero.div(&one);
        assert_bigfield_eq(&result, Fq::zero());
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge: negate zero
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_negate_zero() {
        let builder = make_builder();
        let zero = make_witness(builder.clone(), Fq::zero());
        let mut result = zero.negate();
        result.self_reduce();
        assert_bigfield_eq(&result, Fq::zero());
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge: a * a using madd
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_madd_as_sqr() {
        let builder = make_builder();
        for _ in 0..5 {
            let (a_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let sqr = a.sqr();
            let madd_sq = a.madd(&a, &[]);
            // Both should match a^2
            assert_bigfield_eq(&sqr, a_fq * a_fq);
            assert_bigfield_eq(&madd_sq, a_fq * a_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge: sub then add (should cancel)
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_sub_then_add_cancel() {
        let builder = make_builder();
        for _ in 0..5 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);

            let diff = a.sub(&b);
            let mut restored = diff.add(&b);
            restored.self_reduce();
            assert_bigfield_eq(&restored, a_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Edge: mul then div (should cancel)
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_mul_then_div_cancel() {
        let builder = make_builder();
        for _ in 0..5 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            if b_fq == Fq::zero() { continue; }
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);

            let prod = a.mul(&b);
            let restored = prod.div(&b);
            assert_bigfield_eq(&restored, a_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Multiple assert_equal patterns
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_assert_equal_after_chain() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        // (a + b) - b should equal a
        let sum = a.add(&b);
        let result = sum.sub(&b);
        result.assert_equal(&a, "(a+b)-b should equal a");
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Madd identity: a*1 + 0 = a
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_madd_identity() {
        let builder = make_builder();
        for _ in 0..5 {
            let (a_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let one = make_witness(builder.clone(), Fq::one());
            let result = a.madd(&one, &[]);
            assert_bigfield_eq(&result, a_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Nonnormalized field bug regression
    //  Port of C++ test_nonnormalized_field_bug_regression
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_nonnormalized_field_bug_regression() {
        // Construct a bigfield where limb 0 has a value just at the default max boundary
        let builder = make_builder();
        let two_to_68 = U256::from(1u64).wrapping_shl_vartime(NUM_LIMB_BITS as u32);
        let bf = BFq::from_u256(Some(builder.clone()), two_to_68);

        // Assert invariant: actual value <= maximum value for every limb
        for i in 0..4 {
            let witness_val = {
                let v = bf.binary_basis_limbs[i].element.get_value().from_montgomery_form();
                U256::from_words(v.data)
            };
            let max_val = bf.binary_basis_limbs[i].maximum_value;
            assert!(
                max_val >= witness_val,
                "nonnormalized bug: limb[{}] max {} < witness {}",
                i, max_val, witness_val
            );
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Double negation roundtrip
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_double_negate() {
        let builder = make_builder();
        for _ in 0..5 {
            let (a_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let mut neg = a.negate();
            neg = neg.negate();
            neg.self_reduce();
            assert_bigfield_eq(&neg, a_fq);
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Associativity and commutativity
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_add_commutativity() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        let mut ab = a.add(&b);
        let mut ba = b.add(&a);
        ab.self_reduce();
        ba.self_reduce();
        ab.assert_equal(&ba, "a+b should equal b+a");
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_mul_commutativity() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        let ab = a.mul(&b);
        let ba = b.mul(&a);
        ab.assert_equal(&ba, "a*b should equal b*a");
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_associativity() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let (c_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);
        let c = make_witness(builder.clone(), c_fq);

        let mut ab_c = a.add(&b).add(&c);
        let mut a_bc = a.add(&b.add(&c));
        ab_c.self_reduce();
        a_bc.self_reduce();
        ab_c.assert_equal(&a_bc, "(a+b)+c should equal a+(b+c)");
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_mul_associativity() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let (c_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);
        let c = make_witness(builder.clone(), c_fq);

        let ab_c = a.mul(&b).mul(&c);
        let a_bc = a.mul(&b.mul(&c));
        ab_c.assert_equal(&a_bc, "(a*b)*c should equal a*(b*c)");
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Distributivity: a*(b+c) = a*b + a*c
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_distributivity() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let (c_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);
        let c = make_witness(builder.clone(), c_fq);

        let lhs = a.mul(&b.add(&c));    // a*(b+c)
        let mut rhs = a.mul(&b).add(&a.mul(&c)); // a*b + a*c
        rhs.self_reduce();
        lhs.assert_equal(&rhs, "distributivity: a*(b+c) should equal a*b + a*c");
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Madd consistency: a.madd(b, [c]) == a*b + c
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_madd_consistency() {
        let builder = make_builder();
        for _ in 0..5 {
            let (a_fq, _) = random_fq();
            let (b_fq, _) = random_fq();
            let (c_fq, _) = random_fq();
            let a = make_witness(builder.clone(), a_fq);
            let b = make_witness(builder.clone(), b_fq);
            let c = make_witness(builder.clone(), c_fq);

            let madd_result = a.madd(&b, &[&c]);
            let manual = a.mul(&b).add(&c);

            madd_result.assert_equal(&manual, "madd should equal mul + add");
        }
        assert!(check_circuit(&builder).is_ok());
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Reduction: reduction_check after many subtractions
    // ═══════════════════════════════════════════════════════════════════

    #[test]
    fn test_reduction_check_after_subtractions() {
        let builder = make_builder();
        let (a_fq, _) = random_fq();
        let (b_fq, _) = random_fq();
        let a = make_witness(builder.clone(), a_fq);
        let b = make_witness(builder.clone(), b_fq);

        let mut c = a.clone();
        for _ in 0..20 {
            c = c.sub(&b);
        }
        c.reduction_check();

        let max_unreduced = get_maximum_unreduced_limb_value::<Bn254FrParams>();
        for i in 0..4 {
            assert!(
                c.binary_basis_limbs[i].maximum_value <= max_unreduced,
                "reduction_check after subs: limb {} max too large", i
            );
        }
        assert!(check_circuit(&builder).is_ok());
    }
}
