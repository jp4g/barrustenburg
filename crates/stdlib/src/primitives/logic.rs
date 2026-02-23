//! Circuit logic (bitwise) operations.
//!
//! Port of `barretenberg/stdlib/primitives/logic/logic.{hpp,cpp}`.
//!
//! Provides `create_logic_constraint`, which computes an AND or XOR of two field
//! elements over a specified number of bits, using plookup tables for the 32-bit
//! chunk operations.

use bbrs_ecc::curves::bn254::Bn254FrParams;
use bbrs_ecc::fields::field::Field;
use bbrs_numeric::uint256::U256Ext;
use bbrs_numeric::U256;

use super::field::FieldT;
use super::plookup;
use super::witness::BuilderRef;

use bbrs_circuit_builder::plookup_tables::types::MultiTableId;

type P = Bn254FrParams;
type Fr = Field<P>;

/// Maximum number of bits for logic constraints (BN254 Fr is ~254 bits).
const MAX_LOGIC_BIT_LENGTH: usize = 252;

/// Compute a bitwise AND or XOR of `a` and `b` over `num_bits` bits.
///
/// If either operand exceeds `num_bits`, the result is truncated. The operands
/// are range-constrained via the plookup decomposition: each 32-bit chunk is
/// implicitly range-checked by the lookup, and the final (possibly smaller)
/// chunk gets an explicit range constraint.
///
/// Port of C++ `logic<Builder>::create_logic_constraint`.
pub fn create_logic_constraint(
    a: &FieldT<P>,
    b: &FieldT<P>,
    num_bits: usize,
    is_xor_gate: bool,
) -> FieldT<P> {
    create_logic_constraint_inner(a, b, num_bits, is_xor_gate, default_get_chunk)
}

/// Default chunk extraction: mask both operands to `chunk_size` low bits.
fn default_get_chunk(left: U256, right: U256, chunk_size: usize) -> (u64, u64) {
    let mask = if chunk_size >= 64 {
        u64::MAX
    } else {
        (1u64 << chunk_size) - 1
    };
    let left_chunk = left.as_words()[0] & mask;
    let right_chunk = right.as_words()[0] & mask;
    (left_chunk, right_chunk)
}

/// Inner implementation of `create_logic_constraint` with a configurable chunk
/// extraction function (used in tests to inject malicious witnesses).
fn create_logic_constraint_inner<F>(
    a: &FieldT<P>,
    b: &FieldT<P>,
    num_bits: usize,
    is_xor_gate: bool,
    get_chunk: F,
) -> FieldT<P>
where
    F: Fn(U256, U256, usize) -> (u64, u64),
{
    assert!(
        num_bits <= MAX_LOGIC_BIT_LENGTH,
        "num_bits exceeds maximum logic bit length"
    );
    assert!(num_bits > 0, "num_bits must be positive");

    // ── Constant path ──────────────────────────────────────────────────
    if a.is_constant() && b.is_constant() {
        let a_raw = a.get_value().from_montgomery_form();
        let b_raw = b.get_value().from_montgomery_form();
        let a_u256 = U256::from_words(a_raw.data);
        let b_u256 = U256::from_words(b_raw.data);

        let result = if is_xor_gate {
            a_u256.bitxor(&b_u256)
        } else {
            a_u256.bitand(&b_u256)
        };

        // Truncate to num_bits
        let mask = U256::ONE
            .wrapping_shl_vartime(num_bits as u32)
            .wrapping_sub(&U256::ONE);
        let truncated = result.bitand(&mask);

        let result_field = Fr::from_raw(truncated.limbs()).to_montgomery_form();
        return FieldT::from_field(result_field);
    }

    // ── Convert constant operand to witness ────────────────────────────
    if a.is_constant() && !b.is_constant() {
        let ctx = b.get_context().as_ref().unwrap().clone();
        let a_witness = constant_to_witness(&ctx, a);
        return create_logic_constraint_inner(&a_witness, b, num_bits, is_xor_gate, get_chunk);
    }

    if !a.is_constant() && b.is_constant() {
        let ctx = a.get_context().as_ref().unwrap().clone();
        let b_witness = constant_to_witness(&ctx, b);
        return create_logic_constraint_inner(a, &b_witness, num_bits, is_xor_gate, get_chunk);
    }

    // ── Variable path: both operands are witnesses ─────────────────────
    let ctx = a
        .get_context()
        .as_ref()
        .expect("non-constant must have context")
        .clone();

    let num_chunks = (num_bits + 31) / 32;

    let a_raw = a.get_value().from_montgomery_form();
    let b_raw = b.get_value().from_montgomery_form();
    let mut left = U256::from_words(a_raw.data);
    let mut right = U256::from_words(b_raw.data);

    let mut a_accumulator = FieldT::with_context(ctx.clone()); // zero
    let mut b_accumulator = FieldT::with_context(ctx.clone()); // zero
    let mut res = FieldT::with_context(ctx.clone()); // zero

    let mut scaling = Fr::one(); // 2^(32*i), starts at 1

    for i in 0..num_chunks {
        let chunk_size = if i != num_chunks - 1 {
            32
        } else {
            num_bits - i * 32
        };

        let (left_chunk, right_chunk) = get_chunk(left, right, chunk_size);

        let a_chunk = FieldT::from_witness(ctx.clone(), Fr::from(left_chunk));
        let b_chunk = FieldT::from_witness(ctx.clone(), Fr::from(right_chunk));

        let multi_table_id = if is_xor_gate {
            MultiTableId::UINT32_XOR
        } else {
            MultiTableId::UINT32_AND
        };

        let result_chunk =
            plookup::read_from_2_to_1_table(multi_table_id, &a_chunk, &b_chunk);

        let scaling_ft = FieldT::from_field(scaling);
        a_accumulator = a_accumulator + a_chunk.clone() * scaling_ft.clone();
        b_accumulator = b_accumulator + b_chunk.clone() * scaling_ft.clone();

        if chunk_size != 32 {
            a_chunk.create_range_constraint(
                chunk_size,
                "stdlib logic: bad range on final chunk of left operand",
            );
            b_chunk.create_range_constraint(
                chunk_size,
                "stdlib logic: bad range on final chunk of right operand",
            );
        }

        res = res + result_chunk * scaling_ft;

        left = left.wrapping_shr_vartime(32);
        right = right.wrapping_shr_vartime(32);
        scaling = scaling * Fr::from(1u64 << 32);
    }

    a.assert_equal(
        &a_accumulator,
        "stdlib logic: failed to reconstruct left operand",
    );
    b.assert_equal(
        &b_accumulator,
        "stdlib logic: failed to reconstruct right operand",
    );

    res
}

/// Convert a constant FieldT to a witness by registering it as a constant variable.
fn constant_to_witness(ctx: &BuilderRef<P>, val: &FieldT<P>) -> FieldT<P> {
    let idx = ctx.borrow_mut().put_constant_variable(val.get_value());
    FieldT::from_witness_index(ctx.clone(), idx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use std::cell::RefCell;
    use std::rc::Rc;

    fn make_builder() -> BuilderRef<P> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    /// Port of C++ `LogicTest::TestCorrectLogic`.
    ///
    /// Validates AND and XOR over bit widths 8..248 (step 8) with both
    /// witness and constant inputs.
    #[test]
    fn test_correct_logic() {
        let builder = make_builder();

        for num_bits in (8..=248).step_by(8) {
            let mask = U256::ONE
                .wrapping_shl_vartime(num_bits as u32)
                .wrapping_sub(&U256::ONE);

            // Use deterministic test values derived from num_bits
            let a_u256 = U256::from_words([
                0x1234_5678_9ABC_DEF0u64.wrapping_mul(num_bits as u64),
                0xFEDC_BA98_7654_3210u64.wrapping_mul(num_bits as u64),
                0xAAAA_BBBB_CCCC_DDDDu64.wrapping_mul(num_bits as u64),
                0x1111_2222_3333_4444u64.wrapping_mul(num_bits as u64),
            ])
            .bitand(&mask);
            let b_u256 = U256::from_words([
                0xDEAD_BEEF_CAFE_BABEu64.wrapping_mul(num_bits as u64),
                0x0123_4567_89AB_CDEFu64.wrapping_mul(num_bits as u64),
                0x5555_6666_7777_8888u64.wrapping_mul(num_bits as u64),
                0x9999_0000_AAAA_BBBBu64.wrapping_mul(num_bits as u64),
            ])
            .bitand(&mask);

            let and_expected = a_u256.bitand(&b_u256);
            let xor_expected = a_u256.bitxor(&b_u256);

            let a_fr = Fr::from_raw(a_u256.limbs()).to_montgomery_form();
            let b_fr = Fr::from_raw(b_u256.limbs()).to_montgomery_form();

            let x = FieldT::from_witness(builder.clone(), a_fr);
            let y = FieldT::from_witness(builder.clone(), b_fr);

            let x_const = FieldT::constant_with_context(builder.clone(), a_fr);
            let y_const = FieldT::constant_with_context(builder.clone(), b_fr);

            // Both witnesses
            let and_result = create_logic_constraint(&x, &y, num_bits, false);
            let xor_result = create_logic_constraint(&x, &y, num_bits, true);

            // Left constant
            let and_result_left_const =
                create_logic_constraint(&x_const, &y, num_bits, false);
            let xor_result_left_const =
                create_logic_constraint(&x_const, &y, num_bits, true);

            // Right constant
            let and_result_right_const =
                create_logic_constraint(&x, &y_const, num_bits, false);
            let xor_result_right_const =
                create_logic_constraint(&x, &y_const, num_bits, true);

            // Both constants
            let and_result_both_const =
                create_logic_constraint(&x_const, &y_const, num_bits, false);
            let xor_result_both_const =
                create_logic_constraint(&x_const, &y_const, num_bits, true);

            let and_expected_fr = Fr::from_raw(and_expected.limbs()).to_montgomery_form();
            let xor_expected_fr = Fr::from_raw(xor_expected.limbs()).to_montgomery_form();

            assert_eq!(and_result.get_value(), and_expected_fr, "AND failed for num_bits={num_bits}");
            assert_eq!(and_result_left_const.get_value(), and_expected_fr);
            assert_eq!(and_result_right_const.get_value(), and_expected_fr);
            assert_eq!(and_result_both_const.get_value(), and_expected_fr);

            assert_eq!(xor_result.get_value(), xor_expected_fr, "XOR failed for num_bits={num_bits}");
            assert_eq!(xor_result_left_const.get_value(), xor_expected_fr);
            assert_eq!(xor_result_right_const.get_value(), xor_expected_fr);
            assert_eq!(xor_result_both_const.get_value(), xor_expected_fr);
        }

        UltraCircuitChecker::check(&mut builder.borrow_mut())
            .expect("circuit check failed");
    }

    /// Port of C++ `LogicTest::LargeOperands`.
    ///
    /// Operands are 48-bit but the constraint width is 40 bits, so the
    /// range constraints on the final chunks should fail.
    #[test]
    fn test_large_operands() {
        let builder = make_builder();

        let mask_48 = U256::ONE
            .wrapping_shl_vartime(48)
            .wrapping_sub(&U256::ONE);

        let a_u256 = U256::from_words([0xABCD_1234_5678u64, 0, 0, 0]).bitand(&mask_48);
        let b_u256 = U256::from_words([0xDEAD_BEEF_CAFEu64, 0, 0, 0]).bitand(&mask_48);

        let a_fr = Fr::from_raw(a_u256.limbs()).to_montgomery_form();
        let b_fr = Fr::from_raw(b_u256.limbs()).to_montgomery_form();

        let x = FieldT::from_witness(builder.clone(), a_fr);
        let y = FieldT::from_witness(builder.clone(), b_fr);

        let mask_40 = (1u64 << 40) - 1;
        let xor_expected = (a_u256.as_words()[0] ^ b_u256.as_words()[0]) & mask_40;
        let and_expected = (a_u256.as_words()[0] & b_u256.as_words()[0]) & mask_40;

        let xor_result = create_logic_constraint(&x, &y, 40, true);
        let and_result = create_logic_constraint(&x, &y, 40, false);

        // The values are correct (truncated to 40 bits)...
        assert_eq!(
            xor_result.get_value(),
            Fr::from(xor_expected),
            "XOR value mismatch"
        );
        assert_eq!(
            and_result.get_value(),
            Fr::from(and_expected),
            "AND value mismatch"
        );

        // ...but the circuit should fail because the operands exceed 40 bits
        let result = UltraCircuitChecker::check(&mut builder.borrow_mut());
        assert!(result.is_err(), "circuit should fail with large operands");
    }

    /// Port of C++ `LogicTest::DifferentWitnessSameResult`.
    ///
    /// Injects false witness chunks that produce the same XOR result but
    /// don't match the original operands. The circuit should detect this.
    #[test]
    fn test_different_witness_same_result() {
        let builder = make_builder();

        let a: u64 = 3758096391;
        let b: u64 = 2147483649;
        let xor_expected = a ^ b;

        let x = FieldT::from_witness(builder.clone(), Fr::from(a));
        let y = FieldT::from_witness(builder.clone(), Fr::from(b));

        // Inject bad chunks that produce the same XOR result
        let get_bad_chunk = |_left: U256, _right: U256, _chunk_size: usize| -> (u64, u64) {
            (2684354565u64, 3221225475u64)
        };

        let xor_result =
            create_logic_constraint_inner(&x, &y, 32, true, get_bad_chunk);

        assert_eq!(
            xor_result.get_value(),
            Fr::from(xor_expected),
            "XOR value should match"
        );

        // Circuit should fail because the chunks don't reconstruct to the original operands
        let result = UltraCircuitChecker::check(&mut builder.borrow_mut());
        assert!(
            result.is_err(),
            "circuit should fail with malicious witnesses"
        );
    }

    /// Test that constant-only inputs produce no gates and correct results.
    #[test]
    fn test_constant_logic() {
        let builder = make_builder();

        let a_val: u64 = 0xFF00_FF00;
        let b_val: u64 = 0x0F0F_0F0F;

        let a = FieldT::constant_with_context(builder.clone(), Fr::from(a_val));
        let b = FieldT::constant_with_context(builder.clone(), Fr::from(b_val));

        let and_result = create_logic_constraint(&a, &b, 32, false);
        let xor_result = create_logic_constraint(&a, &b, 32, true);

        assert!(and_result.is_constant());
        assert!(xor_result.is_constant());

        assert_eq!(and_result.get_value(), Fr::from(a_val & b_val));
        assert_eq!(xor_result.get_value(), Fr::from(a_val ^ b_val));
    }
}
