//! Circuit byte-array type.
//!
//! Port of `barretenberg/stdlib/primitives/byte_array/byte_array.hpp` and `byte_array.cpp`.
//!
//! A `ByteArrayT<P>` represents a dynamic array of bytes in-circuit. Each byte
//! is a `FieldT<P>` range-constrained to 8 bits. It supports construction from
//! native values (`Vec<u8>`, `&str`, or `FieldT`) and conversion back to `FieldT`,
//! as well as classical vector operations like slicing, appending, and reversing.
//!
//! Used in hashing primitives (SHA256, Blake2s, Keccak, etc.).

use std::fmt;

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use super::field::FieldT;
use super::witness::{BuilderRef, WitnessT};

/// A dynamic array of bytes in-circuit.
///
/// Each element is a `FieldT<P>` constrained to lie in `[0, 255]`.
pub struct ByteArrayT<P: FieldParams> {
    context: Option<BuilderRef<P>>,
    values: Vec<FieldT<P>>,
}

impl<P: FieldParams> Clone for ByteArrayT<P> {
    fn clone(&self) -> Self {
        Self {
            context: self.context.clone(),
            values: self.values.clone(),
        }
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Constructors
// ════════════════════════════════════════════════════════════════════════

impl<P: FieldParams> ByteArrayT<P> {
    // ── Internal constructors (no constraints) ──────────────────────

    /// Internal: wrap an existing vector of (already-constrained) field
    /// elements. No range constraints are added.
    fn from_values(ctx: Option<BuilderRef<P>>, values: Vec<FieldT<P>>) -> Self {
        Self {
            context: ctx,
            values,
        }
    }

    // ── Public constructors ─────────────────────────────────────────

    /// Create an empty byte array with a context.
    pub fn empty(ctx: BuilderRef<P>) -> Self {
        Self {
            context: Some(ctx),
            values: Vec::new(),
        }
    }

    /// Create a byte array from a slice of native bytes.
    ///
    /// Each byte is instantiated as a **witness** (not a constant) and
    /// range-constrained to 8 bits.
    pub fn from_bytes(ctx: BuilderRef<P>, input: &[u8]) -> Self {
        let mut values = Vec::with_capacity(input.len());
        for &byte in input {
            let value = FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Field::from(byte as u64)));
            value.create_range_constraint(8, "byte_array: vector entry larger than 1 byte.");
            values.push(value);
        }
        Self {
            context: Some(ctx),
            values,
        }
    }

    /// Create a byte array from a string.
    ///
    /// Each character is instantiated as a **witness** (not a constant) and
    /// range-constrained to 8 bits.
    pub fn from_str(ctx: BuilderRef<P>, input: &str) -> Self {
        Self::from_bytes(ctx, input.as_bytes())
    }

    /// Create a byte array from constant values without adding range constraints.
    ///
    /// Safe for padding and other constant data — constants cannot be
    /// manipulated by the prover.
    pub fn from_constants(ctx: BuilderRef<P>, input: &[u8]) -> Self {
        let mut values = Vec::with_capacity(input.len());
        for &byte in input {
            values.push(FieldT::constant_with_context(ctx.clone(), Field::from(byte as u64)));
        }
        Self {
            context: Some(ctx),
            values,
        }
    }

    /// Convenience: create constant zero-padding of the given length.
    pub fn constant_padding(ctx: BuilderRef<P>, num_bytes: usize, value: u8) -> Self {
        Self::from_constants(ctx, &vec![value; num_bytes])
    }

    /// Create a byte array of length `num_bytes` from a field element.
    ///
    /// The field element is decomposed into big-endian bytes. Each byte is
    /// range-constrained to 8 bits. For `num_bytes == 32`, an additional
    /// modular-reduction check ensures the decomposition is unique (i.e. the
    /// value lies in `[0, p-1]`).
    ///
    /// `test_val` may be provided to inject a specific u256 value for testing
    /// overflow scenarios (where the value exceeds the field modulus).
    pub fn from_field(input: &FieldT<P>, num_bytes: usize, test_val: Option<[u64; 4]>) -> Self {
        const MAX_NUM_BYTES: usize = 32;
        const MIDPOINT: usize = MAX_NUM_BYTES / 2;

        assert!(num_bytes <= MAX_NUM_BYTES, "byte_array: num_bytes > 32");

        // Get native value as u256 limbs (little-endian u64 words)
        let value: [u64; 4] = if let Some(tv) = test_val {
            tv
        } else {
            input.get_value().from_montgomery_form().data
        };

        let mut values = vec![FieldT::from_u64(0); num_bytes];

        let context = input.get_context().clone();

        // Accumulators for high and low halves (used for reconstruction check)
        let mut accumulator_lo: Vec<FieldT<P>> = Vec::new();
        let mut accumulator_hi: Vec<FieldT<P>> = Vec::new();

        for i in 0..num_bytes {
            let bit_start = (num_bytes - i - 1) * 8;

            // scaling_factor = 2^bit_start
            let scaling_factor = u256_pow2(bit_start);
            let scaling_field = FieldT::from_field(Field::from_limbs(scaling_factor));

            // Extract the current byte from the u256 value
            let byte_val = u256_slice(&value, bit_start, bit_start + 8);

            // Create the byte element: witness if input is witness, constant otherwise
            let byte = if input.is_constant() {
                FieldT::from_field(Field::from(byte_val as u64))
            } else {
                let ctx = context.as_ref().expect("non-constant must have context");
                FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Field::from(byte_val as u64)))
            };
            byte.create_range_constraint(8, "byte_array: byte extraction failed.");
            values[i] = byte.clone();

            if i < MIDPOINT {
                accumulator_hi.push(&scaling_field * &byte);
            } else {
                accumulator_lo.push(&scaling_field * &byte);
            }
        }

        // Reconstruct the high and low limbs from the byte decomposition
        let reconstructed_lo = FieldT::accumulate(&accumulator_lo);
        let reconstructed_hi = FieldT::accumulate(&accumulator_hi);
        let reconstructed = &reconstructed_hi + &reconstructed_lo;

        // Ensure reconstruction matches input
        input.assert_equal(&reconstructed, "byte_array: reconstruction failed");

        // Handle non-unique decomposition for 32-byte case
        if num_bytes == 32 {
            // modulus - 1, split into low and high 128-bit halves
            let modulus_minus_one = u256_sub(&P::MODULUS, &[1, 0, 0, 0]);
            let s_lo = u256_lo128(&modulus_minus_one);   // (r-1) mod 2^128
            let s_hi = u256_hi128(&modulus_minus_one);   // (r-1) >> 128
            let shift: [u64; 4] = [0, 0, 1, 0]; // 2^128

            // diff_lo = s_lo + 2^128 - reconstructed_lo
            let neg_recon_lo = -reconstructed_lo.clone();
            let diff_lo = &(&neg_recon_lo + &FieldT::from_field(Field::from_limbs(s_lo)))
                + &FieldT::from_field(Field::from_limbs(shift));

            let diff_lo_value = diff_lo.get_value().from_montgomery_form().data;

            // Extract the "borrow" bit: bit 128 of diff_lo
            let diff_lo_hi_value = u256_shr128(&diff_lo_value);
            let diff_lo_hi = if input.is_constant() {
                FieldT::from_field(Field::from_limbs(diff_lo_hi_value))
            } else {
                let ctx = context.as_ref().unwrap();
                FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Field::from_limbs(diff_lo_hi_value)))
            };
            diff_lo_hi.create_range_constraint(1, "byte_array: y_overlap is not a bit");

            // Extract first 128 bits of diff_lo
            let lo_mask: [u64; 4] = [u64::MAX, u64::MAX, 0, 0]; // 2^128 - 1
            let lo = u256_and(&diff_lo_value, &lo_mask);
            let diff_lo_lo = if input.is_constant() {
                FieldT::from_field(Field::from_limbs(lo))
            } else {
                let ctx = context.as_ref().unwrap();
                FieldT::from_witness_t(&WitnessT::new(ctx.clone(), Field::from_limbs(lo)))
            };
            diff_lo_lo.create_range_constraint(128, "byte_array: y_remainder doesn't fit in 128 bits.");

            // Constrain: diff_lo == diff_lo_lo + diff_lo_hi * 2^128
            let shift_field = FieldT::from_field(Field::from_limbs(shift));
            diff_lo.assert_equal(
                &(&diff_lo_lo + &(&diff_lo_hi * &shift_field)),
                "byte_array: diff_lo decomposition failed",
            );

            let overlap = &(-diff_lo_hi.clone()) + &FieldT::from_u64(1);

            // reconstructed_hi is always a multiple of 2^128 by construction.
            // Compute diff_hi = s_hi - reconstructed_hi/shift - overlap
            // Using add_two: (-reconstructed_hi / shift).add_two(s_hi, -overlap)
            let neg_recon_hi = -reconstructed_hi.clone();
            let inv_shift = FieldT::from_field(Field::from_limbs(shift).invert());
            let recon_hi_div_shift = &neg_recon_hi * &inv_shift;
            let diff_hi =
                recon_hi_div_shift.add_two(&FieldT::from_field(Field::from_limbs(s_hi)), &(-overlap));
            diff_hi.create_range_constraint(128, "byte_array: y_hi doesn't fit in 128 bits.");
        }

        Self {
            context,
            values,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Conversion to FieldT
// ════════════════════════════════════════════════════════════════════════

impl<P: FieldParams> ByteArrayT<P> {
    /// Convert this byte array to a field element (big-endian).
    ///
    /// The returned field equals `sum(byte[i] * 256^(n-1-i))`.
    pub fn to_field(&self) -> FieldT<P> {
        let bytes = self.values.len();
        let mut scaled_values = Vec::with_capacity(bytes);

        for i in 0..bytes {
            let bit_start = (bytes - i - 1) * 8;
            let scaling = u256_pow2(bit_start);
            let scaling_field = FieldT::from_field(Field::from_limbs(scaling));
            scaled_values.push(&self.values[i] * &scaling_field);
        }

        FieldT::accumulate(&scaled_values)
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Operations
// ════════════════════════════════════════════════════════════════════════

impl<P: FieldParams> ByteArrayT<P> {
    /// Append the contents of `other` to this byte array.
    pub fn write(&mut self, other: &ByteArrayT<P>) {
        self.values.extend_from_slice(&other.values);
    }

    /// Overwrite bytes starting at `index` with the contents of `other`.
    pub fn write_at(&mut self, other: &ByteArrayT<P>, index: usize) {
        assert!(
            index + other.values.len() <= self.values.len(),
            "byte_array::write_at: out of bounds"
        );
        for i in 0..other.values.len() {
            self.values[i + index] = other.values[i].clone();
        }
    }

    /// Slice bytes from `offset` to the end.
    pub fn slice(&self, offset: usize) -> ByteArrayT<P> {
        assert!(offset <= self.values.len(), "byte_array::slice: offset out of bounds");
        ByteArrayT::from_values(self.context.clone(), self.values[offset..].to_vec())
    }

    /// Slice `length` bytes starting at `offset`.
    pub fn slice_with_length(&self, offset: usize, length: usize) -> ByteArrayT<P> {
        assert!(offset <= self.values.len(), "byte_array::slice: offset out of bounds");
        assert!(
            length <= self.values.len() - offset,
            "byte_array::slice: length out of bounds"
        );
        ByteArrayT::from_values(self.context.clone(), self.values[offset..offset + length].to_vec())
    }

    /// Reverse the order of bytes.
    pub fn reverse(&self) -> ByteArrayT<P> {
        if self.values.is_empty() {
            return self.clone();
        }
        let mut bytes = vec![FieldT::from_u64(0); self.values.len()];
        let mut offset = bytes.len() - 1;
        for i in 0..bytes.len() {
            bytes[offset] = self.values[i].clone();
            if i < bytes.len() - 1 {
                offset -= 1;
            }
        }
        ByteArrayT::from_values(self.context.clone(), bytes)
    }

    /// Number of bytes.
    pub fn size(&self) -> usize {
        self.values.len()
    }

    /// Access the underlying byte elements.
    pub fn bytes(&self) -> &[FieldT<P>] {
        &self.values
    }

    /// Get the builder context.
    pub fn get_context(&self) -> &Option<BuilderRef<P>> {
        &self.context
    }

    /// Index into the byte array.
    pub fn index(&self, i: usize) -> &FieldT<P> {
        assert!(i < self.values.len(), "byte_array: index out of bounds");
        &self.values[i]
    }

    /// Get the native byte values (out-of-circuit).
    pub fn get_value(&self) -> Vec<u8> {
        self.values
            .iter()
            .map(|v| {
                let val = v.get_value().from_montgomery_form().data;
                val[0] as u8
            })
            .collect()
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Display
// ════════════════════════════════════════════════════════════════════════

impl<P: FieldParams> fmt::Display for ByteArrayT<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[")?;
        for byte in self.get_value() {
            write!(f, " {:02x}", byte)?;
        }
        write!(f, " ]")
    }
}

// ════════════════════════════════════════════════════════════════════════
//  u256 helpers (little-endian [u64; 4])
// ════════════════════════════════════════════════════════════════════════

/// Extract bits `[bit_start..bit_end)` from a u256, returning as u64.
/// `bit_end - bit_start` must be <= 64.
fn u256_slice(val: &[u64; 4], bit_start: usize, bit_end: usize) -> u64 {
    assert!(bit_end > bit_start && bit_end - bit_start <= 64);
    let word = bit_start / 64;
    let bit_offset = bit_start % 64;
    let mask = (1u64 << (bit_end - bit_start)) - 1;

    if word >= 4 {
        return 0;
    }

    let lo = val[word] >> bit_offset;
    if bit_offset + (bit_end - bit_start) > 64 && word + 1 < 4 {
        let hi = val[word + 1] << (64 - bit_offset);
        (lo | hi) & mask
    } else {
        lo & mask
    }
}

/// 2^n as a u256
fn u256_pow2(n: usize) -> [u64; 4] {
    assert!(n < 256);
    let word = n / 64;
    let bit = n % 64;
    let mut result = [0u64; 4];
    result[word] = 1u64 << bit;
    result
}

/// Subtract b from a as u256 (a - b), wrapping.
fn u256_sub(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    let mut result = [0u64; 4];
    let mut borrow = 0u64;
    for i in 0..4 {
        let (r, b1) = a[i].overflowing_sub(b[i]);
        let (r2, b2) = r.overflowing_sub(borrow);
        result[i] = r2;
        borrow = (b1 as u64) + (b2 as u64);
    }
    result
}

/// Low 128 bits of a u256 as a u256
fn u256_lo128(val: &[u64; 4]) -> [u64; 4] {
    [val[0], val[1], 0, 0]
}

/// High 128 bits of a u256, shifted right by 128 bits
fn u256_hi128(val: &[u64; 4]) -> [u64; 4] {
    [val[2], val[3], 0, 0]
}

/// Right-shift by 128 bits
fn u256_shr128(val: &[u64; 4]) -> [u64; 4] {
    [val[2], val[3], 0, 0]
}

/// Bitwise AND of two u256s
fn u256_and(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    [a[0] & b[0], a[1] & b[1], a[2] & b[2], a[3] & b[3]]
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
    use std::rc::Rc;

    type Fr = Field<Bn254FrParams>;
    type FrField = FieldT<Bn254FrParams>;
    type FrByteArray = ByteArrayT<Bn254FrParams>;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    fn check_circuit(builder: &BuilderRef<Bn254FrParams>) -> Result<(), String> {
        UltraCircuitChecker::check(&mut builder.borrow_mut())
    }

    /// Slice a Field to fit in n bytes (mask to lower n*8 bits).
    fn slice_to_n_bytes(value: Fr, n: usize) -> Fr {
        let val = value.from_montgomery_form().data;
        // mask = (1 << (8*n)) - 1
        if n >= 32 {
            return value;
        }
        let total_bits = n * 8;
        let mut mask = [0u64; 4];
        for i in 0..4 {
            let word_start = i * 64;
            if word_start >= total_bits {
                break;
            }
            let bits_in_word = std::cmp::min(64, total_bits - word_start);
            if bits_in_word == 64 {
                mask[i] = u64::MAX;
            } else {
                mask[i] = (1u64 << bits_in_word) - 1;
            }
        }
        let masked = [val[0] & mask[0], val[1] & mask[1], val[2] & mask[2], val[3] & mask[3]];
        Fr::from_limbs(masked)
    }

    /// Check that byte decomposition is correct: reconstructing from bytes matches original.
    fn check_byte_decomposition(arr: &FrByteArray, original_val: &FrField) {
        let num_bytes = arr.size();
        let mut reconstructed = Fr::zero();
        for i in 0..num_bytes {
            let byte_val = arr.index(i).get_value();
            let shift = u256_pow2((num_bytes - 1 - i) * 8);
            reconstructed = reconstructed + byte_val * Fr::from_limbs(shift);
        }
        assert_eq!(
            original_val.get_value(),
            reconstructed,
            "byte decomposition mismatch"
        );
    }

    // ── Test: reverse ───────────────────────────────────────────────

    #[test]
    fn test_reverse() {
        let builder = make_builder();

        let expected: Vec<u8> = vec![0x04, 0x03, 0x02, 0x01];
        let arr = FrByteArray::from_bytes(builder, &[0x01, 0x02, 0x03, 0x04]);

        let reversed = arr.reverse();

        assert_eq!(arr.size(), 4);
        assert_eq!(reversed.get_value(), expected);
    }

    // ── Test: byte decomposition < 32 bytes (witness) ───────────────

    #[test]
    fn test_byte_decomposition_less_than_32_bytes() {
        for num_bytes in 1..32 {
            let builder = make_builder();

            let raw_val = Fr::random_element();
            let expected_val = slice_to_n_bytes(raw_val, num_bytes);

            let field = FrField::from_witness(builder.clone(), expected_val);
            let byte_arr = FrByteArray::from_field(&field, num_bytes, None);
            assert_eq!(byte_arr.size(), num_bytes);

            check_byte_decomposition(&byte_arr, &field);

            // Convert back to field
            let reconstructed_field = byte_arr.to_field();
            assert_eq!(reconstructed_field.get_value(), expected_val);

            assert!(check_circuit(&builder).is_ok());
        }
    }

    // ── Test: byte decomposition < 32 bytes (constant) ──────────────

    #[test]
    fn test_byte_decomposition_less_than_32_bytes_const() {
        for num_bytes in 1..32 {
            let builder = make_builder();
            let gates_start = builder.borrow().base.num_gates();

            let raw_val = Fr::random_element();
            let expected_val = slice_to_n_bytes(raw_val, num_bytes);

            let field = FrField::constant_with_context(builder.clone(), expected_val);
            let byte_arr = FrByteArray::from_field(&field, num_bytes, None);
            assert_eq!(byte_arr.size(), num_bytes);

            check_byte_decomposition(&byte_arr, &field);

            // Convert back to field
            let reconstructed_field = byte_arr.to_field();
            assert_eq!(reconstructed_field.get_value(), expected_val);

            // No gates should be added for constants
            assert_eq!(builder.borrow().base.num_gates(), gates_start);
        }
    }

    // ── Test: byte decomposition 32 bytes (witness) ─────────────────

    #[test]
    fn test_byte_decomposition_32_bytes() {
        let builder = make_builder();

        let test_val = FrField::from_witness(builder.clone(), Fr::random_element());
        let arr = FrByteArray::from_field(&test_val, 32, None);

        check_byte_decomposition(&arr, &test_val);
        assert!(check_circuit(&builder).is_ok());

        // Test overflow: value >= modulus should fail circuit check
        {
            let builder2 = make_builder();

            // modulus + 100
            let modulus = Bn254FrParams::MODULUS;
            let overflowing = u256_add(&modulus, &[100, 0, 0, 0]);

            let test_val2 = FrField::from_witness(builder2.clone(), Fr::from_limbs(overflowing));
            let _failure_array = FrByteArray::from_field(&test_val2, 32, Some(overflowing));

            assert!(check_circuit(&builder2).is_err());
        }
    }

    // ── Test: byte decomposition 32 bytes (constant) ────────────────

    #[test]
    fn test_byte_decomposition_32_bytes_const() {
        let builder = make_builder();
        let gates_start = builder.borrow().base.num_gates();

        let test_val = FrField::constant_with_context(builder.clone(), Fr::random_element());
        let arr = FrByteArray::from_field(&test_val, 32, None);

        check_byte_decomposition(&arr, &test_val);

        // No gates should be added for constants
        assert_eq!(builder.borrow().base.num_gates(), gates_start);
    }

    // ── Test: input/output consistency ──────────────────────────────

    #[test]
    fn test_input_output_consistency() {
        let builder = make_builder();

        let a_raw = Fr::random_element();
        let b_raw = Fr::random_element();
        let a_expected = slice_to_n_bytes(a_raw, 31);
        let b_expected = slice_to_n_bytes(b_raw, 31);

        let a = FrField::from_witness(builder.clone(), a_expected);
        let b = FrField::from_witness(builder.clone(), b_expected);

        let a_bytes = FrByteArray::from_field(&a, 31, None);
        let b_bytes = FrByteArray::from_field(&b, 31, None);

        // Build byte_array by writing
        let mut arr = FrByteArray::from_bytes(builder.clone(), &[]);
        arr.write(&a_bytes);
        arr.write(&b_bytes);

        assert_eq!(arr.size(), 62);

        let a_result = arr.slice_with_length(0, 31).to_field();
        let b_result = arr.slice(31).to_field();

        assert_eq!(a_result.get_value(), a_expected);
        assert_eq!(b_result.get_value(), b_expected);

        assert!(check_circuit(&builder).is_ok());
    }

    // ── Test: conversion to field ──────────────────────────────────

    #[test]
    fn test_conversion_to_field() {
        for arr_length in 1..32 {
            let builder = make_builder();

            // Generate random bytes from a random field element
            let random_val = Fr::random_element().from_montgomery_form().data;
            let all_bytes: Vec<u8> = random_val.iter()
                .flat_map(|limb| limb.to_le_bytes())
                .collect();
            let native_bytes: Vec<u8> = all_bytes[..arr_length].to_vec();

            // Create byte_array from vector (witnesses)
            let test_array = FrByteArray::from_bytes(builder.clone(), &native_bytes);

            // Convert to field
            let represented_field = test_array.to_field();

            // Compute expected value manually (big-endian)
            let mut expected_limbs = [0u64; 4];
            for &byte in &native_bytes {
                // Shift left by 8 bits
                expected_limbs = u256_shl8(&expected_limbs);
                expected_limbs[0] |= byte as u64;
            }
            let expected_val = Fr::from_limbs(expected_limbs);

            assert_eq!(represented_field.get_value(), expected_val);
            assert!(check_circuit(&builder).is_ok());
        }
    }

    // ── Test: Display (ostream operator) ────────────────────────────

    #[test]
    fn test_display() {
        let builder = make_builder();

        let arr = FrByteArray::from_bytes(builder, &[0x01, 0x02, 0x03, 0x61]);
        let output = format!("{}", arr);
        assert_eq!(output, "[ 01 02 03 61 ]");
    }

    // ── Test: empty byte array ─────────────────────────────────────

    #[test]
    fn test_empty_byte_array() {
        let builder = make_builder();
        let arr = FrByteArray::empty(builder.clone());
        assert_eq!(arr.size(), 0);
        assert_eq!(arr.get_value(), Vec::<u8>::new());
    }

    // ── Test: write_at ─────────────────────────────────────────────

    #[test]
    fn test_write_at() {
        let builder = make_builder();
        let mut arr = FrByteArray::from_bytes(builder.clone(), &[0x00, 0x00, 0x00, 0x00]);
        let patch = FrByteArray::from_bytes(builder.clone(), &[0xAA, 0xBB]);
        arr.write_at(&patch, 1);

        let values = arr.get_value();
        assert_eq!(values, vec![0x00, 0xAA, 0xBB, 0x00]);
    }

    // ── u256 helper ────────────────────────────────────────────────

    /// Shift left by 8 bits
    fn u256_shl8(val: &[u64; 4]) -> [u64; 4] {
        let mut result = [0u64; 4];
        result[0] = val[0] << 8;
        result[1] = (val[1] << 8) | (val[0] >> 56);
        result[2] = (val[2] << 8) | (val[1] >> 56);
        result[3] = (val[3] << 8) | (val[2] >> 56);
        result
    }

    /// Add two u256 values
    fn u256_add(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
        let mut result = [0u64; 4];
        let mut carry = 0u64;
        for i in 0..4 {
            let (r, c1) = a[i].overflowing_add(b[i]);
            let (r2, c2) = r.overflowing_add(carry);
            result[i] = r2;
            carry = (c1 as u64) + (c2 as u64);
        }
        result
    }
}
