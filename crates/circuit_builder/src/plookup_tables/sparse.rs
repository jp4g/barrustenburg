//! Sparse representation lookup tables for plookup.
//!
//! C++ source: barretenberg/cpp/src/barretenberg/stdlib_circuit_builders/plookup_tables/sparse.hpp
//!
//! Provides functions to generate lookup tables that map binary values into a
//! sparse representation in a given base, optionally with 32-bit rotation.
//! Also provides normalization tables that map sparse-form values back to
//! binary using a caller-supplied normalization table.

use bbrs_ecc::curves::bn254::Fr;
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::curves::bn254::Bn254FrParams;
use bbrs_numeric::bitop::pow64;

use super::types::{BasicTable, BasicTableId, GetValuesFromKey};

// ---------------------------------------------------------------------------
// Helper: sparse form conversion
// ---------------------------------------------------------------------------

/// Map a binary integer into sparse form in the given base.
///
/// Each bit of `input` contributes `base^bit_position` to the output.
/// This mirrors C++ `numeric::map_into_sparse_form<base>(input)`.
///
/// The result can exceed 64 bits for large bases/inputs, so we return [u64; 4]
/// limbs (little-endian) suitable for constructing an Fr via `from_limbs`.
fn map_into_sparse_form(base: u64, input: u64) -> [u64; 4] {
    // Precompute base powers as u128 to handle overflow for up to 32 bits.
    // For base=28, 28^31 overflows u64 but we accumulate in u128 pairs then
    // split into 4 x u64 limbs.
    //
    // Actually, we accumulate into a simple big-integer represented as [u64; 4].
    let mut out = [0u64; 4];

    // Precompute base^i as [u64; 4] for i in 0..32
    let mut base_powers = [[0u64; 4]; 32];
    base_powers[0] = [1, 0, 0, 0];
    for i in 1..32 {
        base_powers[i] = u256_mul_u64(base_powers[i - 1], base);
    }

    for i in 0..32 {
        let sparse_bit = (input >> i) & 1;
        if sparse_bit != 0 {
            out = u256_add(out, base_powers[i]);
        }
    }
    out
}

/// Helper: multiply a little-endian [u64; 4] value by a u64 scalar.
#[inline]
fn u256_mul_u64(a: [u64; 4], b: u64) -> [u64; 4] {
    let mut result = [0u64; 4];
    let mut carry = 0u128;
    for i in 0..4 {
        let prod = (a[i] as u128) * (b as u128) + carry;
        result[i] = prod as u64;
        carry = prod >> 64;
    }
    // carry overflow is ignored (assumed to fit in 256 bits)
    result
}

/// Helper: add two little-endian [u64; 4] values.
#[inline]
fn u256_add(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
    let mut result = [0u64; 4];
    let mut carry = 0u64;
    for i in 0..4 {
        let (sum1, c1) = a[i].overflowing_add(b[i]);
        let (sum2, c2) = sum1.overflowing_add(carry);
        result[i] = sum2;
        carry = (c1 as u64) + (c2 as u64);
    }
    result
}

/// Construct an Fr from the sparse-form [u64; 4] limbs (non-Montgomery).
#[inline]
fn fr_from_sparse(limbs: [u64; 4]) -> Fr {
    Field::<Bn254FrParams>::from_limbs(limbs)
}

// ---------------------------------------------------------------------------
// Helper: 32-bit rotation
// ---------------------------------------------------------------------------

/// Rotate a 32-bit value right by `rotation` bits (with wrap-around).
///
/// Mirrors C++ `numeric::rotate32(value, rotation)`.
#[inline]
fn rotate32(value: u32, rotation: u32) -> u32 {
    if rotation == 0 {
        value
    } else {
        (value >> rotation) | (value << (32 - rotation))
    }
}

// ---------------------------------------------------------------------------
// Helper: sparse_int (mirrors C++ numeric::sparse_int<base, num_bits>)
// ---------------------------------------------------------------------------

/// A sparse integer with `num_bits` limbs in the given base.
///
/// Each limb stores a digit in [0, base). The sparse value is the polynomial
/// evaluation: limbs[0] + limbs[1]*base + limbs[2]*base^2 + ...
struct SparseInt {
    limbs: Vec<u64>,
    base: u64,
}

impl SparseInt {
    /// Construct from a binary integer, decomposing each bit into a limb.
    fn new(base: u64, num_bits: usize, input: u64) -> Self {
        let mut limbs = vec![0u64; num_bits];
        for i in 0..num_bits {
            limbs[i] = (input >> i) & 1;
        }
        SparseInt { limbs, base }
    }

    /// Add another sparse_int in-place with carry propagation.
    fn add_assign(&mut self, other: &SparseInt) {
        let num_bits = self.limbs.len();
        for i in 0..(num_bits - 1) {
            self.limbs[i] += other.limbs[i];
            if self.limbs[i] >= self.base {
                self.limbs[i] -= self.base;
                self.limbs[i + 1] += 1;
            }
        }
        self.limbs[num_bits - 1] += other.limbs[num_bits - 1];
        self.limbs[num_bits - 1] %= self.base;
    }

    /// Compute the sparse polynomial evaluation as a u64.
    ///
    /// result = limbs[0] + limbs[1]*base + limbs[2]*base^2 + ...
    /// Evaluates from the top limb down via Horner's method.
    fn get_sparse_value(&self) -> u64 {
        let num_bits = self.limbs.len();
        let mut result = 0u64;
        // Iterate from most significant limb to least, matching C++:
        //   for (size_t i = num_bits - 1; i < num_bits; --i)
        for i in (0..num_bits).rev() {
            result = result.wrapping_mul(self.base);
            result = result.wrapping_add(self.limbs[i]);
        }
        result
    }
}

// ---------------------------------------------------------------------------
// get_sparse_table_with_rotation_values: runtime-parameterized version
// ---------------------------------------------------------------------------

/// Compute sparse-form values for a key, with optional 32-bit rotation.
///
/// Returns `[sparse(key[0]), sparse(rotate32(key[0], num_rotated_bits))]`.
/// When `num_rotated_bits == 0`, both outputs are the same.
///
/// Mirrors C++ `get_sparse_table_with_rotation_values<base, num_rotated_bits>`.
pub fn get_sparse_table_with_rotation_values(
    base: u64,
    num_rotated_bits: u64,
    key: [u64; 2],
) -> [Fr; 2] {
    let t0 = fr_from_sparse(map_into_sparse_form(base, key[0]));
    let t1 = if num_rotated_bits > 0 {
        fr_from_sparse(map_into_sparse_form(
            base,
            rotate32(key[0] as u32, num_rotated_bits as u32) as u64,
        ))
    } else {
        t0
    };
    [t0, t1]
}

// ---------------------------------------------------------------------------
// Concrete get_values_from_key instantiations
// ---------------------------------------------------------------------------
//
// Since BasicTable.get_values_from_key is `fn([u64; 2]) -> [Fr; 2]` (no
// captured state), we provide concrete functions for each (base, rot)
// combination used in the codebase.

/// base=9, rot=0 (AES_SPARSE_MAP)
pub fn get_sparse_values_base9_rot0(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(9, 0, key)
}

/// base=16, rot=0 (SHA256_WITNESS_SLICE_3, SHA256_BASE16)
pub fn get_sparse_values_base16_rot0(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(16, 0, key)
}

/// base=16, rot=1 (SHA256_WITNESS_SLICE_14_ROTATE_1)
pub fn get_sparse_values_base16_rot1(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(16, 1, key)
}

/// base=16, rot=2 (SHA256_BASE16_ROTATE2)
pub fn get_sparse_values_base16_rot2(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(16, 2, key)
}

/// base=16, rot=4 (SHA256_WITNESS_SLICE_7_ROTATE_4)
pub fn get_sparse_values_base16_rot4(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(16, 4, key)
}

/// base=16, rot=6 (SHA256_BASE16_ROTATE6)
pub fn get_sparse_values_base16_rot6(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(16, 6, key)
}

/// base=16, rot=7 (SHA256_WITNESS_SLICE_8_ROTATE_7, SHA256_BASE16_ROTATE7)
pub fn get_sparse_values_base16_rot7(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(16, 7, key)
}

/// base=16, rot=8 (SHA256_BASE16_ROTATE8)
pub fn get_sparse_values_base16_rot8(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(16, 8, key)
}

/// base=28, rot=0 (SHA256_BASE28)
pub fn get_sparse_values_base28_rot0(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(28, 0, key)
}

/// base=28, rot=3 (SHA256_BASE28_ROTATE3)
pub fn get_sparse_values_base28_rot3(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(28, 3, key)
}

/// base=28, rot=6 (SHA256_BASE28_ROTATE6)
pub fn get_sparse_values_base28_rot6(key: [u64; 2]) -> [Fr; 2] {
    get_sparse_table_with_rotation_values(28, 6, key)
}

// ---------------------------------------------------------------------------
// generate_sparse_table_with_rotation
// ---------------------------------------------------------------------------

/// Select the correct concrete `get_values_from_key` function pointer for the
/// given (base, num_rotated_bits) combination.
///
/// Panics if the combination is not one of the known instantiations.
fn select_sparse_get_values(base: u64, num_rotated_bits: usize) -> GetValuesFromKey {
    match (base, num_rotated_bits) {
        (9, 0) => get_sparse_values_base9_rot0,
        (16, 0) => get_sparse_values_base16_rot0,
        (16, 1) => get_sparse_values_base16_rot1,
        (16, 2) => get_sparse_values_base16_rot2,
        (16, 4) => get_sparse_values_base16_rot4,
        (16, 6) => get_sparse_values_base16_rot6,
        (16, 7) => get_sparse_values_base16_rot7,
        (16, 8) => get_sparse_values_base16_rot8,
        (28, 0) => get_sparse_values_base28_rot0,
        (28, 3) => get_sparse_values_base28_rot3,
        (28, 6) => get_sparse_values_base28_rot6,
        _ => panic!(
            "No concrete get_values_from_key for sparse table with base={}, rot={}",
            base, num_rotated_bits
        ),
    }
}

/// Generate a sparse representation table with optional 32-bit rotation.
///
/// The table maps each integer in `[0, 2^bits_per_slice)` to its sparse-form
/// representation in the given base. If `num_rotated_bits > 0`, column 3
/// stores the sparse form of the rotated value.
///
/// Mirrors C++ `generate_sparse_table_with_rotation<base, bits_per_slice, num_rotated_bits>`.
pub fn generate_sparse_table_with_rotation(
    base: u64,
    bits_per_slice: usize,
    num_rotated_bits: usize,
    id: BasicTableId,
    table_index: usize,
) -> BasicTable {
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    let table_size = 1u64 << bits_per_slice;

    for i in 0..table_size {
        let source = i;
        let target = fr_from_sparse(map_into_sparse_form(base, source));

        table.column_1.push(Fr::from(source));
        table.column_2.push(target);

        if num_rotated_bits > 0 {
            let rotated = fr_from_sparse(map_into_sparse_form(
                base,
                rotate32(source as u32, num_rotated_bits as u32) as u64,
            ));
            table.column_3.push(rotated);
        } else {
            table.column_3.push(target);
        }
    }

    table.get_values_from_key = select_sparse_get_values(base, num_rotated_bits);

    // sparse_step_size = base^bits_per_slice
    // Compute as u128 to avoid overflow, then store as Fr.
    let mut sparse_step_size = [1u64, 0u64, 0u64, 0u64];
    for _ in 0..bits_per_slice {
        sparse_step_size = u256_mul_u64(sparse_step_size, base);
    }
    let sparse_step_fr = Field::<Bn254FrParams>::from_limbs(sparse_step_size);

    // C++ uses a hardcoded column_1_step_size of (1 << 11) regardless of bits_per_slice.
    // This matches the multi-table accumulator convention.
    table.column_1_step_size = Fr::from(1u64 << 11);
    table.column_2_step_size = sparse_step_fr;
    table.column_3_step_size = sparse_step_fr;

    table
}

// ---------------------------------------------------------------------------
// get_sparse_normalization_values
// ---------------------------------------------------------------------------

/// Normalize a sparse-form value by decomposing it into base-ary digits and
/// looking each digit up in the caller-supplied normalization table, then
/// accumulating the looked-up bits.
///
/// Returns `[accumulated_result, 0]`.
///
/// Mirrors C++ `get_sparse_normalization_values<base, base_table>`.
pub fn get_sparse_normalization_values(
    base: u64,
    normalization_table: &[u64],
    key: [u64; 2],
) -> [Fr; 2] {
    let mut accumulator = 0u64;
    let mut input = key[0];
    let mut count = 0u64;
    while input > 0 {
        let slice = input % base;
        let bit = normalization_table[slice as usize];
        accumulator += bit << count;
        input -= slice;
        input /= base;
        count += 1;
    }
    [Fr::from(accumulator), Fr::zero()]
}

// ---------------------------------------------------------------------------
// generate_sparse_normalization_table
// ---------------------------------------------------------------------------

/// Generate a normalization table that maps sparse-form values back to binary
/// using the given normalization table.
///
/// The table enumerates all `base^num_bits` sparse values, and for each one
/// looks up each digit in `normalization_table` to produce the normalized
/// output.
///
/// Mirrors C++ `generate_sparse_normalization_table<base, num_bits, base_table>`.
pub fn generate_sparse_normalization_table(
    base: u64,
    num_bits: usize,
    normalization_table: &[u64],
    id: BasicTableId,
    table_index: usize,
    get_values_from_key: GetValuesFromKey,
) -> BasicTable {
    let mut table = BasicTable::new();
    table.id = id;
    table.table_index = table_index;
    table.use_twin_keys = false;

    let table_size = pow64(base, num_bits as u64);

    let mut accumulator = SparseInt::new(base, num_bits, 0);
    let to_add = SparseInt::new(base, num_bits, 1);

    for _i in 0..table_size {
        let limbs = &accumulator.limbs;
        let mut key = 0u64;
        for j in 0..num_bits {
            let table_idx = limbs[j] as usize;
            key += normalization_table[table_idx] << (j as u64);
        }

        table
            .column_1
            .push(Fr::from(accumulator.get_sparse_value()));
        table.column_2.push(Fr::from(key));
        table.column_3.push(Fr::zero());

        accumulator.add_assign(&to_add);
    }

    table.get_values_from_key = get_values_from_key;

    table.column_1_step_size = Fr::from(table_size);
    table.column_2_step_size = Fr::from(1u64 << (num_bits as u64));
    table.column_3_step_size = Fr::zero();

    table
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_map_into_sparse_form_base16() {
        // 0b0000_0001 => 16^0 = 1
        let result = map_into_sparse_form(16, 1);
        assert_eq!(result, [1, 0, 0, 0]);

        // 0b0000_0011 => 16^0 + 16^1 = 1 + 16 = 17
        let result = map_into_sparse_form(16, 3);
        assert_eq!(result, [17, 0, 0, 0]);

        // 0b0000_0101 => 16^0 + 16^2 = 1 + 256 = 257
        let result = map_into_sparse_form(16, 5);
        assert_eq!(result, [257, 0, 0, 0]);
    }

    #[test]
    fn test_map_into_sparse_form_base9() {
        // 0b01 => 9^0 = 1
        let result = map_into_sparse_form(9, 1);
        assert_eq!(result, [1, 0, 0, 0]);

        // 0b11 => 9^0 + 9^1 = 1 + 9 = 10
        let result = map_into_sparse_form(9, 3);
        assert_eq!(result, [10, 0, 0, 0]);

        // 0xFF => sum of 9^i for i in 0..8 = (9^8 - 1) / 8
        let result = map_into_sparse_form(9, 0xFF);
        let mut expected = 0u64;
        for i in 0..8 {
            expected += pow64(9, i);
        }
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_rotate32() {
        // rotate32(0x80000001, 1) should give 0xC0000000
        assert_eq!(rotate32(0x80000001, 1), 0xC0000000);
        // rotate32(x, 0) == x
        assert_eq!(rotate32(0x12345678, 0), 0x12345678);
        // rotate32(x, 32) wraps: should panic or wrap. In C++ rotate32(x,32) is UB,
        // but we handle rotation=0 case above. For rotation < 32:
        assert_eq!(rotate32(1, 1), 0x80000000);
    }

    #[test]
    fn test_sparse_int_basic() {
        // SparseInt for base=16, num_bits=3, input=5 (binary 101)
        let si = SparseInt::new(16, 3, 5);
        assert_eq!(si.limbs, vec![1, 0, 1]); // bit0=1, bit1=0, bit2=1
        // sparse_value = 1 + 0*16 + 1*256 = 257
        assert_eq!(si.get_sparse_value(), 257);
    }

    #[test]
    fn test_sparse_int_addition() {
        let mut a = SparseInt::new(16, 3, 0);
        let one = SparseInt::new(16, 3, 1);
        // Adding 1 repeatedly increments the lowest digit.
        a.add_assign(&one); // limbs = [1, 0, 0], sparse = 1
        assert_eq!(a.get_sparse_value(), 1);
        a.add_assign(&one); // limbs = [2, 0, 0], sparse = 2
        assert_eq!(a.get_sparse_value(), 2);

        // Add until we hit base (16): should carry into next limb
        for _ in 0..14 {
            a.add_assign(&one);
        }
        // After 16 total additions: limbs = [0, 1, 0], sparse = 0 + 1*16 = 16
        assert_eq!(a.get_sparse_value(), 16);
    }

    #[test]
    fn test_generate_sparse_table_base16_no_rotation() {
        let table = generate_sparse_table_with_rotation(
            16,
            3,
            0,
            BasicTableId::SHA256_WITNESS_SLICE_3,
            0,
        );
        // Table should have 2^3 = 8 entries
        assert_eq!(table.column_1.len(), 8);
        assert_eq!(table.column_2.len(), 8);
        assert_eq!(table.column_3.len(), 8);

        // First entry: source=0, sparse=0
        assert_eq!(table.column_1[0], Fr::from(0u64));
        assert_eq!(table.column_2[0], fr_from_sparse(map_into_sparse_form(16, 0)));

        // Entry for 5: sparse(5) = 16^0 + 16^2 = 257
        assert_eq!(table.column_1[5], Fr::from(5u64));
        assert_eq!(table.column_2[5], Fr::from(257u64));

        // With rot=0, column_3 == column_2
        for i in 0..8 {
            assert_eq!(table.column_2[i], table.column_3[i]);
        }
    }

    #[test]
    fn test_generate_sparse_table_with_rotation() {
        let table = generate_sparse_table_with_rotation(
            16,
            3,
            1,
            BasicTableId::SHA256_WITNESS_SLICE_14_ROTATE_1,
            0,
        );
        assert_eq!(table.column_1.len(), 8);

        // For source=1: sparse(1)=1, rotate32(1,1)=0x80000000,
        // sparse(0x80000000) = 16^31
        let rotated_val = rotate32(1, 1);
        let expected_rot = fr_from_sparse(map_into_sparse_form(16, rotated_val as u64));
        assert_eq!(table.column_3[1], expected_rot);
    }

    #[test]
    fn test_normalization_values() {
        // Simple normalization: identity table for base 16
        // Each digit maps to itself mod 2 (XOR-style)
        let norm_table: Vec<u64> = (0..16).map(|x| x & 1).collect();
        let result = get_sparse_normalization_values(16, &norm_table, [17, 0]);
        // 17 in base 16 = 1*16 + 1 => digits are [1, 1]
        // norm_table[1] = 1 for both digits
        // accumulator = 1 << 0 + 1 << 1 = 3
        assert_eq!(result[0], Fr::from(3u64));
        assert_eq!(result[1], Fr::zero());
    }

    #[test]
    fn test_generate_normalization_table() {
        // Use a trivial normalization table for base=9, num_bits=2
        let norm_table: [u64; 9] = [0, 1, 0, 1, 0, 1, 0, 1, 0];

        // Need a concrete get_values_from_key; use a simple one
        fn test_get_values(key: [u64; 2]) -> [Fr; 2] {
            let norm: [u64; 9] = [0, 1, 0, 1, 0, 1, 0, 1, 0];
            get_sparse_normalization_values(9, &norm, key)
        }

        let table = generate_sparse_normalization_table(
            9,
            2,
            &norm_table,
            BasicTableId::AES_SPARSE_NORMALIZE,
            0,
            test_get_values,
        );

        // 9^2 = 81 entries
        assert_eq!(table.column_1.len(), 81);
        assert_eq!(table.column_2.len(), 81);

        // First entry: all limbs zero => key = 0
        assert_eq!(table.column_1[0], Fr::from(0u64));
        assert_eq!(table.column_2[0], Fr::from(0u64));
    }
}
