use std::marker::PhantomData;

use crate::fields::field_params::FieldParams;
use crate::groups::affine_element::AffineElement;
use crate::groups::curve_params::{BaseField, CurveParams, ScalarField};
use crate::groups::wnaf;

/// An elliptic curve point in Jacobian projective coordinates (X : Y : Z).
///
/// Represents the affine point (X/Z^2, Y/Z^3).
/// Point at infinity uses the same encoding as AffineElement (MSB or modulus).
pub struct Element<C: CurveParams> {
    pub x: BaseField<C>,
    pub y: BaseField<C>,
    pub z: BaseField<C>,
    _phantom: PhantomData<C>,
}

impl<C: CurveParams> Clone for Element<C> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            x: self.x,
            y: self.y,
            z: self.z,
            _phantom: PhantomData,
        }
    }
}

impl<C: CurveParams> Copy for Element<C> {}

impl<C: CurveParams> std::fmt::Debug for Element<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_point_at_infinity() {
            write!(f, "Element(infinity)")
        } else {
            write!(f, "Element({:?}, {:?}, {:?})", self.x, self.y, self.z)
        }
    }
}

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

impl<C: CurveParams> Element<C> {
    /// Construct from (x, y, z) coordinates.
    #[inline]
    pub fn new(x: BaseField<C>, y: BaseField<C>, z: BaseField<C>) -> Self {
        Self {
            x,
            y,
            z,
            _phantom: PhantomData,
        }
    }

    /// Construct from an affine point (sets z = 1).
    #[inline]
    pub fn from_affine(affine: &AffineElement<C>) -> Self {
        Self::new(affine.x, affine.y, BaseField::<C>::one())
    }

    /// The generator point in projective form.
    #[inline]
    pub fn one() -> Self {
        Self::new(C::generator_x(), C::generator_y(), BaseField::<C>::one())
    }

    /// A random curve point, obtained by multiplying the generator by a random scalar.
    pub fn random_element() -> Self {
        let scalar = ScalarField::<C>::random_element();
        Self::one().mul_without_endomorphism(&scalar)
    }

    /// The point at infinity.
    #[inline]
    pub fn infinity() -> Self {
        let mut result = Self::new(
            BaseField::<C>::zero(),
            BaseField::<C>::zero(),
            BaseField::<C>::zero(),
        );
        result.self_set_infinity();
        result
    }
}

// ---------------------------------------------------------------------------
// Infinity flag
// ---------------------------------------------------------------------------

impl<C: CurveParams> Element<C> {
    /// Check if this is the point at infinity.
    #[inline]
    pub fn is_point_at_infinity(&self) -> bool {
        if C::BaseFieldParams::MODULUS_IS_BIG {
            let m = C::BaseFieldParams::MODULUS;
            ((self.x.data[0] ^ m[0])
                | (self.x.data[1] ^ m[1])
                | (self.x.data[2] ^ m[2])
                | (self.x.data[3] ^ m[3]))
                == 0
        } else {
            self.x.is_msb_set()
        }
    }

    /// Set this point to the point at infinity (in-place).
    #[inline]
    pub fn self_set_infinity(&mut self) {
        if C::BaseFieldParams::MODULUS_IS_BIG {
            let m = C::BaseFieldParams::MODULUS;
            self.x.data[0] = m[0];
            self.x.data[1] = m[1];
            self.x.data[2] = m[2];
            self.x.data[3] = m[3];
        } else {
            self.x = BaseField::<C>::zero();
            self.y = BaseField::<C>::zero();
            self.z = BaseField::<C>::zero();
            self.x.self_set_msb();
        }
    }

    /// Return a copy with infinity flag set.
    #[inline]
    pub fn set_infinity(&self) -> Self {
        let mut result = *self;
        result.self_set_infinity();
        result
    }
}

// ---------------------------------------------------------------------------
// Conversion
// ---------------------------------------------------------------------------

impl<C: CurveParams> Element<C> {
    /// Convert to affine coordinates via z-inverse.
    pub fn to_affine(&self) -> AffineElement<C> {
        if self.is_point_at_infinity() {
            return AffineElement::infinity();
        }
        let z_inv = self.z.invert();
        let zz_inv = z_inv.sqr();
        let zzz_inv = zz_inv * z_inv;
        AffineElement::new(self.x * zz_inv, self.y * zzz_inv)
    }

    /// Normalize: convert to affine and back to projective (z = 1).
    pub fn normalize(&self) -> Self {
        let affine = self.to_affine();
        Self::from_affine(&affine)
    }
}

// ---------------------------------------------------------------------------
// On-curve check
// ---------------------------------------------------------------------------

impl<C: CurveParams> Element<C> {
    /// Check if the point lies on the curve: y^2 == x^3 + a*x*z^4 + b*z^6.
    pub fn on_curve(&self) -> bool {
        if self.is_point_at_infinity() {
            return true;
        }
        if self.z.is_zero() {
            return false;
        }
        let zz = self.z.sqr();
        let zzzz = zz.sqr();
        let mut bz_6 = zzzz * zz * C::coeff_b();
        if C::HAS_A {
            bz_6 = bz_6 + (self.x * C::coeff_a()) * zzzz;
        }
        let xxx = self.x.sqr() * self.x + bz_6;
        let yy = self.y.sqr();
        xxx == yy
    }
}

// ---------------------------------------------------------------------------
// Point doubling
// ---------------------------------------------------------------------------

impl<C: CurveParams> Element<C> {
    /// Double this point in-place. Matches C++ `self_dbl()` from element_impl.hpp.
    pub fn self_dbl(&mut self) {
        // Check for infinity
        if C::BaseFieldParams::MODULUS_IS_BIG {
            if self.is_point_at_infinity() {
                return;
            }
        } else {
            if self.x.is_msb_set() {
                return;
            }
        }

        // T0 = x^2
        let t0 = self.x.sqr();
        // T1 = y^2
        let t1 = self.y.sqr();
        // T2 = y^4
        let t2 = t1.sqr();
        // T1 = (y^2 + x)^2
        let t1 = (t1 + self.x).sqr();
        // T3 = x^2 + y^4
        let t3 = t0 + t2;
        // T1 = 2*x*y^2  (= (y^2+x)^2 - x^2 - y^4)
        let t1 = t1 - t3;
        // S = 4*x*y^2
        let t1 = t1 + t1;
        // M = 3*x^2 (+ a*z^4 if has_a)
        let mut t3 = t0 + t0;
        t3 = t3 + t0;
        if C::HAS_A {
            t3 = t3 + (C::coeff_a() * self.z.sqr().sqr());
        }
        // z3 = 2*y*z (must use old y before it's overwritten)
        let new_z = (self.z + self.z) * self.y;
        // 2S
        let two_s = t1 + t1;
        // x3 = M^2 - 2S
        let new_x = t3.sqr() - two_s;
        // 8*y^4
        let t2 = t2 + t2;
        let t2 = t2 + t2;
        let t2 = t2 + t2;
        // y3 = M*(S - x3) - 8*y^4
        let new_y = t3 * (t1 - new_x) - t2;

        self.x = new_x;
        self.y = new_y;
        self.z = new_z;
    }

    /// Return the double of this point.
    #[inline]
    pub fn dbl(&self) -> Self {
        let mut result = *self;
        result.self_dbl();
        result
    }
}

// ---------------------------------------------------------------------------
// Mixed addition (projective += affine)
// ---------------------------------------------------------------------------

impl<C: CurveParams> Element<C> {
    /// Add an affine point to this projective point in-place.
    /// Matches C++ `operator+=(const affine_element&)` from element_impl.hpp.
    pub fn add_assign_affine(&mut self, other: &AffineElement<C>) {
        // Handle infinity edge cases
        if C::BaseFieldParams::MODULUS_IS_BIG {
            if self.is_point_at_infinity() {
                *self = Self::from_affine(other);
                return;
            }
        } else {
            let edge_case = self.x.is_msb_set() || other.x.is_msb_set();
            if edge_case {
                if self.x.is_msb_set() {
                    *self = Self::from_affine(other);
                }
                return;
            }
        }

        // T0 = z1^2
        let t0 = self.z.sqr();
        // T1 = x2*z1^2 - x1 (= H)
        let t1 = other.x * t0 - self.x;
        // T2 = z1^3 * y2 - y1
        let t2 = self.z * t0 * other.y - self.y;

        // Edge case: H == 0 means same x-coordinate
        if t1.is_zero() {
            if t2.is_zero() {
                // Same point: double
                self.self_dbl();
                return;
            }
            // Inverse points: result is infinity
            self.self_set_infinity();
            return;
        }

        // R = 2*(z1^3*y2 - y1)
        let t2 = t2 + t2;
        // z3 = z1 + H (will be squared later)
        self.z = self.z + t1;
        // HH = H^2
        let t3 = t1.sqr();
        // z1^2 + HH (for subtracting from z3^2)
        let t0 = t0 + t3;
        // z3 = (z1 + H)^2 - z1^2 - HH = 2*z1*H
        self.z = self.z.sqr();
        self.z = self.z - t0;
        // 4*HH
        let t3 = t3 + t3;
        let t3 = t3 + t3;
        // 4*HHH
        let t1 = t1 * t3;
        // 4*HH*x1
        let t3 = t3 * self.x;
        // 2*(4*HH*x1) + 4*HHH
        let t0 = t3 + t3;
        let t0 = t0 + t1;
        // x3 = R^2 - (8*HH*x1 + 4*HHH)
        self.x = t2.sqr();
        self.x = self.x - t0;
        // 4*HH*x1 - x3
        let t3 = t3 - self.x;
        // 2*y1 * 4*HHH  (uses old y)
        let t1 = t1 * self.y;
        let t1 = t1 + t1;
        // y3 = R*(4*HH*x1 - x3) - 2*y1*4*HHH
        let t3 = t3 * t2;
        self.y = t3 - t1;
    }

    /// Subtract an affine point from this projective point in-place.
    #[inline]
    pub fn sub_assign_affine(&mut self, other: &AffineElement<C>) {
        let neg_other = AffineElement::new(other.x, -other.y);
        self.add_assign_affine(&neg_other);
    }
}

// ---------------------------------------------------------------------------
// Full projective addition (projective += projective)
// ---------------------------------------------------------------------------

impl<C: CurveParams> Element<C> {
    /// Add another projective point to this one in-place.
    /// Matches C++ `operator+=(const element&)` from element_impl.hpp.
    pub fn add_assign_element(&mut self, other: &Self) {
        // Handle infinity edge cases
        if C::BaseFieldParams::MODULUS_IS_BIG {
            let p1_zero = self.is_point_at_infinity();
            let p2_zero = other.is_point_at_infinity();
            if p1_zero || p2_zero {
                if p1_zero && !p2_zero {
                    *self = *other;
                    return;
                }
                if p2_zero && !p1_zero {
                    return;
                }
                self.self_set_infinity();
                return;
            }
        } else {
            let p1_zero = self.x.is_msb_set();
            let p2_zero = other.x.is_msb_set();
            if p1_zero || p2_zero {
                if p1_zero && !p2_zero {
                    *self = *other;
                    return;
                }
                if p2_zero && !p1_zero {
                    return;
                }
                self.self_set_infinity();
                return;
            }
        }

        let z1z1 = self.z.sqr();
        let z2z2 = other.z.sqr();
        let mut s2 = z1z1 * self.z;
        let u2 = z1z1 * other.x;
        s2 = s2 * other.y;
        let u1 = z2z2 * self.x;
        let mut s1 = z2z2 * other.z;
        s1 = s1 * self.y;

        let f = s2 - s1;
        let h = u2 - u1;

        // Edge case: same x-coordinate in both points
        if h.is_zero() {
            if f.is_zero() {
                // Same point: double
                self.self_dbl();
                return;
            }
            // Inverse points: infinity
            self.self_set_infinity();
            return;
        }

        let f = f + f;
        let mut i = h + h;
        i = i.sqr();
        let j = h * i;
        let u1 = u1 * i;
        let u2_temp = u1 + u1;
        let u2_temp = u2_temp + j;

        self.x = f.sqr();
        self.x = self.x - u2_temp;

        let mut j = j * s1;
        j = j + j;

        self.y = u1 - self.x;
        self.y = self.y * f;
        self.y = self.y - j;

        self.z = self.z + other.z;
        let z1z1_plus_z2z2 = z1z1 + z2z2;
        self.z = self.z.sqr();
        self.z = self.z - z1z1_plus_z2z2;
        self.z = self.z * h;
    }

    /// Subtract another projective point from this one in-place.
    #[inline]
    pub fn sub_assign_element(&mut self, other: &Self) {
        let neg_other = Self::new(other.x, -other.y, other.z);
        self.add_assign_element(&neg_other);
    }
}

// ---------------------------------------------------------------------------
// Scalar multiplication
// ---------------------------------------------------------------------------

impl<C: CurveParams> Element<C> {
    /// Basic double-and-add scalar multiplication (no endomorphism).
    /// Matches C++ `mul_without_endomorphism`.
    pub fn mul_without_endomorphism(&self, scalar: &ScalarField<C>) -> Self {
        let converted = scalar.from_montgomery_form();

        // Check for zero scalar
        if converted.data[0] == 0
            && converted.data[1] == 0
            && converted.data[2] == 0
            && converted.data[3] == 0
        {
            return Self::infinity();
        }

        // Find MSB
        let mut msb = 0u64;
        for i in (0..4).rev() {
            if converted.data[i] != 0 {
                msb = (i as u64) * 64 + (63 - converted.data[i].leading_zeros() as u64);
                break;
            }
        }

        let mut accumulator = *self;
        // Iterate from msb-1 down to 0
        let mut i = msb.wrapping_sub(1);
        while i < msb {
            accumulator.self_dbl();
            let limb_idx = (i / 64) as usize;
            let bit_idx = i % 64;
            if (converted.data[limb_idx] >> bit_idx) & 1 == 1 {
                accumulator.add_assign_element(self);
            }
            i = i.wrapping_sub(1);
        }
        accumulator
    }

    /// Convenience: scalar multiply using raw integer limbs (non-Montgomery).
    /// Constructs a scalar field element from the limbs (converting to Montgomery form)
    /// and then calls mul_without_endomorphism.
    pub fn scalar_mul(&self, raw_scalar: &[u64; 4]) -> Self {
        let scalar = ScalarField::<C>::from_limbs(*raw_scalar);
        self.mul_without_endomorphism(&scalar)
    }

    /// Scalar multiplication using GLV endomorphism.
    ///
    /// Splits the scalar into two ~128-bit half-scalars via the endomorphism,
    /// then uses a 4-bit WNAF with an 8-entry lookup table for efficient
    /// double scalar multiplication (k1*P + k2*endo(P)).
    ///
    /// Matches C++ `mul_with_endomorphism` from element_impl.hpp.
    pub fn mul_with_endomorphism(&self, scalar: &ScalarField<C>) -> Self {
        if self.is_point_at_infinity() {
            return Self::infinity();
        }

        const NUM_ROUNDS: usize = 32;
        const NUM_WNAF_BITS: usize = 4;
        const LOOKUP_SIZE: usize = 8;

        let converted_scalar = scalar.from_montgomery_form();
        if converted_scalar.is_zero() {
            return Self::infinity();
        }

        // Build 8-entry lookup table: table[i] = (2i+1) * P
        let mut lookup_table = [*self; LOOKUP_SIZE];
        let d2 = self.dbl();
        for i in 1..LOOKUP_SIZE {
            lookup_table[i] = lookup_table[i - 1];
            lookup_table[i].add_assign_element(&d2);
        }

        // Split scalar into two half-scalars.
        // C++ pair-returning version only works for 254-bit moduli (static_assert).
        // For those curves, the half-scalars always fit in 128 bits (lower 2 limbs).
        let (k1_full, k2_full) = converted_scalar.split_into_endomorphism_scalars();
        let k1 = [k1_full.data[0], k1_full.data[1]];
        let k2 = [k2_full.data[0], k2_full.data[1]];

        // Compute interleaved WNAF for both half-scalars
        let mut wnaf_table = [0u64; NUM_ROUNDS * 2];
        let mut skew = false;
        let mut endo_skew = false;
        wnaf::fixed_wnaf(&k1, &mut wnaf_table, &mut skew, 0, 2, NUM_WNAF_BITS);
        wnaf::fixed_wnaf(&k2, &mut wnaf_table[1..], &mut endo_skew, 0, 2, NUM_WNAF_BITS);

        let mut accumulator = Self::infinity();
        let beta = BaseField::<C>::cube_root_of_unity();

        // C++ element_impl.hpp:685-702
        for i in 0..(NUM_ROUNDS * 2) {
            let wnaf_entry = wnaf_table[i];
            let index = (wnaf_entry & 0x0fffffff) as usize;
            let sign = ((wnaf_entry >> 31) & 1) != 0;
            let is_odd = (i & 1) == 1;

            let mut to_add = lookup_table[index];
            to_add.y.self_conditional_negate(sign ^ is_odd);
            if is_odd {
                to_add.x = to_add.x * beta;
            }
            accumulator.add_assign_element(&to_add);

            if i != (2 * NUM_ROUNDS - 1) && is_odd {
                for _ in 0..4 {
                    accumulator.self_dbl();
                }
            }
        }

        // Skew correction (C++ element_impl.hpp:704-709)
        if skew {
            let neg = Self::new(lookup_table[0].x, lookup_table[0].y.negate(), lookup_table[0].z);
            accumulator.add_assign_element(&neg);
        }
        if endo_skew {
            let endo_point = Self::new(
                lookup_table[0].x * beta,
                lookup_table[0].y,
                lookup_table[0].z,
            );
            accumulator.add_assign_element(&endo_point);
        }

        accumulator
    }

    /// Scalar multiplication: dispatches to endomorphism or basic path.
    ///
    /// C++ mul_with_endomorphism only compiles for 254-bit scalar fields
    /// (the pair-returning split has a static_assert). For big-modulus curves
    /// like secp256k1, always use the basic path.
    pub fn mul(&self, scalar: &ScalarField<C>) -> Self {
        if C::USE_ENDOMORPHISM && !<C::ScalarFieldParams as FieldParams>::MODULUS_IS_BIG {
            self.mul_with_endomorphism(scalar)
        } else {
            self.mul_without_endomorphism(scalar)
        }
    }

    /// Compute a*self + b*other (double scalar multiplication).
    ///
    /// Used by ECDSA verify and Schnorr verify. Naive approach using two
    /// separate scalar multiplications and addition. (Shamir's trick / interleaved
    /// NAF optimization can be added later when MSM is ported.)
    pub fn double_scalar_mul(
        &self,
        a: &ScalarField<C>,
        other: &Self,
        b: &ScalarField<C>,
    ) -> Self {
        let term1 = self.mul_without_endomorphism(a);
        let term2 = other.mul_without_endomorphism(b);
        term1 + term2
    }
}

// ---------------------------------------------------------------------------
// Batch operations
// ---------------------------------------------------------------------------

impl<C: CurveParams> Element<C> {
    /// Batch normalize: convert N Jacobian points to affine (z=1) using a single inversion.
    ///
    /// Uses Montgomery's trick: accumulates the product of z-coordinates forward,
    /// inverts once, then walks backward to recover individual z-inverses.
    /// Infinity points are skipped in both passes.
    ///
    /// Matches C++ `batch_normalize` from element_impl.hpp.
    pub fn batch_normalize(elements: &mut [Self]) {
        let num_elements = elements.len();
        if num_elements == 0 {
            return;
        }

        let mut temporaries = Vec::with_capacity(num_elements);
        let mut accumulator = BaseField::<C>::one();

        // Forward pass: accumulate z-coordinate products
        for i in 0..num_elements {
            temporaries.push(accumulator);
            if !elements[i].is_point_at_infinity() {
                accumulator = accumulator * elements[i].z;
            }
        }

        // Single inversion of accumulated product
        accumulator = accumulator.invert();

        // Backward pass: compute individual z-inverses and convert to affine
        for i in (0..num_elements).rev() {
            if !elements[i].is_point_at_infinity() {
                let z_inv = accumulator * temporaries[i];
                let zz_inv = z_inv.sqr();
                elements[i].x = elements[i].x * zz_inv;
                elements[i].y = elements[i].y * (zz_inv * z_inv);
                accumulator = accumulator * elements[i].z;
            }
            elements[i].z = BaseField::<C>::one();
        }
    }
}

// ---------------------------------------------------------------------------
// Operator impls
// ---------------------------------------------------------------------------

impl<C: CurveParams> std::ops::Add<AffineElement<C>> for Element<C> {
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: AffineElement<C>) -> Self {
        self.add_assign_affine(&rhs);
        self
    }
}

impl<C: CurveParams> std::ops::AddAssign<AffineElement<C>> for Element<C> {
    #[inline]
    fn add_assign(&mut self, rhs: AffineElement<C>) {
        self.add_assign_affine(&rhs);
    }
}

impl<C: CurveParams> std::ops::Add for Element<C> {
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: Self) -> Self {
        self.add_assign_element(&rhs);
        self
    }
}

impl<C: CurveParams> std::ops::AddAssign for Element<C> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.add_assign_element(&rhs);
    }
}

impl<C: CurveParams> std::ops::Sub<AffineElement<C>> for Element<C> {
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: AffineElement<C>) -> Self {
        self.sub_assign_affine(&rhs);
        self
    }
}

impl<C: CurveParams> std::ops::SubAssign<AffineElement<C>> for Element<C> {
    #[inline]
    fn sub_assign(&mut self, rhs: AffineElement<C>) {
        self.sub_assign_affine(&rhs);
    }
}

impl<C: CurveParams> std::ops::Sub for Element<C> {
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: Self) -> Self {
        self.sub_assign_element(&rhs);
        self
    }
}

impl<C: CurveParams> std::ops::SubAssign for Element<C> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.sub_assign_element(&rhs);
    }
}

impl<C: CurveParams> std::ops::Neg for Element<C> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(self.x, -self.y, self.z)
    }
}

impl<C: CurveParams> PartialEq for Element<C> {
    fn eq(&self, other: &Self) -> bool {
        if !self.on_curve() || !other.on_curve() {
            return false;
        }
        let am_inf = self.is_point_at_infinity();
        let is_inf = other.is_point_at_infinity();
        let both_inf = am_inf && is_inf;
        if !both_inf && (am_inf || is_inf) {
            return false;
        }
        let lhs_zz = self.z.sqr();
        let lhs_zzz = lhs_zz * self.z;
        let rhs_zz = other.z.sqr();
        let rhs_zzz = rhs_zz * other.z;

        let lhs_x = self.x * rhs_zz;
        let lhs_y = self.y * rhs_zzz;
        let rhs_x = other.x * lhs_zz;
        let rhs_y = other.y * lhs_zzz;
        both_inf || (lhs_x == rhs_x && lhs_y == rhs_y)
    }
}

impl<C: CurveParams> Eq for Element<C> {}
