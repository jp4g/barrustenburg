use crate::curves::bn254::{Bn254FqParams, Bn254FrParams, Bn254G1Params, Fq, Fq2, Fr, G2AffineElement, G2Element};
use crate::curves::grumpkin::GrumpkinG1Params;
use crate::curves::secp256k1::{
    Secp256k1FqParams, Secp256k1FrParams, Secp256k1G1Params, Secp256k1Fr,
};
use crate::curves::secp256r1::{Secp256r1FqParams, Secp256r1FrParams, Secp256r1G1Params, Secp256r1Fr};
use crate::fields::field::Field;
use crate::fields::field_params::FieldParams;
use crate::groups::affine_element::AffineElement;
use crate::groups::element::Element;
use crate::scalar_multiplication;
use crate::batched_affine_addition;

// =========================================================================
// Field arithmetic tests
// =========================================================================

#[test]
fn bn254_fq_modulus_matches_known_value() {
    let m = Bn254FqParams::MODULUS;
    assert_eq!(m[0], 0x3C208C16D87CFD47);
    assert_eq!(m[1], 0x97816a916871ca8d);
    assert_eq!(m[2], 0xb85045b68181585d);
    assert_eq!(m[3], 0x30644e72e131a029);
}

#[test]
fn bn254_fr_modulus_matches_known_value() {
    let m = Bn254FrParams::MODULUS;
    assert_eq!(m[0], 0x43E1F593F0000001);
    assert_eq!(m[1], 0x2833E84879B97091);
    assert_eq!(m[2], 0xB85045B68181585D);
    assert_eq!(m[3], 0x30644E72E131A029);
}

#[test]
fn field_one_times_one_is_one() {
    let one = Fq::one();
    assert_eq!(one * one, one);
}

#[test]
fn field_one_times_one_is_one_fr() {
    let one = Fr::one();
    assert_eq!(one * one, one);
}

#[test]
fn field_zero_is_additive_identity() {
    let one = Fq::one();
    let zero = Fq::zero();
    assert_eq!(one + zero, one);
    assert_eq!(zero + one, one);
}

#[test]
fn field_add_small_values() {
    let a = Fq::from(3u64);
    let b = Fq::from(5u64);
    let c = Fq::from(8u64);
    assert_eq!(a + b, c);
}

#[test]
fn field_mul_small_values() {
    let a = Fq::from(3u64);
    let b = Fq::from(5u64);
    let c = Fq::from(15u64);
    assert_eq!(a * b, c);
}

#[test]
fn field_sub_small_values() {
    let a = Fq::from(10u64);
    let b = Fq::from(3u64);
    let c = Fq::from(7u64);
    assert_eq!(a - b, c);
}

#[test]
fn field_negate() {
    let a = Fq::from(5u64);
    let neg_a = -a;
    assert_eq!(a + neg_a, Fq::zero());
}

#[test]
fn field_mul_inverse() {
    let a = Fq::from(7u64);
    let a_inv = a.invert();
    assert_eq!(a * a_inv, Fq::one());
}

#[test]
fn field_mul_inverse_larger() {
    let a = Fq::from(123456789u64);
    let a_inv = a.invert();
    assert_eq!(a * a_inv, Fq::one());
}

#[test]
fn field_montgomery_roundtrip() {
    let raw = [42u64, 0, 0, 0];
    let f = Field::<Bn254FqParams>::from_raw(raw);
    let mont = f.to_montgomery_form();
    let back = mont.from_montgomery_form();
    assert_eq!(back.data, raw);
}

#[test]
fn field_from_creates_montgomery_form() {
    let a = Fq::from(1u64);
    let b = Fq::one();
    assert_eq!(a, b);
}

#[test]
fn field_sqr_equals_mul() {
    let a = Fq::from(42u64);
    assert_eq!(a.sqr(), a * a);
}

#[test]
fn field_pow_small() {
    let a = Fq::from(3u64);
    let a_cubed = a.pow(&[3, 0, 0, 0]);
    assert_eq!(a_cubed, Fq::from(27u64));
}

#[test]
fn field_sqrt_perfect_square() {
    let a = Fq::from(9u64);
    let (is_qr, root) = a.sqrt();
    assert!(is_qr);
    assert_eq!(root.sqr(), a);
}

#[test]
fn field_is_zero() {
    assert!(Fq::zero().is_zero());
    assert!(!Fq::one().is_zero());
}

// =========================================================================
// Big modulus field tests (secp256k1)
// =========================================================================

#[test]
fn secp256k1_fq_one_times_one() {
    type K1Fq = Field<Secp256k1FqParams>;
    let one = K1Fq::one();
    assert_eq!(one * one, one);
}

#[test]
fn secp256k1_fq_add_small() {
    type K1Fq = Field<Secp256k1FqParams>;
    let a = K1Fq::from(3u64);
    let b = K1Fq::from(5u64);
    assert_eq!(a + b, K1Fq::from(8u64));
}

#[test]
fn secp256k1_fq_mul_inverse() {
    type K1Fq = Field<Secp256k1FqParams>;
    let a = K1Fq::from(42u64);
    let a_inv = a.invert();
    assert_eq!(a * a_inv, K1Fq::one());
}

#[test]
fn secp256k1_fq_montgomery_roundtrip() {
    type K1Fq = Field<Secp256k1FqParams>;
    let raw = [12345u64, 0, 0, 0];
    let f = K1Fq::from_raw(raw);
    let mont = f.to_montgomery_form();
    let back = mont.from_montgomery_form();
    assert_eq!(back.data, raw);
}

// =========================================================================
// secp256r1 big modulus tests
// =========================================================================

#[test]
fn secp256r1_fq_one_times_one() {
    type R1Fq = Field<Secp256r1FqParams>;
    let one = R1Fq::one();
    assert_eq!(one * one, one);
}

#[test]
fn secp256r1_fq_mul_inverse() {
    type R1Fq = Field<Secp256r1FqParams>;
    let a = R1Fq::from(99u64);
    let a_inv = a.invert();
    assert_eq!(a * a_inv, R1Fq::one());
}

// =========================================================================
// Group tests — BN254
// =========================================================================

#[test]
fn bn254_generator_on_curve_affine() {
    let g = AffineElement::<Bn254G1Params>::one();
    assert!(g.on_curve(), "BN254 generator should be on the curve");
}

#[test]
fn bn254_generator_on_curve_projective() {
    let g = Element::<Bn254G1Params>::one();
    assert!(g.on_curve(), "BN254 generator should be on the curve (projective)");
}

#[test]
fn bn254_infinity_on_curve() {
    let inf = AffineElement::<Bn254G1Params>::infinity();
    assert!(inf.on_curve());
    assert!(inf.is_point_at_infinity());
}

#[test]
fn bn254_double_equals_add_self() {
    let g = Element::<Bn254G1Params>::one();
    let doubled = g.dbl();
    let added = g + g;
    assert_eq!(doubled, added);
}

#[test]
fn bn254_double_on_curve() {
    let g = Element::<Bn254G1Params>::one();
    let doubled = g.dbl();
    assert!(doubled.on_curve(), "2G should be on the curve");
}

#[test]
fn bn254_triple_on_curve() {
    let g = Element::<Bn254G1Params>::one();
    let tripled = g.dbl() + g;
    assert!(tripled.on_curve(), "3G should be on the curve");
}

#[test]
fn bn254_add_inverse_is_infinity() {
    let g = Element::<Bn254G1Params>::one();
    let neg_g = -g;
    let result = g + neg_g;
    assert!(result.is_point_at_infinity(), "G + (-G) should be infinity");
}

#[test]
fn bn254_mixed_add_matches_projective() {
    let g_proj = Element::<Bn254G1Params>::one();
    let g_affine = AffineElement::<Bn254G1Params>::one();
    let two_g = g_proj.dbl();

    // mixed add: 2G + G(affine)
    let three_g_mixed = two_g + g_affine;
    // projective add: 2G + G(projective)
    let three_g_proj = two_g + g_proj;

    assert_eq!(three_g_mixed, three_g_proj);
}

#[test]
fn bn254_to_affine_roundtrip() {
    let g = Element::<Bn254G1Params>::one();
    let two_g = g.dbl();
    let two_g_affine = two_g.to_affine();
    let two_g_back = Element::from_affine(&two_g_affine);
    assert_eq!(two_g, two_g_back);
}

#[test]
fn bn254_scalar_mul_small() {
    let g = Element::<Bn254G1Params>::one();
    // 3*G via scalar_mul
    let three_g = g.scalar_mul(&[3, 0, 0, 0]);
    // 3*G via add
    let three_g_expected = g.dbl() + g;
    assert_eq!(three_g, three_g_expected);
}

#[test]
fn bn254_scalar_mul_zero_is_infinity() {
    let g = Element::<Bn254G1Params>::one();
    let result = g.scalar_mul(&[0, 0, 0, 0]);
    assert!(result.is_point_at_infinity());
}

#[test]
fn bn254_add_infinity_identity() {
    let g = Element::<Bn254G1Params>::one();
    let inf = Element::<Bn254G1Params>::infinity();
    assert_eq!(g + inf, g);
    assert_eq!(inf + g, g);
}

#[test]
fn bn254_scalar_mul_by_order_is_infinity() {
    let g = Element::<Bn254G1Params>::one();
    // BN254 Fr modulus IS the group order
    let order = Bn254FrParams::MODULUS;
    let result = g.scalar_mul(&order);
    assert!(result.is_point_at_infinity(), "n*G should be infinity");
}

// =========================================================================
// Group tests — Grumpkin
// =========================================================================

#[test]
fn grumpkin_generator_on_curve() {
    let g = AffineElement::<GrumpkinG1Params>::one();
    assert!(g.on_curve(), "Grumpkin generator should be on the curve");
}

#[test]
fn grumpkin_double_on_curve() {
    let g = Element::<GrumpkinG1Params>::one();
    let doubled = g.dbl();
    assert!(doubled.on_curve(), "Grumpkin 2G should be on the curve");
}

#[test]
fn grumpkin_double_equals_add_self() {
    let g = Element::<GrumpkinG1Params>::one();
    assert_eq!(g.dbl(), g + g);
}

#[test]
fn grumpkin_scalar_mul_by_order_is_infinity() {
    let g = Element::<GrumpkinG1Params>::one();
    // Grumpkin's scalar field order = BN254 Fq modulus
    let order = Bn254FqParams::MODULUS;
    let result = g.scalar_mul(&order);
    assert!(result.is_point_at_infinity(), "Grumpkin: n*G should be infinity");
}

// =========================================================================
// Group tests — secp256k1
// =========================================================================

#[test]
fn secp256k1_generator_on_curve() {
    let g = AffineElement::<Secp256k1G1Params>::one();
    assert!(g.on_curve(), "secp256k1 generator should be on the curve");
}

#[test]
fn secp256k1_double_on_curve() {
    let g = Element::<Secp256k1G1Params>::one();
    let doubled = g.dbl();
    assert!(doubled.on_curve(), "secp256k1 2G should be on the curve");
}

#[test]
fn secp256k1_double_equals_add_self() {
    let g = Element::<Secp256k1G1Params>::one();
    assert_eq!(g.dbl(), g + g);
}

#[test]
fn secp256k1_add_inverse_is_infinity() {
    let g = Element::<Secp256k1G1Params>::one();
    let result = g + (-g);
    assert!(result.is_point_at_infinity());
}

// =========================================================================
// Group tests — secp256r1
// =========================================================================

#[test]
fn secp256r1_generator_on_curve() {
    let g = AffineElement::<Secp256r1G1Params>::one();
    assert!(g.on_curve(), "secp256r1 generator should be on the curve");
}

#[test]
fn secp256r1_generator_on_curve_projective() {
    let g = Element::<Secp256r1G1Params>::one();
    assert!(g.on_curve(), "secp256r1 generator should be on the curve (projective)");
}

#[test]
fn secp256r1_coeff_a_is_minus_3() {
    use crate::curves::secp256r1::Secp256r1Fq;
    use crate::groups::curve_params::CurveParams;
    let a = Secp256r1G1Params::coeff_a();
    let three = Secp256r1Fq::from(3u64);
    let neg_three = three.negate();
    assert_eq!(a, neg_three, "secp256r1 a should be -3");
}

#[test]
fn secp256r1_double_on_curve() {
    let g = Element::<Secp256r1G1Params>::one();
    let doubled = g.dbl();
    // Check via affine conversion first (to isolate projective on_curve issues)
    let doubled_affine = doubled.to_affine();
    assert!(doubled_affine.on_curve(), "secp256r1 2G should be on affine curve");
    assert!(doubled.on_curve(), "secp256r1 2G should be on the curve (projective)");
}

#[test]
fn secp256r1_double_equals_add_self() {
    let g = Element::<Secp256r1G1Params>::one();
    let dbl = g.dbl();
    let add = g + g;
    // Check raw field equality (bypassing on_curve check in PartialEq)
    let dbl_a = dbl.to_affine();
    let add_a = add.to_affine();
    assert_eq!(dbl_a.x, add_a.x, "dbl and add should give same x");
    assert_eq!(dbl_a.y, add_a.y, "dbl and add should give same y");
}

#[test]
fn secp256r1_projective_vs_affine_on_curve() {
    let g = Element::<Secp256r1G1Params>::one();
    let doubled = g.dbl();

    // Check projective on_curve
    let proj_on_curve = doubled.on_curve();
    eprintln!("projective on_curve: {}", proj_on_curve);

    // Manual affine check
    let affine = doubled.to_affine();
    let affine_on_curve = affine.on_curve();
    eprintln!("affine on_curve: {}", affine_on_curve);

    // Manual: check x_affine = X*Z^-2, y_affine = Y*Z^-3
    let z_inv = doubled.z.invert();
    let zz_inv = z_inv.sqr();
    let zzz_inv = zz_inv * z_inv;
    let x_aff = doubled.x * zz_inv;
    let y_aff = doubled.y * zzz_inv;
    eprintln!("affine.x: {:?}", affine.x);
    eprintln!("computed x: {:?}", x_aff);
    eprintln!("affine.y: {:?}", affine.y);
    eprintln!("computed y: {:?}", y_aff);
    assert_eq!(affine.x, x_aff);
    assert_eq!(affine.y, y_aff);

    // Now verify the algebra: if Y^2 = X^3 + aXZ^4 + bZ^6, then
    // (Y/Z^3)^2 = (X/Z^2)^3 + a(X/Z^2) + b
    // = X^3/Z^6 + aX/Z^2 + b
    // Multiply both sides by Z^6: Y^2 = X^3 + aXZ^4 + bZ^6
    // So if projective is correct, affine must be correct.
    assert!(proj_on_curve, "projective should be on curve");
    assert!(affine_on_curve, "if projective is on curve, affine must be too");
}

#[test]
fn secp256r1_generator_curve_eq_detailed() {
    use crate::curves::secp256r1::Secp256r1G1Params;
    use crate::groups::curve_params::CurveParams;

    let g = AffineElement::<Secp256r1G1Params>::one();
    let x = g.x;
    let y = g.y;
    let a = Secp256r1G1Params::coeff_a();
    let b = Secp256r1G1Params::coeff_b();

    // Verify x, y roundtrip correctly
    let x_back = x.from_montgomery_form();
    let y_back = y.from_montgomery_form();
    let expected_x = [0xF4A13945D898C296u64, 0x77037D812DEB33A0, 0xF8BCE6E563A440F2, 0x6B17D1F2E12C4247];
    let expected_y = [0xCBB6406837BF51F5u64, 0x2BCE33576B315ECE, 0x8EE7EB4A7C0F9E16, 0x4FE342E2FE1A7F9B];
    assert_eq!(x_back.data, expected_x, "generator x should roundtrip");
    assert_eq!(y_back.data, expected_y, "generator y should roundtrip");

    // Now check y^2 == x^3 + ax + b
    let y_sq = y.sqr();
    let x_sq = x.sqr();
    let x_cubed = x_sq * x;
    let ax = a * x;
    let rhs = x_cubed + ax + b;

    // Convert to standard form for inspection
    eprintln!("y^2 (standard): {:016x?}", y_sq.from_montgomery_form().data);
    eprintln!("x^3+ax+b (standard): {:016x?}", rhs.from_montgomery_form().data);
    eprintln!("x^3 (standard): {:016x?}", x_cubed.from_montgomery_form().data);
    eprintln!("ax (standard): {:016x?}", ax.from_montgomery_form().data);
    eprintln!("b (standard): {:016x?}", b.from_montgomery_form().data);

    assert_eq!(y_sq, rhs, "generator should satisfy curve equation");
}

#[test]
fn secp256r1_montgomery_roundtrip_large() {
    use crate::curves::secp256r1::Secp256r1FqParams;
    type R1Fq = Field<Secp256r1FqParams>;

    // Use the 2G x-coordinate as a test value
    let raw = [0xB1E6AA5FFC0A2A15u64, 0xC90F894FACD572DC, 0x8A52380304B51AC3, 0x7CF27B188D034F7E];
    let f = R1Fq::from_limbs(raw);
    let back = f.from_montgomery_form();
    eprintln!("raw:  {:016x?}", raw);
    eprintln!("mont: {:016x?}", f.data);
    eprintln!("back: {:016x?}", back.data);
    assert_eq!(back.data, raw, "Montgomery roundtrip should preserve large value");
}

#[test]
fn secp256r1_z_inverse_check() {
    let g = Element::<Secp256r1G1Params>::one();
    let doubled = g.dbl();
    let z = doubled.z;
    let z_inv = z.invert();
    let product = z * z_inv;
    assert_eq!(product, Field::<Secp256r1FqParams>::one(), "z * z^-1 should be 1");
}

#[test]
fn secp256r1_double_manual_curve_check() {
    use crate::curves::secp256r1::Secp256r1G1Params;
    use crate::groups::curve_params::CurveParams;

    let g = Element::<Secp256r1G1Params>::one();
    let doubled = g.dbl();
    let p = doubled.to_affine();
    let x = p.x;
    let y = p.y;
    let a = Secp256r1G1Params::coeff_a();
    let b = Secp256r1G1Params::coeff_b();
    let lhs = y.sqr();
    let rhs_no_a = x.sqr() * x + b;
    let rhs_with_a = rhs_no_a + a * x;
    eprintln!("lhs (y^2):      {:?}", lhs);
    eprintln!("rhs (x^3+ax+b): {:?}", rhs_with_a);
    eprintln!("rhs (x^3+b):    {:?}", rhs_no_a);
    eprintln!("a*x:            {:?}", a * x);
    eprintln!("a:              {:?}", a);
    eprintln!("b:              {:?}", b);
    assert_eq!(lhs, rhs_with_a, "doubled point should satisfy curve equation");
}

#[test]
fn secp256r1_coeff_a_is_neg_three() {
    use crate::curves::secp256r1::{Secp256r1FqParams, Secp256r1G1Params};
    use crate::groups::curve_params::CurveParams;
    type R1Fq = Field<Secp256r1FqParams>;

    let a = Secp256r1G1Params::coeff_a();
    let three = R1Fq::from(3u64);
    let sum = a + three;
    eprintln!("a:     {:?}", a);
    eprintln!("three: {:?}", three);
    eprintln!("a+3:   {:?}", sum);
    assert!(sum.is_zero(), "coeff_a + 3 should be zero (a = -3)");

    // Also check from_montgomery_form roundtrip
    let a_standard = a.from_montgomery_form();
    let expected = [0xFFFFFFFFFFFFFFFCu64, 0x00000000FFFFFFFF, 0, 0xFFFFFFFF00000001];
    eprintln!("a (standard form): {:016x} {:016x} {:016x} {:016x}",
        a_standard.data[3], a_standard.data[2], a_standard.data[1], a_standard.data[0]);
    assert_eq!(a_standard.data, expected, "a in standard form should be p-3");
}

#[test]
fn secp256r1_field_distributive() {
    // Test distributive law: a * (b + c) == a*b + a*c
    use crate::curves::secp256r1::Secp256r1FqParams;
    use crate::groups::curve_params::CurveParams;
    type R1Fq = Field<Secp256r1FqParams>;
    let a = R1Fq::from(12345u64);
    let b = R1Fq::from(67890u64);
    let c = R1Fq::from(11111u64);
    let lhs = a * (b + c);
    let rhs = a * b + a * c;
    assert_eq!(lhs, rhs, "distributive law should hold for secp256r1 Fq");

    // Also test with the generator coordinates
    let g = AffineElement::<Secp256r1G1Params>::one();
    let x = g.x;
    let y = g.y;
    // Verify y^2 = x^3 + ax + b
    let a_coeff = Secp256r1G1Params::coeff_a();
    let b_coeff = Secp256r1G1Params::coeff_b();
    let y_sq = y.sqr();
    let x_cubed = x.sqr() * x;
    let rhs = x_cubed + a_coeff * x + b_coeff;
    assert_eq!(y_sq, rhs, "generator should satisfy curve equation manually");
}

// known 2G test removed - the coordinates used were incorrect

#[test]
fn secp256r1_affine_double() {
    use crate::curves::secp256r1::{Secp256r1FqParams, Secp256r1G1Params};
    use crate::groups::curve_params::CurveParams;
    type R1Fq = Field<Secp256r1FqParams>;

    let g = AffineElement::<Secp256r1G1Params>::one();
    let x0 = g.x;
    let y0 = g.y;
    let a = Secp256r1G1Params::coeff_a();

    // Manual affine doubling: lambda = (3*x0^2 + a) / (2*y0)
    let three_x0_sq = x0.sqr() * R1Fq::from(3u64);
    let numerator = three_x0_sq + a;
    let denominator = y0 + y0;
    let lambda = numerator * denominator.invert();
    let x1 = lambda.sqr() - (x0 + x0);
    let y1 = lambda * (x0 - x1) - y0;

    // Check each intermediate step
    assert_eq!(denominator * denominator.invert(), R1Fq::one(), "inverse of 2y check");

    // Verify result via explicit curve equation
    let y1_sq = y1.sqr();
    let x1_cubed = x1.sqr() * x1;
    let rhs = x1_cubed + a * x1 + Secp256r1G1Params::coeff_b();
    eprintln!("x1:         {:?}", x1);
    eprintln!("y1:         {:?}", y1);
    eprintln!("y1^2:       {:?}", y1_sq);
    eprintln!("x1^3+ax+b:  {:?}", rhs);
    eprintln!("x1^3:       {:?}", x1_cubed);
    eprintln!("a*x1:       {:?}", a * x1);
    assert_eq!(y1_sq, rhs, "manually doubled secp256r1 point should satisfy curve eq");
}

#[test]
fn secp256r1_fq_from_consistency() {
    type R1Fq = Field<Secp256r1FqParams>;
    let three = R1Fq::from(3u64);
    let three_add = R1Fq::one() + R1Fq::one() + R1Fq::one();
    assert_eq!(three, three_add, "from(3) should equal 1+1+1");
    assert_eq!(R1Fq::from(3u64) * R1Fq::from(5u64), R1Fq::from(15u64), "3*5=15");
}

#[test]
fn secp256r1_fq_associativity() {
    type R1Fq = Field<Secp256r1FqParams>;
    let x = R1Fq::from_limbs([0xF4A13945D898C296, 0x77037D812DEB33A0, 0xF8BCE6E563A440F2, 0x6B17D1F2E12C4247]);
    let y = R1Fq::from_limbs([0xCBB6406837BF51F5, 0x2BCE33576B315ECE, 0x8EE7EB4A7C0F9E16, 0x4FE342E2FE1A7F9B]);
    let a = R1Fq::from_limbs([0xFFFFFFFFFFFFFFFC, 0x00000000FFFFFFFF, 0, 0xFFFFFFFF00000001]);
    // (x * y) * a == x * (y * a)
    assert_eq!((x * y) * a, x * (y * a), "multiplication should be associative");
    // (x + y) + a == x + (y + a)
    assert_eq!((x + y) + a, x + (y + a), "addition should be associative");
}

#[test]
fn secp256r1_lambda_computation() {
    use crate::curves::secp256r1::{Secp256r1FqParams, Secp256r1G1Params};
    use crate::groups::curve_params::CurveParams;
    type R1Fq = Field<Secp256r1FqParams>;

    let g = AffineElement::<Secp256r1G1Params>::one();
    let x0 = g.x;
    let y0 = g.y;
    let a = Secp256r1G1Params::coeff_a();

    let x0_sq = x0.sqr();
    let three_x0_sq = x0_sq + x0_sq + x0_sq; // via addition (no mul by 3)
    let numerator = three_x0_sq + a;
    let denominator = y0 + y0;
    let lambda = numerator * denominator.invert();

    // lambda should satisfy: lambda^2 - 2*x0 gives a valid x1
    // and y0 = lambda*(x0 - x1) - y1 for some y1 on the curve
    let x1 = lambda.sqr() - (x0 + x0);
    let y1 = lambda * (x0 - x1) - y0;

    // Check the expected 2G coordinates
    let expected_x = R1Fq::from_limbs([0xCC0A48CE78BE3DE4, 0xC90183A1B62FBAFA, 0x8A52380304B51AC3, 0x7CF27B188D034F7E]);
    let expected_y = R1Fq::from_limbs([0x9E04B79D227873D1, 0xBA7DADE63CE98229, 0x293D9AC69F7430DB, 0x07775510DB8ED040]);

    eprintln!("x1 computed:  {:?}", x1);
    eprintln!("x1 expected:  {:?}", expected_x);
    eprintln!("x1 match: {}", x1 == expected_x);
    eprintln!("y1 computed:  {:?}", y1);
    eprintln!("y1 expected:  {:?}", expected_y);
    eprintln!("y1 match: {}", y1 == expected_y);

    // Verify lambda: lambda * denominator == numerator
    assert_eq!(lambda * denominator, numerator, "lambda should satisfy lambda * 2y = 3x^2 + a");

    // Verify x1: lambda^2 == x1 + 2*x0
    assert_eq!(lambda.sqr(), x1 + x0 + x0, "lambda^2 should equal x1 + 2*x0");

    // The key question: do our x1,y1 satisfy the curve equation?
    let y1_sq = y1.sqr();
    let b = Secp256r1G1Params::coeff_b();
    let rhs = x1.sqr() * x1 + a * x1 + b;
    eprintln!("y1^2 = {:?}", y1_sq);
    eprintln!("x1^3+ax+b = {:?}", rhs);
    assert_eq!(y1_sq, rhs, "affine-doubled point should be on curve");
}

#[test]
fn secp256r1_fq_sqr_equals_mul() {
    type R1Fq = Field<Secp256r1FqParams>;
    let x = R1Fq::from_limbs([0xF4A13945D898C296, 0x77037D812DEB33A0, 0xF8BCE6E563A440F2, 0x6B17D1F2E12C4247]);
    assert_eq!(x.sqr(), x * x, "sqr should equal mul for secp256r1 Fq");
}

#[test]
fn secp256r1_fq_mul3_vs_add3() {
    type R1Fq = Field<Secp256r1FqParams>;
    let x = R1Fq::from_limbs([0xF4A13945D898C296, 0x77037D812DEB33A0, 0xF8BCE6E563A440F2, 0x6B17D1F2E12C4247]);
    let x_sq = x.sqr();
    let via_mul = x_sq * R1Fq::from(3u64);
    let via_add = x_sq + x_sq + x_sq;
    eprintln!("x_sq * 3:     {:?}", via_mul);
    eprintln!("x_sq+x_sq+sq: {:?}", via_add);
    eprintln!("mul raw: [{:#018x}, {:#018x}, {:#018x}, {:#018x}]", via_mul.data[0], via_mul.data[1], via_mul.data[2], via_mul.data[3]);
    eprintln!("add raw: [{:#018x}, {:#018x}, {:#018x}, {:#018x}]", via_add.data[0], via_add.data[1], via_add.data[2], via_add.data[3]);
    assert_eq!(via_mul, via_add, "multiplying by 3 should equal adding three times");
    // Also check raw equality
    if via_mul.data != via_add.data {
        eprintln!("WARNING: same field element, different raw representation!");
    }
}

#[test]
fn secp256r1_fq_montgomery_roundtrip_large() {
    type R1Fq = Field<Secp256r1FqParams>;
    let raw = [0xF4A13945D898C296u64, 0x77037D812DEB33A0, 0xF8BCE6E563A440F2, 0x6B17D1F2E12C4247];
    let x = R1Fq::from_limbs(raw);
    let back = x.from_montgomery_form();
    assert_eq!(back.data, raw, "Montgomery roundtrip should preserve value");
}

#[test]
fn secp256r1_fq_mul_then_inverse() {
    type R1Fq = Field<Secp256r1FqParams>;
    let x = R1Fq::from_limbs([0xF4A13945D898C296, 0x77037D812DEB33A0, 0xF8BCE6E563A440F2, 0x6B17D1F2E12C4247]);
    let y = R1Fq::from_limbs([0xCBB6406837BF51F5, 0x2BCE33576B315ECE, 0x8EE7EB4A7C0F9E16, 0x4FE342E2FE1A7F9B]);
    let xy = x * y;
    let x_back = xy * y.invert();
    assert_eq!(x, x_back, "x * y * y^-1 should equal x");
}

#[test]
fn secp256r1_scalar_mul_small() {
    let g = Element::<Secp256r1G1Params>::one();
    let five_g = g.scalar_mul(&[5, 0, 0, 0]);
    let five_g_expected = g.dbl().dbl() + g;
    assert_eq!(five_g, five_g_expected);
}

// =========================================================================
// Field byte serialization tests
// =========================================================================

#[test]
fn field_to_be_bytes_roundtrip() {
    let val = Fq::from(42u64);
    let bytes = val.to_be_bytes();
    let recovered = Fq::from_be_bytes(&bytes);
    assert_eq!(val, recovered, "to_be_bytes / from_be_bytes roundtrip failed");
}

#[test]
fn field_to_be_bytes_known_value() {
    let one = Fq::from(1u64);
    let bytes = one.to_be_bytes();
    // 1 in big-endian 32 bytes: all zeros except last byte = 1
    let mut expected = [0u8; 32];
    expected[31] = 1;
    assert_eq!(bytes, expected, "to_be_bytes(1) should be 0x...01");
}

#[test]
fn field_from_be_bytes_larger_than_modulus() {
    // All 0xFF bytes = 2^256 - 1, should reduce mod p
    let bytes = [0xFFu8; 32];
    let val = Fq::from_be_bytes(&bytes);
    // Should be equivalent to (2^256 - 1) mod p
    let max_limbs = [0xFFFFFFFFFFFFFFFFu64; 4];
    let expected = Fq::from_limbs(max_limbs);
    assert_eq!(val, expected, "from_be_bytes should reduce values >= modulus");
}

#[test]
fn field_to_be_bytes_roundtrip_fr() {
    let val = Fr::from_limbs([0x123456789ABCDEF0, 0xFEDCBA9876543210, 0x1111111111111111, 0x0ABCDEF012345678]);
    let bytes = val.to_be_bytes();
    let recovered = Fr::from_be_bytes(&bytes);
    assert_eq!(val, recovered, "Fr roundtrip failed");
}

// =========================================================================
// Field::get_bit tests
// =========================================================================

#[test]
fn field_get_bit() {
    let val = Fq::from(5u64); // binary: 101
    assert!(val.get_bit(0), "bit 0 should be 1");
    assert!(!val.get_bit(1), "bit 1 should be 0");
    assert!(val.get_bit(2), "bit 2 should be 1");
    assert!(!val.get_bit(3), "bit 3 should be 0");
}

// =========================================================================
// Field::random_element tests
// =========================================================================

#[test]
fn field_random_element_not_zero() {
    let r = Fr::random_element();
    assert!(!r.is_zero(), "random element should not be zero (probabilistically)");
}

#[test]
fn field_random_element_different_values() {
    let r1 = Fr::random_element();
    let r2 = Fr::random_element();
    assert_ne!(r1, r2, "two random elements should differ (probabilistically)");
}

// =========================================================================
// Element::double_scalar_mul tests
// =========================================================================

#[test]
fn element_double_scalar_mul_bn254() {
    let g = Element::<Bn254G1Params>::one();
    let a = Fr::from(3u64);
    let b = Fr::from(5u64);
    let two_g = g.dbl();
    // Compute 3*G + 5*2G = 3G + 10G = 13G
    let result = g.double_scalar_mul(&a, &two_g, &b);
    let expected = g.scalar_mul(&[13, 0, 0, 0]);
    assert_eq!(result, expected, "3*G + 5*2G should equal 13*G");
}

#[test]
fn element_double_scalar_mul_identity() {
    let g = Element::<Bn254G1Params>::one();
    let zero = Fr::zero();
    let one_scalar = Fr::one();
    // 1*G + 0*G = G
    let result = g.double_scalar_mul(&one_scalar, &g, &zero);
    assert_eq!(result, g, "1*G + 0*G should equal G");
}

// =========================================================================
// mul_512 tests
// =========================================================================

#[test]
fn field_mul_512_known_values() {
    // Multiply two known field elements (in raw form, not Montgomery)
    // and verify the 512-bit result
    let a = Fr::from_raw([2, 0, 0, 0]);
    let b = Fr::from_raw([3, 0, 0, 0]);
    let result = a.mul_512(&b);
    assert_eq!(result[0], 6);
    for i in 1..8 {
        assert_eq!(result[i], 0, "limb {} should be 0", i);
    }
}

#[test]
fn field_mul_512_large_values() {
    // Multiply max u64 by itself: (2^64-1)^2 = 2^128 - 2^65 + 1
    let a = Fr::from_raw([u64::MAX, 0, 0, 0]);
    let b = Fr::from_raw([u64::MAX, 0, 0, 0]);
    let result = a.mul_512(&b);
    // (2^64-1)^2 = 0xFFFFFFFFFFFFFFFE_0000000000000001
    assert_eq!(result[0], 1);
    assert_eq!(result[1], 0xFFFFFFFFFFFFFFFE);
    for i in 2..8 {
        assert_eq!(result[i], 0, "limb {} should be 0", i);
    }
}

#[test]
fn field_mul_512_cross_limb() {
    // Verify cross-limb multiplication: a = [0, 1, 0, 0] * b = [0, 1, 0, 0]
    // = 2^64 * 2^64 = 2^128 -> result[2] = 1
    let a = Fr::from_raw([0, 1, 0, 0]);
    let b = Fr::from_raw([0, 1, 0, 0]);
    let result = a.mul_512(&b);
    assert_eq!(result[0], 0);
    assert_eq!(result[1], 0);
    assert_eq!(result[2], 1);
    for i in 3..8 {
        assert_eq!(result[i], 0, "limb {} should be 0", i);
    }
}

// =========================================================================
// split_into_endomorphism_scalars tests
// Matches C++ fr.test.cpp SplitIntoEndomorphismScalars / SplitIntoEndomorphismScalarsSimple
// =========================================================================

/// Helper to verify BN254 Fr endomorphism split: k1 - k2*lambda == k
/// Matches C++ fr.test.cpp SplitIntoEndomorphismScalars exactly:
///   fr::split_into_endomorphism_scalars(k, k1, k2);
///   k1.self_to_montgomery_form();
///   k2.self_to_montgomery_form();
///   fr result = k1 - k2 * lambda;
///   result.self_from_montgomery_form();
///   EXPECT_EQ(result, k);
fn verify_bn254_endo_split(k: Fr) {
    let (k1, k2) = k.split_into_endomorphism_scalars();

    // C++: k1.self_to_montgomery_form(); k2.self_to_montgomery_form();
    let k1_mont = k1.to_montgomery_form();
    let k2_mont = k2.to_montgomery_form();

    let lambda = Fr::from_raw(Bn254FrParams::CUBE_ROOT);
    let result = (k1_mont - k2_mont * lambda).from_montgomery_form();

    assert_eq!(result, k, "k1 - k2*lambda should equal k");
}

#[test]
fn bn254_fr_split_endomorphism_scalars() {
    // C++ SplitIntoEndomorphismScalars: fr k = fr::random_element()
    let k = Fr::random_element();
    verify_bn254_endo_split(k);
}

#[test]
fn bn254_fr_split_endomorphism_scalars_simple() {
    // C++ SplitIntoEndomorphismScalarsSimple: fr input = {1, 0, 0, 0}
    let k = Fr::from_raw([1, 0, 0, 0]);
    verify_bn254_endo_split(k);
}

#[test]
fn secp256k1_fr_split_endomorphism_scalars() {
    // C++ secp256k1::GetEndomorphismScalars: 2048 iterations with random_element()
    // C++: split_into_endomorphism_scalars(k, k1, k2);
    //      if (k2.get_msb() > 200) { k2 = -k2; k2_neg = true; }
    //      EXPECT_LT(k1.get_msb(), 129); EXPECT_LT(k2.get_msb(), 129);
    //      if (k2_neg) k2 = -k2;  // undo
    //      k1.self_to_montgomery_form(); k2.self_to_montgomery_form();
    //      expected = k1 - k2 * beta; expected.self_from_montgomery_form();
    //      EXPECT_EQ(k, expected);
    for _ in 0..2048 {
        let k = Secp256k1Fr::random_element();
        let (k1, k2) = k.split_into_endomorphism_scalars();

        // C++: k1.self_to_montgomery_form(); k2.self_to_montgomery_form();
        let k1_mont = k1.to_montgomery_form();
        let k2_mont = k2.to_montgomery_form();
        let beta = Secp256k1Fr::from_raw(Secp256k1FrParams::CUBE_ROOT);
        let result = (k1_mont - k2_mont * beta).from_montgomery_form();
        assert_eq!(result, k, "secp256k1: k1 - k2*beta should equal k");
    }
}

// =========================================================================
// WNAF tests
// =========================================================================

#[test]
fn wnaf_encode_basic_properties() {
    use crate::groups::wnaf;

    // Test with odd scalar
    let scalar: [u64; 2] = [0xDEADBEEF12345679, 0x1234];
    let wnaf_bits = 4usize;
    let wnaf_entries = (wnaf::SCALAR_BITS + wnaf_bits - 1) / wnaf_bits; // 32
    let mut table = vec![0u64; wnaf_entries + 1];
    let mut skew = false;

    wnaf::fixed_wnaf(&scalar, &mut table, &mut skew, 0, 1, wnaf_bits);

    // Odd scalar -> no skew
    assert!(!skew, "odd scalar should not be skewed");

    // All entries should have abs value < 2^(wnaf_bits-1) = 8
    for i in 0..wnaf_entries {
        let abs_val = table[i] & 0x0fffffff;
        assert!(abs_val < 8, "WNAF entry {} abs value {} >= 8", i, abs_val);
    }
}

#[test]
fn wnaf_even_scalar_skew() {
    use crate::groups::wnaf;

    let scalar: [u64; 2] = [0xDEADBEEF12345678, 0x1234]; // even scalar
    let wnaf_bits = 4usize;
    let wnaf_entries = (wnaf::SCALAR_BITS + wnaf_bits - 1) / wnaf_bits;
    let mut table = vec![0u64; wnaf_entries + 1];
    let mut skew = false;

    wnaf::fixed_wnaf(&scalar, &mut table, &mut skew, 0, 1, wnaf_bits);

    assert!(skew, "even scalar should be skewed");
}

#[test]
fn wnaf_interleaved_two_scalars() {
    use crate::groups::wnaf;

    // Test interleaved encoding with num_points=2 (as used by endomorphism mul)
    let k1: [u64; 2] = [0x123456789ABCDEF1, 0x42];
    let k2: [u64; 2] = [0xFEDCBA9876543211, 0x13];

    let wnaf_entries = (wnaf::SCALAR_BITS + 4 - 1) / 4; // 32
    let mut table = vec![0u64; wnaf_entries * 2];
    let mut skew1 = false;
    let mut skew2 = false;

    wnaf::fixed_wnaf(&k1, &mut table, &mut skew1, 0, 2, 4);
    wnaf::fixed_wnaf(&k2, &mut table[1..], &mut skew2, 0, 2, 4);

    // k1 is odd, k2 is odd
    assert!(!skew1, "k1 should not be skewed");
    assert!(!skew2, "k2 should not be skewed");

    // Even-indexed entries belong to k1, odd-indexed to k2
    for i in (0..wnaf_entries * 2).step_by(2) {
        let abs_val = table[i] & 0x0fffffff;
        assert!(abs_val < 8, "k1 WNAF entry abs value {} >= 8", abs_val);
    }
    for i in (1..wnaf_entries * 2).step_by(2) {
        let abs_val = table[i] & 0x0fffffff;
        assert!(abs_val < 8, "k2 WNAF entry abs value {} >= 8", abs_val);
    }
}

#[test]
fn wnaf_fixed_with_endo_split() {
    // Matches C++ WnafFixedWithEndoSplit:
    //   fr k = engine.get_random_uint256(); k.data[3] &= 0x0fffffffffffffffUL;
    //   split(k) -> k1, k2
    //   fixed_wnaf each with wnaf_bits=5
    //   recover from WNAF
    //   verify k1_recovered - k2_recovered * lambda == k
    use crate::groups::wnaf;

    // C++: k = random, top limb masked to fit < 2^252
    let mut k = Fr::random_element();
    k.data[3] &= 0x0fffffffffffffff;

    let (k1_full, k2_full) = k.split_into_endomorphism_scalars();
    let k1_limbs = [k1_full.data[0], k1_full.data[1]];
    let k2_limbs = [k2_full.data[0], k2_full.data[1]];

    // C++: WNAF_SIZE(5) = (127 + 5 - 1) / 5 = 26
    let wnaf_bits = 5usize;
    let wnaf_entries = (wnaf::SCALAR_BITS + wnaf_bits - 1) / wnaf_bits;
    let mut wnaf_table = vec![0u64; wnaf_entries + 1];
    let mut endo_wnaf_table = vec![0u64; wnaf_entries + 1];
    let mut skew = false;
    let mut endo_skew = false;

    // C++: wnaf::fixed_wnaf<1, 5>(&k1.data[0], wnaf, skew, 0);
    wnaf::fixed_wnaf(&k1_limbs, &mut wnaf_table, &mut skew, 0, 1, wnaf_bits);
    wnaf::fixed_wnaf(&k2_limbs, &mut endo_wnaf_table, &mut endo_skew, 0, 1, wnaf_bits);

    // C++: recover_fixed_wnaf(wnaf, skew, k1_recovered.data[1], k1_recovered.data[0], 5);
    fn recover_wnaf(table: &[u64], skew: bool, num_entries: usize, wnaf_bits: usize) -> [u64; 2] {
        let mut recovered: i128 = 0;
        for i in 0..num_entries {
            let entry = table[i];
            let val = (entry & 0x0fffffff) as i128;
            let sign = ((entry >> 31) & 1) != 0;
            let actual = (2 * val + 1) * if sign { -1 } else { 1 };
            recovered = recovered * (1i128 << wnaf_bits) + actual;
        }
        if skew {
            recovered -= 1;
        }
        [recovered as u64, (recovered >> 64) as u64]
    }

    let k1_recovered = recover_wnaf(&wnaf_table, skew, wnaf_entries, wnaf_bits);
    let k2_recovered = recover_wnaf(&endo_wnaf_table, endo_skew, wnaf_entries, wnaf_bits);

    // C++: result = k2_recovered * lambda; result = k1_recovered - result; EXPECT_EQ(result, k);
    let k1 = Fr::from_limbs([k1_recovered[0], k1_recovered[1], 0, 0]);
    let k2 = Fr::from_limbs([k2_recovered[0], k2_recovered[1], 0, 0]);
    let lambda = Fr::from_raw(Bn254FrParams::CUBE_ROOT);
    let result = (k1 - k2 * lambda).from_montgomery_form();
    assert_eq!(result, k, "WNAF+endo split recovery should match original k");
}

// =========================================================================
// mul_with_endomorphism tests
// Matches C++ affine_element.test.cpp MulWithEndomorphismMatchesMulWithoutEndomorphism:
//   for 100 random elements, verify mul_with_endomorphism == mul_without_endomorphism
// =========================================================================

#[test]
fn bn254_mul_with_endomorphism_matches_basic() {
    // C++ MulWithEndomorphismMatchesMulWithoutEndomorphism:
    //   for (int i = 0; i < 100; i++) {
    //       auto x1 = element(affine_element::random_element());
    //       auto f1 = fr::random_element();
    //       auto r1 = mul_without_endomorphism(x1, f1);
    //       auto r2 = mul_with_endomorphism(x1, f1);
    //       EXPECT_EQ(r1, r2);
    //   }
    let g = Element::<Bn254G1Params>::one();
    for _ in 0..100 {
        let x1 = g.mul_without_endomorphism(&Fr::random_element());
        let f1 = Fr::random_element();
        let r1 = x1.mul_without_endomorphism(&f1);
        let r2 = x1.mul_with_endomorphism(&f1);
        assert_eq!(r1, r2, "endomorphism mul should match basic mul");
    }
}

#[test]
fn bn254_mul_with_endomorphism_by_minus_one() {
    // C++ TestAffineElement BatchEndomorphismByMinusOne: -1 * P == -P
    let g = Element::<Bn254G1Params>::one();
    let minus_one = Fr::one().negate();
    let result = g.mul_with_endomorphism(&minus_one);
    let expected = -g;
    assert_eq!(result, expected, "-1 * G should equal -G");
}

#[test]
fn bn254_mul_dispatch() {
    let g = Element::<Bn254G1Params>::one();
    let scalar = Fr::from(42u64);
    let expected = g.mul_without_endomorphism(&scalar);
    let result = g.mul(&scalar);
    assert_eq!(result, expected, "mul() should match mul_without_endomorphism");
}

// Note: C++ mul_with_endomorphism has a static_assert(modulus_3 < 0x4000...)
// in the pair-returning split, so it only compiles for 254-bit scalar fields.
// secp256k1 (256-bit modulus) does NOT use mul_with_endomorphism in C++.
// The secp256k1 endomorphism split is only used in tests and Pippenger MSM.

#[test]
fn mul_with_endomorphism_zero_returns_infinity() {
    let g = Element::<Bn254G1Params>::one();
    let zero = Fr::zero();
    let result = g.mul_with_endomorphism(&zero);
    assert!(result.is_point_at_infinity(), "mul by zero should be infinity");
}

#[test]
fn mul_with_endomorphism_infinity_returns_infinity() {
    let inf = Element::<Bn254G1Params>::infinity();
    let scalar = Fr::from(42u64);
    let result = inf.mul_with_endomorphism(&scalar);
    assert!(
        result.is_point_at_infinity(),
        "mul on infinity should return infinity"
    );
}

// =========================================================================
// batch_normalize tests
// =========================================================================

// Matches C++ g1.test.cpp / grumpkin.test.cpp BatchNormalize:
//   Creates random points with non-trivial z (a + b where a,b are random),
//   batch normalizes, then verifies x_normalized * z_orig^2 == x_orig
//   and y_normalized * z_orig^3 == y_orig.
#[test]
fn batch_normalize_matches_individual() {
    let g = Element::<Bn254G1Params>::one();
    let num_points = 2; // C++ uses num_points = 2

    // C++: element a = element::random_element(); element b = element::random_element(); points[i] = a + b;
    let mut points: Vec<Element<Bn254G1Params>> = Vec::new();
    let mut normalized: Vec<Element<Bn254G1Params>> = Vec::new();
    for _ in 0..num_points {
        let a = g.mul_without_endomorphism(&Fr::random_element());
        let b = g.mul_without_endomorphism(&Fr::random_element());
        let point = a + b;
        points.push(point);
        normalized.push(point);
    }

    Element::batch_normalize(&mut normalized);

    // C++ verification: result_x = normalized[i].x * zz; EXPECT_EQ(result_x == points[i].x, true);
    for i in 0..num_points {
        let zz = points[i].z.sqr();
        let zzz = points[i].z * zz;
        let result_x = normalized[i].x * zz;
        let result_y = normalized[i].y * zzz;
        assert_eq!(result_x, points[i].x, "batch_normalize x mismatch at {}", i);
        assert_eq!(result_y, points[i].y, "batch_normalize y mismatch at {}", i);
    }
}

#[test]
fn batch_normalize_with_infinity() {
    let g = Element::<Bn254G1Params>::one();
    let point0 = g.mul_without_endomorphism(&Fr::random_element())
        + g.mul_without_endomorphism(&Fr::random_element());
    let point2 = g.mul_without_endomorphism(&Fr::random_element())
        + g.mul_without_endomorphism(&Fr::random_element());

    let orig0 = point0;
    let orig2 = point2;

    let mut points = vec![point0, Element::infinity(), point2];
    Element::batch_normalize(&mut points);

    // Verify non-infinity points
    let zz = orig0.z.sqr();
    let zzz = orig0.z * zz;
    assert_eq!(points[0].x * zz, orig0.x);
    assert_eq!(points[0].y * zzz, orig0.y);

    let zz = orig2.z.sqr();
    let zzz = orig2.z * zz;
    assert_eq!(points[2].x * zz, orig2.x);
    assert_eq!(points[2].y * zzz, orig2.y);
}

#[test]
fn batch_normalize_empty() {
    let mut points: Vec<Element<Bn254G1Params>> = Vec::new();
    Element::batch_normalize(&mut points); // Should not panic
}

#[test]
fn batch_normalize_single() {
    let g = Element::<Bn254G1Params>::one();
    let orig = g.mul_without_endomorphism(&Fr::random_element())
        + g.mul_without_endomorphism(&Fr::random_element());
    let mut points = vec![orig];
    Element::batch_normalize(&mut points);

    let zz = orig.z.sqr();
    let zzz = orig.z * zz;
    assert_eq!(points[0].x * zz, orig.x);
    assert_eq!(points[0].y * zzz, orig.y);
}

// =========================================================================
// Missing field tests — ported from C++ fr.test.cpp / fq.test.cpp
// =========================================================================

#[test]
fn bn254_fr_montgomery_consistency_check() {
    for _ in 0..100 {
        let a = Fr::random_element();
        let roundtrip = a.from_montgomery_form().to_montgomery_form();
        assert_eq!(a, roundtrip, "from_mont(to_mont(a)) should equal a");
    }
}

#[test]
fn bn254_fq_montgomery_consistency_check() {
    for _ in 0..100 {
        let a = Fq::random_element();
        let roundtrip = a.from_montgomery_form().to_montgomery_form();
        assert_eq!(a, roundtrip, "from_mont(to_mont(a)) should equal a");
    }
}

#[test]
fn bn254_fr_add_mul_consistency() {
    for _ in 0..100 {
        let a = Fr::random_element();
        let b = Fr::random_element();
        let c = Fr::random_element();
        let lhs = (a + b) * c;
        let rhs = a * c + b * c;
        assert_eq!(lhs, rhs, "(a+b)*c should equal a*c + b*c");
    }
}

#[test]
fn bn254_fr_sub_mul_consistency() {
    for _ in 0..100 {
        let a = Fr::random_element();
        let b = Fr::random_element();
        let c = Fr::random_element();
        let lhs = (a - b) * c;
        let rhs = a * c - b * c;
        assert_eq!(lhs, rhs, "(a-b)*c should equal a*c - b*c");
    }
}

#[test]
fn bn254_fr_mul_sqr_consistency() {
    for _ in 0..100 {
        let a = Fr::random_element();
        assert_eq!(a * a, a.sqr(), "a*a should equal a.sqr()");
    }
}

#[test]
fn bn254_fq_add_mul_consistency() {
    for _ in 0..100 {
        let a = Fq::random_element();
        let b = Fq::random_element();
        let c = Fq::random_element();
        let lhs = (a + b) * c;
        let rhs = a * c + b * c;
        assert_eq!(lhs, rhs, "(a+b)*c should equal a*c + b*c");
    }
}

#[test]
fn bn254_fq_sub_mul_consistency() {
    for _ in 0..100 {
        let a = Fq::random_element();
        let b = Fq::random_element();
        let c = Fq::random_element();
        let lhs = (a - b) * c;
        let rhs = a * c - b * c;
        assert_eq!(lhs, rhs, "(a-b)*c should equal a*c - b*c");
    }
}

#[test]
fn bn254_fq_mul_sqr_consistency() {
    for _ in 0..100 {
        let a = Fq::random_element();
        assert_eq!(a * a, a.sqr(), "a*a should equal a.sqr()");
    }
}

#[test]
fn bn254_fr_lambda() {
    let cube_root = Fr::cube_root_of_unity();
    let one = Fr::one();
    assert_ne!(cube_root, one, "cube root should not be 1");
    assert_eq!(cube_root * cube_root * cube_root, one, "cube_root^3 should be 1");
}

#[test]
fn bn254_fq_beta() {
    let cube_root = Fq::cube_root_of_unity();
    let one = Fq::one();
    assert_ne!(cube_root, one, "cube root should not be 1");
    assert_eq!(cube_root * cube_root * cube_root, one, "cube_root^3 should be 1");
}

#[test]
fn bn254_fr_invert_one_is_one() {
    assert_eq!(Fr::one().invert(), Fr::one(), "1.invert() should be 1");
}

#[test]
fn bn254_fq_invert_one_is_one() {
    assert_eq!(Fq::one().invert(), Fq::one(), "1.invert() should be 1");
}

#[test]
fn bn254_fr_sqrt_random() {
    for _ in 0..100 {
        let a = Fr::random_element();
        let a_sq = a.sqr();
        let (found, root) = a_sq.sqrt();
        assert!(found, "a^2 should have a sqrt");
        assert!(root == a || root == a.negate(), "sqrt(a^2) should be +/- a");
    }
}

#[test]
fn bn254_fq_sqrt_random() {
    for _ in 0..100 {
        let a = Fq::random_element();
        let a_sq = a.sqr();
        let (found, root) = a_sq.sqrt();
        assert!(found, "a^2 should have a sqrt");
        assert!(root == a || root == a.negate(), "sqrt(a^2) should be +/- a");
    }
}

// --- secp256k1 regression tests ---

#[test]
fn secp256k1_fq_neg_and_self_neg_zero() {
    type K1Fq = Field<Secp256k1FqParams>;
    let a = K1Fq::zero();
    let a_neg = a.negate();
    assert_eq!(a, a_neg, "-0 should equal 0");
}

#[test]
fn secp256k1_fq_montgomery_mul_big_bug() {
    type K1Fq = Field<Secp256k1FqParams>;
    let a = K1Fq::from_limbs([0xfffffffe630dc02f, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff]);
    let a_sqr = a.sqr();
    let expected = K1Fq::from_limbs([0x60381e557e100000, 0, 0, 0]);
    assert_eq!(a_sqr, expected, "secp256k1 Fq MontgomeryMulBigBug regression");
}

// --- secp256r1 regression tests ---

#[test]
fn secp256r1_fr_montgomery_mul_big_bug() {
    type R1Fr = Field<Secp256r1FrParams>;
    // C++ sets a.data[] directly — this is raw Montgomery-form limbs
    let a = R1Fr::from_raw([0xC5BF4F6AFF993D09, 0xA3361BDA67E62E0E, 0xAAAAAAAAAAAAAAAA, 0xFFFFFFFFE38E38E3]);
    let a_sqr = a.sqr();
    let expected = R1Fr::from_limbs([0x57abc6aa0349c084, 0x65b21b232a4cb7a5, 0x5ba781948b0fcd6e, 0xd6e9e0644bda12f7]);
    assert_eq!(a_sqr, expected, "secp256r1 Fr MontgomeryMulBigBug regression");
}

#[test]
fn secp256r1_fq_addition_subtraction_regression() {
    type R1Fq = Field<Secp256r1FqParams>;
    let mut fq1 = R1Fq::from_limbs([0xfffffe0000000200, 0x200fffff9ff, 0xfffffbfffffffe00, 0xfffffbff00000400]);
    let fq2 = R1Fq::from_limbs([0xfffffe0000000200, 0x200fffff9ff, 0xfffffbfffffffe00, 0xfffffbff00000400]);
    // fq1 += (modulus - 2) + 2 = modulus -> wraps to same value
    let mod_minus_2 = {
        let mut m = Secp256r1FqParams::MODULUS;
        // subtract 2 from m[0]
        if m[0] >= 2 {
            m[0] -= 2;
        } else {
            m[0] = m[0].wrapping_sub(2);
            // borrow propagation
            for i in 1..4 {
                if m[i] > 0 { m[i] -= 1; break; }
                m[i] = m[i].wrapping_sub(1);
            }
        }
        R1Fq::from_limbs(m)
    };
    fq1 = fq1 + mod_minus_2;
    fq1 = fq1 + R1Fq::from(2u64);
    // fq1 and fq2 should now be equal since adding modulus wraps
    let lhs = fq1 + fq1;
    let rhs = fq2 + fq2;
    assert_eq!(lhs, rhs, "addition regression: 2p overflow edge case");
    let fq3 = R1Fq::zero() - fq1;
    let fq4 = R1Fq::zero() - fq2;
    assert_eq!(fq3, fq4, "subtraction regression: 2p overflow edge case");
}

// --- BN254 Fq endomorphism split edge case ---

#[test]
fn bn254_fq_split_endomorphism_edge_case() {
    // k = 2^128 (data[2] = 1, rest zero)
    let k = Fq::from_raw([0, 0, 1, 0]);
    let (k1, k2) = k.split_into_endomorphism_scalars();

    // k1 and k2 should fit in 128 bits
    assert_eq!(k1.data[2], 0, "k1 upper limbs should be 0");
    assert_eq!(k1.data[3], 0, "k1 upper limbs should be 0");
    assert_eq!(k2.data[2], 0, "k2 upper limbs should be 0");
    assert_eq!(k2.data[3], 0, "k2 upper limbs should be 0");

    // Verify: k1 - k2*beta == k
    let k1_mont = k1.to_montgomery_form();
    let k2_mont = k2.to_montgomery_form();
    let beta = Fq::cube_root_of_unity();
    let result = (k1_mont - k2_mont * beta).from_montgomery_form();
    assert_eq!(result, k, "endomorphism split edge case: k1 - k2*beta should equal k");
}

// =========================================================================
// Missing group tests — ported from C++ g1.test.cpp / grumpkin.test.cpp
// =========================================================================

#[test]
fn bn254_eq_normalized() {
    let a = Element::<Bn254G1Params>::random_element();
    let a_norm = a.normalize();
    assert_eq!(a, a_norm, "a should equal a.normalize()");
}

#[test]
fn bn254_eq_infinity() {
    let inf1 = Element::<Bn254G1Params>::infinity();
    let inf2 = Element::<Bn254G1Params>::infinity();
    let a = Element::<Bn254G1Params>::random_element();
    assert_eq!(inf1, inf2, "infinity == infinity");
    assert_ne!(a, inf1, "point != infinity");
}

#[test]
fn bn254_add_exception_infinity() {
    let a = Element::<Bn254G1Params>::random_element();
    let neg_a = -a;
    let inf = Element::<Bn254G1Params>::infinity();

    // a + (-a) = infinity
    let result = a + neg_a;
    assert!(result.is_point_at_infinity(), "a + (-a) should be infinity");
    // inf + a = a
    assert_eq!(inf + a, a, "inf + a should be a");
    // a + inf = a
    assert_eq!(a + inf, a, "a + inf should be a");
}

#[test]
fn bn254_add_exception_dbl() {
    let a = Element::<Bn254G1Params>::random_element();
    let result = a + a;
    let expected = a.dbl();
    assert_eq!(result, expected, "a + a should equal a.dbl()");
}

#[test]
fn bn254_add_dbl_consistency() {
    let a = Element::<Bn254G1Params>::random_element();
    // Compute 8P via doubling chain: 2(2(2P))
    let d2 = a.dbl();
    let d4 = d2.dbl();
    let d8 = d4.dbl();
    // Compute 8P via add chain: 4P + 4P
    let d8_add = d4 + d4;
    assert_eq!(d8, d8_add, "doubling chain 8P should match add chain");
}

#[test]
fn bn254_add_dbl_consistency_repeated() {
    let a = Element::<Bn254G1Params>::random_element();
    let b = a.dbl(); // 2a
    let c = b + a;   // 3a
    let d = c + a;   // 4a
    let e = d + a;   // 5a
    let f = e + a;   // 6a
    let g = f + a;   // 7a
    let h = g + a;   // 8a
    let b2 = a.dbl();           // 2a
    let c2 = b2.dbl();          // 4a
    let d2 = c2.dbl();          // 8a
    let e2 = d2 - c2;           // 4a
    let f2 = e2 + b2;           // 6a
    let g2 = f2 - a;            // 5a
    let h2 = g2 - a.dbl().dbl();// 1a
    assert_eq!(b, b2, "2a mismatch");
    assert_eq!(d, c2, "4a mismatch");
    assert_eq!(h, d2, "8a mismatch");
    assert_eq!(d, e2, "4a via sub mismatch");
    assert_eq!(f, f2, "6a mismatch");
    assert_eq!(e, g2, "5a mismatch");
    assert_eq!(a, h2, "1a via sub mismatch");
}

#[test]
fn bn254_mixed_add_exception_infinity() {
    let a = Element::<Bn254G1Params>::random_element();
    let a_aff = a.to_affine();
    let neg_a_aff = (-a).to_affine();
    let inf = Element::<Bn254G1Params>::infinity();

    // a + (-a_affine) = infinity
    let result = a + neg_a_aff;
    assert!(result.is_point_at_infinity(), "a + (-a_affine) should be infinity");
    // inf + a_affine = a (as projective)
    let result2 = inf + a_aff;
    assert_eq!(result2.to_affine(), a_aff, "inf + a_affine should be a");
}

#[test]
fn bn254_mixed_add_exception_dbl() {
    let a = Element::<Bn254G1Params>::random_element();
    let a_aff = a.to_affine();
    let result = a + a_aff;
    let expected = a.dbl();
    assert_eq!(result, expected, "a + a_affine should equal a.dbl()");
}

#[test]
fn bn254_add_mixed_add_consistency_check() {
    let a = Element::<Bn254G1Params>::random_element();
    let b = Element::<Bn254G1Params>::random_element();
    let b_aff = b.to_affine();
    let via_proj = a + b;
    let via_mixed = a + b_aff;
    assert_eq!(via_proj, via_mixed, "projective add should match mixed add");
}

#[test]
fn bn254_group_exponentiation_zero_and_one() {
    let g = Element::<Bn254G1Params>::one();
    let zero = Fr::zero();
    let one_scalar = Fr::one();
    let result_zero = g.mul(&zero);
    let result_one = g.mul(&one_scalar);
    assert!(result_zero.is_point_at_infinity(), "G*0 should be infinity");
    assert_eq!(result_one, g, "G*1 should be G");
}

#[test]
fn bn254_group_exponentiation_consistency_check() {
    for _ in 0..10 {
        let a = Fr::random_element();
        let b = Fr::random_element();
        let g = Element::<Bn254G1Params>::one();
        let ga = g.mul(&a);
        let gab = ga.mul(&b);
        let ab = a * b;
        let g_ab = g.mul(&ab);
        assert_eq!(gab, g_ab, "(G*a)*b should equal G*(a*b)");
    }
}

#[test]
fn bn254_on_curve_random() {
    for _ in 0..100 {
        let p = Element::<Bn254G1Params>::random_element();
        assert!(p.on_curve(), "random element should be on curve");
    }
}

// --- Grumpkin group tests ---

#[test]
fn grumpkin_check_group_modulus() {
    // -1*G + 2G = G
    type GrumpkinFr = Field<Bn254FqParams>;
    let g = Element::<GrumpkinG1Params>::one();
    let minus_one = GrumpkinFr::one().negate();
    let neg_g = g.mul(&minus_one);
    let two_g = g.dbl();
    let result = neg_g + two_g;
    assert_eq!(result, g, "-1*G + 2G should equal G");
}

#[test]
fn grumpkin_check_b() {
    use crate::groups::curve_params::CurveParams;
    type GrumpkinFq = Field<Bn254FrParams>;
    let b = GrumpkinG1Params::coeff_b();
    let neg_17 = GrumpkinFq::from(17u64).negate();
    assert_eq!(b, neg_17, "Grumpkin coeff_b should be -17");
}

#[test]
fn grumpkin_mul_with_endomorphism_matches_basic() {
    type GrumpkinFr = Field<Bn254FqParams>;
    let g = Element::<GrumpkinG1Params>::one();
    for _ in 0..100 {
        let x1 = g.mul_without_endomorphism(&GrumpkinFr::random_element());
        let f1 = GrumpkinFr::random_element();
        let r1 = x1.mul_without_endomorphism(&f1);
        let r2 = x1.mul_with_endomorphism(&f1);
        assert_eq!(r1, r2, "Grumpkin endomorphism mul should match basic mul");
    }
}

#[test]
fn grumpkin_random_element() {
    let p = Element::<GrumpkinG1Params>::random_element();
    assert!(p.on_curve(), "Grumpkin random element should be on curve");
    assert!(!p.is_point_at_infinity(), "random element should not be infinity");
}

#[test]
fn grumpkin_random_affine_element() {
    let p = Element::<GrumpkinG1Params>::random_element().to_affine();
    assert!(p.on_curve(), "Grumpkin random affine element should be on curve");
}

#[test]
fn grumpkin_eq() {
    let a = Element::<GrumpkinG1Params>::random_element();
    let a_norm = a.normalize();
    assert_eq!(a, a_norm, "a should equal a.normalize()");
    let inf1 = Element::<GrumpkinG1Params>::infinity();
    let inf2 = Element::<GrumpkinG1Params>::infinity();
    assert_eq!(inf1, inf2, "infinity == infinity");
    assert_ne!(a, inf1, "point != infinity");
}

#[test]
fn grumpkin_add_exception_infinity() {
    let a = Element::<GrumpkinG1Params>::random_element();
    let neg_a = -a;
    let inf = Element::<GrumpkinG1Params>::infinity();
    let result = a + neg_a;
    assert!(result.is_point_at_infinity(), "a + (-a) should be infinity");
    assert_eq!(inf + a, a, "inf + a should be a");
    assert_eq!(a + inf, a, "a + inf should be a");
}

#[test]
fn grumpkin_add_exception_dbl() {
    let a = Element::<GrumpkinG1Params>::random_element();
    let result = a + a;
    let expected = a.dbl();
    assert_eq!(result, expected, "a + a should equal a.dbl()");
}

#[test]
fn grumpkin_add_dbl_consistency() {
    let a = Element::<GrumpkinG1Params>::random_element();
    let d2 = a.dbl();
    let d4 = d2.dbl();
    let d8 = d4.dbl();
    let d8_add = d4 + d4;
    assert_eq!(d8, d8_add, "Grumpkin doubling chain 8P should match add chain");
}

#[test]
fn grumpkin_add_dbl_consistency_repeated() {
    let a = Element::<GrumpkinG1Params>::random_element();
    let b = a.dbl();
    let c = b + a;
    let d = c + a;
    let e = d + a;
    let f = e + a;
    let g = f + a;
    let h = g + a;
    let b2 = a.dbl();
    let c2 = b2.dbl();
    let d2 = c2.dbl();
    let e2 = d2 - c2;
    let f2 = e2 + b2;
    let g2 = f2 - a;
    let h2 = g2 - a.dbl().dbl();
    assert_eq!(b, b2, "Grumpkin 2a mismatch");
    assert_eq!(d, c2, "Grumpkin 4a mismatch");
    assert_eq!(h, d2, "Grumpkin 8a mismatch");
    assert_eq!(d, e2, "Grumpkin 4a via sub mismatch");
    assert_eq!(f, f2, "Grumpkin 6a mismatch");
    assert_eq!(e, g2, "Grumpkin 5a mismatch");
    assert_eq!(a, h2, "Grumpkin 1a via sub mismatch");
}

#[test]
fn grumpkin_mixed_add_exception_infinity() {
    let a = Element::<GrumpkinG1Params>::random_element();
    let a_aff = a.to_affine();
    let neg_a_aff = (-a).to_affine();
    let inf = Element::<GrumpkinG1Params>::infinity();
    let result = a + neg_a_aff;
    assert!(result.is_point_at_infinity(), "Grumpkin a + (-a_affine) should be infinity");
    let result2 = inf + a_aff;
    assert_eq!(result2.to_affine(), a_aff, "Grumpkin inf + a_affine should be a");
}

#[test]
fn grumpkin_mixed_add_exception_dbl() {
    let a = Element::<GrumpkinG1Params>::random_element();
    let a_aff = a.to_affine();
    let result = a + a_aff;
    let expected = a.dbl();
    assert_eq!(result, expected, "Grumpkin a + a_affine should equal a.dbl()");
}

#[test]
fn grumpkin_add_mixed_add_consistency() {
    let a = Element::<GrumpkinG1Params>::random_element();
    let b = Element::<GrumpkinG1Params>::random_element();
    let b_aff = b.to_affine();
    let via_proj = a + b;
    let via_mixed = a + b_aff;
    assert_eq!(via_proj, via_mixed, "Grumpkin projective add should match mixed add");
}

#[test]
fn grumpkin_on_curve_random() {
    for _ in 0..100 {
        let p = Element::<GrumpkinG1Params>::random_element();
        assert!(p.on_curve(), "Grumpkin random element should be on curve");
    }
}

#[test]
fn grumpkin_batch_normalize() {
    type GrumpkinFr = Field<Bn254FqParams>;
    let g = Element::<GrumpkinG1Params>::one();
    let num_points = 2;
    let mut points: Vec<Element<GrumpkinG1Params>> = Vec::new();
    let mut normalized: Vec<Element<GrumpkinG1Params>> = Vec::new();
    for _ in 0..num_points {
        let a = g.mul_without_endomorphism(&GrumpkinFr::random_element());
        let b = g.mul_without_endomorphism(&GrumpkinFr::random_element());
        let point = a + b;
        points.push(point);
        normalized.push(point);
    }
    Element::batch_normalize(&mut normalized);
    for i in 0..num_points {
        let zz = points[i].z.sqr();
        let zzz = points[i].z * zz;
        let result_x = normalized[i].x * zz;
        let result_y = normalized[i].y * zzz;
        assert_eq!(result_x, points[i].x, "Grumpkin batch_normalize x mismatch at {}", i);
        assert_eq!(result_y, points[i].y, "Grumpkin batch_normalize y mismatch at {}", i);
    }
}

#[test]
fn grumpkin_group_exponentiation_consistency() {
    type GrumpkinFr = Field<Bn254FqParams>;
    for _ in 0..10 {
        let a = GrumpkinFr::random_element();
        let b = GrumpkinFr::random_element();
        let g = Element::<GrumpkinG1Params>::one();
        let ga = g.mul(&a);
        let gab = ga.mul(&b);
        let ab = a * b;
        let g_ab = g.mul(&ab);
        assert_eq!(gab, g_ab, "Grumpkin (G*a)*b should equal G*(a*b)");
    }
}

// --- secp256k1 group tests ---

#[test]
fn secp256k1_on_curve_random() {
    for _ in 0..100 {
        let p = Element::<Secp256k1G1Params>::random_element();
        assert!(p.on_curve(), "secp256k1 random element should be on curve");
    }
}

#[test]
fn secp256k1_add_dbl_consistency_repeated() {
    let a = Element::<Secp256k1G1Params>::random_element();
    let b = a.dbl();
    let c = b + a;
    let d = c + a;
    let e = d + a;
    let b2 = a.dbl();
    let c2 = b2.dbl();
    let d2 = c2 - b2;
    assert_eq!(b, b2, "secp256k1 2a mismatch");
    assert_eq!(d, c2, "secp256k1 4a mismatch");
    assert_eq!(b, d2, "secp256k1 4a-2a mismatch");
    // Chain longer: 5a via adds == dbl+add
    let e2 = c2 + a;
    assert_eq!(e, e2, "secp256k1 5a mismatch");
}

#[test]
fn secp256k1_group_exponentiation_zero_and_one() {
    let g = Element::<Secp256k1G1Params>::one();
    let zero = Secp256k1Fr::zero();
    let one_scalar = Secp256k1Fr::one();
    let result_zero = g.mul(&zero);
    let result_one = g.mul(&one_scalar);
    assert!(result_zero.is_point_at_infinity(), "secp256k1 G*0 should be infinity");
    assert_eq!(result_one, g, "secp256k1 G*1 should be G");
}

#[test]
fn secp256k1_group_exponentiation_consistency_check() {
    for _ in 0..10 {
        let a = Secp256k1Fr::random_element();
        let b = Secp256k1Fr::random_element();
        let g = Element::<Secp256k1G1Params>::one();
        let ga = g.mul(&a);
        let gab = ga.mul(&b);
        let ab = a * b;
        let g_ab = g.mul(&ab);
        assert_eq!(gab, g_ab, "secp256k1 (G*a)*b should equal G*(a*b)");
    }
}

// --- secp256r1 group tests ---

#[test]
fn secp256r1_check_group_modulus() {
    type R1Fr = Field<Secp256r1FrParams>;
    let g = Element::<Secp256r1G1Params>::one();
    let minus_one = R1Fr::one().negate();
    let neg_g = g.mul(&minus_one);
    let two_g = g.dbl();
    let result = neg_g + two_g;
    assert_eq!(result, g, "secp256r1 -1*G + 2G should equal G");
}

#[test]
fn secp256r1_on_curve_random() {
    for _ in 0..100 {
        let p = Element::<Secp256r1G1Params>::random_element();
        assert!(p.on_curve(), "secp256r1 random element should be on curve");
    }
}

#[test]
fn secp256r1_add_dbl_consistency_chain() {
    let a = Element::<Secp256r1G1Params>::random_element();
    let d2 = a.dbl();
    let d4 = d2.dbl();
    let d8 = d4.dbl();
    let d8_add = d4 + d4;
    assert_eq!(d8, d8_add, "secp256r1 doubling chain 8P should match add chain");
}

#[test]
fn secp256r1_add_dbl_consistency_repeated() {
    let a = Element::<Secp256r1G1Params>::random_element();
    let b = a.dbl();
    let c = b + a;
    let d = c + a;
    let e = d + a;
    let b2 = a.dbl();
    let c2 = b2.dbl();
    let e2 = c2 + a;
    assert_eq!(b, b2, "secp256r1 2a mismatch");
    assert_eq!(d, c2, "secp256r1 4a mismatch");
    assert_eq!(e, e2, "secp256r1 5a mismatch");
}

#[test]
fn secp256r1_group_exponentiation_zero_and_one() {
    type R1Fr = Field<Secp256r1FrParams>;
    let g = Element::<Secp256r1G1Params>::one();
    let zero = R1Fr::zero();
    let one_scalar = R1Fr::one();
    let result_zero = g.mul(&zero);
    let result_one = g.mul(&one_scalar);
    assert!(result_zero.is_point_at_infinity(), "secp256r1 G*0 should be infinity");
    assert_eq!(result_one, g, "secp256r1 G*1 should be G");
}

#[test]
fn secp256r1_group_exponentiation_consistency_check() {
    type R1Fr = Field<Secp256r1FrParams>;
    for _ in 0..10 {
        let a = R1Fr::random_element();
        let b = R1Fr::random_element();
        let g = Element::<Secp256r1G1Params>::one();
        let ga = g.mul(&a);
        let gab = ga.mul(&b);
        let ab = a * b;
        let g_ab = g.mul(&ab);
        assert_eq!(gab, g_ab, "secp256r1 (G*a)*b should equal G*(a*b)");
    }
}

// =========================================================================
// WNAF tests — ported from C++ wnaf.test.cpp
// =========================================================================

/// Helper to recover a scalar from WNAF encoding (matches C++ recover_fixed_wnaf)
fn recover_wnaf_scalar(table: &[u64], skew: bool, num_entries: usize, wnaf_bits: usize) -> [u64; 2] {
    let mut recovered: i128 = 0;
    for i in 0..num_entries {
        let entry = table[i];
        let val = (entry & 0x0fffffff) as i128;
        let sign = ((entry >> 31) & 1) != 0;
        let actual = (2 * val + 1) * if sign { -1 } else { 1 };
        recovered = recovered * (1i128 << wnaf_bits) + actual;
    }
    if skew {
        recovered -= 1;
    }
    [recovered as u64, (recovered >> 64) as u64]
}

#[test]
fn wnaf_zero() {
    use crate::groups::wnaf;
    let buffer: [u64; 2] = [0, 0];
    let wnaf_bits = 5;
    let wnaf_entries = (wnaf::SCALAR_BITS + wnaf_bits - 1) / wnaf_bits;
    let mut table = vec![0u64; wnaf_entries + 1];
    let mut skew = false;

    wnaf::fixed_wnaf(&buffer, &mut table, &mut skew, 0, 1, wnaf_bits);
    let recovered = recover_wnaf_scalar(&table, skew, wnaf_entries, wnaf_bits);
    assert_eq!(recovered[0], 0, "WNAF of zero: lo should be 0");
    assert_eq!(recovered[1], 0, "WNAF of zero: hi should be 0");
}

#[test]
fn wnaf_fixed_random() {
    use crate::groups::wnaf;
    let r = Fr::random_element();
    let raw = r.from_montgomery_form();
    let buffer: [u64; 2] = [raw.data[0], raw.data[1] & 0x7fffffffffffffff];
    let wnaf_bits = 5;
    let wnaf_entries = (wnaf::SCALAR_BITS + wnaf_bits - 1) / wnaf_bits;
    let mut table = vec![0u64; wnaf_entries + 1];
    let mut skew = false;

    wnaf::fixed_wnaf(&buffer, &mut table, &mut skew, 0, 1, wnaf_bits);
    let recovered = recover_wnaf_scalar(&table, skew, wnaf_entries, wnaf_bits);
    assert_eq!(recovered[0], buffer[0], "WNAF random: lo mismatch");
    assert_eq!(recovered[1], buffer[1], "WNAF random: hi mismatch");
}

#[test]
fn wnaf_fixed_simple_lo() {
    use crate::groups::wnaf;
    let buffer: [u64; 2] = [1, 0];
    let wnaf_bits = 5;
    let wnaf_entries = (wnaf::SCALAR_BITS + wnaf_bits - 1) / wnaf_bits;
    let mut table = vec![0u64; wnaf_entries + 1];
    let mut skew = false;

    wnaf::fixed_wnaf(&buffer, &mut table, &mut skew, 0, 1, wnaf_bits);
    let recovered = recover_wnaf_scalar(&table, skew, wnaf_entries, wnaf_bits);
    assert_eq!(recovered[0], 1, "WNAF simple lo: lo should be 1");
    assert_eq!(recovered[1], 0, "WNAF simple lo: hi should be 0");
}

#[test]
fn wnaf_fixed_simple_hi() {
    use crate::groups::wnaf;
    let buffer: [u64; 2] = [0, 1];
    let wnaf_bits = 5;
    let wnaf_entries = (wnaf::SCALAR_BITS + wnaf_bits - 1) / wnaf_bits;
    let mut table = vec![0u64; wnaf_entries + 1];
    let mut skew = false;

    wnaf::fixed_wnaf(&buffer, &mut table, &mut skew, 0, 1, wnaf_bits);
    let recovered = recover_wnaf_scalar(&table, skew, wnaf_entries, wnaf_bits);
    assert_eq!(recovered[0], 0, "WNAF simple hi: lo should be 0");
    assert_eq!(recovered[1], 1, "WNAF simple hi: hi should be 1");
}

// =========================================================================
// Pippenger MSM tests
// =========================================================================

#[test]
fn pippenger_empty() {
    let scalars: Vec<Fr> = vec![];
    let points: Vec<AffineElement<Bn254G1Params>> = vec![];
    let result = scalar_multiplication::pippenger_msm::<Bn254G1Params>(&scalars, &points);
    assert!(result.is_point_at_infinity(), "empty MSM should return infinity");
}

#[test]
fn pippenger_all_zeros() {
    let n = 10;
    let scalars: Vec<Fr> = vec![Fr::zero(); n];
    let points: Vec<AffineElement<Bn254G1Params>> = (0..n)
        .map(|_| Element::<Bn254G1Params>::random_element().to_affine())
        .collect();
    let result = scalar_multiplication::pippenger_msm::<Bn254G1Params>(&scalars, &points);
    assert!(result.is_point_at_infinity(), "all-zero scalars should return infinity");
}

#[test]
fn pippenger_single_point() {
    let scalar = Fr::random_element();
    let point = Element::<Bn254G1Params>::random_element().to_affine();
    let expected = Element::from_affine(&point).mul_without_endomorphism(&scalar);
    let result = scalar_multiplication::pippenger_msm::<Bn254G1Params>(&[scalar], &[point]);
    let expected_aff = expected.to_affine();
    let result_aff = result.to_affine();
    assert_eq!(expected_aff.x, result_aff.x, "single point x mismatch");
    assert_eq!(expected_aff.y, result_aff.y, "single point y mismatch");
}

#[test]
fn small_pippenger_correctness() {
    let n = 100;
    let scalars: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let points: Vec<AffineElement<Bn254G1Params>> = (0..n)
        .map(|_| Element::<Bn254G1Params>::random_element().to_affine())
        .collect();

    let expected = scalar_multiplication::naive_msm::<Bn254G1Params>(&scalars, &points);
    let result = scalar_multiplication::pippenger_msm::<Bn254G1Params>(&scalars, &points);
    let expected_aff = expected.to_affine();
    let result_aff = result.to_affine();
    assert_eq!(expected_aff.x, result_aff.x, "small pippenger x mismatch");
    assert_eq!(expected_aff.y, result_aff.y, "small pippenger y mismatch");
}

#[test]
fn pippenger_correctness() {
    let n = 1000;
    let scalars: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let points: Vec<AffineElement<Bn254G1Params>> = (0..n)
        .map(|_| Element::<Bn254G1Params>::random_element().to_affine())
        .collect();

    let expected = scalar_multiplication::naive_msm::<Bn254G1Params>(&scalars, &points);
    let result = scalar_multiplication::pippenger_msm::<Bn254G1Params>(&scalars, &points);
    let expected_aff = expected.to_affine();
    let result_aff = result.to_affine();
    assert_eq!(expected_aff.x, result_aff.x, "pippenger 1000 x mismatch");
    assert_eq!(expected_aff.y, result_aff.y, "pippenger 1000 y mismatch");
}

#[test]
fn pippenger_sparse_scalars() {
    let n = 200;
    let mut scalars: Vec<Fr> = (0..n).map(|_| Fr::random_element()).collect();
    let points: Vec<AffineElement<Bn254G1Params>> = (0..n)
        .map(|_| Element::<Bn254G1Params>::random_element().to_affine())
        .collect();
    // Zero out 80% of scalars
    for i in 0..n {
        if i % 5 != 0 {
            scalars[i] = Fr::zero();
        }
    }

    let expected = scalar_multiplication::naive_msm::<Bn254G1Params>(&scalars, &points);
    let result = scalar_multiplication::pippenger_msm::<Bn254G1Params>(&scalars, &points);
    let expected_aff = expected.to_affine();
    let result_aff = result.to_affine();
    assert_eq!(expected_aff.x, result_aff.x, "sparse pippenger x mismatch");
    assert_eq!(expected_aff.y, result_aff.y, "sparse pippenger y mismatch");
}

#[test]
fn pippenger_small_explicit() {
    // Deterministic 2-point test
    let g = AffineElement::<Bn254G1Params>::one();
    let g2 = (Element::from_affine(&g) + g).to_affine();
    let s1 = Fr::from_limbs([3, 0, 0, 0]);
    let s2 = Fr::from_limbs([7, 0, 0, 0]);
    // Expected: 3*G + 7*2G = 3G + 14G = 17G
    let expected = {
        let s17 = Fr::from_limbs([17, 0, 0, 0]);
        Element::from_affine(&g).mul_without_endomorphism(&s17)
    };

    let result = scalar_multiplication::pippenger_msm::<Bn254G1Params>(&[s1, s2], &[g, g2]);
    let expected_aff = expected.to_affine();
    let result_aff = result.to_affine();
    assert_eq!(expected_aff.x, result_aff.x, "explicit pippenger x mismatch");
    assert_eq!(expected_aff.y, result_aff.y, "explicit pippenger y mismatch");
}

#[test]
fn get_scalar_slice_consistency() {
    // Verify that extracting all slices and reconstructing yields the original scalar.
    // Round 0 = topmost bits, round (num_rounds-1) = lowest bits.
    for _ in 0..100 {
        let r = Fr::random_element();
        let raw = r.from_montgomery_form();
        let scalar = raw.data;
        let num_bits = 254; // BN254 Fr has 254 bits
        let slice_size = 7;
        let num_rounds = (num_bits + slice_size - 1) / slice_size;

        // Reconstruct: start from round 0 (MSB), shift left and OR in each slice
        let mut reconstructed = [0u64; 4];
        for round in 0..num_rounds {
            let slice = scalar_multiplication::get_scalar_slice_test(
                &scalar, round, slice_size, num_bits,
            ) as u64;

            // How many bits does this slice contribute?
            let hi_bit = num_bits - round * slice_size;
            let actual_bits = if hi_bit < slice_size { hi_bit } else { slice_size };

            // Left shift reconstructed by actual_bits
            let mut carry = 0u64;
            for limb in reconstructed.iter_mut() {
                let new_carry = if actual_bits < 64 {
                    *limb >> (64 - actual_bits)
                } else {
                    *limb
                };
                *limb = limb.wrapping_shl(actual_bits as u32) | carry;
                carry = new_carry;
            }
            reconstructed[0] |= slice;
        }

        assert_eq!(reconstructed[0], scalar[0], "slice consistency limb 0");
        assert_eq!(reconstructed[1], scalar[1], "slice consistency limb 1");
        assert_eq!(reconstructed[2], scalar[2], "slice consistency limb 2");
        assert_eq!(reconstructed[3], scalar[3], "slice consistency limb 3");
    }
}

#[test]
fn radix_sort_correctness() {
    let mut schedule: Vec<u64> = vec![
        (5u64 << 32) | 100,
        (3u64 << 32) | 50,
        (1u64 << 32) | 200,
        (7u64 << 32) | 0,
        (2u64 << 32) | 50,
        (4u64 << 32) | 100,
    ];
    let num_zero = scalar_multiplication::radix_sort_test(&mut schedule);

    // Verify sorted by lower 32 bits
    for i in 1..schedule.len() {
        let lo_prev = schedule[i - 1] & 0xFFFFFFFF;
        let lo_curr = schedule[i] & 0xFFFFFFFF;
        assert!(lo_prev <= lo_curr, "radix sort order violated at {}", i);
    }
    // Verify zero count
    assert_eq!(num_zero, 1, "should have exactly one zero-bucket entry");
    // Verify no entries lost
    assert_eq!(schedule.len(), 6, "sort should preserve entry count");
}

// =========================================================================
// Batched affine addition tests
// =========================================================================

#[test]
fn batched_affine_reduce_multiple_sequences() {
    let seq_len = 128;
    let num_seqs = 5;
    let mut points = Vec::new();
    let mut sequence_counts = Vec::new();
    let mut expected_sums = Vec::new();

    for _ in 0..num_seqs {
        let seq_points: Vec<AffineElement<Bn254G1Params>> = (0..seq_len)
            .map(|_| Element::<Bn254G1Params>::random_element().to_affine())
            .collect();

        // Compute expected sum naively
        let mut sum = Element::<Bn254G1Params>::infinity();
        for p in &seq_points {
            sum.add_assign_affine(p);
        }
        expected_sums.push(sum.to_affine());

        points.extend_from_slice(&seq_points);
        sequence_counts.push(seq_len);
    }

    let results = batched_affine_addition::batched_affine_add_in_place::<Bn254G1Params>(
        &mut points,
        &mut sequence_counts,
    );

    assert_eq!(results.len(), num_seqs);
    for (i, (result, expected)) in results.iter().zip(expected_sums.iter()).enumerate() {
        assert_eq!(result.x, expected.x, "seq {} x mismatch", i);
        assert_eq!(result.y, expected.y, "seq {} y mismatch", i);
    }
}

#[test]
fn batched_affine_single_sequence() {
    let seq_len = 64;
    let mut points: Vec<AffineElement<Bn254G1Params>> = (0..seq_len)
        .map(|_| Element::<Bn254G1Params>::random_element().to_affine())
        .collect();
    let mut sequence_counts = vec![seq_len];

    let mut expected = Element::<Bn254G1Params>::infinity();
    for p in points.iter() {
        expected.add_assign_affine(p);
    }
    let expected_aff = expected.to_affine();

    let results = batched_affine_addition::batched_affine_add_in_place::<Bn254G1Params>(
        &mut points,
        &mut sequence_counts,
    );

    assert_eq!(results.len(), 1);
    assert_eq!(results[0].x, expected_aff.x, "single seq x mismatch");
    assert_eq!(results[0].y, expected_aff.y, "single seq y mismatch");
}

#[test]
fn batched_affine_odd_lengths() {
    let lengths = [3, 5, 7, 11, 13];
    let mut points = Vec::new();
    let mut sequence_counts = Vec::new();
    let mut expected_sums = Vec::new();

    for &len in &lengths {
        let seq_points: Vec<AffineElement<Bn254G1Params>> = (0..len)
            .map(|_| Element::<Bn254G1Params>::random_element().to_affine())
            .collect();

        let mut sum = Element::<Bn254G1Params>::infinity();
        for p in &seq_points {
            sum.add_assign_affine(p);
        }
        expected_sums.push(sum.to_affine());

        points.extend_from_slice(&seq_points);
        sequence_counts.push(len);
    }

    let results = batched_affine_addition::batched_affine_add_in_place::<Bn254G1Params>(
        &mut points,
        &mut sequence_counts,
    );

    assert_eq!(results.len(), lengths.len());
    for (i, (result, expected)) in results.iter().zip(expected_sums.iter()).enumerate() {
        assert_eq!(result.x, expected.x, "odd seq {} x mismatch", i);
        assert_eq!(result.y, expected.y, "odd seq {} y mismatch", i);
    }
}

// =========================================================================
// G2 group tests — BN254 G2 over Fq2
// =========================================================================

#[test]
fn g2_random_element() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let mut p = g.dbl();
    p.add_assign_element(&g);
    let aff = p.to_affine();
    assert!(!aff.is_point_at_infinity(), "non-trivial G2 point should not be infinity");
}

#[test]
fn g2_eq() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let a = g.dbl();
    let mut b = g;
    b.add_assign_element(&g);
    assert_eq!(a.to_affine().x, b.to_affine().x, "G2 dbl should equal add self");
    assert_eq!(a.to_affine().y, b.to_affine().y);
}

#[test]
fn g2_add_exception_infinity() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let neg_g = -g;
    let mut result = g;
    result.add_assign_element(&neg_g);
    assert!(result.is_point_at_infinity(), "G2 g + (-g) should be infinity");
}

#[test]
fn g2_add_exception_dbl() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let mut result = g;
    result.add_assign_element(&g);
    let dbl = g.dbl();
    assert_eq!(result.to_affine().x, dbl.to_affine().x, "G2 add self should equal dbl");
    assert_eq!(result.to_affine().y, dbl.to_affine().y);
}

#[test]
fn g2_add_dbl_consistency() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let d2 = g.dbl();
    let d4 = d2.dbl();
    let d8 = d4.dbl();
    let mut d8_add = d4;
    d8_add.add_assign_element(&d4);
    assert_eq!(d8.to_affine().x, d8_add.to_affine().x, "G2 8P dbl vs add");
    assert_eq!(d8.to_affine().y, d8_add.to_affine().y);
}

#[test]
fn g2_add_dbl_consistency_repeated() {
    let a = G2Element::from_affine(&G2AffineElement::generator());
    let b = a.dbl();
    let mut c = b;
    c.add_assign_element(&a);
    let mut d = c;
    d.add_assign_element(&a);
    let b2 = a.dbl();
    let c2 = b2.dbl();
    assert_eq!(b.to_affine().x, b2.to_affine().x, "G2 2a");
    assert_eq!(d.to_affine().x, c2.to_affine().x, "G2 4a");
}

#[test]
fn g2_add_infinity_identity() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let inf = G2Element::infinity();
    let mut r1 = g;
    r1.add_assign_element(&inf);
    assert_eq!(r1.to_affine().x, g.to_affine().x, "G2 g + inf = g");
    let mut r2 = inf;
    r2.add_assign_element(&g);
    assert_eq!(r2.to_affine().x, g.to_affine().x, "G2 inf + g = g");
}

#[test]
fn g2_to_affine_roundtrip() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let doubled = g.dbl();
    let aff = doubled.to_affine();
    let back = G2Element::from_affine(&aff);
    assert_eq!(back.to_affine().x, aff.x, "G2 affine roundtrip x");
    assert_eq!(back.to_affine().y, aff.y, "G2 affine roundtrip y");
}

#[test]
fn g2_infinity_is_infinity() {
    let inf = G2Element::infinity();
    assert!(inf.is_point_at_infinity());
    let inf_aff = G2AffineElement::infinity();
    assert!(inf_aff.is_point_at_infinity());
}

#[test]
fn g2_negation() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let neg_g = -g;
    assert_eq!(neg_g.x, g.x, "negation should preserve x");
    assert_eq!(neg_g.y, -g.y, "negation should negate y");
}

#[test]
fn g2_mul_scalar_small() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let three_g = g.mul_scalar(&Fr::from(3u64));
    let mut expected = g.dbl();
    expected.add_assign_element(&g);
    assert_eq!(three_g.to_affine().x, expected.to_affine().x, "3*G via scalar_mul");
    assert_eq!(three_g.to_affine().y, expected.to_affine().y);
}

#[test]
fn g2_mul_scalar_zero() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let result = g.mul_scalar(&Fr::zero());
    assert!(result.is_point_at_infinity(), "0*G should be infinity");
}

#[test]
fn g2_mul_scalar_one() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let result = g.mul_scalar(&Fr::one());
    assert_eq!(result.to_affine().x, g.to_affine().x, "1*G should be G");
    assert_eq!(result.to_affine().y, g.to_affine().y);
}

#[test]
fn g2_exponentiation_consistency() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let a = Fr::from(5u64);
    let b = Fr::from(7u64);
    let ga = g.mul_scalar(&a);
    let gab = ga.mul_scalar(&b);
    let ab = a * b; // 35
    let g_ab = g.mul_scalar(&ab);
    assert_eq!(gab.to_affine().x, g_ab.to_affine().x, "G2 (G*a)*b = G*(a*b)");
    assert_eq!(gab.to_affine().y, g_ab.to_affine().y);
}

#[test]
fn g2_on_curve_generator() {
    // Verify generator satisfies y^2 = x^3 + b'
    let g2_gen = G2AffineElement::generator();
    let b_twist = Fq2::twist_coeff_b();
    let y_sq = g2_gen.y.sqr();
    let x_cubed = g2_gen.x.sqr() * g2_gen.x;
    let rhs = x_cubed + b_twist;
    assert_eq!(y_sq, rhs, "G2 generator should satisfy curve equation");
}

// =========================================================================
// Additional field tests — Fr, Fq, secp256k1, secp256r1
// =========================================================================

#[test]
fn bn254_fr_multiplicative_generator() {
    // The primitive root (multiplicative generator) should have multiplicative order p-1
    // Verify it's not 0 or 1
    let g = Fr::from_raw(Bn254FrParams::PRIMITIVE_ROOT);
    assert!(!g.is_zero(), "multiplicative generator should not be zero");
    let mont = g.to_montgomery_form();
    assert_ne!(mont, Fr::one(), "multiplicative generator should not be one");
}

#[test]
fn bn254_fr_uint256_conversions() {
    let a = Fr::random_element();
    let raw = a.from_montgomery_form();
    let back = raw.to_montgomery_form();
    assert_eq!(a, back, "Fr to/from raw roundtrip should be consistent");
}

#[test]
fn bn254_fr_equivalent_randomness() {
    // Two random elements should not be equal (probabilistically)
    let a = Fr::random_element();
    let b = Fr::random_element();
    assert_ne!(a, b, "two random Fr elements should differ");
}

#[test]
fn bn254_fq_multiplicative_generator() {
    // Verify Fq cube root is not trivial
    let cube_root = Fq::cube_root_of_unity();
    assert_ne!(cube_root, Fq::one());
    assert_eq!(cube_root * cube_root * cube_root, Fq::one());
}

#[test]
fn bn254_fq_serialize_to_buffer() {
    let a = Fq::random_element();
    let bytes = a.to_be_bytes();
    let recovered = Fq::from_be_bytes(&bytes);
    assert_eq!(a, recovered, "Fq serialization roundtrip");
}

#[test]
fn bn254_fq_equivalent_randomness() {
    let a = Fq::random_element();
    let b = Fq::random_element();
    assert_ne!(a, b, "two random Fq elements should differ");
}

#[test]
fn bn254_fq_neg_and_self_neg_zero() {
    let zero = Fq::zero();
    let neg_zero = zero.negate();
    assert_eq!(zero, neg_zero, "-0 should equal 0 for Fq");
}

#[test]
fn bn254_fq_pow_regression() {
    let a = Fq::from(7u64);
    let a_cubed = a.pow(&[3, 0, 0, 0]);
    assert_eq!(a_cubed, Fq::from(343u64), "7^3 should be 343");
}

#[test]
fn bn254_fq_sqr_regression() {
    let a = Fq::from(12345u64);
    assert_eq!(a.sqr(), a * a, "sqr should match mul for regression value");
}

// --- secp256k1 additional field tests ---

#[test]
fn secp256k1_fq_eq() {
    type K1Fq = Field<Secp256k1FqParams>;
    let a = K1Fq::from(42u64);
    let b = K1Fq::from(42u64);
    assert_eq!(a, b);
    let c = K1Fq::from(43u64);
    assert_ne!(a, c);
}

#[test]
fn secp256k1_fq_is_zero() {
    type K1Fq = Field<Secp256k1FqParams>;
    assert!(K1Fq::zero().is_zero());
    assert!(!K1Fq::one().is_zero());
}

#[test]
fn secp256k1_fq_random_element() {
    type K1Fq = Field<Secp256k1FqParams>;
    let a = K1Fq::random_element();
    let b = K1Fq::random_element();
    assert_ne!(a, b);
}

#[test]
fn secp256k1_fq_sqr_equals_mul() {
    type K1Fq = Field<Secp256k1FqParams>;
    for _ in 0..100 {
        let a = K1Fq::random_element();
        assert_eq!(a.sqr(), a * a);
    }
}

#[test]
fn secp256k1_fq_add_mul_consistency() {
    type K1Fq = Field<Secp256k1FqParams>;
    for _ in 0..100 {
        let a = K1Fq::random_element();
        let b = K1Fq::random_element();
        let c = K1Fq::random_element();
        assert_eq!((a + b) * c, a * c + b * c);
    }
}

#[test]
fn secp256k1_fq_sub_mul_consistency() {
    type K1Fq = Field<Secp256k1FqParams>;
    for _ in 0..100 {
        let a = K1Fq::random_element();
        let b = K1Fq::random_element();
        let c = K1Fq::random_element();
        assert_eq!((a - b) * c, a * c - b * c);
    }
}

#[test]
fn secp256k1_fq_invert() {
    type K1Fq = Field<Secp256k1FqParams>;
    for _ in 0..100 {
        let a = K1Fq::random_element();
        if !a.is_zero() {
            assert_eq!(a * a.invert(), K1Fq::one());
        }
    }
}

#[test]
fn secp256k1_fq_sqrt() {
    type K1Fq = Field<Secp256k1FqParams>;
    for _ in 0..100 {
        let a = K1Fq::random_element();
        let a_sq = a.sqr();
        let (found, root) = a_sq.sqrt();
        assert!(found);
        assert!(root == a || root == a.negate());
    }
}

#[test]
fn secp256k1_fr_eq() {
    let a = Secp256k1Fr::from(42u64);
    let b = Secp256k1Fr::from(42u64);
    assert_eq!(a, b);
}

#[test]
fn secp256k1_fr_add_mul_consistency() {
    for _ in 0..100 {
        let a = Secp256k1Fr::random_element();
        let b = Secp256k1Fr::random_element();
        let c = Secp256k1Fr::random_element();
        assert_eq!((a + b) * c, a * c + b * c);
    }
}

#[test]
fn secp256k1_fr_sub_mul_consistency() {
    for _ in 0..100 {
        let a = Secp256k1Fr::random_element();
        let b = Secp256k1Fr::random_element();
        let c = Secp256k1Fr::random_element();
        assert_eq!((a - b) * c, a * c - b * c);
    }
}

#[test]
fn secp256k1_fr_invert() {
    for _ in 0..100 {
        let a = Secp256k1Fr::random_element();
        if !a.is_zero() {
            assert_eq!(a * a.invert(), Secp256k1Fr::one());
        }
    }
}

// --- secp256r1 additional field tests ---

#[test]
fn secp256r1_fq_eq() {
    type R1Fq = Field<Secp256r1FqParams>;
    let a = R1Fq::from(42u64);
    let b = R1Fq::from(42u64);
    assert_eq!(a, b);
}

#[test]
fn secp256r1_fq_is_zero() {
    type R1Fq = Field<Secp256r1FqParams>;
    assert!(R1Fq::zero().is_zero());
    assert!(!R1Fq::one().is_zero());
}

#[test]
fn secp256r1_fq_random_element() {
    type R1Fq = Field<Secp256r1FqParams>;
    let a = R1Fq::random_element();
    let b = R1Fq::random_element();
    assert_ne!(a, b);
}

#[test]
fn secp256r1_fq_sqr_consistency() {
    type R1Fq = Field<Secp256r1FqParams>;
    for _ in 0..100 {
        let a = R1Fq::random_element();
        assert_eq!(a.sqr(), a * a);
    }
}

#[test]
fn secp256r1_fq_add_mul_consistency() {
    type R1Fq = Field<Secp256r1FqParams>;
    for _ in 0..100 {
        let a = R1Fq::random_element();
        let b = R1Fq::random_element();
        let c = R1Fq::random_element();
        assert_eq!((a + b) * c, a * c + b * c);
    }
}

#[test]
fn secp256r1_fq_sub_mul_consistency() {
    type R1Fq = Field<Secp256r1FqParams>;
    for _ in 0..100 {
        let a = R1Fq::random_element();
        let b = R1Fq::random_element();
        let c = R1Fq::random_element();
        assert_eq!((a - b) * c, a * c - b * c);
    }
}

#[test]
fn secp256r1_fq_invert() {
    type R1Fq = Field<Secp256r1FqParams>;
    for _ in 0..100 {
        let a = R1Fq::random_element();
        if !a.is_zero() {
            assert_eq!(a * a.invert(), R1Fq::one());
        }
    }
}

#[test]
fn secp256r1_fq_sqrt() {
    type R1Fq = Field<Secp256r1FqParams>;
    for _ in 0..100 {
        let a = R1Fq::random_element();
        let a_sq = a.sqr();
        let (found, root) = a_sq.sqrt();
        assert!(found);
        assert!(root == a || root == a.negate());
    }
}

#[test]
fn secp256r1_fr_eq() {
    type R1Fr = Field<Secp256r1FrParams>;
    let a = R1Fr::from(42u64);
    let b = R1Fr::from(42u64);
    assert_eq!(a, b);
}

#[test]
fn secp256r1_fr_add_mul_consistency() {
    type R1Fr = Field<Secp256r1FrParams>;
    for _ in 0..100 {
        let a = R1Fr::random_element();
        let b = R1Fr::random_element();
        let c = R1Fr::random_element();
        assert_eq!((a + b) * c, a * c + b * c);
    }
}

#[test]
fn secp256r1_fr_sub_mul_consistency() {
    type R1Fr = Field<Secp256r1FrParams>;
    for _ in 0..100 {
        let a = R1Fr::random_element();
        let b = R1Fr::random_element();
        let c = R1Fr::random_element();
        assert_eq!((a - b) * c, a * c - b * c);
    }
}

#[test]
fn secp256r1_fr_invert() {
    type R1Fr = Field<Secp256r1FrParams>;
    for _ in 0..100 {
        let a = R1Fr::random_element();
        if !a.is_zero() {
            assert_eq!(a * a.invert(), R1Fr::one());
        }
    }
}

#[test]
fn secp256r1_fr_montgomery_roundtrip() {
    type R1Fr = Field<Secp256r1FrParams>;
    for _ in 0..100 {
        let a = R1Fr::random_element();
        let roundtrip = a.from_montgomery_form().to_montgomery_form();
        assert_eq!(a, roundtrip);
    }
}

// =========================================================================
// Fq field tests — covers MulShortIntegers, CoarseEquivalenceChecks,
// Sqrt (deterministic), OneAndZero, SplitIntoEndomorphismScalarsSimple, RInv
// =========================================================================

#[test]
fn bn254_fq_mul_short_integers() {
    for a in 1u64..=10 {
        for b in 1u64..=10 {
            let fa = Fq::from(a);
            let fb = Fq::from(b);
            let fc = Fq::from(a * b);
            assert_eq!(fa * fb, fc, "Fq: {} * {} should be {}", a, b, a * b);
        }
    }
}

#[test]
fn bn254_fq_coarse_equivalence() {
    for _ in 0..100 {
        let a = Fq::random_element();
        let b = Fq::random_element();
        let c = Fq::random_element();
        // Distributive: a*(b+c) == a*b + a*c
        assert_eq!(a * (b + c), a * b + a * c, "Fq distributive law");
        // Difference of squares: (a-b)*(a+b) == a^2 - b^2
        assert_eq!((a - b) * (a + b), a.sqr() - b.sqr(), "Fq difference of squares");
    }
}

#[test]
fn bn254_fq_sqrt_deterministic() {
    let forty_nine = Fq::from(49u64);
    let (is_qr, root) = forty_nine.sqrt();
    assert!(is_qr, "49 should be a quadratic residue");
    assert_eq!(root.sqr(), forty_nine, "sqrt(49)^2 should be 49");
    // Verify root is either 7 or -7
    let seven = Fq::from(7u64);
    assert!(root == seven || root == seven.negate(), "sqrt(49) should be +/- 7");
}

#[test]
fn bn254_fq_one_and_zero() {
    let x = Fq::random_element();
    assert_eq!(Fq::one() * x, x, "one * x should be x");
    assert_eq!(x * Fq::one(), x, "x * one should be x");
    assert_eq!(Fq::zero() + x, x, "zero + x should be x");
    assert_eq!(x + Fq::zero(), x, "x + zero should be x");
    assert_eq!(x - x, Fq::zero(), "x - x should be zero");
}

#[test]
fn bn254_fq_split_endo_simple() {
    // C++ SplitIntoEndomorphismScalarsSimple: input = {1, 0, 0, 0}
    let k = Fq::from_raw([1, 0, 0, 0]);
    let (k1, k2) = k.split_into_endomorphism_scalars();
    let k1_mont = k1.to_montgomery_form();
    let k2_mont = k2.to_montgomery_form();
    let beta = Fq::cube_root_of_unity();
    let result = (k1_mont - k2_mont * beta).from_montgomery_form();
    assert_eq!(result, k, "Fq endo split: k1 - k2*beta should equal k");
}

#[test]
fn bn254_fq_r_inv() {
    // Montgomery invariant: MODULUS[0] * R_INV + 1 == 0 (mod 2^64)
    let r_inv = Bn254FqParams::R_INV;
    let mod_lo = Bn254FqParams::MODULUS[0];
    let product = mod_lo.wrapping_mul(r_inv);
    assert_eq!(product.wrapping_add(1), 0, "MODULUS[0] * R_INV + 1 should be 0 mod 2^64");
}

// =========================================================================
// G1 group tests — AddAffineTest, OperatorOrdering, Serialize, InitializationCheck
// =========================================================================

#[test]
fn bn254_add_affine_test() {
    let a = Element::<Bn254G1Params>::random_element().to_affine();
    let b = Element::<Bn254G1Params>::random_element().to_affine();
    let result = Element::from_affine(&a) + Element::from_affine(&b);
    assert!(result.on_curve(), "a + b should be on curve");
    assert!(!result.is_point_at_infinity());
}

#[test]
fn bn254_operator_ordering() {
    // Different projective representations of same affine point compare equal after normalization
    let a = Element::<Bn254G1Params>::random_element();
    let a_aff = a.to_affine();
    let b = Element::from_affine(&a_aff);
    // a and b represent the same affine point but may have different Z coords
    let a_norm = a.normalize();
    let b_norm = b.normalize();
    assert_eq!(a_norm.to_affine().x, b_norm.to_affine().x, "normalized x should match");
    assert_eq!(a_norm.to_affine().y, b_norm.to_affine().y, "normalized y should match");
}

#[test]
fn bn254_serialize_affine() {
    // Roundtrip test
    let p = Element::<Bn254G1Params>::random_element().to_affine();
    let buf = p.to_buffer();
    assert_eq!(buf.len(), 64);
    let x_bytes: &[u8; 32] = buf[..32].try_into().unwrap();
    let y_bytes: &[u8; 32] = buf[32..].try_into().unwrap();
    let x_recovered = Fq::from_be_bytes(x_bytes);
    let y_recovered = Fq::from_be_bytes(y_bytes);
    assert_eq!(p.x, x_recovered, "serialized x should roundtrip");
    assert_eq!(p.y, y_recovered, "serialized y should roundtrip");

    // Infinity serialization: all 0xFF
    let inf = AffineElement::<Bn254G1Params>::infinity();
    let buf_inf = inf.to_buffer();
    assert_eq!(buf_inf, [0xFF; 64], "infinity should serialize to all 0xFF");
}

#[test]
fn bn254_initialization_check() {
    use crate::groups::curve_params::CurveParams;
    let g = Element::<Bn254G1Params>::one();
    assert!(g.on_curve(), "generator should be on curve");
    assert!(!g.is_point_at_infinity(), "generator should not be infinity");
    let g_aff = g.to_affine();
    let expected_x = Bn254G1Params::generator_x();
    let expected_y = Bn254G1Params::generator_y();
    assert_eq!(g_aff.x, expected_x, "generator x should match generator_x()");
    assert_eq!(g_aff.y, expected_y, "generator y should match generator_y()");
}

// =========================================================================
// Grumpkin — GroupExponentiationZeroAndOne
// =========================================================================

#[test]
fn grumpkin_group_exponentiation_zero_and_one() {
    type GrumpkinFr = Field<Bn254FqParams>;
    let g = Element::<GrumpkinG1Params>::one();
    let zero = GrumpkinFr::zero();
    let one_scalar = GrumpkinFr::one();
    let result_zero = g.mul(&zero);
    let result_one = g.mul(&one_scalar);
    assert!(result_zero.is_point_at_infinity(), "Grumpkin G*0 should be infinity");
    assert_eq!(result_one, g, "Grumpkin G*1 should be G");
}

// =========================================================================
// G2 tests — RandomAffineElement, DblCheck, ExponentiationZeroAndOne
// =========================================================================

#[test]
fn g2_random_affine_element() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    // Create "random" G2 by scalar mul
    let scalar = Fr::random_element();
    let p = g.mul_scalar(&scalar);
    let aff = p.to_affine();
    assert!(!aff.is_point_at_infinity(), "random G2 should not be infinity");
    // Verify on curve: y^2 == x^3 + b_twist
    let b_twist = Fq2::twist_coeff_b();
    let y_sq = aff.y.sqr();
    let x_cubed = aff.x.sqr() * aff.x;
    let rhs = x_cubed + b_twist;
    assert_eq!(y_sq, rhs, "random G2 affine should be on curve");
}

#[test]
fn g2_dbl_check() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let doubled = g.dbl();
    let g_aff = g.to_affine();
    let d_aff = doubled.to_affine();
    assert_ne!(g_aff.x, d_aff.x, "doubled G2 should differ from generator");
    // Verify doubled point on curve
    let b_twist = Fq2::twist_coeff_b();
    let y_sq = d_aff.y.sqr();
    let x_cubed = d_aff.x.sqr() * d_aff.x;
    let rhs = x_cubed + b_twist;
    assert_eq!(y_sq, rhs, "doubled G2 should be on curve");
}

#[test]
fn g2_exponentiation_zero_and_one() {
    let g = G2Element::from_affine(&G2AffineElement::generator());
    let result_zero = g.mul_scalar(&Fr::zero());
    let result_one = g.mul_scalar(&Fr::one());
    assert!(result_zero.is_point_at_infinity(), "G2 G*0 should be infinity");
    assert_eq!(result_one.to_affine().x, g.to_affine().x, "G2 G*1 should be G (x)");
    assert_eq!(result_one.to_affine().y, g.to_affine().y, "G2 G*1 should be G (y)");
}

// =========================================================================
// secp256k1 group tests — TestArithmetic, CheckGroupModulus,
// AddExceptionInfinity, AddExceptionDbl, MixedAddExceptionInfinity,
// MixedAddExceptionDbl, AddMixedAddConsistency, BatchNormalize
// =========================================================================

#[test]
fn secp256k1_test_arithmetic() {
    let a = Element::<Secp256k1G1Params>::random_element();
    let b = Element::<Secp256k1G1Params>::random_element();
    // a + b - b == a
    let result = (a + b) + (-b);
    assert_eq!(result, a, "secp256k1 a + b - b should equal a");
    // a * 2 == a.dbl()
    assert_eq!(a.scalar_mul(&[2, 0, 0, 0]), a.dbl(), "secp256k1 2*a should equal a.dbl()");
}

#[test]
fn secp256k1_check_group_modulus() {
    // Verify MODULUS matches known secp256k1 order
    let order = Secp256k1FrParams::MODULUS;
    assert_eq!(order[0], 0xBFD25E8CD0364141);
    assert_eq!(order[1], 0xBAAEDCE6AF48A03B);
    assert_eq!(order[2], 0xFFFFFFFFFFFFFFFE);
    assert_eq!(order[3], 0xFFFFFFFFFFFFFFFF);
}

#[test]
fn secp256k1_add_exception_infinity() {
    let a = Element::<Secp256k1G1Params>::random_element();
    let neg_a = -a;
    let inf = Element::<Secp256k1G1Params>::infinity();
    assert!(
        (a + neg_a).is_point_at_infinity(),
        "secp256k1 a + (-a) should be infinity"
    );
    assert_eq!(inf + a, a, "secp256k1 inf + a should be a");
    assert_eq!(a + inf, a, "secp256k1 a + inf should be a");
}

#[test]
fn secp256k1_add_exception_dbl() {
    let a = Element::<Secp256k1G1Params>::random_element();
    assert_eq!(a + a, a.dbl(), "secp256k1 a + a should equal a.dbl()");
}

#[test]
fn secp256k1_mixed_add_exception_infinity() {
    let a = Element::<Secp256k1G1Params>::random_element();
    let a_aff = a.to_affine();
    let neg_a_aff = (-a).to_affine();
    let inf = Element::<Secp256k1G1Params>::infinity();
    assert!(
        (a + neg_a_aff).is_point_at_infinity(),
        "secp256k1 a + (-a_affine) should be infinity"
    );
    assert_eq!(
        (inf + a_aff).to_affine(),
        a_aff,
        "secp256k1 inf + a_affine should be a"
    );
}

#[test]
fn secp256k1_mixed_add_exception_dbl() {
    let a = Element::<Secp256k1G1Params>::random_element();
    let a_aff = a.to_affine();
    assert_eq!(a + a_aff, a.dbl(), "secp256k1 a + a_affine should equal a.dbl()");
}

#[test]
fn secp256k1_add_mixed_add_consistency() {
    let a = Element::<Secp256k1G1Params>::random_element();
    let b = Element::<Secp256k1G1Params>::random_element();
    let b_aff = b.to_affine();
    assert_eq!(a + b, a + b_aff, "secp256k1 projective add should match mixed add");
}

#[test]
fn secp256k1_batch_normalize() {
    let g = Element::<Secp256k1G1Params>::one();
    let num_points = 2;
    let mut points: Vec<Element<Secp256k1G1Params>> = Vec::new();
    let mut normalized: Vec<Element<Secp256k1G1Params>> = Vec::new();
    for _ in 0..num_points {
        let a = g.mul_without_endomorphism(&Secp256k1Fr::random_element());
        let b = g.mul_without_endomorphism(&Secp256k1Fr::random_element());
        let point = a + b;
        points.push(point);
        normalized.push(point);
    }
    Element::batch_normalize(&mut normalized);
    for i in 0..num_points {
        let zz = points[i].z.sqr();
        let zzz = points[i].z * zz;
        assert_eq!(
            normalized[i].x * zz,
            points[i].x,
            "secp256k1 batch_normalize x mismatch at {}",
            i
        );
        assert_eq!(
            normalized[i].y * zzz,
            points[i].y,
            "secp256k1 batch_normalize y mismatch at {}",
            i
        );
    }
}

// =========================================================================
// secp256r1 group tests — TestArithmetic, AddExceptionInfinity,
// AddExceptionDbl, MixedAddExceptionInfinity, MixedAddExceptionDbl,
// AddMixedAddConsistency, BatchNormalize
// =========================================================================

#[test]
fn secp256r1_test_arithmetic() {
    let a = Element::<Secp256r1G1Params>::random_element();
    let b = Element::<Secp256r1G1Params>::random_element();
    // a + b - b == a
    let result = (a + b) + (-b);
    assert_eq!(result, a, "secp256r1 a + b - b should equal a");
    // a * 2 == a.dbl()
    assert_eq!(
        a.scalar_mul(&[2, 0, 0, 0]),
        a.dbl(),
        "secp256r1 2*a should equal a.dbl()"
    );
}

#[test]
fn secp256r1_add_exception_infinity() {
    let a = Element::<Secp256r1G1Params>::random_element();
    let neg_a = -a;
    let inf = Element::<Secp256r1G1Params>::infinity();
    assert!(
        (a + neg_a).is_point_at_infinity(),
        "secp256r1 a + (-a) should be infinity"
    );
    assert_eq!(inf + a, a, "secp256r1 inf + a should be a");
    assert_eq!(a + inf, a, "secp256r1 a + inf should be a");
}

#[test]
fn secp256r1_add_exception_dbl() {
    let a = Element::<Secp256r1G1Params>::random_element();
    assert_eq!(a + a, a.dbl(), "secp256r1 a + a should equal a.dbl()");
}

#[test]
fn secp256r1_mixed_add_exception_infinity() {
    let a = Element::<Secp256r1G1Params>::random_element();
    let a_aff = a.to_affine();
    let neg_a_aff = (-a).to_affine();
    let inf = Element::<Secp256r1G1Params>::infinity();
    assert!(
        (a + neg_a_aff).is_point_at_infinity(),
        "secp256r1 a + (-a_affine) should be infinity"
    );
    assert_eq!(
        (inf + a_aff).to_affine(),
        a_aff,
        "secp256r1 inf + a_affine should be a"
    );
}

#[test]
fn secp256r1_mixed_add_exception_dbl() {
    let a = Element::<Secp256r1G1Params>::random_element();
    let a_aff = a.to_affine();
    assert_eq!(a + a_aff, a.dbl(), "secp256r1 a + a_affine should equal a.dbl()");
}

#[test]
fn secp256r1_add_mixed_add_consistency() {
    let a = Element::<Secp256r1G1Params>::random_element();
    let b = Element::<Secp256r1G1Params>::random_element();
    let b_aff = b.to_affine();
    assert_eq!(a + b, a + b_aff, "secp256r1 projective add should match mixed add");
}

#[test]
fn secp256r1_batch_normalize() {
    let g = Element::<Secp256r1G1Params>::one();
    let num_points = 2;
    let mut points: Vec<Element<Secp256r1G1Params>> = Vec::new();
    let mut normalized: Vec<Element<Secp256r1G1Params>> = Vec::new();
    for _ in 0..num_points {
        let a = g.mul_without_endomorphism(&Secp256r1Fr::random_element());
        let b = g.mul_without_endomorphism(&Secp256r1Fr::random_element());
        let point = a + b;
        points.push(point);
        normalized.push(point);
    }
    Element::batch_normalize(&mut normalized);
    for i in 0..num_points {
        let zz = points[i].z.sqr();
        let zzz = points[i].z * zz;
        assert_eq!(
            normalized[i].x * zz,
            points[i].x,
            "secp256r1 batch_normalize x mismatch at {}",
            i
        );
        assert_eq!(
            normalized[i].y * zzz,
            points[i].y,
            "secp256r1 batch_normalize y mismatch at {}",
            i
        );
    }
}

// =========================================================================
// Other ECC tests — WnafTwoBitWindow, InfinityBatchMulByScalarIsInfinity
// =========================================================================

#[test]
fn wnaf_two_bit_window() {
    use crate::groups::wnaf;
    let r = Fr::random_element();
    let raw = r.from_montgomery_form();
    let buffer: [u64; 2] = [raw.data[0], raw.data[1] & 0x7fffffffffffffff];
    let wnaf_bits = 2usize;
    let wnaf_entries = (wnaf::SCALAR_BITS + wnaf_bits - 1) / wnaf_bits; // 64
    let mut table = vec![0u64; wnaf_entries + 1];
    let mut skew = false;

    wnaf::fixed_wnaf(&buffer, &mut table, &mut skew, 0, 1, wnaf_bits);

    // All entries should have abs value < 2^(wnaf_bits-1) = 2
    for i in 0..wnaf_entries {
        let abs_val = table[i] & 0x0fffffff;
        assert!(abs_val < 2, "WNAF 2-bit entry {} abs value {} >= 2", i, abs_val);
    }

    // Recover scalar and verify
    let recovered = recover_wnaf_scalar(&table, skew, wnaf_entries, wnaf_bits);
    assert_eq!(recovered[0], buffer[0], "WNAF 2-bit: lo mismatch");
    assert_eq!(recovered[1], buffer[1], "WNAF 2-bit: hi mismatch");
}

#[test]
fn infinity_batch_mul_by_scalar_is_infinity() {
    // Multiplying infinity points by scalars should return infinity
    let inf = Element::<Bn254G1Params>::infinity();
    let scalar1 = Fr::random_element();
    let scalar2 = Fr::from(42u64);
    let result1 = inf.mul(&scalar1);
    let result2 = inf.mul(&scalar2);
    assert!(result1.is_point_at_infinity(), "inf * random scalar should be infinity");
    assert!(result2.is_point_at_infinity(), "inf * 42 should be infinity");
}
