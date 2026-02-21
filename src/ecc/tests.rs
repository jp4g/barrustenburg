use crate::ecc::curves::bn254::{Bn254FqParams, Bn254FrParams, Bn254G1Params, Fq, Fr};
use crate::ecc::curves::grumpkin::GrumpkinG1Params;
use crate::ecc::curves::secp256k1::{
    Secp256k1FqParams, Secp256k1FrParams, Secp256k1G1Params, Secp256k1Fr,
};
use crate::ecc::curves::secp256r1::{Secp256r1FqParams, Secp256r1G1Params};
use crate::ecc::fields::field::Field;
use crate::ecc::fields::field_params::FieldParams;
use crate::ecc::groups::affine_element::AffineElement;
use crate::ecc::groups::element::Element;

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
    use crate::ecc::curves::secp256r1::Secp256r1Fq;
    use crate::ecc::groups::curve_params::CurveParams;
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
    use crate::ecc::curves::secp256r1::Secp256r1G1Params;
    use crate::ecc::groups::curve_params::CurveParams;

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
    use crate::ecc::curves::secp256r1::Secp256r1FqParams;
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
    use crate::ecc::curves::secp256r1::Secp256r1G1Params;
    use crate::ecc::groups::curve_params::CurveParams;

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
    use crate::ecc::curves::secp256r1::{Secp256r1FqParams, Secp256r1G1Params};
    use crate::ecc::groups::curve_params::CurveParams;
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
    use crate::ecc::curves::secp256r1::Secp256r1FqParams;
    use crate::ecc::groups::curve_params::CurveParams;
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
    use crate::ecc::curves::secp256r1::{Secp256r1FqParams, Secp256r1G1Params};
    use crate::ecc::groups::curve_params::CurveParams;
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
    use crate::ecc::curves::secp256r1::{Secp256r1FqParams, Secp256r1G1Params};
    use crate::ecc::groups::curve_params::CurveParams;
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
// =========================================================================

#[test]
fn bn254_fr_split_endomorphism_scalars() {
    // Split a known scalar and verify k1 + k2 * lambda == k (mod r)
    let k = Fr::from(12345u64);
    let k_conv = k.from_montgomery_form();
    let (k1_limbs, k2_limbs) = k_conv.split_into_endomorphism_scalars();

    // Reconstruct: k1 + k2 * lambda mod r
    let k1 = Fr::from_limbs([k1_limbs[0], k1_limbs[1], 0, 0]);
    let k2 = Fr::from_limbs([k2_limbs[0], k2_limbs[1], 0, 0]);
    let lambda = Fr::from_raw(Bn254FrParams::CUBE_ROOT);
    let reconstructed = k1 + k2 * lambda;
    assert_eq!(reconstructed, k, "k1 + k2*lambda should equal k (mod r)");
}

#[test]
fn bn254_fr_split_endomorphism_scalars_one() {
    let k = Fr::one();
    let k_conv = k.from_montgomery_form();
    let (k1_limbs, k2_limbs) = k_conv.split_into_endomorphism_scalars();

    let k1 = Fr::from_limbs([k1_limbs[0], k1_limbs[1], 0, 0]);
    let k2 = Fr::from_limbs([k2_limbs[0], k2_limbs[1], 0, 0]);
    let lambda = Fr::from_raw(Bn254FrParams::CUBE_ROOT);
    let reconstructed = k1 + k2 * lambda;
    assert_eq!(reconstructed, k, "split of 1 should reconstruct to 1");
}

#[test]
fn secp256k1_fr_split_endomorphism_scalars() {
    let k = Secp256k1Fr::from(999999u64);
    let k_conv = k.from_montgomery_form();
    let (k1_limbs, k2_limbs) = k_conv.split_into_endomorphism_scalars();

    let k1 = Secp256k1Fr::from_limbs([k1_limbs[0], k1_limbs[1], 0, 0]);
    let k2 = Secp256k1Fr::from_limbs([k2_limbs[0], k2_limbs[1], 0, 0]);
    let lambda = Secp256k1Fr::from_raw(Secp256k1FrParams::CUBE_ROOT);
    let reconstructed = k1 + k2 * lambda;
    assert_eq!(
        reconstructed, k,
        "secp256k1: k1 + k2*lambda should equal k (mod r)"
    );
}

// =========================================================================
// WNAF tests
// =========================================================================

#[test]
fn wnaf_encode_basic_properties() {
    use crate::ecc::groups::wnaf;

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
    use crate::ecc::groups::wnaf;

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
    use crate::ecc::groups::wnaf;

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

// =========================================================================
// mul_with_endomorphism tests
// =========================================================================

#[test]
fn bn254_mul_with_endomorphism_matches_basic() {
    let g = Element::<Bn254G1Params>::one();
    // Test with several scalar values
    for val in [1u64, 2, 7, 42, 12345, 999999] {
        let scalar = Fr::from(val);
        let expected = g.mul_without_endomorphism(&scalar);
        let result = g.mul_with_endomorphism(&scalar);
        assert_eq!(
            result, expected,
            "endomorphism mul should match basic mul for scalar={}",
            val
        );
    }
}

#[test]
fn bn254_mul_with_endomorphism_large_scalar() {
    let g = Element::<Bn254G1Params>::one();
    let scalar = Fr::from_limbs([
        0xDEADBEEF12345678,
        0xCAFEBABE87654321,
        0x1234567890ABCDEF,
        0x0FEDCBA098765432,
    ]);
    let expected = g.mul_without_endomorphism(&scalar);
    let result = g.mul_with_endomorphism(&scalar);
    assert_eq!(result, expected, "endomorphism mul should match for large scalar");
}

#[test]
fn bn254_mul_dispatch() {
    // Test that mul() dispatches to endomorphism path and gives correct result
    let g = Element::<Bn254G1Params>::one();
    let scalar = Fr::from(42u64);
    let expected = g.mul_without_endomorphism(&scalar);
    let result = g.mul(&scalar);
    assert_eq!(result, expected, "mul() should match mul_without_endomorphism");
}

#[test]
fn secp256k1_mul_with_endomorphism_matches_basic() {
    let g = Element::<Secp256k1G1Params>::one();
    for val in [1u64, 2, 7, 42, 12345] {
        let scalar = Secp256k1Fr::from(val);
        let expected = g.mul_without_endomorphism(&scalar);
        let result = g.mul_with_endomorphism(&scalar);
        assert_eq!(
            result, expected,
            "secp256k1: endomorphism mul should match basic mul for scalar={}",
            val
        );
    }
}

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

#[test]
fn batch_normalize_matches_individual() {
    let g = Element::<Bn254G1Params>::one();
    let mut points: Vec<Element<Bn254G1Params>> = Vec::new();

    // Create several points with different z-coordinates
    for i in 1..=8u64 {
        points.push(g.scalar_mul(&[i, 0, 0, 0]));
    }

    // Save expected affine coordinates
    let expected_affines: Vec<_> = points.iter().map(|p| p.to_affine()).collect();

    // Batch normalize
    Element::batch_normalize(&mut points);

    // Verify each normalized point matches
    for (i, point) in points.iter().enumerate() {
        let affine = point.to_affine();
        assert_eq!(
            affine, expected_affines[i],
            "batch_normalize mismatch at index {}",
            i
        );
        // z should be 1 after normalization
        assert_eq!(point.z, Fq::one(), "z should be 1 after batch_normalize");
    }
}

#[test]
fn batch_normalize_with_infinity() {
    let g = Element::<Bn254G1Params>::one();
    let mut points = vec![
        g.scalar_mul(&[1, 0, 0, 0]),
        Element::infinity(),
        g.scalar_mul(&[3, 0, 0, 0]),
    ];

    let expected_0 = points[0].to_affine();
    let expected_2 = points[2].to_affine();

    Element::batch_normalize(&mut points);

    assert_eq!(points[0].to_affine(), expected_0);
    assert!(
        points[1].is_point_at_infinity() || points[1].z == Fq::one(),
        "infinity point should be handled"
    );
    assert_eq!(points[2].to_affine(), expected_2);
}

#[test]
fn batch_normalize_empty() {
    let mut points: Vec<Element<Bn254G1Params>> = Vec::new();
    Element::batch_normalize(&mut points); // Should not panic
}

#[test]
fn batch_normalize_single() {
    let g = Element::<Bn254G1Params>::one();
    let mut points = vec![g.scalar_mul(&[5, 0, 0, 0])];
    let expected = points[0].to_affine();
    Element::batch_normalize(&mut points);
    assert_eq!(points[0].to_affine(), expected);
}
