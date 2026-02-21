// BN254 optimal ate pairing implementation.
//
// C++ source: ecc/curves/bn254/pairing.hpp, pairing_impl.hpp
//
// Implements the reduced ate pairing e: G1 x G2 -> GT (Fq12).

use crate::curves::bn254::{
    Fq, Fq2, Fq12, G1Element, G2AffineElement, G2Element,
};
use crate::fields::field6::Field6;
use crate::fields::field12::{EllCoeffs, Field12};

use crate::curves::bn254::Bn254FqParams;

const LOOP_LENGTH: usize = 64;
const NEG_Z_LOOP_LENGTH: usize = 62;
const PRECOMPUTED_COEFFICIENTS_LENGTH: usize = 87;

/// Miller loop iteration bits for the ate pairing.
/// 0 = no addition, 1 = add Q, 3 = add -Q.
const LOOP_BITS: [u8; LOOP_LENGTH] = [
    1, 0, 1, 0, 0, 0, 3, 0, 3, 0, 0, 0, 3, 0, 1, 0,
    3, 0, 0, 3, 0, 0, 0, 0, 0, 1, 0, 0, 3, 0, 1, 0,
    0, 3, 0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 3, 0, 3, 0,
    0, 1, 0, 0, 0, 3, 0, 0, 3, 0, 1, 0, 1, 0, 0, 0,
];

/// Bits for the final exponentiation exp-by-neg-z step.
const NEG_Z_LOOP_BITS: [bool; NEG_Z_LOOP_LENGTH] = [
    false, false, false, true,  false, false, true,  true,
    true,  false, true,  false, false, true,  true,  false,
    false, true,  false, false, true,  false, true,  false,
    true,  true,  false, true,  false, false, false, true,
    false, false, true,  false, true,  false, false, true,
    true,  false, true,  false, false, true,  false, false,
    false, false, true,  false, false, true,  true,  true,
    true,  true,  false, false, false, true,
];

type EllCoeffsBn254 = EllCoeffs<Bn254FqParams>;

/// Precomputed Miller line evaluations for a fixed G2 point.
pub struct MillerLines {
    pub lines: [EllCoeffsBn254; PRECOMPUTED_COEFFICIENTS_LENGTH],
}

impl MillerLines {
    fn new() -> Self {
        Self {
            lines: [EllCoeffsBn254 {
                o: Fq2::zero(),
                vw: Fq2::zero(),
                vv: Fq2::zero(),
            }; PRECOMPUTED_COEFFICIENTS_LENGTH],
        }
    }
}

/// Frobenius endomorphism on G2: multiply-by-q map.
fn mul_by_q(a: &G2Element, two_inv: Fq) -> G2Element {
    let _ = two_inv; // unused here but kept for API symmetry
    let t0 = a.x.frobenius_map();
    let t1 = a.y.frobenius_map();
    G2Element::new(
        Fq2::twist_mul_by_q_x() * t0,
        Fq2::twist_mul_by_q_y() * t1,
        a.z.frobenius_map(),
    )
}

/// Doubling step for the Miller loop (flipped version).
fn doubling_step(current: &mut G2Element, ell: &mut EllCoeffsBn254, two_inv: Fq) {
    let mut a = current.x.mul_by_fq(two_inv);
    a *= current.y;

    let b = current.y.sqr();
    let c = current.z.sqr();
    let mut d = c + c;
    d = d + c;
    let e = d * Fq2::twist_coeff_b();
    let mut f = e + e;
    f = f + e;

    let g = (b + f).mul_by_fq(two_inv);
    let mut h = current.y + current.z;
    h = h.sqr();
    let i = b + c;
    h = h - i;
    let i_val = e - b;
    let j = current.x.sqr();
    let ee = e.sqr();
    let k = b - f;
    current.x = a * k;

    let mut k = ee + ee;
    k = k + ee;

    let c_val = g.sqr();
    current.y = c_val - k;
    current.z = b * h;

    ell.o = Field6::<Bn254FqParams>::mul_by_non_residue(&i_val);
    ell.vw = -h;
    let mut vv = j + j;
    vv = vv + j;
    ell.vv = vv;
}

/// Mixed addition step for the Miller loop (flipped version).
fn mixed_addition_step(base: &G2Element, q: &mut G2Element, line: &mut EllCoeffsBn254) {
    let mut d = base.x * q.z;
    d = q.x - d;

    let mut e = base.y * q.z;
    e = q.y - e;

    let f = d.sqr();
    let g = e.sqr();
    let h = d * f;
    let i = q.x * f;

    let mut j = q.z * g;
    j = j + h;
    j = j - i;
    j = j - i;

    q.x = d * j;
    let i = i - j;
    let i = i * e;
    let j = q.y * h;
    q.y = i - j;
    q.z = q.z * h;

    let h = e * base.x;
    let i = d * base.y;

    let h = h - i;
    line.o = Field6::<Bn254FqParams>::mul_by_non_residue(&h);
    line.vv = -e;
    line.vw = d;
}

/// Precompute Miller lines for a fixed G2 point.
fn precompute_miller_lines(q: &G2Element, lines: &mut MillerLines, two_inv: Fq) {
    assert!(!q.is_point_at_infinity(), "Cannot compute Miller lines for point at infinity");

    let q_neg = G2Element::new(q.x, -q.y, Fq2::one());
    let mut work_point = *q;

    let mut it = 0;
    for &loop_bit in &LOOP_BITS {
        doubling_step(&mut work_point, &mut lines.lines[it], two_inv);
        it += 1;
        if loop_bit == 1 {
            let q_copy = G2Element::new(q.x, q.y, q.z);
            mixed_addition_step(&q_copy, &mut work_point, &mut lines.lines[it]);
            it += 1;
        } else if loop_bit == 3 {
            let q_neg_copy = G2Element::new(q_neg.x, q_neg.y, q_neg.z);
            mixed_addition_step(&q_neg_copy, &mut work_point, &mut lines.lines[it]);
            it += 1;
        }
    }

    let q1 = mul_by_q(q, two_inv);
    let q2 = -mul_by_q(&q1, two_inv);
    mixed_addition_step(&q1, &mut work_point, &mut lines.lines[it]);
    it += 1;
    mixed_addition_step(&q2, &mut work_point, &mut lines.lines[it]);
}

/// Miller loop: evaluate the pairing product using precomputed lines.
fn miller_loop(p: &G1Element, lines: &MillerLines) -> Fq12 {
    let mut work_scalar = Fq12::one();
    let mut it = 0;

    for &loop_bit in &LOOP_BITS {
        work_scalar = work_scalar.sqr();

        let mut work_line = EllCoeffsBn254 {
            o: lines.lines[it].o,
            vw: lines.lines[it].vw.mul_by_fq(p.y),
            vv: lines.lines[it].vv.mul_by_fq(p.x),
        };
        work_scalar.self_sparse_mul(&work_line);
        it += 1;

        if loop_bit != 0 {
            work_line = EllCoeffsBn254 {
                o: lines.lines[it].o,
                vw: lines.lines[it].vw.mul_by_fq(p.y),
                vv: lines.lines[it].vv.mul_by_fq(p.x),
            };
            work_scalar.self_sparse_mul(&work_line);
            it += 1;
        }
    }

    let mut work_line = EllCoeffsBn254 {
        o: lines.lines[it].o,
        vw: lines.lines[it].vw.mul_by_fq(p.y),
        vv: lines.lines[it].vv.mul_by_fq(p.x),
    };
    work_scalar.self_sparse_mul(&work_line);
    it += 1;

    work_line = EllCoeffsBn254 {
        o: lines.lines[it].o,
        vw: lines.lines[it].vw.mul_by_fq(p.y),
        vv: lines.lines[it].vv.mul_by_fq(p.x),
    };
    work_scalar.self_sparse_mul(&work_line);

    work_scalar
}

/// Easy part of the final exponentiation: f^((p^6 - 1)(p^2 + 1)).
fn final_exponentiation_easy_part(elt: &Fq12) -> Fq12 {
    let a = Field12::new(elt.c0, -elt.c1);
    let a = a * elt.invert();
    a * a.frobenius_map_two()
}

/// Exponentiation by -z (the BN254 parameter) in the cyclotomic subgroup.
fn final_exponentiation_exp_by_neg_z(elt: &Fq12) -> Fq12 {
    let mut r = *elt;
    for &bit in &NEG_Z_LOOP_BITS {
        r = r.cyclotomic_squared();
        if bit {
            r = r * *elt;
        }
    }
    r.unitary_inverse()
}

/// Hard part of the final exponentiation.
fn final_exponentiation_tricky_part(elt: &Fq12) -> Fq12 {
    let a = final_exponentiation_exp_by_neg_z(elt);
    let b = a.cyclotomic_squared();
    let c = b.cyclotomic_squared();
    let d = b * c;
    let e = final_exponentiation_exp_by_neg_z(&d);
    let f = e.cyclotomic_squared();
    let g = final_exponentiation_exp_by_neg_z(&f);
    let h = d.unitary_inverse();
    let ii = g.unitary_inverse();
    let j = ii * e;
    let k = j * h;
    let l = b * k;
    let m = e * k;
    let n = m * *elt;
    let o = l.frobenius_map_one();
    let p = o * n;
    let q = k.frobenius_map_two();
    let r = p * q;
    let s = elt.unitary_inverse();
    let t = l * s;
    let u = t.frobenius_map_three();
    r * u
}

/// Compute the reduced ate pairing e(P, Q) for BN254.
///
/// Returns an element of GT = Fq12.
/// Returns Fq12::one() if either P or Q is the point at infinity.
pub fn reduced_ate_pairing(p_affine: &crate::curves::bn254::G1Affine, q_affine: &G2AffineElement) -> Fq12 {
    use crate::groups::element::Element;

    let p = Element::from_affine(p_affine);
    let q = G2Element::from_affine(q_affine);

    if q.is_point_at_infinity() || p.is_point_at_infinity() {
        return Fq12::one();
    }

    let two_inv = Fq::from(2u64).invert();

    let mut lines = MillerLines::new();
    precompute_miller_lines(&q, &mut lines, two_inv);

    let result = miller_loop(&p, &lines);
    let result = final_exponentiation_easy_part(&result);
    final_exponentiation_tricky_part(&result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curves::bn254::{G1Affine, G2AffineElement, Fq12, Fr};
    use crate::groups::element::Element;

    #[test]
    fn test_pairing_nondegeneracy() {
        // e(G1, G2) should not be the identity
        let p = G1Affine::one();
        let q = G2AffineElement::generator();
        let result = reduced_ate_pairing(&p, &q);
        assert_ne!(result, Fq12::one(), "e(G1, G2) must not be 1");
    }

    #[test]
    fn test_pairing_bilinearity() {
        // e(a*P, Q) == e(P, Q)^a
        // Tested via: e(a*P, Q) == e(P, a*Q) (equivalent by bilinearity)
        //
        // Actually, since we don't have G2 scalar mul, we test:
        // e(a*P, Q) * e(b*P, Q) == e((a+b)*P, Q)
        let p_aff = G1Affine::one();
        let q_aff = G2AffineElement::generator();

        // Use small scalars for speed
        let a = Fr::from(7u64);
        let b = Fr::from(13u64);
        let ab = a + b; // 20

        let p_proj = Element::from_affine(&p_aff);
        let ap = p_proj.mul_without_endomorphism(&a).to_affine();
        let bp = p_proj.mul_without_endomorphism(&b).to_affine();
        let abp = p_proj.mul_without_endomorphism(&ab).to_affine();

        let e_ap_q = reduced_ate_pairing(&ap, &q_aff);
        let e_bp_q = reduced_ate_pairing(&bp, &q_aff);
        let e_abp_q = reduced_ate_pairing(&abp, &q_aff);

        // e(a*P, Q) * e(b*P, Q) should equal e((a+b)*P, Q)
        let product = e_ap_q * e_bp_q;
        assert_eq!(product, e_abp_q, "Pairing bilinearity failed");
    }

    #[test]
    fn test_pairing_infinity() {
        let p_inf = G1Affine::infinity();
        let q = G2AffineElement::generator();
        let result = reduced_ate_pairing(&p_inf, &q);
        assert_eq!(result, Fq12::one(), "e(O, Q) must be 1");

        let p = G1Affine::one();
        let q_inf = G2AffineElement::infinity();
        let result = reduced_ate_pairing(&p, &q_inf);
        assert_eq!(result, Fq12::one(), "e(P, O) must be 1");
    }
}
