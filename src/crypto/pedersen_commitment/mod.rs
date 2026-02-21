use crate::crypto::generators::GeneratorContext;
use crate::ecc::curves::grumpkin::{self, GrumpkinFq};

type GrumpkinAffine = grumpkin::G1Affine;
type GrumpkinElement = grumpkin::G1Element;

/// Pedersen vector commitment over Grumpkin.
///
/// commit(inputs) = Σ inputs[i] * generators[i]
///
/// Inputs are GrumpkinFq (= BN254 Fr) elements. The scalar multiplication
/// uses the raw integer representation (stripped from Montgomery form).
pub fn commit_native(inputs: &[GrumpkinFq], ctx: &GeneratorContext) -> GrumpkinAffine {
    let generators = ctx.get_generators(inputs.len());
    let mut result = GrumpkinElement::infinity();

    for (i, &input) in inputs.iter().enumerate() {
        // Convert from Montgomery to standard form to get raw integer limbs
        let raw = input.from_montgomery_form();
        let gen_proj = GrumpkinElement::from_affine(&generators[i]);
        let term = gen_proj.scalar_mul(&raw.data);
        result += term;
    }

    result.to_affine()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_commit_one_one() {
        let x = GrumpkinFq::one();
        let result = commit_native(&[x, x], &GeneratorContext::default());

        let expected = GrumpkinAffine::new(
            GrumpkinFq::from_limbs([
                0x6f12760266cdfe15,
                0xbf13523c19d7afe3,
                0x82205fb73ee43215,
                0x2f7a8f9a6c969266,
            ]),
            GrumpkinFq::from_limbs([
                0x3878bcc92b6057e6,
                0xec84b46daddf72f4,
                0x10e39b18c1d24b33,
                0x01916b316adbbf0e,
            ]),
        );

        assert_eq!(result, expected, "commit([1, 1]) mismatch");
    }

    #[test]
    fn test_commit_zero_one() {
        let x = GrumpkinFq::zero();
        let y = GrumpkinFq::one();
        let result = commit_native(&[x, y], &GeneratorContext::default());

        // commit([0, 1]) = 0 * g[0] + 1 * g[1] = g[1]
        let expected = GrumpkinAffine::new(
            GrumpkinFq::from_limbs([
                0x8f71df4591bde402,
                0x198e860f5f395026,
                0x525e5bbed6e43ba1,
                0x054aa86a73cb8a34,
            ]),
            GrumpkinFq::from_limbs([
                0xeb621a6287cac126,
                0xf87254afc7407c04,
                0xf6046f44d71ac6fa,
                0x209dcfbf2cfb57f9,
            ]),
        );

        assert_eq!(result, expected, "commit([0, 1]) mismatch");
    }
}
