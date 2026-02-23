//! Port of `barretenberg/stdlib/translator_vm_verifier/` — recursive Translator VM verifier.
//!
//! Implements the in-circuit verification of Translator VM proofs. The Translator
//! proves correct translation between BN254 and Grumpkin curve operations, bridging
//! the ECCVM (which operates on Grumpkin) with the main circuit (which operates on
//! BN254).
//!
//! ## Architecture
//!
//! The recursive verifier follows the same phases as the native verifier:
//! 1. **Fiat-Shamir**: hash VK, receive commitments, compute challenges
//! 2. **Translation data**: decompose non-native field elements into limbs for
//!    relation parameters (evaluation_input_x, batching_challenge_v, accumulated_result)
//! 3. **Sumcheck phase**: verify sumcheck equations with ZK (Libra) over translator
//!    relations (decomposition, delta range, non-native field, permutation, extra)
//! 4. **PCS phase**: Shplemini batch opening + KZG pairing point reduction
//!
//! All operations use circuit types (`FieldT`, `CommitmentT`, `BigFieldT`) so
//! verification logic is encoded as constraints in an outer circuit.
//!
//! ## Current Status
//!
//! **Blocked on B27 (translator_vm native crate).**
//!
//! The translator VM native prover/verifier/flavor types have not yet been ported.
//! This module provides:
//! - Translation data limb decomposition (both native and recursive variants)
//! - Translator-specific constants matching C++ `TranslatorFlavor`
//! - Verifier type skeleton and output types
//! - Tests validating the limb decomposition logic
//!
//! When B27 lands, this module will be completed with full verification logic
//! following the pattern established by `honk_verifier.rs`.
//!
//! ## Differences from Ultra Honk Recursive Verifier
//!
//! - **ZK Sumcheck**: Translator uses Libra masking (3 extra commitments)
//! - **Non-native field**: Translation data uses `BigFieldT` (BN254 Fq over Fr)
//!   with 4 binary limbs + 1 prime basis limb (68-bit limbs)
//! - **Op queue wires**: 4 commitments provided externally (from merge protocol)
//! - **Constant circuit size**: `CONST_TRANSLATOR_LOG_N` (fixed, not variable)
//! - **Concatenated polynomials**: Shifted concat evals reconstructed from groups
//!
//! C++ source: `barretenberg/stdlib/translator_vm_verifier/`
//! C++ verifier: `barretenberg/translator_vm/translator_verifier.{hpp,cpp}`

use bbrs_ecc::curves::bn254::{Bn254FqParams, Bn254FrParams, Fq, Fr};
use bbrs_numeric::U256;
use bbrs_relations::relation_parameters::{
    RelationParameters, NUM_BINARY_LIMBS_IN_GOBLIN_TRANSLATOR, NUM_CHALLENGE_POWERS_IN_GOBLIN_TRANSLATOR,
    NUM_TOTAL_LIMBS,
};

use crate::honk_verifier::{CommitmentT, PairingPointsAccumulator, StdlibTranscript};
use crate::primitives::bigfield::BigFieldT;
use crate::primitives::field::FieldT;
use crate::primitives::witness::BuilderRef;

// ════════════════════════════════════════════════════════════════════════════
//  Constants
// ════════════════════════════════════════════════════════════════════════════

/// Number of bits per binary limb in translator non-native field decomposition.
///
/// Matches C++ `TranslatorFlavor::NUM_LIMB_BITS = 68`.
pub const NUM_LIMB_BITS: u32 = 68;

/// Number of op queue wire polynomials.
///
/// Matches C++ `TranslatorFlavor::NUM_OP_QUEUE_WIRES = 4`.
pub const NUM_OP_QUEUE_WIRES: usize = 4;

/// Number of Libra commitments used in ZK sumcheck.
///
/// Matches C++ `NUM_LIBRA_COMMITMENTS = 3` (concatenation, grand sum, quotient).
pub const NUM_LIBRA_COMMITMENTS: usize = 3;

/// Circuit field element type (BN254 Fr as circuit witness).
type FF = FieldT<Bn254FrParams>;

/// Non-native base field type (BN254 Fq represented over Fr).
type BF = BigFieldT<Bn254FrParams, Bn254FqParams>;

/// Builder context reference type.
type BuilderCtx = BuilderRef<Bn254FrParams>;

// ════════════════════════════════════════════════════════════════════════════
//  Translation data — native limb decomposition
// ════════════════════════════════════════════════════════════════════════════

/// Slice a U256 value into 4 binary limbs of `NUM_LIMB_BITS` each.
///
/// Port of the native C++ `compute_four_limbs` lambda in `translator_verifier.cpp`.
///
/// Returns `[limb_0, limb_1, limb_2, limb_3]` where each limb is a native Fr
/// element containing bits `[i*68, (i+1)*68)` of the input.
fn compute_four_limbs_native(value: &U256) -> [Fr; NUM_BINARY_LIMBS_IN_GOBLIN_TRANSLATOR] {
    let mut result = [Fr::zero(); NUM_BINARY_LIMBS_IN_GOBLIN_TRANSLATOR];
    for i in 0..NUM_BINARY_LIMBS_IN_GOBLIN_TRANSLATOR {
        let start = (i as u32) * NUM_LIMB_BITS;
        let end = start + NUM_LIMB_BITS;
        let sliced = slice_u256(value, start, end);
        result[i] = Fr::from_limbs(*sliced.as_words());
    }
    result
}

/// Slice a U256 value into 5 limbs: 4 binary + 1 prime basis (full value mod Fr).
///
/// Port of the native C++ `compute_five_limbs` lambda in `translator_verifier.cpp`.
///
/// The 5th limb is the full value reduced to the native field, used for
/// `prime_basis_limb` in `TranslatorNonNativeFieldRelation`.
fn compute_five_limbs_native(value: &U256) -> [Fr; NUM_TOTAL_LIMBS] {
    let mut result = [Fr::zero(); NUM_TOTAL_LIMBS];
    let four = compute_four_limbs_native(value);
    result[..4].copy_from_slice(&four);
    // 5th limb: full value as native field element
    result[4] = Fr::from_limbs(*value.as_words());
    result
}

/// Populate relation parameters with translation data from ECCVM verifier (native variant).
///
/// Port of the native C++ `put_translation_data_in_relation_parameters_impl` with
/// `requires(!Flavor::Curve::is_stdlib_type)`.
///
/// Decomposes `evaluation_input_x`, `batching_challenge_v`, and `accumulated_result`
/// from their Fq (non-native) representation into limbed form suitable for
/// translator relations.
///
/// # Arguments
///
/// * `relation_parameters` - Output: populated with limbed translation data
/// * `evaluation_input_x` - Polynomial evaluation challenge (from ECCVM)
/// * `batching_challenge_v` - Batching challenge for translation polynomials (from ECCVM)
/// * `accumulated_result` - Accumulated translation result (from ECCVM)
pub fn put_translation_data_native(
    relation_parameters: &mut RelationParameters<Fr>,
    evaluation_input_x: &Fq,
    batching_challenge_v: &Fq,
    accumulated_result: &Fq,
) {
    // Convert Fq values to U256 for limb slicing
    let x_u256 = fq_to_u256(evaluation_input_x);
    let result_u256 = fq_to_u256(accumulated_result);

    // evaluation_input_x: 5 limbs
    relation_parameters.evaluation_input_x = compute_five_limbs_native(&x_u256);

    // batching_challenge_v: 4 powers × 5 limbs each
    let mut v_power = *batching_challenge_v;
    for i in 0..NUM_CHALLENGE_POWERS_IN_GOBLIN_TRANSLATOR {
        let v_u256 = fq_to_u256(&v_power);
        relation_parameters.batching_challenge_v[i] = compute_five_limbs_native(&v_u256);
        v_power = v_power * *batching_challenge_v;
    }

    // accumulated_result: 4 binary limbs only (sufficient for equality checking)
    relation_parameters.accumulated_result = compute_four_limbs_native(&result_u256);
}

/// Convert an Fq field element to U256 representation.
fn fq_to_u256(fq: &Fq) -> U256 {
    let mont_out = fq.from_montgomery_form();
    U256::from_words(mont_out.data)
}

/// Slice bits `[start, end)` from a U256, returning a U256.
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

// ════════════════════════════════════════════════════════════════════════════
//  Translation data — recursive (in-circuit) limb extraction
// ════════════════════════════════════════════════════════════════════════════

/// Extract 4 binary limbs from a BigFieldT element for in-circuit use.
///
/// Port of the recursive C++ `compute_four_limbs` lambda with
/// `requires(Flavor::Curve::is_stdlib_type)`.
///
/// Uses `BigFieldT::binary_basis_limbs` which are already constrained
/// during bigfield construction.
fn compute_four_limbs_recursive(bf: &BF) -> [FF; NUM_BINARY_LIMBS_IN_GOBLIN_TRANSLATOR] {
    [
        bf.binary_basis_limbs[0].element.clone(),
        bf.binary_basis_limbs[1].element.clone(),
        bf.binary_basis_limbs[2].element.clone(),
        bf.binary_basis_limbs[3].element.clone(),
    ]
}

/// Extract 5 limbs (4 binary + 1 prime) from a BigFieldT for in-circuit use.
///
/// Port of the recursive C++ `compute_five_limbs` lambda with
/// `requires(Flavor::Curve::is_stdlib_type)`.
fn compute_five_limbs_recursive(bf: &BF) -> [FF; NUM_TOTAL_LIMBS] {
    [
        bf.binary_basis_limbs[0].element.clone(),
        bf.binary_basis_limbs[1].element.clone(),
        bf.binary_basis_limbs[2].element.clone(),
        bf.binary_basis_limbs[3].element.clone(),
        bf.prime_basis_limb.clone(),
    ]
}

/// Populate relation parameters with translation data (recursive/in-circuit variant).
///
/// Port of the recursive C++ `put_translation_data_in_relation_parameters_impl` with
/// `requires(Flavor::Curve::is_stdlib_type)`.
///
/// Extracts limbs from BigFieldT elements for use in translator relation constraints.
/// The limbs are already constrained by the bigfield construction, so this simply
/// re-packages them into the relation parameter arrays.
pub fn put_translation_data_recursive(
    relation_parameters: &mut RelationParameters<Fr>,
    evaluation_input_x: &BF,
    batching_challenge_v: &BF,
    accumulated_result: &BF,
) {
    // evaluation_input_x: 5 limbs from bigfield
    let x_limbs = compute_five_limbs_recursive(evaluation_input_x);
    for i in 0..NUM_TOTAL_LIMBS {
        relation_parameters.evaluation_input_x[i] = x_limbs[i].get_value();
    }

    // batching_challenge_v: 4 powers × 5 limbs
    let mut v_power = batching_challenge_v.clone();
    for i in 0..NUM_CHALLENGE_POWERS_IN_GOBLIN_TRANSLATOR {
        let v_limbs = compute_five_limbs_recursive(&v_power);
        for j in 0..NUM_TOTAL_LIMBS {
            relation_parameters.batching_challenge_v[i][j] = v_limbs[j].get_value();
        }
        v_power = &v_power * batching_challenge_v;
    }

    // accumulated_result: 4 binary limbs only
    let result_limbs = compute_four_limbs_recursive(accumulated_result);
    for i in 0..NUM_BINARY_LIMBS_IN_GOBLIN_TRANSLATOR {
        relation_parameters.accumulated_result[i] = result_limbs[i].get_value();
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  TranslatorRecursiveVerifierOutput
// ════════════════════════════════════════════════════════════════════════════

/// Output of the recursive Translator VM verifier.
///
/// Port of C++ `TranslatorVerifier_<TranslatorRecursiveFlavor>::ReductionResult`.
///
/// Contains pairing points for deferred KZG verification and status of internal
/// checks (sumcheck, consistency).
pub struct TranslatorRecursiveVerifierOutput {
    /// Accumulated pairing points for final pairing check.
    pub points_accumulator: PairingPointsAccumulator,
    /// Whether the reduction (sumcheck + consistency checks) succeeded.
    pub reduction_succeeded: bool,
}

// ════════════════════════════════════════════════════════════════════════════
//  TranslatorRecursiveVerifier
// ════════════════════════════════════════════════════════════════════════════

/// Recursive Translator VM verifier for in-circuit proof verification.
///
/// Port of C++ `TranslatorVerifier_<TranslatorRecursiveFlavor>`.
///
/// Constructs a circuit that verifies a Translator VM proof. The translator
/// proves correct translation between BN254 and Grumpkin curve operations
/// (bridging ECCVM and the main circuit).
///
/// ## Verification Inputs (from ECCVM verifier)
///
/// - `evaluation_input_x`: polynomial evaluation challenge (Fq element)
/// - `batching_challenge_v`: batching challenge for translation polynomials (Fq element)
/// - `accumulated_result`: accumulated ECCVM result (Fq element)
/// - `op_queue_commitments`: commitments to 4 op queue wires (from merge protocol)
///
/// ## Status
///
/// **Skeleton only.** Full implementation requires B27 (translator_vm native crate)
/// for the native flavor types, prover, and verification key definitions.
pub struct TranslatorRecursiveVerifier {
    /// Builder for the outer (verifier) circuit.
    pub builder: BuilderCtx,
    /// In-circuit transcript for Fiat-Shamir.
    pub transcript: StdlibTranscript,
    /// Evaluation challenge x (non-native Fq element as BigFieldT).
    pub evaluation_input_x: BF,
    /// Batching challenge v (non-native Fq element as BigFieldT).
    pub batching_challenge_v: BF,
    /// Accumulated result from ECCVM (non-native Fq element as BigFieldT).
    pub accumulated_result: BF,
    /// Op queue wire commitments (from merge protocol).
    pub op_queue_commitments: [CommitmentT; NUM_OP_QUEUE_WIRES],
}

impl TranslatorRecursiveVerifier {
    /// Create a new translator recursive verifier.
    ///
    /// # Arguments
    ///
    /// * `builder` - Builder for the outer circuit
    /// * `evaluation_input_x` - Evaluation challenge from ECCVM (BigFieldT over Fq)
    /// * `batching_challenge_v` - Batching challenge from ECCVM (BigFieldT over Fq)
    /// * `accumulated_result` - Accumulated translation result from ECCVM
    /// * `op_queue_commitments` - 4 op queue wire commitments from merge protocol
    pub fn new(
        builder: BuilderCtx,
        evaluation_input_x: BF,
        batching_challenge_v: BF,
        accumulated_result: BF,
        op_queue_commitments: [CommitmentT; NUM_OP_QUEUE_WIRES],
    ) -> Self {
        let transcript = StdlibTranscript::new(builder.clone());
        Self {
            builder,
            transcript,
            evaluation_input_x,
            batching_challenge_v,
            accumulated_result,
            op_queue_commitments,
        }
    }

    /// Reduce the translator proof to a pairing check (in-circuit).
    ///
    /// Port of C++ `TranslatorVerifier_<TranslatorRecursiveFlavor>::reduce_to_pairing_check`.
    ///
    /// ## Verification Flow
    ///
    /// 1. Load proof into transcript
    /// 2. Fiat-Shamir VK hash
    /// 3. Populate relation parameters with translation data (limb decomposition)
    /// 4. Receive Gemini masking poly commitment (ZK)
    /// 5. Set op queue wire commitments
    /// 6. Receive non-op-queue wire and ordered range constraint commitments
    /// 7. Get beta, gamma challenges → receive z_perm commitment
    /// 8. Get alpha challenge
    /// 9. Run ZK sumcheck with Libra masking
    /// 10. Shplemini batch opening → KZG pairing point reduction
    ///
    /// ## Current Status
    ///
    /// **Returns placeholder output.** Full implementation blocked on B27
    /// (translator_vm native crate) which provides:
    /// - `TranslatorFlavor` (entity definitions, commitment labels, relation types)
    /// - `TranslatorProvingKey` / `TranslatorProver` (for test proof generation)
    /// - Translator-specific relations (5 relation types)
    /// - `VerifierCommitments` and `CommitmentLabels` types
    pub fn reduce_to_pairing_check(self, _proof: &[Fr]) -> TranslatorRecursiveVerifierOutput {
        // TODO(B27): Implement full translator recursive verification.
        //
        // The verification flow will mirror honk_verifier.rs but with:
        // - TranslatorFlavor commitments/labels instead of UltraFlavor
        // - ZK sumcheck with Libra (3 masking commitments)
        // - Translation data in relation parameters (already implemented above)
        // - Concatenated polynomial reconstruction
        // - Constant translator circuit size (CONST_TRANSLATOR_LOG_N)
        //
        // See C++ translator_verifier.cpp reduce_to_pairing_check() for reference.

        TranslatorRecursiveVerifierOutput {
            points_accumulator: PairingPointsAccumulator::new(),
            reduction_succeeded: false,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  Tests
// ════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    use std::cell::RefCell;
    use std::rc::Rc;

    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use bbrs_ecc::curves::bn254::{Bn254FrParams, Fq};

    type P = Bn254FrParams;

    /// Create a builder ref (Rc<RefCell<UltraCircuitBuilder>>).
    fn make_builder_ref() -> BuilderCtx {
        Rc::new(RefCell::new(UltraCircuitBuilder::<P>::new()))
    }

    // ── Test 1: Native translation data limb decomposition ───────────

    /// Verify that native limb decomposition correctly splits Fq elements
    /// into 68-bit limbs and populates relation parameters.
    ///
    /// This tests `put_translation_data_native` which is the native-side
    /// counterpart of the recursive verifier's translation data handling.
    /// Correct limb decomposition is critical for the translator relations
    /// (decomposition, non-native field, accumulator transfer).
    #[test]
    fn test_native_translation_data_limb_decomposition() {
        // Create random Fq elements (simulating ECCVM outputs)
        let evaluation_x = Fq::random_element();
        let batching_v = Fq::random_element();
        let accumulated = Fq::random_element();

        let mut relation_params = RelationParameters::<Fr>::default();
        put_translation_data_native(
            &mut relation_params,
            &evaluation_x,
            &batching_v,
            &accumulated,
        );

        // Verify evaluation_input_x has 5 limbs
        // First 4 should be non-zero for a random element
        let x_u256 = fq_to_u256(&evaluation_x);
        for i in 0..4 {
            let expected_slice = slice_u256(&x_u256, (i as u32) * NUM_LIMB_BITS, ((i as u32) + 1) * NUM_LIMB_BITS);
            let expected_fr = Fr::from_limbs(*expected_slice.as_words());
            assert_eq!(
                relation_params.evaluation_input_x[i],
                expected_fr,
                "evaluation_input_x limb {} mismatch",
                i
            );
        }

        // 5th limb should be the full value mod Fr
        let x_full = Fr::from_limbs(*x_u256.as_words());
        assert_eq!(
            relation_params.evaluation_input_x[4],
            x_full,
            "evaluation_input_x prime basis limb mismatch"
        );

        // Verify accumulated_result has exactly 4 limbs (no prime basis)
        let result_u256 = fq_to_u256(&accumulated);
        for i in 0..4 {
            let expected_slice = slice_u256(
                &result_u256,
                (i as u32) * NUM_LIMB_BITS,
                ((i as u32) + 1) * NUM_LIMB_BITS,
            );
            let expected_fr = Fr::from_limbs(*expected_slice.as_words());
            assert_eq!(
                relation_params.accumulated_result[i],
                expected_fr,
                "accumulated_result limb {} mismatch",
                i
            );
        }

        // Verify batching_challenge_v powers: v^1, v^2, v^3, v^4
        let mut v_power = batching_v;
        for power_idx in 0..NUM_CHALLENGE_POWERS_IN_GOBLIN_TRANSLATOR {
            let v_u256 = fq_to_u256(&v_power);
            for limb_idx in 0..4 {
                let expected_slice = slice_u256(
                    &v_u256,
                    (limb_idx as u32) * NUM_LIMB_BITS,
                    ((limb_idx as u32) + 1) * NUM_LIMB_BITS,
                );
                let expected_fr = Fr::from_limbs(*expected_slice.as_words());
                assert_eq!(
                    relation_params.batching_challenge_v[power_idx][limb_idx],
                    expected_fr,
                    "batching_challenge_v[{}][{}] mismatch",
                    power_idx,
                    limb_idx
                );
            }
            v_power = v_power * batching_v;
        }
    }

    // ── Test 2: Recursive verifier construction and type consistency ──

    /// Verify that the translator recursive verifier can be constructed
    /// and that BigFieldT elements correctly provide limbs for translation
    /// data population.
    ///
    /// This test validates the recursive (in-circuit) path of translation
    /// data handling, which uses `BigFieldT::binary_basis_limbs` and
    /// `prime_basis_limb` instead of U256 slicing.
    #[test]
    fn test_recursive_verifier_construction_and_bigfield_limbs() {
        let builder = make_builder_ref();

        // Create BigFieldT elements from random Fq values
        let eval_x_native = Fq::random_element();
        let batch_v_native = Fq::random_element();
        let result_native = Fq::random_element();

        let eval_x_u256 = fq_to_u256(&eval_x_native);
        let batch_v_u256 = fq_to_u256(&batch_v_native);
        let result_u256 = fq_to_u256(&result_native);

        let eval_x_bf = BF::from_witness(builder.clone(), eval_x_u256);
        let batch_v_bf = BF::from_witness(builder.clone(), batch_v_u256);
        let result_bf = BF::from_witness(builder.clone(), result_u256);

        // Verify BigFieldT limbs match native U256 slicing
        for i in 0..4 {
            let native_slice = slice_u256(
                &eval_x_u256,
                (i as u32) * NUM_LIMB_BITS,
                ((i as u32) + 1) * NUM_LIMB_BITS,
            );
            let native_fr = Fr::from_limbs(*native_slice.as_words());
            let bf_limb_val = eval_x_bf.binary_basis_limbs[i].element.get_value();
            assert_eq!(
                bf_limb_val, native_fr,
                "BigFieldT limb {} should match native slice for evaluation_input_x",
                i
            );
        }

        // Construct the recursive verifier
        let op_queue_commitments = [
            CommitmentT::identity(),
            CommitmentT::identity(),
            CommitmentT::identity(),
            CommitmentT::identity(),
        ];

        let verifier = TranslatorRecursiveVerifier::new(
            builder.clone(),
            eval_x_bf.clone(),
            batch_v_bf.clone(),
            result_bf.clone(),
            op_queue_commitments,
        );

        // Verify verifier holds the correct translation inputs
        for i in 0..4 {
            assert_eq!(
                verifier.evaluation_input_x.binary_basis_limbs[i].element.get_value(),
                eval_x_bf.binary_basis_limbs[i].element.get_value(),
                "Verifier evaluation_input_x limb {} should match input",
                i
            );
        }

        // Verify the recursive limb extraction produces consistent results
        let five_limbs = compute_five_limbs_recursive(&eval_x_bf);
        let four_limbs_native = compute_four_limbs_native(&eval_x_u256);
        for i in 0..4 {
            assert_eq!(
                five_limbs[i].get_value(),
                four_limbs_native[i],
                "Recursive 5-limb decomposition should match native for limb {}",
                i
            );
        }

        // Verify the outer circuit has created witnesses
        let num_vars = builder.borrow().base.get_num_variables();
        assert!(
            num_vars > 1,
            "BigFieldT construction should create witnesses in the outer circuit. Got {} variables",
            num_vars,
        );
    }
}
