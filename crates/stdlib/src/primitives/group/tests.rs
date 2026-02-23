//! Tests for the cycle_group, cycle_scalar, straus_lookup_table, and straus_scalar_slice modules.
//!
//! Port of the C++ cycle_group.test.cpp, cycle_scalar.test.cpp,
//! straus_lookup_table.test.cpp, and straus_scalar_slice.test.cpp.

#[cfg(test)]
mod tests {
    use std::cell::RefCell;
    use std::rc::Rc;

    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use bbrs_ecc::curves::bn254::Bn254FrParams;
    use bbrs_ecc::curves::grumpkin::GrumpkinG1Params;
    use bbrs_ecc::fields::field::Field;
    use bbrs_ecc::groups::affine_element::AffineElement;
    use bbrs_ecc::groups::element::Element;

    use crate::primitives::bool::BoolT;
    use crate::primitives::field::FieldT;
    use crate::primitives::group::cycle_group::CycleGroupT;
    use crate::primitives::group::cycle_scalar::CycleScalarT;
    use crate::primitives::group::straus_lookup_table::StrausLookupTable;
    use crate::primitives::group::straus_scalar_slice::StrausScalarSlices;
    use crate::primitives::witness::{BuilderRef, WitnessT};

    type C = GrumpkinG1Params;
    type Fr = Field<Bn254FrParams>;
    type ScalarField = Field<<GrumpkinG1Params as bbrs_ecc::groups::curve_params::CurveParams>::ScalarFieldParams>;
    type Grumpkin = Element<GrumpkinG1Params>;
    type GrumpkinAffine = AffineElement<GrumpkinG1Params>;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    fn check_circuit(builder: &BuilderRef<Bn254FrParams>) -> Result<(), String> {
        UltraCircuitChecker::check(&mut builder.borrow_mut())
    }

    /// Generate test points by multiplying generator by sequential scalars.
    fn generate_test_points(count: usize) -> Vec<GrumpkinAffine> {
        let grumpkin_gen = Grumpkin::one();
        let mut points = Vec::with_capacity(count);
        for i in 0..count {
            let scalar = ScalarField::from((i as u64 + 7) * 13 + 1);
            let point = grumpkin_gen.mul_without_endomorphism(&scalar);
            points.push(point.to_affine());
        }
        points
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleScalar tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_cycle_scalar_from_witness() {
        let builder = make_builder();
        let scalar_val = ScalarField::random_element();
        let scalar = CycleScalarT::<C>::from_witness(builder.clone(), scalar_val);

        assert_eq!(scalar.get_value(), scalar_val);
        assert!(!scalar.is_constant());

        // Verify lo/hi reconstruction
        let lo_val = scalar.lo().get_value().from_montgomery_form();
        let hi_val = scalar.hi().get_value().from_montgomery_form();
        let mut reconstructed = [0u64; 4];
        reconstructed[0] = lo_val.data[0];
        reconstructed[1] = lo_val.data[1];
        reconstructed[2] = hi_val.data[0];
        reconstructed[3] = hi_val.data[1];
        let reconstructed_scalar = ScalarField::from_limbs(reconstructed);
        assert_eq!(reconstructed_scalar, scalar_val);

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_cycle_scalar_from_native() {
        let scalar_val = ScalarField::random_element();
        let scalar = CycleScalarT::<C>::from_native(scalar_val);

        assert_eq!(scalar.get_value(), scalar_val);
        assert!(scalar.is_constant());
    }

    #[test]
    fn test_cycle_scalar_from_lo_hi() {
        let builder = make_builder();
        let scalar_val = ScalarField::random_element();
        let mont = scalar_val.from_montgomery_form();

        let lo = Field::<Bn254FrParams>::from_limbs([mont.data[0], mont.data[1], 0, 0]);
        let hi = Field::<Bn254FrParams>::from_limbs([mont.data[2], mont.data[3], 0, 0]);

        let lo_field = FieldT::from_witness(builder.clone(), lo);
        let hi_field = FieldT::from_witness(builder.clone(), hi);

        let scalar = CycleScalarT::<C>::from_lo_hi(lo_field, hi_field);
        assert_eq!(scalar.get_value(), scalar_val);

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_cycle_scalar_zero() {
        let builder = make_builder();
        let scalar_val = ScalarField::zero();
        let scalar = CycleScalarT::<C>::from_witness(builder.clone(), scalar_val);

        assert_eq!(scalar.get_value(), ScalarField::zero());
        assert_eq!(scalar.lo().get_value(), Fr::zero());
        assert_eq!(scalar.hi().get_value(), Fr::zero());

        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  StrausScalarSlice tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_straus_scalar_slice_read_and_reconstruction() {
        let builder = make_builder();
        let scalar_val = ScalarField::random_element();
        let scalar = CycleScalarT::<C>::from_witness(builder.clone(), scalar_val);

        let table_bits = 4;
        let slices = StrausScalarSlices::<C>::new(builder.clone(), &scalar, table_bits);

        // Read all slices and verify reconstruction
        let max_slice_val: u64 = (1u64 << table_bits) - 1;
        let mut reconstructed = [0u64; 4];
        for (i, native_val) in slices.slices_native.iter().enumerate() {
            assert!(*native_val <= max_slice_val, "slice {} = {} exceeds max {}", i, native_val, max_slice_val);
            // Reconstruct: each slice contributes `slice_val << (i * table_bits)`
            let bit_offset = i * table_bits;
            let limb_idx = bit_offset / 64;
            let bit_idx = bit_offset % 64;
            if limb_idx < 4 {
                reconstructed[limb_idx] |= native_val << bit_idx;
                if bit_idx + table_bits > 64 && limb_idx + 1 < 4 {
                    reconstructed[limb_idx + 1] |= native_val >> (64 - bit_idx);
                }
            }
        }
        let reconstructed_scalar = ScalarField::from_limbs(reconstructed);
        assert_eq!(reconstructed_scalar, scalar_val);

        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  StrausLookupTable tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_straus_lookup_table_construction() {
        let builder = make_builder();
        let base_point_native = Grumpkin::random_element();
        let offset_gen_native = Grumpkin::random_element();

        let base_point = CycleGroupT::<C>::from_witness(builder.clone(), &base_point_native.to_affine());
        let offset_gen = CycleGroupT::<C>::from_witness(builder.clone(), &offset_gen_native.to_affine());

        let table_bits = 4;
        let _table = StrausLookupTable::<C>::new(
            builder.clone(),
            &base_point,
            &offset_gen,
            table_bits,
            None,
        );

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_straus_lookup_table_read() {
        let builder = make_builder();
        let base_point_native = Grumpkin::random_element();
        let offset_gen_native = Grumpkin::random_element();

        let base_point = CycleGroupT::<C>::from_witness(builder.clone(), &base_point_native.to_affine());
        let offset_gen = CycleGroupT::<C>::from_witness(builder.clone(), &offset_gen_native.to_affine());

        let table_bits = 4;
        let table = StrausLookupTable::<C>::new(
            builder.clone(),
            &base_point,
            &offset_gen,
            table_bits,
            None,
        );

        let table_size = 1usize << table_bits;
        for i in 0..table_size {
            let index = FieldT::from_witness(builder.clone(), Fr::from(i as u64));
            let result = table.read(&index);

            // Expected value: offset_gen + i * base_point
            let mut expected = offset_gen_native.clone();
            for _ in 0..i {
                expected = expected + base_point_native.clone();
            }
            let expected_affine = expected.to_affine();
            assert_eq!(result.get_value(), expected_affine, "Table read at index {} failed", i);
        }

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_straus_lookup_table_with_hints() {
        let builder = make_builder();
        let base_point_native = Grumpkin::random_element();
        let offset_gen_native = Grumpkin::random_element();

        let base_point = CycleGroupT::<C>::from_witness(builder.clone(), &base_point_native.to_affine());
        let offset_gen = CycleGroupT::<C>::from_witness(builder.clone(), &offset_gen_native.to_affine());

        let table_bits = 3;

        // Compute hints
        let hints_elements =
            StrausLookupTable::<C>::compute_native_table(&base_point_native, &offset_gen_native, table_bits);
        let hints_affine: Vec<GrumpkinAffine> = hints_elements[1..].iter().map(|e| e.to_affine()).collect();

        let table = StrausLookupTable::<C>::new(
            builder.clone(),
            &base_point,
            &offset_gen,
            table_bits,
            Some(&hints_affine),
        );

        let index_val = 5;
        let index = FieldT::from_witness(builder.clone(), Fr::from(index_val as u64));
        let result = table.read(&index);

        let mut expected = offset_gen_native.clone();
        for _ in 0..index_val {
            expected = expected + base_point_native.clone();
        }
        assert_eq!(result.get_value(), expected.to_affine());

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_straus_lookup_table_infinity_base_point() {
        let builder = make_builder();
        let offset_gen_native = Grumpkin::random_element();

        let base_point = CycleGroupT::<C>::from_witness(builder.clone(), &GrumpkinAffine::infinity());
        let offset_gen = CycleGroupT::<C>::from_witness(builder.clone(), &offset_gen_native.to_affine());

        let table_bits = 2;
        let table = StrausLookupTable::<C>::new(
            builder.clone(),
            &base_point,
            &offset_gen,
            table_bits,
            None,
        );

        // All entries should equal the offset generator
        let table_size = 1usize << table_bits;
        for i in 0..table_size {
            let index = FieldT::from_witness(builder.clone(), Fr::from(i as u64));
            let result = table.read(&index);
            assert_eq!(
                result.get_value(),
                offset_gen_native.to_affine(),
                "Index {} should be offset generator",
                i
            );
        }

        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: constructor tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_cycle_group_infinity() {
        let inf = CycleGroupT::<C>::infinity(None);
        assert!(inf.is_point_at_infinity().get_value());
        assert!(inf.is_constant());
    }

    #[test]
    fn test_cycle_group_from_affine() {
        let points = generate_test_points(1);
        let cg = CycleGroupT::<C>::from_affine(points[0].clone());
        assert_eq!(cg.get_value(), points[0]);
        assert!(cg.is_constant());
        assert!(!cg.is_point_at_infinity().get_value());
    }

    #[test]
    fn test_cycle_group_from_witness() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let cg = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        assert_eq!(cg.get_value(), points[0]);
        assert!(!cg.is_constant());
        assert!(!cg.is_point_at_infinity().get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_cycle_group_inf_witness_regression() {
        let builder = make_builder();
        let inf_affine = GrumpkinAffine::infinity();
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine);
        assert!(a.is_point_at_infinity().get_value());
        assert!(!builder.borrow().base.failed());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_cycle_group_inf_constant_witness_regression() {
        let builder = make_builder();
        let inf_affine = GrumpkinAffine::infinity();
        let a = CycleGroupT::<C>::from_constant_witness(builder.clone(), &inf_affine);
        assert!(a.is_point_at_infinity().get_value());
        assert!(!builder.borrow().base.failed());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_cycle_group_one() {
        let builder = make_builder();
        let one = CycleGroupT::<C>::one(Some(builder));
        let expected_one = GrumpkinAffine::one();
        let one_native = one.get_value();
        assert_eq!(one_native.x, expected_one.x);
        assert_eq!(one_native.y, expected_one.y);
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: validate_on_curve tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_validate_on_curve_succeed() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let x = FieldT::from_witness(builder.clone(), points[0].x.clone());
        let y = FieldT::from_witness(builder.clone(), points[0].y.clone());
        let is_inf = BoolT::from_witness(&WitnessT::from_bool(builder.clone(), false));

        let _point = CycleGroupT::<C>::from_xy_bool(x, y, is_inf, true);
        assert!(!builder.borrow().base.failed());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_validate_on_curve_infinity_succeed() {
        let builder = make_builder();
        let x = FieldT::from_witness(builder.clone(), Fr::from(1u64));
        let y = FieldT::from_witness(builder.clone(), Fr::from(1u64));

        // Mark as infinity — should pass regardless of (x,y) values
        let _point = CycleGroupT::<C>::from_xy_bool(x, y, BoolT::from_constant(true), true);
        assert!(!builder.borrow().base.failed());
    }

    #[test]
    fn test_validate_on_curve_fail() {
        let builder = make_builder();
        let x = FieldT::from_witness(builder.clone(), Fr::from(1u64));
        let y = FieldT::from_witness(builder.clone(), Fr::from(1u64));

        // (1,1) is not on Grumpkin
        let _point = CycleGroupT::<C>::from_xy_bool(x, y, BoolT::from_constant(false), true);
        assert!(builder.borrow().base.failed());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: dbl tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_dbl_witness() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        let c = a.dbl(None);
        let expected = GrumpkinAffine::from(Grumpkin::from_affine(&points[0]).dbl());
        assert_eq!(c.get_value(), expected);
        let result = check_circuit(&builder);
        assert!(result.is_ok(), "circuit check failed: {:?}", result.err());
    }

    #[test]
    fn test_dbl_constant() {
        let points = generate_test_points(1);
        let a = CycleGroupT::<C>::from_affine(points[0].clone());

        let c = a.dbl(None);
        let expected = GrumpkinAffine::from(Grumpkin::from_affine(&points[0]).dbl());
        assert_eq!(c.get_value(), expected);
        assert!(c.is_constant());
    }

    #[test]
    fn test_dbl_with_hint() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        let doubled = Grumpkin::from_affine(&points[0]).dbl();
        let hint = GrumpkinAffine::from(doubled);

        let result = a.dbl(Some(hint.clone()));
        assert_eq!(result.get_value(), hint);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_dbl_infinity_witness() {
        let builder = make_builder();
        let inf_affine = GrumpkinAffine::infinity();
        let infinity = CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine);

        let result = infinity.dbl(None);
        assert!(result.is_point_at_infinity().get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_dbl_infinity_constant() {
        let infinity = CycleGroupT::<C>::infinity(None);
        let result = infinity.dbl(None);
        assert!(result.is_point_at_infinity().get_value());
        assert!(result.is_constant());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: unconditional_add tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_unconditional_add_witness_witness() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        let c = a.unconditional_add(&b, None);
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) + Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(c.get_value(), expected);
        assert!(!c.is_point_at_infinity().get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_unconditional_add_with_hint() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        let sum = Grumpkin::from_affine(&points[0]) + Grumpkin::from_affine(&points[1]);
        let hint = GrumpkinAffine::from(sum);

        let c = a.unconditional_add(&b, Some(hint.clone()));
        assert_eq!(c.get_value(), hint);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_unconditional_add_mixed_witness_constant() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_affine(points[1].clone()); // constant

        let c = a.unconditional_add(&b, None);
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) + Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(c.get_value(), expected);
        assert!(!c.is_constant());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_unconditional_add_constant_constant() {
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_affine(points[0].clone());
        let b = CycleGroupT::<C>::from_affine(points[1].clone());

        let c = a.unconditional_add(&b, None);
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) + Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(c.get_value(), expected);
        assert!(c.is_constant());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: unconditional_subtract tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_unconditional_subtract_witness_witness() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        let c = a.unconditional_subtract(&b, None);
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) - Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(c.get_value(), expected);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_unconditional_subtract_constant_constant() {
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_affine(points[0].clone());
        let b = CycleGroupT::<C>::from_affine(points[1].clone());

        let c = a.unconditional_subtract(&b, None);
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) - Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(c.get_value(), expected);
        assert!(c.is_constant());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: checked_unconditional_add/sub tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_checked_unconditional_add_succeed() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        let c = a.checked_unconditional_add(&b, None);
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) + Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(c.get_value(), expected);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_checked_unconditional_add_fail() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let neg_point = {
            let mut p = points[0].clone();
            p.y = -p.y;
            p
        };
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &neg_point);

        let _c = a.checked_unconditional_add(&b, None);
        assert!(builder.borrow().base.failed());
    }

    #[test]
    fn test_checked_unconditional_subtract_succeed() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        let c = a.checked_unconditional_subtract(&b, None);
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) - Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(c.get_value(), expected);
        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: complete_add (operator+) tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_add_regular() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        let c = &a + &b;
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) + Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(c.get_value(), expected);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_lhs_infinity() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let inf_affine = GrumpkinAffine::infinity();

        let a = CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        let c = &a + &b;
        assert_eq!(c.get_value(), points[0]);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_rhs_infinity() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let inf_affine = GrumpkinAffine::infinity();

        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine);

        let c = &a + &b;
        assert_eq!(c.get_value(), points[0]);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_both_infinity() {
        let builder = make_builder();
        let inf_affine = GrumpkinAffine::infinity();

        let a = CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine);

        let c = &a + &b;
        assert!(c.is_point_at_infinity().get_value());
        assert!(c.get_value().is_point_at_infinity());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_inverse_points() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let neg_point = {
            let mut p = points[0].clone();
            p.y = -p.y;
            p
        };

        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &neg_point);

        let c = &a + &b;
        assert!(c.is_point_at_infinity().get_value());
        assert!(c.get_value().is_point_at_infinity());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_doubling() {
        let builder = make_builder();
        let points = generate_test_points(1);

        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        let c = &a + &b;
        let expected = GrumpkinAffine::from(Grumpkin::from_affine(&points[0]).dbl());
        assert_eq!(c.get_value(), expected);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_add_constant_points() {
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_affine(points[0].clone());
        let b = CycleGroupT::<C>::from_affine(points[1].clone());

        let result = &a + &b;
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) + Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(result.get_value(), expected);
        assert!(result.is_constant());
    }

    #[test]
    fn test_add_constant_plus_constant_infinity() {
        let builder = make_builder();
        let points = generate_test_points(1);

        let a = CycleGroupT::<C>::from_affine(points[0].clone());
        let b = CycleGroupT::<C>::infinity(Some(builder));

        let result = &a + &b;
        assert_eq!(result.get_value(), points[0]);
        assert!(result.is_constant());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: complete_sub (operator-) tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_subtract_regular() {
        let builder = make_builder();
        let points = generate_test_points(2);

        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        let c = &a - &b;
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) - Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(c.get_value(), expected);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_subtract_lhs_infinity() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let inf_affine = GrumpkinAffine::infinity();

        let a = CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        let c = &a - &b;
        let neg_point = {
            let mut p = points[0].clone();
            p.y = -p.y;
            p
        };
        assert_eq!(c.get_value(), neg_point);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_subtract_rhs_infinity() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let inf_affine = GrumpkinAffine::infinity();

        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine);

        let c = &a - &b;
        assert_eq!(c.get_value(), points[0]);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_subtract_both_infinity() {
        let builder = make_builder();
        let inf_affine = GrumpkinAffine::infinity();

        let a = CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine);

        let c = &a - &b;
        assert!(c.is_point_at_infinity().get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_subtract_same_point() {
        let builder = make_builder();
        let points = generate_test_points(1);

        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        let c = &a - &b;
        assert!(c.is_point_at_infinity().get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_subtract_constant_points() {
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_affine(points[0].clone());
        let b = CycleGroupT::<C>::from_affine(points[1].clone());

        let result = &a - &b;
        let expected = GrumpkinAffine::from(
            Grumpkin::from_affine(&points[0]) - Grumpkin::from_affine(&points[1]),
        );
        assert_eq!(result.get_value(), expected);
        assert!(result.is_constant());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: negate tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_negate() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        let neg_a = a.negate();
        let expected = {
            let mut p = points[0].clone();
            p.y = -p.y;
            p
        };
        assert_eq!(neg_a.get_value(), expected);
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_operator_neg_regression() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let mut b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        b = b.negate();
        let c = a.unconditional_add(&b, None);
        let _ = c;
        assert!(!builder.borrow().base.failed());
        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: witness sum regression
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_witness_sum_regression() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        let c = &a + &b;
        assert!(!c.is_constant());

        let d = &a - &b;
        assert!(!d.is_constant());

        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: conditional_assign tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_conditional_assign() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        // predicate = true → select a
        let pred_true = BoolT::from_witness(&WitnessT::from_bool(builder.clone(), true));
        let result = CycleGroupT::<C>::conditional_assign(&pred_true, &a, &b);
        assert_eq!(result.get_value(), points[0]);

        // predicate = false → select b
        let pred_false = BoolT::from_witness(&WitnessT::from_bool(builder.clone(), false));
        let result = CycleGroupT::<C>::conditional_assign(&pred_false, &a, &b);
        assert_eq!(result.get_value(), points[1]);

        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_conditional_assign_regression() {
        let builder = make_builder();
        let inf_affine = GrumpkinAffine::infinity();
        let c0 = CycleGroupT::<C>::from_affine(inf_affine);
        let pred = BoolT::from_witness(&WitnessT::from_bool(builder.clone(), false));
        let c1 = CycleGroupT::<C>::conditional_assign(&pred, &c0, &c0);
        let w3 = c1.dbl(None);
        let _ = w3;
        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: scalar_mul tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_scalar_mul_witness_point_witness_scalar() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let scalar_val = ScalarField::random_element();

        let point = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let scalar = CycleScalarT::<C>::from_witness(builder.clone(), scalar_val);

        let result = point.scalar_mul(&scalar);

        let expected = Grumpkin::from_affine(&points[0]).mul_without_endomorphism(&scalar_val);
        assert_eq!(result.get_value(), expected.to_affine());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_scalar_mul_constant_point_witness_scalar() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let scalar_val = ScalarField::random_element();

        let point = CycleGroupT::<C>::from_affine(points[0].clone());
        let scalar = CycleScalarT::<C>::from_witness(builder.clone(), scalar_val);

        let result = point.scalar_mul(&scalar);

        let expected = Grumpkin::from_affine(&points[0]).mul_without_endomorphism(&scalar_val);
        assert_eq!(result.get_value(), expected.to_affine());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_scalar_mul_witness_point_constant_scalar() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let scalar_val = ScalarField::random_element();

        let point = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let scalar = CycleScalarT::<C>::from_native(scalar_val);

        let result = point.scalar_mul(&scalar);

        let expected = Grumpkin::from_affine(&points[0]).mul_without_endomorphism(&scalar_val);
        assert_eq!(result.get_value(), expected.to_affine());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_scalar_mul_constant_point_constant_scalar() {
        let points = generate_test_points(1);
        let scalar_val = ScalarField::random_element();

        let point = CycleGroupT::<C>::from_affine(points[0].clone());
        let scalar = CycleScalarT::<C>::from_native(scalar_val);

        let result = point.scalar_mul(&scalar);

        let expected = Grumpkin::from_affine(&points[0]).mul_without_endomorphism(&scalar_val);
        assert_eq!(result.get_value(), expected.to_affine());
    }

    #[test]
    fn test_scalar_mul_by_zero() {
        let builder = make_builder();
        let points = generate_test_points(1);

        let point = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let scalar = CycleScalarT::<C>::from_witness(builder.clone(), ScalarField::zero());

        let result = point.scalar_mul(&scalar);
        assert!(result.is_point_at_infinity().get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: batch_mul tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_batch_mul_general_msm() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let scalar_val = ScalarField::random_element();

        let mut msm_points = Vec::new();
        let mut msm_scalars = Vec::new();
        let mut expected = Grumpkin::infinity();

        // 1: witness point, witness scalar
        expected = expected + Grumpkin::from_affine(&points[0]).mul_without_endomorphism(&scalar_val);
        msm_points.push(CycleGroupT::<C>::from_witness(builder.clone(), &points[0]));
        msm_scalars.push(CycleScalarT::<C>::from_witness(builder.clone(), scalar_val));

        // 2: constant point, witness scalar
        expected = expected + Grumpkin::from_affine(&points[0]).mul_without_endomorphism(&scalar_val);
        msm_points.push(CycleGroupT::<C>::from_affine(points[0].clone()));
        msm_scalars.push(CycleScalarT::<C>::from_witness(builder.clone(), scalar_val));

        // 3: witness point, constant scalar
        expected = expected + Grumpkin::from_affine(&points[0]).mul_without_endomorphism(&scalar_val);
        msm_points.push(CycleGroupT::<C>::from_witness(builder.clone(), &points[0]));
        msm_scalars.push(CycleScalarT::<C>::from_native(scalar_val));

        // 4: constant point, constant scalar
        expected = expected + Grumpkin::from_affine(&points[0]).mul_without_endomorphism(&scalar_val);
        msm_points.push(CycleGroupT::<C>::from_affine(points[0].clone()));
        msm_scalars.push(CycleScalarT::<C>::from_native(scalar_val));

        let result = CycleGroupT::<C>::batch_mul(&msm_points, &msm_scalars);
        assert_eq!(result.get_value(), expected.to_affine());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_batch_mul_produces_infinity() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let scalar_val = ScalarField::random_element();

        let msm_points = vec![
            CycleGroupT::<C>::from_witness(builder.clone(), &points[0]),
            CycleGroupT::<C>::from_witness(builder.clone(), &points[0]),
        ];
        let neg_scalar = -scalar_val;
        let msm_scalars = vec![
            CycleScalarT::<C>::from_witness(builder.clone(), scalar_val),
            CycleScalarT::<C>::from_witness(builder.clone(), neg_scalar),
        ];

        let result = CycleGroupT::<C>::batch_mul(&msm_points, &msm_scalars);
        assert!(result.is_point_at_infinity().get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_batch_mul_multiply_by_zero() {
        let builder = make_builder();
        let points = generate_test_points(1);

        let msm_points = vec![
            CycleGroupT::<C>::from_witness(builder.clone(), &points[0]),
        ];
        let msm_scalars = vec![
            CycleScalarT::<C>::from_witness(builder.clone(), ScalarField::zero()),
        ];

        let result = CycleGroupT::<C>::batch_mul(&msm_points, &msm_scalars);
        assert!(result.is_point_at_infinity().get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_batch_mul_infinity_inputs() {
        let builder = make_builder();
        let scalar_val = ScalarField::random_element();
        let inf_affine = GrumpkinAffine::infinity();

        // Witness infinity
        let msm_points = vec![
            CycleGroupT::<C>::from_witness(builder.clone(), &inf_affine),
            CycleGroupT::<C>::from_affine(inf_affine), // constant infinity
        ];
        let msm_scalars = vec![
            CycleScalarT::<C>::from_witness(builder.clone(), scalar_val),
            CycleScalarT::<C>::from_witness(builder.clone(), scalar_val),
        ];

        let result = CycleGroupT::<C>::batch_mul(&msm_points, &msm_scalars);
        assert!(result.is_point_at_infinity().get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: eq and assert_equal tests
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_eq_same_point() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let mut a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let mut b = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        let result = a.eq(&mut b);
        assert!(result.get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_eq_different_points() {
        let builder = make_builder();
        let points = generate_test_points(2);
        let mut a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let mut b = CycleGroupT::<C>::from_witness(builder.clone(), &points[1]);

        let result = a.eq(&mut b);
        assert!(!result.get_value());
        assert!(check_circuit(&builder).is_ok());
    }

    #[test]
    fn test_assert_equal() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let mut a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);
        let mut b = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        a.assert_equal(&mut b, "test_assert_equal");
        assert!(!builder.borrow().base.failed());
        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: constant/witness mixup regression
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_constant_witness_mixup_regression() {
        let builder = make_builder();
        let points = generate_test_points(1);

        let c1 = CycleGroupT::<C>::from_affine(GrumpkinAffine::one());
        let cw8 = CycleGroupT::<C>::from_constant_witness(builder.clone(), &GrumpkinAffine::infinity());
        let w11 = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        let w9 = &cw8 + &c1;
        let _w26 = &w9 + &w11;

        let w10 = &cw8 - &c1;
        let _w27 = &w10 - &w11;

        assert!(check_circuit(&builder).is_ok());
    }

    // ════════════════════════════════════════════════════════════════════════
    //  CycleGroup: set_public test
    // ════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_set_public() {
        let builder = make_builder();
        let points = generate_test_points(1);
        let a = CycleGroupT::<C>::from_witness(builder.clone(), &points[0]);

        let num_public = a.set_public();
        assert_eq!(num_public, 3); // x, y, is_infinity
        assert!(check_circuit(&builder).is_ok());
    }
}
