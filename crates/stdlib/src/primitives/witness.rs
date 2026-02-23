//! Circuit witness primitives.
//!
//! Port of `barretenberg/stdlib/primitives/witness/witness.hpp`.

use std::cell::RefCell;
use std::rc::Rc;

use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

/// Builder reference type shared across circuit elements.
pub type BuilderRef<P> = Rc<RefCell<UltraCircuitBuilder<P>>>;

/// Sentinel value indicating a constant (no witness).
pub const IS_CONSTANT: u32 = u32::MAX;

/// A witness element in a circuit.
///
/// Wraps a field value and its witness index in the circuit builder.
/// Port of C++ `witness_t<Builder>`.
#[derive(Clone)]
pub struct WitnessT<P: FieldParams> {
    pub witness: Field<P>,
    pub witness_index: u32,
    pub context: BuilderRef<P>,
}

impl<P: FieldParams> WitnessT<P> {
    /// Create a witness from a field element.
    pub fn new(context: BuilderRef<P>, value: Field<P>) -> Self {
        let witness_index = context.borrow_mut().base.add_variable(value);
        Self {
            witness: value,
            witness_index,
            context,
        }
    }

    /// Create a witness from a boolean value.
    pub fn from_bool(context: BuilderRef<P>, value: bool) -> Self {
        let witness = if value { Field::one() } else { Field::zero() };
        Self::new(context, witness)
    }

    /// Create a witness from a u64 value.
    pub fn from_u64(context: BuilderRef<P>, value: u64) -> Self {
        Self::new(context, Field::from(value))
    }

    /// Create a constant witness: the witness value is constrained to equal the
    /// given constant via `assert_equal_constant`.
    pub fn create_constant_witness(context: BuilderRef<P>, value: Field<P>) -> Self {
        let out = Self::new(context.clone(), value);
        context.borrow_mut().assert_equal_constant(
            out.witness_index,
            value,
            "Failed to create constant witness.",
        );
        out
    }

    /// Check whether this witness represents a constant.
    pub fn is_constant(&self) -> bool {
        self.witness_index == IS_CONSTANT
    }
}

/// A public witness element in a circuit.
///
/// Like `WitnessT` but registers the variable as a public input.
/// Port of C++ `public_witness_t<Builder>`.
#[derive(Clone)]
pub struct PublicWitnessT<P: FieldParams> {
    pub witness: Field<P>,
    pub witness_index: u32,
    pub context: BuilderRef<P>,
}

impl<P: FieldParams> PublicWitnessT<P> {
    /// Create a public witness from a field element.
    pub fn new(context: BuilderRef<P>, value: Field<P>) -> Self {
        let witness_index = context.borrow_mut().base.add_public_variable(value);
        Self {
            witness: value,
            witness_index,
            context,
        }
    }

    /// Create a public witness from a boolean value.
    pub fn from_bool(context: BuilderRef<P>, value: bool) -> Self {
        let witness = if value { Field::one() } else { Field::zero() };
        Self::new(context, witness)
    }

    /// Create a public witness from a u64 value.
    pub fn from_u64(context: BuilderRef<P>, value: u64) -> Self {
        Self::new(context, Field::from(value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_ecc::curves::bn254::Bn254FrParams;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    #[test]
    fn test_witness_constructor_from_field() {
        let builder = make_builder();
        let value = Field::from(42u64);
        let w = WitnessT::new(builder.clone(), value);
        assert_eq!(w.witness, value);
        assert!(!w.is_constant());
        assert_eq!(builder.borrow().base.get_variable(w.witness_index), value);
    }

    #[test]
    fn test_witness_constructor_from_bool() {
        let builder = make_builder();
        let w_true = WitnessT::from_bool(builder.clone(), true);
        let w_false = WitnessT::from_bool(builder.clone(), false);
        assert_eq!(w_true.witness, Field::one());
        assert_eq!(w_false.witness, Field::zero());
    }

    #[test]
    fn test_witness_constructor_from_u64() {
        let builder = make_builder();
        let w = WitnessT::from_u64(builder.clone(), 99);
        assert_eq!(w.witness, Field::from(99u64));
        assert_eq!(
            builder.borrow().base.get_variable(w.witness_index),
            Field::from(99u64)
        );
    }

    #[test]
    fn test_witness_create_constant() {
        let builder = make_builder();
        let value = Field::from(7u64);
        let w = WitnessT::create_constant_witness(builder.clone(), value);
        assert_eq!(w.witness, value);
        assert!(!w.is_constant());
    }

    #[test]
    fn test_witness_is_constant_sentinel() {
        // Verify the IS_CONSTANT sentinel
        assert_eq!(IS_CONSTANT, u32::MAX);
    }

    #[test]
    fn test_public_witness_constructor() {
        let builder = make_builder();
        let value = Field::from(123u64);
        let pw = PublicWitnessT::new(builder.clone(), value);
        assert_eq!(pw.witness, value);
        assert_eq!(
            builder.borrow().base.get_variable(pw.witness_index),
            value
        );
        assert_eq!(builder.borrow().base.num_public_inputs(), 1);
    }
}
