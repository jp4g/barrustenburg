//! Runtime-defined read-write memory table.
//!
//! Port of `barretenberg/stdlib/primitives/memory/ram_table.hpp` and `ram_table.cpp`.

use bbrs_ecc::fields::field::Field;
use bbrs_ecc::fields::field_params::FieldParams;

use crate::primitives::field::FieldT;
use crate::primitives::witness::BuilderRef;

/// A runtime-defined read-write memory table.
///
/// Table entries must be initialized before reads. Supports read and write
/// operations via witness indices.
///
/// Port of C++ `ram_table<Builder>`.
pub struct RamTable<P: FieldParams> {
    raw_entries: Vec<FieldT<P>>,
    index_initialized: Vec<bool>,
    length: usize,
    ram_id: usize,
    ram_table_generated_in_builder: bool,
    all_entries_written_to_with_constant_index: bool,
    context: Option<BuilderRef<P>>,
}

impl<P: FieldParams> RamTable<P> {
    /// Create a RAM table from field elements.
    ///
    /// Extracts builder context from entries if available.
    pub fn new(table_entries: Vec<FieldT<P>>) -> Self {
        let mut context = None;
        for entry in &table_entries {
            if entry.get_context().is_some() {
                context = entry.get_context().clone();
                break;
            }
        }
        let length = table_entries.len();
        Self {
            raw_entries: table_entries,
            index_initialized: vec![false; length],
            length,
            ram_id: 0,
            ram_table_generated_in_builder: false,
            all_entries_written_to_with_constant_index: false,
            context,
        }
    }

    /// Create a RAM table with an explicit builder reference.
    pub fn with_builder(ctx: BuilderRef<P>, table_entries: Vec<FieldT<P>>) -> Self {
        let length = table_entries.len();
        Self {
            raw_entries: table_entries,
            index_initialized: vec![false; length],
            length,
            ram_id: 0,
            ram_table_generated_in_builder: false,
            all_entries_written_to_with_constant_index: false,
            context: Some(ctx),
        }
    }

    /// Initialize the internal RAM array in the builder.
    fn initialize_table(&mut self) {
        if self.ram_table_generated_in_builder {
            return;
        }
        let ctx = self
            .context
            .as_ref()
            .expect("ram_table: context required for initialization")
            .clone();

        self.ram_id = ctx.borrow_mut().create_ram_array(self.length);

        if !self.raw_entries.is_empty() {
            for i in 0..self.length {
                if !self.index_initialized[i] {
                    let entry = if self.raw_entries[i].is_constant() {
                        let value = self.raw_entries[i].get_value();
                        let witness_idx = ctx.borrow_mut().put_constant_variable(value);
                        FieldT::from_witness_index(ctx.clone(), witness_idx)
                    } else {
                        self.raw_entries[i].clone()
                    };
                    let witness_idx = entry.get_witness_index();
                    ctx.borrow_mut()
                        .init_ram_element(self.ram_id, i, witness_idx);
                    self.index_initialized[i] = true;
                }
            }
        }

        self.ram_table_generated_in_builder = true;
    }

    /// Check whether all entries have been initialized with constant indices.
    fn check_indices_initialized(&mut self) -> bool {
        if self.all_entries_written_to_with_constant_index {
            return true;
        }
        if self.length == 0 {
            return false;
        }
        let all_init = self.index_initialized.iter().all(|&x| x);
        self.all_entries_written_to_with_constant_index = all_init;
        all_init
    }

    /// Read a field element from the RAM table.
    pub fn read(&mut self, index: &FieldT<P>) -> FieldT<P> {
        if self.context.is_none() {
            self.context = index.get_context().clone();
            assert!(
                self.context.is_some(),
                "ram_table: cannot read without a builder context"
            );
        }

        self.initialize_table();

        let native_index = index.get_value().from_montgomery_form().data[0] as usize;
        if native_index >= self.length {
            self.context
                .as_ref()
                .unwrap()
                .borrow_mut()
                .base
                .failure("ram_table: RAM array access out of bounds".to_string());
        }

        if !self.check_indices_initialized() {
            self.context
                .as_ref()
                .unwrap()
                .borrow_mut()
                .base
                .failure("ram_table must have initialized every RAM entry before reading".to_string());
        }

        let ctx = self.context.as_ref().unwrap().clone();

        // Constant index: convert to fixed witness
        let index_witness = if index.is_constant() {
            let value = index.get_value();
            let idx = ctx.borrow_mut().put_constant_variable(value);
            FieldT::from_witness_index(ctx.clone(), idx)
        } else {
            index.clone()
        };

        let wi = index_witness.get_witness_index();
        let output_idx = ctx.borrow_mut().read_ram_array(self.ram_id, wi);

        FieldT::from_witness_index(ctx, output_idx)
    }

    /// Read from the table with a constant index.
    pub fn read_const(&mut self, index: usize) -> FieldT<P> {
        let field_index = FieldT::from_field(Field::from(index as u64));
        self.read(&field_index)
    }

    /// Write a field element to the RAM table.
    pub fn write(&mut self, index: &FieldT<P>, value: &FieldT<P>) {
        if self.context.is_none() {
            self.context = index.get_context().clone();
            assert!(
                self.context.is_some(),
                "ram_table: cannot write without a builder context"
            );
        }

        self.initialize_table();

        let native_index = index.get_value().from_montgomery_form().data[0] as usize;
        if native_index >= self.length {
            self.context
                .as_ref()
                .unwrap()
                .borrow_mut()
                .base
                .failure("ram_table: RAM array access out of bounds".to_string());
        }

        let ctx = self.context.as_ref().unwrap().clone();

        let index_wire = if index.is_constant() {
            let mut iw = index.clone();
            iw.convert_constant_to_fixed_witness(ctx.clone());
            iw
        } else {
            if !self.check_indices_initialized() {
                ctx.borrow_mut().base.failure(
                    "ram_table must have initialized every RAM entry before a write".to_string(),
                );
            }
            index.clone()
        };

        let value_wire = if value.is_constant() {
            let val = value.get_value();
            let idx = ctx.borrow_mut().put_constant_variable(val);
            FieldT::from_witness_index(ctx.clone(), idx)
        } else {
            value.clone()
        };

        let cast_index = native_index;
        if index.is_constant() && !self.index_initialized[cast_index] {
            let val_wi = value_wire.get_witness_index();
            ctx.borrow_mut()
                .init_ram_element(self.ram_id, cast_index, val_wi);
            self.index_initialized[cast_index] = true;
        } else {
            let idx_wi = index_wire.get_witness_index();
            let val_wi = value_wire.get_witness_index();
            ctx.borrow_mut()
                .write_ram_array(self.ram_id, idx_wi, val_wi);
        }
    }

    /// Write with a constant index.
    pub fn write_const(&mut self, index: usize, value: &FieldT<P>) {
        let field_index = FieldT::from_field(Field::from(index as u64));
        self.write(&field_index, value);
    }

    /// Number of entries in the table.
    pub fn size(&self) -> usize {
        self.length
    }

    /// Get the builder context.
    pub fn get_context(&self) -> &Option<BuilderRef<P>> {
        &self.context
    }
}

impl<P: FieldParams> Clone for RamTable<P> {
    fn clone(&self) -> Self {
        Self {
            raw_entries: self.raw_entries.clone(),
            index_initialized: self.index_initialized.clone(),
            length: self.length,
            ram_id: self.ram_id,
            ram_table_generated_in_builder: self.ram_table_generated_in_builder,
            all_entries_written_to_with_constant_index: self
                .all_entries_written_to_with_constant_index,
            context: self.context.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bbrs_circuit_builder::circuit_checker::UltraCircuitChecker;
    use bbrs_circuit_builder::ultra_builder::UltraCircuitBuilder;
    use bbrs_ecc::curves::bn254::Bn254FrParams;
    use std::cell::RefCell;
    use std::rc::Rc;

    type Fr = Field<Bn254FrParams>;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    #[test]
    fn test_ram_table_init_read_consistency() {
        let ctx = make_builder();
        let table_size = 10;
        let mut table_values = Vec::new();
        let mut expected_values = Vec::new();

        for _ in 0..table_size {
            let val = Fr::random_element();
            expected_values.push(val);
            table_values.push(FieldT::from_witness(ctx.clone(), val));
        }

        let mut table = RamTable::new(table_values);

        let mut result = Fr::zero();
        let mut expected = Fr::zero();

        for i in 0..table_size {
            if i % 2 == 0 {
                // Variable lookup
                let index = FieldT::from_witness(ctx.clone(), Fr::from(i as u64));
                let to_add = table.read(&index);
                result = result + to_add.get_value();
            } else {
                // Constant lookup
                let to_add = table.read_const(i);
                result = result + to_add.get_value();
            }
            expected = expected + expected_values[i];
        }

        assert_eq!(result, expected);

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }

    #[test]
    fn test_ram_table_read_write_consistency() {
        let ctx = make_builder();
        let table_size = 10;

        // Initialize with zeros
        let zeros: Vec<FieldT<Bn254FrParams>> =
            vec![FieldT::from_field(Fr::zero()); table_size];
        let mut table = RamTable::with_builder(ctx.clone(), zeros);

        // Write initial values via constant indices
        for i in 0..table_size {
            let val = FieldT::from_field(Fr::zero());
            table.write_const(i, &val);
        }

        let mut table_values = vec![Fr::zero(); table_size];
        let mut result = Fr::zero();
        let mut expected = Fr::zero();

        // Write+read cycle
        let update = |table: &mut RamTable<Bn254FrParams>,
                      table_values: &mut Vec<Fr>,
                      ctx: &BuilderRef<Bn254FrParams>| {
            for i in 0..table_size / 2 {
                let val1 = Fr::random_element();
                let val2 = Fr::random_element();
                table_values[2 * i] = val1;
                table_values[2 * i + 1] = val2;
                // Constant-index write
                table.write_const(2 * i, &FieldT::from_field(val1));
                // Variable-index write
                table.write_const(2 * i + 1, &FieldT::from_witness(ctx.clone(), val2));
            }
        };

        update(&mut table, &mut table_values, &ctx);

        // Read from table
        for i in 0..table_size {
            let index = FieldT::from_witness(ctx.clone(), Fr::from(i as u64));
            let read_val = table.read(&index);
            result = result + read_val.get_value();
            expected = expected + table_values[i];
        }

        assert_eq!(result, expected);

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }

    #[test]
    fn test_ram_table_write_then_read() {
        let ctx = make_builder();

        let entries: Vec<FieldT<Bn254FrParams>> = vec![
            FieldT::from_witness(ctx.clone(), Fr::from(10u64)),
            FieldT::from_witness(ctx.clone(), Fr::from(20u64)),
            FieldT::from_witness(ctx.clone(), Fr::from(30u64)),
        ];

        let mut table = RamTable::new(entries);

        // Read original
        let idx1 = FieldT::from_witness(ctx.clone(), Fr::from(1u64));
        let val = table.read(&idx1);
        assert_eq!(val.get_value(), Fr::from(20u64));

        // Overwrite entry 1
        let new_val = FieldT::from_witness(ctx.clone(), Fr::from(99u64));
        let idx1_write = FieldT::from_witness(ctx.clone(), Fr::from(1u64));
        table.write(&idx1_write, &new_val);

        // Read back
        let idx1_read = FieldT::from_witness(ctx.clone(), Fr::from(1u64));
        let val2 = table.read(&idx1_read);
        assert_eq!(val2.get_value(), Fr::from(99u64));

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }

    #[test]
    fn test_ram_table_with_builder() {
        let ctx = make_builder();
        let entries: Vec<FieldT<Bn254FrParams>> = vec![
            FieldT::from_witness(ctx.clone(), Fr::from(5u64)),
            FieldT::from_witness(ctx.clone(), Fr::from(15u64)),
            FieldT::from_witness(ctx.clone(), Fr::from(25u64)),
        ];

        let mut table = RamTable::with_builder(ctx.clone(), entries);

        let idx = FieldT::from_witness(ctx.clone(), Fr::from(2u64));
        let val = table.read(&idx);
        assert_eq!(val.get_value(), Fr::from(25u64));

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }

    #[test]
    fn test_ram_table_constant_entries() {
        let ctx = make_builder();
        let entries: Vec<FieldT<Bn254FrParams>> = vec![
            FieldT::from_field(Fr::from(100u64)),
            FieldT::from_field(Fr::from(200u64)),
            FieldT::from_field(Fr::from(300u64)),
        ];

        let mut table = RamTable::with_builder(ctx.clone(), entries);

        let idx = FieldT::from_witness(ctx.clone(), Fr::from(0u64));
        let val = table.read(&idx);
        assert_eq!(val.get_value(), Fr::from(100u64));

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }
}
