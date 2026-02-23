//! Runtime-defined read-only memory table.
//!
//! Port of `barretenberg/stdlib/primitives/memory/rom_table.hpp` and `rom_table.cpp`.

use bbrs_ecc::fields::field_params::FieldParams;

use crate::primitives::field::FieldT;
use crate::primitives::witness::BuilderRef;

/// A runtime-defined read-only memory table.
///
/// Table entries must be initialized in the constructor. Supports both
/// constant lookups (zero gates) and witness-index lookups (adds gates).
///
/// Port of C++ `rom_table<Builder>`.
pub struct RomTable<P: FieldParams> {
    raw_entries: Vec<FieldT<P>>,
    entries: Vec<FieldT<P>>,
    length: usize,
    rom_id: usize,
    initialized: bool,
    context: Option<BuilderRef<P>>,
}

impl<P: FieldParams> RomTable<P> {
    /// Create a ROM table from field elements.
    ///
    /// Extracts builder context from entries if available.
    /// Table initialization is deferred until the first read.
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
            entries: Vec::new(),
            length,
            rom_id: 0,
            initialized: false,
            context,
        }
    }

    /// Create a ROM table with an explicit builder reference.
    pub fn with_builder(ctx: BuilderRef<P>, table_entries: Vec<FieldT<P>>) -> Self {
        let length = table_entries.len();
        Self {
            raw_entries: table_entries,
            entries: Vec::new(),
            length,
            rom_id: 0,
            initialized: false,
            context: Some(ctx),
        }
    }

    /// Initialize the internal ROM array in the builder.
    ///
    /// Converts constant entries to fixed witnesses and registers
    /// all entries with the circuit builder.
    fn initialize_table(&mut self) {
        if self.initialized {
            return;
        }
        let ctx = self
            .context
            .as_ref()
            .expect("rom_table: context required for initialization")
            .clone();

        self.entries.clear();
        for entry in &self.raw_entries {
            if entry.is_constant() {
                let value = entry.get_value();
                let witness_idx = ctx.borrow_mut().put_constant_variable(value);
                self.entries
                    .push(FieldT::from_witness_index(ctx.clone(), witness_idx));
            } else {
                self.entries.push(entry.clone());
            }
        }

        self.rom_id = ctx.borrow_mut().create_rom_array(self.length);

        for i in 0..self.length {
            let witness_idx = self.entries[i].get_witness_index();
            ctx.borrow_mut()
                .set_rom_element(self.rom_id, i, witness_idx);
        }

        self.initialized = true;
    }

    /// Read from the table with a constant index. Does not add any gates.
    pub fn read_const(&mut self, index: usize) -> FieldT<P> {
        if index >= self.length {
            if let Some(ctx) = &self.context {
                ctx.borrow_mut()
                    .base
                    .failure("rom_table: ROM array access out of bounds".to_string());
            }
        }
        // For constant reads, we need initialized entries
        self.initialize_table();
        self.entries[index].clone()
    }

    /// Read from the table with a witness index.
    pub fn read(&mut self, index: &FieldT<P>) -> FieldT<P> {
        if self.context.is_none() {
            self.context = index.get_context().clone();
            assert!(
                self.context.is_some(),
                "rom_table: cannot read without a builder context"
            );
        }

        self.initialize_table();

        if index.is_constant() {
            let idx_val = index.get_value().from_montgomery_form().data[0] as usize;
            return self.read_const(idx_val);
        }

        let native_index = index.get_value().from_montgomery_form().data[0] as usize;
        if native_index >= self.length {
            self.context
                .as_ref()
                .unwrap()
                .borrow_mut()
                .base
                .failure("rom_table: ROM array access out of bounds".to_string());
        }

        let ctx = self.context.as_ref().unwrap().clone();
        let witness_idx = index.get_witness_index();
        let output_idx = ctx
            .borrow_mut()
            .read_rom_array(self.rom_id, witness_idx);

        FieldT::from_witness_index(ctx, output_idx)
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

impl<P: FieldParams> Clone for RomTable<P> {
    fn clone(&self) -> Self {
        Self {
            raw_entries: self.raw_entries.clone(),
            entries: self.entries.clone(),
            length: self.length,
            rom_id: self.rom_id,
            initialized: self.initialized,
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
    use bbrs_ecc::fields::field::Field;
    use std::cell::RefCell;
    use std::rc::Rc;

    type Fr = Field<Bn254FrParams>;

    fn make_builder() -> BuilderRef<Bn254FrParams> {
        Rc::new(RefCell::new(UltraCircuitBuilder::new()))
    }

    #[test]
    fn test_rom_table_init_read_consistency() {
        let ctx = make_builder();
        let table_size = 10;
        let mut table_values = Vec::new();
        let mut expected_values = Vec::new();

        for _ in 0..table_size {
            let val = Fr::random_element();
            expected_values.push(val);
            table_values.push(FieldT::from_witness(ctx.clone(), val));
        }

        let mut table = RomTable::new(table_values);

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
    fn test_rom_table_read_write_consistency() {
        let ctx = make_builder();
        let table_size = 5;
        let mut values = Vec::new();

        for _ in 0..table_size {
            let val = Fr::random_element();
            values.push(val);
        }

        let table_entries: Vec<FieldT<Bn254FrParams>> = values
            .iter()
            .map(|v| FieldT::from_witness(ctx.clone(), *v))
            .collect();

        let mut table = RomTable::new(table_entries);

        // Read and verify each element
        for i in 0..table_size {
            let index = FieldT::from_witness(ctx.clone(), Fr::from(i as u64));
            let read_val = table.read(&index);
            assert_eq!(read_val.get_value(), values[i]);
        }

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }

    #[test]
    fn test_rom_table_copy_read() {
        let ctx = make_builder();
        let table_size = 5;
        let mut table_values = Vec::new();
        let mut expected_values = Vec::new();

        for _ in 0..table_size {
            let val = Fr::random_element();
            expected_values.push(val);
            table_values.push(FieldT::from_witness(ctx.clone(), val));
        }

        let mut table = RomTable::new(table_values);
        let mut copied_table = table.clone();

        let mut result = Fr::zero();
        let mut expected = Fr::zero();

        for i in 0..table_size {
            let index = FieldT::from_witness(ctx.clone(), Fr::from(i as u64));
            let to_add = if i % 2 == 0 {
                copied_table.read(&index)
            } else {
                table.read(&index)
            };
            result = result + to_add.get_value();
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
    fn test_rom_table_with_builder() {
        let ctx = make_builder();
        let values = vec![Fr::from(10u64), Fr::from(20u64), Fr::from(30u64)];
        let entries: Vec<FieldT<Bn254FrParams>> = values
            .iter()
            .map(|v| FieldT::from_witness(ctx.clone(), *v))
            .collect();

        let mut table = RomTable::with_builder(ctx.clone(), entries);

        let idx = FieldT::from_witness(ctx.clone(), Fr::from(1u64));
        let val = table.read(&idx);
        assert_eq!(val.get_value(), Fr::from(20u64));

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        let result = UltraCircuitChecker::check(&mut builder);
        assert!(
            result.is_ok(),
            "Circuit check failed: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_rom_table_constant_entries() {
        let ctx = make_builder();
        // Create a table with constant entries (no witness)
        let entries: Vec<FieldT<Bn254FrParams>> = vec![
            FieldT::from_field(Fr::from(100u64)),
            FieldT::from_field(Fr::from(200u64)),
            FieldT::from_field(Fr::from(300u64)),
        ];

        let mut table = RomTable::with_builder(ctx.clone(), entries);

        // Read with a witness index
        let idx = FieldT::from_witness(ctx.clone(), Fr::from(2u64));
        let val = table.read(&idx);
        assert_eq!(val.get_value(), Fr::from(300u64));

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }
}
