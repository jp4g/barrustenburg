//! Runtime-defined read-only memory table with paired entries.
//!
//! Port of `barretenberg/stdlib/primitives/memory/twin_rom_table.hpp` and `twin_rom_table.cpp`.

use bbrs_ecc::fields::field_params::FieldParams;

use crate::primitives::field::FieldT;
use crate::primitives::witness::BuilderRef;

/// A pair of field elements stored in a twin ROM table.
pub type FieldPair<P> = [FieldT<P>; 2];

/// A runtime-defined read-only memory table where each entry is a pair of values.
///
/// Port of C++ `twin_rom_table<Builder>`.
pub struct TwinRomTable<P: FieldParams> {
    raw_entries: Vec<FieldPair<P>>,
    entries: Vec<FieldPair<P>>,
    length: usize,
    rom_id: usize,
    initialized: bool,
    context: Option<BuilderRef<P>>,
}

impl<P: FieldParams> TwinRomTable<P> {
    /// Create a twin ROM table from pairs of field elements.
    pub fn new(table_entries: Vec<FieldPair<P>>) -> Self {
        let mut context = None;
        'outer: for entry in &table_entries {
            for field in entry {
                if field.get_context().is_some() {
                    context = field.get_context().clone();
                    break 'outer;
                }
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

    /// Initialize the internal ROM array in the builder.
    fn initialize_table(&mut self) {
        if self.initialized {
            return;
        }
        let ctx = self
            .context
            .as_ref()
            .expect("twin_rom_table: context required for initialization")
            .clone();

        self.entries.clear();
        for entry in &self.raw_entries {
            let first = if entry[0].is_constant() {
                let value = entry[0].get_value();
                let witness_idx = ctx.borrow_mut().put_constant_variable(value);
                FieldT::from_witness_index(ctx.clone(), witness_idx)
            } else {
                entry[0].clone()
            };
            let second = if entry[1].is_constant() {
                let value = entry[1].get_value();
                let witness_idx = ctx.borrow_mut().put_constant_variable(value);
                FieldT::from_witness_index(ctx.clone(), witness_idx)
            } else {
                entry[1].clone()
            };
            self.entries.push([first, second]);
        }

        self.rom_id = ctx.borrow_mut().create_rom_array(self.length);

        for i in 0..self.length {
            let w1 = self.entries[i][0].get_witness_index();
            let w2 = self.entries[i][1].get_witness_index();
            ctx.borrow_mut()
                .set_rom_element_pair(self.rom_id, i, [w1, w2]);
        }

        self.initialized = true;
    }

    /// Read from the table with a constant index. Does not add any gates.
    pub fn read_const(&mut self, index: usize) -> FieldPair<P> {
        if index >= self.length {
            if let Some(ctx) = &self.context {
                ctx.borrow_mut()
                    .base
                    .failure("twin_rom_table: ROM array access out of bounds".to_string());
            }
        }
        self.initialize_table();
        self.entries[index].clone()
    }

    /// Read from the table with a witness index.
    pub fn read(&mut self, index: &FieldT<P>) -> FieldPair<P> {
        if index.is_constant() {
            let idx_val = index.get_value().from_montgomery_form().data[0] as usize;
            return self.read_const(idx_val);
        }

        if self.context.is_none() {
            self.context = index.get_context().clone();
        }

        self.initialize_table();

        let native_index = index.get_value().from_montgomery_form().data[0] as usize;
        if native_index >= self.length {
            self.context
                .as_ref()
                .unwrap()
                .borrow_mut()
                .base
                .failure("twin_rom_table: ROM array access out of bounds".to_string());
        }

        let ctx = self.context.as_ref().unwrap().clone();
        let witness_idx = index.get_witness_index();
        let output_indices = ctx
            .borrow_mut()
            .read_rom_array_pair(self.rom_id, witness_idx);

        [
            FieldT::from_witness_index(ctx.clone(), output_indices[0]),
            FieldT::from_witness_index(ctx, output_indices[1]),
        ]
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

impl<P: FieldParams> Clone for TwinRomTable<P> {
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
    fn test_twin_rom_table_read_write_consistency() {
        let ctx = make_builder();
        let table_size = 10;
        let mut table_values: Vec<FieldPair<Bn254FrParams>> = Vec::new();

        for _ in 0..table_size {
            table_values.push([
                FieldT::from_witness(ctx.clone(), Fr::random_element()),
                FieldT::from_witness(ctx.clone(), Fr::random_element()),
            ]);
        }

        let expected_values: Vec<[Fr; 2]> = table_values
            .iter()
            .map(|pair| [pair[0].get_value(), pair[1].get_value()])
            .collect();

        let mut table = TwinRomTable::new(table_values);

        let mut result = [Fr::zero(), Fr::zero()];
        let mut expected = [Fr::zero(), Fr::zero()];

        for i in 0..table_size {
            if i % 2 == 0 {
                // Variable lookup
                let index = FieldT::from_witness(ctx.clone(), Fr::from(i as u64));
                let to_add = table.read(&index);
                result[0] = result[0] + to_add[0].get_value();
                result[1] = result[1] + to_add[1].get_value();
            } else {
                // Constant lookup
                let to_add = table.read_const(i);
                result[0] = result[0] + to_add[0].get_value();
                result[1] = result[1] + to_add[1].get_value();
            }
            expected[0] = expected[0] + expected_values[i][0];
            expected[1] = expected[1] + expected_values[i][1];
        }

        assert_eq!(result[0], expected[0]);
        assert_eq!(result[1], expected[1]);

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }

    #[test]
    fn test_twin_rom_table_constant_entries() {
        let ctx = make_builder();
        let entries: Vec<FieldPair<Bn254FrParams>> = vec![
            [
                FieldT::from_field(Fr::from(10u64)),
                FieldT::from_field(Fr::from(20u64)),
            ],
            [
                FieldT::from_field(Fr::from(30u64)),
                FieldT::from_field(Fr::from(40u64)),
            ],
        ];

        let mut table = TwinRomTable::new(entries);

        // Need context for variable read
        let index = FieldT::from_witness(ctx.clone(), Fr::from(1u64));
        table.context = Some(ctx.clone());
        let pair = table.read(&index);
        assert_eq!(pair[0].get_value(), Fr::from(30u64));
        assert_eq!(pair[1].get_value(), Fr::from(40u64));

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }

    #[test]
    fn test_twin_rom_table_single_entry() {
        let ctx = make_builder();
        let entries: Vec<FieldPair<Bn254FrParams>> = vec![[
            FieldT::from_witness(ctx.clone(), Fr::from(42u64)),
            FieldT::from_witness(ctx.clone(), Fr::from(84u64)),
        ]];

        let mut table = TwinRomTable::new(entries);

        let index = FieldT::from_witness(ctx.clone(), Fr::from(0u64));
        let pair = table.read(&index);
        assert_eq!(pair[0].get_value(), Fr::from(42u64));
        assert_eq!(pair[1].get_value(), Fr::from(84u64));

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }

    #[test]
    fn test_twin_rom_table_all_variable_reads() {
        let ctx = make_builder();
        let table_size = 5;
        let mut table_values: Vec<FieldPair<Bn254FrParams>> = Vec::new();

        for _ in 0..table_size {
            table_values.push([
                FieldT::from_witness(ctx.clone(), Fr::random_element()),
                FieldT::from_witness(ctx.clone(), Fr::random_element()),
            ]);
        }

        let expected: Vec<[Fr; 2]> = table_values
            .iter()
            .map(|p| [p[0].get_value(), p[1].get_value()])
            .collect();

        let mut table = TwinRomTable::new(table_values);

        for i in 0..table_size {
            let index = FieldT::from_witness(ctx.clone(), Fr::from(i as u64));
            let pair = table.read(&index);
            assert_eq!(pair[0].get_value(), expected[i][0]);
            assert_eq!(pair[1].get_value(), expected[i][1]);
        }

        let mut builder = ctx.borrow_mut();
        builder.finalize_circuit(false);
        assert!(
            UltraCircuitChecker::check(&mut builder).is_ok(),
            "Circuit check failed"
        );
    }
}
