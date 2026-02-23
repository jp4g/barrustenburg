//! Circuit memory table primitives: ROM, RAM, and twin ROM tables.
//!
//! Port of `barretenberg/stdlib/primitives/memory/`.

pub mod ram_table;
pub mod rom_table;
pub mod twin_rom_table;

pub use ram_table::RamTable;
pub use rom_table::RomTable;
pub use twin_rom_table::TwinRomTable;
