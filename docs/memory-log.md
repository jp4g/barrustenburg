# Barrustenburg Memory Log

General learnings and patterns discovered during the Barretenberg C++ to Rust port.

## Project structure

- Rust module tree mirrors BB's C++ directory layout: `src/crypto/keccak/` maps to `barretenberg/cpp/src/barretenberg/crypto/keccak/`
- Project is a library (`src/lib.rs`) with a binary target (`src/main.rs`)

## Leaf crypto modules — all use standard algorithms

All 5 Layer 0 crypto modules in BB implement standard algorithms from scratch (no external crypto libs). None have algorithmic modifications for ZK — the ZK-specific code lives in the `stdlib_*` circuit gadget modules higher up the dependency graph.

| Module | Standard | Rust crate replacement |
|--------|----------|----------------------|
| `crypto_aes128` | AES-128-CBC | `aes` + `cbc` |
| `crypto_blake2s` | BLAKE2s-256 | `blake2` |
| `crypto_blake3s_full` | BLAKE3 (SIMD disabled) | `blake3` with `default-features = false` |
| `crypto_keccak` | Ethereum Keccak-256 (0x01 padding) | `tiny-keccak` |
| `crypto_sha256` | SHA-256 (FIPS 180-4) | `sha2` with `features = ["compress"]` |

### Key nuances

- **Keccak uses 0x01 padding** (Ethereum variant), NOT SHA-3's 0x06. Use `Keccak::v256()`, not `Sha3`.
- **SHA-256 exposes `sha256_block()`** — the raw compression function — because ACIR's `Sha256Compression` opcode calls it directly. The `sha2` crate exposes this as `compress256()` behind the `compress` feature flag (pin to v0.10.x; v0.11 moves it to `block_api`).
- **BLAKE3 disables SIMD** in BB for deterministic cross-platform behavior. Use `default-features = false` in the Rust crate.
- **Keccak has field element glue** (`hash_field_elements`) that serializes `[u64; 4]` limbs to big-endian bytes before hashing. This is the only "custom" code; the core hash is standard.

## Porting patterns

- **Prefer crate delegation over hand-porting** for standard cryptographic primitives. Only write glue code for BB-specific API surfaces (field element serialization, exposed compression functions).
- **Test against known constants first** (e.g. empty-input hash digests), then cross-validate against C++ test vectors.
- **Higher-level wrappers** (like `Keccak::hash` returning `bb::fr`) depend on types from other modules (`ecc`, `numeric`). Port the wrapper when those modules arrive; don't introduce placeholder types.

## Dependency graph insights

- `common` is the most foundational module (8+ dependents) — port it early
- `stdlib_circuit_builders` is the most depended-on module (15+ dependents) — it's the critical path bottleneck
- `dsl` and `chonk` have a circular dependency — will need a shared types crate or feature flags
- The AVM (`vm2`) is 31% of the codebase (~127K LOC) and can be skipped if not needed
