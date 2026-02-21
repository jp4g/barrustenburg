# Barretenberg C++ to Rust Porting Log

## Port #1: `crypto_keccak` (Layer 0 leaf)

**Date**: 2026-02-20
**C++ source**: `barretenberg/cpp/src/barretenberg/crypto/keccak/`
**Rust module**: `src/crypto/keccak/mod.rs`

### What was ported

The `crypto_keccak` module provides:
1. `ethash_keccak256(data, size)` — standard Ethereum Keccak-256 hash
2. `hash_field_elements(limbs, num_elements)` — serializes BN254 field elements (4 x u64 limbs) to big-endian bytes, then hashes with Keccak-256
3. `hash_field_element(limb)` — convenience wrapper for a single element

### Approach

- **Core primitive**: Delegated entirely to `tiny-keccak` crate (v2.0, `Keccak::v256()`). The BB C++ implementation is standard Ethereum Keccak-256 (0x01 padding, not SHA-3 0x06) with no algorithmic modifications.
- **Glue code**: Only the field element serialization (`hash_field_elements`, `hash_field_element`) was hand-written. This is ~15 lines of big-endian byte serialization matching C++ lines 114-134 of `keccak.cpp`.
- **`Keccak` wrapper class** (keccak.hpp lines 52-77): NOT ported yet. It depends on `bb::fr` (BN254 field element type) and `uint256_t`, which live in the `ecc` and `numeric` modules respectively. Will be added when those modules are ported.

### What was NOT ported

- `keccakf1600.cpp` — the raw Keccak-f[1600] permutation. Handled internally by `tiny-keccak`.
- `Keccak::hash(vector<uint256_t>)` — the higher-level wrapper that returns `bb::fr`. Blocked on `ecc`/`numeric` modules.

### Output type

Introduced `Keccak256` struct with `[u64; 4]` internal representation, mirroring C++ `struct keccak256 { uint64_t word64s[4]; }`.

### Tests

4 tests:
- `empty_hash_matches_known_keccak256` — validates against the well-known Keccak-256 empty-input digest
- `field_element_serialization_is_big_endian` — confirms limb serialization matches C++ byte extraction order
- `hash_multiple_field_elements` — multi-element hashing
- `rejects_invalid_limb_count` — panics on non-multiple-of-4 input

### Cross-validation status

- [x] Empty hash matches known Keccak-256 constant
- [x] Field element serialization order verified by construction
- [ ] TODO: Generate test vectors from C++ (`hash_field_elements` with known BN254 field elements) and validate in Rust

### Dependencies added

- `tiny-keccak = { version = "2.0", features = ["keccak"] }`

---

## Port #2: `crypto_aes128` (Layer 0 leaf) — crate only

**Date**: 2026-02-20
**Decision**: No custom code needed. Standard AES-128-CBC — use `aes` + `cbc` crates when a downstream module needs it. No glue code required.

---

## Port #3: `crypto_blake2s` (Layer 0 leaf) — crate only

**Date**: 2026-02-20
**Decision**: No custom code needed. Standard BLAKE2s-256 — use `blake2` crate when a downstream module needs it. No glue code required.

---

## Port #4: `crypto_blake3s_full` (Layer 0 leaf) — crate only

**Date**: 2026-02-20
**Decision**: No custom code needed. Standard BLAKE3 — use `blake3` with `default-features = false` (SIMD disabled for determinism) when a downstream module needs it. No glue code required.

---

## Port #5: `common` (Layer 1)

**Date**: 2026-02-20
**C++ source**: `barretenberg/cpp/src/barretenberg/common/` (~55 files)

### Assessment

`common` is a grab-bag utility module. Most of it maps to Rust built-ins or standard crates. Rather than porting it monolithically, we port only what downstream modules need as we go.

### Per-file porting decisions

#### Rust built-ins (no code needed)

| File | C++ purpose | Rust equivalent |
|------|-------------|-----------------|
| `assert.hpp` | `BB_ASSERT`, `BB_ASSERT_EQ`, debug-only variants | `assert!`, `assert_eq!`, `debug_assert!` |
| `net.hpp` | `is_little_endian()`, `htonl`/`ntohl` byte swaps | `u32::to_be()`, `u64::to_be_bytes()`, `cfg(target_endian)` |
| `throw_or_abort.hpp` | Throws exception or calls `std::abort` | `panic!` / `Result<T, E>` |
| `compiler_hints.hpp` | `BB_INLINE`, `BB_LIKELY`, `BB_UNLIKELY` | `#[inline(always)]`, `likely()`/`unlikely()` (nightly) |
| `constexpr_utils.hpp` | Compile-time `for` loop via templates | `const` generics, macros, or regular loops (LLVM unrolls) |
| `streams.hpp` | `operator<<` for vectors, arrays, tuples | `Debug`/`Display` trait impls |
| `container.hpp` | Deduplication helpers for sorted containers | `Vec::dedup()`, `BTreeSet` |
| `ref_array.hpp`, `ref_span.hpp`, `ref_vector.hpp` | Non-owning reference wrappers | `&[T]`, `&[&T]` — native Rust slices |
| `std_array.hpp`, `std_vector.hpp`, `std_string.hpp` | Stream operators for std types | `Debug`/`Display` |
| `zip_view.hpp` | Zip iterator over multiple containers | `Iterator::zip()` |
| `tuple.hpp`, `tuplet.hpp` | Tuple utilities | Native Rust tuples |
| `map.hpp` | `parallel_for` over map entries | `rayon::par_iter()` on `HashMap` |
| `timer.hpp` | Wall-clock timing | `std::time::Instant` |
| `printf.hpp` | Format string wrapper | `format!` macro |

#### Needs crate (add when needed)

| File | C++ purpose | Rust crate | When needed |
|------|-------------|------------|-------------|
| `serialize.hpp` | Big-endian binary serialization to/from buffers | `serde` + `bincode` or custom | `numeric` (Layer 1) |
| `log.hpp` / `log.cpp` / `debug_log.hpp` | Logging with levels, benchmark metrics | `tracing` or `log` + `env_logger` | Throughout, add early |
| `base64.hpp` / `base64.cpp` | Base64 encode/decode | `base64` crate | CLI/API layers |
| `msgpack_to_json.hpp` | Msgpack to JSON conversion | `rmp-serde` + `serde_json` | DSL/API layers |

#### TODO later (not needed until higher layers)

| File | C++ purpose | Notes |
|------|-------------|-------|
| `thread.hpp` / `thread.cpp` | Thread pool management | Use `rayon`. Needed at Layer 3+ (polynomials, MSM) |
| `parallel_for_*.cpp` (5 variants) | Parallel iteration with cost heuristics | Use `rayon::par_iter()`. Needed at Layer 3+ |
| `thread_pool.hpp` / `thread_pool.cpp` | Custom work-stealing thread pool | Use `rayon::ThreadPool`. Needed at Layer 3+ |
| `mem.hpp` | Aligned alloc (32/64-byte), Tracy malloc wrappers | `#[repr(align(N))]` + `std::alloc::Layout`. Needed at Layer 2+ (ecc) |
| `bbmalloc.hpp` / `bbmalloc.cpp` | Custom allocator with Tracy integration | Skip Tracy for now. Revisit for profiling |
| `tracy_mem/` | Tracy memory profiler integration | Skip entirely until profiling phase |
| `wasm_export.hpp` | `WASM_EXPORT`/`WASM_IMPORT` macros | `wasm-bindgen`. Needed when targeting WASM |
| `c_bind.hpp` / `c_bind.cpp` | C FFI bindings for WASM | `extern "C"` + `#[no_mangle]`. Needed for WASM |
| `flock.hpp` | File locking | `fs2` crate or `flock` syscall. Needed for SRS file access |
| `fuzzer.hpp` / `fuzzer_constants.hpp` | Fuzzer harness support | `cargo-fuzz`. Needed in testing phase |
| `bb_bench.hpp` / `bb_bench.cpp` / `google_bb_bench.hpp` | Benchmark harness | `criterion` crate. Needed in benchmarking phase |
| `test.hpp` / `benchmark.hpp` | Test/bench utilities | `#[cfg(test)]` + `criterion`. Needed in testing phase |
| `get_bytecode.hpp` / `get_bytecode.cpp` | Load bytecode from base64 | Needed at DSL layer |
| `utils.hpp` / `utils.cpp` | Hex conversion, tuple hashing | `hex` crate. Needed by `numeric` |
| `version.hpp` / `version.cpp` | Build version info | `env!("CARGO_PKG_VERSION")`. Needed at API layer |
| `named_union.hpp` | Tagged union helper | Rust `enum`. No port needed |
| `try_catch_shim.hpp` | Exception shim for no-exceptions builds | Not applicable in Rust |

### Impact on blockers

With this analysis, `common` does NOT need to be ported as a standalone module to unblock Layer 1:
- **`crypto_sha256`** needs only `assert.hpp` + `net.hpp` → Rust built-ins
- **`crypto_blake3s`** needs only `assert.hpp` + `constexpr_utils.hpp` → Rust built-ins
- **`numeric`** needs `assert.hpp` + `serialize.hpp` + `throw_or_abort.hpp` + `utils.hpp` → mostly built-ins, serialize is the only real work

`serialize.hpp` is the one piece that needs actual porting for `numeric` to work — it defines how uint128/uint256 are written to binary buffers.

---

## Port #6: `crypto_sha256` (Layer 1) — crate only

**Date**: 2026-02-20
**Decision**: Use `sha2` crate with `features = ["compress"]` (pin v0.10.x). The `compress256()` function maps directly to BB's `sha256_block()`. Only needs a thin wrapper to match BB's `fn(h_init, input) -> [u32; 8]` signature (the crate mutates in place). The `common` dependency is only `assert.hpp` + `net.hpp`, both of which are Rust built-ins.

---

## Port #7: `crypto_blake3s` (Layer 1) — crate only

**Date**: 2026-02-20
**Decision**: Same as `crypto_blake3s_full` — use `blake3` with `default-features = false`. The `common` dependency is only `assert.hpp` + `constexpr_utils.hpp`, both of which are Rust built-ins. Need to investigate what distinguishes this from `crypto_blake3s_full`.
