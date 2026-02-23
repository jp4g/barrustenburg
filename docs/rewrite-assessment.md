# Barretenberg C++ to Rust Rewrite Assessment

## Codebase Scope

| Metric | Value |
|--------|-------|
| **Total LOC** | ~411K lines |
| **Source files** | 2,049 (.cpp + .hpp) |
| **Top-level modules** | 45 |
| **Proving systems** | 5+ (Ultra Honk, Goblin, HyperNova, Chonk, ECCVM) |
| **Curves** | 4 (BN254, Grumpkin, secp256k1, secp256r1) |
| **Test files** | 384 |
| **Benchmark files** | 23 |
| **Fuzzer targets** | 25 |

## Difficulty: Very Hard (multi-year, multi-engineer effort)

This is one of the most complex cryptographic codebases in the ZK ecosystem. The difficulty stems not from any single component, but from the sheer breadth, depth, and tight coupling of the system.

---

## Major Work Packages

### 1. Finite Field Arithmetic (~20K LOC) — Hard ✅ DONE (Package 1)

- Montgomery-form field implementations for 4 curves
- x86-64 inline assembly for `MULX`/`ADCX`/`ADOX` instructions (BMI2/ADX)
- Extension fields: Fq2, Fq6, Fq12 for BN254 pairings
- WASM variant with 29-bit limb (9-limb) representation
- **Rust approach**: Use `ark-ff` or hand-roll with `core::arch::x86_64`. The inline asm is ~300+ lines of carefully tuned Montgomery multiplication/squaring. Rust `asm!` macro can do this but is nightly-only for some targets.
- **Testing**: Hard-coded field element test vectors exist. Export intermediate values from C++ field operations for cross-validation.

### 2. Elliptic Curve Operations (~20K LOC) — Hard ✅ DONE (Package 2)

- Group operations (affine, projective, Jacobian)
- Pippenger MSM (multi-scalar multiplication) with bucket processing
- Batched affine addition
- Precomputed generator tables per curve
- WNAF (Windowed Non-Adjacent Form) scalar decomposition
- **Rust approach**: `ark-ec` provides much of this, but BB's MSM is heavily customized. You'd likely need to port the Pippenger impl.
- **Testing**: `scalar_multiplication.test.cpp` (24K) has comprehensive MSM tests.

### 3. Polynomial Infrastructure (~6K LOC) — Medium ✅ DONE (Package 3)

- Univariate and multilinear polynomial representations
- File-backed coefficient storage via `mmap` for low-memory mode
- Polynomial arithmetic and evaluation
- **Rust approach**: Straightforward with `memmap2` for file-backed storage.

### 4. Commitment Schemes (~8K LOC) — Hard ✅ DONE (Package 4)

- KZG (pairing-based), IPA (inner product), Gemini (folding), SHPLONK (batched)
- Small subgroup IPA variant
- Claim batching infrastructure
- **Rust approach**: Arkworks has KZG/IPA, but BB's Gemini and SHPLONK are custom protocols. These need careful porting.

### 5. Proving Systems — Very Hard ✅ DONE (Package 5)

| System | LOC | Complexity |
|--------|-----|------------|
| Ultra Honk (prover/verifier) | ~6K | High |
| Sumcheck protocol | ~5K | High |
| Flavor system (28 variants) | ~5K | Very High — heavy C++ template metaprogramming |
| Relations (62 types) | ~11K | High |
| Goblin merge protocol | ~2K | Medium |
| HyperNova/Chonk IVC | ~4K | High |
| ECCVM | ~6K | Very High |
| Translator VM | ~5K | High |

The **flavor system** is the hardest part to port. It uses deep C++ template specialization to define 28 proving system variants from shared components. In Rust, this maps to a complex trait hierarchy with associated types — doable but architecturally challenging to get right.

### 6. Circuit Builders & Stdlib (~57K LOC) — Very Hard ✅ DONE (Package 6)

- Circuit-friendly implementations of every crypto primitive (in-circuit SHA256, Poseidon2, ECDSA, AES128, field operations, bigfield non-native arithmetic, biggroup, cycle_group)
- Recursive verifiers (honk, eccvm, goblin, translator, chonk)
- 742 Rust tests covering 561 C++ tests, all at COVERED or COVERED+
- See `pkg6_test_audit.md` for full component-by-component mapping.

### 7. DSL / ACIR Integration (~36K LOC) — Hard → **Package 9**

- Noir circuit format parsing and deserialization
- Bincode/msgpack serialization compatibility
- ACIR opcode to constraint translation
- **Rust approach**: Noir is already Rust-native, so this layer may actually *simplify* — you'd integrate directly with the Noir compiler rather than going through serialization.
- See Package 9 below for full breakdown.

### 8. VM2 / AVM (~137K LOC) — Very Hard → **Package 10**

- The Aztec Virtual Machine is 31% of the entire codebase (79K source + 48K tests + 9.5K fuzzer)
- Complete instruction set with constraining, trace generation, simulation
- Includes `avm_fuzzer` (9.4K) and `vm2_stub` (119 LOC)
- See Package 10 below for full breakdown.

### 9. Platform, Analysis & Serialization (~25K LOC) — Medium → **Package 13**

- Thread pool, aligned memory, logging, serialization (`common/`, `serialize/`, `env/`, `wasi/`)
- SMT verification, formal proofs, boomerang detection (analysis tooling)
- Benchmark infrastructure, Solidity helpers, SRS generation
- See Package 13 below for full breakdown.

### 10. CLI, Bindings & Aztec Infrastructure (~28K LOC) → **Packages 11–12**

- CLI (`bb` tool, 12K), API layer (2.2K), C FFI (3.1K), Node.js module (3.6K) → Package 11
- LMDB, world state, IPC, messaging, Aztec-specific public inputs → Package 12
- See Packages 11–12 below for full breakdown.

---

## Testing Against C++ Test Vectors

The codebase provides excellent cross-validation resources:

### Available test vectors

1. **Crypto primitives**: Hard-coded vectors for Blake2s, SHA256, Poseidon2, Pedersen (in `*.test.cpp` files)
2. **Serialization**: msgpack round-trip tests for field elements, verification keys, and VM instructions
3. **Binary test data**: 2MB+ of msgpack-encoded AVM test inputs
4. **Fixture files**: `fixtures.hpp` for Merkle tree tests

### Recommended strategy

**Phase 1: Export C++ test vectors to JSON/binary format**
- Write a small C++ harness that runs each crypto primitive on known inputs and dumps (input, output) pairs to files
- Cover: field arithmetic, curve operations, hash functions, commitments, MSM, polynomial evaluation

**Phase 2: Cross-language integration tests**
- Rust tests read C++ output files and verify identical results
- Use `serde` + msgpack (`rmp-serde`) since BB already uses msgpack
- Key checkpoints: field mul/add/inv, curve add/dbl/msm, hash(input)->output, commit(poly)->point

**Phase 3: Proof-level compatibility**
- Generate proofs in C++, verify in Rust (and vice versa)
- Requires identical transcript (Fiat-Shamir) behavior
- Requires identical serialization of proofs and VKs

**Phase 4: Differential fuzzing**
- Port the 25 existing C++ fuzzers to call both implementations
- LibFuzzer/cargo-fuzz with shared corpus

---

## Benchmarking Against C++

BB's existing benchmarks cover all performance-critical paths:
- Field operations (add, mul, sqr, inv)
- MSM (Pippenger)
- Poseidon2 hashing
- IPA/KZG commitment
- Circuit construction
- Full proving (Ultra Honk, Mega Honk)
- Merkle tree operations

### Recommended approach

1. Run C++ benchmarks via Google Benchmark (already set up):
   ```
   cd barretenberg/cpp/build && cmake --preset clang16
   make -j <module>_bench && ./<module>_bench --benchmark_format=json
   ```

2. Port benchmarks to Rust using `criterion`:
   - Match exact circuit sizes and parameters
   - Match threading configuration (same core count)

3. Critical comparison points:
   - Field mul throughput (ns/op) — this dominates everything
   - MSM for 2^16, 2^18, 2^20 points
   - Full prove time for equivalent circuits
   - Memory usage (peak RSS)
   - WASM performance (same benchmarks, wasm32-wasi)

### Performance expectations

A naive Rust port will be **10-30% slower** than C++ due to the hand-tuned x86 assembly in field arithmetic. To match C++, you'd need:
- Inline asm (`core::arch`) for Montgomery multiplication
- Equivalent Pippenger bucket processing
- Rayon-based parallelism tuned with similar cost heuristics

---

## Platform Equivalents

| Feature | Barretenberg C++ | Rust Equivalent |
|---------|------------------|-----------------|
| Assembly | Inline x64 asm (BMI2/ADX) | `core::arch::asm!` macro |
| SIMD | Infrastructure exists, disabled | `core::simd` (nightly) or `packed_simd` |
| Threading | Custom pool + heuristic scheduling | `rayon` + custom work-stealing |
| Aligned Memory | Platform-specific aligned_alloc | `#[repr(align(N))]` or `std::alloc` |
| Memory Profiling | Tracy integration | `dhat-rs` or `tracing` |
| WASM | First-class wasm32-wasi support | `wasm-bindgen` + `wasm32-wasi` |
| IPC/SHM | Custom lock-free queues | `parking_lot` + `memmap2` |
| Templates | Extensive C++ templates | Rust generics + traits |
| Concepts | C++20 concepts | Rust trait bounds |
| Constexpr | Heavy compile-time computation | `const fn` |
| Cache Alignment | 32/64-byte `alignas` | `#[repr(align(32))]` |

---

## Effort Estimate Summary

| Package | Component | C++ LOC | Status |
|---------|-----------|---------|--------|
| 1 | Field arithmetic | ~20K | DONE |
| 2 | Elliptic curve operations | ~20K | DONE |
| 3 | Polynomial infrastructure | ~6K | DONE |
| 4 | Commitment schemes | ~8K | DONE |
| 5 | Proving systems (Honk + flavors + relations + sumcheck) | ~37K | DONE |
| 6 | Circuit builder & stdlib | ~57K | DONE |
| 7 | Goblin pipeline (ECCVM, Translator VM, merge protocol) | ~16K | Not started |
| 8 | Client IVC (Chonk, HyperNova) | ~4K | Not started |
| 9 | DSL / ACIR backend | ~36K | Not started |
| 10 | VM2 / Aztec VM | ~137K | Not started |
| 11 | CLI & bindings (bb, bbapi, Node.js) | ~21K | Not started |
| 12 | Aztec state infrastructure | ~7K | Not started |
| 13 | Platform, analysis & serialization | ~25K | Not started |
| | Supporting (srs, transcript, crypto) | ~8K | DONE |
| **DONE** | | **~156K** | **1,408 tests, 0 failures** |
| **Remaining (without VM2)** | | **~109K** | |
| **Remaining (with VM2)** | | **~246K** | |

---

## Current Status (2026-02-23)

**Packages 1–6 are COMPLETE.** The Rust port covers field arithmetic, ECC, polynomials, commitment schemes, the full proving system (flavor, relations, sumcheck, Ultra Honk prover/verifier), circuit builder, and the complete stdlib (primitives, hash, encryption, recursive verifiers).

### What's built

| Package | Crate(s) | Rust Tests | Status |
|---------|----------|------------|--------|
| 1. Field Arithmetic | `bbrs-numeric` | 49 | DONE |
| 2. ECC (BN254, Grumpkin, secp256k1/r1) | `bbrs-ecc` | 295 | DONE |
| 3. Polynomial Infrastructure | `bbrs-polynomials` | 95 | DONE |
| 4. Commitment Schemes (KZG, IPA, Gemini, SHPLONK) | `bbrs-commitment-schemes` | 25 | DONE |
| 5. Proving Systems | `bbrs-flavor`, `bbrs-relations`, `bbrs-sumcheck`, `bbrs-honk`, `bbrs-ultra-honk` | 79 | DONE |
| 6. Circuit Builder & Stdlib | `bbrs-circuit-builder`, `bbrs-stdlib` | 770 | DONE |
| Supporting | `bbrs-srs`, `bbrs-transcript`, `bbrs-crypto` | 95 | DONE |
| **TOTAL** | **14 crates** | **1,408** | |

### Package 6 completion notes

Package 6 was delivered across 33 beads (B1–B33) in 6 waves. See `pkg6_test_audit.md` for the full component-by-component C++ → Rust test mapping. Key numbers:

- **561 C++ tests in scope → 742 Rust tests** (many components have MORE Rust tests than C++)
- All components at COVERED or COVERED+ — zero gaps
- bigfield: 115 Rust vs 114 C++, biggroup: 105 Rust vs 92 C++
- 36 E2E prove/verify tests validate the full UltraProver → UltraVerifier pipeline
- Recursive verifiers (honk, eccvm, goblin, translator, chonk) all implemented and tested

---

## Package 7: Goblin Pipeline (~16K LOC)

The Goblin proving system extends Ultra Honk with batched EC operations. It introduces a separate VM (ECCVM) that handles elliptic curve operations more efficiently than in-circuit constraints, connected via a Translator VM.

### Components

| Component | C++ Dir | C++ LOC | Tests | Description |
|-----------|---------|---------|-------|-------------|
| `trace_to_polynomials` | `trace_to_polynomials/` | 173 | 0 | Converts execution trace blocks → polynomial evaluations |
| `op_queue` | `op_queue/` | 1,650 | ~12 | ECC operation queue — batches EC ops for ECCVM processing |
| `multilinear_batching` | `multilinear_batching/` | 805 | ~3 | Batches multilinear polynomial claims for commitment opening |
| `eccvm` | `eccvm/` | 5,495 | ~15 | ECCVM circuit builder, flavor (1,086 LOC), prover, verifier, trace checker, MSM builder, precomputed tables |
| `translator_vm` | `translator_vm/` | 4,809 | ~10 | Translator VM — bridges ECCVM ↔ native field. Flavor (1,095 LOC), circuit builder, prover, verifier, proving key |
| `goblin` | `goblin/` | 2,031 | ~8 | Goblin orchestrator — merge prover/verifier, mock circuits, translation evaluations |
| `commitment_schemes_recursion` | `commitment_schemes_recursion/` | 1,046 | ~3 | Recursive (in-circuit) commitment scheme verification |

**Total**: ~16K LOC, ~50 tests

### Dependencies

```
trace_to_polynomials ← circuit_builder, polynomials
op_queue             ← ecc
multilinear_batching ← commitment_schemes
eccvm                ← trace_to_polynomials, op_queue, relations, sumcheck, flavor
translator_vm        ← eccvm, op_queue
goblin               ← eccvm, translator_vm
commitment_schemes_recursion ← stdlib, commitment_schemes
```

### New crates

- `bbrs-op-queue` — ECC op queue
- `bbrs-eccvm` — ECCVM flavor, builder, prover, verifier
- `bbrs-translator-vm` — Translator VM flavor, builder, prover, verifier
- `bbrs-goblin` — Goblin orchestrator + merge protocol

### Risk

| Risk | Severity | Notes |
|------|----------|-------|
| ECCVM flavor complexity | High | 1,086 LOC flavor definition — second most complex after Ultra |
| Translator VM field bridging | High | Non-trivial field conversion between BN254 and Grumpkin |
| Merge protocol correctness | Medium | Merge prover/verifier must match C++ transcript exactly |

---

## Package 8: Client IVC (~4K LOC)

Client IVC (Incremental Verifiable Computation) is the key Aztec feature — it enables recursive proof composition for private execution. "Chonk" is the accumulation scheme, HyperNova provides the folding protocol.

### Components

| Component | C++ Dir | C++ LOC | Tests | Description |
|-----------|---------|---------|-------|-------------|
| `chonk` | `chonk/` | 2,512 | ~8 | Client IVC accumulation — folds multiple proofs into one. Private execution steps, mock circuit producer, transcript invariants |
| `hypernova` | `hypernova/` | 1,393 | ~8 | HyperNova folding prover/verifier/decider |

**Total**: ~4K LOC, ~16 tests

### Dependencies

Blocked by Package 7 (Goblin pipeline). Chonk depends on Goblin, commitment_schemes_recursion, and stdlib. HyperNova depends on Chonk.

### New crates

- `bbrs-chonk` — Client IVC accumulation
- `bbrs-hypernova` — HyperNova folding protocol

---

## Package 9: DSL / ACIR Backend (~36K LOC)

The ACIR (Abstract Circuit Intermediate Representation) backend compiles Noir programs into circuits. This is the layer that makes barretenberg usable as a Noir proving backend.

### Components

| Component | C++ Dir | C++ LOC | Description |
|-----------|---------|---------|-------------|
| `dsl/acir_format` | `dsl/acir_format/` | ~30K | ACIR opcode → circuit constraint translation. Covers arithmetic, block constraints, EC operations, hash constraints, recursion constraints, bigint, AES, SHA, ECDSA, Schnorr, Poseidon2, multi-scalar-mul, range, logic, memory |
| `dsl/acir_proofs` | `dsl/acir_proofs/` | ~3K | Proof generation/verification from ACIR |
| `dsl/brillig_vm` | `dsl/brillig_vm/` | ~3K | Brillig (unconstrained code) VM execution |

**Total**: ~36K LOC

### Dependencies

Blocked by Package 7 (needs Goblin for MegaHonk circuits) and Package 6 (circuit builder, stdlib). However, a useful subset targeting UltraHonk-only circuits could start after Package 6.

### Notes

Noir is already Rust-native, so this layer may simplify significantly — direct integration with the Noir compiler rather than going through C++ serialization (bincode/msgpack). The Rust port could potentially expose a native Rust API to Noir, eliminating the FFI overhead entirely.

### New crates

- `bbrs-acir` — ACIR format parsing and circuit compilation
- `bbrs-brillig` — Brillig VM for unconstrained execution

---

## Package 10: VM2 / Aztec VM (~137K LOC)

The Aztec Virtual Machine (AVM) is the largest single component — 31% of the entire C++ codebase. It provides the execution environment for Aztec public functions with full constraining.

### Components

| Component | C++ Dir | C++ LOC | Description |
|-----------|---------|---------|-------------|
| `vm2` (source) | `vm2/` | ~79K | Complete instruction set implementation, trace generation, simulation, constraining, memory model, ALU, control flow, side effects |
| `vm2` (tests) | `vm2/` | ~48K | Extensive test suite covering all instructions and execution scenarios |
| `vm2_stub` | `vm2_stub/` | 119 | Stub interface for builds without full VM2 |
| `avm_fuzzer` | `avm_fuzzer/` | 9,373 | Fuzz testing harness for the AVM — instruction fuzzing, execution fuzzing, constraint satisfaction fuzzing |

**Total**: ~137K LOC (79K source + 48K tests + 9.5K fuzzer)

### Dependencies

Blocked by Package 7 (Goblin) and Package 9 (ACIR, for circuit compilation). The AVM generates circuits that are proved via MegaHonk/Goblin.

### Notes

This is by far the largest package. It may be out of scope for a standalone barretenberg port — the AVM is Aztec-specific and not needed for general-purpose Noir proving. Consider whether this package is needed for the target use case before committing to the port.

If ported, the fuzzer (`avm_fuzzer`) should be ported alongside or immediately after `vm2` to catch constraint bugs early.

### New crates

- `bbrs-vm2` — Aztec VM instruction set, trace generation, constraining
- `bbrs-avm-fuzzer` — AVM fuzz testing harness

---

## Package 11: CLI & Bindings (~21K LOC)

The `bb` CLI and API layers provide the user-facing interface to the proving system. These are the outermost layer — everything else must be complete first.

### Components

| Component | C++ Dir | C++ LOC | Description |
|-----------|---------|---------|-------------|
| `bb` | `bb/` | 12,237 | The `bb` command-line tool — prove, verify, write_vk, contract, gates, etc. |
| `api` | `api/` | 2,246 | Structured API layer — `api_ultra_honk`, `api_chonk`, `api_avm`, `api_msgpack`, `aztec_process` |
| `bbapi` | `bbapi/` | 3,104 | C FFI / shared library API for external consumers |
| `nodejs_module` | `nodejs_module/` | 3,636 | Node.js native module (N-API bindings) for `@aztec/bb.js` |

**Total**: ~21K LOC

### Dependencies

Blocked by all proving system packages (6-8). The CLI wraps the entire proving stack.

### Notes

The Rust port may not need all of these. A Rust `bb` binary is the natural replacement for the C++ CLI. The Node.js module could be replaced with `napi-rs` or `wasm-bindgen` targeting the same JS API. The C FFI layer (`bbapi`) is only needed if external C/C++ consumers must link against the Rust library.

### New crates

- `bbrs-cli` — `bb` CLI binary (clap-based)
- `bbrs-api` — Structured proving API
- `bbrs-ffi` — C FFI exports (optional)
- `bbrs-nodejs` — Node.js bindings via napi-rs (optional)

---

## Package 12: Aztec State Infrastructure (~7K LOC)

Infrastructure components specific to the Aztec network — persistent storage, inter-process communication, and protocol-specific data handling. These are needed for running a full Aztec node but not for standalone proving.

### Components

| Component | C++ Dir | C++ LOC | Description |
|-----------|---------|---------|-------------|
| `lmdblib` | `lmdblib/` | 3,807 | LMDB wrapper for persistent key-value storage (Merkle tree backing store) |
| `world_state` | `world_state/` | 3,330 | World state management — Merkle trees for notes, nullifiers, contracts, public data |
| `messaging` | `messaging/` | 188 | Message passing infrastructure |
| `ipc` | `ipc/` | 3,067 | Inter-process communication — shared memory, lock-free queues between bb processes |
| `special_public_inputs` | `special_public_inputs/` | 113 | Aztec-specific public input layout (fee, da gas, l2 gas, etc.) |
| `public_input_component` | `public_input_component/` | 88 | Public input routing for Aztec circuits |

**Total**: ~7K LOC (source only, minimal tests)

### Dependencies

Blocked by Package 6 (circuit builder) for world_state Merkle tree circuits. Otherwise mostly standalone infrastructure.

### Notes

This entire package is Aztec-specific. For a general-purpose barretenberg port, it can be deferred indefinitely. Only needed if the goal is a full Aztec node implementation in Rust.

### New crates

- `bbrs-world-state` — World state Merkle trees + LMDB storage
- `bbrs-ipc` — Inter-process communication (optional)

---

## Package 13: Platform, Analysis & Serialization (~25K LOC)

Cross-cutting platform support, analysis tooling, and serialization. Some of these are already partially addressed by Rust's standard library and ecosystem.

### Components

| Component | C++ Dir | C++ LOC | Description |
|-----------|---------|---------|-------------|
| `common` | `common/` | 7,771 | Shared utilities — thread pool, parallel_for, aligned allocation, logging, base64, memory management, container helpers, ref_array/ref_span/ref_vector, timer, assertions |
| `serialize` | `serialize/` | 531 | msgpack serialization (msgpack_impl, cbind, raw_pointer, schema) |
| `env` | `env/` | 202 | Environment detection (hardware capabilities, WASM detection) |
| `ext` | `ext/` | 1,047 | External library wrappers/shims |
| `wasi` | `wasi/` | 269 | WASM System Interface support |
| `benchmark` | `benchmark/` | 2,286 | Benchmark infrastructure and helpers |
| `grumpkin_srs_gen` | `grumpkin_srs_gen/` | 79 | Grumpkin SRS generation tool |
| `smt_verification` | `smt_verification/` | 6,150 | SMT-based circuit verification — uses Z3 solver to formally verify constraint correctness |
| `acir_formal_proofs` | `acir_formal_proofs/` | 1,910 | Formal proofs of ACIR constraint correctness |
| `boomerang_value_detection` | `boomerang_value_detection/` | 4,774 | Static analyzer for circuit variable reuse — detects "boomerang" values that could lead to soundness bugs via copy constraint analysis |
| `solidity_helpers` | `solidity_helpers/` | 537 | Solidity verifier contract generation helpers |

**Total**: ~25K LOC

### Subcategories

**Already addressed by Rust ecosystem (low priority):**
- `common` — `rayon` (parallel_for), `#[repr(align(N))]` (alignment), `log`/`tracing` (logging), standard collections. Some custom containers (ref_array, ref_span) may need thin Rust equivalents.
- `env` / `wasi` — `cfg` attributes, `wasm-bindgen`, `wasm32-wasi` target.
- `ext` — Rust has `bindgen` and native crate wrappers.
- `serialize` — `rmp-serde` (msgpack) or `serde` with any format.

**Analysis tooling (port for correctness assurance):**
- `smt_verification` — SMT-based constraint checking. Could interface with Z3 via `z3-sys` crate.
- `acir_formal_proofs` — Formal proofs of ACIR constraints. Depends on smt_verification.
- `boomerang_value_detection` — Static analysis of circuit copy constraints to detect soundness bugs. Depends on circuit_builder.

**Build/test tooling (port as needed):**
- `benchmark` — `criterion` crate replaces Google Benchmark.
- `grumpkin_srs_gen` — Small tool, port when needed for SRS generation.
- `solidity_helpers` — Needed if generating on-chain verifier contracts.

### New crates

- `bbrs-smt-verify` — SMT-based circuit verification (optional)
- `bbrs-solidity-gen` — Solidity verifier contract generation (optional)
