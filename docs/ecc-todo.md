# ECC Module — Remaining Work

## Done

| Component | File(s) | Tests | Notes |
|-----------|---------|-------|-------|
| Field arithmetic | `fields/field.rs`, `field_params.rs` | 30+ | Montgomery form, u128 generic path |
| BN254 Fq/Fr | `curves/bn254.rs` | included | Small modulus path |
| Grumpkin | `curves/grumpkin.rs` | included | Field-swapped BN254 |
| secp256k1 Fq/Fr | `curves/secp256k1.rs` | included | Big modulus path |
| secp256r1 Fq/Fr | `curves/secp256r1.rs` | included | Big modulus, `has_a = true` |
| CurveParams trait | `groups/curve_params.rs` | — | 4 impls (G1 for each curve) |
| AffineElement | `groups/affine_element.rs` | included | Infinity encoding, on_curve |
| Element (Jacobian) | `groups/element.rs` | included | dbl, mixed add, full add |
| Scalar mul (basic) | `groups/element.rs` | included | Double-and-add, no endomorphism |

---

## Required API — Not Yet Implemented

These are needed by higher layers (commitments, prover) for correctness.

### derive_generators (hash-to-curve)
- **Stub:** `groups/group.rs`
- **C++ source:** `ecc/groups/group.hpp` lines 56-109
- **Blocked on:** `crypto_blake3s` crate integration
- **Used by:** Pedersen commitments, IPA, any protocol needing independent generators
- **Effort:** Medium — algorithm is straightforward, mostly BLAKE3 + field ops

### Extension fields (Fq2, Fq6, Fq12)
- **Stubs:** `fields/field2.rs`, `fields/field6.rs`, `fields/field12.rs`
- **C++ source:** `ecc/fields/field2.hpp`, `field6.hpp`, `field12.hpp`
- **Used by:** BN254 G2 operations, pairings
- **Effort:** Medium-large — standard tower construction but lots of operations
- **Note:** Only needed if pairings are in scope. Many Aztec circuits don't use pairings directly.

---

## Performance Optimizations — Deferred

These use the basic double-and-add scalar mul as a correct fallback. Only needed for prover performance.

### Pippenger MSM
- **Stub:** `scalar_multiplication.rs`
- **C++ source:** `ecc/scalar_multiplication/scalar_multiplication.hpp` + `process_buckets.hpp` + `bitvector.hpp`
- **Depends on:** WNAF, batched affine addition
- **Impact:** 10-100x speedup for multi-scalar multiplication (prover inner loop)
- **Effort:** Large — complex algorithm with threading, dynamic bucket sizing

### WNAF encoding
- **Stub:** `groups/wnaf.rs`
- **C++ source:** `ecc/groups/wnaf.hpp`
- **Used by:** Pippenger MSM
- **Effort:** Small-medium — bit manipulation, lookup tables

### Batched affine addition
- **Stub:** `batched_affine_addition.rs`
- **C++ source:** `ecc/batched_affine_addition/batched_affine_addition.hpp`
- **Used by:** Pippenger MSM (bucket accumulation)
- **Effort:** Medium — Montgomery batch inversion trick + threading

### GLV endomorphism
- **Stub:** None (flag exists: `CurveParams::USE_ENDOMORPHISM`)
- **C++ source:** Endomorphism params in field_params, splitting logic in scalar_mul
- **Used by:** Faster scalar multiplication on BN254/secp256k1 (curves with efficiently computable endomorphism)
- **Effort:** Medium — scalar decomposition + parallel double-and-add

### x86-64 ASM paths
- **C++ source:** `ecc/fields/asm_macros.hpp`, `ecc/groups/group_impl_asm.tcc`
- **Used by:** All field ops on x86-64
- **Impact:** ~2-3x speedup over generic u128 path
- **Effort:** Large — inline assembly via `core::arch`
- **Note:** Can defer indefinitely; Rust u128 path is correct

---

## Dependency Order for Remaining Work

```
crypto_blake3s ──► derive_generators
                        │
WNAF ──► Pippenger MSM ◄── Batched Affine Addition
              │
         GLV Endomorphism (optional speedup)

field2 ──► field6 ──► field12 ──► pairings (if needed)
```
