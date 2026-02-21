// Group utilities (derive_generators, hash-to-curve)
//
// The actual implementation is in bbrs_crypto::generators.
// This module previously re-exported for convenience but now that
// ecc and crypto are separate crates, callers should use
// bbrs_crypto::generators::derive_generators directly.
