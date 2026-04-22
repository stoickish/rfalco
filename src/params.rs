//! FN-DSA parameter sets (round-3 / expected FIPS 206).

/// NTT modulus. `q - 1 = 12288 = 2^12 · 3` supports NTT for `n` up to 1024.
pub const Q: u32 = 12289;

/// Bit width of the modulus.
pub const Q_BITS: u32 = 14;

/// Ring dimension for FN-DSA-512.
pub const N_512: usize = 512;

/// Ring dimension for FN-DSA-1024.
pub const N_1024: usize = 1024;

/// FN-DSA-512 serialized signature length (draft-ietf-cose-falcon-04 Table 1).
pub const SIG_512_BYTES: usize = 666;

/// FN-DSA-1024 serialized signature length.
pub const SIG_1024_BYTES: usize = 1280;

/// Fixed nonce / salt length in bytes (Falcon round-3 §3.8).
pub const NONCE_BYTES: usize = 40;

/// FN-DSA-512 serialized public-key length.
pub const PUB_512_BYTES: usize = 897;

/// FN-DSA-512 serialized private-key length.
pub const PRIV_512_BYTES: usize = 1281;

/// FN-DSA-1024 serialized public-key length.
pub const PUB_1024_BYTES: usize = 1793;

/// FN-DSA-1024 serialized private-key length.
pub const PRIV_1024_BYTES: usize = 2305;

/// Round-3 signature norm bound (squared) for FN-DSA-512.
pub const SIG_BOUND_SQ_512: u64 = 34_034_726;

/// Round-3 signature norm bound (squared) for FN-DSA-1024.
pub const SIG_BOUND_SQ_1024: u64 = 70_265_242;
