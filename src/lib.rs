//! FN-DSA (Falcon) pure-Rust implementation.
//!
//! See `CLAUDE.md` in the crate root for scope, KAT strategy, and
//! constant-time / formal-verification posture.

pub mod codec;
pub mod cose;
pub mod drbg;
pub mod fft;
pub mod field;
pub mod keygen;
pub mod ntt;
pub mod params;
pub mod sampler;
pub mod shake;
pub mod sign;
pub mod verify;
