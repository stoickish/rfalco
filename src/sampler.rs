//! Constant-time integer Gaussian sampler (Karney / cumulative-table rejection).
//!
//! Implementation pending — tracked as task #9. Consequence: sign-KAT bit-exact
//! match to the reference FP sampler is not achievable; KAT validation uses
//! `verify(pk, msg, kat_sig) == accept` instead. See `CLAUDE.md`.
