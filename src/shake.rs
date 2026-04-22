//! SHAKE256 extendable-output function (FIPS 202).
//!
//! Zero-dependency pure-Rust implementation of Keccak-f[1600] + SHAKE256
//! sponge. Used by FN-DSA for:
//!
//! * deterministic PRNG expansion of the round-3 reference
//!   `inner_shake256_context` (seed → unbounded byte stream),
//! * hash-to-point during signing and verification,
//! * feeding entropy into the constant-time Gaussian sampler.
//!
//! The public type [`Shake256`] exposes `absorb` / `finalize` / `squeeze`.
//! Compatibility aliases [`Shake256::init`], [`Shake256::inject`],
//! [`Shake256::flip`], [`Shake256::extract`] mirror the reference C API so
//! that the rest of the port reads naturally against the Falcon round-3
//! source.
//!
//! # Constant-time
//!
//! The permutation and sponge are branch-free on the value bits of their
//! inputs. The only data-dependent control flow is on *lengths* — which are
//! public (a byte-count of produced output, not secret material).
//!
//! # TODO
//! * `Drop`-time zeroization of the state (see task #13 / memory hygiene).

// Pedantic lints tuned for this module: round constants are deliberately
// multi-digit; Keccak's natural expression uses short, single-letter
// coordinate names; and byte truncation on the u64-to-u8 read path is
// exactly what the spec demands.
#![allow(
    clippy::cast_possible_truncation,
    clippy::unreadable_literal,
    clippy::similar_names,
    clippy::needless_range_loop,
    clippy::many_single_char_names
)]

/// Number of Keccak-f[1600] rounds.
const KECCAK_ROUNDS: usize = 24;

/// SHAKE256 rate in bytes: (1600 − 2·256) / 8 = 136.
pub const SHAKE256_RATE_BYTES: usize = 136;

/// Iota round constants.
const RC: [u64; KECCAK_ROUNDS] = [
    0x0000000000000001,
    0x0000000000008082,
    0x800000000000808a,
    0x8000000080008000,
    0x000000000000808b,
    0x0000000080000001,
    0x8000000080008081,
    0x8000000000008009,
    0x000000000000008a,
    0x0000000000000088,
    0x0000000080008009,
    0x000000008000000a,
    0x000000008000808b,
    0x800000000000008b,
    0x8000000000008089,
    0x8000000000008003,
    0x8000000000008002,
    0x8000000000000080,
    0x000000000000800a,
    0x800000008000000a,
    0x8000000080008081,
    0x8000000000008080,
    0x0000000080000001,
    0x8000000080008008,
];

/// Rho rotation offsets, indexed `ROT[x][y]`.
const ROT: [[u32; 5]; 5] = [
    [0, 36, 3, 41, 18],
    [1, 44, 10, 45, 2],
    [62, 6, 43, 15, 61],
    [28, 55, 25, 21, 56],
    [27, 20, 39, 8, 14],
];

/// Keccak-f[1600] permutation over a 25-lane state.
fn keccak_f1600(s: &mut [u64; 25]) {
    for round in 0..KECCAK_ROUNDS {
        // θ
        let mut c = [0u64; 5];
        for x in 0..5 {
            c[x] = s[x] ^ s[x + 5] ^ s[x + 10] ^ s[x + 15] ^ s[x + 20];
        }
        let mut d = [0u64; 5];
        for x in 0..5 {
            d[x] = c[(x + 4) % 5] ^ c[(x + 1) % 5].rotate_left(1);
        }
        for y in 0..5 {
            for x in 0..5 {
                s[x + 5 * y] ^= d[x];
            }
        }

        // ρ + π, combined into a single lane permutation.
        let mut b = [0u64; 25];
        for x in 0..5 {
            for y in 0..5 {
                b[y + 5 * ((2 * x + 3 * y) % 5)] = s[x + 5 * y].rotate_left(ROT[x][y]);
            }
        }

        // χ
        for y in 0..5 {
            for x in 0..5 {
                s[x + 5 * y] = b[x + 5 * y] ^ ((!b[(x + 1) % 5 + 5 * y]) & b[(x + 2) % 5 + 5 * y]);
            }
        }

        // ι
        s[0] ^= RC[round];
    }
}

/// Xor `data` into the state starting at byte offset `byte_off` within the
/// rate block (lane `byte_off / 8`, little-endian within each lane).
fn xor_bytes_into_state(state: &mut [u64; 25], byte_off: usize, data: &[u8]) {
    for (i, &b) in data.iter().enumerate() {
        let off = byte_off + i;
        let lane = off / 8;
        let shift = (off % 8) * 8;
        state[lane] ^= u64::from(b) << shift;
    }
}

/// Extract `out.len()` bytes out of the state starting at byte offset
/// `byte_off`. `byte_off + out.len()` must not exceed [`SHAKE256_RATE_BYTES`].
fn read_bytes_from_state(state: &[u64; 25], byte_off: usize, out: &mut [u8]) {
    for (i, slot) in out.iter_mut().enumerate() {
        let off = byte_off + i;
        let lane = off / 8;
        let shift = (off % 8) * 8;
        *slot = (state[lane] >> shift) as u8;
    }
}

/// Sponge phase.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Phase {
    /// Bytes may still be absorbed into the state.
    Absorbing,
    /// The padding has been applied; only `squeeze` is legal.
    Squeezing,
}

/// SHAKE256 extendable-output function.
///
/// Instantiate with [`Shake256::new`], feed bytes via [`Shake256::absorb`],
/// transition with [`Shake256::finalize`], then pull an arbitrary number of
/// output bytes via one or more [`Shake256::squeeze`] calls.
#[derive(Clone, Debug)]
pub struct Shake256 {
    state: [u64; 25],
    /// Current byte position within the rate block.
    rate_pos: usize,
    phase: Phase,
}

impl Default for Shake256 {
    fn default() -> Self {
        Self::new()
    }
}

impl Shake256 {
    /// Create a fresh SHAKE256 context in the absorbing phase.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            state: [0; 25],
            rate_pos: 0,
            phase: Phase::Absorbing,
        }
    }

    /// Reset to the absorbing phase with a zero state.
    pub fn reset(&mut self) {
        self.state = [0; 25];
        self.rate_pos = 0;
        self.phase = Phase::Absorbing;
    }

    /// Absorb `data` into the sponge. May be called any number of times
    /// before [`finalize`](Self::finalize).
    ///
    /// # Panics
    /// Panics if called after [`finalize`](Self::finalize).
    pub fn absorb(&mut self, data: &[u8]) {
        assert!(self.phase == Phase::Absorbing, "absorb after finalize");
        let mut pos = 0;
        while pos < data.len() {
            let can = SHAKE256_RATE_BYTES - self.rate_pos;
            let take = can.min(data.len() - pos);
            xor_bytes_into_state(&mut self.state, self.rate_pos, &data[pos..pos + take]);
            self.rate_pos += take;
            pos += take;
            if self.rate_pos == SHAKE256_RATE_BYTES {
                keccak_f1600(&mut self.state);
                self.rate_pos = 0;
            }
        }
    }

    /// Apply SHAKE domain-separation padding (`0x1F … 0x80`) and transition
    /// to the squeezing phase.
    ///
    /// # Panics
    /// Panics if called more than once.
    pub fn finalize(&mut self) {
        assert!(self.phase == Phase::Absorbing, "double finalize");
        xor_bytes_into_state(&mut self.state, self.rate_pos, &[0x1F]);
        xor_bytes_into_state(&mut self.state, SHAKE256_RATE_BYTES - 1, &[0x80]);
        keccak_f1600(&mut self.state);
        self.rate_pos = 0;
        self.phase = Phase::Squeezing;
    }

    /// Squeeze `out.len()` bytes from the sponge. May be called any number
    /// of times after [`finalize`](Self::finalize).
    ///
    /// # Panics
    /// Panics if called before [`finalize`](Self::finalize).
    pub fn squeeze(&mut self, out: &mut [u8]) {
        assert!(self.phase == Phase::Squeezing, "squeeze before finalize");
        let mut pos = 0;
        while pos < out.len() {
            if self.rate_pos == SHAKE256_RATE_BYTES {
                keccak_f1600(&mut self.state);
                self.rate_pos = 0;
            }
            let can = SHAKE256_RATE_BYTES - self.rate_pos;
            let take = can.min(out.len() - pos);
            read_bytes_from_state(&self.state, self.rate_pos, &mut out[pos..pos + take]);
            self.rate_pos += take;
            pos += take;
        }
    }

    // --- Compatibility aliases matching Falcon round-3 reference naming. ---

    /// Alias for [`reset`](Self::reset) — matches `inner_shake256_init`.
    pub fn init(&mut self) {
        self.reset();
    }

    /// Alias for [`absorb`](Self::absorb) — matches `inner_shake256_inject`.
    pub fn inject(&mut self, data: &[u8]) {
        self.absorb(data);
    }

    /// Alias for [`finalize`](Self::finalize) — matches `inner_shake256_flip`.
    pub fn flip(&mut self) {
        self.finalize();
    }

    /// Alias for [`squeeze`](Self::squeeze) — matches `inner_shake256_extract`.
    pub fn extract(&mut self, out: &mut [u8]) {
        self.squeeze(out);
    }
}

/// Type alias mirroring the Falcon round-3 reference type name.
pub type InnerShake256Context = Shake256;

/// One-shot convenience: absorb `data`, finalize, and squeeze `out.len()`
/// bytes. Equivalent to `let mut s = Shake256::new(); s.absorb(data);
/// s.finalize(); s.squeeze(out);`.
pub fn shake256(data: &[u8], out: &mut [u8]) {
    let mut s = Shake256::new();
    s.absorb(data);
    s.finalize();
    s.squeeze(out);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reset_returns_to_absorb_phase() {
        let mut s = Shake256::new();
        s.absorb(b"hello");
        s.finalize();
        let mut buf = [0u8; 8];
        s.squeeze(&mut buf);
        s.reset();
        assert_eq!(s.phase, Phase::Absorbing);
        assert_eq!(s.rate_pos, 0);
        assert_eq!(s.state, [0u64; 25]);
    }

    #[test]
    fn incremental_absorb_equivalent_to_single() {
        let mut a = Shake256::new();
        a.absorb(b"hello ");
        a.absorb(b"world");
        a.finalize();
        let mut ao = [0u8; 32];
        a.squeeze(&mut ao);

        let mut b = Shake256::new();
        b.absorb(b"hello world");
        b.finalize();
        let mut bo = [0u8; 32];
        b.squeeze(&mut bo);

        assert_eq!(ao, bo);
    }

    #[test]
    fn incremental_squeeze_equivalent_to_single() {
        let mut a = Shake256::new();
        a.absorb(b"rfalco");
        a.finalize();
        let mut ao = [0u8; 300];
        a.squeeze(&mut ao);

        let mut b = Shake256::new();
        b.absorb(b"rfalco");
        b.finalize();
        let mut bo = [0u8; 300];
        // Non-aligned chunk sizes exercise the rate-boundary logic.
        let chunks = [7usize, 1, 128, 1, 135, 1, 27];
        let mut pos = 0;
        for &c in &chunks {
            b.squeeze(&mut bo[pos..pos + c]);
            pos += c;
        }
        assert_eq!(ao, bo);
    }

    #[test]
    #[should_panic(expected = "absorb after finalize")]
    fn absorb_after_finalize_panics() {
        let mut s = Shake256::new();
        s.finalize();
        s.absorb(b"x");
    }

    #[test]
    #[should_panic(expected = "squeeze before finalize")]
    fn squeeze_before_finalize_panics() {
        let mut s = Shake256::new();
        let mut out = [0u8; 1];
        s.squeeze(&mut out);
    }

    #[test]
    #[should_panic(expected = "double finalize")]
    fn double_finalize_panics() {
        let mut s = Shake256::new();
        s.finalize();
        s.finalize();
    }
}
