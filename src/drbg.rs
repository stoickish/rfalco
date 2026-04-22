//! NIST SP 800-90A AES-256 `CTR_DRBG`, no derivation function.
//!
//! Reproduces the byte-exact output of the NIST PQC KAT harness's
//! `randombytes()` (`rng.c` in the round-3 submission package). Seed
//! with a 48-byte entropy input (and optional 48-byte personalization
//! string) and call [`CtrDrbg::fill`] to draw bytes deterministically.
//!
//! # Scope
//! This DRBG exists to reproduce NIST PQC KAT files and drive
//! deterministic tests. It is **not** a production RNG: the AES
//! implementation here is a plain byte-oriented S-box version and is
//! therefore not constant-time against cache-timing adversaries. In
//! production, callers should supply their own CSPRNG. Seeded test
//! harnesses and KAT reproduction are the only intended callers.
//!
//! # Algorithm
//! Per SP 800-90A §10.2.1:
//! - `CTR_DRBG_Update(provided[48], Key, V)`: generate three AES-256(V+1..V+3)
//!   blocks, XOR with `provided`, split into new `(Key, V)`.
//! - `Instantiate(entropy[48], personalization)`: start with `(Key, V) = 0`,
//!   `Update(entropy XOR personalization)`.
//! - `Generate(out, len)`: counter-mode AES-256(V+i) into `out`, then
//!   `Update(zero[48])` to reseed internal state.
//!
//! `reseed_counter` and `additional_input` are tracked but not
//! exposed; the PQC harness does not use either.

#![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]

/// AES S-box (FIPS 197 Figure 7).
const SBOX: [u8; 256] = [
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16,
];

/// Round constants for AES-256 key expansion. Index `i` is used when
/// `word_index / Nk == i + 1` (Nk = 8 for AES-256).
const RCON: [u8; 7] = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40];

/// AES block size in bytes.
const BLOCK: usize = 16;
/// AES-256 round count.
const NR: usize = 14;
/// Number of round keys (one per round plus the initial whitening).
const NUM_RK: usize = NR + 1;

/// GF(2^8) multiplication by x: shift left with conditional XOR of 0x1B.
#[inline]
const fn xtime(x: u8) -> u8 {
    let hi = x >> 7;
    (x << 1) ^ (hi * 0x1b)
}

/// AES-256 key schedule: expand a 32-byte key into 15 × 16-byte round keys.
fn key_expand_256(key: &[u8; 32]) -> [[u8; BLOCK]; NUM_RK] {
    let mut w = [[0u8; 4]; 4 * NUM_RK];
    for i in 0..8 {
        w[i].copy_from_slice(&key[4 * i..4 * i + 4]);
    }
    for i in 8..4 * NUM_RK {
        let mut temp = w[i - 1];
        if i % 8 == 0 {
            let t0 = temp[0];
            temp[0] = SBOX[temp[1] as usize] ^ RCON[i / 8 - 1];
            temp[1] = SBOX[temp[2] as usize];
            temp[2] = SBOX[temp[3] as usize];
            temp[3] = SBOX[t0 as usize];
        } else if i % 8 == 4 {
            for b in &mut temp {
                *b = SBOX[*b as usize];
            }
        }
        for j in 0..4 {
            w[i][j] = w[i - 8][j] ^ temp[j];
        }
    }
    let mut rk = [[0u8; BLOCK]; NUM_RK];
    for r in 0..NUM_RK {
        for c in 0..4 {
            rk[r][4 * c..4 * c + 4].copy_from_slice(&w[4 * r + c]);
        }
    }
    rk
}

/// AES `SubBytes`: apply the S-box byte-wise.
#[inline]
fn sub_bytes(s: &mut [u8; BLOCK]) {
    for b in s.iter_mut() {
        *b = SBOX[*b as usize];
    }
}

/// AES `ShiftRows`: row `r` cyclically shifts left by `r`. State is laid
/// out column-major, so byte index `r + 4*c` holds row `r`, column `c`.
#[inline]
fn shift_rows(s: &mut [u8; BLOCK]) {
    let mut t = [0u8; BLOCK];
    for c in 0..4 {
        for r in 0..4 {
            t[r + 4 * c] = s[r + 4 * ((c + r) % 4)];
        }
    }
    *s = t;
}

/// AES `MixColumns`: mix each column by the fixed MDS matrix in GF(2^8).
#[inline]
fn mix_columns(s: &mut [u8; BLOCK]) {
    for c in 0..4 {
        let a0 = s[4 * c];
        let a1 = s[4 * c + 1];
        let a2 = s[4 * c + 2];
        let a3 = s[4 * c + 3];
        let t = a0 ^ a1 ^ a2 ^ a3;
        s[4 * c] ^= t ^ xtime(a0 ^ a1);
        s[4 * c + 1] ^= t ^ xtime(a1 ^ a2);
        s[4 * c + 2] ^= t ^ xtime(a2 ^ a3);
        s[4 * c + 3] ^= t ^ xtime(a3 ^ a0);
    }
}

/// `AddRoundKey`: XOR the round key into the state.
#[inline]
fn add_round_key(s: &mut [u8; BLOCK], rk: &[u8; BLOCK]) {
    for i in 0..BLOCK {
        s[i] ^= rk[i];
    }
}

/// AES-256 single-block encrypt.
fn aes256_encrypt(pt: &[u8; BLOCK], rk: &[[u8; BLOCK]; NUM_RK]) -> [u8; BLOCK] {
    let mut s = *pt;
    add_round_key(&mut s, &rk[0]);
    for round_key in rk.iter().take(NR).skip(1) {
        sub_bytes(&mut s);
        shift_rows(&mut s);
        mix_columns(&mut s);
        add_round_key(&mut s, round_key);
    }
    sub_bytes(&mut s);
    shift_rows(&mut s);
    add_round_key(&mut s, &rk[NR]);
    s
}

/// Big-endian in-place increment of a 16-byte counter. Wraps at `2^128`;
/// SP 800-90A permits wrap because the `reseed_counter` bounds the
/// number of draws before reseeding.
#[inline]
fn inc_v(v: &mut [u8; BLOCK]) {
    for byte in v.iter_mut().rev() {
        *byte = byte.wrapping_add(1);
        if *byte != 0 {
            return;
        }
    }
}

/// `CTR_DRBG_Update(provided, Key, V)` per SP 800-90A §10.2.1.2.
fn update(provided: &[u8; 48], key: &mut [u8; 32], v: &mut [u8; BLOCK]) {
    let rk = key_expand_256(key);
    let mut tmp = [0u8; 48];
    for blk in 0..3 {
        inc_v(v);
        let ct = aes256_encrypt(v, &rk);
        tmp[blk * BLOCK..(blk + 1) * BLOCK].copy_from_slice(&ct);
    }
    for i in 0..48 {
        tmp[i] ^= provided[i];
    }
    key.copy_from_slice(&tmp[..32]);
    v.copy_from_slice(&tmp[32..48]);
}

/// Seeded AES-256 `CTR_DRBG` matching the NIST PQC KAT harness.
///
/// See the module doc for the threat-model caveat: not for production.
pub struct CtrDrbg {
    /// 256-bit AES key of the DRBG state.
    key: [u8; 32],
    /// 128-bit counter of the DRBG state.
    v: [u8; BLOCK],
    /// Number of `Generate` calls since instantiation (capped by
    /// SP 800-90A §10.2.1.5 at `2^48`; we track but do not enforce).
    reseed_counter: u64,
}

impl CtrDrbg {
    /// Instantiate a DRBG with a 48-byte entropy input and an optional
    /// personalization string. If `personalization` is shorter than 48
    /// bytes it is zero-padded; longer inputs are rejected by
    /// `debug_assert!`.
    #[must_use]
    pub fn new(entropy: &[u8; 48], personalization: Option<&[u8]>) -> Self {
        let mut seed = *entropy;
        if let Some(p) = personalization {
            debug_assert!(p.len() <= 48, "personalization string too long");
            for (s, b) in seed.iter_mut().zip(p.iter()) {
                *s ^= *b;
            }
        }
        let mut key = [0u8; 32];
        let mut v = [0u8; BLOCK];
        update(&seed, &mut key, &mut v);
        Self {
            key,
            v,
            reseed_counter: 1,
        }
    }

    /// Fill `out` with pseudorandom bytes. Matches the reference
    /// `randombytes(out, outlen)` byte-for-byte.
    pub fn fill(&mut self, out: &mut [u8]) {
        let rk = key_expand_256(&self.key);
        let mut i = 0;
        while i < out.len() {
            inc_v(&mut self.v);
            let ct = aes256_encrypt(&self.v, &rk);
            let take = core::cmp::min(BLOCK, out.len() - i);
            out[i..i + take].copy_from_slice(&ct[..take]);
            i += take;
        }
        let zeros = [0u8; 48];
        update(&zeros, &mut self.key, &mut self.v);
        self.reseed_counter += 1;
    }

    /// Current reseed counter value (test-only accessor).
    #[cfg(test)]
    pub(crate) const fn reseed_counter(&self) -> u64 {
        self.reseed_counter
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// FIPS 197 Appendix C.3 — AES-256 known answer.
    #[test]
    fn aes256_fips197_appendix_c3() {
        let mut key = [0u8; 32];
        for (i, b) in key.iter_mut().enumerate() {
            *b = i as u8;
        }
        let rk = key_expand_256(&key);
        let pt: [u8; 16] = [
            0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd,
            0xee, 0xff,
        ];
        let ct = aes256_encrypt(&pt, &rk);
        let expected: [u8; 16] = [
            0x8e, 0xa2, 0xb7, 0xca, 0x51, 0x67, 0x45, 0xbf, 0xea, 0xfc, 0x49, 0x90, 0x4b, 0x49,
            0x60, 0x89,
        ];
        assert_eq!(ct, expected);
    }

    /// Big-endian counter increment wraps through the boundary.
    #[test]
    fn counter_increment_big_endian() {
        let mut v = [0u8; BLOCK];
        inc_v(&mut v);
        assert_eq!(v[15], 1);
        assert_eq!(v[..15], [0u8; 15]);

        let mut v = [0u8; BLOCK];
        v[15] = 0xff;
        inc_v(&mut v);
        assert_eq!(v[14], 1);
        assert_eq!(v[15], 0);
    }

    /// NIST PQC KAT harness reproduction. Entropy input = `0..48`,
    /// no personalization: the first 48-byte `randombytes` draw must
    /// match the first `seed = ...` line of the round-3 KAT files
    /// (value shared by every NIST PQC submission using `rng.c`).
    #[test]
    fn pqc_kat_harness_first_draw() {
        let mut entropy = [0u8; 48];
        for (i, b) in entropy.iter_mut().enumerate() {
            *b = i as u8;
        }
        let mut drbg = CtrDrbg::new(&entropy, None);
        let mut seed = [0u8; 48];
        drbg.fill(&mut seed);
        let expected: [u8; 48] = hex48(
            "061550234D158C5EC95595FE04EF7A25767F2E24CC2BC479D09D86DC9ABCFDE7056A8C266F9EF97ED08541DBD2E1FFA1",
        );
        assert_eq!(seed, expected);
        assert_eq!(drbg.reseed_counter(), 2);
    }

    /// Second draw continues the same state: confirms the internal
    /// reseed `Update(0^48)` step after each `fill`.
    #[test]
    fn pqc_kat_harness_second_draw() {
        let mut entropy = [0u8; 48];
        for (i, b) in entropy.iter_mut().enumerate() {
            *b = i as u8;
        }
        let mut drbg = CtrDrbg::new(&entropy, None);
        let mut first = [0u8; 48];
        drbg.fill(&mut first);
        let mut second = [0u8; 48];
        drbg.fill(&mut second);
        let expected: [u8; 48] = hex48(
            "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC81ADDE6AEEB4A5A875C3BFCADFA958F",
        );
        assert_eq!(second, expected);
        assert_eq!(drbg.reseed_counter(), 3);
    }

    /// Non-block-aligned draws equal the prefix of a single block-aligned draw.
    #[test]
    fn partial_block_draw_is_prefix_of_full_block() {
        let mut entropy = [0u8; 48];
        for (i, b) in entropy.iter_mut().enumerate() {
            *b = i as u8;
        }
        let mut a = CtrDrbg::new(&entropy, None);
        let mut b = CtrDrbg::new(&entropy, None);
        let mut big = [0u8; 48];
        a.fill(&mut big);
        let mut small = [0u8; 17];
        b.fill(&mut small);
        assert_eq!(small[..], big[..17]);
    }

    /// Hex decoder for 48-byte fixture strings. Panics on bad input.
    fn hex48(s: &str) -> [u8; 48] {
        assert_eq!(s.len(), 96);
        let mut out = [0u8; 48];
        for (i, b) in out.iter_mut().enumerate() {
            let hi = hex_digit(s.as_bytes()[2 * i]);
            let lo = hex_digit(s.as_bytes()[2 * i + 1]);
            *b = (hi << 4) | lo;
        }
        out
    }

    fn hex_digit(c: u8) -> u8 {
        match c {
            b'0'..=b'9' => c - b'0',
            b'a'..=b'f' => c - b'a' + 10,
            b'A'..=b'F' => c - b'A' + 10,
            _ => panic!("bad hex digit: {c}"),
        }
    }
}
