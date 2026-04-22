//! Bit-level serialization of FN-DSA public keys, signatures, and private
//! keys, matching the Falcon round-3 reference (`codec.c`).
//!
//! Layout (round-3 §2.7):
//!
//! * **Public key**: 1 header byte `0x00 | logn`, then `n` coefficients
//!   packed as 14-bit big-endian fields (MSB-first within each byte).
//! * **Signature**: 1 header byte `0x30 | logn`, 40-byte nonce (salt),
//!   then compressed `s2` — per-coefficient 1 sign bit, 7 low bits, and
//!   a unary-encoded high tail. Zero-padded to a fixed length per logn.
//! * **Private key**: 1 header byte `0x50 | logn`, then packed `f`, `g`,
//!   `F` with widths determined by logn (6/6/8 for n=512, 5/5/8 for
//!   n=1024). `G` is recomputable from the NTRU equation and omitted.
//!
//! # Constant time
//! Control flow and buffer indexing depend only on `logn` and the fixed
//! output lengths — both public. The inner unary loop in signature
//! compression *does* iterate a data-dependent number of times; this
//! matches the round-3 reference, and the signature bytes themselves are
//! public output, not secret material.

#![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]

use crate::params::{
    NONCE_BYTES, N_1024, N_512, PRIV_1024_BYTES, PRIV_512_BYTES, PUB_1024_BYTES, PUB_512_BYTES, Q,
    SIG_1024_BYTES, SIG_512_BYTES,
};

const HDR_PUBKEY: u8 = 0x00;
const HDR_SIG: u8 = 0x30;
const HDR_SECKEY: u8 = 0x50;
const HDR_MASK: u8 = 0xF0;
const LOGN_MASK: u8 = 0x0F;

/// Coefficient bound for `s2` accepted by the compressed signature encoder
/// (round-3 reference: `|s2[i]| <= 2047`).
const S2_ABS_MAX: i32 = 2047;

/// Errors returned by the codec. Every variant reflects a public, structural
/// property of the input — never a secret.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CodecError {
    /// An input or output buffer has the wrong length, or `logn` is not
    /// one of the supported values (9 or 10).
    BadLength,
    /// The header byte does not match the expected type tag.
    BadHeader,
    /// A coefficient falls outside its declared range (e.g. `h[i] >= Q`,
    /// `|s2[i]| > 2047`, or an encoding of `-0`).
    OutOfRange,
    /// Non-zero bits appear in the trailing padding of a fixed-length
    /// container.
    BadPadding,
    /// The compressed stream demanded more bytes than were available.
    Truncated,
}

/// Parsed signature: header `logn`, 40-byte salt, and `s2` coefficients.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedSignature {
    /// Log-base-2 of the ring dimension recovered from the header.
    pub logn: u32,
    /// Random salt / nonce bound into the hash-to-point input.
    pub salt: [u8; NONCE_BYTES],
    /// Decompressed signature coefficients, `|s2[i]| <= 2047`.
    pub s2: Vec<i16>,
}

/// Parsed private key: `(f, g, F)`. `G` is not stored; it is recomputed from
/// the NTRU equation `fG - gF = q` when needed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedSecKey {
    /// Log-base-2 of the ring dimension recovered from the header.
    pub logn: u32,
    /// Small polynomial `f`, coefficients trimmed to `fg_bits(logn)` width.
    pub f: Vec<i8>,
    /// Small polynomial `g`, same width as `f`.
    pub g: Vec<i8>,
    /// Larger polynomial `F`, 8-bit signed coefficients.
    pub big_f: Vec<i8>,
}

/// Parsed public key: header `logn` and the NTT-domain coefficients `h`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedPubKey {
    /// Log-base-2 of the ring dimension recovered from the header.
    pub logn: u32,
    /// Public key coefficients in `[0, Q)`, natural order.
    pub h: Vec<u16>,
}

#[inline]
const fn check_logn(logn: u32) -> Result<(), CodecError> {
    if logn == 9 || logn == 10 {
        Ok(())
    } else {
        Err(CodecError::BadLength)
    }
}

#[inline]
const fn pubkey_bytes(logn: u32) -> usize {
    if logn == 9 {
        PUB_512_BYTES
    } else {
        PUB_1024_BYTES
    }
}

#[inline]
const fn sig_bytes(logn: u32) -> usize {
    if logn == 9 {
        SIG_512_BYTES
    } else {
        SIG_1024_BYTES
    }
}

#[inline]
const fn seckey_bytes(logn: u32) -> usize {
    if logn == 9 {
        PRIV_512_BYTES
    } else {
        PRIV_1024_BYTES
    }
}

#[inline]
const fn n_of(logn: u32) -> usize {
    if logn == 9 {
        N_512
    } else {
        N_1024
    }
}

/// Width (in bits) of each `f`/`g` coefficient in the encoded private key.
#[inline]
const fn fg_bits(logn: u32) -> u32 {
    if logn == 9 {
        6
    } else {
        5
    }
}

/// Width (in bits) of each `F` coefficient in the encoded private key.
#[inline]
const fn big_f_bits(_logn: u32) -> u32 {
    8
}

// ---- bit stream helpers -----------------------------------------------

struct BitWriter<'a> {
    buf: &'a mut [u8],
    pos: usize,
    acc: u64,
    acc_len: u32,
}

impl<'a> BitWriter<'a> {
    fn new(buf: &'a mut [u8]) -> Self {
        Self {
            buf,
            pos: 0,
            acc: 0,
            acc_len: 0,
        }
    }

    fn push(&mut self, val: u64, bits: u32) -> Result<(), CodecError> {
        debug_assert!(bits <= 32);
        let mask = (1u64 << bits) - 1;
        self.acc = (self.acc << bits) | (val & mask);
        self.acc_len += bits;
        while self.acc_len >= 8 {
            self.acc_len -= 8;
            if self.pos >= self.buf.len() {
                return Err(CodecError::BadLength);
            }
            self.buf[self.pos] = (self.acc >> self.acc_len) as u8;
            self.pos += 1;
        }
        Ok(())
    }

    /// Append a unary-coded value: `val` zero bits followed by a single `1`.
    fn push_unary(&mut self, val: u32) -> Result<(), CodecError> {
        // Emit `val + 1` bits of which only the LSB is 1.
        // For small `val` we can do this in one shot (bits ≤ 32).
        let bits = val + 1;
        debug_assert!(bits <= 32);
        self.push(1, bits)
    }

    fn finish(mut self) -> Result<usize, CodecError> {
        if self.acc_len > 0 {
            if self.pos >= self.buf.len() {
                return Err(CodecError::BadLength);
            }
            let shift = 8 - self.acc_len;
            self.buf[self.pos] = (self.acc << shift) as u8;
            self.pos += 1;
            self.acc_len = 0;
        }
        Ok(self.pos)
    }
}

struct BitReader<'a> {
    buf: &'a [u8],
    pos: usize,
    acc: u64,
    acc_len: u32,
}

impl<'a> BitReader<'a> {
    const fn new(buf: &'a [u8]) -> Self {
        Self {
            buf,
            pos: 0,
            acc: 0,
            acc_len: 0,
        }
    }

    fn pop(&mut self, bits: u32) -> Result<u64, CodecError> {
        debug_assert!(bits <= 32);
        while self.acc_len < bits {
            if self.pos >= self.buf.len() {
                return Err(CodecError::Truncated);
            }
            self.acc = (self.acc << 8) | u64::from(self.buf[self.pos]);
            self.pos += 1;
            self.acc_len += 8;
        }
        self.acc_len -= bits;
        let mask = (1u64 << bits) - 1;
        Ok((self.acc >> self.acc_len) & mask)
    }

    fn pop_bit(&mut self) -> Result<u8, CodecError> {
        Ok(self.pop(1)? as u8)
    }

    /// Verify that all remaining bits in `acc` and all remaining bytes in
    /// `buf` are zero. Used to enforce canonical padding.
    fn drain_zero(&self) -> Result<(), CodecError> {
        if self.acc_len > 0 && (self.acc & ((1u64 << self.acc_len) - 1)) != 0 {
            return Err(CodecError::BadPadding);
        }
        for &b in &self.buf[self.pos..] {
            if b != 0 {
                return Err(CodecError::BadPadding);
            }
        }
        Ok(())
    }
}

// ---- public key -------------------------------------------------------

/// Encode public key `h` (coefficients in `[0, Q)`) into `out`. Returns the
/// number of bytes written (always `pubkey_bytes(logn)`).
///
/// # Errors
/// * `BadLength` — `logn` unsupported, or `h.len() != 1 << logn`, or `out`
///   is not exactly `pubkey_bytes(logn)` bytes.
/// * `OutOfRange` — some `h[i] >= Q`.
pub fn encode_pubkey(logn: u32, h: &[u16], out: &mut [u8]) -> Result<usize, CodecError> {
    check_logn(logn)?;
    if h.len() != n_of(logn) || out.len() != pubkey_bytes(logn) {
        return Err(CodecError::BadLength);
    }
    for &x in h {
        if u32::from(x) >= Q {
            return Err(CodecError::OutOfRange);
        }
    }
    out[0] = HDR_PUBKEY | (logn as u8);
    let mut w = BitWriter::new(&mut out[1..]);
    for &x in h {
        w.push(u64::from(x), 14)?;
    }
    let body = w.finish()?;
    Ok(1 + body)
}

/// Decode a public key from `buf`. Returns the parsed `logn` and
/// coefficient vector.
///
/// # Errors
/// * `BadHeader` — header byte tag mismatch.
/// * `BadLength` — unsupported `logn` or wrong byte length.
/// * `OutOfRange` — decoded coefficient `>= Q`.
/// * `BadPadding` — trailing bits are non-zero.
pub fn decode_pubkey(buf: &[u8]) -> Result<DecodedPubKey, CodecError> {
    if buf.is_empty() {
        return Err(CodecError::BadLength);
    }
    if (buf[0] & HDR_MASK) != HDR_PUBKEY {
        return Err(CodecError::BadHeader);
    }
    let logn = u32::from(buf[0] & LOGN_MASK);
    check_logn(logn)?;
    if buf.len() != pubkey_bytes(logn) {
        return Err(CodecError::BadLength);
    }
    let n = n_of(logn);
    let mut h = vec![0u16; n];
    let mut r = BitReader::new(&buf[1..]);
    for slot in &mut h {
        let x = r.pop(14)? as u16;
        if u32::from(x) >= Q {
            return Err(CodecError::OutOfRange);
        }
        *slot = x;
    }
    r.drain_zero()?;
    Ok(DecodedPubKey { logn, h })
}

// ---- signature --------------------------------------------------------

/// Encode a signature: header + salt + compressed `s2` zero-padded to a
/// fixed length. Returns `sig_bytes(logn)` on success.
///
/// # Errors
/// * `BadLength` — wrong buffer sizes, or the compressed stream did not fit
///   in the allotted body length (caller should re-sign with a new salt).
/// * `OutOfRange` — some `|s2[i]| > 2047`.
pub fn encode_signature(
    logn: u32,
    salt: &[u8; NONCE_BYTES],
    s2: &[i16],
    out: &mut [u8],
) -> Result<usize, CodecError> {
    check_logn(logn)?;
    if s2.len() != n_of(logn) || out.len() != sig_bytes(logn) {
        return Err(CodecError::BadLength);
    }
    for &x in s2 {
        if !(-S2_ABS_MAX..=S2_ABS_MAX).contains(&i32::from(x)) {
            return Err(CodecError::OutOfRange);
        }
    }
    out[0] = HDR_SIG | (logn as u8);
    out[1..=NONCE_BYTES].copy_from_slice(salt);
    let body = &mut out[1 + NONCE_BYTES..];
    let body_len = body.len();
    let mut w = BitWriter::new(body);
    for &t in s2 {
        let (sign_bit, abs_val) = if t < 0 {
            (1u64, (-i32::from(t)) as u32)
        } else {
            (0u64, i32::from(t) as u32)
        };
        // 1 sign bit + 7 low bits, then unary for the remaining high bits.
        w.push(sign_bit, 1)?;
        w.push(u64::from(abs_val) & 0x7F, 7)?;
        w.push_unary(abs_val >> 7)?;
    }
    let written = w.finish()?;
    // Zero-pad any unused tail.
    for b in &mut body[written..body_len] {
        *b = 0;
    }
    Ok(1 + NONCE_BYTES + body_len)
}

/// Decode a signature. Returns `logn`, salt, and `s2` coefficients.
///
/// # Errors
/// * `BadHeader` — header byte tag mismatch.
/// * `BadLength` — unsupported `logn` or wrong byte length.
/// * `OutOfRange` — compressed `|s2|` exceeds 2047, or `-0` appears.
/// * `Truncated` — stream ended before decoding all `n` coefficients.
/// * `BadPadding` — trailing body bytes are non-zero.
pub fn decode_signature(buf: &[u8]) -> Result<DecodedSignature, CodecError> {
    if buf.is_empty() {
        return Err(CodecError::BadLength);
    }
    if (buf[0] & HDR_MASK) != HDR_SIG {
        return Err(CodecError::BadHeader);
    }
    let logn = u32::from(buf[0] & LOGN_MASK);
    check_logn(logn)?;
    if buf.len() != sig_bytes(logn) {
        return Err(CodecError::BadLength);
    }
    let n = n_of(logn);
    let mut salt = [0u8; NONCE_BYTES];
    salt.copy_from_slice(&buf[1..=NONCE_BYTES]);
    let body = &buf[1 + NONCE_BYTES..];
    let mut r = BitReader::new(body);
    let mut s2 = vec![0i16; n];
    for slot in &mut s2 {
        let sign = r.pop_bit()?;
        let mut mag = r.pop(7)? as u32;
        // Unary high bits: zeros until a terminating 1.
        loop {
            let bit = r.pop_bit()?;
            if bit == 1 {
                break;
            }
            mag += 128;
            if mag > S2_ABS_MAX as u32 {
                return Err(CodecError::OutOfRange);
            }
        }
        if sign == 1 && mag == 0 {
            return Err(CodecError::OutOfRange);
        }
        let mag_i = i32::try_from(mag).map_err(|_| CodecError::OutOfRange)?;
        *slot = if sign == 1 {
            -mag_i as i16
        } else {
            mag_i as i16
        };
    }
    r.drain_zero()?;
    Ok(DecodedSignature { logn, salt, s2 })
}

// ---- secret key -------------------------------------------------------

#[inline]
const fn fg_range(logn: u32) -> i32 {
    (1i32 << (fg_bits(logn) - 1)) - 1
}

#[inline]
const fn big_f_range(logn: u32) -> i32 {
    (1i32 << (big_f_bits(logn) - 1)) - 1
}

fn pack_signed(w: &mut BitWriter<'_>, src: &[i8], bits: u32, bound: i32) -> Result<(), CodecError> {
    for &x in src {
        let xi = i32::from(x);
        if xi < -bound || xi > bound {
            return Err(CodecError::OutOfRange);
        }
        let mask = (1u64 << bits) - 1;
        w.push(u64::from(x as u8) & mask, bits)?;
    }
    Ok(())
}

fn unpack_signed(
    r: &mut BitReader<'_>,
    dst: &mut [i8],
    bits: u32,
    bound: i32,
) -> Result<(), CodecError> {
    let sign_bit = 1u64 << (bits - 1);
    let forbidden = -(1i32 << (bits - 1));
    for slot in dst.iter_mut() {
        let raw = r.pop(bits)?;
        // Sign-extend from `bits` bits to i32.
        let value = if (raw & sign_bit) != 0 {
            (raw as i32) - (1i32 << bits)
        } else {
            raw as i32
        };
        if value == forbidden || value < -bound || value > bound {
            return Err(CodecError::OutOfRange);
        }
        *slot = value as i8;
    }
    Ok(())
}

/// Encode a private key as (header, `f`, `g`, `F`).
///
/// # Errors
/// * `BadLength` — unsupported `logn`, wrong slice or output lengths.
/// * `OutOfRange` — some coefficient falls outside the declared field width.
pub fn encode_seckey(
    logn: u32,
    f: &[i8],
    g: &[i8],
    big_f: &[i8],
    out: &mut [u8],
) -> Result<usize, CodecError> {
    check_logn(logn)?;
    let n = n_of(logn);
    if f.len() != n || g.len() != n || big_f.len() != n || out.len() != seckey_bytes(logn) {
        return Err(CodecError::BadLength);
    }
    out[0] = HDR_SECKEY | (logn as u8);
    let mut w = BitWriter::new(&mut out[1..]);
    pack_signed(&mut w, f, fg_bits(logn), fg_range(logn))?;
    pack_signed(&mut w, g, fg_bits(logn), fg_range(logn))?;
    pack_signed(&mut w, big_f, big_f_bits(logn), big_f_range(logn))?;
    let body = w.finish()?;
    Ok(1 + body)
}

/// Decode a private key.
///
/// # Errors
/// * `BadHeader`, `BadLength`, `OutOfRange`, `BadPadding` — see variants.
pub fn decode_seckey(buf: &[u8]) -> Result<DecodedSecKey, CodecError> {
    if buf.is_empty() {
        return Err(CodecError::BadLength);
    }
    if (buf[0] & HDR_MASK) != HDR_SECKEY {
        return Err(CodecError::BadHeader);
    }
    let logn = u32::from(buf[0] & LOGN_MASK);
    check_logn(logn)?;
    if buf.len() != seckey_bytes(logn) {
        return Err(CodecError::BadLength);
    }
    let n = n_of(logn);
    let mut f = vec![0i8; n];
    let mut g = vec![0i8; n];
    let mut big_f = vec![0i8; n];
    let mut r = BitReader::new(&buf[1..]);
    unpack_signed(&mut r, &mut f, fg_bits(logn), fg_range(logn))?;
    unpack_signed(&mut r, &mut g, fg_bits(logn), fg_range(logn))?;
    unpack_signed(&mut r, &mut big_f, big_f_bits(logn), big_f_range(logn))?;
    r.drain_zero()?;
    Ok(DecodedSecKey { logn, f, g, big_f })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand_u16_range(seed: u64, n: usize, hi: u16) -> Vec<u16> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = s;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                (z % u64::from(hi)) as u16
            })
            .collect()
    }

    fn rand_i8_range(seed: u64, n: usize, bound: i32) -> Vec<i8> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = s;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z ^= z >> 27;
                let span = (2 * bound + 1) as u64;
                ((z % span) as i32 - bound) as i8
            })
            .collect()
    }

    fn rand_i16_range(seed: u64, n: usize, bound: i32) -> Vec<i16> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = s;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z ^= z >> 27;
                let span = (2 * bound + 1) as u64;
                ((z % span) as i32 - bound) as i16
            })
            .collect()
    }

    #[test]
    fn pubkey_roundtrip_512() {
        let h = rand_u16_range(1, N_512, Q as u16);
        let mut buf = vec![0u8; PUB_512_BYTES];
        assert_eq!(encode_pubkey(9, &h, &mut buf).unwrap(), PUB_512_BYTES);
        assert_eq!(buf[0], 0x09);
        let dec = decode_pubkey(&buf).unwrap();
        assert_eq!(dec.logn, 9);
        assert_eq!(dec.h, h);
    }

    #[test]
    fn pubkey_roundtrip_1024() {
        let h = rand_u16_range(2, N_1024, Q as u16);
        let mut buf = vec![0u8; PUB_1024_BYTES];
        assert_eq!(encode_pubkey(10, &h, &mut buf).unwrap(), PUB_1024_BYTES);
        assert_eq!(buf[0], 0x0A);
        let dec = decode_pubkey(&buf).unwrap();
        assert_eq!(dec.logn, 10);
        assert_eq!(dec.h, h);
    }

    #[test]
    fn pubkey_rejects_coef_ge_q() {
        let mut h = vec![0u16; N_512];
        h[0] = Q as u16;
        let mut buf = vec![0u8; PUB_512_BYTES];
        assert_eq!(encode_pubkey(9, &h, &mut buf), Err(CodecError::OutOfRange));
    }

    #[test]
    fn pubkey_rejects_bad_header() {
        let mut buf = vec![0u8; PUB_512_BYTES];
        buf[0] = 0x39; // signature header, not pubkey
        assert_eq!(decode_pubkey(&buf), Err(CodecError::BadHeader));
    }

    #[test]
    fn pubkey_rejects_bad_length() {
        // Valid header byte, wrong total length.
        let mut buf = vec![0u8; PUB_512_BYTES - 1];
        buf[0] = 0x09;
        assert_eq!(decode_pubkey(&buf), Err(CodecError::BadLength));
        // Unsupported logn.
        let buf = vec![0u8; 1];
        assert_eq!(decode_pubkey(&buf), Err(CodecError::BadLength));
    }

    #[test]
    fn pubkey_rejects_out_of_range_coef() {
        // Craft a pubkey whose first 14-bit field encodes exactly Q.
        let mut buf = vec![0u8; PUB_512_BYTES];
        buf[0] = 0x09;
        let q = Q as u16;
        buf[1] = (q >> 6) as u8;
        buf[2] = ((q << 2) & 0xFC) as u8;
        assert_eq!(decode_pubkey(&buf), Err(CodecError::OutOfRange));
    }

    #[test]
    fn sig_roundtrip_512_small_coefs() {
        let salt = [0x42u8; NONCE_BYTES];
        let s2 = rand_i16_range(3, N_512, 100);
        let mut buf = vec![0u8; SIG_512_BYTES];
        assert_eq!(
            encode_signature(9, &salt, &s2, &mut buf).unwrap(),
            SIG_512_BYTES
        );
        assert_eq!(buf[0], 0x39);
        let dec = decode_signature(&buf).unwrap();
        assert_eq!(dec.logn, 9);
        assert_eq!(dec.salt, salt);
        assert_eq!(dec.s2, s2);
    }

    #[test]
    fn sig_roundtrip_extreme_coefs() {
        // Mix of magnitudes up to the allowed bound.
        let salt = [0u8; NONCE_BYTES];
        let mut s2 = vec![0i16; N_512];
        for (i, slot) in s2.iter_mut().enumerate() {
            *slot = match i % 4 {
                0 => 0,
                1 => 1,
                2 => -1,
                _ => 127,
            };
        }
        s2[0] = 2047;
        s2[1] = -2047;
        let mut buf = vec![0u8; SIG_512_BYTES];
        encode_signature(9, &salt, &s2, &mut buf).unwrap();
        let dec = decode_signature(&buf).unwrap();
        assert_eq!(dec.s2, s2);
    }

    #[test]
    fn sig_rejects_out_of_range() {
        let salt = [0u8; NONCE_BYTES];
        let mut s2 = vec![0i16; N_512];
        s2[0] = 2048;
        let mut buf = vec![0u8; SIG_512_BYTES];
        assert_eq!(
            encode_signature(9, &salt, &s2, &mut buf),
            Err(CodecError::OutOfRange)
        );
    }

    #[test]
    fn sig_rejects_bad_header() {
        let mut buf = vec![0u8; SIG_512_BYTES];
        buf[0] = 0x09;
        assert_eq!(decode_signature(&buf), Err(CodecError::BadHeader));
    }

    #[test]
    fn sig_rejects_minus_zero() {
        // Hand-crafted sig body: sign=1, low7=0, unary=1. That's 0b10000000_1
        // = the bit string "1 0000000 1" — 9 bits → first 8 bits = 0x80,
        // then the trailing 1 bit starts byte 2.
        let mut buf = vec![0u8; SIG_512_BYTES];
        buf[0] = 0x39;
        // salt stays zero.
        let body = &mut buf[1 + NONCE_BYTES..];
        body[0] = 0x80; // sign=1, low7=0000000
        body[1] = 0x80; // unary high: first bit=1 (stop) then 7 pad bits = 0
                        // The remaining 511 coefficients are zero — bits of "0 0000000 1"
                        // = 0x01 per coefficient... too fiddly to craft here. Instead rely
                        // on encoder to validate by building a full bad stream via encode
                        // path isn't possible; trust the first coef path.
                        // Truncate this test: we just check that an explicitly malformed
                        // body byte triggers the guard when parsing stops at the first
                        // coefficient — we need a fully formed sig for the full decode, so
                        // this test asserts the coefficient builder path rejects -0 by
                        // constructing a minimal valid sig and flipping one coef.
        let salt = [0u8; NONCE_BYTES];
        let mut s2 = vec![0i16; N_512];
        s2[0] = 1;
        let mut ok = vec![0u8; SIG_512_BYTES];
        encode_signature(9, &salt, &s2, &mut ok).unwrap();
        // Flip sign bit on the encoded first coefficient: its bit layout
        // starts at the first body byte. Byte pattern for `1`: sign=0,
        // low7=0000001, unary=1 → bits "0 0000001 1" = 0x01, 0x80.
        // Flipping sign gives "1 0000001 1" = 0x81, 0x80 → decodes as -1,
        // not -0, so not a good test. Instead craft explicit -0 stream:
        // bits "1 0000000 1" = 0x80, 0x80.
        let body = &mut ok[1 + NONCE_BYTES..];
        body[0] = 0x80;
        body[1] = 0x80;
        // Remaining bytes: each coef 0 encodes as "0 0000000 1" = 0x01.
        // Fill body with a valid pattern for all remaining coefs: each
        // coef adds 9 bits, so after coef 0 (bits 0..8 byte0 + bit0 byte1)
        // the cursor is at bit 1 of byte 1. Easier: expect Truncated or
        // OutOfRange; guaranteed we hit OutOfRange on coef 0.
        assert_eq!(decode_signature(&ok), Err(CodecError::OutOfRange));
    }

    #[test]
    fn seckey_roundtrip_512() {
        let fg_b = fg_range(9);
        let fb = big_f_range(9);
        let f = rand_i8_range(10, N_512, fg_b);
        let g = rand_i8_range(11, N_512, fg_b);
        let big_f = rand_i8_range(12, N_512, fb);
        let mut buf = vec![0u8; PRIV_512_BYTES];
        assert_eq!(
            encode_seckey(9, &f, &g, &big_f, &mut buf).unwrap(),
            PRIV_512_BYTES
        );
        assert_eq!(buf[0], 0x59);
        let dec = decode_seckey(&buf).unwrap();
        assert_eq!(dec.logn, 9);
        assert_eq!(dec.f, f);
        assert_eq!(dec.g, g);
        assert_eq!(dec.big_f, big_f);
    }

    #[test]
    fn seckey_roundtrip_1024() {
        let fg_b = fg_range(10);
        let fb = big_f_range(10);
        let f = rand_i8_range(20, N_1024, fg_b);
        let g = rand_i8_range(21, N_1024, fg_b);
        let big_f = rand_i8_range(22, N_1024, fb);
        let mut buf = vec![0u8; PRIV_1024_BYTES];
        encode_seckey(10, &f, &g, &big_f, &mut buf).unwrap();
        assert_eq!(buf[0], 0x5A);
        let dec = decode_seckey(&buf).unwrap();
        assert_eq!(dec.f, f);
        assert_eq!(dec.g, g);
        assert_eq!(dec.big_f, big_f);
    }

    #[test]
    fn seckey_rejects_out_of_range() {
        let mut f = vec![0i8; N_512];
        f[0] = 32; // fg_range(9) = 31
        let g = vec![0i8; N_512];
        let big_f = vec![0i8; N_512];
        let mut buf = vec![0u8; PRIV_512_BYTES];
        assert_eq!(
            encode_seckey(9, &f, &g, &big_f, &mut buf),
            Err(CodecError::OutOfRange)
        );
    }

    #[test]
    fn seckey_rejects_bad_header() {
        let mut buf = vec![0u8; PRIV_512_BYTES];
        buf[0] = 0x09;
        assert_eq!(decode_seckey(&buf), Err(CodecError::BadHeader));
    }
}
