//! Negacyclic Number-Theoretic Transform over `Z_q[x] / (x^n + 1)` with
//! `q = 12289` and `n ∈ {512, 1024}`.
//!
//! The transform evaluates a polynomial at the `n` primitive `2n`-th roots
//! of unity `{ψ^(2k+1) : k ∈ [0, n)}`, which allows a length-`n` point-wise
//! multiplication to implement negacyclic convolution (i.e. multiplication
//! modulo `x^n + 1`).
//!
//! # Algorithm
//!
//! Forward transform: in-place Cooley-Tukey (decimation-in-time) with
//! twiddles pre-stored in bit-reversed order, producing bit-reversed output.
//! Inverse transform: in-place Gentleman-Sande (decimation-in-frequency)
//! with inverse twiddles, consuming bit-reversed input and producing
//! natural-order output; a final multiplication by `n^{-1} mod q` closes
//! the inverse.
//!
//! A polynomial multiplication therefore reads: forward NTT on both
//! operands, point-wise multiply, inverse NTT — no explicit bit-reversal
//! step is needed because the pattern cancels.
//!
//! # Constant time
//!
//! All twiddle indexing is derived solely from the loop counters, which
//! depend on the public dimension `n` and not on polynomial coefficients.
//! The butterfly itself delegates to [`crate::field`], which is branch-free.
//! No early exits, no data-dependent memory access.

// All arithmetic below is bounded by `Q = 12289 < 2^14`, so every u32→u16
// truncation is provably safe. Range-based loops read more naturally than
// iterator chains when we need paired indices `(i, j)` within the
// butterflies.
#![allow(
    clippy::cast_possible_truncation,
    clippy::needless_range_loop,
    clippy::many_single_char_names
)]

mod tables;

use crate::field;
use crate::params::{N_1024, N_512, Q};

/// Forward negacyclic NTT for `n = 512`. Output is in bit-reversed order.
pub fn ntt_512(a: &mut [u16; N_512]) {
    ntt_generic(a, &tables::GM_512);
}

/// Forward negacyclic NTT for `n = 1024`. Output is in bit-reversed order.
pub fn ntt_1024(a: &mut [u16; N_1024]) {
    ntt_generic(a, &tables::GM_1024);
}

/// Inverse negacyclic NTT for `n = 512`. Input must be in bit-reversed order.
pub fn intt_512(a: &mut [u16; N_512]) {
    intt_generic(a, &tables::IGM_512, tables::N_INV_512);
}

/// Inverse negacyclic NTT for `n = 1024`. Input must be in bit-reversed order.
pub fn intt_1024(a: &mut [u16; N_1024]) {
    intt_generic(a, &tables::IGM_1024, tables::N_INV_1024);
}

/// Point-wise multiplication modulo `Q`, in place into `a`.
///
/// # Panics
/// Panics if `a.len() != b.len()`.
pub fn pointwise_mul(a: &mut [u16], b: &[u16]) {
    assert_eq!(a.len(), b.len(), "pointwise_mul: length mismatch");
    for (ai, bi) in a.iter_mut().zip(b.iter()) {
        *ai = mul_u16(*ai, *bi);
    }
}

/// Negacyclic polynomial multiplication `c = a · b mod (x^n + 1, Q)` for
/// `n = 512`. Consumes both operands, returns the product in natural order.
#[must_use]
pub fn negacyclic_mul_512(a: &[u16; N_512], b: &[u16; N_512]) -> [u16; N_512] {
    let mut ac = *a;
    let mut bc = *b;
    ntt_512(&mut ac);
    ntt_512(&mut bc);
    pointwise_mul(&mut ac, &bc);
    intt_512(&mut ac);
    ac
}

/// Negacyclic polynomial multiplication for `n = 1024`.
#[must_use]
pub fn negacyclic_mul_1024(a: &[u16; N_1024], b: &[u16; N_1024]) -> [u16; N_1024] {
    let mut ac = *a;
    let mut bc = *b;
    ntt_1024(&mut ac);
    ntt_1024(&mut bc);
    pointwise_mul(&mut ac, &bc);
    intt_1024(&mut ac);
    ac
}

// ---- internals ---------------------------------------------------------

#[inline]
fn mul_u16(a: u16, b: u16) -> u16 {
    let r = field::mul(u32::from(a), u32::from(b));
    debug_assert!(r < Q);
    r as u16
}

#[inline]
fn add_u16(a: u16, b: u16) -> u16 {
    field::add(u32::from(a), u32::from(b)) as u16
}

#[inline]
fn sub_u16(a: u16, b: u16) -> u16 {
    field::sub(u32::from(a), u32::from(b)) as u16
}

/// Cooley-Tukey forward NTT driver, parameterised by the twiddle table.
///
/// The table `gm[k]` contains `ψ^{brv_{log n}(k)}` in canonical form for
/// `k ∈ [1, n)`; `gm[0]` is unused.
fn ntt_generic(a: &mut [u16], gm: &[u16]) {
    let n = a.len();
    debug_assert_eq!(n, gm.len());
    debug_assert!(n.is_power_of_two());

    let mut t = n;
    let mut m = 1;
    while m < n {
        t >>= 1;
        let mut j1 = 0;
        for i in 0..m {
            let j2 = j1 + t;
            let s = gm[m + i];
            for j in j1..j2 {
                let u = a[j];
                let v = mul_u16(a[j + t], s);
                a[j] = add_u16(u, v);
                a[j + t] = sub_u16(u, v);
            }
            j1 += t << 1;
        }
        m <<= 1;
    }
}

/// Gentleman-Sande inverse NTT driver.
fn intt_generic(a: &mut [u16], igm: &[u16], n_inv: u16) {
    let n = a.len();
    debug_assert_eq!(n, igm.len());
    debug_assert!(n.is_power_of_two());

    let mut t = 1;
    let mut m = n;
    while m > 1 {
        let hm = m >> 1;
        let dt = t << 1;
        let mut j1 = 0;
        for i in 0..hm {
            let j2 = j1 + t;
            let s = igm[hm + i];
            for j in j1..j2 {
                let u = a[j];
                let v = a[j + t];
                a[j] = add_u16(u, v);
                a[j + t] = mul_u16(sub_u16(u, v), s);
            }
            j1 += dt;
        }
        t = dt;
        m = hm;
    }
    for x in a.iter_mut() {
        *x = mul_u16(*x, n_inv);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn schoolbook_512(a: &[u16; N_512], b: &[u16; N_512]) -> [u16; N_512] {
        let mut c = [0u16; N_512];
        let q64 = u64::from(Q);
        for i in 0..N_512 {
            for j in 0..N_512 {
                let prod = (u64::from(a[i]) * u64::from(b[j])) % q64;
                let k = (i + j) % N_512;
                let sign_neg = (i + j) >= N_512;
                let ck = u64::from(c[k]);
                let updated = if sign_neg {
                    (ck + q64 - prod) % q64
                } else {
                    (ck + prod) % q64
                };
                c[k] = u16::try_from(updated).expect("< q");
            }
        }
        c
    }

    fn small_rand(seed: u64, n: usize) -> Vec<u16> {
        // SplitMix64 — deterministic, dep-free. Good enough for tests.
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = s;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                u16::try_from(z % u64::from(Q)).expect("< q")
            })
            .collect()
    }

    #[test]
    fn roundtrip_512_zero() {
        let mut a = [0u16; N_512];
        ntt_512(&mut a);
        intt_512(&mut a);
        assert!(a.iter().all(|&x| x == 0));
    }

    #[test]
    fn roundtrip_1024_zero() {
        let mut a = [0u16; N_1024];
        ntt_1024(&mut a);
        intt_1024(&mut a);
        assert!(a.iter().all(|&x| x == 0));
    }

    #[test]
    fn roundtrip_512_random() {
        let src: Vec<u16> = small_rand(0xDEAD_BEEF, N_512);
        let mut a = [0u16; N_512];
        a.copy_from_slice(&src);
        ntt_512(&mut a);
        intt_512(&mut a);
        assert_eq!(a.as_slice(), src.as_slice());
    }

    #[test]
    fn roundtrip_1024_random() {
        let src: Vec<u16> = small_rand(0xCAFE_F00D, N_1024);
        let mut a = [0u16; N_1024];
        a.copy_from_slice(&src);
        ntt_1024(&mut a);
        intt_1024(&mut a);
        assert_eq!(a.as_slice(), src.as_slice());
    }

    #[test]
    fn negacyclic_mul_matches_schoolbook_512() {
        let av = small_rand(1, N_512);
        let bv = small_rand(2, N_512);
        let mut a = [0u16; N_512];
        let mut b = [0u16; N_512];
        a.copy_from_slice(&av);
        b.copy_from_slice(&bv);
        let fast = negacyclic_mul_512(&a, &b);
        let slow = schoolbook_512(&a, &b);
        assert_eq!(fast, slow);
    }

    #[test]
    fn linearity_512() {
        let av = small_rand(10, N_512);
        let bv = small_rand(20, N_512);
        let mut a = [0u16; N_512];
        let mut b = [0u16; N_512];
        let mut s = [0u16; N_512];
        for i in 0..N_512 {
            a[i] = av[i];
            b[i] = bv[i];
            s[i] = add_u16(a[i], b[i]);
        }
        ntt_512(&mut a);
        ntt_512(&mut b);
        ntt_512(&mut s);
        for i in 0..N_512 {
            assert_eq!(s[i], add_u16(a[i], b[i]), "linearity at {i}");
        }
    }

    #[test]
    fn constant_poly_evaluates_to_constant() {
        // a(x) = c (a constant) — NTT evaluates at n distinct roots, so all
        // outputs are c.
        let c: u16 = 7;
        let mut a = [0u16; N_512];
        a[0] = c;
        ntt_512(&mut a);
        assert!(a.iter().all(|&x| x == c), "all entries should equal {c}");
    }

    #[test]
    fn identity_poly_times_x_plus_1() {
        // (x + 1)*(x + 1) = x^2 + 2x + 1 mod (x^n + 1) for small n embedded.
        let mut a = [0u16; N_512];
        a[0] = 1;
        a[1] = 1;
        let prod = negacyclic_mul_512(&a, &a);
        assert_eq!(prod[0], 1);
        assert_eq!(prod[1], 2);
        assert_eq!(prod[2], 1);
        for i in 3..N_512 {
            assert_eq!(prod[i], 0, "coef {i}");
        }
    }
}
