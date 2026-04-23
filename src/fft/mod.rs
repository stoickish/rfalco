//! FP tower-of-rings FFT for the Falcon sign path.
//!
//! Ports the round-3 reference FFT and polynomial-tower operations
//! (`PQClean` `falcon-512/clean/fft.c`, MIT) with matching op order so
//! signatures can be bit-exact against round-3 `PQCgenKAT_sign.rsp`.
//!
//! # Representation
//! Polynomials live in `R = R[x] / (x^n + 1)` with `n = 2^logn`. After
//! [`fft`], an `n`-element `f64` slice holds `n/2` complex evaluations
//! `f(w_j)` at the primitive `2n`-th roots of unity `w_j = w^(2j+1)`.
//! Only half of the roots are stored because the other half are
//! complex conjugates of the first.
//!
//! Layout, for a slice `f` of length `n`:
//! - real part of the `u`-th stored evaluation: `f[u]`
//! - imaginary part: `f[u + n/2]`
//!
//! The indexing is bit-reversed per the round-3 reference; the same
//! permutation is used by [`poly_split`] / [`poly_merge`], so callers
//! that only compose these public functions never see the reversal.
//!
//! # Arithmetic posture
//! IEEE 754 binary64 round-to-nearest-even, native `f64` ops.
//! - No `fast-math`.
//! - No FMA contraction — `a * b + c * d` is two multiplies and one
//!   add, never a fused `mul_add`.
//! - No reassociation of additions.
//!
//! This matches the reference's "no fused ops" posture. Bit-exact
//! cross-platform behavior is expected on `x86_64/ARM64` with default
//! rounding (no x87 80-bit reg use, no flush-to-zero).
//!
//! # Constant time
//! FP operations on secret data follow isochronous sequences — no
//! branches on secret magnitudes, no data-dependent memory access.
//! This matches the reference's CT stance; microarchitectural FP
//! timing quirks are acknowledged and out of scope.

// `suboptimal_flops` would rewrite `a*b + c*d` to `mul_add`, which is
// FMA contraction — forbidden by the Falcon round-3 reference FFT for
// bit-exact KAT reproduction. Disabled crate-wide for this module.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::module_name_repetitions,
    clippy::similar_names,
    clippy::suboptimal_flops,
    clippy::needless_range_loop,
    clippy::many_single_char_names,
    clippy::too_many_lines,
    clippy::too_long_first_doc_paragraph,
    clippy::float_cmp
)]

mod tables;

use tables::{GM_TAB, P2_TAB};

/// Reinterpret the GM table entry at index `k` as a pair `(re, im) = w^rev(k)`.
#[inline]
fn gm(k: usize) -> (f64, f64) {
    (
        f64::from_bits(GM_TAB[2 * k]),
        f64::from_bits(GM_TAB[2 * k + 1]),
    )
}

/// Complex add: `(dr, di) = (ar, ai) + (br, bi)`.
#[inline]
fn fpc_add(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    (ar + br, ai + bi)
}

/// Complex sub: `(dr, di) = (ar, ai) - (br, bi)`.
#[inline]
fn fpc_sub(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    (ar - br, ai - bi)
}

/// Complex mul: `(dr, di) = (ar, ai) * (br, bi)`. Matches the
/// reference `FPC_MUL` macro op order (no FMA).
#[inline]
fn fpc_mul(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    let dr = ar * br - ai * bi;
    let di = ar * bi + ai * br;
    (dr, di)
}

/// Complex div: `(dr, di) = (ar, ai) / (br, bi)`. Matches the reference
/// `FPC_DIV` macro op order: first invert the denominator via
/// `1 / (br² + bi²)`, then multiply by the numerator with the sign of
/// the imaginary part flipped.
#[inline]
fn fpc_div(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    let m = 1.0 / (br * br + bi * bi);
    let br = br * m;
    let bi = -bi * m;
    let dr = ar * br - ai * bi;
    let di = ar * bi + ai * br;
    (dr, di)
}

/// Forward FFT of a real polynomial in `R[x]/(x^n+1)`, `n = 2^logn`,
/// `1 <= logn <= 10`. Input `f` holds `n` real coefficients; output
/// holds the Falcon FFT representation (see module doc).
///
/// # Panics
/// Panics (via `debug_assert`) if `f.len() != 1 << logn` or if `logn`
/// is out of range.
pub fn fft(f: &mut [f64], logn: u32) {
    debug_assert!((1..=10).contains(&logn));
    let n = 1usize << logn;
    debug_assert_eq!(f.len(), n);
    let hn = n >> 1;
    let mut t = hn;
    let mut m: usize = 2;
    let mut u: u32 = 1;
    while u < logn {
        let ht = t >> 1;
        let hm = m >> 1;
        let mut j1: usize = 0;
        for i1 in 0..hm {
            let j2 = j1 + ht;
            let (s_re, s_im) = gm(m + i1);
            for j in j1..j2 {
                let x_re = f[j];
                let x_im = f[j + hn];
                let y_re0 = f[j + ht];
                let y_im0 = f[j + ht + hn];
                let (y_re, y_im) = fpc_mul(y_re0, y_im0, s_re, s_im);
                let (sum_re, sum_im) = fpc_add(x_re, x_im, y_re, y_im);
                f[j] = sum_re;
                f[j + hn] = sum_im;
                let (dif_re, dif_im) = fpc_sub(x_re, x_im, y_re, y_im);
                f[j + ht] = dif_re;
                f[j + ht + hn] = dif_im;
            }
            j1 += t;
        }
        t = ht;
        m <<= 1;
        u += 1;
    }
}

/// Inverse FFT. Undoes [`fft`] up to `f64` rounding.
///
/// # Panics
/// Panics (via `debug_assert`) if `f.len() != 1 << logn` or if `logn`
/// is out of range.
pub fn ifft(f: &mut [f64], logn: u32) {
    debug_assert!((1..=10).contains(&logn));
    let n = 1usize << logn;
    debug_assert_eq!(f.len(), n);
    let hn = n >> 1;
    let mut t: usize = 1;
    let mut m: usize = n;
    let mut u = logn;
    while u > 1 {
        let hm = m >> 1;
        let dt = t << 1;
        let mut j1: usize = 0;
        let mut i1: usize = 0;
        while j1 < hn {
            let j2 = j1 + t;
            let (s_re0, s_im0) = gm(hm + i1);
            let s_re = s_re0;
            let s_im = -s_im0;
            for j in j1..j2 {
                let x_re = f[j];
                let x_im = f[j + hn];
                let y_re = f[j + t];
                let y_im = f[j + t + hn];
                let (sum_re, sum_im) = fpc_add(x_re, x_im, y_re, y_im);
                f[j] = sum_re;
                f[j + hn] = sum_im;
                let (dif_re, dif_im) = fpc_sub(x_re, x_im, y_re, y_im);
                let (res_re, res_im) = fpc_mul(dif_re, dif_im, s_re, s_im);
                f[j + t] = res_re;
                f[j + t + hn] = res_im;
            }
            i1 += 1;
            j1 += dt;
        }
        t = dt;
        m = hm;
        u -= 1;
    }
    // Final division by N/2 folded with the skipped last butterfly
    // layer — see reference comment.
    if logn > 0 {
        let ni = f64::from_bits(P2_TAB[logn as usize]);
        for slot in f.iter_mut() {
            *slot *= ni;
        }
    }
}

/// `a += b` in coefficient or FFT domain.
pub fn poly_add(a: &mut [f64], b: &[f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    debug_assert_eq!(b.len(), n);
    for u in 0..n {
        a[u] += b[u];
    }
}

/// `a -= b`.
pub fn poly_sub(a: &mut [f64], b: &[f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    debug_assert_eq!(b.len(), n);
    for u in 0..n {
        a[u] -= b[u];
    }
}

/// `a = -a`.
pub fn poly_neg(a: &mut [f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    for slot in a.iter_mut().take(n) {
        *slot = -*slot;
    }
}

/// Complex conjugate of an FFT-domain polynomial: negate the
/// imaginary halves of every stored evaluation.
pub fn poly_adj(a: &mut [f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    for u in (n >> 1)..n {
        a[u] = -a[u];
    }
}

/// Pointwise complex multiply in FFT domain: `a *= b`.
pub fn poly_mul(a: &mut [f64], b: &[f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    debug_assert_eq!(b.len(), n);
    let hn = n >> 1;
    for u in 0..hn {
        let (a_re, a_im) = (a[u], a[u + hn]);
        let (b_re, b_im) = (b[u], b[u + hn]);
        let (dr, di) = fpc_mul(a_re, a_im, b_re, b_im);
        a[u] = dr;
        a[u + hn] = di;
    }
}

/// Pointwise `a *= adj(b)` in FFT domain.
pub fn poly_muladj(a: &mut [f64], b: &[f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    debug_assert_eq!(b.len(), n);
    let hn = n >> 1;
    for u in 0..hn {
        let (a_re, a_im) = (a[u], a[u + hn]);
        let b_re = b[u];
        let b_im = -b[u + hn];
        let (dr, di) = fpc_mul(a_re, a_im, b_re, b_im);
        a[u] = dr;
        a[u + hn] = di;
    }
}

/// Pointwise `a *= adj(a)` — result is real-valued; the imaginary
/// half is zeroed.
pub fn poly_mulselfadj(a: &mut [f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    let hn = n >> 1;
    for u in 0..hn {
        let (a_re, a_im) = (a[u], a[u + hn]);
        a[u] = a_re * a_re + a_im * a_im;
        a[u + hn] = 0.0;
    }
}

/// Scale every coefficient by `x`.
pub fn poly_mulconst(a: &mut [f64], x: f64, logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    for slot in a.iter_mut().take(n) {
        *slot *= x;
    }
}

/// Pointwise complex divide in FFT domain: `a /= b`.
pub fn poly_div(a: &mut [f64], b: &[f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    debug_assert_eq!(b.len(), n);
    let hn = n >> 1;
    for u in 0..hn {
        let (a_re, a_im) = (a[u], a[u + hn]);
        let (b_re, b_im) = (b[u], b[u + hn]);
        let (dr, di) = fpc_div(a_re, a_im, b_re, b_im);
        a[u] = dr;
        a[u + hn] = di;
    }
}

/// Write `d[u] = 1 / (|a(w_u)|² + |b(w_u)|²)` — real-valued output.
pub fn poly_invnorm2(d: &mut [f64], a: &[f64], b: &[f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(d.len(), n);
    debug_assert_eq!(a.len(), n);
    debug_assert_eq!(b.len(), n);
    let hn = n >> 1;
    for u in 0..hn {
        let (a_re, a_im) = (a[u], a[u + hn]);
        let (b_re, b_im) = (b[u], b[u + hn]);
        d[u] = 1.0 / ((a_re * a_re + a_im * a_im) + (b_re * b_re + b_im * b_im));
    }
}

/// `d = F * adj(f) + G * adj(g)` pointwise in FFT domain.
pub fn poly_add_muladj(
    d: &mut [f64],
    big_f: &[f64],
    big_g: &[f64],
    f: &[f64],
    g: &[f64],
    logn: u32,
) {
    let n = 1usize << logn;
    debug_assert_eq!(d.len(), n);
    let hn = n >> 1;
    for u in 0..hn {
        let big_f_re = big_f[u];
        let big_f_im = big_f[u + hn];
        let big_g_re = big_g[u];
        let big_g_im = big_g[u + hn];
        let f_re = f[u];
        let f_im = f[u + hn];
        let g_re = g[u];
        let g_im = g[u + hn];
        let (a_re, a_im) = fpc_mul(big_f_re, big_f_im, f_re, -f_im);
        let (b_re, b_im) = fpc_mul(big_g_re, big_g_im, g_re, -g_im);
        d[u] = a_re + b_re;
        d[u + hn] = a_im + b_im;
    }
}

/// `a *= b` where `b` is self-adjoint (real-valued, imaginary halves
/// unused). Both real and imaginary halves of `a` are scaled by
/// `b[u]`.
pub fn poly_mul_autoadj(a: &mut [f64], b: &[f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    let hn = n >> 1;
    for u in 0..hn {
        let bu = b[u];
        a[u] *= bu;
        a[u + hn] *= bu;
    }
}

/// `a /= b` where `b` is self-adjoint (real-valued).
pub fn poly_div_autoadj(a: &mut [f64], b: &[f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(a.len(), n);
    let hn = n >> 1;
    for u in 0..hn {
        let ib = 1.0 / b[u];
        a[u] *= ib;
        a[u + hn] *= ib;
    }
}

/// LDL decomposition of a self-adjoint `2x2` matrix
/// `[[g00, g01], [adj(g01), g11]]`. Overwrites `g01` with the
/// off-diagonal `L` factor `L_10 = g01 / g00`, and updates `g11` to
/// the Schur complement `g11 - L_10 * adj(g01_original)`. `g00` is
/// left unchanged (it already is `D_0`).
pub fn poly_ldl(g00: &[f64], g01: &mut [f64], g11: &mut [f64], logn: u32) {
    let n = 1usize << logn;
    debug_assert_eq!(g00.len(), n);
    debug_assert_eq!(g01.len(), n);
    debug_assert_eq!(g11.len(), n);
    let hn = n >> 1;
    for u in 0..hn {
        let g00_re = g00[u];
        let g00_im = g00[u + hn];
        let g01_re = g01[u];
        let g01_im = g01[u + hn];
        let g11_re = g11[u];
        let g11_im = g11[u + hn];
        let (mu_re, mu_im) = fpc_div(g01_re, g01_im, g00_re, g00_im);
        let (p_re, p_im) = fpc_mul(mu_re, mu_im, g01_re, -g01_im);
        let (r_re, r_im) = fpc_sub(g11_re, g11_im, p_re, p_im);
        g11[u] = r_re;
        g11[u + hn] = r_im;
        g01[u] = mu_re;
        g01[u + hn] = -mu_im;
    }
}

/// LDL with separate output buffers: leaves `g00`, `g01`, `g11`
/// untouched and writes the updated diagonal to `d11` and the
/// `L_10` factor to `l10`.
pub fn poly_ldl_mv(
    d11: &mut [f64],
    l10: &mut [f64],
    g00: &[f64],
    g01: &[f64],
    g11: &[f64],
    logn: u32,
) {
    let n = 1usize << logn;
    let hn = n >> 1;
    for u in 0..hn {
        let g00_re = g00[u];
        let g00_im = g00[u + hn];
        let g01_re = g01[u];
        let g01_im = g01[u + hn];
        let g11_re = g11[u];
        let g11_im = g11[u + hn];
        let (mu_re, mu_im) = fpc_div(g01_re, g01_im, g00_re, g00_im);
        let (p_re, p_im) = fpc_mul(mu_re, mu_im, g01_re, -g01_im);
        let (r_re, r_im) = fpc_sub(g11_re, g11_im, p_re, p_im);
        d11[u] = r_re;
        d11[u + hn] = r_im;
        l10[u] = mu_re;
        l10[u + hn] = -mu_im;
    }
}

/// Split `f ∈ R_n` in FFT representation into `(f0, f1) ∈ R_{n/2}^2`
/// such that `f(x) = f0(x²) + x * f1(x²)`. `logn >= 1`.
pub fn poly_split(f0: &mut [f64], f1: &mut [f64], f: &[f64], logn: u32) {
    debug_assert!(logn >= 1);
    let n = 1usize << logn;
    let hn = n >> 1;
    let qn = hn >> 1;
    debug_assert_eq!(f.len(), n);
    debug_assert_eq!(f0.len(), hn);
    debug_assert_eq!(f1.len(), hn);

    // logn=1 special case: only one stored evaluation, the inner
    // loop below runs zero times.
    f0[0] = f[0];
    f1[0] = f[hn];

    for u in 0..qn {
        let a_re = f[2 * u];
        let a_im = f[2 * u + hn];
        let b_re = f[2 * u + 1];
        let b_im = f[2 * u + 1 + hn];
        let (t_re, t_im) = fpc_add(a_re, a_im, b_re, b_im);
        f0[u] = 0.5 * t_re;
        f0[u + qn] = 0.5 * t_im;
        let (d_re, d_im) = fpc_sub(a_re, a_im, b_re, b_im);
        let (s_re, s_im) = gm(u + hn);
        let (m_re, m_im) = fpc_mul(d_re, d_im, s_re, -s_im);
        f1[u] = 0.5 * m_re;
        f1[u + qn] = 0.5 * m_im;
    }
}

/// Inverse of [`poly_split`]: merge `(f0, f1) ∈ R_{n/2}^2` back into
/// `f ∈ R_n`.
pub fn poly_merge(f: &mut [f64], f0: &[f64], f1: &[f64], logn: u32) {
    debug_assert!(logn >= 1);
    let n = 1usize << logn;
    let hn = n >> 1;
    let qn = hn >> 1;
    debug_assert_eq!(f.len(), n);
    debug_assert_eq!(f0.len(), hn);
    debug_assert_eq!(f1.len(), hn);

    f[0] = f0[0];
    f[hn] = f1[0];

    for u in 0..qn {
        let a_re = f0[u];
        let a_im = f0[u + qn];
        let (s_re, s_im) = gm(u + hn);
        let (b_re, b_im) = fpc_mul(f1[u], f1[u + qn], s_re, s_im);
        let (t_re, t_im) = fpc_add(a_re, a_im, b_re, b_im);
        f[2 * u] = t_re;
        f[2 * u + hn] = t_im;
        let (d_re, d_im) = fpc_sub(a_re, a_im, b_re, b_im);
        f[2 * u + 1] = d_re;
        f[2 * u + 1 + hn] = d_im;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// FFT then iFFT recovers the input within tight `f64` tolerance.
    #[test]
    fn fft_round_trip_small() {
        for &logn in &[1u32, 2, 3, 4, 8, 9, 10] {
            let n = 1usize << logn;
            let mut a: Vec<f64> = (0..n).map(|i| (i as f64).sin() * 7.0 - 3.0).collect();
            let orig = a.clone();
            fft(&mut a, logn);
            ifft(&mut a, logn);
            for i in 0..n {
                assert!(
                    (a[i] - orig[i]).abs() < 1e-9,
                    "logn={logn} i={i} got={} want={}",
                    a[i],
                    orig[i]
                );
            }
        }
    }

    /// FFT-domain multiplication equals negacyclic schoolbook in
    /// coefficient domain, within tight `f64` tolerance.
    #[test]
    fn fft_multiply_matches_schoolbook() {
        for &logn in &[1u32, 2, 3, 4, 6] {
            let n = 1usize << logn;
            let a: Vec<f64> = (0..n).map(|i| ((i * 3 + 1) % 17) as f64 - 8.0).collect();
            let b: Vec<f64> = (0..n).map(|i| ((i * 5 + 2) % 19) as f64 - 9.0).collect();

            // Schoolbook negacyclic product mod x^n+1.
            let mut c = vec![0.0f64; n];
            for i in 0..n {
                for j in 0..n {
                    let k = i + j;
                    if k < n {
                        c[k] += a[i] * b[j];
                    } else {
                        c[k - n] -= a[i] * b[j];
                    }
                }
            }

            // FFT-domain product.
            let mut af = a.clone();
            let mut bf = b.clone();
            fft(&mut af, logn);
            fft(&mut bf, logn);
            poly_mul(&mut af, &bf, logn);
            ifft(&mut af, logn);

            for i in 0..n {
                assert!(
                    (af[i] - c[i]).abs() < 1e-9,
                    "logn={logn} i={i} fft={} school={}",
                    af[i],
                    c[i]
                );
            }
        }
    }

    /// `poly_split` then `poly_merge` is identity up to `f64` tolerance.
    #[test]
    fn split_merge_round_trip() {
        for &logn in &[2u32, 3, 4, 8] {
            let n = 1usize << logn;
            let hn = n >> 1;
            let mut f: Vec<f64> = (0..n).map(|i| ((i * 13 + 7) % 23) as f64 - 11.0).collect();
            fft(&mut f, logn);
            let orig = f.clone();
            let mut f0 = vec![0.0f64; hn];
            let mut f1 = vec![0.0f64; hn];
            poly_split(&mut f0, &mut f1, &f, logn);
            let mut g = vec![0.0f64; n];
            poly_merge(&mut g, &f0, &f1, logn);
            for i in 0..n {
                assert!(
                    (g[i] - orig[i]).abs() < 1e-9,
                    "logn={logn} i={i} got={} want={}",
                    g[i],
                    orig[i]
                );
            }
        }
    }

    /// `poly_mulselfadj` yields the real sum of squared moduli.
    #[test]
    fn mulselfadj_is_real_sum_of_squares() {
        let logn: u32 = 4;
        let n = 1usize << logn;
        let hn = n >> 1;
        let mut a: Vec<f64> = (0..n).map(|i| (i as f64).cos()).collect();
        fft(&mut a, logn);
        let snapshot = a.clone();
        poly_mulselfadj(&mut a, logn);
        for u in 0..hn {
            let want = snapshot[u] * snapshot[u] + snapshot[u + hn] * snapshot[u + hn];
            assert!((a[u] - want).abs() < 1e-12);
            assert_eq!(a[u + hn], 0.0);
        }
    }

    /// `poly_invnorm2(d, a, b)` is `1 / (|a|² + |b|²)` pointwise.
    #[test]
    fn invnorm2_is_reciprocal_of_squared_norms() {
        let logn: u32 = 3;
        let n = 1usize << logn;
        let hn = n >> 1;
        let mut a: Vec<f64> = (0..n).map(|i| 1.0 + i as f64).collect();
        let mut b: Vec<f64> = (0..n).map(|i| 2.0 - i as f64).collect();
        fft(&mut a, logn);
        fft(&mut b, logn);
        let mut d = vec![0.0f64; n];
        poly_invnorm2(&mut d, &a, &b, logn);
        for u in 0..hn {
            let want = 1.0
                / ((a[u] * a[u] + a[u + hn] * a[u + hn]) + (b[u] * b[u] + b[u + hn] * b[u + hn]));
            assert!((d[u] - want).abs() < 1e-12);
        }
    }

    /// `poly_mul_autoadj` scales both halves by the real factor `b[u]`.
    #[test]
    fn mul_autoadj_scales_both_halves() {
        let logn: u32 = 4;
        let n = 1usize << logn;
        let hn = n >> 1;
        let mut a: Vec<f64> = (0..n).map(|i| (i as f64) - 5.0).collect();
        let b: Vec<f64> = (0..n).map(|i| 0.5 + i as f64 * 0.25).collect();
        let snap = a.clone();
        poly_mul_autoadj(&mut a, &b, logn);
        for u in 0..hn {
            assert!((a[u] - snap[u] * b[u]).abs() < 1e-12);
            assert!((a[u + hn] - snap[u + hn] * b[u]).abs() < 1e-12);
        }
    }

    /// `poly_div_autoadj` is the inverse of `poly_mul_autoadj`.
    #[test]
    fn div_autoadj_inverts_mul_autoadj() {
        let logn: u32 = 4;
        let n = 1usize << logn;
        let mut a: Vec<f64> = (0..n).map(|i| (i as f64) - 7.5).collect();
        let b: Vec<f64> = (0..n).map(|i| 1.25 + i as f64 * 0.5).collect();
        let snap = a.clone();
        poly_mul_autoadj(&mut a, &b, logn);
        poly_div_autoadj(&mut a, &b, logn);
        for i in 0..n {
            assert!((a[i] - snap[i]).abs() < 1e-10);
        }
    }

    /// LDL yields `g11_new + L_10 * adj(L_10) * g00 == g11_old`.
    #[test]
    fn ldl_satisfies_schur_identity() {
        let logn: u32 = 4;
        let n = 1usize << logn;
        let hn = n >> 1;
        // Build a positive-definite 2x2 self-adjoint block:
        //   g00 = |x|² + 1,  g11 = |y|² + 1,  g01 = x * adj(y)
        let mut x: Vec<f64> = (0..n).map(|i| 0.5 + i as f64 * 0.1).collect();
        let mut y: Vec<f64> = (0..n).map(|i| 1.0 - i as f64 * 0.05).collect();
        fft(&mut x, logn);
        fft(&mut y, logn);

        let mut g00 = x.clone();
        poly_mulselfadj(&mut g00, logn);
        for u in 0..hn {
            g00[u] += 1.0;
        }

        let mut g11 = y.clone();
        poly_mulselfadj(&mut g11, logn);
        for u in 0..hn {
            g11[u] += 1.0;
        }

        let mut g01 = x.clone();
        poly_muladj(&mut g01, &y, logn);
        let g01_snapshot = g01.clone();
        let g11_snapshot = g11.clone();

        poly_ldl(&g00, &mut g01, &mut g11, logn);

        // Reconstruct g11_old = g11_new + L_10 * adj(g01_original).
        // After LDL: g01 now holds adj(L_10) (reference stores conj).
        // So L_10 = adj(g01).
        let mut l10 = g01.clone();
        poly_adj(&mut l10, logn);

        let mut check = l10.clone();
        poly_muladj(&mut check, &g01_snapshot, logn); // L_10 * adj(g01_orig)
        for u in 0..n {
            check[u] += g11[u];
        }
        for u in 0..n {
            assert!(
                (check[u] - g11_snapshot[u]).abs() < 1e-9,
                "u={u} got={} want={}",
                check[u],
                g11_snapshot[u]
            );
        }
    }
}
