//! Constant-time arithmetic in the prime field GF(Q) with Q = 12289.
//!
//! All public operations are branch-free on the value bits of their inputs.
//! Inputs must already be reduced modulo `Q`; `debug_assert!` guards are
//! erased in release builds to preserve constant-time behaviour.

use crate::params::Q;

/// Precomputed Barrett multiplier: floor(2^32 / Q).
const Q_BAR: u64 = (1u64 << 32) / Q as u64;

/// Constant-time conditional subtraction: returns `s - m` when `s >= m`,
/// otherwise `s`. Branch-free on the value bits.
#[inline]
const fn csub(s: u32, m: u32) -> u32 {
    let t = s.wrapping_sub(m);
    // `t >> 31` is 1 iff `s < m` (underflow), else 0.
    // `mask` is all-ones when no underflow (so we pick `t`), zeros otherwise.
    let mask = (t >> 31).wrapping_sub(1);
    (t & mask) | (s & !mask)
}

/// Barrett reduction of `x` into `[0, Q)`. Correct for `x < Q^2 + Q`.
#[must_use]
#[inline]
pub const fn reduce(x: u32) -> u32 {
    let q_hat = ((x as u64).wrapping_mul(Q_BAR) >> 32) as u32;
    let r = x.wrapping_sub(q_hat.wrapping_mul(Q));
    csub(r, Q)
}

/// Modular addition in GF(Q).
#[must_use]
#[inline]
pub const fn add(a: u32, b: u32) -> u32 {
    debug_assert!(a < Q);
    debug_assert!(b < Q);
    csub(a + b, Q)
}

/// Modular subtraction in GF(Q).
#[must_use]
#[inline]
pub const fn sub(a: u32, b: u32) -> u32 {
    debug_assert!(a < Q);
    debug_assert!(b < Q);
    csub(a + Q - b, Q)
}

/// Modular multiplication in GF(Q).
#[must_use]
#[inline]
pub const fn mul(a: u32, b: u32) -> u32 {
    debug_assert!(a < Q);
    debug_assert!(b < Q);
    // `a * b < Q^2 ≈ 2^27.2`, so `u32` holds the product without overflow.
    reduce(a * b)
}

/// Modular negation in GF(Q).
#[must_use]
#[inline]
pub const fn neg(a: u32) -> u32 {
    debug_assert!(a < Q);
    // Branchless: when `a == 0` we must also return 0, not Q.
    let r = Q - a;
    csub(r, Q)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_wraps_at_q() {
        assert_eq!(add(Q - 1, 1), 0);
        assert_eq!(add(Q - 1, 2), 1);
        assert_eq!(add(0, 0), 0);
        assert_eq!(add(5000, 5000), 10_000);
        assert_eq!(add(6144, 6145), 0); // 6144 + 6145 = 12289 = Q
    }

    #[test]
    fn sub_borrows() {
        assert_eq!(sub(0, 1), Q - 1);
        assert_eq!(sub(Q - 1, Q - 1), 0);
        assert_eq!(sub(10, 3), 7);
        assert_eq!(sub(0, 0), 0);
    }

    #[test]
    fn neg_zero_is_zero() {
        assert_eq!(neg(0), 0);
        assert_eq!(neg(1), Q - 1);
        assert_eq!(neg(Q - 1), 1);
    }

    fn mod_ref(a: u32, b: u32) -> u32 {
        u32::try_from((u64::from(a) * u64::from(b)) % u64::from(Q)).expect("< Q fits u32")
    }

    fn reduce_ref(x: u32) -> u32 {
        u32::try_from(u64::from(x) % u64::from(Q)).expect("< Q fits u32")
    }

    #[test]
    fn mul_matches_reference() {
        for a in [0, 1, 2, Q - 1, 6144, 12_000] {
            for b in [0, 1, 2, Q - 1, 4096, 2048, 11_111] {
                assert_eq!(mul(a, b), mod_ref(a, b), "mul({a},{b})");
            }
        }
    }

    #[test]
    fn reduce_full_range_sample() {
        let qsq = (Q - 1) * (Q - 1);
        for x in [
            0u32,
            1,
            Q - 1,
            Q,
            Q + 1,
            2 * Q,
            qsq,
            1_000_000,
            50_000_000,
            150_000_000,
        ] {
            assert_eq!(reduce(x), reduce_ref(x), "reduce({x})");
        }
    }

    #[test]
    fn exhaustive_add_sub_roundtrip_sample() {
        // Exhaustive-ish spot check; full exhaustive left to a release-gated test.
        for a in (0..Q).step_by(97) {
            for b in (0..Q).step_by(131) {
                let s = add(a, b);
                assert_eq!(sub(s, b), a);
                assert_eq!(sub(s, a), b);
            }
        }
    }

    #[test]
    fn exhaustive_mul_against_u64_sparse() {
        for a in (0..Q).step_by(101) {
            for b in (0..Q).step_by(149) {
                assert_eq!(mul(a, b), mod_ref(a, b));
            }
        }
    }
}
