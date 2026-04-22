//! NTT known-answer tests, backed by Python-generated ground truth.
//!
//! `A` and `B` are random polynomials in `Z_q[x]/(x^n+1)` with n=512.
//! `A_MUL_B` is their negacyclic product computed by schoolbook convolution
//! in Python. `A_NTT` is the forward NTT of `A`, letting us cross-check the
//! transform itself independently of the inverse.

#[path = "common/ntt_testvec.rs"]
mod tv;

use rfalco::ntt::{negacyclic_mul_512, ntt_512};
use rfalco::params::N_512;

fn to_arr(slice: &[u16]) -> [u16; N_512] {
    let mut a = [0u16; N_512];
    a.copy_from_slice(slice);
    a
}

#[test]
fn kat_negacyclic_product_matches_python_schoolbook() {
    let a = to_arr(&tv::A);
    let b = to_arr(&tv::B);
    let got = negacyclic_mul_512(&a, &b);
    assert_eq!(got.as_slice(), tv::A_MUL_B.as_slice());
}

#[test]
fn kat_forward_ntt_matches_python_reference() {
    let mut a = to_arr(&tv::A);
    ntt_512(&mut a);
    assert_eq!(a.as_slice(), tv::A_NTT.as_slice());
}
