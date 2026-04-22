//! SHAKE256 known-answer tests.
//!
//! Ground truth generated via `CPython` `hashlib.shake_256` (commit-time
//! snapshot). Each vector exercises a specific boundary in the sponge:
//! empty input, inputs shorter than/equal to/longer than the 136-byte rate,
//! and output lengths that span one, exactly-one, and multiple rate blocks.

use rfalco::shake::{shake256, Shake256, SHAKE256_RATE_BYTES};

fn hex(s: &str) -> Vec<u8> {
    assert!(s.len() % 2 == 0, "hex string must have even length");
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).expect("hex digit"))
        .collect()
}

fn shake(input: &[u8], n: usize) -> Vec<u8> {
    let mut out = vec![0u8; n];
    shake256(input, &mut out);
    out
}

#[test]
fn kat_empty_64() {
    let exp = hex(
        "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f\
         d75dc4ddd8c0f200cb05019d67b592f6fc821c49479ab48640292eacb3b7c4be",
    );
    assert_eq!(shake(b"", 64), exp);
}

#[test]
fn kat_abc_32() {
    let exp = hex("483366601360a8771c6863080cc4114d8db44530f8f1e1ee4f94ea37e78b5739");
    assert_eq!(shake(b"abc", 32), exp);
}

#[test]
fn kat_abc_136_exact_one_rate_block() {
    let exp = hex(
        "483366601360a8771c6863080cc4114d8db44530f8f1e1ee4f94ea37e78b5739\
         d5a15bef186a5386c75744c0527e1faa9f8726e462a12a4feb06bd8801e751e4\
         1385141204f329979fd3047a13c5657724ada64d2470157b3cdc288620944d78\
         dbcddbd912993f0913f164fb2ce95131a2d09a3e6d51cbfc622720d7a75c6334\
         e8a2d7ec71a7cc29",
    );
    assert_eq!(shake(b"abc", 136), exp);
}

#[test]
fn kat_abc_137_just_past_first_rate() {
    let exp = hex(
        "483366601360a8771c6863080cc4114d8db44530f8f1e1ee4f94ea37e78b5739\
         d5a15bef186a5386c75744c0527e1faa9f8726e462a12a4feb06bd8801e751e4\
         1385141204f329979fd3047a13c5657724ada64d2470157b3cdc288620944d78\
         dbcddbd912993f0913f164fb2ce95131a2d09a3e6d51cbfc622720d7a75c6334\
         e8a2d7ec71a7cc29cf",
    );
    assert_eq!(shake(b"abc", 137), exp);
}

#[test]
fn kat_abc_272_two_rate_blocks() {
    let exp = hex(
        "483366601360a8771c6863080cc4114d8db44530f8f1e1ee4f94ea37e78b5739\
         d5a15bef186a5386c75744c0527e1faa9f8726e462a12a4feb06bd8801e751e4\
         1385141204f329979fd3047a13c5657724ada64d2470157b3cdc288620944d78\
         dbcddbd912993f0913f164fb2ce95131a2d09a3e6d51cbfc622720d7a75c6334\
         e8a2d7ec71a7cc29cf0ea610eeff1a588290a53000faa79932becec0bd3cd0b3\
         3a7e5d397fed1ada9442b99903f4dcfd8559ed3950faf40fe6f3b5d710ed3b67\
         7513771af6bfe11934817e8762d9896ba579d88d84ba7aa3cdc7055f6796f195\
         bd9ae788f2f5bb96100d6bbaff7fbc6eea24d4449a2477d172a5507dcc931412\
         fc346b1bb39b878330e026b12ddf384a",
    );
    assert_eq!(shake(b"abc", 272), exp);
}

#[test]
fn kat_200_bytes_a3_64() {
    let exp = hex(
        "cd8a920ed141aa0407a22d59288652e9d9f1a7ee0c1e7c1ca699424da84a904d\
         2d700caae7396ece96604440577da4f3aa22aeb8857f961c4cd8e06f0ae6610b",
    );
    let input = [0xA3u8; 200];
    assert_eq!(shake(&input, 64), exp);
}

#[test]
fn kat_134_zeros_absorb_just_below_rate() {
    let exp = hex("6b580f0f88ae529dd7ab023e55a82bf71558d4021ce408ccee5cb79d1c20fd85");
    assert_eq!(shake(&[0u8; 134], 32), exp);
}

#[test]
fn kat_135_zeros_absorb_one_below_rate() {
    let exp = hex("4a6c0970c326babfaeef17f91988d1b4c5e95ed584c21b55b9f92e0d3671ddf9");
    assert_eq!(shake(&[0u8; 135], 32), exp);
}

#[test]
fn kat_136_zeros_absorb_exactly_one_rate() {
    let exp = hex("ea947b835fec1f9b0a7eabba901deb7881fd9999a1cbd5ccbb5a9afab7f6fe70");
    assert_eq!(shake(&[0u8; 136], 32), exp);
}

#[test]
fn kat_137_zeros_absorb_just_past_rate() {
    let exp = hex("60691a6b6b79c4abf99438b3f7a6455f2ce44fed8c8546cc90c218fe37ba5466");
    assert_eq!(shake(&[0u8; 137], 32), exp);
}

#[test]
fn kat_48_zero_seed_64() {
    // 48-byte all-zeros seed matches the NIST KAT framework seed width;
    // this vector anchors future Falcon KAT plumbing.
    let exp = hex(
        "eda313c95591a023a5b37f361c07a5753a92d3d0427459f34c7895d727d62816\
         b3aa2224eb9d823127d4f9f8a30fd7a1a02c6483d9c0f1fd41957b9ae4dfc63a",
    );
    assert_eq!(shake(&[0u8; 48], 64), exp);
}

#[test]
fn kat_48_ff_seed_40() {
    let exp = hex(
        "4c6bd2ee197ac771ab00b3040623880193377a21fe44f90b720d9d191f06db11\
         3da64afe09fbd90b",
    );
    assert_eq!(shake(&[0xffu8; 48], 40), exp);
}

#[test]
fn byte_by_byte_absorb_matches_bulk() {
    let msg = [0xA3u8; 200];
    let mut a = Shake256::new();
    for b in &msg {
        a.absorb(core::slice::from_ref(b));
    }
    a.finalize();
    let mut ao = [0u8; 64];
    a.squeeze(&mut ao);

    let mut b = Shake256::new();
    b.absorb(&msg);
    b.finalize();
    let mut bo = [0u8; 64];
    b.squeeze(&mut bo);

    assert_eq!(ao, bo);
}

#[test]
fn squeeze_spanning_many_blocks_agrees_with_single_block_steps() {
    let mut a = Shake256::new();
    a.absorb(b"long squeeze test");
    a.finalize();
    let mut big = vec![0u8; 5 * SHAKE256_RATE_BYTES + 7];
    a.squeeze(&mut big);

    let mut b = Shake256::new();
    b.absorb(b"long squeeze test");
    b.finalize();
    let mut stepped = vec![0u8; big.len()];
    let mut pos = 0;
    while pos < stepped.len() {
        let n = (SHAKE256_RATE_BYTES).min(stepped.len() - pos);
        b.squeeze(&mut stepped[pos..pos + n]);
        pos += n;
    }
    assert_eq!(big, stepped);
}
