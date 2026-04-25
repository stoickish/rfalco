# rfalco

Pure-Rust implementation of the **FN-DSA** (Falcon) digital signature
algorithm. Research-grade; not yet production-ready.

## Status

Work in progress. The pipeline is being built in small, test-verified
increments. As of the latest commit:

| Module    | Status      | Notes                                             |
|-----------|-------------|---------------------------------------------------|
| `params`  | done        | constants, key/sig sizes, norm bounds             |
| `field`   | done        | GF(12289) constant-time arithmetic (Barrett)      |
| `shake`   | done        | SHAKE256 sponge, zero-dep, matches round-3 API    |
| `ntt`     | done        | Cooley-Tukey fwd + Gentleman-Sande inv, n=512,1024|
| `codec`   | done        | pk/sig/sk bit packing per round-3 layout          |
| `drbg`    | done        | AES-256 CTR_DRBG for NIST PQC KAT reproduction    |
| `verify`  | done        | hash-to-point, NTT multiply, norm check           |
| `sampler` | done        | FP `SamplerZ`, bit-exact vs PQClean reference     |
| `fft`     | done        | tower-of-rings FFT, `f64` bit-exact vs round-3    |
| `keygen`  | pending     | NTRUGen (integer, CT)                             |
| `sign`    | pending     | hash-to-point, ffSampling, norm check             |
| `cose`    | pending     | COSE/JOSE serialization per draft-04              |

## Design goals

1. **Spec fidelity**: target NIST FIPS 206 (Falcon round-3 working
   reference until FIPS 206 is published). Parameter set: `q = 12289`,
   `n ∈ {512, 1024}`.
2. **Constant-time on secrets (integer path)**: every integer operation
   touching secret data is branch-free in control flow and memory access.
   FP operations on secret data (the ffSampling tree and `SamplerZ`)
   follow the round-3 reference's isochronous sequences — same posture as
   upstream, microarchitectural FP timing quirks acknowledged.
3. **Side-channel countermeasures**: Boolean/arithmetic masking on the
   integer secret polynomial state in the sign path; randomized loop
   ordering where data-independence allows.
4. **Formal verification**: ESBMC over all code. Per-module bounded
   verification first; whole-pipeline verification gated behind module
   proofs. Harnesses live in `tests/esbmc/`.
5. **No `unsafe`**: forbidden at the crate root.

## KAT strategy

Two sources of ground truth:

- **Primary** — NIST PQC round-3 KAT files (`PQCgenKAT_sign.rsp`):
  **bit-exact for keygen, sign, and verify**.
- **Secondary** — draft-ietf-cose-falcon-04 Appendix A.2 COSE hex:
  serialization round-trip only.

Sign bit-exactness requires reproducing the reference's IEEE 754
binary64 operation sequence: no `fast-math`, no reassociation, no FMA
contraction unless the reference performs one at that call site. The
`sampler` and `fft` modules use `f64` directly and follow the reference
op order line-for-line.

## Build & test

```sh
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test --release
```

All three are required to pass before any commit. `cargo test` is
always run in `--release`; debug builds are never a sign-off target.

Toolchain pinned via `rust-toolchain.toml`; formatter pinned via
`rustfmt.toml`.

## Layout

```
src/
  params.rs    constants
  field.rs     GF(12289) arithmetic
  shake.rs     SHAKE256
  ntt/         negacyclic NTT + twiddle tables
  fft.rs       tower-of-rings FFT (pending)
  sampler.rs   FP discrete Gaussian (done, bit-exact)
  keygen.rs    NTRUGen (pending)
  sign.rs      signing (pending)
  verify.rs    verification (pending)
  codec.rs     bit-level serialization
  cose.rs      COSE/JOSE (pending)
tests/
  common/      shared test vectors
  *_kat.rs     known-answer tests per module
```

## License

Dual-licensed under MIT or Apache-2.0, at your option.
