# rfalco — project context

Pure-Rust FN-DSA (Falcon) signature library.

## Build discipline

- `cargo fmt --check` must pass before every commit.
- `cargo clippy --all-targets -- -D warnings` must pass (all + pedantic denied in `Cargo.toml`).
- `cargo test --release` must pass. Never run `cargo test` in debug for sign-off.
- RED-GREEN-BLUE cycle per task: write failing test first, make it pass, then refactor.
- Never relax tests to make them pass. Fix the implementation.
- No `unsafe` (forbidden at crate root).

## Specification

- Target: **NIST FIPS 206 / Falcon round-3 submission** (q = 12289).
- FIPS 206 is not yet published (expected late 2026/2027); round-3 package is the working reference.
- Serialization follow-on: **draft-ietf-cose-falcon-04** (in repo root). COSE Appendix A.2 hex is complete; JWK Appendix A.1 is truncated and **not usable as a KAT**.

## KATs

- Primary: NIST PQC round-3 KAT files (`PQCgenKAT_sign.rsp`) — bit-exact keygen + verify; sign is deterministic under reference FP sampler.
- Secondary: draft-04 Appendix A.2 COSE hex — serialization round-trip only.
- **Sign-KAT reconciliation with CT sampler**: the constant-time integer Gaussian sampler will NOT produce bit-exact signatures matching the FP reference. Sign KATs are verified via `verify(pk, msg, kat_sig) == accept`, not byte-equality. Keygen and verify KATs remain bit-exact.

## Side-channel posture

- Constant-time: all operations on secret data. Verified via integration tests (dudect-style timing + a no-secret-branch harness).
- Gaussian sampler: CT integer variant (Karney / cumulative-table rejection). Departs from reference FP sampler; see KAT note above.
- Masking: Boolean/arithmetic masking on secret polynomial state in sign path.
- Shuffling: randomized loop ordering where data-independent correctness permits.

## Formal verification

- Target: ESBMC over all code. This is research-grade work; expect to carve around floats and deep recursion.
- Harnesses in `tests/esbmc/` (C harnesses over `cbindgen`-generated headers, or direct `extern "C"` exports).
- Per-module bounded verification first; whole-pipeline verification gated behind module-level proofs.

## Module layout (`src/`)

| module    | role                                                 |
|-----------|------------------------------------------------------|
| `params`  | constants (Q, N, sizes, norm bounds)                 |
| `field`   | GF(12289) constant-time arithmetic                   |
| `shake`   | SHAKE256 XOF / PRNG matching round-3 seed expansion  |
| `ntt`     | NTT / INTT mod 12289                                 |
| `fft`     | tower-of-rings FFT for sign path                     |
| `sampler` | CT integer Gaussian sampler                          |
| `keygen`  | NTRUGen (integer, CT)                                |
| `sign`    | hash-to-point, ffSampling, norm check                |
| `verify`  | hash-to-point, NTT multiply, norm check              |
| `codec`   | sig/pk/sk bit-packing per round-3 spec               |
| `cose`    | COSE/JOSE serialization per draft-04                 |

## Open decisions

- ESBMC toolchain version + harness style (C shim vs direct) — settle before task #14.
- Masking order (first-order sufficient for this project; higher-order out of scope unless raised).

## User-facing protocol

- Caveman mode active globally; commits, code, docs remain normal prose.
- Auto mode: proceed on low-risk scaffolding + implementation; pause for spec-divergence or destructive actions.
