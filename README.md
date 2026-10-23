# Navier–Stokes ternary shell model — invariant 2

This repository adds the explanation that the [first evaluation](https://github.com/JesseBrown1980/Checked-by-GPT-sol-6.1-Solved-Navier-Stokes-Global-Regularity-Finite-Time-Singularity-Invariant) omitted: **the exact all-time law of the unchanged 27-mode ternary shell update, with viscosity included.**

The update in [`Solved-Navier-Stokes…`](https://github.com/JesseBrown1980/Solved-Navier-Stokes-Global-Regularity-Finite-Time-Singularity-Invariant-Solved/tree/5147e98b7539b079c3af23b6ea1fda5ce78bc89f) at `5147e98` satisfies, for every state, every viscosity $\nu$ and every step $t$:

$$
Q(u_t) \equiv Q(u_0) - t\,\nu\,D \pmod P,\qquad D=\sum_{n=0}^{26}9^n \equiv 883\,973,
$$

so $I = Q + t\nu D$ is conserved, where $Q=\sum_n u_n$ is the DC lane of the canonical prism of the state. With $\nu = 0$, $Q$ itself is conserved. With $\nu \ne 0$, $Q$ cycles through every residue with exact period $P$. The proof, corollaries and scope are in **[PROOF.md](PROOF.md)**.

## Results (Rust 1.81.0, integer only, offline, WSL)

| Gate | Result |
|---|---|
| Clippy, all targets, `-D warnings` | PASS ([log](logs/clippy.log)) |
| Original upstream tests, unchanged ([`tests/upstream.rs`](tests/upstream.rs)) | **8 passed, 0 failed** |
| Canonical tribit kernel unit tests (vendored, unchanged) | **11 passed, 0 failed** |
| New contract tests ([`tests/contract.rs`](tests/contract.rs)) | **13 passed, 0 failed** |
| Release verifier: step closed form | 200,000 random (state, $\nu$) pairs, 0 mismatches |
| Release verifier: full period | 6 orbits × 1,000,081 upstream steps, 0 violations ([receipt](receipts/VERIFY.hbp)) |

Four of the contract tests record what the theorem does **not** give: zero is not a fixed point when $\nu > 0$; linear conservation is not quadratic energy cancellation; the README's Theorem 1 formula sums to −17; and the modular enstrophy bound is a range fact. The continuum Clay conclusion remains unverified by this argument (see [PROOF.md § Relation to the Clay problem](PROOF.md#relation-to-the-clay-problem)).

## Sources

Upstream code is vendored byte-for-byte and not modified. Hashes are compared with the git blobs in [`receipts/SOURCES.hbp`](receipts/SOURCES.hbp).

- `vendor/navier/`: Navier–Stokes crate at `5147e98b7539b079c3af23b6ea1fda5ce78bc89f` (Jesse Daniel Brown, Rayssa Chiqueto).
- `vendor/canonical/tribit/mod.rs`: canonical Asolaria tribit kernel from [ASI-OS](https://github.com/JesseBrown1980/Asolaria-Argon-Antigravity-Unification-ASI-OS) at `c0dfee6da17a4811964ebb4b9261d854075deac1`.
- `receipts/source-view/`: canonical ternary GitHub 3-D view of both source commits.

## Reproduce

```sh
rustup toolchain install 1.81.0 --profile minimal --component clippy
bash scripts/gate.sh        # clippy + all tests + full-period verifier; exit 0 = all pass
```
