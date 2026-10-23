# All-time theorem for the unchanged ternary shell update

This document covers the implemented 27-mode finite-field system in [`navier-stokes-regularity-spectral`](https://github.com/JesseBrown1980/Solved-Navier-Stokes-Global-Regularity-Finite-Time-Singularity-Invariant-Solved/tree/5147e98b7539b079c3af23b6ea1fda5ce78bc89f) at commit `5147e98b7539b079c3af23b6ea1fda5ce78bc89f`. Its source is vendored byte-for-byte in [`vendor/navier`](vendor/navier/src/lib.rs). The canonical Asolaria tribit kernel comes from ASI-OS commit `c0dfee6da17a4811964ebb4b9261d854075deac1` (`canonical/asolaria-os/kernel/core/src/tribit/mod.rs`) and is included unchanged as [`vendor/canonical/tribit/mod.rs`](vendor/canonical/tribit/mod.rs). Neither file is edited. SHA-256 values are listed in [`receipts/SOURCES.hbp`](receipts/SOURCES.hbp).

The first evaluation showed that the transfer sums to zero. It did not say what that fact implies once viscosity is present, which is the case every upstream time-stepping test runs. The answer is a single exact identity that holds for every time step, every state and every viscosity value.

## Setting

- $P = 1\,000\,081$ (prime), $K = 27$, state $u \in \mathbb F_P^{27}$, viscosity $\nu \in \mathbb Z_{\ge 0}$. Upstream accepts raw `u64` lanes; every operation reduces them mod $P$.
- Shell total: $Q(u) = \sum_{n=0}^{26} u_n \bmod P$.
- Drift constant: $D = \sum_{n=0}^{26} 9^n \bmod P = 883\,973$. Since $8D \equiv 9^{27} - 1 \pmod P$ and $9^{27} \equiv 71\,218 \ne 1$, we have $D \neq 0$.

## Lemma 1 — transfer (all states)

For $1 \le n \le 25$, upstream `non_linear_transfer` computes $NL_n = u_n(u_{n+1} - u_{n-1}) \bmod P$. Let $S$ be the residue of $\sum_{n=1}^{25} NL_n$. It then sets $NL_0 = -\lfloor S/2 \rfloor$ and $NL_{26} = -S + \lfloor S/2 \rfloor$. Therefore $\sum_{n=0}^{26} NL_n \equiv 0$. Each floor term appears once with each sign, so they cancel.

## Lemma 2 — step in closed form (all states, all $\nu$)

Upstream `step` sets $u_n' = (u_n + NL_n) \bmod P$ and $\mathrm{diss}_n = (\nu \cdot 9^n) \bmod P$. It then takes `u' - diss` when `u' > diss`, and `(u' + P - diss) mod P` otherwise. Both values lie in $[0, P)$, so both branches compute the same residue. Hence

$$
u_n^{(t+1)} \equiv u_n^{(t)} + NL_n(u^{(t)}) - \nu\,9^n \pmod P .
$$

The subtracted term does not contain $u_n$.

## Theorem — exact affine drift and its conserved quantity

Every orbit $u^{(t+1)} = \mathrm{step}(u^{(t)})$ satisfies, for every $t \ge 0$:

$$
Q\big(u^{(t)}\big) \equiv Q\big(u^{(0)}\big) - t\,\nu\,D \pmod P,
\qquad
I(t,u) := Q(u) + t\,\nu\,D \ \text{ is constant along the orbit.}
$$

**Proof.** Sum Lemma 2 over $n$ and apply Lemma 1: $Q(u^{(t+1)}) = Q(u^{(t)}) + 0 - \nu D$. Induction on $t$ completes the proof. $\blacksquare$

**Corollaries.**

1. **Inviscid case.** If $\nu \equiv 0 \pmod P$, then $Q$ itself is conserved for all time.
2. **Viscous case.** If $\nu \not\equiv 0$, then $Q$ is *not* conserved. Because $\nu D$ is a unit of $\mathbb F_P$, the map $t \mapsto Q(u^{(t)})$ is a bijection of $\mathbb Z/P$ on each window of $P$ consecutive steps. So $Q$ has exact period $P$ and takes every residue once per period.
3. **Spectral reading.** $Q(u)$ is lane 0 (the DC glyph) of the canonical 27-point prism of the **state**. The DC lane of the **transfer** is always 0, which is the upstream Law 0 test. The DC lane of the state drifts by exactly $-\nu D$ per step. "Uncomputed centre = 0" therefore describes the transfer, not the evolving field.
4. **Finiteness.** Every lane is a residue in $[0, P)$ after every step, and the state space has $P^{27}$ elements, so every orbit is eventually periodic. Upstream assertions of the form `< P` check this, and only this.

## What the tests establish

| Statement | Where | How |
|---|---|---|
| $D = 883\,973$, $D \ne 0$, $8D \equiv 9^{27}-1$ | `drift_constant_matches_its_definition_and_is_nonzero` | three independent computations |
| Lemma 1 | `upstream_transfer_equals_closed_form_and_sums_to_zero` | upstream vs independent restatement: 24,000+ states (edge, uniform, raw `u64`, Host8-style trit rooms) |
| Lemma 2 | `upstream_step_equals_closed_form_for_every_tested_nu` and the `STEPFORM` receipt line | 3,500+ states × 9 viscosities, including $0, P-1, P, P+1, 2^{64}-1$; 200,000 random (state, $\nu$) pairs in release mode |
| Theorem, per step | `one_step_changes_q_by_exactly_minus_nu_d` | 5,000 random (state, $\nu$) pairs |
| Theorem along orbits | `invariant_holds_along_*` | the four upstream test initial conditions for 20,000 steps each; random, room and zero orbits |
| Corollary 1 | `inviscid_case_conserves_q_exactly` | 50 orbits × 2,000 steps |
| Corollary 2 | `drift_has_exact_period_p_in_time`, `ORBIT` receipt lines | **one full period of $P$ upstream steps** per orbit: every residue visited exactly once, $Q(P) = Q(0)$ |
| Corollary 3 | `q_is_dc_lane_of_state_and_transfer_dc_is_zero_canonical_prism` | upstream prism = canonical tribit prism; exact inverses |

The proof covers every time step by induction. The tests and the full-period receipt check that the shipped binary implements the operators the proof assumes.

## What the theorem does not give (also tested)

These are asserted as passing tests so they cannot silently drift:

- `scope_zero_state_is_not_fixed_under_positive_viscosity`: $u = 0$, $\nu = 1$ gives $u_0' = P - 1$, with all 27 lanes nonzero and $Q = -D$. The PDE term $-\nu |k|^2 \hat u_k$ maps 0 to 0. The implemented term does not.
- `scope_linear_conservation_is_not_quadratic_energy_cancellation`: for $u_0 = u_1 = 1$, $\sum NL_n = 0$ but $\sum u_n NL_n \equiv 500\,040 \ne 0$. Lemma 1 is a linear identity. It is not kinetic-energy conservation.
- `scope_readme_theorem_1_formula_does_not_sum_to_zero`: the README's Theorem 1 formula with $u_1 = u_2 = 1$ sums to $-17 \pmod P$. The code implements a different transfer, and the theorem above concerns the code.
- `scope_modular_enstrophy_bound_is_a_range_fact`: `total_enstrophy` is unchanged when $P$ is added to a lane. A residue below $P$ bounds no real quantity.

## Relation to the Clay problem

The theorem is exact and holds for all time, but it is a statement about a linear functional of a 27-lane state over $\mathbb F_P$. The [official problem statement](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf) asks for smooth real solutions of the 3-D incompressible PDE for arbitrary smooth divergence-free data. Connecting the two would still require everything listed in the [first evaluation, §7](https://github.com/JesseBrown1980/Checked-by-GPT-sol-6.1-Solved-Navier-Stokes-Global-Regularity-Finite-Time-Singularity-Invariant/blob/main/EXPLANATION.md): an encoding of real fields, operators consistent with the PDE, a passage to the continuum, and real-norm estimates controlling vorticity. Nothing in the vendored sources, the canonical tribit kernel, the Host8 room store or the theorem above supplies that bridge. The continuum conclusion therefore remains **unverified by this argument**. What *is* verified is that the implemented system has a precise, exactly computable all-time law, and that law is $I = Q + t\nu D$.
