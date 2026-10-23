//! The explicit mathematical contract of the unchanged ternary shell dynamics.
//!
//! The upstream crate (`navier-stokes-regularity-spectral`, commit `5147e98`) is vendored
//! byte-for-byte and is not modified. This library restates its update in closed form and
//! exposes the quantities needed to state one exact, all-time theorem about it:
//!
//! ```text
//! Q(u)      = sum_{n=0}^{26} u_n                      (mod P)   -- shell-coordinate total
//! D         = sum_{n=0}^{26} 9^n  = 883_973           (mod P)   -- viscous drift per unit nu
//! Q(u_t)    = Q(u_0) - t * nu * D                     (mod P)   for every t >= 0
//! I(t, u_t) = Q(u_t) + t * nu * D = Q(u_0)            (mod P)   conserved for every t >= 0
//! ```
//!
//! `Q(u)` is exactly lane 0 (the DC glyph) of the canonical 27-point prism of the *state*.
//! Q is a linear shell-coordinate total. No identification with a real fluid norm or
//! continuum velocity is assumed by this library.

#![forbid(unsafe_code)]
#![deny(
    clippy::float_arithmetic,
    clippy::float_cmp,
    clippy::cast_precision_loss
)]

pub use navier_stokes_regularity_spectral as upstream;
use upstream::{ShellModel27, K, P};

/// The canonical Asolaria tribit kernel module (ASI-OS commit `c0dfee6`,
/// `canonical/asolaria-os/kernel/core/src/tribit/mod.rs`), included unchanged so the
/// prism used for cross-checks is the canonical one rather than a re-implementation.
#[path = "../vendor/canonical/tribit/mod.rs"]
pub mod canonical_tribit;

/// `D = sum_{n=0}^{26} 9^n mod P`. Checked against [`drift_constant`] in the tests.
pub const DRIFT: u64 = 883_973;

#[inline]
const fn mul(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % P as u128) as u64
}

#[inline]
const fn add(a: u64, b: u64) -> u64 {
    ((a as u128 + b as u128) % P as u128) as u64
}

/// `a - b mod P` for any `a, b` (both reduced first).
#[inline]
pub const fn sub_mod(a: u64, b: u64) -> u64 {
    let a = a % P;
    let b = b % P;
    if a >= b {
        a - b
    } else {
        a + P - b
    }
}

/// `9^n mod P` — the coefficient `k_n^2` with `k_n = 3^n` used by the upstream `step`.
pub const fn dissipation_weight(n: usize) -> u64 {
    let mut acc = 1u64;
    let mut i = 0;
    while i < n {
        acc = mul(acc, 9);
        i += 1;
    }
    acc
}

/// `D = sum_{n=0}^{26} 9^n mod P`, computed from its definition.
pub const fn drift_constant() -> u64 {
    let mut acc = 0u64;
    let mut n = 0;
    while n < K {
        acc = add(acc, dissipation_weight(n));
        n += 1;
    }
    acc
}

/// `Q(u) = sum_n u_n mod P`.
pub fn shell_total(u: &[u64; K]) -> u64 {
    u.iter().fold(0u64, |acc, &x| add(acc, x % P))
}

/// Independent restatement of the upstream transfer from its formula:
/// `NL_n = u_n (u_{n+1} - u_{n-1})` for `1 <= n <= 25`, with the interior residue
/// `S` split as `NL_0 = -floor(S/2)`, `NL_26 = -S + floor(S/2)`.
pub fn reference_transfer(u: &[u64; K]) -> [u64; K] {
    let mut nl = [0u64; K];
    for n in 1..K - 1 {
        nl[n] = mul(u[n], sub_mod(u[n + 1], u[n - 1]));
    }
    let s = nl[1..K - 1].iter().fold(0u64, |acc, &x| add(acc, x));
    let half = s / 2;
    nl[0] = sub_mod(0, half);
    nl[K - 1] = add(sub_mod(0, s), half);
    nl
}

/// Closed form of one upstream `step`: `u_n + NL_n - nu * 9^n (mod P)`.
pub fn predicted_step(u: &[u64; K], nu: u64) -> [u64; K] {
    let nl = reference_transfer(u);
    let nu = nu % P;
    let mut out = [0u64; K];
    for n in 0..K {
        out[n] = sub_mod(add(u[n] % P, nl[n]), mul(nu, dissipation_weight(n)));
    }
    out
}

/// `t * nu * D mod P` — the exact accumulated drift of `Q` after `t` steps.
pub fn accumulated_drift(nu: u64, t: u64) -> u64 {
    mul(mul(t % P, nu % P), DRIFT)
}

/// `Q(u_0) - t * nu * D mod P` — the predicted shell total after `t` steps.
pub fn predicted_total(q0: u64, nu: u64, t: u64) -> u64 {
    sub_mod(q0, accumulated_drift(nu, t))
}

/// `I(t, u) = Q(u) + t * nu * D mod P` — constant along every forward orbit.
pub fn conserved_invariant(u: &[u64; K], nu: u64, t: u64) -> u64 {
    add(shell_total(u), accumulated_drift(nu, t))
}

/// Runs the *unchanged upstream* `ShellModel27::step` once and returns the new state.
pub fn upstream_step(u: &[u64; K], nu: u64) -> [u64; K] {
    let mut m = ShellModel27 { u: *u, nu };
    m.step();
    m.u
}

/// Deterministic integer generator (SplitMix64) for reproducible state families.
pub struct SplitMix64(pub u64);

impl SplitMix64 {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A state with every lane uniformly reduced into `[0, P)`.
    pub fn state(&mut self) -> [u64; K] {
        let mut u = [0u64; K];
        for slot in &mut u {
            *slot = self.next_u64() % P;
        }
        u
    }

    /// A state with raw, unreduced `u64` lanes (the upstream API accepts these).
    pub fn raw_state(&mut self) -> [u64; K] {
        let mut u = [0u64; K];
        for slot in &mut u {
            *slot = self.next_u64();
        }
        u
    }
}

/// Unpacks a 27-trit room address into lanes with the same amplitude map as the
/// upstream Host8 test, using the canonical `Zero` type.
pub fn room_state(stub: u64) -> [u64; K] {
    use canonical_tribit::Zero;
    let mut u = [0u64; K];
    let mut v = stub;
    for (i, slot) in u.iter_mut().enumerate() {
        let z = Zero::from_digit((v % 3) as u8);
        v /= 3;
        let amp = match z {
            Zero::Neg => P - 333_360,
            Zero::Nil => 0,
            Zero::Pos => 333_360,
        };
        *slot = (amp + stub.wrapping_add(i as u64) % 10_000) % P;
    }
    u
}
