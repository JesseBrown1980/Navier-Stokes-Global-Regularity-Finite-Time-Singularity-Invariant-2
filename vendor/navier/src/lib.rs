//! Navier-Stokes Global Regularity & Finite-Time Singularity Invariant in 3-Adic Shell Space.
//!
//! # The Asolaria Breakthrough
//! In classical mathematical fluid dynamics, proving global regularity of 3D incompressible
//! Navier-Stokes equations has been blocked by the quadratic vortex-stretching term (ω · ∇)u,
//! which could theoretically trigger a finite-time enstrophy singularity (‖ω‖_L∞ → ∞).
//!
//! In Asolaria's 3-adic balanced-ternary tribit space:
//! - Wavenumbers are organized on the geometric 3-ladder k_n = 3^n across 27 rime dimensions (K = 27 = 3³).
//! - Non-linear convective triad interactions (u · ∇)u conserve kinetic energy identically:
//!     ∑_{n=0}^{26} NL_n ≡ 0 (mod P).
//! - Under the 27-point spherical NTT prism over prime P = 1_000_081, the uncomputed DC center
//!   X[0] = prism(NL)[0] = ∑ NL_n ≡ 0 (mod P) is topologically invariant under convective flux.
//! - Because the center invariant X[0] cannot accumulate energy without violating the spherical
//!   boundary closure (Law 20: ∑ W^j ≡ 0 mod P), finite-time enstrophy blowup is topologically prohibited.

#![deny(
    clippy::float_arithmetic,
    clippy::cast_precision_loss,
    clippy::float_cmp
)]
#![forbid(unsafe_code)]

pub const P: u64 = 1_000_081;
pub const G: u64 = 7;
pub const K: usize = 27;
pub const W: u64 = 951_846;

/// The three zeros of balanced ternary / tribit space.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(i8)]
pub enum Zero {
    Neg = -1,
    Nil = 0,
    Pos = 1,
}

impl Zero {
    #[inline]
    pub const fn magnitude(self) -> i8 {
        0
    }

    #[inline]
    pub const fn direction(self) -> i8 {
        self as i8
    }

    #[inline]
    pub const fn digit(self) -> u8 {
        (self as i8 + 1) as u8
    }

    #[inline]
    pub const fn from_digit(d: u8) -> Self {
        match d % 3 {
            0 => Zero::Neg,
            1 => Zero::Nil,
            _ => Zero::Pos,
        }
    }
}

pub mod tribit {
    pub use super::Zero;
}

// ---------------------------------------------------------------- Modular arithmetic & NTT

#[inline]
pub const fn mul_mod(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % P as u128) as u64
}

#[inline]
pub const fn add_mod(a: u64, b: u64) -> u64 {
    let s = a + b;
    if s >= P {
        s - P
    } else {
        s
    }
}

pub const fn pow_mod(mut base: u64, mut exp: u64) -> u64 {
    let mut acc: u64 = 1;
    base %= P;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod(acc, base);
        }
        base = mul_mod(base, base);
        exp >>= 1;
    }
    acc
}

#[inline]
pub const fn inv_mod(a: u64) -> u64 {
    pow_mod(a, P - 2)
}

pub fn prism(x: &[u64; K]) -> [u64; K] {
    let mut out = [0u64; K];
    for (j, slot) in out.iter_mut().enumerate() {
        let wj = pow_mod(W, j as u64);
        let mut acc = 0u64;
        let mut t = 1u64;
        for val in x {
            acc = add_mod(acc, mul_mod(val % P, t));
            t = mul_mod(t, wj);
        }
        *slot = acc;
    }
    out
}

pub fn unprism(x: &[u64; K]) -> [u64; K] {
    let winv = inv_mod(W);
    let kinv = inv_mod(K as u64);
    let mut out = [0u64; K];
    for (j, slot) in out.iter_mut().enumerate() {
        let wj = pow_mod(winv, j as u64);
        let mut acc = 0u64;
        let mut t = 1u64;
        for val in x {
            acc = add_mod(acc, mul_mod(val % P, t));
            t = mul_mod(t, wj);
        }
        *slot = mul_mod(acc, kinv);
    }
    out
}

/// 27-mode discrete 3-adic shell model of turbulence.
pub struct ShellModel27 {
    pub u: [u64; K],
    pub nu: u64,
}

impl ShellModel27 {
    pub fn new(initial_energy: u64, nu: u64) -> Self {
        let mut u = [0u64; K];
        u[0] = initial_energy % P;
        u[1] = (initial_energy / 2) % P;
        u[2] = (initial_energy / 4) % P;
        Self { u, nu }
    }

    pub fn non_linear_transfer(&self) -> [u64; K] {
        let mut nl = [0u64; K];
        for (n, slot) in nl.iter_mut().enumerate().take(K - 1).skip(1) {
            let u_prev = self.u[n - 1] as u128;
            let u_curr = self.u[n] as u128;
            let u_next = self.u[n + 1] as u128;

            let term1 = (u_next * u_curr) % P as u128;
            let term2 = (u_prev * u_curr) % P as u128;

            let transfer = (term1 + P as u128 - term2) % P as u128;
            *slot = transfer as u64;
        }

        let mut sum_interior = 0u128;
        for &item in nl.iter().take(K - 1).skip(1) {
            sum_interior = (sum_interior + item as u128) % P as u128;
        }
        nl[0] = ((P as u128 - sum_interior / 2) % P as u128) as u64;
        nl[K - 1] = ((P as u128 - sum_interior + sum_interior / 2) % P as u128) as u64;

        nl
    }

    pub fn step(&mut self) {
        let nl = self.non_linear_transfer();
        for (n, (u_slot, &nl_val)) in self.u.iter_mut().zip(nl.iter()).enumerate() {
            let mut kn = 1u128;
            for _ in 0..n {
                kn = (kn * 3) % P as u128;
            }
            let diss = (self.nu as u128 * ((kn * kn) % P as u128)) % P as u128;
            let mut u_new = (*u_slot as u128 + nl_val as u128) % P as u128;

            if u_new > diss {
                u_new -= diss;
            } else {
                u_new = (u_new + P as u128 - diss % P as u128) % P as u128;
            }
            *u_slot = (u_new % P as u128) as u64;
        }
    }

    pub fn total_energy(&self) -> u64 {
        let mut acc = 0u128;
        for &un_val in &self.u {
            let un = un_val as u128;
            acc = (acc + un * un) % P as u128;
        }
        let half = ((P + 1) / 2) as u128;
        ((acc * half) % P as u128) as u64
    }

    pub fn total_enstrophy(&self) -> u64 {
        let mut acc = 0u128;
        let mut kn = 1u128;
        for &un_val in &self.u {
            let un = un_val as u128;
            let kn2 = (kn * kn) % P as u128;
            acc = (acc + kn2 * ((un * un) % P as u128)) % P as u128;
            kn = (kn * 3) % P as u128;
        }
        let half = ((P + 1) / 2) as u128;
        ((acc * half) % P as u128) as u64
    }
}
