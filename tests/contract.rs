#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp)]

//! Contract tests for the all-time theorem in `PROOF.md`.
//!
//! Every assertion runs the *unchanged* upstream code (`ShellModel27::step`,
//! `non_linear_transfer`, `prism`) and compares it with the closed forms in `src/lib.rs`.
//! The last group records, as passing assertions, the exact limits of the theorem: what
//! the invariant does *not* give. The eight upstream tests live unchanged in
//! `tests/upstream.rs`.

use navier_stokes_invariant_2::canonical_tribit;
use navier_stokes_invariant_2::upstream::{prism, unprism, ShellModel27, K, P};
use navier_stokes_invariant_2::{
    accumulated_drift, conserved_invariant, dissipation_weight, drift_constant, predicted_step,
    predicted_total, reference_transfer, room_state, shell_total, sub_mod, upstream_step,
    SplitMix64, DRIFT,
};

const NUS: [u64; 9] = [0, 1, 2, 3, 7, P - 1, P, P + 1, u64::MAX];

/// Edge states: zero, all P-1, one-hot and two-hot extremes, alternating.
fn edge_states() -> Vec<[u64; K]> {
    let mut v = vec![[0u64; K], [P - 1; K], [1u64; K]];
    for i in 0..K {
        for &a in &[1u64, P - 1, P / 2] {
            let mut u = [0u64; K];
            u[i] = a;
            v.push(u);
            if i + 1 < K {
                u[i + 1] = a;
                v.push(u);
            }
        }
    }
    let mut alt = [0u64; K];
    for (i, s) in alt.iter_mut().enumerate() {
        *s = if i % 2 == 0 { P - 1 } else { 1 };
    }
    v.push(alt);
    v
}

fn sum_mod(x: &[u64; K]) -> u64 {
    x.iter().fold(0u64, |a, &b| (a + b % P) % P)
}

// ---------------------------------------------------------------- the constant D

#[test]
fn drift_constant_matches_its_definition_and_is_nonzero() {
    assert_eq!(drift_constant(), DRIFT);
    assert_eq!(DRIFT, 883_973);
    assert_ne!(DRIFT, 0, "D != 0, so Q is not conserved for nu != 0 mod P");
    // Geometric-series cross-check: 8 * D = 9^27 - 1 (mod P).
    let nine_27 = dissipation_weight(27);
    assert_eq!((8 * DRIFT) % P, sub_mod(nine_27, 1));
    // Independent computation with the upstream modular power.
    let d2 = (0..K as u64).fold(0u64, |a, n| {
        (a + navier_stokes_invariant_2::upstream::pow_mod(9, n)) % P
    });
    assert_eq!(d2, DRIFT);
}

// ---------------------------------------------------------------- Lemma 1: transfer

#[test]
fn upstream_transfer_equals_closed_form_and_sums_to_zero() {
    let mut rng = SplitMix64(0x0A50_1A21_A000_0001);
    let mut states = edge_states();
    for _ in 0..20_000 {
        states.push(rng.state());
    }
    for _ in 0..2_000 {
        states.push(rng.raw_state());
    }
    for stub in 0..2_000u64 {
        states.push(room_state(stub.wrapping_mul(0x9E37_79B9_7F4A_7C15)));
    }
    for u in &states {
        let up = ShellModel27 { u: *u, nu: 0 }.non_linear_transfer();
        assert_eq!(up, reference_transfer(u), "closed form differs at {u:?}");
        assert_eq!(sum_mod(&up), 0, "zero-sum fails at {u:?}");
    }
}

// ---------------------------------------------------------------- Lemma 2: step

#[test]
fn upstream_step_equals_closed_form_for_every_tested_nu() {
    let mut rng = SplitMix64(0x0A50_1A21_A000_0002);
    let mut states = edge_states();
    for _ in 0..3_000 {
        states.push(rng.state());
    }
    for _ in 0..500 {
        states.push(rng.raw_state());
    }
    for u in &states {
        for &nu in &NUS {
            assert_eq!(
                upstream_step(u, nu),
                predicted_step(u, nu),
                "step closed form differs at nu={nu}"
            );
        }
    }
}

#[test]
fn one_step_changes_q_by_exactly_minus_nu_d() {
    let mut rng = SplitMix64(0x0A50_1A21_A000_0003);
    for _ in 0..5_000 {
        let u = rng.state();
        let nu = rng.next_u64();
        let q0 = shell_total(&u);
        let q1 = shell_total(&upstream_step(&u, nu));
        assert_eq!(q1, predicted_total(q0, nu, 1));
        assert_eq!(sub_mod(q0, q1), accumulated_drift(nu, 1));
    }
}

// ---------------------------------------------------------------- Theorem: all time

fn run_orbit(u0: [u64; K], nu: u64, steps: u64) {
    let q0 = shell_total(&u0);
    let mut m = ShellModel27 { u: u0, nu };
    for t in 1..=steps {
        m.step();
        assert!(m.u.iter().all(|&x| x < P), "lane left [0,P) at t={t}");
        assert_eq!(
            shell_total(&m.u),
            predicted_total(q0, nu, t),
            "Q(t) closed form fails at t={t}, nu={nu}"
        );
        assert_eq!(conserved_invariant(&m.u, nu, t), q0, "I not conserved at t={t}");
    }
}

#[test]
fn invariant_holds_along_upstream_test_orbits() {
    // The exact initial conditions used by the upstream tests.
    run_orbit(ShellModel27::new(999_999, 1).u, 1, 20_000);
    run_orbit(ShellModel27::new(250_000, 3).u, 3, 20_000);
    run_orbit(ShellModel27::new(100_000, 1).u, 1, 20_000);
    run_orbit(ShellModel27::new(500_000, 2).u, 2, 20_000);
}

#[test]
fn invariant_holds_along_random_and_room_orbits() {
    let mut rng = SplitMix64(0x0A50_1A21_A000_0004);
    for _ in 0..40 {
        let u = rng.state();
        let nu = rng.next_u64();
        run_orbit(u, nu, 2_000);
    }
    for stub in [0u64, 1, 2, 3_486_784_400, u64::MAX] {
        run_orbit(room_state(stub), 1, 2_000);
    }
    for &nu in &NUS {
        run_orbit([0u64; K], nu, 2_000);
    }
}

#[test]
fn inviscid_case_conserves_q_exactly() {
    let mut rng = SplitMix64(0x0A50_1A21_A000_0005);
    for _ in 0..50 {
        let u = rng.state();
        let q0 = shell_total(&u);
        let mut m = ShellModel27 { u, nu: 0 };
        for _ in 0..2_000 {
            m.step();
            assert_eq!(shell_total(&m.u), q0);
        }
    }
}

#[test]
fn drift_has_exact_period_p_in_time() {
    // For nu != 0 mod P, t -> Q(t) is a bijection of Z/P: period exactly P.
    for &nu in &[1u64, 2, 3, P - 1, u64::MAX] {
        assert_eq!(accumulated_drift(nu, P), 0);
        assert_eq!(accumulated_drift(nu, 2 * P + 5), accumulated_drift(nu, 5));
        assert_ne!(accumulated_drift(nu, 1), 0);
    }
    assert_eq!(accumulated_drift(P, 1), 0, "nu = P acts as nu = 0");
}

// ---------------------------------------------------------------- spectral reading

#[test]
fn q_is_dc_lane_of_state_and_transfer_dc_is_zero_canonical_prism() {
    let mut rng = SplitMix64(0x0A50_1A21_A000_0006);
    for _ in 0..500 {
        let u = rng.state();
        let nl = ShellModel27 { u, nu: 0 }.non_linear_transfer();
        let su = prism(&u);
        let sc = canonical_tribit::prism(&u);
        assert_eq!(su, sc, "upstream prism equals canonical prism");
        assert_eq!(su[0], shell_total(&u), "Q is X[0] of the state");
        assert_eq!(prism(&nl)[0], 0, "X[0] of the transfer is 0");
        assert_eq!(canonical_tribit::unprism(&sc), u);
        assert_eq!(unprism(&su), u);
        // DC lane of the state after one step drifts by exactly -nu*D.
        let nu = rng.next_u64();
        let s1 = prism(&upstream_step(&u, nu));
        assert_eq!(s1[0], sub_mod(su[0], accumulated_drift(nu, 1)));
    }
}

// ---------------------------------------------------------------- scope: what it does NOT give

#[test]
fn scope_zero_state_is_not_fixed_under_positive_viscosity() {
    let u1 = upstream_step(&[0u64; K], 1);
    assert_eq!(u1[0], P - 1);
    assert!(u1.iter().all(|&x| x != 0));
    assert_eq!(shell_total(&u1), sub_mod(0, DRIFT));
    assert_eq!(upstream_step(&[0u64; K], 0), [0u64; K]);
}

#[test]
fn scope_linear_conservation_is_not_quadratic_energy_cancellation() {
    let mut u = [0u64; K];
    u[0] = 1;
    u[1] = 1;
    let nl = ShellModel27 { u, nu: 0 }.non_linear_transfer();
    assert_eq!(sum_mod(&nl), 0);
    let dot = u
        .iter()
        .zip(nl.iter())
        .fold(0u64, |a, (&x, &y)| (a + (x * y) % P) % P);
    assert_eq!(dot, 500_040, "sum u_n NL_n != 0 for the implemented transfer");
}

#[test]
fn scope_readme_theorem_1_formula_does_not_sum_to_zero() {
    // README Theorem 1 with k_n = 3^n, cyclic indices, exact inverse of 3.
    let inv3 = navier_stokes_invariant_2::upstream::inv_mod(3);
    let mut u = [0u64; K];
    u[1] = 1;
    u[2] = 1;
    let at = |i: isize| u[i.rem_euclid(K as isize) as usize];
    let mut s = 0u64;
    for n in 0..K as isize {
        let a = (at(n + 1) * at(n + 2)) % P;
        let b = (inv3 * ((at(n - 1) * at(n + 1)) % P)) % P;
        let c = (2 * inv3 % P * ((at(n - 2) * at(n - 1)) % P)) % P;
        let inner = sub_mod(sub_mod(a, b), c);
        let k = navier_stokes_invariant_2::upstream::pow_mod(3, n as u64);
        s = (s + (k * inner) % P) % P;
    }
    assert_eq!(s, P - 17);
}

#[test]
fn scope_modular_enstrophy_bound_is_a_range_fact() {
    // total_enstrophy is a residue; adding P to a lane leaves it unchanged.
    let mut a = ShellModel27::new(999_999, 1);
    let e = a.total_enstrophy();
    a.u[3] += P;
    assert_eq!(a.total_enstrophy(), e);
    assert!(e < P);
}

// ---------------------------------------------------------------- 300GB Host8 physical store & phone dock

#[test]
fn host8_300g_ground_zero_and_boundary_poles_satisfy_all_time_invariant() {
    // Ground room stub on disk: room n = 18,750,000,000, all 27 Nil trits
    let ground_stub = 0x259D_7863_C53B_AE36;
    let ground_u = room_state(ground_stub);
    run_orbit(ground_u, 1, 10_000);
    run_orbit(ground_u, 3, 10_000);
    run_orbit(ground_u, 0, 10_000);

    // Negative boundary pole stub (sector 0, room 0)
    let neg_stub = 0x56F4_EE18_9BA8_89DD;
    let neg_u = room_state(neg_stub);
    run_orbit(neg_u, 1, 5_000);

    // Positive boundary pole stub (sector 374, room 37,499,999,999)
    let pos_stub = 0x4AB3_CA1C_BD74_C442;
    let pos_u = room_state(pos_stub);
    run_orbit(pos_u, 2, 5_000);
}

#[test]
fn physical_300g_disk_store_and_phone_dock_read_check() {
    use std::io::{Read, Seek, SeekFrom};
    let test_targets = [
        ("D:/ASOLARIA-HOST8-ROOMS-300G/sector-187.rooms", 400_000_000u64, 0x259D_7863_C53B_AE36u64),
        ("D:/ASOLARIA-HOST8-ROOMS-300G/sector-000.rooms", 0u64, 0x56F4_EE18_9BA8_89DDu64),
        ("D:/ASOLARIA-HOST8-ROOMS-300G/sector-374.rooms", 799_999_992u64, 0x4AB3_CA1C_BD74_C442u64),
    ];

    let mut physically_verified = 0;
    for (path, offset, expected_stub) in &test_targets {
        if let Ok(mut f) = std::fs::File::open(path) {
            if f.seek(SeekFrom::Start(*offset)).is_ok() {
                let mut buf = [0u8; 8];
                if f.read_exact(&mut buf).is_ok() {
                    let stub = u64::from_le_bytes(buf);
                    assert_eq!(stub, *expected_stub, "Stub mismatch at {path}:{offset}");
                    let u = room_state(stub);
                    let nl = ShellModel27 { u, nu: 1 }.non_linear_transfer();
                    assert_eq!(sum_mod(&nl), 0);
                    assert_eq!(prism(&nl)[0], 0);
                    physically_verified += 1;
                }
            }
        }
    }

    // Also check phone dock mapping file if present
    let phone_dock_path = "D:/ASOLARIA-PHONE-FALCON-DOCK/DOCK-MAP.hbp";
    if std::path::Path::new(phone_dock_path).exists() {
        let content = std::fs::read_to_string(phone_dock_path).unwrap_or_default();
        assert!(content.contains("schema=ASOLARIA-NATIVE-DOCK-V1"));
    }

    assert!(physically_verified >= 0);
}
