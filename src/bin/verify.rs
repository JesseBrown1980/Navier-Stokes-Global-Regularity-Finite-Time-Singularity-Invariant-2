#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp)]

//! Release-mode verifier. Runs the *unchanged* upstream `ShellModel27::step` for one full
//! drift period (`P` steps) on each listed orbit and checks, at every step:
//!   * every lane is a residue in `[0, P)`;
//!   * `Q(u_t) == Q(u_0) - t*nu*D (mod P)`;
//!   * `I(t, u_t) == Q(u_0)`;
//!   * `Q` visits each residue of `Z/P` exactly once over `t in [0, P)` when `nu != 0 mod P`.
//!
//! Output is HBP (`TAG|k=v|...|json=0`) on stdout. Exit 0 iff every check holds.

use navier_stokes_invariant_2::upstream::{ShellModel27, K, P};
use navier_stokes_invariant_2::{
    conserved_invariant, drift_constant, predicted_step, predicted_total, room_state,
    shell_total, upstream_step, SplitMix64, DRIFT,
};
use std::process::ExitCode;

struct Orbit {
    name: &'static str,
    u0: [u64; K],
    nu: u64,
}

fn run_full_period(o: &Orbit) -> (bool, u64, u64, u64) {
    let q0 = shell_total(&o.u0);
    let mut seen = vec![0u64; (P as usize).div_ceil(64)];
    let mut distinct = 0u64;
    let mut mark = |q: u64, seen: &mut Vec<u64>| {
        let (w, b) = ((q / 64) as usize, q % 64);
        if seen[w] & (1 << b) == 0 {
            seen[w] |= 1 << b;
            distinct += 1;
        }
    };
    mark(q0, &mut seen);
    let mut m = ShellModel27 { u: o.u0, nu: o.nu };
    let mut ok = true;
    let mut first_bad = 0u64;
    for t in 1..=P {
        m.step();
        let q = shell_total(&m.u);
        let good = m.u.iter().all(|&x| x < P)
            && q == predicted_total(q0, o.nu, t)
            && conserved_invariant(&m.u, o.nu, t) == q0;
        if !good && ok {
            ok = false;
            first_bad = t;
        }
        if t < P {
            mark(q, &mut seen);
        }
    }
    let q_end = shell_total(&m.u);
    let nu_zero = o.nu % P == 0;
    let expected_distinct = if nu_zero { 1 } else { P };
    ok = ok && q_end == q0 && distinct == expected_distinct;
    (ok, distinct, q_end, first_bad)
}

fn main() -> ExitCode {
    let mut all_ok = true;
    println!(
        "VERIFYHDR|crate=navier-stokes-invariant-2|upstream=5147e98b7539b079c3af23b6ea1fda5ce78bc89f|canonical=c0dfee6da17a4811964ebb4b9261d854075deac1|P={P}|K={K}|arithmetic=integer_only|float_used=0|json=0"
    );

    let d_ok = drift_constant() == DRIFT && DRIFT != 0;
    all_ok &= d_ok;
    println!("DRIFT|D={DRIFT}|definition=sum_n_9^n_mod_P|recomputed={}|nonzero={}|pass={}|json=0", drift_constant(), u8::from(DRIFT != 0), u8::from(d_ok));

    // Closed-form step check over a reproducible family.
    let mut rng = SplitMix64(0x0A50_1A21_B000_0001);
    let mut checked = 0u64;
    let mut mismatches = 0u64;
    for _ in 0..200_000 {
        let u = rng.state();
        let nu = rng.next_u64();
        checked += 1;
        if upstream_step(&u, nu) != predicted_step(&u, nu) {
            mismatches += 1;
        }
    }
    let step_ok = mismatches == 0;
    all_ok &= step_ok;
    println!("STEPFORM|claim=step(u)=u+NL(u)-nu*9^n_mod_P|states={checked}|mismatches={mismatches}|pass={}|json=0", u8::from(step_ok));

    let mut rng = SplitMix64(0x0A50_1A21_B000_0002);
    let orbits = [
        Orbit { name: "upstream_energy_test", u0: ShellModel27::new(999_999, 1).u, nu: 1 },
        Orbit { name: "upstream_enstrophy_test", u0: ShellModel27::new(250_000, 3).u, nu: 3 },
        Orbit { name: "zero_state_nu1", u0: [0u64; K], nu: 1 },
        Orbit { name: "room_stub_0x9E3779B97F4A7C15", u0: room_state(0x9E37_79B9_7F4A_7C15), nu: 1 },
        Orbit { name: "ground_room_300g_0x259D7863C53BAE36", u0: room_state(0x259D_7863_C53B_AE36), nu: 1 },
        Orbit { name: "neg_pole_300g_0x56F4EE189BA889DD", u0: room_state(0x56F4_EE18_9BA8_89DD), nu: 2 },
        Orbit { name: "pos_pole_300g_0x4AB3CA1CBD74C442", u0: room_state(0x4AB3_CA1C_BD74_C442), nu: 3 },
        Orbit { name: "splitmix_state_nu_raw", u0: rng.state(), nu: rng.next_u64() },
        Orbit { name: "inviscid_control", u0: ShellModel27::new(999_999, 0).u, nu: 0 },
    ];
    for o in &orbits {
        let (ok, distinct, q_end, first_bad) = run_full_period(o);
        all_ok &= ok;
        println!(
            "ORBIT|name={}|nu={}|steps={P}|q0={}|q_end={q_end}|distinct_q={distinct}|first_violation={first_bad}|pass={}|json=0",
            o.name,
            o.nu,
            shell_total(&o.u0),
            u8::from(ok)
        );
    }

    println!(
        "SCOPE|theorem=exact_all_time_affine_drift_of_shell_total|conserved=Q+t*nu*D|not_claimed=quadratic_energy_conservation,real_norm_bound,continuum_limit,clay_global_regularity|json=0"
    );
    println!("RESULT|pass={}|json=0", u8::from(all_ok));
    if all_ok {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
