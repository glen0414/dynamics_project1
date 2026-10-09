"""Experiment 2: numerical energy error without damping, VV vs RK4.

Nonlinear pendulum, no damping, large amplitude. The energy is exactly
conserved by the physics, so any change in E is numerical.
"""

import matplotlib.pyplot as plt
import numpy as np

from common import (COLORS, DT, Summary, accel_rk4, accel_vv, energy,
                    exact_period_nonlinear, measure_periods, rk4, save_fig,
                    time_grid, velocity_verlet)

THETA0 = np.radians(120)
T_END = 2000.0      # long run so the RK4 drift is measurable at DT
DT_COARSE = 0.1     # deliberately coarse step to make the qualitative behavior visible
T_END_COARSE = 2000.0


def run_both(dt, t_end):
    t = time_grid(t_end, dt)
    th_rk, om_rk = rk4(accel_rk4("nonlinear"), t, THETA0, 0.0)
    th_vv, om_vv = velocity_verlet(accel_vv("nonlinear"), t, THETA0, 0.0)
    return t, {"RK4": (th_rk, om_rk), "velocity Verlet": (th_vv, om_vv)}


def rel_energy_error(theta, omega):
    e0 = energy(THETA0, 0.0, "nonlinear")
    return energy(theta, omega, "nonlinear") / e0 - 1


def main():
    s = Summary("exp2", "Experiment 2 - energy error without damping (VV vs RK4)")
    t_exact = exact_period_nonlinear(THETA0)
    s(f"theta0 = 120 deg, no damping, t_end = {T_END:.0f} s, exact period = {t_exact:.6f} s")
    s()

    stats = {}
    runs = {}
    for dt in (DT, DT / 2):
        t, res = run_both(dt, T_END)
        runs[dt] = (t, res)
        for name, (th, om) in res.items():
            de = rel_energy_error(th, om)
            # RK4 drift rate: least-squares slope of dE/E0 vs t
            slope = np.polyfit(t, de, 1)[0]
            period = measure_periods(t, th)[0].mean()
            stats[(name, dt)] = dict(max_abs=np.abs(de).max(), end=de[-1], slope=slope,
                                     period_err=period - t_exact)

    s(f"{'method':>16} {'dt':>7} {'max|dE/E0|':>11} {'dE/E0 at end':>13} "
      f"{'drift rate [1/s]':>17} {'period error [s]':>17}")
    for (name, dt), d in stats.items():
        s(f"{name:>16} {dt:>7.4f} {d['max_abs']:>11.2e} {d['end']:>13.2e} "
          f"{d['slope']:>17.2e} {d['period_err']:>17.2e}")
    s()
    s("Ratios when dt is halved (A1 section 4 check):")
    for name in ("velocity Verlet", "RK4"):
        a, b = stats[(name, DT)], stats[(name, DT / 2)]
        s(f"  {name:>16}: max|dE| ratio = {a['max_abs'] / b['max_abs']:6.1f}, "
          f"period-error ratio = {a['period_err'] / b['period_err']:6.1f}")
    rk_rate = -stats[("RK4", DT)]["slope"]
    vv_bound = stats[("velocity Verlet", DT)]["max_abs"]
    s()
    s(f"Crossover estimate at dt = {DT}: RK4 drift reaches the VV error band after about")
    s(f"  t* = {vv_bound:.2e} / {rk_rate:.2e} = {vv_bound / rk_rate:.2e} s "
      f"(~{vv_bound / rk_rate / t_exact:.1e} periods)")

    # ---- figure 1: per-period envelope of |dE/E0| vs time, DT and DT/2 (log scale)
    # VV's error oscillates within every swing; the envelope (max over one exact
    # period) shows the long-term trend of both methods on the same plot.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, dt in zip(axes, (DT, DT / 2)):
        t, res = runs[dt]
        win = int(round(t_exact / dt))
        n_win = len(t) // win
        t_mid = t[: n_win * win].reshape(n_win, win).mean(axis=1)
        for color, z, (name, (th, om)) in zip(COLORS, (3, 2), res.items()):
            de = np.abs(rel_energy_error(th, om))[: n_win * win].reshape(n_win, win).max(axis=1)
            ax.plot(t_mid, de, color=color, lw=1.6, zorder=z, label=name)
        ax.set_yscale("log")
        ax.set_title(f"dt = {dt:g} s", loc="left", fontsize=10)
        ax.set_xlabel("t (s)")
    axes[0].set_ylabel(r"max $|\Delta E / E_0|$ over each period")
    axes[0].legend(loc="center right")
    fig.suptitle(r"Exp 2: energy error, nonlinear pendulum, $\theta_0=120^\circ$, no damping")
    save_fig(fig, "exp2_energy_error")

    # ---- figure 2: coarse step, signed energy error + phase space
    t, res = run_both(DT_COARSE, T_END_COARSE)
    fig, (ax_e, ax_p) = plt.subplots(1, 2, figsize=(10, 4.2))
    for color, (name, (th, om)) in zip(COLORS, res.items()):
        de = rel_energy_error(th, om)
        ax_e.plot(t, de, color=color, lw=0.8, label=name)
        ax_p.plot(np.degrees(th), om, color=color, lw=0.4, alpha=0.8, label=name)
        stats[(name, DT_COARSE)] = de
    ax_e.axhline(0, color="#6b6a66", lw=0.8, ls="--")
    ax_e.set_xlabel("t (s)")
    ax_e.set_ylabel(r"$\Delta E / E_0$")
    ax_e.set_title(f"energy error, dt = {DT_COARSE:g} s", loc="left", fontsize=10)
    ax_e.legend(loc="lower left")
    ax_p.set_xlabel(r"$\theta$ (deg)")
    ax_p.set_ylabel(r"$\omega$ (rad/s)")
    ax_p.set_title(f"phase space, dt = {DT_COARSE:g} s, t = 0-{T_END_COARSE:.0f} s",
                   loc="left", fontsize=10)
    fig.suptitle("Exp 2 (coarse step): RK4 spirals inward, VV stays on a closed band")
    save_fig(fig, "exp2_phase_space_coarse")

    s()
    s(f"Coarse run dt = {DT_COARSE}: dE/E0 at t = {T_END_COARSE:.0f} s: "
      f"RK4 {stats[('RK4', DT_COARSE)][-1]:.2e}, "
      f"VV {stats[('velocity Verlet', DT_COARSE)][-1]:.2e} "
      f"(VV max |dE/E0| {np.abs(stats[('velocity Verlet', DT_COARSE)]).max():.2e})")
    s()
    s("Reading:")
    s("- VV: energy error oscillates in a bounded band (no secular drift); band shrinks ~4x per dt halving.")
    s("- RK4: energy decreases monotonically (numerical dissipation), linear in t.")
    s("  Its drift rate shrinks ~32x per dt halving (h^5, faster than the h^4 assumed in A1).")
    s("- At dt = 0.01 the RK4 drift stays far below the VV band for any practical run time;")
    s("  the inward phase-space spiral predicted in A1 is only visible with a coarse step.")
    s.save()


if __name__ == "__main__":
    main()
