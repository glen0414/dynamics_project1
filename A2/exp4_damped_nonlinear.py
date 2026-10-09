"""Experiment 4: small-angle vs nonlinear equation with damping (RK4).

Large initial angle with damping. As the amplitude decays, the nonlinear
period should approach the damped small-angle period.
"""

import matplotlib.pyplot as plt
import numpy as np

from common import (COLORS, DT, GAMMA, OMEGA0, Summary, accel_rk4, measure_periods,
                    rk4, save_fig, time_grid)

THETA0 = np.radians(120)
T_END = 30.0


def run(dt, model):
    t = time_grid(T_END, dt)
    theta, _ = rk4(accel_rk4(model, GAMMA), t, THETA0, 0.0)
    return t, theta


def main():
    s = Summary("exp4", "Experiment 4 - damped pendulum, small-angle vs nonlinear (RK4)")
    td = 2 * np.pi / np.sqrt(OMEGA0**2 - GAMMA**2 / 4)
    s(f"theta0 = 120 deg, gamma = {GAMMA} 1/s, damped small-angle period Td = {td:.6f} s")
    s()

    runs = {m: {dt: run(dt, m) for dt in (DT, DT / 2)} for m in ("linear", "nonlinear")}

    fig, (ax_th, ax_T) = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
    for color, model in zip(COLORS, ("linear", "nonlinear")):
        t, th = runs[model][DT]
        label = "small-angle (linear)" if model == "linear" else "nonlinear"
        ax_th.plot(t, np.degrees(th), color=color, label=label)
        periods, mid = measure_periods(t, th)
        ax_T.plot(mid, periods, "o-", color=color, ms=4, label=label)
        # dt-halving check on the trajectory itself
        t2, th2 = runs[model][DT / 2]
        diff = np.abs(th - th2[::2]).max()
        s(f"{model:>9}: max |theta(dt) - theta(dt/2)| = {diff:.2e} rad")
        s(f"{'':>9}  successive periods [s]: " + ", ".join(f"{p:.4f}" for p in periods[:6])
          + " ... " + ", ".join(f"{p:.4f}" for p in periods[-3:]))
    ax_T.axhline(td, color="#6b6a66", ls="--", lw=1, label=r"damped small-angle $T_d$")
    ax_th.set_ylabel(r"$\theta$ (deg)")
    ax_th.legend(loc="upper right")
    ax_T.set_ylabel("period of each swing (s)")
    ax_T.set_xlabel("t (s)")
    ax_T.legend(loc="upper right")
    fig.suptitle(rf"Exp 4: damped pendulum, $\theta_0 = 120^\circ$, $\gamma = {GAMMA}$ 1/s (RK4)")
    save_fig(fig, "exp4_theta_t_periods")

    # phase offset between the two models after the amplitude has decayed
    c_lin = measure_periods(*runs["linear"][DT])[1]
    c_non = measure_periods(*runs["nonlinear"][DT])[1]
    n = min(len(c_lin), len(c_non))
    lag = (c_non[:n] - c_lin[:n])[-1]
    s()
    s(f"Time lag of nonlinear behind linear at the last common swing: {lag:.3f} s "
      f"({lag / td * 360:.0f} deg of phase)")
    s()
    s("Reading:")
    s("- Nonlinear period starts longer than Td and converges to Td as the amplitude decays")
    s("  (the frequencies synchronise) - this is physics: unchanged when dt is halved.")
    s("- The phase lag built up during the large-amplitude swings does NOT disappear,")
    s("  so the two curves keep a constant offset rather than coinciding exactly.")
    s.save()


if __name__ == "__main__":
    main()
