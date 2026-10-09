"""Experiment 1: physical effect of amplitude on the period (no damping, RK4).

Compare the small-angle (linear) and full nonlinear equations at theta0 = 5 deg
and 60 deg. Run at DT and DT/2 to confirm the difference is physical.
"""

import matplotlib.pyplot as plt
import numpy as np

from common import (COLORS, DT, T0, Summary, accel_rk4, exact_period_nonlinear,
                    measure_periods, rk4, save_fig, time_grid)

ANGLES_DEG = [5, 60]
N_PERIODS = 10


def run(dt, theta0, model):
    t = time_grid(N_PERIODS * T0 * 1.4, dt)
    theta, omega = rk4(accel_rk4(model), t, theta0, 0.0)
    return t, theta


def main():
    s = Summary("exp1", "Experiment 1 - amplitude vs period (no damping, RK4)")
    s(f"small-angle period T0 = 2 pi sqrt(L/g) = {T0:.6f} s")
    s()
    s(f"{'theta0':>7} {'model':>9} {'T(dt) [s]':>12} {'T(dt/2) [s]':>12} "
      f"{'|change|':>10} {'T exact [s]':>12} {'T/T0':>8}")

    fig, axes = plt.subplots(len(ANGLES_DEG), 1, figsize=(8, 5.5), sharex=True)
    for ax, deg in zip(axes, ANGLES_DEG):
        theta0 = np.radians(deg)
        for color, model in zip(COLORS, ["linear", "nonlinear"]):
            periods = {}
            for step in (DT, DT / 2):
                t, theta = run(step, theta0, model)
                periods[step] = measure_periods(t, theta)[0].mean()
                if step == DT:
                    label = "small-angle (linear)" if model == "linear" else "nonlinear"
                    ax.plot(t, np.degrees(theta), color=color, label=label)
            exact = T0 if model == "linear" else exact_period_nonlinear(theta0)
            s(f"{deg:>6}° {model:>9} {periods[DT]:>12.6f} {periods[DT / 2]:>12.6f} "
              f"{abs(periods[DT] - periods[DT / 2]):>10.1e} {exact:>12.6f} "
              f"{periods[DT] / T0:>8.4f}")
        ax.set_ylabel(r"$\theta$ (deg)")
        ax.set_title(rf"$\theta_0 = {deg}^\circ$", loc="left", fontsize=10)
        ax.set_xlim(0, N_PERIODS * T0)
    axes[0].legend(loc="upper right", ncol=2)
    axes[-1].set_xlabel("t (s)")
    fig.suptitle("Exp 1: small-angle vs nonlinear pendulum (RK4, no damping)")
    save_fig(fig, "exp1_theta_t")

    s()
    s("Reading: the nonlinear/linear period difference is set by theta0 (physics);")
    s("halving dt changes the measured period by orders of magnitude less than that")
    s("difference, and RK4 matches the exact AGM/elliptic-integral period.")
    s.save()


if __name__ == "__main__":
    main()
