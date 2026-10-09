"""Experiment 3: accuracy with linear damping, against the analytic solution.

Small-angle equation with damping gamma * omega. Includes the finding that the
toolkit's velocity Verlet cannot take a velocity-dependent force, and compares
two modified VV schemes with RK4.
"""

import matplotlib.pyplot as plt
import numpy as np

from common import (COLORS, DT, GAMMA, T0, Summary, accel_rk4, exact_damped_linear,
                    rk4, save_fig, time_grid, velocity_verlet, vv_damped_implicit,
                    vv_damped_lagged)

THETA0 = np.radians(10)
T_END = 30.0
CONV_DTS = [0.04, 0.02, 0.01, 0.005, 0.0025]

METHODS = {
    "RK4 (toolkit)": lambda t: rk4(accel_rk4("linear", GAMMA), t, THETA0, 0.0),
    "VV, implicit damping (ours)": lambda t: vv_damped_implicit("linear", GAMMA, t, THETA0, 0.0),
    "VV, lagged damping (ours)": lambda t: vv_damped_lagged("linear", GAMMA, t, THETA0, 0.0),
}


def toolkit_limit_demo():
    """Try to pass the damped acceleration to the toolkit's velocity Verlet."""
    t = time_grid(1.0, DT)
    try:
        velocity_verlet(accel_rk4("linear", GAMMA), t, THETA0, 0.0)
    except TypeError as err:
        return f"TypeError: {err}"
    return "no error (unexpected)"


def max_error(dt):
    t = time_grid(T_END, dt)
    exact, _ = exact_damped_linear(t, THETA0)
    return {name: np.abs(f(t)[0] - exact).max() for name, f in METHODS.items()}


def main():
    s = Summary("exp3", "Experiment 3 - damped small-angle pendulum vs analytic solution")
    s(f"theta0 = 10 deg, gamma = {GAMMA} 1/s, t_end = {T_END:.0f} s")
    s()
    s("Finding (toolkit limit): passing a(t, theta, omega) with damping to the toolkit's")
    s("velocity_verlet, which calls acceleration(t, x), gives:")
    s(f"  {toolkit_limit_demo()}")
    s("The standard VV interface cannot represent -gamma*omega, so VV must be modified.")
    s("Two modifications are compared (see common.py): lagged and implicit damping.")
    s()

    # ---- figure 1: error vs time at DT and DT/2
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    errs = {}
    for ax, dt in zip(axes, (DT, DT / 2)):
        t = time_grid(T_END, dt)
        exact, _ = exact_damped_linear(t, THETA0)
        for color, (name, f) in zip(COLORS, METHODS.items()):
            err = np.abs(f(t)[0] - exact)
            errs[(name, dt)] = err.max()
            ax.plot(t, np.maximum(err, 1e-16), color=color, lw=0.8, label=name)
        ax.set_yscale("log")
        ax.set_xlabel("t (s)")
        ax.set_title(f"dt = {dt:g} s", loc="left", fontsize=10)
    axes[0].set_ylabel(r"$|\theta_{num} - \theta_{exact}|$ (rad)")
    axes[0].legend(loc="lower right", fontsize=8)
    fig.suptitle(rf"Exp 3: error vs analytic solution, small-angle, $\gamma = {GAMMA}$ 1/s")
    save_fig(fig, "exp3_error_vs_time")

    s(f"{'method':>28} {'max err (dt)':>13} {'max err (dt/2)':>15} {'ratio':>7}")
    for name in METHODS:
        a, b = errs[(name, DT)], errs[(name, DT / 2)]
        s(f"{name:>28} {a:>13.2e} {b:>15.2e} {a / b:>7.1f}")
    s()

    # ---- figure 2: convergence
    conv = {name: [] for name in METHODS}
    for dt in CONV_DTS:
        for name, e in max_error(dt).items():
            conv[name].append(e)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    s("Convergence (max error over 0-30 s) and fitted order p (error ~ dt^p):")
    for color, (name, e) in zip(COLORS, conv.items()):
        p = np.polyfit(np.log(CONV_DTS), np.log(e), 1)[0]
        ax.loglog(CONV_DTS, e, "o-", color=color, ms=6, label=f"{name}, p = {p:.2f}")
        s(f"  {name:>28}: p = {p:.2f}   errors = " + ", ".join(f"{x:.1e}" for x in e))
    ax.set_xlabel("dt (s)")
    ax.set_ylabel("max error (rad)")
    ax.legend(fontsize=8)
    ax.set_title("Exp 3: convergence with dt (damped, small-angle)", loc="left", fontsize=10)
    save_fig(fig, "exp3_convergence")

    # ---- figure 3: trajectories, for context
    t = time_grid(T_END, DT)
    exact, _ = exact_damped_linear(t, THETA0)
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.plot(t, np.degrees(METHODS["RK4 (toolkit)"](t)[0]), color=COLORS[0], label="RK4")
    ax.plot(t, np.degrees(exact), "--", color="#6b6a66", label="analytic")
    ax.set_xlabel("t (s)")
    ax.set_ylabel(r"$\theta$ (deg)")
    ax.legend(loc="upper right")
    ax.set_title("Exp 3: damped small-angle pendulum", loc="left", fontsize=10)
    save_fig(fig, "exp3_theta_t")

    s()
    s(f"Period scale for reference: T0 = {T0:.3f} s, decay time 2/gamma = {2 / GAMMA:.2f} s")
    s("Reading:")
    s("- Lagged-damping VV is first order (error halves with dt): this is the A1-predicted")
    s("  'velocity lag' error, and it dominates.")
    s("- Implicit-damping VV recovers second order (~4x per halving).")
    s("- RK4 handles a(t, x, v) natively and is fourth order (~16x per halving).")
    s.save()


if __name__ == "__main__":
    main()
