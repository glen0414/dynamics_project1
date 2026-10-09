"""Shared parameters, models and tools for the Group 22 pendulum experiments.

State convention follows the course toolkit: x = theta (rad), v = omega (rad/s).
All quantities are per unit mass.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
TOOLKIT_DIR = HERE.parent / "numerical_tools_student_release" / "student_release"
sys.path.insert(0, str(TOOLKIT_DIR))

from mechanics_integrators import _prepare_inputs, rk4, velocity_verlet  # noqa: E402,F401

FIG_DIR = HERE / "figures"
RESULT_DIR = HERE / "results"
FIG_DIR.mkdir(exist_ok=True)
RESULT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------- parameters
G = 9.8                    # m/s^2
L = 1.0                    # m
OMEGA0 = np.sqrt(G / L)    # small-angle angular frequency
T0 = 2 * np.pi / OMEGA0    # small-angle period, about 2.007 s
GAMMA = 0.3                # damping rate (1/s): theta'' = ... - GAMMA * omega
DT = 0.01                  # default time step (s); every experiment also runs DT/2

# ---------------------------------------------------------------- plot style
COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]   # fixed categorical order
REF_COLOR = "#6b6a66"                        # exact / analytic reference (dashed)
plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 150,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#e4e3df",
    "grid.linewidth": 0.6,
    "lines.linewidth": 1.6,
    "legend.frameon": False,
})


# ---------------------------------------------------------------- models
def restoring(theta, model):
    """Angular acceleration from gravity alone."""
    if model == "nonlinear":
        return -OMEGA0**2 * np.sin(theta)
    if model == "linear":
        return -OMEGA0**2 * theta
    raise ValueError(f"unknown model {model!r}")


def accel_rk4(model, gamma=0.0):
    """Acceleration a(t, x, v) for the toolkit's RK4."""
    def a(t, theta, omega):
        return restoring(theta, model) - gamma * omega
    return a


def accel_vv(model):
    """Position-only acceleration a(t, x) for the toolkit's velocity Verlet."""
    def a(t, theta):
        return restoring(theta, model)
    return a


def energy(theta, omega, model):
    """Mechanical energy per unit mass."""
    kinetic = 0.5 * L**2 * omega**2
    if model == "nonlinear":
        return kinetic + G * L * (1 - np.cos(theta))
    return kinetic + 0.5 * G * L * theta**2


def time_grid(t_end, dt):
    n = int(round(t_end / dt))
    return np.linspace(0.0, n * dt, n + 1)


# ---------------------------------------------------------------- exact results
def exact_period_nonlinear(theta0):
    """Exact period of the undamped nonlinear pendulum.

    T = 4 K(k) / omega0 with k = sin(theta0/2), and K(k) = pi / (2 AGM(1, cos(theta0/2))).
    """
    a, b = 1.0, np.cos(theta0 / 2)
    for _ in range(60):
        a, b = 0.5 * (a + b), np.sqrt(a * b)
        if abs(a - b) < 1e-15:
            break
    return 2 * np.pi / (OMEGA0 * a)


def exact_damped_linear(t, theta0, omega_init=0.0, gamma=GAMMA):
    """Analytic solution of theta'' + gamma theta' + omega0^2 theta = 0 (underdamped)."""
    wd = np.sqrt(OMEGA0**2 - gamma**2 / 4)
    b = (omega_init + 0.5 * gamma * theta0) / wd
    decay = np.exp(-0.5 * gamma * t)
    theta = decay * (theta0 * np.cos(wd * t) + b * np.sin(wd * t))
    omega = decay * ((b * wd - 0.5 * gamma * theta0) * np.cos(wd * t)
                     - (theta0 * wd + 0.5 * gamma * b) * np.sin(wd * t))
    return theta, omega


# ---------------------------------------------------------------- periods
def upward_crossings(t, theta):
    """Times where theta crosses zero going upward (linear interpolation)."""
    idx = np.where((theta[:-1] < 0) & (theta[1:] >= 0))[0]
    frac = -theta[idx] / (theta[idx + 1] - theta[idx])
    return t[idx] + frac * (t[idx + 1] - t[idx])


def measure_periods(t, theta):
    """Successive periods (s) and the time at the middle of each period."""
    c = upward_crossings(t, theta)
    return np.diff(c), 0.5 * (c[1:] + c[:-1])


# ---------------------------------------------------------------- damped VV variants
# The toolkit's velocity_verlet only accepts a(t, x). A velocity-dependent force
# therefore needs a modified scheme; these two are the group's own versions.
def vv_damped_lagged(model, gamma, t, x0, v0):
    """Velocity Verlet with damping evaluated at the OLD velocity in a(n+1).

    This is the naive modification: a_{n+1} ~ f(x_{n+1}) - gamma v_n.
    The lag makes the scheme only first-order accurate in dt.
    """
    t, x, v = _prepare_inputs(t, x0, v0)
    dt = t[1] - t[0]
    for n in range(len(t) - 1):
        a_n = restoring(x[n], model) - gamma * v[n]
        x[n + 1] = x[n] + v[n] * dt + 0.5 * a_n * dt**2
        a_next = restoring(x[n + 1], model) - gamma * v[n]
        v[n + 1] = v[n] + 0.5 * (a_n + a_next) * dt
    return x, v


def vv_damped_implicit(model, gamma, t, x0, v0):
    """Velocity Verlet with the damping term at v_{n+1} solved exactly.

    v_{n+1} = v_n + dt/2 [a_n + f(x_{n+1}) - gamma v_{n+1}]
      =>  v_{n+1} = (v_n + dt/2 [a_n + f(x_{n+1})]) / (1 + gamma dt / 2)
    Possible without iteration because the damping is linear in v; second order.
    """
    t, x, v = _prepare_inputs(t, x0, v0)
    dt = t[1] - t[0]
    for n in range(len(t) - 1):
        a_n = restoring(x[n], model) - gamma * v[n]
        x[n + 1] = x[n] + v[n] * dt + 0.5 * a_n * dt**2
        v[n + 1] = (v[n] + 0.5 * dt * (a_n + restoring(x[n + 1], model))) / (1 + 0.5 * gamma * dt)
    return x, v


# ---------------------------------------------------------------- output helpers
class Summary:
    """Collects lines, prints them and writes results/<name>_summary.txt."""

    def __init__(self, name, title):
        self.path = RESULT_DIR / f"{name}_summary.txt"
        self.lines = [title, "=" * len(title)]

    def __call__(self, line=""):
        self.lines.append(line)

    def save(self):
        text = "\n".join(self.lines) + "\n"
        self.path.write_text(text, encoding="utf-8")
        print(text)


def save_fig(fig, name):
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path
