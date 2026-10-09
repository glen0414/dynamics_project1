"""Minimal example: projectile motion with linear drag."""

import numpy as np
import matplotlib.pyplot as plt

from mechanics_integrators import semi_implicit_euler, rk4

m = 1.0
b = 0.15
g = 9.8


def acceleration(t, x, v):
    # Drag is written as a force proportional to velocity, so the returned
    # acceleration depends on the current state as well as gravity.
    # F_drag = -b v
    return np.array([-b * v[0] / m, -g - b * v[1] / m])


t = np.linspace(0, 5, 1001)
x0 = np.array([0.0, 0.0])
v0 = np.array([10.0, 15.0])

x_euler, v_euler = semi_implicit_euler(acceleration, t, x0, v0)
x_rk4, v_rk4 = rk4(acceleration, t, x0, v0)

# Both methods use the same physical model and time grid, so the plotted
# difference comes from the numerical integrator rather than from the setup.
plt.plot(x_euler[:, 0], x_euler[:, 1], label="semi-implicit Euler")
plt.plot(x_rk4[:, 0], x_rk4[:, 1], label="RK4")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.tight_layout()
plt.show()
