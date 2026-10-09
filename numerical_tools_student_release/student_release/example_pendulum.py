"""Minimal example: nonlinear pendulum without damping."""

import numpy as np
import matplotlib.pyplot as plt

from mechanics_integrators import rk4, velocity_verlet

g = 9.8
L = 1.0


def acceleration(t, theta, omega):
    # This is the nonlinear pendulum equation written as an angular
    # acceleration function for the integrator.
    return -g / L * np.sin(theta)


def position_acceleration(t, theta):
    return -g / L * np.sin(theta)


t = np.linspace(0, 20, 2001)
theta0 = 1.0
omega0 = 0.0

theta_rk4, omega_rk4 = rk4(acceleration, t, theta0, omega0)
theta_vv, omega_vv = velocity_verlet(position_acceleration, t, theta0, omega0)

# The two methods see the same physical model and initial conditions; this
# makes it easy to compare how they evolve the angle over time.
plt.plot(t, theta_rk4, label="RK4")
plt.plot(t, theta_vv, label="velocity Verlet")
plt.xlabel("t")
plt.ylabel(r"$\theta$")
plt.legend()
plt.tight_layout()
plt.show()
