"""Getting started example for the mechanics numerical toolkit.

This example shows a one-dimensional free-fall problem with constant
acceleration:

    x'' = -g

It is intentionally simple and is meant only as an onboarding example.
"""

import numpy as np
import matplotlib.pyplot as plt

from mechanics_integrators import rk4


g = 9.8


def acceleration(time, x, v):
    # The integrator only needs the acceleration; the free-fall model is kept
    # in this function so the physical assumptions are explicit.
    return -g


t = np.linspace(0, 2.0, 201)
x0 = 1.5
v0 = 4.0

x, v = rk4(acceleration, t, x0, v0)

# The analytic solution is a simple reference check for this onboarding case.
x_exact = x0 + v0 * t - 0.5 * g * t**2
v_exact = v0 - g * t

plt.figure(figsize=(7, 4))
plt.plot(t, x, label="RK4 numerical x(t)")
plt.plot(t, x_exact, "--", label="analytic x(t)")
plt.xlabel("t")
plt.ylabel("x")
plt.title("Free fall with constant acceleration")
plt.legend()
plt.tight_layout()
plt.show()

# The printed values show that the numerical result matches the reference.
print(f"final numerical position: {x[-1]:.6f}")
print(f"final analytic position:  {x_exact[-1]:.6f}")
print(f"final numerical velocity: {v[-1]:.6f}")
print(f"final analytic velocity:  {v_exact[-1]:.6f}")