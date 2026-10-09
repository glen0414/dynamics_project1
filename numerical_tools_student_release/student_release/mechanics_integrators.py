"""Numerical integrators for undergraduate Classical Mechanics.

The library provides four methods used in the Mechanics course:
    - semi-implicit Euler
    - Verlet
    - velocity Verlet
    - RK4 (fourth-order Runge-Kutta)

State convention
----------------
The equations of motion are written as

    dx/dt = v
    dv/dt = a(t, x, v)

where x and v may be scalars or NumPy arrays.

Verlet and velocity Verlet are provided in their standard form for
problems whose acceleration is independent of velocity:

    a = a(t, x)

This means the standard Verlet and velocity Verlet interfaces are
intended for position-based acceleration models rather than explicit
velocity-dependent acceleration.
"""

from __future__ import annotations

from typing import Callable
import numpy as np


Array = np.ndarray
Acceleration = Callable[[float, Array, Array], Array]
PositionAcceleration = Callable[[float, Array], Array]


def _validated_acceleration(acceleration_value, state_shape) -> Array:
    """Convert an acceleration value to an array and validate its shape.

    The integrators expect the returned acceleration to match the shape of
    the state they are updating, so broadcasting does not silently change the
    physics of a vector problem.
    """
    acceleration_array = np.asarray(acceleration_value, dtype=float)
    if acceleration_array.shape != state_shape:
        raise ValueError(
            "Acceleration shape is incompatible with the state shape: "
            f"got {acceleration_array.shape}, expected {state_shape}."
        )
    return acceleration_array


def _prepare_inputs(t: Array, x0, v0) -> tuple[Array, Array, Array]:
    """Validate and normalize common inputs."""
    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or len(t) < 3:
        raise ValueError("t must be a one-dimensional array with at least three points.")
    dt = np.diff(t)
    if not np.allclose(dt, dt[0]):
        raise ValueError("The time array must be uniformly spaced.")
    if dt[0] <= 0:
        raise ValueError("The time array must be strictly increasing.")

    x0 = np.asarray(x0, dtype=float)
    v0 = np.asarray(v0, dtype=float)
    if x0.shape != v0.shape:
        raise ValueError("x0 and v0 must have the same shape.")

    x = np.empty((len(t),) + x0.shape, dtype=float)
    v = np.empty_like(x)
    x[0] = x0
    v[0] = v0
    return t, x, v


def semi_implicit_euler(
    acceleration: Acceleration, t: Array, x0, v0
) -> tuple[Array, Array]:
    """Integrate using semi-implicit Euler.

    Parameters
    ----------
    acceleration : callable
        Function ``acceleration(t, x, v)`` returning the acceleration.
    t : array-like
        Uniform time grid.
    x0, v0 : array-like or float
        Initial position and velocity.

    Returns
    -------
    x, v : ndarray
        Position and velocity at each time point.
    """
    t, x, v = _prepare_inputs(t, x0, v0)
    dt = t[1] - t[0]

    for n in range(len(t) - 1):
        a = _validated_acceleration(acceleration(t[n], x[n], v[n]), x[n].shape)
        v[n + 1] = v[n] + a * dt
        x[n + 1] = x[n] + v[n + 1] * dt

    return x, v


def rk4(
    acceleration: Acceleration, t: Array, x0, v0
) -> tuple[Array, Array]:
    """Integrate the equations of motion using fourth-order Runge-Kutta.

    The acceleration may depend on time, position, and velocity.
    """
    t, x, v = _prepare_inputs(t, x0, v0)
    dt = t[1] - t[0]

    def rhs(time, position, velocity):
        a = _validated_acceleration(acceleration(time, position, velocity), position.shape)
        return velocity, a

    for n in range(len(t) - 1):
        tn = t[n]
        xn, vn = x[n], v[n]

        k1x, k1v = rhs(tn, xn, vn)
        k2x, k2v = rhs(tn + dt / 2, xn + dt * k1x / 2, vn + dt * k1v / 2)
        k3x, k3v = rhs(tn + dt / 2, xn + dt * k2x / 2, vn + dt * k2v / 2)
        k4x, k4v = rhs(tn + dt, xn + dt * k3x, vn + dt * k3v)

        x[n + 1] = xn + dt * (k1x + 2 * k2x + 2 * k3x + k4x) / 6
        v[n + 1] = vn + dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6

    return x, v


def verlet(
    acceleration: PositionAcceleration, t: Array, x0, v0
) -> tuple[Array, Array]:
    """Integrate using the standard position-Verlet method.

    ``acceleration`` must have the form ``acceleration(t, x)`` and must
    not depend on velocity. The initial velocity is used to construct
    the first position step. Velocities at later times are reconstructed
    from neighboring positions.
    """
    t, x, v = _prepare_inputs(t, x0, v0)
    dt = t[1] - t[0]

    a0 = _validated_acceleration(acceleration(t[0], x[0]), x[0].shape)
    x[1] = x[0] + v[0] * dt + 0.5 * a0 * dt**2

    for n in range(1, len(t) - 1):
        a = _validated_acceleration(acceleration(t[n], x[n]), x[n].shape)
        x[n + 1] = 2 * x[n] - x[n - 1] + a * dt**2

    # Centered difference for velocity; one-sided second-order estimates
    # at the endpoints.
    v[0] = v0
    v[-1] = (3 * x[-1] - 4 * x[-2] + x[-3]) / (2 * dt)
    if len(t) > 2:
        v[1:-1] = (x[2:] - x[:-2]) / (2 * dt)

    return x, v


def velocity_verlet(
    acceleration: PositionAcceleration, t: Array, x0, v0
) -> tuple[Array, Array]:
    """Integrate using the standard velocity-Verlet method.

    ``acceleration`` must have the form ``acceleration(t, x)`` and must
    not depend on velocity.
    """
    t, x, v = _prepare_inputs(t, x0, v0)
    dt = t[1] - t[0]

    v[0] = v0
    for n in range(len(t) - 1):
        a_n = _validated_acceleration(acceleration(t[n], x[n]), x[n].shape)
        x[n + 1] = x[n] + v[n] * dt + 0.5 * a_n * dt**2
        a_next = _validated_acceleration(acceleration(t[n + 1], x[n + 1]), x[n + 1].shape)
        v[n + 1] = v[n] + 0.5 * (a_n + a_next) * dt

    return x, v


METHODS = {
    "semi-implicit Euler": semi_implicit_euler,
    "Verlet": verlet,
    "velocity Verlet": velocity_verlet,
    "RK4": rk4,
}


def integrate(method: str, acceleration, t: Array, x0, v0):
    """Convenience wrapper for selecting one of the four methods.

    For semi-implicit Euler and RK4, ``acceleration`` must have the
    signature ``acceleration(t, x, v)``.

    For Verlet and velocity Verlet, ``acceleration`` must have the
    signature ``acceleration(t, x)`` and must be independent of velocity.
    """
    key = method.strip()
    if key not in METHODS:
        available = ", ".join(METHODS)
        raise ValueError(f"Unknown method {method!r}. Available methods: {available}")
    return METHODS[key](acceleration, t, x0, v0)
