"""Low-level geometric and circular-angle helpers."""

from __future__ import annotations

import numpy as np

TAU = 2.0 * np.pi


def wrap_angle(theta: np.ndarray | float) -> np.ndarray | float:
    """Wrap angle(s) into [0, 2pi)."""
    return np.mod(theta, TAU)


def cw_unwrap(a: float, b: float) -> float:
    """CW arc length from a to b in [0, 2pi)."""
    d = b - a
    if d < 0:
        d += TAU
    return d


def cw_is_between(theta: np.ndarray, start: float, end: float) -> np.ndarray:
    """True where theta lies on CW arc [start, end), with wrapping."""
    t = wrap_angle(theta)
    s = wrap_angle(start)
    e = wrap_angle(end)
    if s <= e:
        return (t >= s) & (t < e)
    return (t >= s) | (t < e)


def circular_mean(theta: np.ndarray) -> float:
    """Circular mean angle in [0, 2pi)."""
    s = np.sin(theta).mean()
    c = np.cos(theta).mean()
    return float(wrap_angle(np.arctan2(s, c)))


def inertia_tensor(points: np.ndarray, center: np.ndarray | None = None) -> np.ndarray:
    """Unweighted inertia tensor around center (3x3)."""
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 3:
        raise ValueError("points must have shape (N,3)")
    if center is None:
        center = p.mean(axis=0)
    r = p - center
    x, y, z = r[:, 0], r[:, 1], r[:, 2]
    ixx = np.sum(y * y + z * z)
    iyy = np.sum(x * x + z * z)
    izz = np.sum(x * x + y * y)
    ixy = -np.sum(x * y)
    ixz = -np.sum(x * z)
    iyz = -np.sum(y * z)
    return np.array(
        [[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]],
        dtype=float,
    )
