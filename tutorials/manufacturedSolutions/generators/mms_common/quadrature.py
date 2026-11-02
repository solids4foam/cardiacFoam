"""Gauss-Legendre quadratures for the manufactured-solution checks.

Integrands take the points as one array of shape (dim, ...) and return an
array of shape (...).
"""

import numpy as np


def _nodes(a, b, n):
    """Gauss-Legendre nodes and weights on [a, b]."""
    t, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * t + 0.5 * (b + a), 0.5 * (b - a) * w


def gauss_legendre_box(f, lo, hi, n):
    """Tensor-product Gauss-Legendre integral of f over the box [lo, hi]."""
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    axes, weights = zip(*(_nodes(a, b, n) for a, b in zip(lo, hi)))
    grid = np.meshgrid(*axes, indexing="ij")
    wgrid = np.ones_like(grid[0])
    for i, w in enumerate(weights):
        shape = [1] * len(weights)
        shape[i] = n
        wgrid = wgrid * w.reshape(shape)
    return float(np.sum(wgrid * f(np.array(grid))))


def _ray_length(direction, centre, radius, box_lo, box_hi):
    """Length of each ray from the centre, cut by the ball and the box."""
    with np.errstate(divide="ignore", invalid="ignore"):
        to_hi = np.where(direction > 0, (box_hi[:, None] - centre[:, None]) / direction, np.inf)
        to_lo = np.where(direction < 0, (box_lo[:, None] - centre[:, None]) / direction, np.inf)
    return np.minimum(radius, np.minimum(to_hi, to_lo).min(axis=0))


def half_ball_kernel_integral(f, centre, radius, box_lo, box_hi, n):
    """Integrals of w*f and w over the ball intersected with the box.

    w = max(0, 1 - |x - centre| / radius).  Spherical coordinates about the
    centre with the polar axis along the first coordinate: the polar cosine
    mu is split where the first-axis faces start to cut the ball, so a ball
    centred on a first-axis face (a half ball) is integrated with smooth
    integrands on every piece; n points per direction and per piece.  Cuts by
    the other faces enter through the ray length and converge algebraically.
    """
    centre = np.asarray(centre, dtype=float)
    box_lo = np.asarray(box_lo, dtype=float)
    box_hi = np.asarray(box_hi, dtype=float)
    breaks = {-1.0, 1.0}
    for d in (centre[0] - box_lo[0], box_hi[0] - centre[0]):
        if 0.0 <= d < radius:
            breaks.update((d / radius, -d / radius))
    breaks = sorted(breaks)

    phi, w_phi = _nodes(0.0, 2.0 * np.pi, n)
    t_r, w_r = np.polynomial.legendre.leggauss(n)
    total_f = 0.0
    total_w = 0.0
    for a, b in zip(breaks[:-1], breaks[1:]):
        if b - a < 1e-14:
            continue
        mu, w_mu = _nodes(a, b, n)
        MU, PHI = np.meshgrid(mu, phi, indexing="ij")
        W2 = w_mu[:, None] * w_phi[None, :]
        sin_t = np.sqrt(np.maximum(0.0, 1.0 - MU**2))
        direction = np.array([MU, sin_t * np.cos(PHI), sin_t * np.sin(PHI)]).reshape(3, -1)
        length = _ray_length(direction, centre, radius, box_lo, box_hi)
        r = 0.5 * length[None, :] * (t_r[:, None] + 1.0)
        wr = 0.5 * length[None, :] * w_r[:, None]
        pts = centre[:, None, None] + direction[:, None, :] * r[None]
        kernel = np.maximum(0.0, 1.0 - r / radius)
        weight = (W2.reshape(1, -1) * wr * r**2) * kernel
        total_f += float(np.sum(weight * f(pts)))
        total_w += float(np.sum(weight))
    return total_f, total_w
