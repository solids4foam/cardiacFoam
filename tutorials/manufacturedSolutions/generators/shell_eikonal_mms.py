#!/usr/bin/env python3
"""Heart-only eikonal manufactured solution on the spherical shell.

The heart is the shell r1 < |x| < r2; both spheres are sealed walls.  The
eikonal operator of the solver is

    -div(M grad psi) + c0 sqrt(grad psi . M grad psi) = 1 + S,   M = G/(chi cm),

with G the full conductivity of shell_monodomain_mms (three distinct
eigenvalues, fibre and sheet rotating through the wall).  psi = hhat + omega*m
with hhat = sqrt((x - xs)^T Ms^-1 (x - xs) + delta^2)/c0, the smoothed
anisotropic distance from the stimulus point xs (Ms = M(xs)), omega zero on
both spheres and m the blend of the per-sphere corrections, satisfies
n.M grad(psi) = 0 exactly on each sphere.  S is the exact residual of the
operator, including the gradient of M.  The sealed problem fixes psi up to a
constant, so the verifier pins the cells within rPin of xs to the exact psi.

psi scales as 1/c0 (hhat and m both carry 1/c0), so the critical points and
their types do not depend on c0.  --c0 VALUE runs --report and --self-check at
that c0; the header (--write, --check) and the emitted reference carry c0 as a
runtime parameter and are written at the default c0 = 1.

--emit writes reference/shell_eikonal.py (psi, gradPsi, source and sigma as
NumPy functions of x, y, z and the parameters, with PARAMS and field()) and the
LaTeX building blocks latex/shell_eikonal_*.tex: hhat, M, x_s, m_j, m, psi and
the source operator.  The geometry fragments (lambda, omega, fibre frame,
sigma, wall normals, d_j) are those of shell_monodomain_mms.py.  The expanded
source is not printed; it is source of the reference module.

psi has critical points on the walls (where the wave collides behind the
cavity); there the metric sqrt(grad psi . M grad psi) is not differentiable,
which is the feature of the solution.

    python shell_eikonal_mms.py --write ../../../src/verificationModels/eikonalVerification/manufacturedShellEikonalReference.H
    python shell_eikonal_mms.py --check ../../../src/verificationModels/eikonalVerification/manufacturedShellEikonalReference.H
    python shell_eikonal_mms.py --report [--c0 VALUE]
    python shell_eikonal_mms.py --self-check [--c0 VALUE]
    python shell_eikonal_mms.py --emit
"""

import argparse
import functools
import math
import os
import sys

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.optimize import least_squares

import mms_codegen as codegen
import shell_monodomain_mms as smm
from mms_common.latex_emit import emit_latex_fragment

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "shell_eikonal"
HEADER_NAME = "manufacturedShellEikonalReference.H"
HEADER_PATH = "../../../src/verificationModels/eikonalVerification/" + HEADER_NAME
SAMPLE_DPS = 30

SCALAR_PARAMS = ("r1", "r2", "gamma1", "gamma2", "gil", "git", "gin", "theta0", "theta1",
                 "chi", "cm", "c0", "delta", "rPin")
TAIL_PARAMS = ()
ARRAY_PARAMS = (("xs", 3), ("Msi", 6))
ANGLE_PARAMS = ("gamma1", "gamma2", "theta0", "theta1")
PARAM_ORDER = list(SCALAR_PARAMS) + ["xs%d" % i for i in range(3)] + ["Msi%d" % i for i in range(6)]
P = {name: sp.Symbol(name, real=True) for name in PARAM_ORDER}

x, y, z = sp.symbols("x y z", real=True)
ARGS = [x, y, z] + [P[n] for n in PARAM_ORDER]
XYZ = sp.Matrix([x, y, z])

SOLUTIONS = ("shell",)
HEART_ZONE = 1
ZONE_NAMES = {1: "myocardium"}


def _stimulus(params, direction, radius):
    """xs on the ray along direction at the given radius, and the six components
    (xx xy xz yy yz zz) of M(xs)^-1, inverted at 50 digits."""
    norm = math.sqrt(math.fsum(c*c for c in direction))
    params["xs"] = [radius*c/norm for c in direction]
    params["Msi"] = [0.0]*6
    with mp.workdps(50):
        G6 = evaluate("G", params["xs"], params, "mpmath")
        scale = mp.mpf(params["chi"])*mp.mpf(params["cm"])
        Ms = mp.matrix([[G6[0], G6[1], G6[2]], [G6[1], G6[3], G6[4]],
                        [G6[2], G6[4], G6[5]]])/scale
        inv = Ms**-1
        params["Msi"] = [float(inv[i, j]) for i, j in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))]
    return params


STIMULUS = {"default": ((1.0, 0.3, 0.2), 0.75), "alternate": ((0.9, 0.35, 0.25), 0.7)}


def default_params():
    p = {"r1": 0.5, "r2": 1.0, "gamma1": 20.0, "gamma2": 50.0,
         "gil": 0.17, "git": 0.019, "gin": 0.01, "theta0": -20.0, "theta1": 90.0,
         "chi": 3.0, "cm": 2.0, "c0": 1.0, "delta": 0.05, "rPin": 0.15}
    return _stimulus(p, *STIMULUS["default"])


def alternate_params():
    p = default_params()
    p.update({"theta1": 70.0, "delta": 0.06, "rPin": 0.12})
    return _stimulus(p, *STIMULUS["alternate"])


def param_values(params):
    """Flat {symbol name: float}; angles in radians as Foam::degToRad does."""
    v = {}
    for name in SCALAR_PARAMS:
        val = float(params[name])
        v[name] = val*math.pi/180.0 if name in ANGLE_PARAMS else val
    for name, n in ARRAY_PARAMS:
        assert len(params[name]) == n, name
        for i in range(n):
            v["%s%d" % (name, i)] = float(params[name][i])
    return v


def _grad(u):
    return sp.Matrix([sp.diff(u, v) for v in (x, y, z)])


def _div(J):
    return sum(sp.diff(J[i], v) for i, v in enumerate((x, y, z)))


def _dot(a, b):
    return (a.T*b)[0]


def _components(T):
    return [T[0, 0], T[0, 1], T[0, 2], T[1, 1], T[1, 2], T[2, 2]]


def _metric_tensor(G):
    return G/(P["chi"]*P["cm"])


def _hhat(quadratic):
    """sqrt(q + delta^2)/c0 for the quadratic form q = (x - xs).Ms^-1 (x - xs)."""
    return sp.sqrt(quadratic + P["delta"]**2)/P["c0"]


def _operator(div_flux, quadratic):
    """-div(M grad psi) + c0 sqrt(grad psi.M grad psi) - 1."""
    return -div_flux + P["c0"]*sp.sqrt(quadratic) - 1


@functools.lru_cache(maxsize=None)
def expressions():
    r1, r2 = P["r1"], P["r2"]
    lam, om = smm._radial(x**2 + y**2 + z**2, r1, r2)
    G = smm._conductivity(P["theta0"] + P["theta1"]*lam)
    M = _metric_tensor(G)
    Msi = sp.Matrix([[P["Msi0"], P["Msi1"], P["Msi2"]],
                     [P["Msi1"], P["Msi3"], P["Msi4"]],
                     [P["Msi2"], P["Msi4"], P["Msi5"]]])
    dvec = XYZ - sp.Matrix([P["xs0"], P["xs1"], P["xs2"]])
    hhat = _hhat(_dot(dvec, Msi*dvec))
    n, d = smm._walls(XYZ, r1, r2)
    mj = {j: smm._wall_correction(_dot(n[j], M*_grad(hhat)), _dot(n[j], M*n[j]), d[j])
          for j in (1, 2)}
    m = smm._blend(lam, mj[1], mj[2])
    psi = hhat + om*m
    gradPsi = _grad(psi)
    flux = M*gradPsi
    source = _operator(_div(flux), _dot(gradPsi, flux))
    return {"psi": psi, "gradPsi": list(gradPsi), "source": source,
            "G": _components(G), "Gmat": G, "hhat": hhat,
            "denominators": [_dot(n[j], G*n[j]) for j in (1, 2)]}


def _expr(solution, zone, name):
    return expressions()[name]


FUNCTIONS = (
    (0, "psi", "heart", "psi", "scalar"),
    (1, "gradPsi", "heart", "gradPsi", "vector"),
    (2, "source", "heart", "source", "scalar"),
    (3, "sigma", "heart", "G", "symmTensor"),
)
SIGNATURES = {"heart": "const point& pt, const Params& p"}
RETURN_TYPE = {"scalar": "scalar", "vector": "vector", "symmTensor": "symmTensor"}
ZONES_OF_KIND = {"heart": (1,)}
N_COMPONENTS = {"scalar": 1, "vector": 3, "symmTensor": 6}
PREAMBLE = (
    (x, (), "const scalar x = pt.x();"),
    (y, (), "const scalar y = pt.y();"),
    (z, (), "const scalar z = pt.z();"),
)
PREAMBLE_NEEDS = {}


@functools.lru_cache(maxsize=None)
def _compiled(name, lib):
    return sp.lambdify(ARGS, expressions()[name], modules=lib, cse=True)


def evaluate(name, pts, params=None, lib="numpy"):
    pv = param_values(params or default_params())
    if lib == "numpy":
        pts = np.asarray(pts, dtype=float)
        out = _compiled(name, "numpy")(pts[..., 0], pts[..., 1], pts[..., 2],
                                       *[pv[n] for n in PARAM_ORDER])
        if isinstance(out, (list, tuple)):
            return np.stack([np.broadcast_to(np.asarray(c, dtype=float), pts.shape[:-1])
                             for c in out], -1)
        return np.broadcast_to(np.asarray(out, dtype=float), pts.shape[:-1])
    pvm = [mp.mpf(pv[n]) for n in PARAM_ORDER]
    out = _compiled(name, "mpmath")(*[mp.mpf(c) for c in pts], *pvm)
    return [mp.mpf(c) for c in out] if isinstance(out, (list, tuple)) else mp.mpf(out)


def _full(T6):
    """symmTensor components (xx xy xz yy yz zz) as 3x3 matrices, stacked last."""
    return np.stack([np.stack([T6[..., 0], T6[..., 1], T6[..., 2]], -1),
                     np.stack([T6[..., 1], T6[..., 3], T6[..., 4]], -1),
                     np.stack([T6[..., 2], T6[..., 4], T6[..., 5]], -1)], -2)


def _M(pts, params=None):
    p = params or default_params()
    return _full(evaluate("G", pts, p))/(p["chi"]*p["cm"])


def _fibonacci_sphere(n, radius):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2*i/n)
    th = math.pi*(1 + 5**0.5)*i
    return radius*np.stack([np.cos(th)*np.sin(phi), np.sin(th)*np.sin(phi), np.cos(phi)], -1)


def _shell_points(n, band, seed=1, params=None):
    """Random points with r1 - band <= |x| <= r2 + band."""
    p = params or default_params()
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    r = rng.uniform(p["r1"] - band, p["r2"] + band, size=n)
    return u*r[:, None]


def max_sealing_residual(n=200, params=None):
    """max |n.G grad(psi)| on both spheres, mpmath at SAMPLE_DPS digits."""
    mp.mp.dps = SAMPLE_DPS
    p = params or default_params()
    worst = mp.mpf(0)
    for radius, sign in ((p["r1"], -1), (p["r2"], 1)):
        for pt in _fibonacci_sphere(n, radius):
            # project onto the sphere at working precision: the sealing identity
            # holds on |x| = radius exactly, not on its double-precision image
            ptm = [mp.mpf(float(c)) for c in pt]
            r = mp.sqrt(sum(c*c for c in ptm))
            ptm = [mp.mpf(radius)*c/r for c in ptm]
            nrm = [sign*c/mp.mpf(radius) for c in ptm]
            gP = evaluate("gradPsi", ptm, p, "mpmath")
            Gc = evaluate("G", ptm, p, "mpmath")
            Gm = [[Gc[0], Gc[1], Gc[2]], [Gc[1], Gc[3], Gc[4]], [Gc[2], Gc[4], Gc[5]]]
            res = sum(nrm[i]*Gm[i][j]*gP[j] for i in range(3) for j in range(3))
            worst = max(worst, abs(res))
    return float(worst)


def max_gradient_fd_error(h=1e-3, n=40, params=None):
    """gradPsi against a fourth-order central difference of psi."""
    p = params or default_params()

    def psi(q):
        return float(evaluate("psi", q, p))
    err = 0.0
    for pt in _shell_points(n, 0.0, seed=7, params=p):
        fd = np.zeros(3)
        for a in range(3):
            e = np.zeros(3)
            e[a] = h
            fd[a] = (psi(pt - 2*e) - 8*psi(pt - e) + 8*psi(pt + e) - psi(pt + 2*e))/(12*h)
        err = max(err, np.abs(fd - evaluate("gradPsi", pt, p)).max())
    return err


def _flux(q, params=None):
    return np.einsum("...ij,...j->...i", _M(q, params), evaluate("gradPsi", q, params))


def max_source_fd_error(h=1e-5, n=40, params=None):
    """source against central differences of the flux M grad(psi).

    The metric term needs no derivative, so this holds at the critical points
    of psi as well as away from them; the error is the O(h^2) truncation of
    the difference.
    """
    p = params or default_params()
    err = 0.0
    for pt in _shell_points(n, 0.0, seed=2, params=p):
        fd = 0.0
        for a in range(3):
            e = np.zeros(3)
            e[a] = h
            fd += (_flux(pt + e, p)[a] - _flux(pt - e, p)[a])/(2*h)
        gP = evaluate("gradPsi", pt, p)
        metric = math.sqrt(float(gP @ _M(pt, p) @ gP))
        exact = -fd + p["c0"]*metric - 1
        err = max(err, abs(exact - float(evaluate("source", pt, p))))
    return err


def min_denominator(band, params=None):
    """Smallest n.G n of both spheres over a band around the shell."""
    p = params or default_params()
    return float(evaluate("denominators", _shell_points(4000, band, seed=3, params=p), p).min())


def _volume_points(n, seed, params=None):
    """Points uniform in the volume of the shell."""
    p = params or default_params()
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    r = (p["r1"]**3 + rng.random(n)*(p["r2"]**3 - p["r1"]**3))**(1/3)
    return u*r[:, None]


def source_statistics(n=200000, params=None):
    """Finiteness and size of the source over the shell and a band around it, and the
    range of the metric term 1 + S outside the pinned ball."""
    p = params or default_params()
    pts = _shell_points(40000, 0.05, seed=8, params=p)
    s = evaluate("source", pts, p)
    inner = _volume_points(n, seed=11, params=p)
    t = evaluate("source", inner, p)
    outside = np.linalg.norm(inner - np.array(p["xs"]), axis=1) >= p["rPin"]
    t = 1 + t[outside]
    return {"finite": bool(np.isfinite(s).all() and np.isfinite(t).all()),
            "maxAbs": float(np.abs(s).max()), "median1pS": float(np.median(t)),
            "min1pS": float(t.min()), "fractionNegative": float(np.mean(t < 0))}


def _gradient(q, params=None):
    return np.asarray(evaluate("gradPsi", q, params), dtype=float)


def _hessian(q, h=1e-5, params=None):
    H = np.zeros((3, 3))
    for a in range(3):
        e = np.zeros(3)
        e[a] = h
        H[:, a] = (_gradient(q + e, params) - _gradient(q - e, params))/(2*h)
    return 0.5*(H + H.T)


def _tangent_basis(n):
    a = np.cross(n, [1.0, 0.0, 0.3])
    a /= np.linalg.norm(a)
    return np.array([a, np.cross(n, a)])


def critical_points(params=None):
    """The critical points of psi in the closed shell, grad(psi) = 0.

    On a wall the sealing n.M grad(psi) = 0 makes grad(psi) = 0 equivalent to a
    vanishing tangential gradient, so the wall points are the critical points of psi
    on the sphere; they are found from the local minima of the tangential gradient on
    a Fibonacci grid, refined by Levenberg-Marquardt, and are complete when the
    Morse count (minima - saddles + maxima) is 2 on each sphere.  The interior
    points are refined from xs and from the smallest-gradient random samples.
    Each point carries its position, wall (1, 2, or 0 in the interior), psi, type
    (from the Hessian tangential to the wall, or the full Hessian inside), the
    Hessian eigenvalues of that classification, the normal second derivative on a
    wall and whether it lies in the pinned ball.
    """
    p = params or default_params()
    gradient = functools.partial(_gradient, params=p)
    xs = np.array(p["xs"])
    kw = dict(method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15)
    found = []

    def add(q, wall):
        if np.linalg.norm(gradient(q)) > 1e-9 or any(np.linalg.norm(q - c["x"]) < 1e-6 for c in found):
            return
        H = _hessian(q, params=p)
        if wall:
            n = q/np.linalg.norm(q)
            T = _tangent_basis(n)
            eig = np.linalg.eigvalsh(T @ H @ T.T)
            normal = float(n @ H @ n)
        else:
            eig = np.linalg.eigvalsh(H)
            normal = float("nan")
        found.append({"x": q, "r": float(np.linalg.norm(q)), "wall": wall,
                      "psi": float(evaluate("psi", q, p)), "eig": eig, "normal": normal,
                      "type": {0: "maximum", 1: "saddle", 2: "minimum"}[int((eig > 0).sum())] if wall
                      else {0: "maximum", 1: "saddle", 2: "saddle", 3: "minimum"}[int((eig > 0).sum())],
                      "inPin": bool(np.linalg.norm(q - xs) < p["rPin"])})

    from scipy.spatial import cKDTree
    for wall, radius in ((1, p["r1"]), (2, p["r2"])):
        pts = _fibonacci_sphere(40000, radius)
        nrm = pts/radius
        g = gradient(pts)
        gt = np.linalg.norm(g - np.sum(g*nrm, -1)[:, None]*nrm, axis=1)
        nbr = cKDTree(pts).query(pts, k=11)[1][:, 1:]
        for i in np.where((gt[:, None] <= gt[nbr]).all(1))[0]:
            def residual(u):
                q = radius*u/np.linalg.norm(u)
                gq = gradient(q)
                nq = q/radius
                return gq - (gq @ nq)*nq
            sol = least_squares(residual, pts[i]/radius, **kw)
            add(radius*sol.x/np.linalg.norm(sol.x), wall)

    pts = _volume_points(100000, seed=10, params=p)
    gm = np.linalg.norm(gradient(pts), axis=1)
    seeds = [xs] + [pts[i] for i in np.argsort(gm)[:2000:50]]
    for s in seeds:
        sol = least_squares(gradient, s, **kw)
        r = np.linalg.norm(sol.x)
        if p["r1"] + 1e-7 < r < p["r2"] - 1e-7:
            add(sol.x, 0)

    return sorted(found, key=lambda c: (c["wall"], c["inPin"], -c["psi"]))


def morse_count(cps, wall):
    """minima - saddles + maxima of psi on the given wall (2 on a sphere)."""
    sign = {"minimum": 1, "saddle": -1, "maximum": 1}
    return sum(sign[c["type"]] for c in cps if c["wall"] == wall)


def format_critical_points(cps):
    lines = ["critical points of psi in the closed shell (grad psi = 0); on a wall the type and the "
             "two eigenvalues are those of the Hessian tangential to the sphere, then the normal "
             "second derivative, inside the three Hessian eigenvalues:",
             "  wall  inPin  type     x                          r        psi       eigenvalues"]
    for c in cps:
        if c["wall"]:
            tail = "(%9.3f %9.3f)  n.H.n %9.3f" % (*c["eig"], c["normal"])
        else:
            tail = "(%9.3f %9.3f %9.3f)" % tuple(c["eig"])
        lines.append("  %-4d  %-5s  %-8s (%8.5f %8.5f %8.5f)  %.5f  %8.4f  %s" % (
            c["wall"], "yes" if c["inPin"] else "no", c["type"], *c["x"], c["r"], c["psi"], tail))
    lines.append("Morse count (minima - saddles + maxima): inner wall %d, outer wall %d"
                 % (morse_count(cps, 1), morse_count(cps, 2)))
    return "\n".join(lines)


def min_nondegeneracy(cps):
    """Smallest ratio of the smallest to the largest absolute Hessian eigenvalue of a critical
    point outside the pinned ball (zero would be a degenerate point)."""
    return min(np.abs(c["eig"]).min()/np.abs(c["eig"]).max() for c in cps if not c["inPin"])


def smallest_value_outside_pin(cps):
    """Smallest psi among the local minima outside the pinned ball."""
    return min(c["psi"] for c in cps if c["type"] == "minimum" and not c["inPin"])


def pinned_minimum(cps):
    """The interior minimum of psi, which lies in the pinned ball and is the only point
    where psi can be negative."""
    return min(c["psi"] for c in cps if c["inPin"])


def report(params=None):
    p = params or default_params()
    cps = critical_points(p)
    stats = source_statistics(params=p)
    print("c0 = %g; stimulus xs = (%.6f, %.6f, %.6f), pinned ball radius %.3f"
          % (p["c0"], *p["xs"], p["rPin"]))
    print("psi: minimum in the pinned ball %.5f, smallest local minimum outside it %.5f, "
          "largest critical value %.5f"
          % (pinned_minimum(cps), smallest_value_outside_pin(cps), max(c["psi"] for c in cps)))
    print("sealing residual (30 digits): %.3e" % max_sealing_residual(params=p))
    print("gradient against finite differences of psi: %.3e" % max_gradient_fd_error(params=p))
    print("source against finite differences of the flux: %.3e" % max_source_fd_error(params=p))
    print("smallest n.G n over the band: %.5f" % min_denominator(0.05, p))
    print("source: finite %s, max |S| %.3f; outside the pinned ball median 1+S %.3f, min %.3f, "
          "fraction 1+S<0 %.3f" % (stats["finite"], stats["maxAbs"], stats["median1pS"],
                                    stats["min1pS"], stats["fractionNegative"]))
    print(format_critical_points(cps))
    return cps


def run_self_checks(params=None):
    """Assert every self-check at params (default_params(), whose c0 may be changed);
    return the critical points."""
    p = params or default_params()
    scale = p["chi"]*p["cm"]
    G6 = evaluate("G", p["xs"], p, "mpmath")
    Ms = np.array([[G6[0], G6[1], G6[2]], [G6[1], G6[3], G6[4]], [G6[2], G6[4], G6[5]]],
                  dtype=float)/scale
    Msi = np.array([[p["Msi"][0], p["Msi"][1], p["Msi"][2]], [p["Msi"][1], p["Msi"][3], p["Msi"][4]],
                    [p["Msi"][2], p["Msi"][4], p["Msi"][5]]])
    assert np.abs(Msi @ Ms - np.eye(3)).max() < 1e-12, "stimulus metric"
    assert abs(math.hypot(*p["xs"]) - STIMULUS["default"][1]) < 1e-15, "stimulus radius"
    assert max_sealing_residual(params=p) < 1e-25, "sealing"
    assert max_gradient_fd_error(params=p) < 1e-5, "gradient"
    assert max_source_fd_error(params=p) < 1e-6, "source"
    assert min_denominator(0.05, p) > 0, "denominator"
    assert source_statistics(params=p)["finite"], "source finite"
    cps = critical_points(p)
    assert morse_count(cps, 1) == 2 and morse_count(cps, 2) == 2, "critical points of the walls"
    assert all(c["wall"] or c["inPin"] for c in cps), "interior critical point outside the pinned ball"
    assert min_nondegeneracy(cps) > 1e-3, "nondegenerate critical points"
    assert smallest_value_outside_pin(cps) > 0 and pinned_minimum(cps) > -1e-2, "psi positive"
    for c in cps:
        near = c["x"] + 1e-6*np.array([0.6, -0.5, 0.4])
        assert np.isfinite(float(evaluate("source", near, p))), "source at critical point"
    return cps


def printed_check_context():
    pv = param_values(alternate_params())
    ctx = {"x": 0.61, "y": -0.29, "z": 0.43}
    ctx.update(pv)
    return ctx


SAMPLE_POINTS = [(0.36, 0.27, -0.23), (0.55, -0.47, 0.36), (-0.12, 0.64, 0.71)]


def self_test_samples():
    mp.mp.dps = SAMPLE_DPS
    rows = []
    for idx, cname, kind, ename, vkind in FUNCTIONS:
        for param, params in ((0, default_params()), (1, alternate_params())):
            for pt in SAMPLE_POINTS:
                val = evaluate(ename, [mp.mpf(float(c)) for c in pt], params, "mpmath")
                vals = val if isinstance(val, list) else [val]
                for comp in range(N_COMPONENTS[vkind]):
                    rows.append((idx, param, 0, 1, comp, pt[0], pt[1], pt[2], 0.0, float(vals[comp])))
    return rows


def latex_fragments():
    q, div_flux, quad, Gs = sp.symbols("q divflux quad Gs", real=True)
    names = smm.shell_symbol_names([
        (q, r"(\mathbf x-\mathbf x_s)\cdot\mathbf M_s^{-1}(\mathbf x-\mathbf x_s)"),
        (div_flux, r"\nabla\!\cdot\!\left(\mathbf M\nabla\psi\right)"),
        (quad, r"\nabla\psi\cdot\mathbf M\nabla\psi"), (Gs, r"\boldsymbol\sigma")])
    direction, radius = STIMULUS["default"]
    d = sp.Matrix([sp.nsimplify(c) for c in direction])
    out = {
        "shell_eikonal_hhat.tex": emit_latex_fragment(r"\hat h", _hhat(q), names),
        "shell_eikonal_metric.tex": emit_latex_fragment(r"\mathbf M", _metric_tensor(Gs), names),
        "shell_eikonal_source.tex": emit_latex_fragment(r"S", _operator(div_flux, quad), names),
        "shell_eikonal_xs.tex": emit_latex_fragment(
            r"\mathbf x_s", sp.nsimplify(radius)*d/sp.sqrt((d.T*d)[0]), names),
        "shell_eikonal_Ms.tex": emit_latex_fragment(
            r"\mathbf M_s", sp.Symbol(r"M(x_s)", real=True),
            {sp.Symbol(r"M(x_s)", real=True): r"\mathbf M(\mathbf x_s)"}),
    }
    out.update(smm.blend_fragments("shell_eikonal", r"\mathbf M", r"\psi"))
    return {"latex/" + k: v for k, v in out.items()}


def numpy_module():
    e = expressions()
    functions = {"psi": e["psi"], "gradPsi": sp.Matrix(e["gradPsi"]), "source": e["source"],
                 "sigma": e["Gmat"]}
    doc = ("Heart-only eikonal manufactured solution on the spherical shell r1 < |x| < r2: "
           "psi, gradPsi (shape (3, ...)), the source S of -div(M grad psi) + c0 sqrt(grad psi."
           "M grad psi) = 1 + S and sigma (shape (3, 3, ...)), each a function of (x, y, z, "
           "*PARAM_ORDER) with angles in radians and M = sigma/(chi cm).  PARAMS holds the "
           "default and alternate parameters; field(name, x, y, z, **overrides) evaluates with "
           "them.  c0 may be overridden freely; xs and Msi = M(xs)^-1 belong together.  "
           "Generated by shell_eikonal_mms.py.")
    return smm.reference_module(_THIS, functions, doc)


def render_outputs():
    """{path relative to this folder: text} of every emitted file."""
    out = {"reference/%s.py" % REFERENCE_NAME: numpy_module()}
    out.update(latex_fragments())
    return out


def write_outputs(root=HERE):
    return smm.write_outputs(render_outputs, root)


def outputs_match(root=HERE):
    return smm.outputs_match(render_outputs, root)


def _positive_number(text):
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError("not a number: %r" % text)
    if not (math.isfinite(value) and value > 0):
        raise argparse.ArgumentTypeError("must be a positive finite number: %r" % text)
    return value


LICENSE = """\
/*---------------------------------------------------------------------------*\\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Exact activation time psi, its gradient, the manufactured eikonal source
    and the conductivity sigma of the heart-only eikonal manufactured
    solution on the spherical shell, sealed on both spheres, with the
    stimulus point xs, the pinned radius rPin and self-test samples.

    Regenerate from generator/:
        python shell_eikonal_mms.py --write %s

    sympy %s
\\*---------------------------------------------------------------------------*/
"""


def generate_header():
    guard = "manufacturedShellEikonalReference_H"
    L = [LICENSE % (HEADER_PATH, sp.__version__)]
    L += ["#ifndef " + guard, "#define " + guard, ""]
    L += codegen.INCLUDES
    L += ["namespace Foam", "{", "namespace shellEikonalMMS", "{", ""]
    L += ["struct Params", "{", "    scalar r1, r2, gamma1, gamma2;",
          "    scalar gil, git, gin, theta0, theta1;",
          "    scalar chi, cm, c0, delta, rPin;",
          "    scalar xs[3];", "    scalar Msi[6];",
          "    label solution;", "};", ""]
    L += codegen.params_function(_THIS, "defaultParams", default_params(), 0) + [""]
    L += codegen.params_function(_THIS, "alternateParams", alternate_params(), 0) + [""]
    for spec in FUNCTIONS:
        L += codegen.emit_function(_THIS, *spec) + [""]
    L += codegen.self_test_table(_rows_cache())
    L += ["} // End namespace shellEikonalMMS", "} // End namespace Foam", "",
          "#endif", ""]
    return "\n".join(L)


def header_matches(path, text=None):
    return codegen.header_matches(__file__, generate_header, path, text)


@functools.lru_cache(maxsize=None)
def _rows_cache():
    return tuple(self_test_samples())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--write", metavar="PATH")
    grp.add_argument("--check", metavar="PATH")
    grp.add_argument("--report", action="store_true")
    grp.add_argument("--self-check", action="store_true")
    grp.add_argument("--emit", action="store_true")
    ap.add_argument("--c0", type=_positive_number, metavar="VALUE",
                    help="c0 of --report and --self-check (default 1)")
    args = ap.parse_args(argv)
    if args.c0 is not None and not (args.report or args.self_check):
        ap.error("--c0 applies to --report and --self-check; the header and the emitted "
                 "reference carry c0 as a runtime parameter")
    params = default_params()
    if args.c0 is not None:
        params["c0"] = args.c0
    if args.report:
        report(params)
        return 0
    if args.self_check or args.emit:
        try:
            run_self_checks(params)
        except AssertionError as err:
            print("self-check failed at c0 = %g: %s" % (params["c0"], err))
            return 1
        if args.self_check:
            print("self-checks passed at c0 = %g" % params["c0"])
            return 0
        for rel in write_outputs():
            print("wrote %s" % rel)
        return 0
    if args.write:
        try:
            run_self_checks()
        except AssertionError as err:
            print("self-check failed, header not written: %s" % err)
            return 1
        report()
        text = generate_header()
        path = codegen.resolve(__file__, args.write)
        with open(path, "w") as fh:
            fh.write(text)
        print("wrote %s (%d bytes, %d self-test samples)" % (path, len(text.encode()), len(_rows_cache())))
        return 0
    if header_matches(args.check):
        print("%s is up to date" % codegen.resolve(__file__, args.check))
        return 0
    print("%s differs from the generator output" % codegen.resolve(__file__, args.check))
    return 1


_THIS = sys.modules[__name__]

if __name__ == "__main__":
    sys.exit(main())
