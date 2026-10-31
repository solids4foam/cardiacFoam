#!/usr/bin/env python3
"""Heart-only bidomain manufactured solution on the spherical shell.

The heart is the shell r1 < |x| < r2; both spheres are sealed walls.  The
conductivities are the intracellular Gi and extracellular Ge tensors
(unequal eigenvalue ratios, one fibre frame, fibre and sheet rotating
through the wall), so Gi + Ge is not proportional to Gi and both are full and
oblique to both walls.  Fi = hhat_i + omega*m_i satisfies n.Gi grad(Fi) = 0 and
Fe = hhat_e + omega*m_e satisfies n.Ge grad(Fe) = 0 exactly on each sphere,
with omega zero on both spheres and m the blend of the per-sphere corrections.
Fi is the field of shell_monodomain_mms; hhat_e is a smooth field built from
polynomial, trigonometric and exponential terms.  The verifier uses
phiI = sqrt(1 + t) Fi, phiE = sqrt(1 + t) Fe, Vm = phiI - phiE and the FDA
ionic closure with F = Fi - Fe and the positive function g.  The phiE equation
carries the source sqrt(1 + t) (div(Gi grad Fi) + div(Ge grad Fe)), whose
integral over the shell is zero.

--emit writes reference/shell_bidomain.py (Fi, gradFi, Fe, gradFe,
divGiGradFi, divGeGradFe, Gi, Ge and g as NumPy functions of x, y, z and the
parameters, with PARAMS and field()) and the LaTeX fragments
latex/shell_bidomain_*.tex: the tensors Gi and Ge, hhat_i, hhat_e, the wall
corrections and blends of both potentials, Fi, Fe, phi_i, phi_e, V_m, the
gauge constant c, u_1, u_2, the applied current of the V_m equation, the
phi_e source and the parameters.  The shared geometry fragments
latex/shell_*.tex are written by shell_monodomain_mms.py.

    python shell_bidomain_mms.py --emit
    python shell_bidomain_mms.py --self-check
    python shell_bidomain_mms.py --write ../../../src/verificationModels/bidomainVerification/manufacturedShellBidomainReference.H
    python shell_bidomain_mms.py --check ../../../src/verificationModels/bidomainVerification/manufacturedShellBidomainReference.H
"""

import argparse
import functools
import math
import os
import sys

import mpmath as mp
import numpy as np
import sympy as sp

import fda_monodomain_mms as fmm
import mms_codegen as codegen
import shell_monodomain_mms as smm
from mms_common.latex_emit import emit_latex_fragment, emit_parameter_fragment

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "shell_bidomain"
# phiERefPoint of the tutorial manufacturedSolutions/bidomain/sphericalShell
PHIE_REF_POINT = (0.6, 0.3, 0.2)
PHIE_REF_VALUE = 0

HEADER_NAME = "manufacturedShellBidomainReference.H"
HEADER_PATH = "../../../src/verificationModels/bidomainVerification/" + HEADER_NAME
SAMPLE_DPS = 30

SCALAR_PARAMS = ("r1", "r2", "gamma1", "gamma2", "k",
                 "gil", "git", "gin", "gel", "get", "gen", "theta0", "theta1")
TAIL_PARAMS = ()
ARRAY_PARAMS = (("fA", 4), ("hA", 3))
ANGLE_PARAMS = ("gamma1", "gamma2", "theta0", "theta1")
PARAM_ORDER = list(SCALAR_PARAMS) + ["fA%d" % i for i in range(4)] + ["hA%d" % i for i in range(3)]
P = {name: sp.Symbol(name, real=True) for name in PARAM_ORDER}

x, y, z = sp.symbols("x y z", real=True)
ARGS = [x, y, z] + [P[n] for n in PARAM_ORDER]
XYZ = sp.Matrix([x, y, z])

SOLUTIONS = ("shell",)
HEART_ZONE = 1
ZONE_NAMES = {1: "myocardium"}


def default_params():
    return {"r1": 0.5, "r2": 1.0, "gamma1": 20.0, "gamma2": 50.0, "k": math.pi,
            "gil": 0.17, "git": 0.019, "gin": 0.01,
            "gel": 0.62, "get": 0.12, "gen": 0.08,
            "theta0": -20.0, "theta1": 90.0,
            "fA": [0.3, 0.2, 0.4, 0.15], "hA": [-0.1, 0.3, 0.25]}


def alternate_params():
    p = default_params()
    p["fA"] = [0.25, 0.25, 0.35, 0.2]
    p["hA"] = [-0.15, 0.35, 0.15]
    p["theta1"] = 70.0
    return p


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


def _frame():
    """Fixed orthonormal frame (e_a, e_b, n_s): columns of Rz(gamma1) Rx(gamma2)."""
    c1, s1 = sp.cos(P["gamma1"]), sp.sin(P["gamma1"])
    c2, s2 = sp.cos(P["gamma2"]), sp.sin(P["gamma2"])
    Rz = sp.Matrix([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]])
    Rx = sp.Matrix([[1, 0, 0], [0, c2, -s2], [0, s2, c2]])
    Q = Rz*Rx
    return Q[:, 0], Q[:, 1], Q[:, 2]


def _tensor(l, tr, nn, theta):
    """l f f + tr s s + nn n_s n_s with f, s rotated by theta in the (e_a, e_b) plane."""
    ea, eb, ec = _frame()
    f = sp.cos(theta)*ea + sp.sin(theta)*eb
    s = -sp.sin(theta)*ea + sp.cos(theta)*eb
    return l*f*f.T + tr*s*s.T + nn*ec*ec.T


def _ionic_G():
    k = P["k"]
    return 1 + sp.Rational(1, 4)*sp.sin(k*x + sp.Rational(1, 5))*sp.cos(k*y - sp.Rational(1, 3))*sp.cos(k*z)


def _components(T):
    return [T[0, 0], T[0, 1], T[0, 2], T[1, 1], T[1, 2], T[2, 2]]


@functools.lru_cache(maxsize=None)
def expressions():
    r1, r2, k = P["r1"], P["r2"], P["k"]
    rho = x**2 + y**2 + z**2
    lam = (rho - r1**2)/(r2**2 - r1**2)
    om = (rho - r1**2)*(rho - r2**2)/((r2**2 - r1**2)*(r1 + r2))
    theta = P["theta0"] + P["theta1"]*lam
    Gi = smm.expressions()["Gmat"]
    Ge = _tensor(P["gel"], P["get"], P["gen"], theta)
    a = [P["fA%d" % i] for i in range(4)]
    ehat = (a[0]*x + a[1]*y*z
            + a[2]*sp.sin(k*x)*(1 + sp.Rational(1, 2)*y)
            + a[3]*sp.cos(k*y + sp.Rational(3, 10))*sp.exp(sp.Rational(2, 5)*z))
    # unit normals pointing out of the heart; d_j = grad(omega).n_j on sphere j
    n = {1: -XYZ/r1, 2: XYZ/r2}
    d = {1: 2*r1/(r1 + r2), 2: 2*r2/(r1 + r2)}
    mj = {j: -_dot(n[j], Ge*_grad(ehat))/(_dot(n[j], Ge*n[j])*d[j]) for j in (1, 2)}
    m = (1 - lam)*mj[1] + lam*mj[2]
    Fi = smm.expressions()["F"]
    Fe = ehat + om*m
    gradFi, gradFe = _grad(Fi), _grad(Fe)
    g = _ionic_G()
    return {"Fi": Fi, "gradFi": list(gradFi), "Fe": Fe, "gradFe": list(gradFe),
            "divGiGradFi": _div(Gi*gradFi), "divGeGradFe": _div(Ge*gradFe),
            "Gi": _components(Gi), "Ge": _components(Ge), "g": g,
            "ehat": ehat, "n1": n[1], "n2": n[2], "Gimat": Gi, "Gemat": Ge,
            "denominators": [_dot(n[j], T*n[j]) for T in (Gi, Ge) for j in (1, 2)]}


def _expr(solution, zone, name):
    return expressions()[name]


FUNCTIONS = (
    (0, "Fi", "heart", "Fi", "scalar"),
    (1, "gradFi", "heart", "gradFi", "vector"),
    (2, "Fe", "heart", "Fe", "scalar"),
    (3, "gradFe", "heart", "gradFe", "vector"),
    (4, "divGiGradFi", "heart", "divGiGradFi", "scalar"),
    (5, "divGeGradFe", "heart", "divGeGradFe", "scalar"),
    (6, "Gi", "heart", "Gi", "symmTensor"),
    (7, "Ge", "heart", "Ge", "symmTensor"),
    (8, "g", "heart", "g", "scalar"),
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


def _fibonacci_sphere(n, radius):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2*i/n)
    th = math.pi*(1 + 5**0.5)*i
    return radius*np.stack([np.cos(th)*np.sin(phi), np.sin(th)*np.sin(phi), np.cos(phi)], -1)


def _shell_points(n, band, seed=1):
    """Random points with r1 - band <= |x| <= r2 + band."""
    p = default_params()
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    r = rng.uniform(p["r1"] - band, p["r2"] + band, size=n)
    return u*r[:, None]


# (potential, its gradient, its conductivity, its divergence of flux)
PAIRS = (("Fi", "gradFi", "Gi", "divGiGradFi"), ("Fe", "gradFe", "Ge", "divGeGradFe"))


def max_sealing_residual(n=200, params=None):
    """max |n.G grad(F)| on both spheres for (Fi, Gi) and (Fe, Ge), mpmath at SAMPLE_DPS digits."""
    mp.mp.dps = SAMPLE_DPS
    p = params or default_params()
    worst = mp.mpf(0)
    for _, gname, tname, _ in PAIRS:
        for radius, sign in ((p["r1"], -1), (p["r2"], 1)):
            for pt in _fibonacci_sphere(n, radius):
                # project onto the sphere at working precision: the sealing identity
                # holds on |x| = radius exactly, not on its double-precision image
                ptm = [mp.mpf(float(c)) for c in pt]
                r = mp.sqrt(sum(c*c for c in ptm))
                ptm = [mp.mpf(radius)*c/r for c in ptm]
                nrm = [sign*c/mp.mpf(radius) for c in ptm]
                gF = evaluate(gname, ptm, p, "mpmath")
                Gc = evaluate(tname, ptm, p, "mpmath")
                Gm = [[Gc[0], Gc[1], Gc[2]], [Gc[1], Gc[3], Gc[4]], [Gc[2], Gc[4], Gc[5]]]
                res = sum(nrm[i]*Gm[i][j]*gF[j] for i in range(3) for j in range(3))
                worst = max(worst, abs(res))
    return float(worst)


def max_source_fd_error(h=1e-5, n=40):
    """divGiGradFi and divGeGradFe against central differences of the fluxes.

    The error is the O(h^2) truncation of the difference (it falls by 100 per
    decade of h), set by the large third derivatives of the flux (m carries
    1/(n.G n) with gin small), so h is chosen small enough for 1e-6.
    """
    pts = _shell_points(n, 0.0, seed=2)
    err = 0.0
    for _, gname, tname, dname in PAIRS:
        for pt in pts:
            fd = 0.0
            for a in range(3):
                e = np.zeros(3)
                e[a] = h
                fd += (_flux(gname, tname, pt + e)[a] - _flux(gname, tname, pt - e)[a])/(2*h)
            err = max(err, abs(fd - float(evaluate(dname, pt))))
    return err


def _flux(gname, tname, q):
    return _full(evaluate(tname, q)) @ evaluate(gname, q)


def min_denominator(band):
    """Smallest n.Gi n and n.Ge n of both spheres over a band around the shell."""
    pts = _shell_points(4000, band, seed=3)
    return float(evaluate("denominators", pts).min())


def min_g(band):
    return float(evaluate("g", _shell_points(4000, band, seed=4)).min())


def wall_richness(n=400):
    """Smallest over both walls of the wall mean of the tangential parts of Gi n/gil, Ge n/gel,
    grad Fi, grad Fe and grad F = grad Fi - grad Fe.

    The tangential part of G n is a tangent vector field on a sphere, so it
    vanishes at isolated points (hairy-ball theorem) and its pointwise minimum
    is zero in the continuum; the wall mean measures how oblique G is.
    """
    p = default_params()
    out = {"GiN": np.inf, "GeN": np.inf, "gradFi": np.inf, "gradFe": np.inf, "gradF": np.inf}
    for radius, sign in ((p["r1"], -1), (p["r2"], 1)):
        pts = _fibonacci_sphere(n, radius)
        nrm = sign*pts/radius
        tang = lambda v: v - np.sum(v*nrm, -1)[:, None]*nrm
        mean = lambda v: float(np.linalg.norm(tang(v), axis=1).mean())
        giN = np.einsum("pij,pj->pi", _full(evaluate("Gi", pts)), nrm)
        geN = np.einsum("pij,pj->pi", _full(evaluate("Ge", pts)), nrm)
        gFi, gFe = evaluate("gradFi", pts), evaluate("gradFe", pts)
        for key, val in (("GiN", mean(giN)/p["gil"]), ("GeN", mean(geN)/p["gel"]),
                         ("gradFi", mean(gFi)), ("gradFe", mean(gFe)),
                         ("gradF", mean(gFi - gFe))):
            out[key] = min(out[key], val)
    return out


def phiE_source_mean(nr=60, nc=60, nphi=120):
    """Integral of divGiGradFi + divGeGradFe over the shell and of its absolute value.

    Gauss-Legendre in r and cos(theta), uniform in phi; the integral is zero
    in the continuum because both fluxes vanish on both spheres.
    """
    p = default_params()
    xr, wr = np.polynomial.legendre.leggauss(nr)
    xc, wc = np.polynomial.legendre.leggauss(nc)
    r = 0.5*(p["r2"] - p["r1"])*xr + 0.5*(p["r2"] + p["r1"])
    wr = 0.5*(p["r2"] - p["r1"])*wr*r**2
    ph = 2*math.pi*(np.arange(nphi) + 0.5)/nphi
    R, C, PH = np.meshgrid(r, xc, ph, indexing="ij")
    S = np.sqrt(1 - C**2)
    pts = np.stack([R*S*np.cos(PH), R*S*np.sin(PH), R*C], -1)
    W = wr[:, None, None]*wc[None, :, None]*(2*math.pi/nphi)
    f = evaluate("divGiGradFi", pts) + evaluate("divGeGradFe", pts)
    return float((W*f).sum()), float((W*np.abs(f)).sum())


def max_difference_to_monodomain_F(n=200):
    """Fi against the field of shell_monodomain_mms."""
    pts = _shell_points(n, 0.0, seed=5)
    return float(np.abs(smm.evaluate("F", pts) - evaluate("Fi", pts)).max())


def frame_check(n=500):
    """Gi and Ge share the eigenvectors (f, s, n_s); their eigenvalue ratios
    differ, so Gi + Ge is not a multiple of Gi.  Returns the largest commutator
    [Gi, Ge], the smallest relative distance of Ge from the multiples of Gi, and
    the largest difference between Gi and sigma of shell_monodomain_mms."""
    pts = _shell_points(n, 0.0, seed=7)
    Gi = _full(evaluate("Gi", pts))
    Ge = _full(evaluate("Ge", pts))
    comm = np.abs(Gi @ Ge - Ge @ Gi).max()
    c = np.einsum("pij,pij->p", Ge, Gi)/np.einsum("pij,pij->p", Gi, Gi)
    off = np.linalg.norm(Ge - c[:, None, None]*Gi, axis=(1, 2))/np.linalg.norm(Ge, axis=(1, 2))
    mono = np.abs(Gi - _full(smm.evaluate("G", pts))).max()
    return float(comm), float(off.min()), float(mono)


def vm_equation_residual():
    """Residuals of the V_m equation and of the two ionic rates on the exact
    fields, with the ionic closure of monodomainFDAManufactured and the
    applied current of the verifier, the spatial factors written as symbols."""
    t = sp.Symbol("t", positive=True)
    Fi, Fe, g, b, dGi = sp.symbols("Fi Fe g beta dGi", real=True)
    chi, cm = sp.symbols("chi C_m", positive=True)
    S = sp.sqrt(1 + t)
    V = S*(Fi - Fe)
    a1 = (1 + t)*g + V
    a2 = 1/((1 + t)*sp.sqrt(g))
    J = -sp.Rational(1, 2)*cm*a2**2*V*(a1 - V) + b*V/chi
    Is = S*(b*(Fi - Fe) - dGi)
    res = chi*cm*sp.diff(V, t) - (S*dGi - chi*cm*J/cm + Is)
    rate1 = a2**2*(a1 - V)**2 + sp.Rational(1, 2)*(a1 - V)*a2**2*V
    rate2 = -(a1 - V)*a2**3
    return [sp.simplify(res), sp.simplify(sp.diff(a1, t) - rate1),
            sp.simplify(sp.diff(a2, t) - rate2)]


SOURCE_MEAN_TOLERANCE = 1e-12


def run_self_checks():
    assert max_sealing_residual() < 1e-25, "sealing"
    assert max_source_fd_error() < 1e-6, "source"
    assert min_denominator(0.05) > 0, "denominator"
    assert min_g(0.05) > 0.5, "g"
    assert all(v > 1e-2 for v in wall_richness().values()), "richness"
    mean, absolute = phiE_source_mean()
    assert abs(mean) < SOURCE_MEAN_TOLERANCE*absolute, "phiE source mean"
    assert max_difference_to_monodomain_F() < 1e-13, "monodomain F"
    comm, off, mono = frame_check()
    assert comm < 1e-14, "Gi and Ge share the frame"
    assert off > 0.1, "Ge is not proportional to Gi"
    assert mono < 1e-14, "Gi equals sigma of the monodomain shell"
    assert all(v == 0 for v in vm_equation_residual()), "V_m equation"
    r = math.sqrt(sum(c*c for c in PHIE_REF_POINT))
    p = default_params()
    assert p["r1"] < r < p["r2"], "reference point inside the shell"


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


# ------------------------------------------------------------------ LaTeX ----

# LaTeX names of the shell parameters: those of the monodomain shell, with the
# eigenvalues of Gi and Ge and the coefficients of hhat_e
SHELL_TEX = dict(smm.SHELL_SYMBOLS)
SHELL_TEX.update({"gil": r"\sigma_{i,f}", "git": r"\sigma_{i,s}", "gin": r"\sigma_{i,n}",
                  "gel": r"\sigma_{e,f}", "get": r"\sigma_{e,s}", "gen": r"\sigma_{e,n}",
                  "fA0": "b_1", "fA1": "b_2", "fA2": "b_3", "fA3": "b_4"})


def _shell_names(extra=()):
    names = {P[n]: tex for n, tex in SHELL_TEX.items() if n in P}
    names[sp.Symbol("k", real=True)] = "k"
    names.update(dict(extra))
    return names


def latex_fragments():
    out = {}
    ff, ss, nn = sp.symbols("ff ss nn", real=True)
    nm = _shell_names([(ff, r"\mathbf f\otimes\mathbf f"), (ss, r"\mathbf s\otimes\mathbf s"),
                       (nn, r"\mathbf n_s\otimes\mathbf n_s")])
    for tag, (l, tr, n) in (("Gi", ("gil", "git", "gin")), ("Ge", ("gel", "get", "gen"))):
        out["shell_bidomain_%s.tex" % tag] = emit_latex_fragment(
            r"\boldsymbol\sigma_%s" % tag[1], P[l]*ff + P[tr]*ss + P[n]*nn, nm)
    names = _shell_names()
    out["shell_bidomain_hhat_e.tex"] = emit_latex_fragment(
        r"\hat h_e", expressions()["ehat"], names)
    out["shell_bidomain_hhat_i.tex"] = emit_latex_fragment(r"\hat h_i", smm._hhat(), names)
    for tag, tens in (("i", r"\boldsymbol\sigma_i"), ("e", r"\boldsymbol\sigma_e")):
        a, b, dj, lam, m1, m2, hh, om, m = sp.symbols("a b dj lam m1 m2 hh om m", real=True)
        nmj = _shell_names([
            (a, r"\left(\mathbf n_j\cdot %s\nabla\hat h_%s\right)" % (tens, tag)),
            (b, r"\left(\mathbf n_j\cdot %s\,\mathbf n_j\right)" % tens),
            (dj, "d_j"), (lam, r"\lambda"), (m1, "m_{%s,1}" % tag), (m2, "m_{%s,2}" % tag),
            (hh, r"\hat h_%s" % tag), (om, r"\omega"), (m, "m_%s" % tag)])
        out["shell_bidomain_mj_%s.tex" % tag] = emit_latex_fragment(
            "m_{%s,j}" % tag, smm._wall_correction(a, b, dj), nmj)
        out["shell_bidomain_m_%s.tex" % tag] = emit_latex_fragment(
            "m_%s" % tag, smm._blend(lam, m1, m2), nmj)
        out["shell_bidomain_F%s.tex" % tag] = emit_latex_fragment(
            "F_%s" % tag, hh + om*m, nmj)
    t, Fi, Fe, g, b, Vm, dGi, dGe, shift = sp.symbols(
        "t Fi Fe g beta Vm dGi dGe shift", real=True)
    nf = _shell_names([(Fi, "F_i"), (Fe, "F_e"), (g, "g"), (b, r"\beta"), (Vm, "V_m"),
                       (dGi, r"\nabla\!\cdot\!\left(\boldsymbol\sigma_i\nabla F_i\right)"),
                       (dGe, r"\nabla\!\cdot\!\left(\boldsymbol\sigma_e\nabla F_e\right)"),
                       (shift, "c")])
    S = sp.sqrt(1 + t)
    out["shell_bidomain_phiI.tex"] = emit_latex_fragment(r"\phi_i", S*Fi + shift, nf)
    xr, ve = sp.symbols("xr ve", real=True)
    nf.update({xr: r"F_e(\mathbf x_{\mathrm{ref}})", ve: r"\phi_{e,\mathrm{ref}}"})
    out["shell_bidomain_c.tex"] = emit_latex_fragment("c", ve - S*xr, nf)
    out["shell_bidomain_phiE.tex"] = emit_latex_fragment(r"\phi_e", S*Fe + shift, nf)
    out["shell_bidomain_Vm.tex"] = emit_latex_fragment("V_m", S*(Fi - Fe), nf)
    out["shell_bidomain_u1.tex"] = emit_latex_fragment("u_1", (1 + t)*g + Vm, nf)
    out["shell_bidomain_u2.tex"] = emit_latex_fragment("u_2", g**sp.Rational(-1, 2)/(1 + t), nf)
    out["shell_bidomain_Is.tex"] = emit_latex_fragment(
        r"I_{\mathrm{app}}", S*(b*(Fi - Fe) - dGi), nf)
    out["shell_bidomain_se.tex"] = emit_latex_fragment(r"s_e", S*(dGi + dGe), nf)
    p = default_params()
    pairs = []
    for name in SCALAR_PARAMS:
        v = p[name]
        if name in ANGLE_PARAMS:
            v = r"%g^\circ" % v
        pairs.append((SHELL_TEX[name], v))
    pairs += [("a_%d" % (i + 1), v) for i, v in enumerate(p["hA"])]
    pairs += [("b_%d" % (i + 1), v) for i, v in enumerate(p["fA"])]
    out["shell_bidomain_params.tex"] = emit_parameter_fragment(pairs, per_row=5)
    out["shell_bidomain_ref.tex"] = emit_parameter_fragment(
        [(r"\mathbf x_{\mathrm{ref}}", sp.Matrix(PHIE_REF_POINT)),
         (r"\phi_{e,\mathrm{ref}}", PHIE_REF_VALUE), (r"\beta", fmm.BETA[3])])
    return {"latex/" + k: v for k, v in out.items()}


def numpy_module():
    e = expressions()
    functions = {"Fi": e["Fi"], "gradFi": sp.Matrix(e["gradFi"]), "Fe": e["Fe"],
                 "gradFe": sp.Matrix(e["gradFe"]), "divGiGradFi": e["divGiGradFi"],
                 "divGeGradFe": e["divGeGradFe"], "Gi": e["Gimat"], "Ge": e["Gemat"],
                 "g": e["g"]}
    doc = ("Heart-only bidomain manufactured solution on the spherical shell r1 < |x| < r2: "
           "Fi, Fe, their gradients (shape (3, ...)), divGiGradFi = div(Gi grad Fi), "
           "divGeGradFe = div(Ge grad Fe), Gi and Ge (shape (3, 3, ...)) and the ionic "
           "function g, each a function of (x, y, z, *PARAM_ORDER) with angles in radians; "
           "phi_i = sqrt(1 + t) Fi + c, phi_e = sqrt(1 + t) Fe + c, V_m = sqrt(1 + t) "
           "(Fi - Fe).  PARAMS holds the default and alternate parameters; "
           "field(name, x, y, z, **overrides) evaluates with them.  Generated by "
           "shell_bidomain_mms.py.")
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
    Exact fields Fi and Fe, their gradients, div(Gi grad Fi), div(Ge grad Fe),
    the conductivities Gi and Ge and the positive ionic function g of the
    heart-only bidomain manufactured solution on the spherical shell, sealed
    on both spheres, with self-test samples.

    Regenerate from generator/:
        python shell_bidomain_mms.py --write %s

    sympy %s
\\*---------------------------------------------------------------------------*/
"""


def generate_header():
    guard = "manufacturedShellBidomainReference_H"
    L = [LICENSE % (HEADER_PATH, sp.__version__)]
    L += ["#ifndef " + guard, "#define " + guard, ""]
    L += codegen.INCLUDES
    L += ["namespace Foam", "{", "namespace shellBidomainMMS", "{", ""]
    L += ["struct Params", "{", "    scalar r1, r2, gamma1, gamma2, k;",
          "    scalar gil, git, gin, gel, get, gen, theta0, theta1;",
          "    scalar fA[4];", "    scalar hA[3];",
          "    label solution;", "};", ""]
    L += codegen.params_function(_THIS, "defaultParams", default_params(), 0) + [""]
    L += codegen.params_function(_THIS, "alternateParams", alternate_params(), 0) + [""]
    for spec in FUNCTIONS:
        L += codegen.emit_function(_THIS, *spec) + [""]
    L += codegen.self_test_table(_rows_cache())
    L += ["} // End namespace shellBidomainMMS", "} // End namespace Foam", "",
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
    grp.add_argument("--emit", action="store_true")
    grp.add_argument("--self-check", action="store_true")
    args = ap.parse_args(argv)
    if args.self_check:
        try:
            run_self_checks()
        except AssertionError as err:
            print("self-check failed: %s" % err)
            return 1
        print("self-checks passed")
        return 0
    if args.emit:
        try:
            run_self_checks()
        except AssertionError as err:
            print("self-check failed, nothing emitted: %s" % err)
            return 1
        for rel in write_outputs():
            print("wrote %s" % rel)
        return 0
    if args.write:
        try:
            run_self_checks()
        except AssertionError as err:
            print("self-check failed, header not written: %s" % err)
            return 1
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
