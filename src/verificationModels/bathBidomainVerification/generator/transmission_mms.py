#!/usr/bin/env python3
"""Bath-bidomain transmission manufactured solution on the oblique slab.

A heart slab s1 < s < s2 lies between the bath slabs bathLow (s < s1) and
bathHigh (s > s2), with s = n.x, n = (cos alpha, sin alpha, 0), tau = t.x,
t = (-sin alpha, cos alpha, 0), and every field periodic in tau.  The exact
fields satisfy phiE continuity, extracellular flux continuity and intracellular
sealing on both interfaces.  The module builds the fields and sources with
sympy, runs the self-checks used by test_transmission_mms.py, and writes the
C++ reference header (namespace Foam::transmissionMMS) with CSE-compacted code
and mpmath self-test samples.

    python transmission_mms.py --write ../manufacturedTransmissionBathBidomainReference.H
    python transmission_mms.py --check ../manufacturedTransmissionBathBidomainReference.H
"""

import argparse
import copy
import functools
import math
import os
import re
import sys

import mpmath as mp
import numpy as np
import sympy as sp
from sympy.printing.c import C99CodePrinter

HEADER_NAME = "manufacturedTransmissionBathBidomainReference.H"
SAMPLE_DPS = 30

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

SCALAR_PARAMS = (
    "alpha", "s1", "s2", "period",
    "gil", "git", "gin", "gel", "get", "gen", "theta0", "theta1",
    "sigma0Low", "sigma0High", "sigmaEps",
)
TAIL_PARAMS = ("cmAmplitude", "cmOmega", "P0", "P1", "C")
ARRAY_PARAMS = (("fA", 4), ("hA", 3))
ANGLE_PARAMS = ("alpha", "theta0", "theta1")

PARAM_ORDER = (
    list(SCALAR_PARAMS)
    + ["fA%d" % i for i in range(4)]
    + ["hA%d" % i for i in range(3)]
    + list(TAIL_PARAMS)
)
P = {name: sp.Symbol(name, real=True) for name in PARAM_ORDER}

# Coordinates: x, y, z, t; the fields are built in the interface frame
# s = n.x, tau = t.x.  St = sqrt(1 + t), Tt = cmAmplitude sin(cmOmega t),
# thC = theta0 + theta1 clip(lambda, 0, 1).
x, y, z, t = sp.symbols("x y z t", real=True)
s_, tau_ = sp.symbols("s tau", real=True)
St_, Tt_, thC_ = sp.symbols("St Tt thC", real=True)
ARGS = [s_, tau_, t, St_, Tt_, thC_] + [P[n] for n in PARAM_ORDER]

ZONES = (0, 1, 2)
ZONE_NAMES = {0: "bathLow", 1: "myocardium", 2: "bathHigh"}
DOMAIN_S = (0.0, 2.0)  # slab extent in s (mesh)


def default_params():
    return {
        "alpha": 25.0, "s1": 0.5, "s2": 1.5, "period": 1.0,
        "gil": 0.17, "git": 0.019, "gin": 0.019,
        "gel": 0.62, "get": 0.12, "gen": 0.12,
        "theta0": -20.0, "theta1": 90.0,
        "sigma0Low": 0.30, "sigma0High": 0.45, "sigmaEps": 0.2,
        "fA": [0.3, 0.2, 0.4, 0.15],
        "hA": [-0.2, 0.3, 0.1],
        "cmAmplitude": 0.0, "cmOmega": 0.0,
        "P0": 0.1, "P1": 0.4, "C": 0.3,
    }


def alternate_params():
    p = default_params()
    p["cmAmplitude"] = 0.3
    p["cmOmega"] = 20.0
    p["fA"] = [0.25, 0.25, 0.35, 0.2]
    p["hA"] = [-0.15, 0.35, 0.15]
    return p


def param_values(params):
    """Flat {symbol name: float}; angles in radians as Foam::degToRad does."""
    v = {}
    for name in SCALAR_PARAMS + TAIL_PARAMS:
        val = float(params[name])
        if name in ANGLE_PARAMS:
            val = val * math.pi / 180.0
        v[name] = val
    for name, n in ARRAY_PARAMS:
        assert len(params[name]) == n, name
        for i in range(n):
            v["%s%d" % (name, i)] = float(params[name][i])
    return v


# ---------------------------------------------------------------------------
# Symbolic construction
# ---------------------------------------------------------------------------

def _grad(u):
    return sp.Matrix([sp.diff(u, s_), sp.diff(u, tau_), 0])


def _div(J):
    return sp.diff(J[0], s_) + sp.diff(J[1], tau_)


def _dot(a, b):
    return (a.T * b)[0]


def _frame_tensor(l, tr, nn, theta):
    """R(theta) diag(l, tr) R^T (+) nn expressed in the (s, tau, z) frame."""
    ph = theta - P["alpha"]
    c, sn = sp.cos(ph), sp.sin(ph)
    return sp.Matrix([
        [l*c**2 + tr*sn**2, (l - tr)*c*sn, 0],
        [(l - tr)*c*sn, l*sn**2 + tr*c**2, 0],
        [0, 0, nn],
    ])


def _xyz_tensor(l, tr, nn, theta):
    """R(theta) diag(l, tr) R^T (+) nn, components xx xy xz yy yz zz."""
    c, sn = sp.cos(theta), sp.sin(theta)
    return [l*c**2 + tr*sn**2, (l - tr)*c*sn, sp.Integer(0),
            l*sn**2 + tr*c**2, sp.Integer(0), nn]


def _to_xyz(v):
    a = P["alpha"]
    return [v[0]*sp.cos(a) - v[1]*sp.sin(a),
            v[0]*sp.sin(a) + v[1]*sp.cos(a),
            sp.Integer(0)]


def _tensors(theta):
    i = (P["gil"], P["git"], P["gin"])
    e = (P["gel"], P["get"], P["gen"])
    return _frame_tensor(*i, theta), _frame_tensor(*e, theta), \
        _xyz_tensor(*i, theta), _xyz_tensor(*e, theta)


def _ionic_G():
    k = 2*sp.pi/P["period"]
    return 1 + sp.Rational(1, 4)*sp.sin(k*tau_)*sp.cos(sp.pi*s_)


def _heart_exports(phiE, phiI, Gi, Ge):
    Vm = phiI - phiE
    G = _ionic_G()
    divI = _div(Gi*_grad(phiI))
    return {
        "phiE": phiE, "phiI": phiI, "Vm": Vm,
        "gradPhiE": _to_xyz(_grad(phiE)),
        "gradPhiI": _to_xyz(_grad(phiI)),
        "gradVm": _to_xyz(_grad(Vm)),
        "divGiGradPhiI": divI,
        "phiESource": divI + _div(Ge*_grad(phiE)),
        "J": _to_xyz(Ge*_grad(phiE)),
        "G": G,
        "u1": (1 + t)*G + Vm,
        "u2": 1/(sp.sqrt(G)*(1 + t)),
    }


def _bath_exports(phiE, sigma):
    return {
        "phiE": phiE,
        "gradPhiE": _to_xyz(_grad(phiE)),
        "phiESource": _div(sigma*_grad(phiE)),
        "sigmaB": sigma,
    }


@functools.lru_cache(maxsize=None)
def expressions(solution):
    """Per-zone expressions in (s, tau, t, St, Tt, thC, params)."""
    p = P
    s1, s2 = p["s1"], p["s2"]
    th0, th1 = p["theta0"], p["theta1"]
    E = {0: {}, 1: {}, 2: {}}
    if solution == "affine":
        Gi, Ge, Gix, Gex = _tensors(th0)
        nGen = Ge[0, 0]
        QL = nGen*p["P1"]/p["sigma0Low"]
        QH = nGen*p["P1"]/p["sigma0High"]
        phiEH = St_*(p["P0"] + p["P1"]*(s_ - s1))
        phiEL = St_*(p["P0"] + QL*(s_ - s1))
        phiER = St_*(p["P0"] + p["P1"]*(s2 - s1) + QH*(s_ - s2))
        phiI = St_*p["C"]
        E[1].update(_heart_exports(phiEH, phiI, Gi, Ge))
        E[0].update(_bath_exports(phiEL, p["sigma0Low"]))
        E[2].update(_bath_exports(phiER, p["sigma0High"]))
        for zn in ZONES:
            E[zn]["GiZ"], E[zn]["GeZ"] = Gix, Gex
        E["Gi"], E["Ge"] = Gix, Gex
        return E
    if solution != "levelSet":
        raise ValueError("solution must be 'affine' or 'levelSet'")

    k = 2*sp.pi/p["period"]
    lam = (s_ - s1)/(s2 - s1)
    om = (s_ - s1)*(s_ - s2)/(s2 - s1)
    thH = th0 + th1*lam
    GiH, GeH, GixH, GexH = _tensors(thH)
    GiL, GeL, GixL, GexL = _tensors(th0)
    GiR, GeR, GixR, GexR = _tensors(th0 + th1)
    sighat = sp.sin(k*tau_)*(1 + sp.Rational(1, 2)*sp.sin(2*sp.pi*s_))
    sig = {0: p["sigma0Low"]*(1 + p["sigmaEps"]*sighat),
           2: p["sigma0High"]*(1 + p["sigmaEps"]*sighat)}
    a = [p["fA%d" % i] for i in range(4)]
    b = [p["hA%d" % i] for i in range(3)]
    fhat = (a[0]*s_ + a[1]*s_**2
            + a[2]*sp.sin(k*tau_)*(1 + sp.Rational(1, 2)*s_)
            + a[3]*sp.cos(2*k*tau_ + sp.Rational(3, 10))*sp.exp(sp.Rational(2, 5)*s_))
    hhat = (b[0]*s_ + b[1]*sp.cos(k*tau_)*(1 - sp.Rational(3, 10)*s_**2)
            + b[2]*sp.sin(2*k*tau_)*s_)
    f, h = St_*fhat, St_*hhat
    f2 = sp.Rational(1, 5)*sp.cos(k*tau_) + sp.Rational(1, 10)*s_
    n = {1: sp.Matrix([-1, 0, 0]), 2: sp.Matrix([1, 0, 0])}
    sj = {1: s1, 2: s2}
    sigj = {1: sig[0], 2: sig[2]}
    d = {j: sp.simplify(_dot(n[j], _grad(om)).subs(s_, sj[j])) for j in (1, 2)}
    gj, mj, qj = {}, {}, {}
    for j in (1, 2):
        nj = n[j]
        gj[j] = _dot(nj, (sigj[j]*sp.eye(3) - GeH)*_grad(f))/(_dot(nj, GeH*nj)*d[j])
        mj[j] = -_dot(nj, GiH*_grad(h))/(_dot(nj, GiH*nj)*d[j])
        qj[j] = -_dot(nj, GiH*_grad(f2))/(_dot(nj, GiH*nj)*d[j])
    g = (1 - lam)*gj[1] + lam*gj[2]
    m = (1 - lam)*mj[1] + lam*mj[2]
    q = (1 - lam)*qj[1] + lam*qj[2]
    pcm = f2 + om*q
    phiEH = f + om*g + Tt_*pcm
    phiI = h + om*m + Tt_*pcm
    E[1].update(_heart_exports(phiEH, phiI, GiH, GeH))
    E[1].update({"g": g, "m": m, "q": q})

    def r_closed(nt, Gi, Ge, sigma, e):
        gf2 = _grad(f2)
        return (((_dot(nt, Ge*gf2) - _dot(nt, Gi*gf2)*_dot(nt, Ge*nt)/_dot(nt, Gi*nt))
                 / sigma - _dot(nt, gf2))/e)

    beta = {0: -(s_ - s1), 2: s_ - s2}
    for zn, j, Gi, Ge in ((0, 1, GiL, GeL), (2, 2, GiR, GeR)):
        e = sp.simplify(_dot(n[j], _grad(beta[zn])).subs(s_, sj[j]))
        r = r_closed(n[j], Gi, Ge, sig[zn], e)
        E[zn].update(_bath_exports(f + Tt_*(f2 + beta[zn]*r), sig[zn]))
        E[zn]["betaR"] = beta[zn]*r
    E[1]["GiZ"], E[1]["GeZ"] = GixH, GexH
    E[0]["GiZ"], E[0]["GeZ"] = GixL, GexL
    E[2]["GiZ"], E[2]["GeZ"] = GixR, GexR
    for zn in ZONES:
        E[zn].update({"fhat": fhat, "hhat": hhat, "f2": f2, "sighat": sighat})
    E["Gi"] = _xyz_tensor(P["gil"], P["git"], P["gin"], thC_)
    E["Ge"] = _xyz_tensor(P["gel"], P["get"], P["gen"], thC_)
    return E


def _expr(solution, zone, name):
    E = expressions(solution)
    if name in ("Gi", "Ge"):
        return E[name]
    return E[zone][name]


# ---------------------------------------------------------------------------
# Numerical evaluation
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _compiled(solution, zone, name, lib):
    e = _expr(solution, zone, name)
    if lib == "numpy":
        return sp.lambdify(ARGS, e, modules="numpy", cse=True)
    return sp.lambdify(ARGS, e, modules="mpmath", cse=True)


def _context(sv, tv, tt, pv, lib):
    """Values of the preamble symbols, as the C++ functions compute them."""
    if lib == "numpy":
        sqrt, sin = np.sqrt, np.sin
        lam = (sv - pv["s1"])/(pv["s2"] - pv["s1"])
        clip = np.clip(lam, 0.0, 1.0)
    else:
        sqrt, sin = mp.sqrt, mp.sin
        lam = (sv - pv["s1"])/(pv["s2"] - pv["s1"])
        clip = min(max(lam, 0), 1)
    St = sqrt(1 + tt)
    Tt = pv["cmAmplitude"]*sin(pv["cmOmega"]*tt)
    thC = pv["theta0"] + pv["theta1"]*clip
    return [sv, tv, tt, St, Tt, thC] + [pv[n] for n in PARAM_ORDER]


def _st_from_xyz(pt, alpha, lib):
    if lib == "numpy":
        pt = np.asarray(pt, dtype=float)
        ca, sa = np.cos(alpha), np.sin(alpha)
        return ca*pt[..., 0] + sa*pt[..., 1], -sa*pt[..., 0] + ca*pt[..., 1]
    ca, sa = mp.cos(alpha), mp.sin(alpha)
    return ca*pt[0] + sa*pt[1], -sa*pt[0] + ca*pt[1]


def _full(T6):
    """symmTensor components (xx xy xz yy yz zz) as 3x3 matrices, stacked last."""
    return np.stack([np.stack([T6[..., 0], T6[..., 1], T6[..., 2]], -1),
                     np.stack([T6[..., 1], T6[..., 3], T6[..., 4]], -1),
                     np.stack([T6[..., 2], T6[..., 4], T6[..., 5]], -1)], -2)


class Model:
    def __init__(self, solution, params):
        self.solution = solution
        self.params = copy.deepcopy(params)
        self.pv = param_values(params)
        self.E = expressions(solution)
        a = self.pv["alpha"]
        self.nhat = np.array([math.cos(a), math.sin(a), 0.0])
        self.that = np.array([-math.sin(a), math.cos(a), 0.0])

    # -- evaluation helpers -------------------------------------------------
    def eval_st(self, zone, name, sv, tv, tt):
        """numpy evaluation on (s, tau) arrays; vectors/tensors stacked last."""
        sv = np.asarray(sv, dtype=float)
        tv = np.asarray(tv, dtype=float)
        sv, tv = np.broadcast_arrays(sv, tv)
        out = _compiled(self.solution, zone, name, "numpy")(
            *_context(sv, tv, float(tt), self.pv, "numpy"))
        if isinstance(out, (list, tuple)):
            return np.stack([np.broadcast_to(np.asarray(c, dtype=float), sv.shape)
                             for c in out], axis=-1)
        return np.broadcast_to(np.asarray(out, dtype=float), sv.shape)

    def eval_xyz(self, zone, name, pts, tt):
        sv, tv = _st_from_xyz(pts, self.pv["alpha"], "numpy")
        return self.eval_st(zone, name, sv, tv, tt)

    def eval_mp(self, zone, name, pt, tt, pv=None):
        """mpmath evaluation at a point (x, y, z) given as floats."""
        pv = {k: mp.mpf(v) for k, v in (pv or self.pv).items()}
        ptm = [mp.mpf(c) for c in pt]
        sv, tv = _st_from_xyz(ptm, pv["alpha"], "mpmath")
        out = _compiled(self.solution, zone, name, "mpmath")(
            *_context(sv, tv, mp.mpf(tt), pv, "mpmath"))
        if isinstance(out, (list, tuple)):
            return [mp.mpf(c) for c in out]
        return mp.mpf(out)

    def zone_of(self, pt):
        sv = float(np.dot(self.nhat, pt))
        if sv < self.pv["s1"]:
            return 0
        if sv > self.pv["s2"]:
            return 2
        return 1

    def point(self, sv, tv, zv=0.05):
        return sv*self.nhat + tv*self.that + np.array([0.0, 0.0, zv])

    # -- interface conditions -----------------------------------------------
    def interface_residuals(self, j, tau, tt):
        mp.mp.dps = SAMPLE_DPS
        zb = 0 if j == 1 else 2
        sj = mp.mpf(self.pv["s1"] if j == 1 else self.pv["s2"])
        a = mp.mpf(self.pv["alpha"])
        nh = [mp.cos(a), mp.sin(a), mp.mpf(0)]
        th = [-mp.sin(a), mp.cos(a), mp.mpf(0)]
        tv = mp.mpf(float(tau))
        pt = [sj*nh[i] + tv*th[i] for i in range(2)] + [mp.mpf("0.05")]
        n = [-c for c in nh] if j == 1 else nh

        def tdot(T6, v):
            M = [[T6[0], T6[1], T6[2]], [T6[1], T6[3], T6[4]], [T6[2], T6[4], T6[5]]]
            return [sum(M[i][k]*v[k] for k in range(3)) for i in range(3)]

        def vdot(a, b):
            return sum(a[i]*b[i] for i in range(3))

        Ge = self.eval_mp(1, "GeZ", pt, tt)
        Gi = self.eval_mp(1, "GiZ", pt, tt)
        cont = self.eval_mp(1, "phiE", pt, tt) - self.eval_mp(zb, "phiE", pt, tt)
        flux = (vdot(n, tdot(Ge, self.eval_mp(1, "gradPhiE", pt, tt)))
                - self.eval_mp(zb, "sigmaB", pt, tt)*vdot(n, self.eval_mp(zb, "gradPhiE", pt, tt)))
        seal = vdot(n, tdot(Gi, self.eval_mp(1, "gradPhiI", pt, tt)))
        return {"continuity": float(cont), "extracellularFlux": float(flux),
                "sealing": float(seal)}

    # -- sources versus finite differences ------------------------------------
    def sample_points(self, seed, n, min_distance_to_interface):
        rng = np.random.default_rng(seed)
        pts = []
        while len(pts) < n:
            sv = rng.uniform(*DOMAIN_S)
            if min(abs(sv - self.pv["s1"]), abs(sv - self.pv["s2"])) < min_distance_to_interface:
                continue
            pts.append(self.point(sv, rng.uniform(0.0, 1.0)))
        return pts

    def source_fd_error(self, pt, tt, h):
        """Max relative error of the exact gradients/sources against 4th-order
        central differences of the fields; the stencil uses the smooth formula
        of the zone that contains pt."""
        zone = self.zone_of(pt)
        pt = np.asarray(pt, dtype=float)
        D = ((-2, 1.0/12.0), (-1, -8.0/12.0), (1, 8.0/12.0), (2, -1.0/12.0))
        outer = [(ax, m, c) for ax in range(3) for m, c in D]
        Q = np.array([pt + m*h*np.eye(3)[ax] for ax, m, c in outer])

        def grad_fd(name, X):
            g = np.zeros((len(X), 3))
            for ax in range(3):
                for m, c in D:
                    g[:, ax] += c*self.eval_xyz(zone, name, X + m*h*np.eye(3)[ax], tt)/h
            return g

        def div_fd(Tname, name):
            if Tname == "sigma":
                T = self.eval_xyz(zone, "sigmaB", Q, tt)[:, None, None]*np.eye(3)
            else:
                T = _full(self.eval_xyz(zone, Tname, Q, tt))
            flux = np.einsum("nij,nj->ni", T, grad_fd(name, Q))
            return sum(c*flux[i, ax]/h for i, (ax, m, c) in enumerate(outer))

        def rel(exact, approx):
            exact, approx = np.atleast_1d(exact), np.atleast_1d(approx)
            return float(np.max(np.abs(exact - approx)/np.maximum(1.0, np.abs(exact))))

        P0 = pt[None, :]
        ex = lambda name: self.eval_xyz(zone, name, P0, tt)[0]  # noqa: E731
        errs = [rel(ex("gradPhiE"), grad_fd("phiE", P0)[0])]
        # the clipped tensors equal the zone formula inside each zone
        errs.append(rel(self.eval_xyz(zone, "Gi", P0, tt)[0], ex("GiZ")))
        errs.append(rel(self.eval_xyz(zone, "Ge", P0, tt)[0], ex("GeZ")))
        if zone == 1:
            divI = div_fd("GiZ", "phiI")
            errs.append(rel(ex("divGiGradPhiI"), divI))
            errs.append(rel(ex("phiESource"), divI + div_fd("GeZ", "phiE")))
            errs.append(rel(ex("gradVm"), grad_fd("Vm", P0)[0]))
            errs.append(rel(ex("J"), _full(ex("GeZ")) @ grad_fd("phiE", P0)[0]))
        else:
            errs.append(rel(ex("phiESource"), div_fd("sigma", "phiE")))
        return max(errs)

    # -- boundedness, periodicity, richness -----------------------------------
    def _names(self, zone):
        if zone == 1:
            names = ["phiE", "phiI", "Vm", "gradPhiE", "gradPhiI", "gradVm",
                     "divGiGradPhiI", "phiESource", "J", "G", "u1", "u2", "GiZ", "GeZ"]
            aux = ["g", "m", "q"]
        else:
            names = ["phiE", "gradPhiE", "phiESource", "sigmaB", "GiZ", "GeZ"]
            aux = ["betaR"]
        return names + [a for a in aux if a in self.E[zone]]

    def _zone_s_range(self, zone):
        return {0: (DOMAIN_S[0], self.pv["s1"]), 1: (self.pv["s1"], self.pv["s2"]),
                2: (self.pv["s2"], DOMAIN_S[1])}[zone]

    def dense_grid_values(self, n=200, tt=0.2):
        vals = []
        for zone in ZONES:
            lo, hi = self._zone_s_range(zone)
            S, TA = np.meshgrid(np.linspace(lo, hi, n), np.linspace(0.0, self.pv["period"], n))
            for name in self._names(zone):
                vals.append(np.ravel(self.eval_st(zone, name, S, TA, tt)))
        return np.concatenate(vals)

    def periodicity_error(self):
        err = 0.0
        L = self.pv["period"]
        for zone in ZONES:
            lo, hi = self._zone_s_range(zone)
            S, TA = np.meshgrid(np.linspace(lo, hi, 9)[1:-1], np.linspace(0.0, L, 11)[:-1])
            names = self._names(zone) + [nm for nm in ("fhat", "hhat", "f2", "sighat")
                                         if nm in self.E[zone]]
            for tt in (0.0, 0.13):
                for name in names:
                    a = self.eval_st(zone, name, S, TA, tt)
                    b = self.eval_st(zone, name, S, TA + L, tt)
                    err = max(err, float(np.max(np.abs(a - b))))
        return err

    def richness(self, nq=400, tt=0.0):
        """Per interface (Gamma1, Gamma2) and per tensor (Ge, Gi): RMS over a
        uniform Gamma quadrature of |(G n)_t . grad phi| / |n.Ge grad phiE|
        (Ge) or / |Gi grad phiI| (Gi), and the fraction with ratio >= 0.2."""
        rms, frac = [], []
        taus = (np.arange(nq) + 0.5)/nq*self.pv["period"]
        for j in (1, 2):
            sj = self.pv["s1"] if j == 1 else self.pv["s2"]
            n = -self.nhat if j == 1 else self.nhat
            pts = np.array([self.point(sj, tv) for tv in taus])

            Ge = _full(self.eval_xyz(1, "GeZ", pts, tt))
            Gi = _full(self.eval_xyz(1, "GiZ", pts, tt))
            ge = self.eval_xyz(1, "gradPhiE", pts, tt)
            gi = self.eval_xyz(1, "gradPhiI", pts, tt)
            for G, gr, kind in ((Ge, ge, "e"), (Gi, gi, "i")):
                Gn = G @ n
                Gn_t = Gn - (Gn @ n)[:, None]*n[None, :]
                num = np.abs(np.einsum("ni,ni->n", Gn_t, gr))
                Gg = np.einsum("nij,nj->ni", G, gr)
                den = np.abs(Gg @ n) if kind == "e" else np.linalg.norm(Gg, axis=1)
                ratio = num/den
                rms.append(float(np.sqrt(np.mean(ratio**2))))
                frac.append(float(np.mean(ratio >= 0.2)))
        return rms, frac

    # -- affine checks ----------------------------------------------------------
    def max_second_derivative(self):
        worst = 0.0
        for zone in ZONES:
            for name in ("phiE", "phiI", "Vm"):
                if name not in self.E[zone]:
                    continue
                e = self.E[zone][name]
                for a in (s_, tau_):
                    for b in (s_, tau_):
                        d2 = sp.simplify(sp.diff(e, a, b))
                        if d2 != 0:
                            f = sp.lambdify(ARGS, d2, "numpy")
                            for sv in np.linspace(*self._zone_s_range(zone), 5):
                                worst = max(worst, abs(float(f(*_context(
                                    sv, 0.3, 0.1, self.pv, "numpy")))))
        return worst

    def tensors_constant(self):
        space = {s_, tau_, thC_}
        exprs = list(self.E["Gi"]) + list(self.E["Ge"])
        for zone in ZONES:
            exprs += list(self.E[zone]["GiZ"]) + list(self.E[zone]["GeZ"])
            if "sigmaB" in self.E[zone]:
                exprs.append(self.E[zone]["sigmaB"])
        return all(not (sp.sympify(e).free_symbols & space) for e in exprs)

    # -- ionic closure ----------------------------------------------------------
    def ionic_residual(self):
        """ODE closure of the FDA ionic model (Cm = 2, Chi = 3) with the exact
        states, and the Vm-equation balance with the generated source, at
        chi*cm != Chi*Cm."""
        mp.mp.dps = SAMPLE_DPS
        Cm, Chi = 2, 3
        chi, cm = sp.Rational(19, 10), sp.Rational(7, 10)
        Beta = -sp.Rational(11, 10)*(1 - 1/sp.sqrt(2))
        sub = {St_: sp.sqrt(1 + t), Tt_: P["cmAmplitude"]*sp.sin(P["cmOmega"]*t)}
        V = self.E[1]["Vm"].subs(sub)
        u1 = self.E[1]["u1"].subs(sub)
        u2 = self.E[1]["u2"].subs(sub)
        u3 = 0
        divI = self.E[1]["divGiGradPhiI"].subs(sub)
        w = u1 + u3 - V
        rate1 = w**2*u2**2 + sp.Rational(1, 2)*w*u2**2*(V - u3)
        rate2 = -w*u2**3
        Iion = -sp.Rational(1, 2)*Cm*w*u2**2*(V - u3) + Beta*(V - u3)/Chi
        Im = Iion/Cm
        S_Vm = (chi*cm)/(Chi*Cm)*Beta*V - divI
        res = [sp.diff(u1, t) - rate1,
               sp.diff(u2, t) - rate2,
               Cm*sp.diff(V, t) + Iion - Beta*(V - u3)/Chi,
               chi*cm*sp.diff(V, t) - (divI - chi*cm*Im + S_Vm)]
        f = sp.lambdify(ARGS, res, modules="mpmath")
        pv = {k: mp.mpf(v) for k, v in self.pv.items()}
        worst = mp.mpf(0)
        for pt in self.sample_points(seed=7, n=60, min_distance_to_interface=0.0):
            if self.zone_of(pt) != 1:
                continue
            ptm = [mp.mpf(c) for c in pt]
            sv, tv = _st_from_xyz(ptm, pv["alpha"], "mpmath")
            for tt in (0.0, 0.07, 0.13, 0.5):
                for r in f(*_context(sv, tv, mp.mpf(tt), pv, "mpmath")):
                    worst = max(worst, abs(r))
        return float(worst)


def build(solution, params):
    return Model(solution, params)


def vm_common_mode_difference(params_default, params_alternate):
    """Vm, grad Vm, u1, u2 with and without the common mode: the default set
    with the alternate common mode switched on, and the alternate set with it
    switched off."""
    mp.mp.dps = SAMPLE_DPS
    on = copy.deepcopy(params_default)
    on["cmAmplitude"] = params_alternate["cmAmplitude"]
    on["cmOmega"] = params_alternate["cmOmega"]
    off = copy.deepcopy(params_alternate)
    off["cmAmplitude"] = 0.0
    worst = 0.0
    for pa, pb in ((params_default, on), (params_alternate, off)):
        ma, mb = build("levelSet", pa), build("levelSet", pb)
        for pt in ma.sample_points(seed=3, n=30, min_distance_to_interface=0.0):
            for tt in (0.0, 0.07, 0.137):
                for name in ("Vm", "gradVm", "u1", "u2"):
                    a, b = ma.eval_mp(1, name, pt, tt), mb.eval_mp(1, name, pt, tt)
                    a, b = (a, b) if isinstance(a, list) else ([a], [b])
                    worst = max(worst, max(float(abs(u - v)) for u, v in zip(a, b)))
    return worst


# ---------------------------------------------------------------------------
# Self-checks (run before any header is written; also used by the tests)
# ---------------------------------------------------------------------------

def _require(cond, *info):
    if not cond:
        raise AssertionError(" ".join(str(i) for i in info) or "check failed")


def check_interfaces(m):
    for j in (1, 2):
        for tau in np.linspace(0.0, 1.0, 13):
            for tt in (0.0, 0.07, 0.2):
                r = m.interface_residuals(j, tau, tt)
                _require(all(abs(v) < 1e-12 for v in r.values()),
                         "interface", m.solution, j, tau, tt, r)


def check_all_interfaces():
    for sol in ("affine", "levelSet"):
        check_interfaces(build(sol, default_params()))
    check_interfaces(build("levelSet", alternate_params()))


def check_sources_match_fd():
    for sol, p in (("affine", default_params()), ("levelSet", default_params()),
                   ("levelSet", alternate_params())):
        m = build(sol, p)
        for pt in m.sample_points(seed=1, n=40, min_distance_to_interface=2e-3):
            for tt in (0.0, 0.13):
                _require(m.source_fd_error(pt, tt, h=1e-3) < 1e-6,
                         "source/gradient FD", sol, pt, tt)


def check_bounded_and_finite():
    v = build("levelSet", alternate_params()).dense_grid_values()
    _require(np.all(np.isfinite(v)) and np.max(np.abs(v)) < 1e3, "bounded/finite")


def check_periodic():
    err = build("levelSet", alternate_params()).periodicity_error()
    _require(err < 1e-12, "periodicity", err)


def check_richness():
    rms, frac = build("levelSet", default_params()).richness()
    _require(min(rms) >= 0.2 and min(frac) >= 0.8, "richness", rms, frac)


def check_common_mode_keeps_vm():
    d = vm_common_mode_difference(default_params(), alternate_params())
    _require(d < 1e-14, "common mode changes Vm", d)


def check_affine():
    m = build("affine", default_params())
    _require(m.max_second_derivative() == 0 and m.tensors_constant(), "affine")


def check_ionic_closure():
    r = build("levelSet", default_params()).ionic_residual()
    _require(r < 1e-13, "ionic closure", r)


SELF_CHECKS = (check_all_interfaces, check_sources_match_fd, check_bounded_and_finite,
               check_periodic, check_richness, check_common_mode_keeps_vm,
               check_affine, check_ionic_closure)


def run_self_checks():
    """Every generator check except the header comparison; raises AssertionError."""
    for check in SELF_CHECKS:
        check()


# ---------------------------------------------------------------------------
# C++ emission
# ---------------------------------------------------------------------------

class _CxxPrinter(C99CodePrinter):
    """C++ for OpenFOAM scalars: no integer-exponent pow(), pi from Foam."""

    def _print_Pi(self, expr):
        return "constant::mathematical::pi"

    def _print_Pow(self, expr):
        # Every result is a single parenthesised or call operand, so it is
        # safe after '/' in the Mul printer.
        b, e = expr.base, expr.exp
        bs = self._print(b) if b.is_Symbol else "(%s)" % self._print(b)
        if e.is_Integer and 1 <= abs(int(e)) <= 4:
            prod = "*".join([bs]*abs(int(e)))
            return "(%s)" % prod if e > 0 else "(1.0/(%s))" % prod
        if e == sp.Rational(1, 2):
            return "sqrt(%s)" % self._print(b)
        if e == -sp.Rational(1, 2):
            return "(1.0/sqrt(%s))" % self._print(b)
        return "pow(%s, %s)" % (self._print(b), self._print(sp.Float(e, 17)))


_PRINTER = _CxxPrinter()


def _cxx(expr):
    return _PRINTER.doprint(expr)


# fn index, C++ name, kind, expression name, value kind
FUNCTIONS = (
    (0, "Vm", "heart", "Vm", "scalar"),
    (1, "phiE", "zone", "phiE", "scalar"),
    (2, "phiI", "heart", "phiI", "scalar"),
    (3, "gradPhiE", "zone", "gradPhiE", "vector"),
    (4, "gradVm", "heart", "gradVm", "vector"),
    (5, "u1", "heart", "u1", "scalar"),
    (6, "u2", "heart", "u2", "scalar"),
    (7, "Gi", "tensor", "Gi", "symmTensor"),
    (8, "Ge", "tensor", "Ge", "symmTensor"),
    (9, "sigmaB", "bath", "sigmaB", "scalar"),
    (10, "divGiGradPhiI", "heart", "divGiGradPhiI", "scalar"),
    (11, "phiESource", "zone", "phiESource", "scalar"),
    (12, "interfaceFlux", "flux", "J", "flux"),
    (13, "G", "heartNoTime", "G", "scalar"),
)

SIGNATURES = {
    "heart": "const point& pt, const scalar t, const Params& p",
    "heartNoTime": "const point& pt, const Params& p",
    "flux": "const point& pt, const vector& n, const scalar t, const Params& p",
    "zone": "const point& pt, const scalar t, const label zone, const Params& p",
    "bath": "const point& pt, const label zone, const Params& p",
    "tensor": "const point& pt, const Params& p",
}
RETURN_TYPE = {"scalar": "scalar", "vector": "vector", "symmTensor": "symmTensor",
               "flux": "scalar"}
ZONES_OF_KIND = {"heart": (1,), "heartNoTime": (1,), "flux": (1,), "zone": (0, 1, 2),
                 "bath": (0, 2), "tensor": (0, 1, 2)}
N_COMPONENTS = {"scalar": 1, "vector": 3, "symmTensor": 6, "flux": 3}


def _as_list(e):
    return list(e) if isinstance(e, (list, tuple, sp.MatrixBase)) else [e]


def _check_printed(exprs, repl, red):
    """Evaluate the printed C++ arithmetic as Python at one point and compare
    with sympy, so a printer precedence error cannot reach the header."""
    pv = param_values(alternate_params())
    env = {"sin": math.sin, "cos": math.cos, "exp": math.exp, "sqrt": math.sqrt,
           "pow": math.pow}
    ctx = dict(zip([str(a) for a in ARGS], _context(0.83, 0.29, 0.11, pv, "numpy")))
    ctx = {k: float(v) for k, v in ctx.items()}
    env.update(ctx)

    def py(code):
        if re.search(r"(?<![\w.])\d+/\d+(?![\w.])", code):
            raise RuntimeError("integer division in emitted code: " + code)
        return eval(code.replace("constant::mathematical::pi", "%r" % math.pi), {}, env)

    for sym, val in repl:
        env[str(sym)] = py(_cxx(val))
    ref = sp.lambdify(ARGS, exprs, modules="math")(*[ctx[str(a)] for a in ARGS])
    for e, r in zip(red, ref):
        v = py(_cxx(e))
        if abs(v - r) > 1e-12*max(1.0, abs(r)):
            raise RuntimeError("printed C++ differs from sympy: %r vs %r" % (v, r))


def _branch_code(exprs, vkind, indent):
    exprs = [sp.sympify(e) for e in exprs]
    repl, red = sp.cse(exprs, symbols=sp.numbered_symbols("w"), order="canonical")
    _check_printed(exprs, repl, red)
    pad = " "*indent
    lines = ["%sconst scalar %s = %s;" % (pad, _cxx(sym), _cxx(val)) for sym, val in repl]
    r = [_cxx(e) for e in red]
    if vkind == "scalar":
        lines.append("%sreturn %s;" % (pad, r[0]))
    elif vkind == "vector":
        lines.append("%sreturn vector(%s, %s, %s);" % (pad, *r))
    elif vkind == "symmTensor":
        lines.append("%sreturn symmTensor(%s, %s, %s, %s, %s, %s);" % (pad, *r))
    else:
        lines.append("%sreturn n.x()*(%s) + n.y()*(%s) + n.z()*(%s);" % (pad, *r))
    return lines


PREAMBLE = (
    (s_, ("alpha",), "const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();"),
    (tau_, ("alpha",), "const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();"),
    (St_, (), "const scalar St = sqrt(1 + t);"),
    (Tt_, ("cmAmplitude", "cmOmega"), "const scalar Tt = cmAmplitude*sin(cmOmega*t);"),
    (thC_, ("theta0", "theta1", "s1", "s2"),
     "const scalar thC = theta0 + theta1*min(max((s - s1)/(s2 - s1), scalar(0)), scalar(1));"),
)


# preamble locals that other preamble locals are computed from
PREAMBLE_NEEDS = {thC_: (s_,)}


def _param_access(name):
    for arr, _ in ARRAY_PARAMS:
        if name.startswith(arr) and name[len(arr):].isdigit():
            return "p.%s[%s]" % (arr, name[len(arr):])
    return "p." + name


def _emit_function(idx, cname, kind, ename, vkind):
    zones = ZONES_OF_KIND[kind]
    branches = {}
    used = set()
    for sol in ("affine", "levelSet"):
        for zn in zones:
            ex = _as_list(_expr(sol, zn, ename))
            branches[(sol, zn)] = ex
            for e in ex:
                used |= sp.sympify(e).free_symbols
    for sym, needs in PREAMBLE_NEEDS.items():
        if sym in used:
            used |= set(needs)
    pre = [(sym, deps, code) for sym, deps, code in PREAMBLE if sym in used]
    params = [n for n in PARAM_ORDER
              if P[n] in used or any(n in deps for _, deps, _ in pre)]
    out = ["inline %s %s(%s)" % (RETURN_TYPE[vkind], cname, SIGNATURES[kind]), "{"]
    out += ["    const scalar %s = %s;" % (n, _param_access(n)) for n in params]
    out += ["    " + code for _, _, code in pre]
    if kind in ("zone", "bath"):
        out.append("    checkZone(zone, %s);" % ("true" if kind == "zone" else "false"))
    out.append("")

    def zone_block(sol, indent):
        lines = []
        if len(zones) == 1 or kind == "tensor":
            return _branch_code(branches[(sol, zones[0])], vkind, indent)
        pad = " "*indent
        for zn in zones:
            if zn == 1:
                continue
            lines.append("%sif (zone == %s)" % (pad, ZONE_NAMES[zn]))
            lines.append(pad + "{")
            lines += _branch_code(branches[(sol, zn)], vkind, indent + 4)
            lines.append(pad + "}")
        if 1 in zones:
            lines += _branch_code(branches[(sol, 1)], vkind, indent)
        else:
            lines.append("%sreturn 0;" % pad)
        return lines

    out.append("    if (p.solution == 0)")
    out.append("    {")
    out += zone_block("affine", 8)
    out.append("    }")
    out.append("")
    out += zone_block("levelSet", 4)
    out.append("}")
    return out


def _sample_points_st():
    return {
        0: [(0.07, 0.13), (0.31, 0.62), (0.4973, -0.41)],
        1: [(0.5031, 0.27), (0.94, 0.81), (1.4962, -0.18)],
        2: [(1.5027, 0.44), (1.71, -0.07), (1.93, 0.71)],
    }


SAMPLE_TIMES = (0.0, 0.137)
# functions without a time argument are sampled once, at the first time
TIME_INDEPENDENT_KINDS = ("heartNoTime", "bath", "tensor")


def self_test_samples():
    mp.mp.dps = SAMPLE_DPS
    sets = ((0, "levelSet", default_params()), (1, "levelSet", alternate_params()),
            (0, "affine", default_params()))
    rows = []
    for idx, cname, kind, ename, vkind in FUNCTIONS:
        for param, sol, params in sets:
            m = build(sol, params)
            a = mp.mpf(m.pv["alpha"])
            for zn in ZONES_OF_KIND[kind]:
                for sv, tv in _sample_points_st()[zn]:
                    pt = (float(sv*mp.cos(a) - tv*mp.sin(a)),
                          float(sv*mp.sin(a) + tv*mp.cos(a)), 0.05)
                    for tt in (SAMPLE_TIMES if kind not in TIME_INDEPENDENT_KINDS
                               else SAMPLE_TIMES[:1]):
                        val = m.eval_mp(zn, ename, pt, tt)
                        vals = val if isinstance(val, list) else [val]
                        for comp in range(N_COMPONENTS[vkind]):
                            rows.append((idx, param, 0 if sol == "affine" else 1, zn, comp,
                                         pt[0], pt[1], pt[2], tt, float(vals[comp])))
    return rows


def _g(v):
    return "%.17g" % v


def _params_function(name, params):
    v = param_values(params)
    out = ["inline Params %s()" % name, "{", "    Params p;"]
    for n in SCALAR_PARAMS + TAIL_PARAMS:
        if n in ANGLE_PARAMS:
            out.append("    p.%s = degToRad(%r);" % (n, float(params[n])))
        else:
            out.append("    p.%s = %r;" % (n, v[n]))
    for arr, cnt in ARRAY_PARAMS:
        for i in range(cnt):
            out.append("    p.%s[%d] = %r;" % (arr, i, float(params[arr][i])))
    out.append("    p.solution = 1;")
    out += ["    return p;", "}"]
    return out


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
    Exact fields, conductivities, sources and self-test samples of the
    bath-bidomain transmission manufactured solution on the oblique slab.

    Regenerate from generator/:
        python transmission_mms.py --write ../%s

    sympy %s
\\*---------------------------------------------------------------------------*/
"""


def generate_header():
    guard = "manufacturedTransmissionBathBidomainReference_H"
    L = [LICENSE % (HEADER_NAME, sp.__version__)]
    L += ["#ifndef " + guard, "#define " + guard, ""]
    L += ['#include "error.H"', '#include "label.H"', '#include "mathematicalConstants.H"',
          '#include "point.H"', '#include "scalar.H"', '#include "symmTensor.H"',
          '#include "unitConversion.H"', '#include "vector.H"', ""]
    L += ["namespace Foam", "{", "namespace transmissionMMS", "{", ""]
    L += ["struct Params", "{", "    scalar alpha;", "    scalar s1, s2, period;",
          "    scalar gil, git, gin, gel, get, gen, theta0, theta1;",
          "    scalar sigma0Low, sigma0High, sigmaEps;",
          "    scalar fA[4], hA[3];", "    scalar cmAmplitude, cmOmega;",
          "    scalar P0, P1, C;", "    label solution;", "};", ""]
    L += ["enum Zone : label { bathLow = 0, myocardium = 1, bathHigh = 2 };", ""]
    L += _params_function("defaultParams", default_params()) + [""]
    L += _params_function("alternateParams", alternate_params()) + [""]
    L += ["inline label zoneOf(const point& pt, const Params& p)", "{",
          "    const scalar s = cos(p.alpha)*pt.x() + sin(p.alpha)*pt.y();",
          "    if (s < p.s1)", "    {", "        return bathLow;", "    }",
          "    if (s > p.s2)", "    {", "        return bathHigh;", "    }",
          "    return myocardium;", "}", ""]
    L += ["inline void checkZone(const label zone, const bool heartAllowed)", "{",
          "    if (zone < bathLow || zone > bathHigh || (!heartAllowed && zone == myocardium))",
          "    {", "        FatalErrorInFunction",
          "            << \"Zone \" << zone << \" is not valid here\" << exit(FatalError);",
          "    }", "}", ""]
    for spec in FUNCTIONS:
        L += _emit_function(*spec) + [""]
    rows = _rows_cache()
    L += ["struct SelfTestSample", "{", "    label fn;", "    label param;",
          "    label solution;", "    label zone;", "    label component;",
          "    scalar x, y, z, t, value;", "};", ""]
    L.append("inline const SelfTestSample selfTestSamples[] =")
    L.append("{")
    for r in rows:
        L.append("    {%d, %d, %d, %d, %d, %s, %s, %s, %s, %s}," % (
            r[0], r[1], r[2], r[3], r[4], _g(r[5]), _g(r[6]), _g(r[7]), _g(r[8]), _g(r[9])))
    L.append("};")
    L += ["", "inline const label nSelfTestSamples = %d;" % len(rows), ""]
    L += ["} // End namespace transmissionMMS", "} // End namespace Foam", "",
          "#endif", ""]
    return "\n".join(L)


def _resolve(path):
    """Relative paths are taken against this directory, for --write and --check
    alike, so the result does not depend on the working directory."""
    if os.path.isabs(path):
        return path
    return os.path.normpath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), path))


def header_matches(path, text=None):
    path = _resolve(path)
    if not os.path.exists(path):
        return False
    with open(path) as fh:
        current = fh.read()
    return current == (text if text is not None else generate_header())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--write", metavar="PATH")
    g.add_argument("--check", metavar="PATH")
    args = ap.parse_args(argv)
    if args.write:
        try:
            run_self_checks()
        except AssertionError as err:
            print("self-check failed, header not written: %s" % err)
            return 1
        text = generate_header()
        path = _resolve(args.write)
        with open(path, "w") as fh:
            fh.write(text)
        print("wrote %s (%d bytes, %d self-test samples)"
              % (path, len(text.encode()), len(_rows_cache())))
        return 0
    if header_matches(args.check):
        print("%s is up to date" % _resolve(args.check))
        return 0
    print("%s differs from the generator output" % _resolve(args.check))
    return 1


@functools.lru_cache(maxsize=None)
def _rows_cache():
    return tuple(self_test_samples())


if __name__ == "__main__":
    sys.exit(main())
