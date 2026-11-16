#!/usr/bin/env python3
"""Bath bidomain manufactured solution on the spherical shell (the bath shell D).

The heart shell r1 < |x| < r2 lies between the bath ball bloodPool (|x| < r1)
and the outer bath r2 < |x| < r3; bloodPool and bath are one extracellular
conductor with the heart, coupled at every step (torsoCoupling bath).  The
geometry, the fibre frame, the tensors Gi and Ge (unequal eigenvalue ratios,
oblique to both spheres), omega, the base fields f = S fhat and h = S hhat,
the bath conductivity, f2 and the ionic closure are those of
shell_bidomain_mms.py and shell_mms.py; what this generator adds is the
coupling of the bath shell:

    heart     phi_e = f + omega g + T p,   phi_i = h + omega m + T p,
              p = f2 + omega q,            Vm = phi_i - phi_e
    bath      phi_e = Phi + (outer closure, bath only),
              Phi = f + T (f2 + beta_j r_j)

with S = sqrt(1 + t) and the common mode T = cmAmplitude sin(cmOmega t).
On both spheres the potential is continuous, the extracellular current is
continuous through Ge on the heart side (n.Ge grad phi_e = sigma_b n.grad
phi_e: the bath rule of passiveConductor), and the intracellular space is
sealed (n.Gi grad phi_i = 0), so n.Gi grad Vm = -n.Gi grad phi_e is the
virtual-electrode current.  The corrections are

    g_j = n_j.(sigma_b,j I - Ge) grad f/((n_j.Ge n_j) d_j)
    m_j = -n_j.Gi grad h/((n_j.Gi n_j) d_j),  q_j likewise with f2
    r_j = n_j.Ge grad p/sigma_b - n_j.grad f2

blended through the wall as (1 - lambda) X_1 + lambda X_2.  The normals
n_1 = -x/r1, n_2 = x/r2 are the unit normals on their own sphere and are
polynomial, so every correction is smooth at the origin.  The outer sphere r3
is closed with constant data only, as PR 62 imposes:

    insulated (solution 0): phi_e = Phi + eta3 c3,
        eta3 = (r - r2)^2 (r - r3)/(r3 - r2)^2, c3 = -xhat.grad Phi,
        so d phi_e/dr = 0 on r3 (pure Neumann; the reference cell)
    grounded (solution 1):  phi_e = Phi + eta3g (c0 - Phi),
        eta3g = ((r - r2)/(r3 - r2))^2, so phi_e = c0 on r3

Both closures leave the value and the normal derivative on r2 unchanged.
The common mode leaves Vm, and so the ionic closure, unchanged.

--emit writes the LaTeX fragments latex/shell_bath_*.tex; --write writes the C++
reference header (namespace Foam::shellBathBidomainMMS) with mpmath
self-test samples.

    python bath_shell_bidomain_mms.py --emit
    python bath_shell_bidomain_mms.py --self-check
    python bath_shell_bidomain_mms.py --write ../../../src/verificationModels/bathBidomainVerification/manufacturedShellBathBidomainReference.H
    python bath_shell_bidomain_mms.py --check ../../../src/verificationModels/bathBidomainVerification/manufacturedShellBathBidomainReference.H
"""

import argparse
import copy
import functools
import math
import os
import sys

import mpmath as mp
import numpy as np
import sympy as sp

import mms_codegen as codegen
from mms_common.latex_emit import (emit_cases_fragment, emit_latex_fragment,
                                   emit_parameter_fragment, format_value)
import shell_mms as stm
from shell_bidomain_mms import _components, _full, _ionic_G, _tensor
from shell_monodomain_mms import (_blend, _div, _dot, _fibonacci_sphere, _grad, _hhat,
                                  _radial, _wall_correction, _walls)

HERE = os.path.dirname(os.path.abspath(__file__))
HEADER_NAME = "manufacturedShellBathBidomainReference.H"
HEADER_PATH = "../../../src/verificationModels/bathBidomainVerification/" + HEADER_NAME
SAMPLE_DPS = 30

# phiERefPoint of the insulated case: in the outer bath
PHIE_REF_POINT = (1.2, 0.3, -0.2)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

SCALAR_PARAMS = (
    "r1", "r2", "r3", "gamma1", "gamma2", "k",
    "gil", "git", "gin", "gel", "get", "gen", "theta0", "theta1",
    "sigma0In", "sigma0Out", "sigmaEps", "c0",
)
TAIL_PARAMS = ("cmAmplitude", "cmOmega")
ARRAY_PARAMS = (("fA", 4), ("hA", 3))
ANGLE_PARAMS = ("gamma1", "gamma2", "theta0", "theta1")

PARAM_ORDER = (
    list(SCALAR_PARAMS)
    + ["fA%d" % i for i in range(4)]
    + ["hA%d" % i for i in range(3)]
    + list(TAIL_PARAMS)
)
P = {name: sp.Symbol(name, real=True) for name in PARAM_ORDER}

# St = sqrt(1 + t), Tt = cmAmplitude sin(cmOmega t)
x, y, z, t = sp.symbols("x y z t", real=True)
St_, Tt_ = sp.symbols("St Tt", real=True)
ARGS = [x, y, z, t, St_, Tt_] + [P[n] for n in PARAM_ORDER]
XYZ = sp.Matrix([x, y, z])

ZONES = (0, 1, 2)
ZONE_NAMES = {0: "bloodPool", 1: "myocardium", 2: "bath"}
HEART_ZONE = 1
# outer closure: 0 insulated (O1), 1 grounded (O2)
SOLUTIONS = ("insulated", "grounded")


def default_params():
    """The shells' geometry, tensors and base fields, the curved conductor's
    bath conductivity, the outer sphere r3 = 1.5 and the ground value c0."""
    return {
        "r1": 0.5, "r2": 1.0, "r3": 1.5, "gamma1": 20.0, "gamma2": 50.0, "k": math.pi,
        "gil": 0.17, "git": 0.019, "gin": 0.01,
        "gel": 0.62, "get": 0.12, "gen": 0.08,
        "theta0": -20.0, "theta1": 90.0,
        "sigma0In": 0.6, "sigma0Out": 0.25, "sigmaEps": 0.2, "c0": 0.1,
        "fA": [0.3, 0.2, 0.4, 0.15],
        "hA": [-0.1, 0.3, 0.25],
        "cmAmplitude": 0.0, "cmOmega": 0.0,
    }


def alternate_params():
    """The common mode on, other base fields and a fibre rotation of 70 degrees."""
    p = default_params()
    p["cmAmplitude"] = 0.3
    p["cmOmega"] = 20.0
    p["fA"] = [0.25, 0.25, 0.35, 0.2]
    p["hA"] = [-0.15, 0.35, 0.15]
    p["theta1"] = 70.0
    p["c0"] = -0.05
    return p


def param_values(params):
    """Flat {symbol name: float}; angles in radians as Foam::degToRad does."""
    v = {}
    for name in SCALAR_PARAMS + TAIL_PARAMS:
        val = float(params[name])
        if name in ANGLE_PARAMS:
            val = val*math.pi/180.0
        v[name] = val
    for name, n in ARRAY_PARAMS:
        assert len(params[name]) == n, name
        for i in range(n):
            v["%s%d" % (name, i)] = float(params[name][i])
    return v


# ---------------------------------------------------------------------------
# Symbolic construction
# ---------------------------------------------------------------------------

def flux_correction(n_jump_grad_f, n_Ge_n, d):
    """g_j = n_j.(sigma_b,j I - Ge) grad f/((n_j.Ge n_j) d_j): the jump of the
    extracellular flux, Ge on the heart side, that omega g_j removes on sphere j."""
    return n_jump_grad_f/(n_Ge_n*d)


def mode_closure(n_Ge_grad_p, n_grad_f2, sigma):
    """r_j = n_j.Ge grad p/sigma_b - n_j.grad f2: the bath part of the common mode
    that carries its extracellular current across sphere j."""
    return n_Ge_grad_p/sigma - n_grad_f2


def pool_denominator_term(rho):
    """kappa (r1^2 - |x|^2)^2, kappa = gin/r1^4: zero on r1, gin at the origin."""
    r1 = P["r1"]
    return P["gin"]/r1**4*(r1**2 - rho)**2


def insulated_closure(Phi, r, r2, r3):
    """Phi + eta3 c3 with d/dr = 0 on r3 and value and d/dr unchanged on r2."""
    eta3 = (r - r2)**2*(r - r3)/(r3 - r2)**2
    c3 = -_dot(XYZ/r, _grad(Phi))
    return Phi + eta3*c3


def grounded_closure(Phi, r, r2, r3, c0):
    """Phi + eta3g (c0 - Phi), equal to c0 on r3, value and d/dr unchanged on r2."""
    eta3g = ((r - r2)/(r3 - r2))**2
    return Phi + eta3g*(c0 - Phi)


def _heart_exports(phiE, phiI, Gi, Ge):
    Vm = phiI - phiE
    G = _ionic_G()
    divI = _div(Gi*_grad(phiI))
    divE = _div(Ge*_grad(phiE))
    return {
        "phiE": phiE, "phiI": phiI, "Vm": Vm,
        "gradPhiE": list(_grad(phiE)),
        "gradPhiI": list(_grad(phiI)),
        "gradVm": list(_grad(Vm)),
        "divGiGradPhiI": divI,
        "divGeGradPhiE": divE,
        "phiESource": divI + divE,
        "J": list(Ge*_grad(phiE)),
        "G": G,
        "u1": (1 + t)*G + Vm,
        "u2": 1/(sp.sqrt(G)*(1 + t)),
    }


def _bath_exports(phiE, sigma):
    return {
        "phiE": phiE,
        "gradPhiE": list(_grad(phiE)),
        "phiESource": _div(sigma*_grad(phiE)),
        "sigmaB": sigma,
    }


@functools.lru_cache(maxsize=None)
def expressions(solution):
    """Per-zone expressions in (x, y, z, t, St, Tt, params)."""
    if solution not in SOLUTIONS:
        raise ValueError("solution must be one of %s" % (SOLUTIONS,))
    p = P
    r1, r2, r3 = p["r1"], p["r2"], p["r3"]
    th0, th1 = p["theta0"], p["theta1"]
    rho = x**2 + y**2 + z**2
    E = {0: {}, 1: {}, 2: {}}

    lam, om = _radial(rho, r1, r2)
    GiH, GeH = stm._tensors(th0 + th1*lam)
    # the heart tensors at the wall angle of each sphere, extended into its bath zone
    wallT = {1: stm._tensors(th0), 2: stm._tensors(th0 + th1)}
    sighat = stm._sigma_hat()
    sig = {0: p["sigma0In"]*(1 + p["sigmaEps"]*sighat),
           2: p["sigma0Out"]*(1 + p["sigmaEps"]*sighat)}
    sigj = {1: sig[0], 2: sig[2]}
    fhat, hhat = stm._fhat(), _hhat()
    f, h = St_*fhat, St_*hhat
    f2 = stm._f2()
    n, d = _walls(XYZ, r1, r2)

    # heart corrections, with the heart tensors
    gj, mj, qj = {}, {}, {}
    for j in (1, 2):
        nj = n[j]
        gj[j] = flux_correction(_dot(nj, (sigj[j]*sp.eye(3) - GeH)*_grad(f)),
                                _dot(nj, GeH*nj), d[j])
        mj[j] = _wall_correction(_dot(nj, GiH*_grad(h)), _dot(nj, GiH*nj), d[j])
        qj[j] = _wall_correction(_dot(nj, GiH*_grad(f2)), _dot(nj, GiH*nj), d[j])
    g = _blend(lam, gj[1], gj[2])
    m = _blend(lam, mj[1], mj[2])
    q = _blend(lam, qj[1], qj[2])
    pmode = f2 + om*q
    phiEH = f + om*g + Tt_*pmode
    phiIH = h + om*m + Tt_*pmode
    E[1].update(_heart_exports(phiEH, phiIH, GiH, GeH))
    E[1].update({"g": g, "m": m, "q": q, "GiZ": _components(GiH),
                 "GeZ": _components(GeH)})

    # bath zones: the common mode closed by beta_j r_j (beta_j = 0 and
    # n_j.grad beta_j = 1 on sphere j); grad p on sphere j is grad f2 + q_j d_j n_j,
    # with q_j at the wall angle.  Only the value of r_j on sphere j matters, so
    # its denominator n_j.Gi n_j, which the polynomial n_1 = -x/r1 makes vanish
    # at the origin, carries kappa (r1^2 - |x|^2)^2 in the pool: unchanged on
    # r1 and positive in the whole ball
    beta = {0: (r1**2 - rho)/(2*r1), 2: (rho - r2**2)/(2*r2)}
    poolDen = pool_denominator_term(rho)
    r = sp.sqrt(rho)
    for zn, j in ((0, 1), (2, 2)):
        Gi, Ge = wallT[j]
        nj = n[j]
        den = _dot(nj, Gi*nj) + (poolDen if zn == 0 else 0)
        qw = _wall_correction(_dot(nj, Gi*_grad(f2)), den, d[j])
        rj = mode_closure(_dot(nj, Ge*(_grad(f2) + qw*d[j]*nj)), _dot(nj, _grad(f2)),
                          sig[zn])
        Phi = f + Tt_*(f2 + beta[zn]*rj)
        if zn == 2:
            Phi = (insulated_closure(Phi, r, r2, r3) if solution == "insulated"
                   else grounded_closure(Phi, r, r2, r3, p["c0"]))
        E[zn].update(_bath_exports(Phi, sig[zn]))
        E[zn]["betaR"] = beta[zn]*rj
        E[zn]["GiZ"], E[zn]["GeZ"] = _components(Gi), _components(Ge)
    for zn in ZONES:
        E[zn].update({"fhat": fhat, "hhat": hhat, "f2": f2, "sighat": sighat})
    # the cell tensors are the heart formula everywhere, so a heart cell whose centre lies
    # beyond a polyhedral interface carries the tensor its sources use
    E["Gi"], E["Ge"] = E[1]["GiZ"], E[1]["GeZ"]
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
    return sp.lambdify(ARGS, _expr(solution, zone, name), modules=lib, cse=True)


def _context(pts, tt, pv, lib):
    """Values of ARGS at points (..., 3), as the C++ functions compute them."""
    if lib == "numpy":
        pts = np.asarray(pts, dtype=float)
        X, Y, Z = pts[..., 0], pts[..., 1], pts[..., 2]
        St, Tt = np.sqrt(1 + tt), pv["cmAmplitude"]*np.sin(pv["cmOmega"]*tt)
    else:
        X, Y, Z = pts
        St, Tt = mp.sqrt(1 + tt), pv["cmAmplitude"]*mp.sin(pv["cmOmega"]*tt)
    return [X, Y, Z, tt, St, Tt] + [pv[n] for n in PARAM_ORDER]


def _tdot(T6, v):
    M = [[T6[0], T6[1], T6[2]], [T6[1], T6[3], T6[4]], [T6[2], T6[4], T6[5]]]
    return [sum(M[i][kk]*v[kk] for kk in range(3)) for i in range(3)]


def _vdot(a, b):
    return sum(a[i]*b[i] for i in range(3))


class Model:
    def __init__(self, solution, params):
        self.solution = solution
        self.params = copy.deepcopy(params)
        self.pv = param_values(params)
        self.E = expressions(solution)

    def eval_xyz(self, zone, name, pts, tt):
        """numpy evaluation at points (..., 3); vectors/tensors stacked last."""
        pts = np.asarray(pts, dtype=float)
        shape = pts.shape[:-1]
        out = _compiled(self.solution, zone, name, "numpy")(*_context(pts, float(tt), self.pv, "numpy"))
        if isinstance(out, (list, tuple)):
            return np.stack([np.broadcast_to(np.asarray(c, dtype=float), shape) for c in out], -1)
        return np.broadcast_to(np.asarray(out, dtype=float), shape)

    def eval_mp(self, zone, name, pt, tt):
        """mpmath evaluation at a point (x, y, z)."""
        pv = {kk: mp.mpf(v) for kk, v in self.pv.items()}
        ptm = [mp.mpf(c) for c in pt]
        out = _compiled(self.solution, zone, name, "mpmath")(*_context(ptm, mp.mpf(tt), pv, "mpmath"))
        if isinstance(out, (list, tuple)):
            return [mp.mpf(c) for c in out]
        return mp.mpf(out)

    def zone_of(self, pt):
        rho = float(np.dot(pt, pt))
        if rho < self.pv["r1"]**2:
            return 0
        if rho > self.pv["r2"]**2:
            return 2
        return 1

    def radius(self, j):
        return self.pv["r%d" % j]

    def _point(self, j, u):
        um = [mp.mpf(float(c)) for c in u]
        nu = mp.sqrt(sum(c*c for c in um))
        um = [c/nu for c in um]
        return [mp.mpf(self.radius(j))*c for c in um], um

    # -- interface and boundary conditions ------------------------------------
    def interface_residuals(self, j, u, tt):
        """Residuals on sphere j (1, 2) at the unit direction u."""
        mp.mp.dps = SAMPLE_DPS
        zb = 0 if j == 1 else 2
        pt, um = self._point(j, u)
        n = [-c for c in um] if j == 1 else um
        ev = lambda zone, name: self.eval_mp(zone, name, pt, tt)  # noqa: E731
        Gi, Ge = ev(1, "GiZ"), ev(1, "GeZ")
        cont = ev(1, "phiE") - ev(zb, "phiE")
        flux = (_vdot(n, _tdot(Ge, ev(1, "gradPhiE")))
                - ev(zb, "sigmaB")*_vdot(n, ev(zb, "gradPhiE")))
        seal = _vdot(n, _tdot(Gi, ev(1, "gradPhiI")))
        jflux = _vdot(n, ev(1, "J")) - _vdot(n, _tdot(Ge, ev(1, "gradPhiE")))
        return {"continuity": float(cont), "extracellularFlux": float(flux),
                "intracellularSealing": float(seal), "interfaceFluxExport": float(jflux)}

    def outer_residual(self, u, tt):
        """d phi_e/dr (insulated) or phi_e - c0 (grounded) on r3."""
        mp.mp.dps = SAMPLE_DPS
        pt, um = self._point(3, u)
        if self.solution == "insulated":
            return float(_vdot(um, self.eval_mp(2, "gradPhiE", pt, tt)))
        return float(self.eval_mp(2, "phiE", pt, tt) - self.pv["c0"])

    def r2_unchanged_by_closure(self, u, tt):
        """phi_e and its gradient on r2 with and without the outer closure."""
        mp.mp.dps = SAMPLE_DPS
        pt, _ = self._point(2, u)
        other = Model("grounded" if self.solution == "insulated" else "insulated", self.params)
        a = [self.eval_mp(2, "phiE", pt, tt)] + self.eval_mp(2, "gradPhiE", pt, tt)
        b = [other.eval_mp(2, "phiE", pt, tt)] + other.eval_mp(2, "gradPhiE", pt, tt)
        return max(float(abs(u_ - v_)) for u_, v_ in zip(a, b))

    # -- sources versus finite differences ------------------------------------
    def sample_points(self, seed, n, min_distance_to_interface):
        rng = np.random.default_rng(seed)
        r1, r2, r3 = self.pv["r1"], self.pv["r2"], self.pv["r3"]
        pts = []
        while len(pts) < n:
            pt = rng.uniform(-r3, r3, 3)
            r = float(np.linalg.norm(pt))
            if r > r3 - min_distance_to_interface:
                continue
            if min(abs(r - r1), abs(r - r2)) < min_distance_to_interface:
                continue
            pts.append(pt)
        return pts

    def source_fd_error(self, pt, tt, h, zone=None):
        """Max relative error of the exact gradients/sources against 4th-order
        central differences of the fields; the stencil uses the formula of
        `zone` (default: the zone containing pt)."""
        zone = self.zone_of(pt) if zone is None else zone
        pt = np.asarray(pt, dtype=float)
        Dst = ((-2, 1.0/12.0), (-1, -8.0/12.0), (1, 8.0/12.0), (2, -1.0/12.0))
        outer = [(ax, mm, c) for ax in range(3) for mm, c in Dst]
        Q = np.array([pt + mm*h*np.eye(3)[ax] for ax, mm, c in outer])

        def grad_fd(name, X):
            gr = np.zeros((len(X), 3))
            for ax in range(3):
                for mm, c in Dst:
                    gr[:, ax] += c*self.eval_xyz(zone, name, X + mm*h*np.eye(3)[ax], tt)/h
            return gr

        def div_fd(Tname, name):
            if Tname == "sigma":
                T = self.eval_xyz(zone, "sigmaB", Q, tt)[:, None, None]*np.eye(3)
            else:
                T = _full(self.eval_xyz(zone, Tname, Q, tt))
            flux = np.einsum("nij,nj->ni", T, grad_fd(name, Q))
            return sum(c*flux[i, ax]/h for i, (ax, mm, c) in enumerate(outer))

        def rel(exact, approx):
            exact, approx = np.atleast_1d(exact), np.atleast_1d(approx)
            return float(np.max(np.abs(exact - approx)/np.maximum(1.0, np.abs(exact))))

        P0 = pt[None, :]
        ex = lambda name: self.eval_xyz(zone, name, P0, tt)[0]  # noqa: E731
        errs = [rel(ex("gradPhiE"), grad_fd("phiE", P0)[0])]
        if zone == 1:
            errs.append(rel(self.eval_xyz(zone, "Gi", P0, tt)[0], ex("GiZ")))
            errs.append(rel(self.eval_xyz(zone, "Ge", P0, tt)[0], ex("GeZ")))
            divI = div_fd("GiZ", "phiI")
            divE = div_fd("GeZ", "phiE")
            errs.append(rel(ex("divGiGradPhiI"), divI))
            errs.append(rel(ex("divGeGradPhiE"), divE))
            errs.append(rel(ex("phiESource"), divI + divE))
            errs.append(rel(ex("gradPhiI"), grad_fd("phiI", P0)[0]))
            errs.append(rel(ex("gradVm"), grad_fd("Vm", P0)[0]))
            errs.append(rel(ex("J"), _full(ex("GeZ")) @ grad_fd("phiE", P0)[0]))
        else:
            errs.append(rel(ex("phiESource"), div_fd("sigma", "phiE")))
        return max(errs)

    # -- boundedness, positivity, richness --------------------------------------
    def _names(self, zone):
        if zone == 1:
            return ["phiE", "phiI", "Vm", "gradPhiE", "gradPhiI", "gradVm",
                    "divGiGradPhiI", "divGeGradPhiE", "phiESource", "J", "G", "u1",
                    "u2", "GiZ", "GeZ", "g", "m", "q"]
        return ["phiE", "gradPhiE", "phiESource", "sigmaB", "GiZ", "GeZ", "betaR"]

    def zone_grid(self, zone, n=41, band=0.0):
        """Grid points of the ball r < r3 in `zone`, widened by `band` across its
        interface spheres (not across r3)."""
        r3 = self.pv["r3"]
        s = np.linspace(-r3, r3, n)
        G = np.stack(np.meshgrid(s, s, s, indexing="ij"), -1).reshape(-1, 3)
        r = np.linalg.norm(G, axis=1)
        r1, r2 = self.pv["r1"], self.pv["r2"]
        keep = {0: r < r1 + band, 1: (r > r1 - band) & (r < r2 + band),
                2: (r > r2 - band) & (r <= r3)}[zone]
        return G[keep]

    def dense_grid_values(self, tt=0.2, band=0.0):
        vals = []
        for zone in ZONES:
            G = self.zone_grid(zone, band=band)
            for name in self._names(zone):
                vals.append(np.ravel(self.eval_xyz(zone, name, G, tt)))
        return np.concatenate(vals)

    def min_denominators(self, band):
        """Smallest n_j.K n_j of the corrections (Gi for m and q, Ge for g) over the
        heart zone widened by band, and of the wall-angle tensors over each bath
        zone, and the smallest sigma_b."""
        out = {}
        G = self.zone_grid(1, band=band)
        for j in (1, 2):
            nj = (-G/self.pv["r1"]) if j == 1 else (G/self.pv["r2"])
            for name in ("GiZ", "GeZ"):
                T = _full(self.eval_xyz(1, name, G, 0.0))
                out["heart %s n%d" % (name, j)] = float(np.min(np.einsum("ni,nij,nj->n", nj, T, nj)))
        for zn, j in ((0, 1), (2, 2)):
            B = self.zone_grid(zn, band=band)
            B = B[np.linalg.norm(B, axis=1) > 1e-9]
            nj = (-B/self.pv["r1"]) if j == 1 else (B/self.pv["r2"])
            unit = nj/np.linalg.norm(nj, axis=1)[:, None]
            for name in ("GiZ", "GeZ"):
                T = _full(self.eval_xyz(zn, name, B, 0.0))
                out["zone %d %s" % (zn, name)] = float(np.min(np.einsum("ni,nij,nj->n", unit, T, unit)))
            out["zone %d sigmaB" % zn] = float(np.min(self.eval_xyz(zn, "sigmaB", B, 0.0)))
        # the pool denominator of q_1 over the whole ball r < r1 + band, origin included
        B = self.zone_grid(0, band=band)
        nj = -B/self.pv["r1"]
        Gi = _full(self.eval_xyz(0, "GiZ", B, 0.0))
        rho = np.einsum("ni,ni->n", B, B)
        kappa = self.pv["gin"]/self.pv["r1"]**4
        out["pool q denominator"] = float(np.min(np.einsum("ni,nij,nj->n", nj, Gi, nj)
                                                 + kappa*(self.pv["r1"]**2 - rho)**2))
        return out

    def richness(self, nq=800, tt=0.0):
        """Per sphere, over an equal-area quadrature: RMS ratios of
        - the tangential part of grad phi_e on the heart side to |grad phi_e|;
        - the tangential part of Ge n to n.Ge n (obliqueness of the flux);
        - the tangential part of grad Vm to |grad Vm|;
        - the jump of the normal derivative of phi_e to its heart-side value (kink);
        - n.Gi grad Vm, the virtual-electrode current, to |Gi grad Vm|."""
        out = {}
        for j in (1, 2):
            zb = 0 if j == 1 else 2
            pts = _fibonacci_sphere(nq, self.radius(j))
            n = (-1 if j == 1 else 1)*pts/self.radius(j)
            ev = lambda zone, name: self.eval_xyz(zone, name, pts, tt)  # noqa: E731
            Ge, Gi = _full(ev(1, "GeZ")), _full(ev(1, "GiZ"))
            gpH, gpB, gv = ev(1, "gradPhiE"), ev(zb, "gradPhiE"), ev(1, "gradVm")

            def tang(v):
                return v - np.einsum("ni,ni->n", v, n)[:, None]*n

            def rms_ratio(a, b):
                return float(np.sqrt(np.mean(a**2)/np.mean(b**2)))

            Gen = np.einsum("nij,nj->ni", Ge, n)
            dnH, dnB = np.einsum("ni,ni->n", gpH, n), np.einsum("ni,ni->n", gpB, n)
            Gigv = np.einsum("nij,nj->ni", Gi, gv)
            out[j] = {
                "tangential grad phiE": rms_ratio(np.linalg.norm(tang(gpH), axis=1),
                                                  np.linalg.norm(gpH, axis=1)),
                "obliqueness Ge n": rms_ratio(np.linalg.norm(tang(Gen), axis=1),
                                              np.einsum("ni,ni->n", Gen, n)),
                "tangential grad Vm": rms_ratio(np.linalg.norm(tang(gv), axis=1),
                                                np.linalg.norm(gv, axis=1)),
                "kink": rms_ratio(dnB - dnH, dnH),
                "virtual electrode": rms_ratio(np.einsum("ni,ni->n", Gigv, n),
                                               np.linalg.norm(Gigv, axis=1)),
            }
        return out

    # -- compatibility ------------------------------------------------------------
    def source_integral(self, tt=0.13, nr=60, nmu=64, nphi=128):
        """Gauss-Legendre integrals of the phi_e source and of its modulus over the
        three zones (r, cos(polar angle)) with the periodic trapezoid in the azimuth."""
        r1, r2, r3 = self.pv["r1"], self.pv["r2"], self.pv["r3"]
        xr, wr = np.polynomial.legendre.leggauss(nr)
        xm, wm = np.polynomial.legendre.leggauss(nmu)
        ph = 2*np.pi*np.arange(nphi)/nphi
        wph = 2*np.pi/nphi
        total, total_abs = 0.0, 0.0
        for zn, (a, b) in ((0, (0.0, r1)), (1, (r1, r2)), (2, (r2, r3))):
            rr = 0.5*(b - a)*xr + 0.5*(b + a)
            wrr = 0.5*(b - a)*wr*rr**2
            R, MU, PH = np.meshgrid(rr, xm, ph, indexing="ij")
            W = np.einsum("i,j->ij", wrr, wm)[:, :, None]*wph
            S = np.sqrt(1 - MU**2)
            pts = np.stack([R*S*np.cos(PH), R*S*np.sin(PH), R*MU], -1)
            s = self.eval_xyz(zn, "phiESource", pts, tt)
            total += float(np.sum(W*s))
            total_abs += float(np.sum(W*np.abs(s)))
        return total, total_abs

    # -- ionic closure ----------------------------------------------------------
    def ionic_residual(self):
        """ODE closure of the FDA ionic model (Cm = 2, Chi = 3) with the exact
        states, and the Vm-equation balance with the generated source, at
        chi*cm != Chi*Cm (as shell_mms.py)."""
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
        f = sp.lambdify(ARGS, res, modules="mpmath", cse=True)
        pv = {kk: mp.mpf(v) for kk, v in self.pv.items()}
        worst = mp.mpf(0)
        for pt in self.sample_points(seed=7, n=80, min_distance_to_interface=0.0):
            if self.zone_of(pt) != 1:
                continue
            ptm = [mp.mpf(c) for c in pt]
            for tt in (0.0, 0.07, 0.13, 0.5):
                for rr in f(*_context(ptm, mp.mpf(tt), pv, "mpmath")):
                    worst = max(worst, abs(rr))
        return float(worst)


def build(solution, params):
    return Model(solution, params)


def vm_common_mode_difference(solution, params_default, params_alternate):
    """Vm, grad Vm, u1, u2 with and without the common mode."""
    mp.mp.dps = SAMPLE_DPS
    on = copy.deepcopy(params_default)
    on["cmAmplitude"] = params_alternate["cmAmplitude"]
    on["cmOmega"] = params_alternate["cmOmega"]
    off = copy.deepcopy(params_alternate)
    off["cmAmplitude"] = 0.0
    worst = 0.0
    for pa, pb in ((params_default, on), (params_alternate, off)):
        ma, mb = build(solution, pa), build(solution, pb)
        for pt in ma.sample_points(seed=3, n=30, min_distance_to_interface=0.0):
            if ma.zone_of(pt) != 1:
                continue
            for tt in (0.0, 0.07, 0.137):
                for name in ("Vm", "gradVm", "u1", "u2"):
                    a, b = ma.eval_mp(1, name, pt, tt), mb.eval_mp(1, name, pt, tt)
                    a, b = (a, b) if isinstance(a, list) else ([a], [b])
                    worst = max(worst, max(float(abs(u - v)) for u, v in zip(a, b)))
    return worst


# ---------------------------------------------------------------------------
# Self-checks (run before any header is written; also used by the tests)
# ---------------------------------------------------------------------------

EXTENSION_BAND = 0.05
INTERFACE_TOL = 1e-12
COMPATIBILITY_TOL = 1e-12
# thresholds of the richness ratios (RMS over each sphere)
RICHNESS_MIN = {"tangential grad phiE": 1e-2, "obliqueness Ge n": 1e-2,
                "tangential grad Vm": 1e-2, "kink": 0.1, "virtual electrode": 0.1}


def _require(cond, *info):
    if not cond:
        raise AssertionError(" ".join(str(i) for i in info) or "check failed")


def _all_models():
    return [build(sol, p) for sol in SOLUTIONS for p in (default_params(), alternate_params())]


def check_interfaces():
    rng = np.random.default_rng(11)
    dirs = rng.normal(size=(24, 3))
    for m in _all_models():
        for j in (1, 2):
            for u in dirs:
                for tt in (0.0, 0.07, 0.2):
                    r = m.interface_residuals(j, u, tt)
                    _require(all(abs(v) < INTERFACE_TOL for v in r.values()),
                             "interface", m.solution, j, u, tt, r)


def check_outer_closure():
    rng = np.random.default_rng(13)
    dirs = rng.normal(size=(24, 3))
    for m in _all_models():
        for u in dirs:
            for tt in (0.0, 0.07, 0.2):
                _require(abs(m.outer_residual(u, tt)) < INTERFACE_TOL,
                         "outer closure", m.solution, u, tt)
                _require(m.r2_unchanged_by_closure(u, tt) < INTERFACE_TOL,
                         "closure changes r2", m.solution, u, tt)


def check_sources_match_fd():
    for m in _all_models():
        pts = m.sample_points(seed=1, n=40, min_distance_to_interface=2e-3)
        pts += [np.array([1e-3, -2e-3, 1.5e-3]), np.array([0.0, 0.0, 0.0])]
        for pt in pts:
            for tt in (0.0, 0.13):
                _require(m.source_fd_error(pt, tt, h=1e-3) < 1e-6, "source/gradient FD",
                         m.solution, pt, tt)
        # zone formulas extended across their sphere, as polyhedral faces need
        for j, (inner, outer) in ((1, (0, 1)), (2, (1, 2))):
            for u in np.random.default_rng(5).normal(size=(6, 3)):
                u = u/np.linalg.norm(u)
                for dr in (-0.5*EXTENSION_BAND, 0.5*EXTENSION_BAND):
                    pt = (m.radius(j) + dr)*u
                    for zone in (inner, outer):
                        _require(m.source_fd_error(pt, 0.13, h=1e-3, zone=zone) < 1e-6,
                                 "extension FD", m.solution, j, zone, pt)


def check_bounded_and_finite():
    for sol in SOLUTIONS:
        v = build(sol, alternate_params()).dense_grid_values(band=EXTENSION_BAND)
        _require(np.all(np.isfinite(v)) and np.max(np.abs(v)) < 1e3, "bounded/finite",
                 sol, np.max(np.abs(v)))


def check_denominators():
    for p in (default_params(), alternate_params()):
        d = build("insulated", p).min_denominators(band=EXTENSION_BAND)
        _require(min(d.values()) > 1e-3, "denominators", d)


def check_richness():
    for p in (default_params(), alternate_params()):
        rich = build("insulated", p).richness()
        for j, vals in rich.items():
            for name, v in vals.items():
                _require(v >= RICHNESS_MIN[name], "richness", j, name, v)


def check_compatibility():
    """Insulated outer surface: the phi_e source integrates to zero over the
    three zones; grounded: it does not have to (and does not)."""
    for p in (default_params(), alternate_params()):
        total, total_abs = build("insulated", p).source_integral()
        _require(abs(total) < COMPATIBILITY_TOL*total_abs, "compatibility", total, total_abs)
        total, total_abs = build("grounded", p).source_integral()
        _require(abs(total) > 1e-6*total_abs, "grounded source unexpectedly balanced",
                 total, total_abs)


def check_common_mode_keeps_vm():
    for sol in SOLUTIONS:
        dd = vm_common_mode_difference(sol, default_params(), alternate_params())
        _require(dd < 1e-14, "common mode changes Vm", sol, dd)


def check_ionic_closure():
    r = build("insulated", alternate_params()).ionic_residual()
    _require(r < 1e-13, "ionic closure", r)


def check_reference_point():
    m = build("insulated", default_params())
    r = float(np.linalg.norm(PHIE_REF_POINT))
    _require(m.pv["r2"] + EXTENSION_BAND < r < m.pv["r3"] - EXTENSION_BAND,
             "reference point outside the outer bath", r)


def check_shared_construction():
    """The heart-only bidomain shell is the sigma_b = 0, T = 0 member of this
    construction: with sigma_b = 0 the extracellular correction g_j is the
    sealing correction of shell_bidomain_mms.py, so phi_e = Fe and phi_i = Fi
    of that generator at t = 0 (S = 1)."""
    import shell_bidomain_mms as sbm
    p = default_params()
    p["sigma0In"] = p["sigma0Out"] = 0.0
    m = build("insulated", p)
    pts = np.array([(0.61, -0.29, 0.43), (-0.2, 0.55, 0.47), (0.33, 0.3, -0.62),
                    (0.5, 0.1, 0.05), (0.0, -0.8, 0.55)])
    for ours, theirs in (("phiE", "Fe"), ("phiI", "Fi")):
        a = m.eval_xyz(1, ours, pts, 0.0)
        b = sbm.evaluate(theirs, pts, sbm.default_params())
        _require(np.max(np.abs(a - b)) < 1e-13, "sigma_b = 0", ours, "differs from", theirs,
                 np.max(np.abs(a - b)))


SELF_CHECKS = (check_shared_construction, check_interfaces, check_outer_closure,
               check_sources_match_fd, check_bounded_and_finite, check_denominators,
               check_richness, check_compatibility, check_common_mode_keeps_vm,
               check_ionic_closure, check_reference_point)


def run_self_checks():
    """Every generator check except the header comparison; raises AssertionError."""
    for check in SELF_CHECKS:
        check()


# ---------------------------------------------------------------------------
# C++ emission
# ---------------------------------------------------------------------------

_THIS = sys.modules[__name__]

# fn index, C++ name, kind, expression name, value kind
FUNCTIONS = (
    (0, "Vm", "heart", "Vm", "scalar"),
    (1, "phiE", "zone", "phiE", "scalar"),
    (2, "phiI", "heart", "phiI", "scalar"),
    (3, "gradPhiE", "zone", "gradPhiE", "vector"),
    (4, "gradPhiI", "heart", "gradPhiI", "vector"),
    (5, "gradVm", "heart", "gradVm", "vector"),
    (6, "u1", "heart", "u1", "scalar"),
    (7, "u2", "heart", "u2", "scalar"),
    (8, "Gi", "tensor", "Gi", "symmTensor"),
    (9, "Ge", "tensor", "Ge", "symmTensor"),
    (10, "sigmaB", "bath", "sigmaB", "scalar"),
    (11, "divGiGradPhiI", "heart", "divGiGradPhiI", "scalar"),
    (12, "divGeGradPhiE", "heart", "divGeGradPhiE", "scalar"),
    (13, "phiESource", "zone", "phiESource", "scalar"),
    (14, "interfaceFlux", "flux", "J", "flux"),
    (15, "G", "heartNoTime", "G", "scalar"),
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


def printed_check_context():
    pv = param_values(alternate_params())
    ctx = dict(zip([str(a) for a in ARGS], _context(np.array([0.61, -0.29, 0.43]), 0.11, pv, "numpy")))
    return {kk: float(v) for kk, v in ctx.items()}


PREAMBLE = (
    (x, (), "const scalar x = pt.x();"),
    (y, (), "const scalar y = pt.y();"),
    (z, (), "const scalar z = pt.z();"),
    (St_, (), "const scalar St = sqrt(1 + t);"),
    (Tt_, ("cmAmplitude", "cmOmega"), "const scalar Tt = cmAmplitude*sin(cmOmega*t);"),
)

PREAMBLE_NEEDS = {}


def _sample_points():
    return {
        0: [(0.0, 0.0, 0.0), (0.07, 0.13, -0.05), (0.21, -0.3, 0.25), (-0.28, 0.24, -0.33)],
        1: [(0.36, 0.27, -0.23), (0.55, -0.47, 0.36), (-0.12, 0.64, 0.71)],
        2: [(0.62, 0.58, -0.59), (-1.1, 0.4, 0.83), (0.83, -0.81, 0.79)],
    }


SAMPLE_TIMES = (0.0, 0.137)
# functions without a time argument are sampled once, at the first time
TIME_INDEPENDENT_KINDS = ("heartNoTime", "bath", "tensor")


def self_test_samples():
    mp.mp.dps = SAMPLE_DPS
    sets = ((0, default_params()), (1, alternate_params()))
    rows = []
    for idx, cname, kind, ename, vkind in FUNCTIONS:
        for si, sol in enumerate(SOLUTIONS):
            for param, params in sets:
                m = build(sol, params)
                for zn in ZONES_OF_KIND[kind]:
                    for pt in _sample_points()[zn]:
                        assert m.zone_of(np.array(pt)) == zn, (zn, pt)
                        assert np.linalg.norm(pt) < m.pv["r3"], pt
                        for tt in (SAMPLE_TIMES if kind not in TIME_INDEPENDENT_KINDS
                                   else SAMPLE_TIMES[:1]):
                            val = m.eval_mp(zn, ename, pt, tt)
                            vals = val if isinstance(val, list) else [val]
                            for comp in range(N_COMPONENTS[vkind]):
                                rows.append((idx, param, si, zn, comp, pt[0], pt[1], pt[2],
                                             tt, float(vals[comp])))
    return rows


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
    Exact fields, conductivities, sources and self-test samples of the bath
    bidomain manufactured solution on the spherical shell: the heart shell
    r1 < |x| < r2 coupled every step to the bath ball bloodPool and the
    outer bath r2 < |x| < r3 (torsoCoupling bath), phi_i sealed and the
    extracellular current continuous through Ge on both spheres; the outer
    sphere insulated (solution 0) or grounded at c0 (solution 1).

    Regenerate from generators/:
        python bath_shell_bidomain_mms.py --write ../../../src/verificationModels/bathBidomainVerification/%s

    sympy %s
\\*---------------------------------------------------------------------------*/
"""


def generate_header():
    guard = "manufacturedShellBathBidomainReference_H"
    L = [LICENSE % (HEADER_NAME, sp.__version__)]
    L += ["#ifndef " + guard, "#define " + guard, ""]
    L += codegen.INCLUDES
    L += ["namespace Foam", "{", "namespace shellBathBidomainMMS", "{", ""]
    L += ["struct Params", "{", "    scalar r1, r2, r3, gamma1, gamma2, k;",
          "    scalar gil, git, gin, gel, get, gen, theta0, theta1;",
          "    scalar sigma0In, sigma0Out, sigmaEps, c0;",
          "    scalar fA[4], hA[3];", "    scalar cmAmplitude, cmOmega;",
          "    label solution;", "};", ""]
    L += ["enum Zone : label { bloodPool = 0, myocardium = 1, bath = 2 };", ""]
    L += ["//- Outer closure (Params::solution)",
          "enum Solution : label { insulated = 0, grounded = 1 };", ""]
    L += codegen.params_function(_THIS, "defaultParams", default_params(), 0) + [""]
    L += codegen.params_function(_THIS, "alternateParams", alternate_params(), 0) + [""]
    L += ["//- phiERefPoint of the insulated case (in the outer bath)",
          "inline point phiERefPoint()", "{",
          "    return point(%r, %r, %r);" % PHIE_REF_POINT, "}", ""]
    L += ["inline label zoneOf(const point& pt, const Params& p)", "{",
          "    const scalar rho = magSqr(pt);",
          "    if (rho < p.r1*p.r1)", "    {", "        return bloodPool;", "    }",
          "    if (rho > p.r2*p.r2)", "    {", "        return bath;", "    }",
          "    return myocardium;", "}", ""]
    L += codegen.check_zone_function("bloodPool", "bath", "myocardium")
    for spec in FUNCTIONS:
        L += codegen.emit_function(_THIS, *spec) + [""]
    L += codegen.self_test_table(_rows_cache())
    L += ["} // End namespace shellBathBidomainMMS", "} // End namespace Foam", "",
          "#endif", ""]
    return "\n".join(L)


def header_matches(path, text=None):
    return codegen.header_matches(__file__, generate_header, path, text)


# ---------------------------------------------------------------------------
# LaTeX fragments of the bath document
# ---------------------------------------------------------------------------

SHELL_TEX = dict(stm.SHELL_TEX)
SHELL_TEX.update({"r3": "r_3", "c0": "c_0"})
LATEX_PREFIX = "shell_bath_"


def _names(extra=()):
    names = {P[n]: tex for n, tex in SHELL_TEX.items()}
    names.update(dict(extra))
    return names


def latex_fragments():
    """The fragments the bath document adds to the shared shell and curved
    conductor fragments (latex/shell_*.tex, latex/shell_transmission_*.tex)."""
    out = {}
    lam, om, rho, r, S, T = sp.symbols("lam om rho r S T", real=True)
    base = [(lam, r"\lambda"), (om, r"\omega"), (rho, r"\lvert\mathbf x\rvert^2"),
            (r, "r"), (S, "S"), (T, "T")]
    a, b, dj = sp.symbols("a b dj", real=True)
    nf = _names(base + [(dj, "d_j"),
                        (a, r"\mathbf n_j\cdot\left(\sigma_{b,j}\mathbf I-\boldsymbol\sigma_e\right)\nabla f"),
                        (b, r"\left(\mathbf n_j\cdot\boldsymbol\sigma_e\mathbf n_j\right)")])
    out["gj"] = emit_latex_fragment("g_j", flux_correction(a, b, dj), nf)
    ns = _names(base + [(dj, "d_j"),
                        (a, r"\left(\mathbf n_j\cdot\boldsymbol\sigma_i\nabla h\right)"),
                        (b, r"\left(\mathbf n_j\cdot\boldsymbol\sigma_i\mathbf n_j\right)")])
    out["mj"] = emit_latex_fragment("m_j", _wall_correction(a, b, dj), ns)
    nq = _names(base + [(dj, "d_j"),
                        (a, r"\left(\mathbf n_j\cdot\boldsymbol\sigma_i\nabla f_2\right)"),
                        (b, r"\left(\mathbf n_j\cdot\boldsymbol\sigma_i\mathbf n_j\right)")])
    out["qj"] = emit_latex_fragment("q_j", _wall_correction(a, b, dj), nq)
    c1, c2, sj = sp.symbols("c1 c2 sj", real=True)
    nr = _names(base + [(c1, r"\left(\mathbf n_j\cdot\boldsymbol\sigma_e\nabla p\right)"),
                        (c2, r"\left(\mathbf n_j\cdot\nabla f_2\right)"),
                        (sj, r"\sigma_{b,j}")])
    out["rj"] = emit_latex_fragment("r_j", mode_closure(c1, c2, sj), nr)

    f, h, f2, g, m, q, beta_j, r_j, Phi = sp.symbols("f h f2 g m q beta_j r_j Phi", real=True)
    nz = _names(base + [(f, "f"), (h, "h"), (f2, "f_2"), (g, "g"), (m, "m"), (q, "q"),
                        (beta_j, r"\beta_j"), (r_j, "r_j"), (Phi, r"\Phi")])
    out["p"] = emit_latex_fragment("p", f2 + om*q, nz)
    pp = sp.Symbol("p", real=True)
    nz[pp] = "p"
    out["phiE"] = emit_cases_fragment(r"\phi_e", [
        (f + om*g + T*pp, r"r_1\le\lvert\mathbf x\rvert\le r_2"),
        (Phi, r"\text{bath zones, before the outer closure}")], nz)
    out["phiI"] = emit_latex_fragment(r"\phi_i", h + om*m + T*pp, nz)
    out["Phi"] = emit_latex_fragment(r"\Phi", f + T*(f2 + beta_j*r_j), nz)
    r2, r3, c0 = P["r2"], P["r3"], P["c0"]
    out["eta3"] = emit_latex_fragment(r"\eta_3", (r - r2)**2*(r - r3)/(r3 - r2)**2, nz)
    out["eta3g"] = emit_latex_fragment(r"\tilde\eta_3", ((r - r2)/(r3 - r2))**2, nz)

    p = default_params()
    pairs = []
    for name in SCALAR_PARAMS + TAIL_PARAMS:
        v = p[name]
        if name in ANGLE_PARAMS:
            v = r"%s^\circ" % format_value(float(v))
        pairs.append((SHELL_TEX[name], v))
    pairs += [(SHELL_TEX["hA%d" % i], v) for i, v in enumerate(p["hA"])]
    pairs += [(SHELL_TEX["fA%d" % i], v) for i, v in enumerate(p["fA"])]
    out["params"] = emit_parameter_fragment(pairs, per_row=5)
    return {"latex/%s%s.tex" % (LATEX_PREFIX, k_): v for k_, v in out.items()}


def render_outputs():
    """{path relative to this folder: text} of every emitted file."""
    return latex_fragments()


def write_outputs(root=HERE):
    outputs = render_outputs()
    for rel, text in outputs.items():
        path = os.path.join(root, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)
    return sorted(outputs)


def outputs_match(root=HERE):
    for rel, text in render_outputs().items():
        path = os.path.join(root, rel)
        if not os.path.exists(path):
            return False
        with open(path) as fh:
            if fh.read() != text:
                return False
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--write", metavar="PATH")
    grp.add_argument("--check", metavar="PATH")
    grp.add_argument("--emit", action="store_true")
    grp.add_argument("--self-check", action="store_true")
    args = ap.parse_args(argv)
    if args.self_check or args.emit:
        try:
            run_self_checks()
        except AssertionError as err:
            print("self-check failed%s: %s" % (", nothing emitted" if args.emit else "", err))
            return 1
        if args.emit:
            for rel in write_outputs():
                print("wrote %s" % rel)
        else:
            print("self-checks passed")
        return 0
    if args.write:
        try:
            run_self_checks()
        except AssertionError as err:
            print("self-check failed, header not written: %s" % err)
            return 1
        text = generate_header()
        path = codegen.resolve(__file__, args.write)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)
        print("wrote %s (%d bytes, %d self-test samples)"
              % (path, len(text.encode()), len(_rows_cache())))
        return 0
    if header_matches(args.check):
        print("%s is up to date" % codegen.resolve(__file__, args.check))
        return 0
    print("%s differs from the generator output" % codegen.resolve(__file__, args.check))
    return 1


@functools.lru_cache(maxsize=None)
def _rows_cache():
    return tuple(self_test_samples())


if __name__ == "__main__":
    sys.exit(main())
