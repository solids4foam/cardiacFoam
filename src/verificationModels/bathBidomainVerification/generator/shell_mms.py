#!/usr/bin/env python3
"""Bath-bidomain transmission manufactured solution on the spherical shell.

A heart shell r1 < |x| < r2 lies between the bath ball bloodPool (|x| < r1)
and the outer bath (|x| > r2).  Every field is a polynomial, trigonometric or
exponential function of x, y, z, or a ratio whose denominator stays positive in
its zone and in a band beyond it, so each zone formula is analytic and extends
smoothly across the polyhedral interfaces of a mesh.  The exact fields satisfy
phiE continuity, extracellular flux continuity and intracellular sealing on
both spheres.  The module builds the fields and sources with sympy, runs the
self-checks used by test_shell_mms.py, and writes the C++ reference header
(namespace Foam::shellTransmissionMMS) with the function names and signatures
of the oblique-slab header.

    python shell_mms.py --write ../manufacturedShellTransmissionBathBidomainReference.H
    python shell_mms.py --check ../manufacturedShellTransmissionBathBidomainReference.H
"""

import argparse
import copy
import functools
import math
import sys

import mpmath as mp
import numpy as np
import sympy as sp

import mms_codegen as codegen

HEADER_NAME = "manufacturedShellTransmissionBathBidomainReference.H"
SAMPLE_DPS = 30

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

SCALAR_PARAMS = (
    "r1", "r2", "gamma1", "gamma2", "k",
    "gil", "git", "gin", "gel", "get", "gen", "theta0", "theta1",
    "sigma0In", "sigma0Out", "sigmaEps",
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
SOLUTIONS = ("levelSet",)
BOX_HALF_WIDTH = 1.5


def default_params():
    return {
        "r1": 0.5, "r2": 1.0, "gamma1": 20.0, "gamma2": 50.0, "k": math.pi,
        "gil": 0.17, "git": 0.019, "gin": 0.01,
        "gel": 0.62, "get": 0.12, "gen": 0.08,
        "theta0": -20.0, "theta1": 90.0,
        "sigma0In": 0.6, "sigma0Out": 0.25, "sigmaEps": 0.2,
        "fA": [0.3, 0.2, 0.4, 0.15],
        "hA": [-0.1, 0.3, 0.25],
        "cmAmplitude": 0.0, "cmOmega": 0.0,
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

def _grad(u):
    return sp.Matrix([sp.diff(u, v) for v in (x, y, z)])


def _div(J):
    return sum(sp.diff(J[i], v) for i, v in enumerate((x, y, z)))


def _dot(a, b):
    return (a.T*b)[0]


def _frame():
    """Fixed orthonormal frame (e_a, e_b, e_c): columns of Rz(gamma1) Rx(gamma2)."""
    c1, s1 = sp.cos(P["gamma1"]), sp.sin(P["gamma1"])
    c2, s2 = sp.cos(P["gamma2"]), sp.sin(P["gamma2"])
    Rz = sp.Matrix([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]])
    Rx = sp.Matrix([[1, 0, 0], [0, c2, -s2], [0, s2, c2]])
    Q = Rz*Rx
    return Q[:, 0], Q[:, 1], Q[:, 2]


def _tensor(l, tr, nn, theta):
    """l f f + tr s s + nn e_c e_c with f, s rotated by theta in the (e_a, e_b) plane."""
    ea, eb, ec = _frame()
    f = sp.cos(theta)*ea + sp.sin(theta)*eb
    s = -sp.sin(theta)*ea + sp.cos(theta)*eb
    return l*f*f.T + tr*s*s.T + nn*ec*ec.T


def _components(T):
    return [T[0, 0], T[0, 1], T[0, 2], T[1, 1], T[1, 2], T[2, 2]]


def _tensors(theta):
    i = (P["gil"], P["git"], P["gin"])
    e = (P["gel"], P["get"], P["gen"])
    return _tensor(*i, theta), _tensor(*e, theta)


def _ionic_G():
    k = P["k"]
    return 1 + sp.Rational(1, 4)*sp.sin(k*x + sp.Rational(1, 5))*sp.cos(k*y - sp.Rational(1, 3))*sp.cos(k*z)


def _heart_exports(phiE, phiI, Gi, Ge):
    Vm = phiI - phiE
    G = _ionic_G()
    divI = _div(Gi*_grad(phiI))
    return {
        "phiE": phiE, "phiI": phiI, "Vm": Vm,
        "gradPhiE": list(_grad(phiE)),
        "gradPhiI": list(_grad(phiI)),
        "gradVm": list(_grad(Vm)),
        "divGiGradPhiI": divI,
        "phiESource": divI + _div(Ge*_grad(phiE)),
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
    if solution != "levelSet":
        raise ValueError("solution must be 'levelSet'")
    p = P
    r1, r2, k = p["r1"], p["r2"], p["k"]
    th0, th1 = p["theta0"], p["theta1"]
    rho = x**2 + y**2 + z**2
    E = {0: {}, 1: {}, 2: {}}

    lam = (rho - r1**2)/(r2**2 - r1**2)
    om = (rho - r1**2)*(rho - r2**2)/((r2**2 - r1**2)*(r1 + r2))
    GiH, GeH = _tensors(th0 + th1*lam)
    GiL, GeL = _tensors(th0)
    GiR, GeR = _tensors(th0 + th1)
    sighat = (sp.sin(k*x + sp.Rational(1, 5))*sp.cos(k*y)
              * (1 + sp.Rational(1, 2)*sp.sin(k*z)))
    sig = {0: p["sigma0In"]*(1 + p["sigmaEps"]*sighat),
           2: p["sigma0Out"]*(1 + p["sigmaEps"]*sighat)}
    a = [p["fA%d" % i] for i in range(4)]
    b = [p["hA%d" % i] for i in range(3)]
    fhat = (a[0]*x + a[1]*y*z
            + a[2]*sp.sin(k*x)*(1 + sp.Rational(1, 2)*y)
            + a[3]*sp.cos(k*y + sp.Rational(3, 10))*sp.exp(sp.Rational(2, 5)*z))
    hhat = (b[0]*z + b[1]*sp.cos(k*y)*(1 - sp.Rational(3, 10)*x**2)
            + b[2]*sp.sin(k*z)*x)
    f, h = St_*fhat, St_*hhat
    f2 = sp.Rational(1, 5)*sp.cos(k*x) + sp.Rational(1, 10)*y + sp.Rational(1, 10)*z**2

    # unit normals on their own sphere, pointing from the heart into the bath zone
    n = {1: -XYZ/r1, 2: XYZ/r2}
    sigj = {1: sig[0], 2: sig[2]}
    d = {1: 2*r1/(r1 + r2), 2: 2*r2/(r1 + r2)}
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

    # n.Gi n on its sphere; on bloodPool kappa (r1^2 - rho), kappa below every eigenvalue of
    # Gi, keeps it >= kappa everywhere, the origin included
    kappa = 1/(1/p["gil"] + 1/p["git"] + 1/p["gin"])
    D = {1: (_dot(XYZ, GiL*XYZ) + kappa*(r1**2 - rho))/r1**2,
         2: _dot(XYZ, GiR*XYZ)/r2**2}

    def r_closed(nt, Gi, Ge, Dj, sigma):
        gf2 = _grad(f2)
        return (_dot(nt, Ge*gf2) - _dot(nt, Gi*gf2)*_dot(nt, Ge*nt)/Dj
                - sigma*_dot(nt, gf2))/sigma

    beta = {0: (r1**2 - rho)/(2*r1), 2: (rho - r2**2)/(2*r2)}
    for zn, j, Gi, Ge in ((0, 1, GiL, GeL), (2, 2, GiR, GeR)):
        r = r_closed(n[j], Gi, Ge, D[j], sig[zn])
        E[zn].update(_bath_exports(f + Tt_*(f2 + beta[zn]*r), sig[zn]))
        E[zn]["betaR"] = beta[zn]*r
        E[zn]["D"] = D[j]
    E[1]["GiZ"], E[1]["GeZ"] = _components(GiH), _components(GeH)
    E[0]["GiZ"], E[0]["GeZ"] = _components(GiL), _components(GeL)
    E[2]["GiZ"], E[2]["GeZ"] = _components(GiR), _components(GeR)
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
    e = _expr(solution, zone, name)
    return sp.lambdify(ARGS, e, modules=lib, cse=True)


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


def _full(T6):
    """symmTensor components (xx xy xz yy yz zz) as 3x3 matrices, stacked last."""
    return np.stack([np.stack([T6[..., 0], T6[..., 1], T6[..., 2]], -1),
                     np.stack([T6[..., 1], T6[..., 3], T6[..., 4]], -1),
                     np.stack([T6[..., 2], T6[..., 4], T6[..., 5]], -1)], -2)


def _fibonacci_sphere(n, radius):
    """n points of equal area on the sphere of the given radius."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2*i/n)
    th = math.pi*(1 + 5**0.5)*i
    return radius*np.stack([np.cos(th)*np.sin(phi), np.sin(th)*np.sin(phi), np.cos(phi)], -1)


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

    def eval_mp(self, zone, name, pt, tt, pv=None):
        """mpmath evaluation at a point (x, y, z)."""
        pv = {kk: mp.mpf(v) for kk, v in (pv or self.pv).items()}
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
        return self.pv["r1"] if j == 1 else self.pv["r2"]

    # -- interface conditions -----------------------------------------------
    def interface_residuals(self, j, u, tt):
        """Residuals at the point of sphere j in the unit direction u."""
        mp.mp.dps = SAMPLE_DPS
        zb = 0 if j == 1 else 2
        um = [mp.mpf(float(c)) for c in u]
        nu = mp.sqrt(sum(c*c for c in um))
        um = [c/nu for c in um]
        pt = [mp.mpf(self.radius(j))*c for c in um]
        n = [-c for c in um] if j == 1 else um

        def tdot(T6, v):
            M = [[T6[0], T6[1], T6[2]], [T6[1], T6[3], T6[4]], [T6[2], T6[4], T6[5]]]
            return [sum(M[i][kk]*v[kk] for kk in range(3)) for i in range(3)]

        def vdot(a, b):
            return sum(a[i]*b[i] for i in range(3))

        ev = lambda zone, name: self.eval_mp(zone, name, pt, tt)  # noqa: E731
        Ge, Gi = ev(1, "GeZ"), ev(1, "GiZ")
        cont = ev(1, "phiE") - ev(zb, "phiE")
        flux = vdot(n, tdot(Ge, ev(1, "gradPhiE"))) - ev(zb, "sigmaB")*vdot(n, ev(zb, "gradPhiE"))
        seal = vdot(n, tdot(Gi, ev(1, "gradPhiI")))
        return {"continuity": float(cont), "extracellularFlux": float(flux),
                "sealing": float(seal)}

    # -- sources versus finite differences ------------------------------------
    def sample_points(self, seed, n, min_distance_to_interface):
        rng = np.random.default_rng(seed)
        r1, r2 = self.pv["r1"], self.pv["r2"]
        pts = []
        while len(pts) < n:
            pt = rng.uniform(-BOX_HALF_WIDTH, BOX_HALF_WIDTH, 3)
            r = float(np.linalg.norm(pt))
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
            g = np.zeros((len(X), 3))
            for ax in range(3):
                for mm, c in Dst:
                    g[:, ax] += c*self.eval_xyz(zone, name, X + mm*h*np.eye(3)[ax], tt)/h
            return g

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
            errs.append(rel(ex("divGiGradPhiI"), divI))
            errs.append(rel(ex("phiESource"), divI + div_fd("GeZ", "phiE")))
            errs.append(rel(ex("gradVm"), grad_fd("Vm", P0)[0]))
            errs.append(rel(ex("J"), _full(ex("GeZ")) @ grad_fd("phiE", P0)[0]))
        else:
            errs.append(rel(ex("phiESource"), div_fd("sigma", "phiE")))
        return max(errs)

    # -- boundedness, positivity, richness --------------------------------------
    def _names(self, zone):
        if zone == 1:
            return ["phiE", "phiI", "Vm", "gradPhiE", "gradPhiI", "gradVm",
                    "divGiGradPhiI", "phiESource", "J", "G", "u1", "u2", "GiZ", "GeZ",
                    "g", "m", "q"]
        return ["phiE", "gradPhiE", "phiESource", "sigmaB", "GiZ", "GeZ", "betaR", "D"]

    def zone_grid(self, zone, n=41, band=0.0):
        """Grid points of the box in `zone`, widened by `band` across its spheres."""
        s = np.linspace(-BOX_HALF_WIDTH, BOX_HALF_WIDTH, n)
        G = np.stack(np.meshgrid(s, s, s, indexing="ij"), -1).reshape(-1, 3)
        r = np.linalg.norm(G, axis=1)
        r1, r2 = self.pv["r1"], self.pv["r2"]
        keep = {0: r < r1 + band, 1: (r > r1 - band) & (r < r2 + band), 2: r > r2 - band}[zone]
        return G[keep]

    def dense_grid_values(self, tt=0.2, band=0.0):
        vals = []
        for zone in ZONES:
            G = self.zone_grid(zone, band=band)
            for name in self._names(zone):
                vals.append(np.ravel(self.eval_xyz(zone, name, G, tt)))
        return np.concatenate(vals)

    def min_denominators(self, band):
        """Smallest D (bath zones) and n_j.G n_j (heart) over each zone widened by band."""
        out = {}
        for zone in (0, 2):
            out[ZONE_NAMES[zone]] = float(np.min(self.eval_xyz(zone, "D", self.zone_grid(zone, band=band), 0.0)))
        G = self.zone_grid(1, band=band)
        for j in (1, 2):
            nj = (-G/self.pv["r1"]) if j == 1 else (G/self.pv["r2"])
            for name in ("GiZ", "GeZ"):
                T = _full(self.eval_xyz(1, name, G, 0.0))
                out["heart %s n%d" % (name, j)] = float(np.min(np.einsum("ni,nij,nj->n", nj, T, nj)))
        return out

    def richness(self, nq=800, tt=0.0):
        """Per sphere (Gamma1, Gamma2) and per tensor (Ge, Gi), over an equal-area
        sphere quadrature: RMS of |(G n)_t . grad phi| over RMS of |n.Ge grad phiE|
        (Ge) or |Gi grad phiI| (Gi), and the area fraction where the pointwise
        ratio is >= 0.2."""
        rms, frac = [], []
        for j in (1, 2):
            pts = _fibonacci_sphere(nq, self.radius(j))
            n = (-1 if j == 1 else 1)*pts/self.radius(j)
            Ge = _full(self.eval_xyz(1, "GeZ", pts, tt))
            Gi = _full(self.eval_xyz(1, "GiZ", pts, tt))
            ge = self.eval_xyz(1, "gradPhiE", pts, tt)
            gi = self.eval_xyz(1, "gradPhiI", pts, tt)
            for G, gr, kind in ((Ge, ge, "e"), (Gi, gi, "i")):
                Gn = np.einsum("nij,nj->ni", G, n)
                Gn_t = Gn - np.einsum("ni,ni->n", Gn, n)[:, None]*n
                num = np.abs(np.einsum("ni,ni->n", Gn_t, gr))
                Gg = np.einsum("nij,nj->ni", G, gr)
                den = np.abs(np.einsum("ni,ni->n", Gg, n)) if kind == "e" else np.linalg.norm(Gg, axis=1)
                rms.append(float(np.sqrt(np.mean(num**2)/np.mean(den**2))))
                frac.append(float(np.mean(num >= 0.2*den)))
        return rms, frac

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
        f = sp.lambdify(ARGS, res, modules="mpmath", cse=True)
        pv = {kk: mp.mpf(v) for kk, v in self.pv.items()}
        worst = mp.mpf(0)
        for pt in self.sample_points(seed=7, n=80, min_distance_to_interface=0.0):
            if self.zone_of(pt) != 1:
                continue
            ptm = [mp.mpf(c) for c in pt]
            for tt in (0.0, 0.07, 0.13, 0.5):
                for r in f(*_context(ptm, mp.mpf(tt), pv, "mpmath")):
                    worst = max(worst, abs(r))
        return float(worst)


def build(solution, params):
    return Model(solution, params)


def vm_common_mode_difference(params_default, params_alternate):
    """Vm, grad Vm, u1, u2 with and without the common mode."""
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

EXTENSION_BAND = 0.05


def _require(cond, *info):
    if not cond:
        raise AssertionError(" ".join(str(i) for i in info) or "check failed")


def check_interfaces(m):
    rng = np.random.default_rng(11)
    for j in (1, 2):
        for u in rng.normal(size=(24, 3)):
            for tt in (0.0, 0.07, 0.2):
                r = m.interface_residuals(j, u, tt)
                _require(all(abs(v) < 1e-12 for v in r.values()), "interface", j, u, tt, r)


def check_all_interfaces():
    check_interfaces(build("levelSet", default_params()))
    check_interfaces(build("levelSet", alternate_params()))


def check_sources_match_fd():
    for p in (default_params(), alternate_params()):
        m = build("levelSet", p)
        pts = m.sample_points(seed=1, n=40, min_distance_to_interface=2e-3)
        pts += [np.array([1e-3, -2e-3, 1.5e-3]), np.array([0.0, 0.0, 0.0])]
        for pt in pts:
            for tt in (0.0, 0.13):
                _require(m.source_fd_error(pt, tt, h=1e-3) < 1e-6, "source/gradient FD", pt, tt)
        # zone formulas extended across their sphere, as polyhedral faces need
        for j, (inner, outer) in ((1, (0, 1)), (2, (1, 2))):
            for u in np.random.default_rng(5).normal(size=(6, 3)):
                u = u/np.linalg.norm(u)
                for dr in (-0.5*EXTENSION_BAND, 0.5*EXTENSION_BAND):
                    pt = (m.radius(j) + dr)*u
                    for zone in (inner, outer):
                        _require(m.source_fd_error(pt, 0.13, h=1e-3, zone=zone) < 1e-6,
                                 "extension FD", j, zone, pt)


def check_bounded_and_finite():
    v = build("levelSet", alternate_params()).dense_grid_values(band=EXTENSION_BAND)
    _require(np.all(np.isfinite(v)) and np.max(np.abs(v)) < 1e3, "bounded/finite", np.max(np.abs(v)))


def check_denominators():
    d = build("levelSet", default_params()).min_denominators(band=EXTENSION_BAND)
    _require(min(d.values()) > 1e-3, "denominators", d)


def check_richness():
    rms, frac = build("levelSet", default_params()).richness()
    _require(min(rms) >= 0.2 and min(frac) >= RICHNESS_FRACTION, "richness", rms, frac)


def check_common_mode_keeps_vm():
    dd = vm_common_mode_difference(default_params(), alternate_params())
    _require(dd < 1e-14, "common mode changes Vm", dd)


def check_ionic_closure():
    r = build("levelSet", default_params()).ionic_residual()
    _require(r < 1e-13, "ionic closure", r)


# (G n)_t vanishes where the sphere normal is an eigenvector of G, so a sphere cannot reach
# the slab's 0.8
RICHNESS_FRACTION = 0.7

SELF_CHECKS = (check_all_interfaces, check_sources_match_fd, check_bounded_and_finite,
               check_denominators, check_richness, check_common_mode_keeps_vm,
               check_ionic_closure)


def run_self_checks():
    """Every generator check except the header comparison; raises AssertionError."""
    for check in SELF_CHECKS:
        check()


# ---------------------------------------------------------------------------
# C++ emission
# ---------------------------------------------------------------------------

_THIS = sys.modules[__name__]

# fn index, C++ name, kind, expression name, value kind (as the slab header)
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
        2: [(0.62, 0.58, -0.59), (-1.1, 0.4, 0.83), (1.47, -1.41, 1.38)],
    }


SAMPLE_TIMES = (0.0, 0.137)
# functions without a time argument are sampled once, at the first time
TIME_INDEPENDENT_KINDS = ("heartNoTime", "bath", "tensor")


def self_test_samples():
    mp.mp.dps = SAMPLE_DPS
    sets = ((0, default_params()), (1, alternate_params()))
    rows = []
    for idx, cname, kind, ename, vkind in FUNCTIONS:
        for param, params in sets:
            m = build("levelSet", params)
            for zn in ZONES_OF_KIND[kind]:
                for pt in _sample_points()[zn]:
                    assert m.zone_of(np.array(pt)) == zn, (zn, pt)
                    for tt in (SAMPLE_TIMES if kind not in TIME_INDEPENDENT_KINDS
                               else SAMPLE_TIMES[:1]):
                        val = m.eval_mp(zn, ename, pt, tt)
                        vals = val if isinstance(val, list) else [val]
                        for comp in range(N_COMPONENTS[vkind]):
                            rows.append((idx, param, 1, zn, comp, pt[0], pt[1], pt[2], tt,
                                         float(vals[comp])))
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
    Exact fields, conductivities, sources and self-test samples of the
    bath-bidomain transmission manufactured solution on the spherical shell.

    Regenerate from generator/:
        python shell_mms.py --write ../%s

    sympy %s
\\*---------------------------------------------------------------------------*/
"""


def generate_header():
    guard = "manufacturedShellTransmissionBathBidomainReference_H"
    L = [LICENSE % (HEADER_NAME, sp.__version__)]
    L += ["#ifndef " + guard, "#define " + guard, ""]
    L += codegen.INCLUDES
    L += ["namespace Foam", "{", "namespace shellTransmissionMMS", "{", ""]
    L += ["struct Params", "{", "    scalar r1, r2, gamma1, gamma2, k;",
          "    scalar gil, git, gin, gel, get, gen, theta0, theta1;",
          "    scalar sigma0In, sigma0Out, sigmaEps;",
          "    scalar fA[4], hA[3];", "    scalar cmAmplitude, cmOmega;",
          "    label solution;", "};", ""]
    L += ["enum Zone : label { bloodPool = 0, myocardium = 1, bath = 2 };", ""]
    L += codegen.params_function(_THIS, "defaultParams", default_params(), 1) + [""]
    L += codegen.params_function(_THIS, "alternateParams", alternate_params(), 1) + [""]
    L += ["inline label zoneOf(const point& pt, const Params& p)", "{",
          "    const scalar rho = magSqr(pt);",
          "    if (rho < p.r1*p.r1)", "    {", "        return bloodPool;", "    }",
          "    if (rho > p.r2*p.r2)", "    {", "        return bath;", "    }",
          "    return myocardium;", "}", ""]
    L += codegen.check_zone_function("bloodPool", "bath", "myocardium")
    for spec in FUNCTIONS:
        L += codegen.emit_function(_THIS, *spec) + [""]
    L += codegen.self_test_table(_rows_cache())
    L += ["} // End namespace shellTransmissionMMS", "} // End namespace Foam", "",
          "#endif", ""]
    return "\n".join(L)


def header_matches(path, text=None):
    return codegen.header_matches(__file__, generate_header, path, text)


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
        path = codegen.resolve(__file__, args.write)
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
