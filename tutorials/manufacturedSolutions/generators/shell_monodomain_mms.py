#!/usr/bin/env python3
"""Heart-only monodomain manufactured solution on the spherical shell.

The heart is the shell r1 < |x| < r2; both spheres are sealed walls.  The
conductivity G has three distinct eigenvalues, with fibre and sheet rotating
through the wall in a fixed frame, so
it is full and oblique to both walls.  F = hhat + omega*m, with omega zero on
both spheres and m the blend of the per-sphere corrections, satisfies
n.G grad(F) = 0 exactly on each sphere.  The verifier uses V = sqrt(1 + t) F
and the FDA ionic closure with the positive function g.

--emit writes reference/shell_monodomain.py (F, gradF, divSigmaGradF, sigma
and g as NumPy functions of x, y, z and the parameters, with PARAMS and
field()) and the LaTeX building blocks latex/shell_*.tex and
latex/shell_monodomain_*.tex: lambda, omega, the fibre frame and sigma, the
wall normals and d_j, hhat, m_j, m, F, V_m and g.  div(sigma grad F) is not
printed; it is divSigmaGradF of the reference module.

    python shell_monodomain_mms.py --write ../../../src/verificationModels/monodomainVerification/manufacturedShellMonodomainReference.H
    python shell_monodomain_mms.py --check ../../../src/verificationModels/monodomainVerification/manufacturedShellMonodomainReference.H
    python shell_monodomain_mms.py --emit
    python shell_monodomain_mms.py --self-check
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
from mms_common.latex_emit import (emit_latex_fragment, emit_parameter_fragment,
                                   format_value, paper_symbol_names)
from mms_common.numpy_emit import emit_numpy_module

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "shell_monodomain"

HEADER_NAME = "manufacturedShellMonodomainReference.H"
HEADER_PATH = "../../../src/verificationModels/monodomainVerification/" + HEADER_NAME
SAMPLE_DPS = 30

SCALAR_PARAMS = ("r1", "r2", "gamma1", "gamma2", "k",
                 "gil", "git", "gin", "theta0", "theta1")
TAIL_PARAMS = ()
ARRAY_PARAMS = (("hA", 3),)
ANGLE_PARAMS = ("gamma1", "gamma2", "theta0", "theta1")
PARAM_ORDER = list(SCALAR_PARAMS) + ["hA%d" % i for i in range(3)]
P = {name: sp.Symbol(name, real=True) for name in PARAM_ORDER}

x, y, z = sp.symbols("x y z", real=True)
ARGS = [x, y, z] + [P[n] for n in PARAM_ORDER]
XYZ = sp.Matrix([x, y, z])

SOLUTIONS = ("shell",)
HEART_ZONE = 1
ZONE_NAMES = {1: "myocardium"}


def default_params():
    return {"r1": 0.5, "r2": 1.0, "gamma1": 20.0, "gamma2": 50.0, "k": math.pi,
            "gil": 0.17, "git": 0.019, "gin": 0.01, "theta0": -20.0, "theta1": 90.0,
            "hA": [-0.1, 0.3, 0.25]}


def alternate_params():
    p = default_params()
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
        for i in range(n):
            v["%s%d" % (name, i)] = float(params[name][i])
    return v


def _grad(u):
    return sp.Matrix([sp.diff(u, v) for v in (x, y, z)])


def _div(J):
    return sum(sp.diff(J[i], v) for i, v in enumerate((x, y, z)))


def _dot(a, b):
    return (a.T*b)[0]


def _conductivity(theta):
    """gil f f + git s s + gin e_c e_c, f and s rotated by theta in the (e_a, e_b) plane."""
    Rz, Rx = _rotations()
    Q = Rz*Rx
    ea, eb, ec = Q[:, 0], Q[:, 1], Q[:, 2]
    f, s = _fibre_sheet(theta, ea, eb)
    return _assemble(f*f.T, s*s.T, ec*ec.T)


def _rotations():
    """R_z(gamma1) and R_x(gamma2); the columns of their product are e_a, e_b, e_c."""
    c1, s1 = sp.cos(P["gamma1"]), sp.sin(P["gamma1"])
    c2, s2 = sp.cos(P["gamma2"]), sp.sin(P["gamma2"])
    return (sp.Matrix([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]]),
            sp.Matrix([[1, 0, 0], [0, c2, -s2], [0, s2, c2]]))


def _fibre_sheet(theta, ea, eb):
    """Fibre f and sheet s, rotated by theta in the (e_a, e_b) plane."""
    return sp.cos(theta)*ea + sp.sin(theta)*eb, -sp.sin(theta)*ea + sp.cos(theta)*eb


def _assemble(ff, ss, nn):
    """gil ff + git ss + gin nn from the outer products of fibre, sheet and e_c."""
    return P["gil"]*ff + P["git"]*ss + P["gin"]*nn


def _radial(rho, r1, r2):
    """lambda (0 on the inner, 1 on the outer sphere) and omega (zero on both) of rho = |x|^2."""
    lam = (rho - r1**2)/(r2**2 - r1**2)
    om = (rho - r1**2)*(rho - r2**2)/((r2**2 - r1**2)*(r1 + r2))
    return lam, om


def _walls(X, r1, r2):
    """Unit normals pointing out of the heart and d_j = grad(omega).n_j on sphere j."""
    return {1: -X/r1, 2: X/r2}, {1: 2*r1/(r1 + r2), 2: 2*r2/(r1 + r2)}


def _wall_correction(n_K_grad_h, n_K_n, d):
    """m_j = -(n_j.K grad hhat)/((n_j.K n_j) d_j) for the operator tensor K."""
    return -n_K_grad_h/(n_K_n*d)


def _blend(lam, m1, m2):
    return (1 - lam)*m1 + lam*m2


def _hhat():
    b = [P["hA%d" % i] for i in range(3)]
    return (b[0]*z + b[1]*sp.cos(P["k"]*y)*(1 - sp.Rational(3, 10)*x**2)
            + b[2]*sp.sin(P["k"]*z)*x)


def _g():
    k = P["k"]
    return 1 + sp.Rational(1, 4)*sp.sin(k*x + sp.Rational(1, 5))*sp.cos(k*y - sp.Rational(1, 3))*sp.cos(k*z)


@functools.lru_cache(maxsize=None)
def expressions():
    r1, r2 = P["r1"], P["r2"]
    lam, om = _radial(x**2 + y**2 + z**2, r1, r2)
    G = _conductivity(P["theta0"] + P["theta1"]*lam)
    hhat = _hhat()
    n, d = _walls(XYZ, r1, r2)
    mj = {j: _wall_correction(_dot(n[j], G*_grad(hhat)), _dot(n[j], G*n[j]), d[j]) for j in (1, 2)}
    m = _blend(lam, mj[1], mj[2])
    F = hhat + om*m
    gradF = _grad(F)
    g = _g()
    return {"F": F, "gradF": list(gradF), "divGGradF": _div(G*gradF),
            "G": [G[0, 0], G[0, 1], G[0, 2], G[1, 1], G[1, 2], G[2, 2]], "g": g,
            "n1": n[1], "n2": n[2], "Gmat": G, "denominators": [_dot(n[j], G*n[j]) for j in (1, 2)]}


def _expr(solution, zone, name):
    return expressions()[name]


FUNCTIONS = (
    (0, "F", "heart", "F", "scalar"),
    (1, "gradF", "heart", "gradF", "vector"),
    (2, "divGGradF", "heart", "divGGradF", "scalar"),
    (3, "G", "heart", "G", "symmTensor"),
    (4, "g", "heart", "g", "scalar"),
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


def max_sealing_residual(n=200, params=None):
    """max |n.G grad(F)| on both spheres, mpmath at SAMPLE_DPS digits."""
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
            gF = evaluate("gradF", ptm, p, "mpmath")
            Gc = evaluate("G", ptm, p, "mpmath")
            Gm = [[Gc[0], Gc[1], Gc[2]], [Gc[1], Gc[3], Gc[4]], [Gc[2], Gc[4], Gc[5]]]
            res = sum(nrm[i]*Gm[i][j]*gF[j] for i in range(3) for j in range(3))
            worst = max(worst, abs(res))
    return float(worst)


def max_source_fd_error(h=1e-5, n=40):
    """divGGradF against central differences of the flux G grad(F).

    The error is the O(h^2) truncation of the difference (it falls by 100 per
    decade of h), set by the large third derivatives of the flux (m carries
    1/(n.G n) with gin small), so h is chosen small enough for 1e-6.
    """
    pts = _shell_points(n, 0.0, seed=2)
    err = 0.0
    for pt in pts:
        fd = 0.0
        for a in range(3):
            e = np.zeros(3)
            e[a] = h
            Gp = evaluate("G", pt + e)
            Gm = evaluate("G", pt - e)

            def flux(G6, q):
                Gm3 = np.array([[G6[0], G6[1], G6[2]], [G6[1], G6[3], G6[4]], [G6[2], G6[4], G6[5]]])
                return Gm3 @ evaluate("gradF", q)
            fd += (flux(Gp, pt + e)[a] - flux(Gm, pt - e)[a])/(2*h)
        err = max(err, abs(fd - float(evaluate("divGGradF", pt))))
    return err


def min_denominator(band):
    pts = _shell_points(4000, band, seed=3)
    d = np.stack([evaluate("denominators", pts)[..., j] for j in (0, 1)], -1)
    return float(d.min())


def min_g(band):
    return float(evaluate("g", _shell_points(4000, band, seed=4)).min())


def wall_richness(n=400):
    """Over the two walls, the smaller wall mean of |tangential part of G n|/gil and of
    |tangential grad F|.

    The tangential part of G n is a tangent vector field on a sphere, so it
    vanishes at isolated points (hairy-ball theorem) and its pointwise minimum
    is zero in the continuum; the wall mean measures how oblique G is.
    """
    p = default_params()
    tc, tg = np.inf, np.inf
    for radius, sign in ((p["r1"], -1), (p["r2"], 1)):
        pts = _fibonacci_sphere(n, radius)
        nrm = sign*pts/radius
        G6 = evaluate("G", pts)
        Gm = np.stack([np.stack([G6[:, 0], G6[:, 1], G6[:, 2]], -1),
                       np.stack([G6[:, 1], G6[:, 3], G6[:, 4]], -1),
                       np.stack([G6[:, 2], G6[:, 4], G6[:, 5]], -1)], -2)
        Gn = np.einsum("pij,pj->pi", Gm, nrm)
        gF = evaluate("gradF", pts)
        tang = lambda v: v - np.sum(v*nrm, -1)[:, None]*nrm
        tc = min(tc, np.linalg.norm(tang(Gn), axis=1).mean()/p["gil"])
        tg = min(tg, np.linalg.norm(tang(gF), axis=1).mean())
    return float(tc), float(tg)


def run_self_checks():
    assert max_sealing_residual() < 1e-25, "sealing"
    assert max_source_fd_error() < 1e-6, "source"
    assert min_denominator(0.05) > 0, "denominator"
    assert min_g(0.05) > 0.5, "g"
    tc, tg = wall_richness()
    assert tc > 1e-2 and tg > 1e-2, "richness"


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


RUNTIME = '''

PARAM_ORDER = %r


def field(name, x, y, z, params=None, **overrides):
    """The function name at (x, y, z) with PARAMS["default"] (or params), updated by
    overrides; angles in radians."""
    values = dict(PARAMS["default"] if params is None else params)
    unknown = set(overrides) - set(PARAM_ORDER)
    if unknown:
        raise KeyError("unknown parameters %%s" %% sorted(unknown))
    values.update(overrides)
    return globals()[name](x, y, z, *[values[n] for n in PARAM_ORDER])
'''


def reference_module(module, functions, doc):
    """Emitted NumPy module of a shell generator: functions of (x, y, z, *PARAM_ORDER),
    PARAMS with the default and alternate parameters, and field()."""
    args = tuple(module.ARGS)
    emitted = {name: (args, expr) for name, expr in functions.items()}
    params = {"default": module.param_values(module.default_params()),
              "alternate": module.param_values(module.alternate_params())}
    text = emit_numpy_module(module.REFERENCE_NAME, emitted, params, doc)
    return text.rstrip("\n") + "\n" + RUNTIME % (list(module.PARAM_ORDER),)


# LaTeX names of the shell symbols, in Paper I notation where it has one.
SHELL_SYMBOLS = {"r1": "r_1", "r2": "r_2", "gamma1": r"\gamma_1", "gamma2": r"\gamma_2",
                 "gil": r"\sigma_f", "git": r"\sigma_s", "gin": r"\sigma_n",
                 "theta0": r"\theta_0", "theta1": r"\theta_1", "k": "k",
                 "hA0": "a_1", "hA1": "a_2", "hA2": "a_3", "delta": r"\delta",
                 "xs0": r"x_{s,1}", "xs1": r"x_{s,2}", "xs2": r"x_{s,3}"}


def shell_symbol_names(extra=()):
    """symbol_names for the shell parameters, the paper symbols and the given
    (symbol, LaTeX) pairs."""
    names = {P[n] if n in P else sp.Symbol(n, real=True): tex for n, tex in SHELL_SYMBOLS.items()}
    names.update(paper_symbol_names([sp.Symbol(n, real=True) for n in ("chi", "cm", "c0")]))
    names.update(dict(extra))
    return names


def geometry_fragments():
    """lambda, omega, theta, the fibre frame, f, s, sigma, the wall normals and d_j."""
    rho, lam, theta, xv, ea, eb = sp.symbols("rho lam theta xv ea eb", real=True)
    ff, ss, nn = sp.symbols("ff ss nn", real=True)
    names = shell_symbol_names([(rho, r"\lvert\mathbf x\rvert^2"), (lam, r"\lambda"),
                                (theta, r"\theta"), (xv, r"\mathbf x"), (ea, r"\mathbf e_a"),
                                (eb, r"\mathbf e_b"), (ff, r"\mathbf f\otimes\mathbf f"),
                                (ss, r"\mathbf s\otimes\mathbf s"),
                                (nn, r"\mathbf n_s\otimes\mathbf n_s")])
    r1, r2 = P["r1"], P["r2"]
    lam_e, om_e = _radial(rho, r1, r2)
    Rz, Rx = _rotations()
    f, s = _fibre_sheet(theta, ea, eb)
    n, d = _walls(xv, r1, r2)
    out = {
        "shell_lambda.tex": emit_latex_fragment(r"\lambda", lam_e, names),
        "shell_omega.tex": emit_latex_fragment(r"\omega", om_e, names),
        "shell_theta.tex": emit_latex_fragment(r"\theta", P["theta0"] + P["theta1"]*lam, names),
        "shell_frame.tex": emit_latex_fragment(
            r"\left[\mathbf e_a\;\mathbf e_b\;\mathbf n_s\right]",
            sp.MatMul(Rz, Rx, evaluate=False), names),
        "shell_fibre.tex": emit_latex_fragment(r"\mathbf f", f, names),
        "shell_sheet.tex": emit_latex_fragment(r"\mathbf s", s, names),
        "shell_sigma.tex": emit_latex_fragment(r"\boldsymbol\sigma", _assemble(ff, ss, nn), names),
    }
    for j in (1, 2):
        out["shell_n%d.tex" % j] = emit_latex_fragment(r"\mathbf n_%d" % j, n[j], names)
        out["shell_d%d.tex" % j] = emit_latex_fragment(r"d_%d" % j, d[j], names)
    return out


def blend_fragments(prefix, tensor_tex, field_tex):
    """m_j for the operator tensor tensor_tex, the blend m, and the field hhat + omega m."""
    a, b, dj, lam, m1, m2, hh, om, m = sp.symbols("a b dj lam m1 m2 hh om m", real=True)
    names = shell_symbol_names([
        (a, r"\left(\mathbf n_j\cdot %s\nabla\hat h\right)" % tensor_tex),
        (b, r"\left(\mathbf n_j\cdot %s\,\mathbf n_j\right)" % tensor_tex),
        (dj, "d_j"), (lam, r"\lambda"), (m1, "m_1"), (m2, "m_2"), (hh, r"\hat h"),
        (om, r"\omega"), (m, "m")])
    return {
        prefix + "_mj.tex": emit_latex_fragment("m_j", _wall_correction(a, b, dj), names),
        prefix + "_m.tex": emit_latex_fragment("m", _blend(lam, m1, m2), names),
        prefix + "_field.tex": emit_latex_fragment(field_tex, hh + om*m, names),
    }


def latex_fragments():
    t, Fs = sp.symbols("t F", real=True)
    names = shell_symbol_names([(Fs, "F")])
    out = geometry_fragments()
    out.update(blend_fragments("shell_monodomain", r"\boldsymbol\sigma", "F"))
    out.update({
        "shell_monodomain_hhat.tex": emit_latex_fragment(r"\hat h", _hhat(), names),
        "shell_monodomain_Vm.tex": emit_latex_fragment("V_m", sp.sqrt(1 + t)*Fs, names),
        "shell_monodomain_g.tex": emit_latex_fragment("g", _g(), names),
    })
    Vm, gs, beta, dsig = sp.symbols("Vm gs beta dsig", real=True)
    snames = shell_symbol_names([(Fs, "F"), (Vm, "V_m"), (gs, "g"), (beta, r"\beta"),
                                 (dsig, r"\nabla\!\cdot\!\left(\boldsymbol\sigma\nabla F\right)")])
    out.update({
        "shell_monodomain_u1.tex": emit_latex_fragment("u_1", (1 + t)*gs + Vm, snames),
        "shell_monodomain_u2.tex": emit_latex_fragment(
            "u_2", gs**sp.Rational(-1, 2)/(1 + t), snames),
        "shell_monodomain_source.tex": emit_latex_fragment(
            "S", sp.sqrt(1 + t)*(beta*Fs - dsig), snames),
        "shell_monodomain_params.tex": parameter_fragment(default_params()),
    })
    return {"latex/" + k: v for k, v in out.items()}


def parameter_fragment(params, extra=()):
    """The shell parameters in Paper I notation (angles in degrees)."""
    pairs = []
    for name in SCALAR_PARAMS:
        if name not in SHELL_SYMBOLS or name not in params:
            continue
        value = params[name]
        if name in ANGLE_PARAMS:
            value = r"%s^\circ" % format_value(float(value))
        pairs.append((SHELL_SYMBOLS[name], value))
    if "hA" in params:
        pairs += [(SHELL_SYMBOLS["hA%d" % i], params["hA"][i]) for i in range(3)]
    return emit_parameter_fragment(pairs + list(extra), per_row=5)


def numpy_module():
    e = expressions()
    functions = {"F": e["F"], "gradF": sp.Matrix(e["gradF"]), "divSigmaGradF": e["divGGradF"],
                 "sigma": e["Gmat"], "g": e["g"]}
    doc = ("Heart-only monodomain manufactured solution on the spherical shell r1 < |x| < r2: "
           "F, gradF (shape (3, ...)), divSigmaGradF = div(sigma grad F), sigma (shape "
           "(3, 3, ...)) and the ionic function g, each a function of (x, y, z, *PARAM_ORDER) "
           "with angles in radians; V_m = sqrt(1 + t) F.  PARAMS holds the default and "
           "alternate parameters; field(name, x, y, z, **overrides) evaluates with them.  "
           "Generated by shell_monodomain_mms.py.")
    return reference_module(_THIS, functions, doc)


def render_outputs():
    """{path relative to this folder: text} of every emitted file."""
    out = {"reference/%s.py" % REFERENCE_NAME: numpy_module()}
    out.update(latex_fragments())
    return out


def write_outputs(render=render_outputs, root=HERE):
    outputs = render()
    for rel, text in outputs.items():
        path = os.path.join(root, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)
    return sorted(outputs)


def outputs_match(render=render_outputs, root=HERE):
    for rel, text in render().items():
        path = os.path.join(root, rel)
        if not os.path.exists(path):
            return False
        with open(path) as fh:
            if fh.read() != text:
                return False
    return True


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
    Exact field F, its gradient, div(G grad F), the conductivity G and the
    positive ionic function g of the heart-only monodomain manufactured
    solution on the spherical shell, sealed on both spheres, with self-test
    samples.

    Regenerate from generator/:
        python shell_monodomain_mms.py --write %s

    sympy %s
\\*---------------------------------------------------------------------------*/
"""


def generate_header():
    guard = "manufacturedShellMonodomainReference_H"
    L = [LICENSE % (HEADER_PATH, sp.__version__)]
    L += ["#ifndef " + guard, "#define " + guard, ""]
    L += codegen.INCLUDES
    L += ["namespace Foam", "{", "namespace shellMonodomainMMS", "{", ""]
    L += ["struct Params", "{", "    scalar r1, r2, gamma1, gamma2, k;",
          "    scalar gil, git, gin, theta0, theta1;", "    scalar hA[3];",
          "    label solution;", "};", ""]
    L += codegen.params_function(_THIS, "defaultParams", default_params(), 0) + [""]
    L += codegen.params_function(_THIS, "alternateParams", alternate_params(), 0) + [""]
    for spec in FUNCTIONS:
        L += codegen.emit_function(_THIS, *spec) + [""]
    L += codegen.self_test_table(_rows_cache())
    L += ["} // End namespace shellMonodomainMMS", "} // End namespace Foam", "",
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
