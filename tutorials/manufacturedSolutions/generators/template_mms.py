#!/usr/bin/env python3
"""Template manufactured solution: passive monodomain on the unit square.

The smallest complete generator of this folder, to be copied when a new
manufactured solution is created.  Fields on [0, 1]^2, with t in units of 1 s
and x in units of 1 m:

    V_m = sqrt(1 + t) cos(pi x) cos(pi y)
    S   = chi C_m dV_m/dt - div(sigma grad V_m)

with the diagonal conductivity sigma = diag(sigma_x, sigma_y) and no ionic
current.  S is the applied current that makes V_m exact for

    chi C_m dV_m/dt = div(sigma grad V_m) + S,   n . sigma grad V_m = 0.

The x-derivative of cos(pi x) vanishes at x = 0 and x = 1 (likewise in y), so
the zero-flux condition holds by construction of the field.

Checks (--self-check): the residual of the equation vanishes on the exact
fields; the conormal flux vanishes on the four sides; S agrees with a
finite-difference evaluation of the operator applied to V_m; 1 + t, the
denominator of dV_m/dt, and sigma are positive; the printed constants
round-trip; the printed form of S equals the derived one.

This generator owns no C++ header, so it has no --write or --check.

    python template_mms.py --emit
    python template_mms.py --self-check

Each step to change in a copy is marked "REPLACE".
"""

import math
import os
import random
import sys

import sympy as sp

import mms_codegen as codegen
from mms_common.latex_emit import (emit_latex_fragment, emit_parameter_fragment,
                                   paper_symbol_names, printed_constant)
from mms_common.numpy_emit import emit_numpy_module

HERE = os.path.dirname(os.path.abspath(__file__))

# REPLACE: the name of reference/<name>.py and the prefix of latex/<name>_*.tex.
REFERENCE_NAME = "template"

# REPLACE: the coordinates of the domain (here the unit square) and time.
x, y, t = sp.symbols("x y t", real=True)
AXES = (x, y)
T_END = 1

# REPLACE: the constants of the problem instance, as exact rationals so that
# the residual vanishes exactly; floats would leave a rounding residual.
CHI = sp.Integer(3)
CM = sp.Integer(2)
SIGMA = sp.diag(sp.Rational("0.13"), sp.Rational("0.0176"))


def exact_field():
    """REPLACE: the exact field, chosen so that the boundary condition holds
    by construction (here d/dn cos(pi x) cos(pi y) = 0 on every side)."""
    return sp.sqrt(1 + t)*sp.cos(sp.pi*x)*sp.cos(sp.pi*y)


def operator(V, chi, cm, sigma):
    """REPLACE: the operator the solver implements, with its sign and scaling
    conventions: chi C_m dV/dt - div(sigma grad V), the left-hand side minus
    the diffusion term, so that the source is its value on the exact field."""
    grad = sp.Matrix([sp.diff(V, a) for a in AXES])
    flux = sigma*grad
    div = sum(sp.diff(flux[i], AXES[i]) for i in range(len(AXES)))
    return chi*cm*sp.diff(V, t) - div


def fields(chi=CHI, cm=CM, sigma=SIGMA):
    """{Vm, S} as sympy expressions in x, y, t; S is the residual of the
    operator on the exact field."""
    V = exact_field()
    S = sp.simplify(operator(V, chi, cm, sigma))
    return {"Vm": V, "S": S}


def printed_source(chi, cm, sx, sy):
    """REPLACE: S in the grouped form the document prints; a self-check
    asserts that it equals the derived S."""
    F = sp.cos(sp.pi*x)*sp.cos(sp.pi*y)
    return F*(chi*cm/(2*sp.sqrt(1 + t)) + sp.pi**2*(sx + sy)*sp.sqrt(1 + t))


# ---------------------------------------------------------------------------
# Self-checks.  REPLACE or extend them: every property the construction relies
# on gets a check that fails loudly before anything is written.
# ---------------------------------------------------------------------------

FD_STEP = 1e-4
FD_RTOL = 1e-6


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def residual():
    """chi C_m dV/dt - div(sigma grad V) - S on the exact fields."""
    f = fields()
    return sp.simplify(operator(f["Vm"], CHI, CM, SIGMA) - f["S"])


def conormal_fluxes():
    """{(axis, side): n . sigma grad V_m} on each side of the square."""
    V = exact_field()
    grad = sp.Matrix([sp.diff(V, a) for a in AXES])
    flux = SIGMA*grad
    return {(str(a), side): sp.simplify(flux[i].subs(a, side))
            for i, a in enumerate(AXES) for side in (0, 1)}


def max_source_fd_error(samples=20, seed=1):
    """Largest relative difference between S and the operator applied to V_m
    with second-order central differences, at random points."""
    f = fields()
    V = sp.lambdify((x, y, t), f["Vm"], "math")
    S = sp.lambdify((x, y, t), f["S"], "math")
    sx, sy = float(SIGMA[0, 0]), float(SIGMA[1, 1])
    h = FD_STEP
    rng = random.Random(seed)
    worst = 0.0
    for _ in range(samples):
        px, py, pt = rng.uniform(0, 1), rng.uniform(0, 1), rng.uniform(0, T_END)
        dt = (V(px, py, pt + h) - V(px, py, pt - h))/(2*h)
        dxx = (V(px + h, py, pt) - 2*V(px, py, pt) + V(px - h, py, pt))/h**2
        dyy = (V(px, py + h, pt) - 2*V(px, py, pt) + V(px, py - h, pt))/h**2
        terms = (float(CHI*CM)*dt, sx*dxx, sy*dyy)
        fd = math.fsum((terms[0], -terms[1], -terms[2]))
        # Scale by the size of the terms, not of their sum, which can cancel.
        scale = math.fsum(abs(v) for v in terms)
        worst = max(worst, abs(S(px, py, pt) - fd)/scale)
    return worst


def run_self_checks():
    _check(residual() == 0, "the residual does not vanish on the exact fields")
    for (axis, side), value in conormal_fluxes().items():
        _check(value == 0, "n . sigma grad V_m does not vanish on %s=%d" % (axis, side))
    printed = printed_source(CHI, CM, SIGMA[0, 0], SIGMA[1, 1])
    _check(sp.simplify(printed - fields()["S"]) == 0, "the printed S differs from S")
    error = max_source_fd_error()
    _check(error < FD_RTOL, "S differs from finite differences by %.3e" % error)
    # REPLACE: every denominator of the fields and sources must stay positive
    # on the domain and time interval (here 1 + t from d sqrt(1 + t)/dt).
    _check(sp.minimum(1 + t, t, sp.Interval(0, T_END)) > 0, "1 + t is not positive")
    _check(all(SIGMA[i, i] > 0 for i in range(2)), "sigma is not positive definite")
    for v in (SIGMA[0, 0], SIGMA[1, 1]):
        printed_constant(v, 3)


# ---------------------------------------------------------------------------
# Emission.  REPLACE the fragments and functions; keep the structure.
# ---------------------------------------------------------------------------

def latex_fragments():
    """{latex/<name>_*.tex: text}.  Symbols named as in PAPER_SYMBOLS print in
    the documents' notation; the rest are named here."""
    s = {n: sp.Symbol(n, real=True) for n in "chi cm sx sy".split()}
    names = paper_symbol_names(s.values())
    names.update({s["sx"]: r"\sigma_x", s["sy"]: r"\sigma_y"})
    out = {
        "template_Vm.tex": emit_latex_fragment("V_m", exact_field(), names),
        "template_S.tex": emit_latex_fragment(
            "S", printed_source(s["chi"], s["cm"], s["sx"], s["sy"]), names),
        "template_params.tex": emit_parameter_fragment(
            [(r"\chi", CHI), ("C_m", CM),
             (r"\sigma_x", printed_constant(SIGMA[0, 0], 3)),
             (r"\sigma_y", printed_constant(SIGMA[1, 1], 3))]),
    }
    return {"latex/" + k: v for k, v in out.items()}


def numpy_module():
    f = fields()
    functions = {name: ((x, y, t), f[name]) for name in ("Vm", "S")}
    doc = ("Template passive monodomain manufactured fields on the unit square; "
           "<name>(x, y, t) for Vm and S.  Generated by template_mms.py.")
    params = {"CHI": float(CHI), "CM": float(CM),
              "SIGMA": [float(SIGMA[0, 0]), float(SIGMA[1, 1])]}
    return emit_numpy_module(REFERENCE_NAME, functions, params, doc)


def render_outputs():
    """{path relative to this folder: text} of every emitted file."""
    out = {"reference/%s.py" % REFERENCE_NAME: numpy_module()}
    out.update(latex_fragments())
    return out


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


# REPLACE: a generator that owns a C++ header maps its file name to a function
# returning {region name: C++ lines}, as rotated_monodomain_mms.py does; the
# command line then gains --write PATH and --check PATH.
HEADERS = {}


def main(argv=None):
    return codegen.generator_main(argv, __doc__, __file__, run_self_checks, write_outputs,
                                  HEADERS)


if __name__ == "__main__":
    sys.exit(main())
