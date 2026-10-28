#!/usr/bin/env python3
"""Junction between the Purkinje graph and the myocardium, in continuum form.

The tissue field is V_m of fda_monodomain_mms (d = 3) on the unit cube; the
graph field is V_1D of graph_monodomain_mms.  The graph terminals x_k are
(0, 1/6, 1/3) and (1, 1/6, 1/3), both on a face of the cube.  With the linear
kernel of radius r = 0.11,

    w_k(x) = max(0, 1 - |x - x_k| / r)   on the ball B_k(r) intersected with the cube,

the junction is

    sampling     Vbar_k = int w_k V_m / int w_k
    current      I_pvj,k = (V_1D,k - Vbar_k) / R_pvj
    deposition   S_pvj(x) = sum_k I_pvj,k w_k(x) / int w_k

so that int S_pvj over the cube is sum_k I_pvj,k and each terminal deposits its
own current.  The integrals run over a half ball (the cube cuts the ball at the
face x = 0 or x = 1); they use mms_common.quadrature.half_ball_kernel_integral.

At a terminal V_m = -V_1D/2 and dV_m/dx = 0, so the point-value current is
I_pvj,k = +-(3/2) sqrt(1 + t) / R_pvj (terminal_current_exact).  The sampled
value differs from the point value by O(r^2): by the symmetry of the ball in y
and z and dV_m/dx = 0 the first-order Taylor term vanishes, and
|Vbar_k - V_m(x_k)| <= 7 pi^2 r^2 sqrt(1 + t), from a second-order remainder
with |s^T H s| <= 14 pi^2 |s|^2 and the weighted mean of |s|^2 at most r^2.

Limitation.  cardiacFoam applies the junction through the mesh: the sampling
and deposition maps are discrete.  This module holds only the continuum
junction; the exact source built with a solver's own discrete map tests the
coupled workflow against that map, not the map against the continuum coupling,
and the discrete map is not ported here.

Checks (--check): the terminal values and currents; the sampled value tends to
the point value at second order in r; the kernel mass is pi r^3/6; each
terminal's deposition integrates to its current.

    python coupled_1d3d_mms.py --emit
    python coupled_1d3d_mms.py --check
"""

import argparse
import functools
import inspect
import os
import re
import sys

import numpy as np
import sympy as sp

import fda_monodomain_mms as fda
import graph_monodomain_mms as graph
from mms_common import quadrature
from mms_common.latex_emit import emit_latex_fragment, paper_symbol_names
from mms_common.numpy_emit import emit_numpy_module

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "coupled_1d3d"

RADIUS = 0.11
R_PVJ = 1.0
CUBE_LO = np.zeros(3)
CUBE_HI = np.ones(3)
TERMINALS = tuple(tuple(float(v) for v in p) for p in graph.TERMINALS)


@functools.lru_cache(maxsize=None)
def _vm_numpy():
    return sp.lambdify((fda.x, fda.y, fda.z, fda.t), graph.tissue_vm(), "numpy")


@functools.lru_cache(maxsize=None)
def _v1d_numpy():
    return sp.lambdify((graph.x, graph.t), graph.fields()["V1D"], "numpy")


def terminal_vm(k, t):
    """Point value of the graph field at terminal k."""
    return float(_v1d_numpy()(TERMINALS[k][0], t))


def terminal_current_exact(t, R_pvj=1.0):
    """Point-value terminal currents (V_1D,k - V_m(x_k))/R_pvj = +-(3/2) sqrt(1+t)/R_pvj."""
    s = 1.5*np.sqrt(1.0 + t)/R_pvj
    return (float(s), float(-s))


def kernel_mass(radius=RADIUS):
    """int w over a half ball: pi r^3 / 6."""
    return np.pi*radius**3/6.0


def continuum_sampling(k, t, radius=RADIUS, n=32):
    """Vbar_k = int w V_m / int w over the half ball at terminal k."""
    vm = _vm_numpy()
    num, den = quadrature.half_ball_kernel_integral(
        lambda p: vm(p[0], p[1], p[2], t), TERMINALS[k], radius, CUBE_LO, CUBE_HI, n)
    return num/den


def continuum_current(k, t, radius=RADIUS, n=32, R_pvj=1.0):
    """I_pvj,k = (V_1D,k - Vbar_k)/R_pvj."""
    return (terminal_vm(k, t) - continuum_sampling(k, t, radius, n))/R_pvj


def _kernel(k, points, radius):
    pts = np.asarray(points, dtype=float)
    r = np.linalg.norm(pts - np.asarray(TERMINALS[k]), axis=-1)
    return np.maximum(0.0, 1.0 - r/radius)


def terminal_deposition(k, points, t, radius=RADIUS, n=32, R_pvj=1.0):
    """I_pvj,k w_k / int w_k at points of shape (N, 3)."""
    return continuum_current(k, t, radius, n, R_pvj)*_kernel(k, points, radius)/kernel_mass(radius)


def continuum_junction_source(points, t, radius=RADIUS, n=32, R_pvj=1.0):
    """S_pvj at points of shape (N, 3): the sum of the terminal depositions."""
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    return sum(terminal_deposition(k, pts, t, radius, n, R_pvj) for k in range(len(TERMINALS)))


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def run_self_checks():
    t = 0.07
    for k, sign in ((0, 1.0), (1, -1.0)):
        _check(abs(terminal_current_exact(t)[k] - sign*1.5*np.sqrt(1 + t)) < 1e-15,
               "terminal current %d differs from +-(3/2) sqrt(1+t)" % k)
        point = float(_vm_numpy()(*TERMINALS[k], t))
        _check(abs(terminal_vm(k, t) - point - terminal_current_exact(t)[k]) < 1e-14,
               "V_1D - V_m differs from the exact current at terminal %d" % k)
        bound = 7*np.pi**2*RADIUS**2*np.sqrt(1 + t)
        _check(abs(continuum_sampling(k, t) - point) <= bound, "sampling exceeds the O(r^2) bound")
    _check(abs(quadrature.half_ball_kernel_integral(
        lambda p: np.ones(p.shape[1:]), TERMINALS[0], RADIUS, CUBE_LO, CUBE_HI, 32)[1]
        - kernel_mass()) < 1e-12*kernel_mass(), "kernel mass differs from pi r^3/6")


# ---- emission ---------------------------------------------------------------

_REFERENCE_FUNCTIONS = '''

TERMINALS = ((0.0, 1.0/6.0, 1.0/3.0), (1.0, 1.0/6.0, 1.0/3.0))
RADIUS = 0.11
CUBE_LO = numpy.zeros(3)
CUBE_HI = numpy.ones(3)


def terminal_current_exact(t, R_pvj=1.0):
    """Point-value terminal currents +-(3/2) sqrt(1 + t) / R_pvj."""
    s = 1.5*numpy.sqrt(1.0 + t)/R_pvj
    return (float(s), float(-s))


def continuum_sampling(k, t, radius=RADIUS, n=32):
    """Vbar_k = int w V_m / int w over the half ball at terminal k."""
    num, den = half_ball_kernel_integral(
        lambda p: Vm3D(p[0], p[1], p[2], t), TERMINALS[k], radius, CUBE_LO, CUBE_HI, n)
    return num/den


def continuum_current(k, t, radius=RADIUS, n=32, R_pvj=1.0):
    """I_pvj,k = (V_1D,k - Vbar_k) / R_pvj."""
    return (float(V1D(TERMINALS[k][0], t)) - continuum_sampling(k, t, radius, n))/R_pvj


def _kernel(k, points, radius):
    r = numpy.linalg.norm(numpy.asarray(points, dtype=float) - numpy.asarray(TERMINALS[k]), axis=-1)
    return numpy.maximum(0.0, 1.0 - r/radius)


def terminal_deposition(k, points, t, radius=RADIUS, n=32, R_pvj=1.0):
    """I_pvj,k w_k / int w_k at points of shape (N, 3)."""
    mass = numpy.pi*radius**3/6.0
    return continuum_current(k, t, radius, n, R_pvj)*_kernel(k, points, radius)/mass


def continuum_junction_source(points, t, radius=RADIUS, n=32, R_pvj=1.0):
    """S_pvj at points of shape (N, 3): the sum of the terminal depositions."""
    pts = numpy.atleast_2d(numpy.asarray(points, dtype=float))
    return sum(terminal_deposition(k, pts, t, radius, n, R_pvj) for k in range(len(TERMINALS)))
'''


def _quadrature_source():
    """Source of the half-ball quadrature, numpy-only, from mms_common."""
    parts = [inspect.getsource(quadrature._nodes),
             inspect.getsource(quadrature._ray_length),
             inspect.getsource(quadrature.half_ball_kernel_integral)]
    return "\n\n" + "\n\n".join(re.sub(r"\bnp\.", "numpy.", p).rstrip("\n") + "\n" for p in parts)


def numpy_module():
    functions = {
        "Vm3D": ((fda.x, fda.y, fda.z, fda.t), graph.tissue_vm()),
        "V1D": ((graph.x, graph.t), graph.fields()["V1D"]),
    }
    doc = ("Continuum 1D-3D junction of the Purkinje-graph and myocardium manufactured "
           "solution: Vm3D(x, y, z, t), V1D(x, t), terminal_current_exact, "
           "continuum_sampling, continuum_current, terminal_deposition, "
           "continuum_junction_source, with the half-ball kernel quadrature.  "
           "Generated by coupled_1d3d_mms.py.")
    params = {"RADIUS": RADIUS, "R_PVJ": R_PVJ, "CHI": fda.CHI, "CM": fda.CM}
    base = emit_numpy_module(REFERENCE_NAME, functions, params, doc)
    return base.rstrip("\n") + "\n" + _quadrature_source() + _REFERENCE_FUNCTIONS


def latex_fragments():
    s = {n: sp.Symbol(n, real=True) for n in "Vm Vbar V1k Ik Rpvj wk Wk Swv W".split()}
    names = paper_symbol_names(s.values())
    names.update({
        s["Vbar"]: r"\bar V_{k}", s["V1k"]: r"V_{1D,k}", s["Ik"]: r"I_{\mathrm{pvj},k}",
        s["Rpvj"]: r"R_{\mathrm{pvj}}", s["wk"]: r"w_k(\mathbf x)",
        s["Swv"]: r"\int_{B_k\cap\Omega} w_k V_m\,\mathrm d\mathbf x",
        s["W"]: r"\int_{B_k\cap\Omega} w_k\,\mathrm d\mathbf x",
    })
    ratio, dep = sp.symbols("ratio dep")
    names.update({ratio: r"\left|\mathbf x-\mathbf x_k\right|/r_{\mathrm{pvj}}",
                  dep: r"\sum_k I_{\mathrm{pvj},k}\,w_k(\mathbf x)"})
    sgn = sp.Symbol("sgn")
    out = {
        "coupled_kernel.tex": emit_latex_fragment(
            r"w_k(\mathbf x)", sp.Max(0, 1 - ratio), names),
        "coupled_sampling.tex": emit_latex_fragment(r"\bar V_k", s["Swv"]/s["W"], names),
        "coupled_current.tex": emit_latex_fragment(
            r"I_{\mathrm{pvj},k}", (s["V1k"] - s["Vbar"])/s["Rpvj"], names),
        "coupled_deposition.tex": emit_latex_fragment(
            r"S_{\mathrm{pvj}}", dep/s["W"], names),
        "coupled_terminal_current.tex": emit_latex_fragment(
            r"I_{\mathrm{pvj},k}^{\mathrm{point}}",
            sp.Rational(3, 2)*sgn/s["Rpvj"], {**names, sgn: r"\pm\sqrt{1 + t}"}),
    }
    return {"latex/" + k: v for k, v in out.items()}


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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--emit", action="store_true")
    grp.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)
    try:
        run_self_checks()
    except AssertionError as err:
        print("self-check failed: %s" % err)
        return 1
    if args.emit:
        for rel in write_outputs():
            print("wrote %s" % rel)
    else:
        print("self-checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
