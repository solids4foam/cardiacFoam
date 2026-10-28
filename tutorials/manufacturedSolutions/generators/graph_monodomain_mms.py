#!/usr/bin/env python3
"""Reaction-diffusion Purkinje-graph manufactured solution (one dimension).

On a graph edge with arc coordinate x in [0, 1] (t in units of 1 s, x in units
of 1 m) the field is the d = 1 field of fda_monodomain_mms:

    V_1D = sqrt(1 + t) cos(pi x)
    u_1 = (1 + t)(1 + x) + V_1D,  u_2 = (1 + x)^(-1/2) / (1 + t),  u_3 = 0

with the graph conductivity g = G_GRAPH = 0.111453302 (the first diagonal entry
of the FDA tensor), chi = 3, C_m = 2 and the FDA ionic closure with beta_1.  The
graph equation is

    chi C_m dV_1D/dt = d/dx(g dV_1D/dx) - chi C_m I_ion,1D,

where I_ion,1D = J_ion/C_m is the solver's ionic current.  The solver applies
no bulk source on the graph: the FDA closure makes the residual of this
equation vanish up to the rounding-level (g pi^2 + beta_1) V_1D, with
|g pi^2 + beta_1| <= 1e-8 |beta_1| because g is tuned to the eigenvalue of
cos(pi x).  That residual is kept only as a check (R_1D, residual_1D), not as
an applied term.

The only term added on the graph belongs to the coupled problem: in
bidirectional transfer the terminal node k receives -I_pvj,k/(pi rho^2 L), with
rho the fibre radius and L the edge length (node volume pi rho^2 L); in
unidirectional transfer it is zero.

The graph has a root (1/2, 1/6, 1/3) and two straight branches to the terminals
(0, 1/6, 1/3) and (1, 1/6, 1/3).  At a terminal the tissue field
V_m = sqrt(1 + t) cos(pi x) cos(2 pi y) cos(3 pi z) takes the value
-V_1D/2 and has dV_m/dx = 0.

Checks (--check): the graph equation holds without a source up to the residual
(g pi^2 + beta_1) V_1D, which is rounding-level; the ionic rates hold; V_m = -V_1D/2 and dV_m/dx = 0 at both
terminals.

    python graph_monodomain_mms.py --emit
    python graph_monodomain_mms.py --check
"""

import argparse
import functools
import os
import sys

import sympy as sp

import fda_monodomain_mms as fda
from mms_common.latex_emit import emit_latex_fragment, paper_symbol_names
from mms_common.numpy_emit import emit_numpy_module

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "graph_monodomain"

x, t = fda.x, fda.t
G_GRAPH = 0.111453302
ROOT = (sp.Rational(1, 2), sp.Rational(1, 6), sp.Rational(1, 3))
TERMINALS = ((sp.Integer(0), sp.Rational(1, 6), sp.Rational(1, 3)),
             (sp.Integer(1), sp.Rational(1, 6), sp.Rational(1, 3)))


def graph_conductivity():
    return sp.Rational(str(G_GRAPH))


@functools.lru_cache(maxsize=None)
def _fields():
    f = fda.fields(1)
    V = f["Vm"]
    g = graph_conductivity()
    # chi C_m dV/dt - d/dx(g dV/dx) + chi C_m (J_ion / C_m)
    residual = (fda.CHI*fda.CM*sp.diff(V, t) - g*sp.diff(V, x, 2)
                + fda.CHI*fda.CM*(f["Iion"]/fda.CM))
    return {"V1D": V, "u1": f["u1"], "u2": f["u2"], "u3": f["u3"],
            "Iion": f["Iion"], "residual_1D": sp.simplify(residual)}


def fields():
    """{V1D, u1, u2, u3, Iion, residual_1D} in x, t; Iion is J_ion (solver
    current Iion/C_m), residual_1D the graph-equation residual (no source is
    applied)."""
    return dict(_fields())


def graph_residual():
    """chi C_m dV/dt - d/dx(g dV/dx) + chi C_m I_ion,1D, which is residual_1D."""
    f = fields()
    g = graph_conductivity()
    return (fda.CHI*fda.CM*sp.diff(f["V1D"], t) - g*sp.diff(f["V1D"], x, 2)
            + fda.CHI*fda.CM*(f["Iion"]/fda.CM))


def tissue_vm():
    """The 3D field V_m (FDA d = 3), in x, y, z, t."""
    return fda.fields(3)["Vm"]


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def run_self_checks():
    f = fields()
    _check(sp.simplify(graph_residual() - f["residual_1D"]) == 0,
           "graph_residual differs from residual_1D")
    k1 = graph_conductivity()*sp.pi**2
    _check(sp.simplify(f["residual_1D"] - (k1 + sp.Rational(str(fda.BETA[1])))*f["V1D"]) == 0,
           "the residual differs from (g pi^2 + beta_1) V_1D")
    _check(abs(float(k1) + fda.BETA[1]) <= 1e-8*abs(fda.BETA[1]),
           "g pi^2 + beta_1 is not rounding-level")
    _check(G_GRAPH == float(fda.SIGMA_DIAG[0]), "G_GRAPH differs from sigma_xx")
    rates = fda.ionic_rates(f["V1D"], f["u1"], f["u2"], f["u3"])
    for name in ("u1", "u2"):
        _check(sp.simplify(sp.diff(f[name], t) - rates[name]) == 0,
               "d%s/dt differs from its ionic rate" % name)
    Vm = tissue_vm()
    for tx, ty, tz in TERMINALS:
        at = {fda.x: tx, fda.y: ty, fda.z: tz}
        _check(sp.simplify(Vm.subs(at) + f["V1D"].subs(x, tx)/2) == 0,
               "V_m != -V_1D/2 at a terminal")
        _check(sp.simplify(sp.diff(Vm, fda.x).subs(at)) == 0,
               "dV_m/dx != 0 at a terminal")


def latex_fragments():
    s = {n: sp.Symbol(n, real=True) for n in "V1 beta g Ik vol".split()}
    names = paper_symbol_names(s.values())
    names.update({s["V1"]: r"V_{1D}", s["g"]: r"g",
                  s["Ik"]: r"I_{\mathrm{pvj},k}", s["vol"]: r"\pi\rho^{2}L"})
    V1 = sp.sqrt(1 + t)*sp.cos(sp.pi*x)
    out = {
        "graph_V1D.tex": emit_latex_fragment(r"V_{1D}", V1, names),
        "graph_residual.tex": emit_latex_fragment(
            r"R_{1D}", (s["g"]*sp.pi**2 + s["beta"])*s["V1"], names),
        "graph_terminal_source.tex": emit_latex_fragment(
            r"I_{\mathrm{app},1D}(\mathbf x_k)", -s["Ik"]/s["vol"], names),
        "graph_terminal_Vm.tex": emit_latex_fragment(r"V_m(\mathbf x_k)", -s["V1"]/2, names),
    }
    return {"latex/" + k: v for k, v in out.items()}


def numpy_module():
    f = fields()
    functions = {name: ((x, t), f[name]) for name in ("V1D", "u1", "u2", "u3", "Iion", "residual_1D")}
    doc = ("Purkinje-graph manufactured fields; <name>(x, t) for V1D, u1, u2, u3, "
           "Iion (J_ion), residual_1D (graph-equation residual; no source is applied).  Generated by "
           "graph_monodomain_mms.py.")
    params = {"G_GRAPH": G_GRAPH, "CM": fda.CM, "CHI": fda.CHI, "BETA": fda.BETA[1],
              "ROOT": [float(v) for v in ROOT],
              "TERMINALS": [[float(v) for v in p] for p in TERMINALS]}
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
