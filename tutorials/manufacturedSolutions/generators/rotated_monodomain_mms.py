#!/usr/bin/env python3
"""Monodomain manufactured solution on the unit cube with a full conductivity tensor.

Fields on [0, 1]^3, with t in units of 1 s and x in units of 1 m:

    F   = sin^2(pi x) sin^2(pi y) sin^2(pi z)
    V_m = sqrt(1 + t) F
    u_1 = (1 + t) G + V_m,  u_2 = G^(-1/2) / (1 + t),  u_3 = 0,
    G   = 1 + x y^2 z^3
    S   = sqrt(1 + t) (beta F - div(K grad F))

with the constant symmetric tensor K of the tutorial (K_CASE, off-diagonal
entries nonzero) and beta = -8.6, the three-dimensional constant of the FDA
ionic closure.  The ionic closure is that of fda_monodomain_mms; S is the
applied-current source the solver adds on the right-hand side, so that the
residual of the monodomain equation vanishes.

The gradient of F vanishes on every face of the cube, so the conormal flux
K grad V_m . n is zero because of the field and not because of the operator:
the boundary flux reconstruction is not exercised by this case.

Checks (--check): K is symmetric positive definite with nonzero off-diagonal
entries; grad F vanishes on every face; the source equals the six-term Hessian
sum of the C++ reference; the ionic rates of u_1 and u_2 hold on the exact
fields.

    python rotated_monodomain_mms.py --emit
    python rotated_monodomain_mms.py --check
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
REFERENCE_NAME = "rotated_monodomain"

x, y, z, t = fda.x, fda.y, fda.z, fda.t
AXES = (x, y, z)
BETA = fda.BETA[3]

_K = [sp.Rational(s) for s in
      ("0.106875", "-0.0084931", "-0.022561", "0.116682", "-0.0130256", "0.0398782")]
# symmTensor order xx xy xz yy yz zz
K_CASE = sp.Matrix([[_K[0], _K[1], _K[2]],
                    [_K[1], _K[3], _K[4]],
                    [_K[2], _K[4], _K[5]]])


def _F():
    return sp.sin(sp.pi*x)**2*sp.sin(sp.pi*y)**2*sp.sin(sp.pi*z)**2


def _G():
    return 1 + x*y**2*z**3


def div_K_grad(F, K):
    """div(K grad F) for a constant tensor K."""
    grad = sp.Matrix([sp.diff(F, a) for a in AXES])
    flux = K*grad
    return sum(sp.diff(flux[i], AXES[i]) for i in range(3))


@functools.lru_cache(maxsize=None)
def _fields(K):
    F, G = _F(), _G()
    V = sp.sqrt(1 + t)*F
    u1 = (1 + t)*G + V
    u2 = G**sp.Rational(-1, 2)/(1 + t)
    S = sp.sqrt(1 + t)*(sp.Rational(str(BETA))*F - div_K_grad(F, sp.Matrix(K)))
    return {"Vm": V, "u1": u1, "u2": u2, "u3": sp.Integer(0), "F": F, "S": S}


def fields(K):
    """{Vm, u1, u2, u3, F, S} as sympy expressions in x, y, z, t for a
    constant symmetric 3x3 tensor K."""
    K = sp.ImmutableMatrix(K)
    if K.shape != (3, 3) or K != K.T:
        raise ValueError("K must be a symmetric 3x3 matrix")
    return dict(_fields(K))


def cxx_source(K):
    """The source in the form of the C++ reference: six-term Hessian sum."""
    pi = sp.pi
    sx2, sy2, sz2 = (sp.sin(pi*a)**2 for a in AXES)
    dx, dy, dz = (pi*sp.sin(2*pi*a) for a in AXES)
    dxx, dyy, dzz = (2*pi**2*sp.cos(2*pi*a) for a in AXES)
    div = (K[0, 0]*dxx*sy2*sz2 + K[1, 1]*sx2*dyy*sz2 + K[2, 2]*sx2*sy2*dzz
           + (K[0, 1] + K[1, 0])*dx*dy*sz2
           + (K[0, 2] + K[2, 0])*dx*sy2*dz
           + (K[1, 2] + K[2, 1])*sx2*dy*dz)
    return sp.sqrt(1 + t)*(sp.Rational(str(BETA))*sx2*sy2*sz2 - div)


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def run_self_checks():
    K = K_CASE
    _check(K == K.T and all(K[i, j] != 0 for i in range(3) for j in range(i)),
           "K is not symmetric with nonzero off-diagonal entries")
    _check(all(K[:k, :k].det() > 0 for k in (1, 2, 3)), "K is not positive definite")
    f = fields(K)
    for a in AXES:
        for face in (0, 1):
            g = sp.diff(f["F"], a).subs(a, face)
            _check(sp.simplify(g) == 0, "d F/d%s does not vanish on %s=%d" % (a, a, face))
    _check(sp.simplify(f["S"] - cxx_source(K)) == 0,
           "source differs from the six-term Hessian sum")
    rates = fda.ionic_rates(f["Vm"], f["u1"], f["u2"], f["u3"])
    for name in ("u1", "u2"):
        _check(sp.simplify(sp.diff(f[name], t) - rates[name]) == 0,
               "d%s/dt differs from its ionic rate" % name)


def latex_fragments():
    s = {n: sp.Symbol(n, real=True) for n in "Vm F Gfun beta Dsig".split()}
    names = paper_symbol_names(s.values())
    names.update({s["F"]: "F", s["Gfun"]: "G_3",
                  s["Dsig"]: r"\nabla\!\cdot\!\left(\boldsymbol\sigma\nabla F\right)"})
    Gf = s["Gfun"]
    out = {
        "rotated_F.tex": emit_latex_fragment("F", _F(), names),
        "rotated_G.tex": emit_latex_fragment("G_3", _G(), names),
        "rotated_Vm.tex": emit_latex_fragment("V_m", sp.sqrt(1 + t)*s["F"], names),
        "rotated_u1.tex": emit_latex_fragment("u_1", (1 + t)*Gf + s["Vm"], names),
        "rotated_u2.tex": emit_latex_fragment("u_2", Gf**sp.Rational(-1, 2)/(1 + t), names),
        "rotated_S.tex": emit_latex_fragment(
            "S", sp.sqrt(1 + t)*(s["beta"]*s["F"] - s["Dsig"]), names),
        "rotated_K.tex": emit_latex_fragment(r"\boldsymbol\sigma", K_CASE, names),
    }
    return {"latex/" + k: v for k, v in out.items()}


def numpy_module():
    f = fields(K_CASE)
    functions = {name: ((x, y, z, t), f[name]) for name in ("Vm", "u1", "u2", "u3", "F", "S")}
    doc = ("Rotated-anisotropy monodomain manufactured fields on the unit cube; "
           "<name>(x, y, z, t) for Vm, u1, u2, u3, F, S with the tensor K in PARAMS.  "
           "Generated by rotated_monodomain_mms.py.")
    params = {"BETA": BETA, "K": [[float(K_CASE[i, j]) for j in range(3)] for i in range(3)]}
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
