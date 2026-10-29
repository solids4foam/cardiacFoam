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

Checks (--self-check): K is symmetric positive definite with nonzero
off-diagonal entries; grad F vanishes on every face; the source equals the
six-term Hessian sum (cxx_source); the ionic rates of u_1 and u_2 hold on the
exact fields.

--write PATH rewrites, and --check PATH compares, the generated regions of the
hand-written manufacturedAnisotropicMonodomainReference.H: F, V_m, u_1, u_2,
u_3 and S, with K and beta the runtime arguments of the C++ (K need not be
symmetric there: div(K grad F) is written for a general constant tensor).

    python rotated_monodomain_mms.py --emit
    python rotated_monodomain_mms.py --self-check
    python rotated_monodomain_mms.py --write ../../../src/verificationModels/monodomainVerification/manufacturedAnisotropicMonodomainReference.H
    python rotated_monodomain_mms.py --check ../../../src/verificationModels/monodomainVerification/manufacturedAnisotropicMonodomainReference.H
"""

import functools
import os
import sys

import sympy as sp

import fda_monodomain_mms as fda
import mms_codegen as codegen
from mms_common.latex_emit import emit_latex_fragment, paper_symbol_names, printed_constant
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


def source_of(F, div_K_grad_F, beta):
    """S = sqrt(1 + t) (beta F - div(K grad F))."""
    return fda.vm_of(beta*F - div_K_grad_F)


def pde_source(F, K, beta):
    """S for a constant tensor K."""
    return source_of(F, div_K_grad(F, sp.Matrix(K)), beta)


@functools.lru_cache(maxsize=None)
def _fields(K):
    F, G = _F(), _G()
    V = fda.vm_of(F)
    u1, u2, u3 = fda.u_of(F, G)
    S = pde_source(F, K, sp.Rational(str(BETA)))
    return {"Vm": V, "u1": u1, "u2": u2, "u3": u3, "F": F, "S": S}


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
        "rotated_K.tex": emit_latex_fragment(
            r"\boldsymbol\sigma", K_CASE.applyfunc(lambda v: printed_constant(v, 6)), names),
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


def _sym(name):
    return sp.Symbol(name, real=True)


def _at_centre(expr):
    """expr with x, y, z replaced by the cell centre c."""
    return expr.subs({x: _sym("c.x()"), y: _sym("c.y()"), z: _sym("c.z()")},
                     simultaneous=True)


# The tensor argument of the C++ source, a general (not necessarily
# symmetric) constant tensor.
K_CXX = sp.Matrix(3, 3, lambda i, j: _sym("K.%s%s()" % ("xyz"[i], "xyz"[j])))


def anisotropic_header_regions(text):
    """{region name: C++ lines} of manufacturedAnisotropicMonodomainReference.H."""
    ind = lambda name: codegen.region_indent(text, name)
    out = {}
    name = "computeAnisotropicManufacturedF"
    out[name] = codegen.assignment("F[cellI]", _at_centre(_F()), ind(name))
    F = _sym("F[cellI]")
    name = "computeAnisotropicManufacturedV"
    out[name] = codegen.assignment("V[cellI]", fda.vm_of(F), ind(name))
    name = "computeAnisotropicManufacturedU"
    G = _sym("G")
    u = fda.u_of(F, G)
    out[name] = codegen.assignment("G", _at_centre(_G()), ind(name), decl="const scalar")
    for m in range(3):
        out[name] += codegen.assignment("u%d[cellI]" % (m + 1), u[m], ind(name))
    name = "computeAnisotropicManufacturedPdeSource"
    Fc, D = _at_centre(_F()), _at_centre(div_K_grad(_F(), K_CXX))
    lines, (Fr, Dr) = codegen.cse_lines([Fc, D], ind(name))
    lines += codegen.assignment("F", Fr, ind(name), decl="const scalar")
    lines += codegen.assignment("divKGradF", Dr, ind(name), decl="const scalar")
    lines += codegen.assignment(
        "source[cellI]", source_of(_sym("F"), _sym("divKGradF"), _sym("beta")), ind(name))
    out[name] = lines
    return out


HEADER_PATH = ("../../../src/verificationModels/monodomainVerification/"
               "manufacturedAnisotropicMonodomainReference.H")
HEADER_PATHS = {"manufacturedAnisotropicMonodomainReference.H": HEADER_PATH}
HEADERS = {"manufacturedAnisotropicMonodomainReference.H": anisotropic_header_regions}


def header_matches(path):
    return codegen.regions_match(__file__, path, HEADERS[os.path.basename(path)])


def main(argv=None):
    return codegen.generator_main(argv, __doc__, __file__, run_self_checks, write_outputs,
                                  HEADERS)


if __name__ == "__main__":
    sys.exit(main())
