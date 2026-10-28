#!/usr/bin/env python3
"""FDA monodomain manufactured solution on the unit cube, d = 1, 2, 3.

Fields on [0, 1]^d, with t in units of 1 s and x in units of 1 m:

    F_1 = cos(pi x)   F_2 = F_1 cos(2 pi y)   F_3 = F_2 cos(3 pi z)
    G_1 = 1 + x       G_2 = 1 + x y^2         G_3 = 1 + x y^2 z^3
    V_m = sqrt(1 + t) F_d
    u_1 = (1 + t) G_d + V_m,  u_2 = G_d^(-1/2) / (1 + t),  u_3 = 0

The ionic current of the manufactured model is
J_ion = -(C_m/2)(u_1 + u_3 - V_m) u_2^2 (V_m - u_3) + (beta/chi)(V_m - u_3),
the quantity the ionic model writes as Iion; the myocardium solver receives
I_ion = J_ion/C_m and solves
chi C_m dV_m/dt = div(sigma grad V_m) - chi C_m I_ion, with C_m = 2 and chi = 3.

F_d is an eigenfunction of the diagonal operator,
div(sigma grad F_d) = -k_d F_d with k_d = sum_{m<=d} sigma_mm (m pi)^2, and
(u_1 - V_m) u_2^2 = 1/(1 + t), so the PDE residual is exactly (k_d + beta_d) V_m:
beta_d is -k_d rounded to two significant digits, and no applied-current source
is needed.  sigma = diag(0.111453302, 0.121585420, 0.030396355) is given to nine
significant digits and is tuned so that k_d equals -beta_d to that precision: the
residual is nonzero but below RESIDUAL_RTOL times |div(sigma grad V_m)|.

Checks (--self-check): the ionic rates of u_1 and u_2 hold on the exact fields;
the residual equals (k_d + beta_d) V_m symbolically and |k_d + beta_d|/k_d is at
most RESIDUAL_RTOL; n.grad V_m vanishes on every face of the cube.

--write PATH rewrites, and --check PATH compares, the generated regions of the
hand-written manufacturedFDAReference.H: F_d (computeF and manufacturedFDAFAt),
G_d (computeG), V_m (computeManufacturedV), u_1, u_2, u_3
(computeManufacturedU) and J_ion (computeIion).  Each sits between
`// BEGIN GENERATED <name>` and `// END GENERATED <name>`; the rest of the
file is left as it is.

    python fda_monodomain_mms.py --emit
    python fda_monodomain_mms.py --self-check
    python fda_monodomain_mms.py --write ../../../src/verificationModels/monodomainVerification/manufacturedFDAReference.H
    python fda_monodomain_mms.py --check ../../../src/verificationModels/monodomainVerification/manufacturedFDAReference.H
"""

import functools
import os
import sys

import sympy as sp

import mms_codegen as codegen
from mms_common.latex_emit import emit_latex_fragment, paper_symbol_names
from mms_common.numpy_emit import emit_numpy_module

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "fda_monodomain"
DIMENSIONS = (1, 2, 3)

BETA = {1: -1.1, 2: -5.9, 3: -8.6}
# sigma carries 9 significant digits, so k_d + beta_d is rounding-level (<= 3e-9 k_d)
RESIDUAL_RTOL = 1e-8
CM = 2
CHI = 3
SIGMA_DIAG = tuple(sp.Rational(s) for s in ("0.111453302", "0.121585420", "0.030396355"))

x, y, z, t = sp.symbols("x y z t", real=True)
AXES = (x, y, z)


def _beta(d):
    return sp.Rational(str(BETA[d]))


def eigenvalue(d):
    """k_d with div(sigma grad F_d) = -k_d F_d."""
    return sum(SIGMA_DIAG[m]*(m + 1)**2*sp.pi**2 for m in range(d))


def _F(d):
    f = sp.cos(sp.pi*x)
    for m, ax in ((2, y), (3, z)):
        if d >= m:
            f *= sp.cos(m*sp.pi*ax)
    return f


def _G(d):
    return {1: 1 + x, 2: 1 + x*y**2, 3: 1 + x*y**2*z**3}[d]


def ionic_rates(V, u1, u2, u3):
    """Right-hand sides of the auxiliary-state equations (V_m held fixed)."""
    a = u1 + u3 - V
    return {"u1": a**2*u2**2 + sp.Rational(1, 2)*a*u2**2*(V - u3),
            "u2": -a*u2**3,
            "u3": sp.Integer(0)}


def ionic_current(V, u1, u2, u3, cm, beta, chi):
    """J_ion, the model's Iion (not yet divided by C_m)."""
    return -(cm/2)*(u1 + u3 - V)*u2**2*(V - u3) + (beta/chi)*(V - u3)


def vm_of(F):
    """V_m = sqrt(1 + t) F."""
    return sp.sqrt(1 + t)*F


def u_of(F, G):
    """(u_1, u_2, u_3) of the spatial factors F and G."""
    return (1 + t)*G + vm_of(F), G**sp.Rational(-1, 2)/(1 + t), sp.Integer(0)


@functools.lru_cache(maxsize=None)
def _fields(d):
    F, G = _F(d), _G(d)
    V = vm_of(F)
    u1, u2, u3 = u_of(F, G)
    Iion = ionic_current(V, u1, u2, u3, sp.Integer(CM), _beta(d), sp.Integer(CHI))
    return {"Vm": V, "u1": u1, "u2": u2, "u3": u3, "Iion": Iion, "F": F, "G": G}


def fields(d):
    """{Vm, u1, u2, u3, Iion, F, G} as sympy expressions in x, y, z, t.

    Iion is J_ion; the solver's current is Iion/C_m.
    """
    if d not in DIMENSIONS:
        raise ValueError("d must be 1, 2 or 3, got %r" % (d,))
    return dict(_fields(d))


def pde_residual(d):
    """chi C_m dV/dt - div(sigma grad V) + chi C_m (J_ion/C_m)."""
    f = fields(d)
    V = f["Vm"]
    div = sum(SIGMA_DIAG[m]*sp.diff(V, AXES[m], 2) for m in range(d))
    return CHI*CM*sp.diff(V, t) - div + CHI*CM*(f["Iion"]/CM)


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def run_self_checks():
    for d in DIMENSIONS:
        f = fields(d)
        V, u1, u2, u3 = f["Vm"], f["u1"], f["u2"], f["u3"]
        rates = ionic_rates(V, u1, u2, u3)
        for name, u in (("u1", u1), ("u2", u2), ("u3", u3)):
            _check(sp.simplify(sp.diff(u, t) - rates[name]) == 0,
                   "d=%d: d%s/dt differs from its ionic rate" % (d, name))
        k = eigenvalue(d)
        _check(sp.simplify(pde_residual(d) - (k + _beta(d))*V) == 0,
               "d=%d: residual is not (k + beta) V" % d)
        gap = abs(float(k) + BETA[d])/float(k)
        _check(gap <= RESIDUAL_RTOL,
               "d=%d: |k + beta|/k = %.3e exceeds %.1e" % (d, gap, RESIDUAL_RTOL))
        for m in range(d):
            for face in (0, 1):
                flux = sp.diff(V, AXES[m]).subs(AXES[m], face)
                _check(sp.simplify(flux) == 0,
                       "d=%d: n.grad V_m does not vanish on %s=%d" % (d, AXES[m], face))


def _symbols():
    names = "Vm u1 u2 u3 cm chi beta Fd Gd Jion Iion".split()
    return {n: sp.Symbol(n, real=True) for n in names}


def _latex_names(s):
    names = paper_symbol_names(s.values())
    names.update({s["Fd"]: "F_d", s["Gd"]: "G_d",
                  s["Jion"]: r"J_{\mathrm{ion}}", s["Iion"]: r"I_{\mathrm{ion}}"})
    return names


def latex_fragments():
    s = _symbols()
    names = _latex_names(s)
    out = {}
    for d in DIMENSIONS:
        out["fda_F_d%d.tex" % d] = emit_latex_fragment("F_%d" % d, _F(d), names)
        out["fda_G_d%d.tex" % d] = emit_latex_fragment("G_%d" % d, _G(d), names)
        out["fda_beta_d%d.tex" % d] = emit_latex_fragment(
            r"\beta_%d" % d, _beta(d), names)
    out["fda_Vm.tex"] = emit_latex_fragment("V_m", sp.sqrt(1 + t)*s["Fd"], names)
    out["fda_u1.tex"] = emit_latex_fragment("u_1", (1 + t)*s["Gd"] + s["Vm"], names)
    out["fda_u2.tex"] = emit_latex_fragment("u_2", s["Gd"]**sp.Rational(-1, 2)/(1 + t), names)
    out["fda_u3.tex"] = emit_latex_fragment("u_3", sp.Integer(0), names)
    rates = ionic_rates(s["Vm"], s["u1"], s["u2"], s["u3"])
    out["fda_rate_u1.tex"] = emit_latex_fragment(r"\dot u_1", rates["u1"], names)
    out["fda_rate_u2.tex"] = emit_latex_fragment(r"\dot u_2", rates["u2"], names)
    out["fda_Jion.tex"] = emit_latex_fragment(
        r"J_{\mathrm{ion}}",
        ionic_current(s["Vm"], s["u1"], s["u2"], s["u3"], s["cm"], s["beta"], s["chi"]), names)
    out["fda_Iion.tex"] = emit_latex_fragment(
        r"I_{\mathrm{ion}}", s["Jion"]/s["cm"], names)
    return {"latex/" + k: v for k, v in out.items()}


def numpy_module():
    functions = {}
    for d in DIMENSIONS:
        for name, expr in fields(d).items():
            functions["%s_d%d" % (name, d)] = ((x, y, z, t), expr)
    doc = ("FDA monodomain manufactured fields on the unit cube; <name>_d<d>(x, y, z, t) "
           "for Vm, u1, u2, u3, Iion (J_ion, the model's Iion), F, G.  Generated by "
           "fda_monodomain_mms.py.")
    params = {"BETA": {d: BETA[d] for d in DIMENSIONS}, "CM": CM, "CHI": CHI,
              "SIGMA_DIAG": [float(v) for v in SIGMA_DIAG]}
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


def _at(expr, xs, ys, zs):
    """expr with x, y, z replaced by C++ symbols."""
    return expr.subs({x: _sym(xs), y: _sym(ys), z: _sym(zs)}, simultaneous=True)


def fda_header_regions(text):
    """{region name: C++ lines} of manufacturedFDAReference.H."""
    ind = functools.partial(codegen.region_indent, text)
    out = {}
    for d in DIMENSIONS:
        name = "computeF_d%d" % d
        out[name] = codegen.assignment("F[i]", _at(_F(d), "x[i]", "y[i]", "z[i]"), ind(name))
        name = "manufacturedFDAFAt_d%d" % d
        out[name] = codegen.statement(
            "return", _at(_F(d), "p.x()", "p.y()", "p.z()"), ind(name))
        name = "computeG_d%d" % d
        out[name] = codegen.assignment("G[i]", _at(_G(d), "x[i]", "y[i]", "z[i]"), ind(name))
    F, G = _sym("Fx[i]"), _sym("Gx[i]")
    out["computeManufacturedV"] = codegen.assignment(
        "V[i]", vm_of(F), ind("computeManufacturedV"))
    u = u_of(F, G)
    out["computeManufacturedU"] = sum(
        (codegen.assignment("u%d[i]" % (m + 1), u[m], ind("computeManufacturedU"))
         for m in range(3)), [])
    V, u1, u2, u3 = (_sym(n) for n in ("V[i]", "u1[i]", "u2[i]", "u3[i]"))
    out["computeIion"] = codegen.assignment(
        "Iion[i]", ionic_current(V, u1, u2, u3, _sym("Cm"), _sym("beta"), _sym("chi")),
        ind("computeIion"))
    return out


HEADER_PATH = ("../../../src/verificationModels/monodomainVerification/"
               "manufacturedFDAReference.H")
HEADER_PATHS = {"manufacturedFDAReference.H": HEADER_PATH}
HEADERS = {"manufacturedFDAReference.H": fda_header_regions}


def header_matches(path):
    return codegen.regions_match(__file__, path, HEADERS[os.path.basename(path)])


def main(argv=None):
    return codegen.generator_main(argv, __doc__, __file__, run_self_checks, write_outputs,
                                  HEADERS)


if __name__ == "__main__":
    sys.exit(main())
