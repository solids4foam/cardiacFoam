#!/usr/bin/env python3
"""Lead-field ECG manufactured solutions of the insulated slab and the layered torso.

The box [-1, 2] x [0, 1]^2 holds the heart slab 0 <= x <= 1 between two torso
slabs; every field depends on x and t only, with t in units of 1 s and x in
units of 1 m.  With S = sqrt(1 + t) and k = Gi_xx/(Gi_xx + Ge_xx):

    V_m = S cos(pi x)                                   (heart)
    phi*_slab = -k S + C              x <= 0
              = -k S cos(pi x) + C    0 < x <= 1        C = k S
              =  k S + C              x > 1
    phi*_layered = phi*_slab + alpha R(x),  R' = 1/sigma(x),  R(0) = 0,

with sigma(x) = sigma_2, sigma_1, (Gi + Ge)_xx, sigma_1, sigma_2 on the
layers split at x_L, 0, 1, x_R.  The slab torso is insulated; the layered
torso carries the surface currents n.sigma grad phi = -alpha on x = -1 and
+alpha on x = 2.  The verifier gauges phi* at the centroid of the cell that
holds the reference point, so the constant C cancels.

Checks (--self-check): in the heart (Gi + Ge)_xx phi'' = -Gi_xx V_m'' with
k (Gi + Ge)_xx = Gi_xx, and V_m' vanishes on both walls; phi* is continuous
at x = 0 and x = 1 and linear in each torso slab; R is continuous at x_L, 0, 1,
x_R with sigma R' = 1 in every layer, so the normal current is continuous;
the surface currents are -alpha and +alpha; the slab potential gauged in the
left torso is 0 there and 2 k S in the right torso; with alpha = 0 the
layered potential equals the slab's; k equals 1/sqrt(2) to K_RTOL for both
conductivity pairs of the tutorial; the printed constants round-trip.

The cardiacFoam functions of these solutions (manufacturedLeadFieldECGReference.H
and manufacturedFDABathPhiEAt of manufacturedFDABathBidomainReference.H) carry
no generated regions, so this generator has no --write or --check;
test_slab_leadfield_mms.py compares the emitted NumPy module with
docs/leadField/leadfield_reference.py, which follows those headers branch by
branch.

    python slab_leadfield_mms.py --emit
    python slab_leadfield_mms.py --self-check
"""

import functools
import os
import sys

import sympy as sp

import fda_bidomain_mms as fbm
import mms_codegen as codegen
from mms_common.latex_emit import (emit_cases_fragment, emit_latex_fragment,
                                   emit_parameter_fragment, printed_constant)
from mms_common.numpy_emit import emit_numpy_module

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "slab_leadfield"

# Heart conductivities (xx entries, S/m) of the ECG domain of the tutorial
# manufacturedSolutions/leadFieldECG: the monodomain pair (2 + sqrt 2) G_m and
# sqrt 2 G_m, and the bidomain pair of fda_bidomain_mms.py
MONODOMAIN_PAIR = (sp.Rational("0.380525375260"), sp.Rational("0.157618771260"))
BIDOMAIN_PAIR = (fbm.GI_DIAG[0], fbm.GE_DIAG[0])
K_RTOL = 1e-8

SIGMA_SLAB = sp.Rational("0.0230827345")
X_LEFT = sp.Rational(-1, 2)
X_RIGHT = sp.Rational(3, 2)
SIGMA_INNER = sp.Rational(1, 10)
SIGMA_OUTER = sp.Rational(1, 50)
ALPHA = sp.Rational(1, 100)
BOX = ((-1, 2), (0, 1), (0, 1))
PHI_REF_POINT = tuple(sp.Rational(s) for s in ("-0.9013", "0.0437", "0.0437"))
PHI_REF_VALUE = sp.Integer(0)

x, t, k, alpha = sp.symbols("x t k alpha", real=True)
LAYER_NAMES = ("xLeft", "xRight", "sigmaOuterLeft", "sigmaInnerLeft", "sigmaHeart",
               "sigmaInnerRight", "sigmaOuterRight")
L = {n: sp.Symbol(n, real=True) for n in LAYER_NAMES}
LAYER_ARGS = tuple(L[n] for n in LAYER_NAMES)
S = sp.sqrt(1 + t)


def vm():
    """V_m = S cos(pi x) in the heart."""
    return S*sp.cos(sp.pi*x)


def slab_constant():
    """C = k S, the constant of the FDA bath slab with zero bath current."""
    return k*S


def slab_branches(C=None):
    """[(phi*_slab, upper end of its interval)] from left to right; C defaults
    to slab_constant()."""
    C = slab_constant() if C is None else C
    return [(-k*S + C, 0), (-k*S*sp.cos(sp.pi*x) + C, 1), (k*S + C, None)]


def resistance_branches():
    """[(R, upper end of its interval)] from left to right."""
    xl, xr = L["xLeft"], L["xRight"]
    heart = 1/L["sigmaHeart"]
    return [(xl/L["sigmaInnerLeft"] + (x - xl)/L["sigmaOuterLeft"], xl),
            (x/L["sigmaInnerLeft"], 0),
            (x/L["sigmaHeart"], 1),
            (heart + (x - 1)/L["sigmaInnerRight"], xr),
            (heart + (xr - 1)/L["sigmaInnerRight"] + (x - xr)/L["sigmaOuterRight"], None)]


def layer_conductivities():
    """sigma of each branch of resistance_branches()."""
    return [L["sigmaOuterLeft"], L["sigmaInnerLeft"], L["sigmaHeart"], L["sigmaInnerRight"],
            L["sigmaOuterRight"]]


def _piecewise(branches):
    return sp.Piecewise(*[(e, True if b is None else x <= b) for e, b in branches])


@functools.lru_cache(maxsize=None)
def fields():
    """{Vm, phiSlab, resistance, phiLayered} as sympy expressions."""
    phi_slab = _piecewise(slab_branches())
    R = _piecewise(resistance_branches())
    return {"Vm": vm(), "phiSlab": phi_slab, "resistance": R,
            "phiLayered": phi_slab + alpha*R}


def case_layers():
    """{layer symbol: value} of the layered problem of the tutorial."""
    pair = MONODOMAIN_PAIR
    return {L["xLeft"]: X_LEFT, L["xRight"]: X_RIGHT,
            L["sigmaOuterLeft"]: SIGMA_OUTER, L["sigmaInnerLeft"]: SIGMA_INNER,
            L["sigmaHeart"]: pair[0] + pair[1],
            L["sigmaInnerRight"]: SIGMA_INNER, L["sigmaOuterRight"]: SIGMA_OUTER}


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def _zero(expr):
    return sp.simplify(expr) == 0


def run_self_checks():
    gi, ge = sp.symbols("gi ge", positive=True)
    kk = gi/(gi + ge)
    slab = [e.subs(k, kk) for e, _ in slab_branches()]
    _check(_zero((gi + ge)*sp.diff(slab[1], x, 2) + gi*sp.diff(vm(), x, 2)),
           "heart equation")
    for wall in (0, 1):
        _check(_zero(sp.diff(vm(), x).subs(x, wall)), "V_m' on x=%d" % wall)
    for (a, end), (b, _) in zip(slab_branches()[:-1], slab_branches()[1:]):
        _check(_zero((a - b).subs(x, end)), "phi* continuity at x=%s" % end)
    for side in (0, 2):
        _check(_zero(sp.diff(slab_branches()[side][0], x, 2)), "torso slab not linear")
    branches = resistance_branches()
    for (a, end), (b, _) in zip(branches[:-1], branches[1:]):
        _check(_zero((a - b).subs(x, end)), "R continuity at x=%s" % end)
    for (R, _), sigma in zip(branches, layer_conductivities()):
        _check(_zero(sigma*sp.diff(R, x) - 1), "sigma R' = 1")
    _check(_zero(branches[1][0].subs(x, 0)), "R(0) = 0")
    # outward normals -e_x on x = -1 and +e_x on x = 2
    left = -L["sigmaOuterLeft"]*sp.diff(alpha*branches[0][0], x)
    right = L["sigmaOuterRight"]*sp.diff(alpha*branches[-1][0], x)
    _check(_zero(left + alpha) and _zero(right - alpha), "surface currents")
    left, right = slab_branches()[0][0], slab_branches()[2][0]
    _check(_zero(left.subs(x, PHI_REF_POINT[0]) - left)
           and _zero(right - left - 2*k*S), "gauged slab potential")
    f = fields()
    _check(_zero(f["phiLayered"].subs(alpha, 0) - f["phiSlab"]), "alpha = 0")
    for name, (a, b) in (("monodomain", MONODOMAIN_PAIR), ("bidomain", BIDOMAIN_PAIR)):
        value = float(a/(a + b))
        _check(abs(value - 2**-0.5) <= K_RTOL*2**-0.5,
               "k of the %s pair differs from 1/sqrt(2)" % name)
    _check(BOX[0][0] < PHI_REF_POINT[0] < X_LEFT, "reference point in the outer left layer")
    for v in (SIGMA_SLAB, SIGMA_INNER, SIGMA_OUTER, ALPHA) + PHI_REF_POINT:
        printed_constant(v, 10)


# ------------------------------------------------------------------ LaTeX ----

def _names():
    sigma1, sigma2, sh, xl, xr = (sp.Symbol(n, real=True)
                                  for n in ("sigma1", "sigma2", "sH", "xL", "xR"))
    return ({k: "k", alpha: r"\alpha", sigma1: r"\sigma_1", sigma2: r"\sigma_2",
             sh: r"(G_i+G_e)_{xx}", xl: "x_L", xr: "x_R"},
            {L["sigmaOuterLeft"]: sigma2, L["sigmaInnerLeft"]: sigma1,
             L["sigmaHeart"]: sh, L["sigmaInnerRight"]: sigma1,
             L["sigmaOuterRight"]: sigma2, L["xLeft"]: xl, L["xRight"]: xr})


def _intervals(branches, labels):
    """Condition LaTeX of each branch from the upper ends of the intervals."""
    out, lower = [], None
    for (_, upper), label in zip(branches, labels):
        if lower is None:
            out.append(r"x\le %s" % label[1])
        elif upper is None:
            out.append(r"x>%s" % label[0])
        else:
            out.append(r"%s<x\le %s" % label)
        lower = upper
    return out


def latex_fragments():
    names, symmetric = _names()
    Ss = sp.Symbol("S", positive=True)
    names[Ss] = "S"
    sub = {S: Ss}
    gxx, exx = sp.symbols("gxx exx", real=True)
    names.update({gxx: r"G_{i,xx}", exx: r"G_{e,xx}"})
    out = {
        "leadfield_S.tex": emit_latex_fragment("S", S, names),
        "leadfield_k.tex": emit_latex_fragment("k", gxx/(gxx + exx), names),
        "leadfield_slab_Vm.tex": emit_latex_fragment("V_m", vm().subs(sub), names),
        "leadfield_slab_C.tex": emit_latex_fragment("C", slab_constant().subs(sub), names),
    }
    Cs = sp.Symbol("C", real=True)
    names[Cs] = "C"
    slab = slab_branches(Cs)
    conditions = _intervals(slab, [(None, "0"), ("0", "1"), ("1", None)])
    out["leadfield_slab_phi.tex"] = emit_cases_fragment(
        r"\phi^\ast_{\mathrm{slab}}",
        [(e.subs(sub), c) for (e, _), c in zip(slab, conditions)], names)
    # gauged in the left torso; the middle branch printed as k S (1 - cos(pi x))
    gauged = [sp.Integer(0), k*Ss*(1 - sp.cos(sp.pi*x)), 2*k*Ss]
    for g, (e, _) in zip(gauged, slab):
        _check(_zero(g - (e - slab[0][0]).subs(sub)), "printed gauged slab potential")
    out["leadfield_slab_phi_gauged.tex"] = emit_cases_fragment(
        r"\phi_{\mathrm{slab}}", list(zip(gauged, conditions)), names)
    res = resistance_branches()
    out["leadfield_layered_R.tex"] = emit_cases_fragment(
        "R", [(e.subs(symmetric), c) for (e, _), c in zip(res, _intervals(
            res, [(None, "x_L"), ("x_L", "0"), ("0", "1"), ("1", "x_R"), ("x_R", None)]))],
        names)
    phis, R = sp.symbols("phis R", real=True)
    names.update({phis: r"\phi^\ast_{\mathrm{slab}}", R: "R"})
    out["leadfield_layered_phi.tex"] = emit_latex_fragment(
        r"\phi^\ast_{\mathrm{layered}}", phis + alpha*R, names)
    out["leadfield_slab_params.tex"] = emit_parameter_fragment(
        [(r"\sigma", printed_constant(SIGMA_SLAB, 9)),
         (r"\mathbf r_{\mathrm{ref}}",
          sp.Matrix([printed_constant(v, 4) for v in PHI_REF_POINT])),
         (r"\phi_{\mathrm{ref}}", int(PHI_REF_VALUE))], per_row=3)
    out["leadfield_layered_params.tex"] = emit_parameter_fragment(
        [("x_L", printed_constant(X_LEFT, 1)), ("x_R", printed_constant(X_RIGHT, 2)),
         (r"\sigma_1", printed_constant(SIGMA_INNER, 1)),
         (r"\sigma_2", printed_constant(SIGMA_OUTER, 1)),
         (r"\alpha", printed_constant(ALPHA, 1))], per_row=5)
    return {"latex/" + k_: v for k_, v in out.items()}


# ------------------------------------------------------------------ NumPy ----

def numpy_module():
    f = fields()
    functions = {"Vm": ((x, t), f["Vm"]),
                 "phi_slab": ((x, t, k), f["phiSlab"]),
                 "resistance": ((x,) + LAYER_ARGS, f["resistance"]),
                 "phi_layered": ((x, t, k, alpha) + LAYER_ARGS, f["phiLayered"])}
    doc = ("Lead-field ECG slab and layered manufactured potentials, functions of x and t: "
           "Vm(x, t) in the heart; phi_slab(x, t, k) and phi_layered(x, t, k, alpha, "
           "*LAYER_NAMES) before the gauge (phi(x) - phi(x_ref) + PHI_REF_VALUE); "
           "resistance(x, *LAYER_NAMES), R with R' = 1/sigma and R(0) = 0.  PARAMS holds "
           "the values of the tutorial.  Generated by slab_leadfield_mms.py.")
    layers = {str(s): float(v) for s, v in case_layers().items()}
    params = {"K": float(MONODOMAIN_PAIR[0]/(MONODOMAIN_PAIR[0] + MONODOMAIN_PAIR[1])),
              "ALPHA": float(ALPHA), "SIGMA_SLAB": float(SIGMA_SLAB), "LAYERS": layers,
              "LAYER_NAMES": list(LAYER_NAMES),
              "PHI_REF_POINT": [float(v) for v in PHI_REF_POINT],
              "PHI_REF_VALUE": float(PHI_REF_VALUE)}
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
    return codegen.generator_main(argv, __doc__, __file__, run_self_checks, write_outputs, {})


if __name__ == "__main__":
    sys.exit(main())
