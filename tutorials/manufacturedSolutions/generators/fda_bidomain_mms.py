#!/usr/bin/env python3
"""FDA bidomain manufactured solution on the unit cube, d = 1, 2, 3.

The transmembrane potential and the auxiliary states are those of the FDA
monodomain problem (fda_monodomain_mms.py), with t in units of 1 s and x in
units of 1 m:

    V_m = sqrt(1 + t) F_d,  u_1 = (1 + t) G_d + V_m,  u_2 = G_d^(-1/2)/(1 + t),
    u_3 = 0,
    phi_e = -kappa V_m + c(t),  phi_i = V_m + phi_e,
    kappa = Gi_xx/(Gi_xx + Ge_xx),
    c(t) = phi_e,ref + kappa V_m(x_ref, t),

so that phi_e(x_ref, t) = phi_e,ref.  Gi and Ge are diagonal, with Ge given to
the nine decimals of Gi and Gi/(Gi + Ge) equal to kappa on every axis to about
1e-9, so the elliptic equation div(Gi grad V_m) + div((Gi + Ge) grad phi_e) =
s_e needs no source.  The ionic model bidomainFDAManufactured has the current
of the FDA monodomain model with beta_d = -sum_{m<=d} H_mm (m pi)^2,
H_mm = Gi_mm Ge_mm/(Gi_mm + Ge_mm), which cancels the diffusion term of the
V_m equation chi C_m dV_m/dt = div(Gi grad phi_i) - chi C_m I_ion with
I_ion = J_ion/C_m, C_m = 2 and chi = 3.

Checks (--self-check), in exact rational arithmetic on the decimals of the
conductivities: the ionic rates of u_1 and u_2 hold on the exact fields; the
residuals of both equations are rho_p V_m and rho_e V_m with constant rho_p,
rho_e and |rho|/k_d at most RESIDUAL_RTOL, k_d = sum_{m<=d} Gi_mm (m pi)^2;
n.Gi grad V_m, n.Ge grad phi_e and n.Gi grad phi_i vanish on every face of
the cube; phi_e(x_ref, t) = phi_e,ref; the ratios Gi_mm/(Gi_mm + Ge_mm)
agree with kappa to KAPPA_RTOL, the tolerance of the ionic model; the
printed constants round-trip.

The cardiacFoam functions of this solution, manufacturedFDABidomainReference.H
and the constants of the ionic model bidomainFDAManufactured, carry no
generated regions, so this generator has no --write or --check.

    python fda_bidomain_mms.py --emit
    python fda_bidomain_mms.py --self-check
"""

import functools
import os
import sys

import sympy as sp

import fda_monodomain_mms as fmm
import mms_codegen as codegen
from mms_common.latex_emit import (emit_latex_fragment, emit_parameter_fragment,
                                   printed_constant)
from mms_common.numpy_emit import emit_numpy_module

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "fda_bidomain"
DIMENSIONS = fmm.DIMENSIONS

# Diagonals of the conductivities of the tutorial manufacturedSolutions/bidomain
GI_DIAG = tuple(sp.Rational(s) for s in ("0.111453302", "0.121585420", "0.030396355"))
GE_DIAG = tuple(sp.Rational(s) for s in ("0.046165469", "0.050362330", "0.012590582"))
CM = 2
CHI = 3
PHIE_REF_POINT = tuple(sp.Rational(s) for s in ("0.125", "0.05", "0.05"))
PHIE_REF_VALUE = sp.Integer(0)
# the ratios Gi_mm/(Gi_mm + Ge_mm) agree with kappa to about 1e-9 (input rounding)
RESIDUAL_RTOL = 1e-8
# tolerance of bidomainFDAManufactured on the ratios Gi_mm/(Gi_mm + Ge_mm)
KAPPA_RTOL = 1e-6

x, y, z, t = fmm.x, fmm.y, fmm.z, fmm.t
AXES = fmm.AXES


def kappa():
    """Gi_xx/(Gi_xx + Ge_xx), exact."""
    return GI_DIAG[0]/(GI_DIAG[0] + GE_DIAG[0])


def harmonic(m):
    """H_mm = Gi_mm Ge_mm/(Gi_mm + Ge_mm)."""
    return GI_DIAG[m]*GE_DIAG[m]/(GI_DIAG[m] + GE_DIAG[m])


def beta(d):
    """beta_d = -sum_{m<=d} H_mm (m pi)^2, the constant of the ionic model."""
    return -sum(harmonic(m)*((m + 1)*sp.pi)**2 for m in range(d))


def eigenvalue(d):
    """k_d with div(Gi grad F_d) = -k_d F_d."""
    return sum(GI_DIAG[m]*((m + 1)*sp.pi)**2 for m in range(d))


def _at_reference(expr):
    return expr.subs(dict(zip(AXES, PHIE_REF_POINT)), simultaneous=True)


def shift(d):
    """c(t) = phi_e,ref + kappa V_m(x_ref, t)."""
    return PHIE_REF_VALUE + kappa()*_at_reference(fmm.vm_of(fmm._F(d)))


@functools.lru_cache(maxsize=None)
def _fields(d):
    F, G = fmm._F(d), fmm._G(d)
    V = fmm.vm_of(F)
    phiE = -kappa()*V + shift(d)
    u1, u2, u3 = fmm.u_of(F, G)
    Iion = fmm.ionic_current(V, u1, u2, u3, sp.Integer(CM), beta(d), sp.Integer(CHI))
    return {"Vm": V, "phiE": phiE, "phiI": V + phiE, "shift": shift(d),
            "u1": u1, "u2": u2, "u3": u3, "Iion": Iion, "F": F, "G": G}


def fields(d):
    """{Vm, phiE, phiI, shift, u1, u2, u3, Iion, F, G} as sympy expressions in
    x, y, z, t.  Iion is J_ion; the solver's current is Iion/C_m."""
    if d not in DIMENSIONS:
        raise ValueError("d must be 1, 2 or 3, got %r" % (d,))
    return dict(_fields(d))


def _div(diag, u, d):
    return sum(diag[m]*sp.diff(u, AXES[m], 2) for m in range(d))


def residuals(d):
    """(residual of the V_m equation, residual of the phi_e equation), each
    the left minus the right side on the exact fields."""
    f = fields(d)
    V, phiE, phiI = f["Vm"], f["phiE"], f["phiI"]
    total = [GI_DIAG[m] + GE_DIAG[m] for m in range(3)]
    res_p = CHI*CM*sp.diff(V, t) - (_div(GI_DIAG, phiI, d) - CHI*CM*(f["Iion"]/CM))
    res_e = _div(GI_DIAG, V, d) + _div(total, phiE, d)
    return res_p, res_e


@functools.lru_cache(maxsize=None)
def residual_coefficients(d):
    """(rho_p, rho_e): the residuals divided by V_m, constants."""
    V = fields(d)["Vm"]
    return tuple(sp.simplify(r/V) for r in residuals(d))


def relative_residuals(d):
    """(|rho_p|/k_d, |rho_e|/k_d)."""
    k = eigenvalue(d)
    return tuple(float(abs(r)/k) for r in residual_coefficients(d))


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def run_self_checks():
    kap = kappa()
    for m in range(3):
        ratio = GI_DIAG[m]/(GI_DIAG[m] + GE_DIAG[m])
        _check(abs(float(ratio - kap)) <= KAPPA_RTOL*float(kap),
               "Gi/(Gi + Ge) on axis %d differs from kappa" % m)
    for d in DIMENSIONS:
        f = fields(d)
        V, u1, u2, u3 = f["Vm"], f["u1"], f["u2"], f["u3"]
        rates = fmm.ionic_rates(V, u1, u2, u3)
        for name, u in (("u1", u1), ("u2", u2), ("u3", u3)):
            _check(sp.simplify(sp.diff(u, t) - rates[name]) == 0,
                   "d=%d: d%s/dt differs from its ionic rate" % (d, name))
        for rho, label in zip(residual_coefficients(d), ("V_m", "phi_e")):
            _check(rho.free_symbols == set(),
                   "d=%d: the %s residual is not a constant multiple of V_m" % (d, label))
        for rel, label in zip(relative_residuals(d), ("V_m", "phi_e")):
            _check(rel <= RESIDUAL_RTOL,
                   "d=%d: |rho|/k_d = %.3e of the %s equation exceeds %.1e"
                   % (d, rel, label, RESIDUAL_RTOL))
        for m in range(d):
            for face in (0, 1):
                for u, diag, name in ((V, GI_DIAG, "Gi grad V_m"),
                                      (f["phiE"], GE_DIAG, "Ge grad phi_e"),
                                      (f["phiI"], GI_DIAG, "Gi grad phi_i")):
                    flux = (diag[m]*sp.diff(u, AXES[m])).subs(AXES[m], face)
                    _check(sp.simplify(flux) == 0,
                           "d=%d: n.%s does not vanish on %s=%d" % (d, name, AXES[m], face))
        _check(sp.simplify(_at_reference(f["phiE"]) - PHIE_REF_VALUE) == 0,
               "d=%d: phi_e(x_ref, t) differs from phi_e,ref" % d)
    for v in GI_DIAG + GE_DIAG + PHIE_REF_POINT:
        printed_constant(v, 9)


# ------------------------------------------------------------------ LaTeX ----

def _sci(v):
    if v == 0:
        return "0"
    mantissa, exponent = ("%.1e" % v).split("e")
    return r"%s\times10^{%d}" % (mantissa, int(exponent))


def latex_fragments():
    out = {}
    out["fda_bidomain_Gi.tex"] = emit_latex_fragment(
        r"\mathbf G_i", sp.diag(*[printed_constant(v, 9) for v in GI_DIAG]), {})
    out["fda_bidomain_Ge.tex"] = emit_latex_fragment(
        r"\mathbf G_e", sp.diag(*[printed_constant(v, 8) for v in GE_DIAG]), {})
    gxx, exx, V, c, kap, xr, ve = sp.symbols("gxx exx V c kappa xr ve", real=True)
    names = {gxx: r"G_{i,xx}", exx: r"G_{e,xx}", V: "V_m", c: "c", kap: r"\kappa",
             xr: r"V_m(\mathbf x_{\mathrm{ref}},t)", ve: r"\phi_{e,\mathrm{ref}}"}
    out["fda_bidomain_kappa.tex"] = emit_latex_fragment(r"\kappa", gxx/(gxx + exx), names)
    out["fda_bidomain_phiE.tex"] = emit_latex_fragment(r"\phi_e", -kap*V + c, names)
    out["fda_bidomain_phiI.tex"] = emit_latex_fragment(r"\phi_i", V*(1 - kap) + c, names)
    out["fda_bidomain_c.tex"] = emit_latex_fragment("c", ve + kap*xr, names)
    out["fda_bidomain_beta.tex"] = (
        "$$\n\\beta_d = -\\sum_{m=1}^{d} H_{mm}\\,(m\\pi)^2,\\qquad "
        "H_{mm} = \\frac{G_{i,mm}\\,G_{e,mm}}{G_{i,mm}+G_{e,mm}}\n$$\n")
    out["fda_bidomain_kd.tex"] = (
        "$$\nk_d = \\sum_{m=1}^{d} G_{i,mm}\\,(m\\pi)^2\n$$\n")
    pairs = [(r"C_m", CM), (r"\chi", CHI), (r"\kappa", sp.Float(float(kappa()), 12))]
    pairs += [(r"\beta_%d" % d, sp.Float(float(beta(d)), 10)) for d in DIMENSIONS]
    out["fda_bidomain_params.tex"] = emit_parameter_fragment(pairs, per_row=3)
    out["fda_bidomain_ref.tex"] = emit_parameter_fragment(
        [(r"\mathbf x_{\mathrm{ref}}",
          sp.Matrix([printed_constant(v, 9) for v in PHIE_REF_POINT])),
         (r"\phi_{e,\mathrm{ref}}", int(PHIE_REF_VALUE))])
    rows = [r"d=%d:\ \ \frac{|\rho_p|}{k_d} = %s,\quad \frac{|\rho_e|}{k_d} = %s"
            % ((d,) + tuple(_sci(v) for v in relative_residuals(d))) for d in DIMENSIONS]
    out["fda_bidomain_residual_relative.tex"] = (
        "$$\n\\begin{gathered}\n" + ",\\\\\n".join(rows) + "\n\\end{gathered}\n$$\n")
    return {"latex/" + k: v for k, v in out.items()}


# ------------------------------------------------------------------ NumPy ----

def numpy_module():
    functions = {}
    for d in DIMENSIONS:
        for name, expr in fields(d).items():
            functions["%s_d%d" % (name, d)] = ((x, y, z, t), expr)
    doc = ("FDA bidomain manufactured fields on the unit cube; <name>_d<d>(x, y, z, t) "
           "for Vm, phiE, phiI, shift (the constant c of phi_e = -kappa Vm + c), u1, u2, "
           "u3, Iion (J_ion, the model's Iion), F, G.  phiE is gauged to PHIE_REF_VALUE "
           "at PHIE_REF_POINT.  Generated by fda_bidomain_mms.py.")
    params = {"GI_DIAG": [float(v) for v in GI_DIAG], "GE_DIAG": [float(v) for v in GE_DIAG],
              "CM": CM, "CHI": CHI, "KAPPA": float(kappa()),
              "BETA": {d: float(beta(d)) for d in DIMENSIONS},
              "PHIE_REF_POINT": [float(v) for v in PHIE_REF_POINT],
              "PHIE_REF_VALUE": float(PHIE_REF_VALUE)}
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
