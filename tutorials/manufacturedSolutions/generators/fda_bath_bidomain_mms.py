#!/usr/bin/env python3
"""FDA bath bidomain manufactured solution on the box [-1, 2] x [0, 1]^2.

The heart slab 0 <= x <= 1 lies between two bath slabs; V_m and phi_e depend
on x and t only, the auxiliary states also on y and z through G_d of
fda_monodomain_mms.py.  With S = sqrt(1 + t), F = cos(pi x), s_e = Ge_xx and
sigma_b = s_e/2 (t in units of 1 s, x in units of 1 m):

    V_m = S F - alpha x/s_e,
    u_1 = (1 + t) G_d + S F,  u_2 = G_d^(-1/2)/(1 + t),  u_3 = -alpha x/s_e,
    phi_e = -k S + alpha x/sigma_b + C                          x <= 0
          = -k S F + alpha x/s_e + C                            0 < x <= 1
          =  k S + alpha/s_e + alpha (x - 1)/sigma_b + C        x > 1
    C = k S + alpha/sigma_b,  phi_i = V_m + phi_e,

so that phi_e = 0 on the grounded face x = -1 and the bath carries the
current alpha out through x = 2.  The ionic model bathBidomainFDAManufactured
(hand-written) has the rates and the current of the FDA monodomain model with
beta = -1.1 (1 - k), C_m = 2 and chi = 3; k = 1/sqrt(2) and alpha = 1/100 are
the tutorial's.

Checks (--self-check): the rates of u_1, u_2, u_3 hold on the exact fields;
the V_m equation chi C_m dV_m/dt = div(Gi grad phi_i) - chi C_m I_ion and the
elliptic equation div(Gi grad V_m) + div((Gi + Ge) grad phi_e) = 0 hold in
the heart to RESIDUAL_RTOL of Gi_xx pi^2 at the tutorial's k (the nine-decimal
input makes Gi_xx pi^2 = 1.1 and Gi/(Gi + Ge) = k to rounding); phi_e is
linear in the bath, continuous at x = 0 and x = 1 with s_e phi_e' (heart) =
sigma_b phi_e' (bath) there, and n.Gi grad phi_i vanishes on both walls;
phi_e(-1) = 0 and sigma_b phi_e'(2) = alpha; sigma_b of the tutorial is
Ge_xx/2; with alpha = 0 phi_e is the slab potential of slab_leadfield_mms.py;
the printed constants round-trip.

--emit writes reference/fda_bath_bidomain.py (no LaTeX: no document prints
this solution).  --write PATH rewrites, and --check PATH compares, the
generated regions of manufacturedFDABathBidomainReference.H.

    python fda_bath_bidomain_mms.py --emit
    python fda_bath_bidomain_mms.py --self-check
    python fda_bath_bidomain_mms.py --write ../../../src/verificationModels/bathBidomainVerification/manufacturedFDABathBidomainReference.H
    python fda_bath_bidomain_mms.py --check ../../../src/verificationModels/bathBidomainVerification/manufacturedFDABathBidomainReference.H
"""

import functools
import os
import sys

import sympy as sp

import fda_bidomain_mms as fbm
import fda_monodomain_mms as fmm
import mms_codegen as codegen
from mms_common.latex_emit import printed_constant
from mms_common.numpy_emit import emit_numpy_module

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "fda_bath_bidomain"
DIMENSIONS = fmm.DIMENSIONS

# Tutorial manufacturedSolutions/bathBidomain: constant/electroProperties
# (conductivities, chi, cm, k, alpha), system/setTorsoOrganConductivityFieldDict
# (bath conductivity) and the ionic model bathBidomainFDAManufactured_2014.H
# (betaScale)
GI_DIAG = fbm.GI_DIAG
GE_DIAG = fbm.GE_DIAG
CM = 2
CHI = 3
K = 1/sp.sqrt(2)
ALPHA = sp.Rational(1, 100)
SIGMA_B = sp.Rational("0.0230827345")
BETA_SCALE = sp.Rational(-11, 10)
# Gi_xx pi^2 = 1.1 and Gi/(Gi + Ge) = k to the rounding of the nine-decimal input
RESIDUAL_RTOL = 1e-8

x, y, z, t = fmm.x, fmm.y, fmm.z, fmm.t
AXES = fmm.AXES
k, alpha = sp.symbols("k alpha", real=True)
se = sp.Symbol("se", positive=True)
# sigma_b and C as the C++ locals of the header
SB = sp.Symbol("sigmaB", positive=True)
C_S = sp.Symbol("C", real=True)
S = sp.sqrt(1 + t)


def sigma_b_of(s_e):
    """sigma_b = s_e/2."""
    return s_e/2


def beta_of(kk):
    """beta = -1.1 (1 - k), the constant of the ionic model."""
    return BETA_SCALE*(1 - kk)


def constant_of():
    """C = k S + alpha/sigma_b, so that phi_e(-1) = 0."""
    return k*S + alpha/SB


def spatial_factor():
    """F = cos(pi x)."""
    return sp.cos(sp.pi*x)


def vm_of(F):
    """V_m = S F - alpha x/s_e."""
    return S*F - alpha*x/se


def u_of(F, G):
    """(u_1, u_2, u_3)."""
    return (1 + t)*G + S*F, G**sp.Rational(-1, 2)/(1 + t), -alpha*x/se


def phie_branches(C=None):
    """[(phi_e, upper end of its interval)] from left to right, in x, t, k,
    alpha, s_e, sigma_b (symbol SB) and C (defaults to constant_of())."""
    C = constant_of() if C is None else C
    return [(-k*S + alpha*x/SB + C, 0),
            (-k*S*spatial_factor() + alpha*x/se + C, 1),
            (k*S + alpha/se + alpha*(x - 1)/SB + C, None)]


def phii_of(V, phiE):
    """phi_i = V_m + phi_e."""
    return V + phiE


def _physical(e):
    return e.subs(SB, sigma_b_of(se))


@functools.lru_cache(maxsize=None)
def _fields(d):
    F, G = spatial_factor(), fmm._G(d)
    V = vm_of(F)
    u1, u2, u3 = u_of(F, G)
    Iion = fmm.ionic_current(V, u1, u2, u3, sp.Integer(CM), beta_of(k), sp.Integer(CHI))
    phiE = _physical(sp.Piecewise(*[(e, True if b is None else x <= b)
                                    for e, b in phie_branches()]))
    return {"Vm": V, "phiE": phiE, "phiI": phii_of(V, phiE), "u1": u1, "u2": u2, "u3": u3,
            "Iion": Iion, "F": F, "G": G, "C": _physical(constant_of()),
            "sigmaB": sigma_b_of(se)}


def fields(d):
    """{Vm, phiE, phiI, u1, u2, u3, Iion, F, G, C, sigmaB} in x, y, z, t, k,
    alpha, se.  Iion is J_ion; the solver's current is Iion/C_m."""
    if d not in DIMENSIONS:
        raise ValueError("d must be 1, 2 or 3, got %r" % (d,))
    return dict(_fields(d))


def _div(diag, u, d):
    return sum(diag[m]*sp.diff(u, AXES[m], 2) for m in range(d))


def residuals(d):
    """(V_m equation, elliptic equation) in the heart, left minus right side."""
    f = fields(d)
    V = f["Vm"]
    phiE = _physical(phie_branches()[1][0])
    phiI = phii_of(V, phiE)
    total = [GI_DIAG[m] + GE_DIAG[m] for m in range(3)]
    res_p = CHI*CM*sp.diff(V, t) - (_div(GI_DIAG, phiI, d) - CHI*CM*(f["Iion"]/CM))
    res_e = _div(GI_DIAG, V, d) + _div(total, phiE, d)
    return res_p, res_e


def relative_residuals(d):
    """|residual|/(Gi_xx pi^2 S F) at the tutorial's k, s_e = Ge_xx; constants."""
    scale = GI_DIAG[0]*sp.pi**2
    out = []
    for r in residuals(d):
        rho = sp.simplify((r/(S*spatial_factor())).subs({k: K, se: GE_DIAG[0]}))
        if rho.free_symbols:
            raise AssertionError("d=%d: residual is not a multiple of S F" % d)
        out.append(float(abs(rho)/scale))
    return tuple(out)


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def _zero(e):
    return sp.simplify(e) == 0


def run_self_checks():
    for d in DIMENSIONS:
        f = fields(d)
        rates = fmm.ionic_rates(f["Vm"], f["u1"], f["u2"], f["u3"])
        for name in ("u1", "u2", "u3"):
            _check(_zero(sp.diff(f[name], t) - rates[name]),
                   "d=%d: d%s/dt differs from its ionic rate" % (d, name))
        for rel, label in zip(relative_residuals(d), ("V_m", "elliptic")):
            _check(rel <= RESIDUAL_RTOL, "d=%d: %s residual %.3e" % (d, label, rel))
    sb = sigma_b_of(se)
    left, heart, right = (_physical(e) for e, _ in phie_branches())
    phiI = phii_of(vm_of(spatial_factor()), heart)
    for wall, bath in ((0, left), (1, right)):
        _check(_zero((bath - heart).subs(x, wall)), "phi_e continuity at x=%d" % wall)
        _check(_zero((se*sp.diff(heart, x) - sb*sp.diff(bath, x)).subs(x, wall)),
               "current continuity at x=%d" % wall)
        _check(_zero(sp.diff(phiI, x).subs(x, wall)), "n.Gi grad phi_i at x=%d" % wall)
        _check(_zero(sp.diff(bath, x, 2)), "bath potential not linear")
        _check(_zero(sp.diff(bath, x) - alpha/sb), "bath gradient alpha/sigma_b")
    _check(_zero(left.subs(x, -1)), "phi_e(-1) = 0")
    _check(_zero(sb*sp.diff(right, x) - alpha), "surface current on x = 2")
    _check(SIGMA_B == GE_DIAG[0]/2, "sigma_b of the tutorial is not Ge_xx/2")
    import slab_leadfield_mms as slm
    for (a, _), (b, _) in zip(phie_branches(), slm.slab_branches()):
        _check(_zero(_physical(a).subs(alpha, 0) - b.subs(slm.k, k)),
               "alpha = 0 differs from the slab potential")
    for v in (SIGMA_B, ALPHA):
        printed_constant(v, 10)


# ------------------------------------------------------------------ NumPy ----

def numpy_module():
    args = (x, y, z, t, k, alpha, se)
    functions = {}
    for d in DIMENSIONS:
        for name in ("Vm", "phiE", "phiI", "u1", "u2", "u3", "Iion"):
            functions["%s_d%d" % (name, d)] = (args, fields(d)[name])
    doc = ("FDA bath bidomain manufactured fields on [-1, 2] x [0, 1]^2; "
           "<name>_d<d>(x, y, z, t, k, alpha, se) for Vm, phiE (all three slabs), "
           "phiI = Vm + phiE, u1, u2, u3, Iion (J_ion, the model's Iion), with "
           "sigma_b = se/2.  PARAMS holds the tutorial's values.  Generated by "
           "fda_bath_bidomain_mms.py.")
    params = {"K": float(K), "ALPHA": float(ALPHA), "SE": float(GE_DIAG[0]),
              "SIGMA_B": float(SIGMA_B), "CM": CM, "CHI": CHI,
              "BETA": float(beta_of(K))}
    return emit_numpy_module(REFERENCE_NAME, functions, params, doc)


def render_outputs():
    """{path relative to this folder: text} of every emitted file."""
    return {"reference/%s.py" % REFERENCE_NAME: numpy_module()}


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


# -------------------------------------------------------------------- C++ ----

PHIE_REGIONS = ("bathPhiE_leftBath", "bathPhiE_heart", "bathPhiE_rightBath")


def _sym(name):
    return sp.Symbol(name, real=True)


def bath_header_regions(text):
    """{region name: C++ lines} of manufacturedFDABathBidomainReference.H.
    sigma_b and C enter the potentials as the locals sigmaB and C."""
    ind = functools.partial(codegen.region_indent, text)
    xi, Fi, Gi = _sym("x[i]"), _sym("F[i]"), _sym("G[i]")
    out = {"manufacturedFDABathSigmaB": codegen.statement(
        "return", sigma_b_of(se), ind("manufacturedFDABathSigmaB"))}
    out["manufacturedFDABathC"] = codegen.statement(
        "return", constant_of(), ind("manufacturedFDABathC"))
    out["computeManufacturedFDABathF"] = codegen.assignment(
        "F[i]", spatial_factor().subs(x, xi), ind("computeManufacturedFDABathF"))
    out["computeManufacturedFDABathV"] = codegen.assignment(
        "V[i]", vm_of(Fi).subs(x, xi), ind("computeManufacturedFDABathV"))
    names = ("u1[i]", "u2[i]", "u3[i]")
    out["computeManufacturedFDABathU"] = sum(
        (codegen.assignment(n, u.subs(x, xi), ind("computeManufacturedFDABathU"))
         for n, u in zip(names, u_of(Fi, Gi))), [])
    for name, (e, _) in zip(PHIE_REGIONS, phie_branches(C_S)):
        out[name] = codegen.statement("return", e, ind(name))
    out["computeManufacturedFDABathPhiI"] = codegen.assignment(
        "phiI[i]", phii_of(_sym("V[i]"), _sym("phiE[i]")),
        ind("computeManufacturedFDABathPhiI"))
    V = vm_of(spatial_factor())
    out["manufacturedFDABathVmAt"] = codegen.statement(
        "return", V, ind("manufacturedFDABathVmAt"))
    out["manufacturedFDABathPhiIAt"] = codegen.statement(
        "return", phii_of(_sym("V"), _sym("phiE")), ind("manufacturedFDABathPhiIAt"))
    zero = sp.Integer(0)
    grads = (("manufacturedFDABathGradVm", sp.diff(V, x)),
             ("manufacturedFDABathGradPhiEHeart", sp.diff(phie_branches()[1][0], x)),
             ("manufacturedFDABathGradPhiIHeart", sp.diff(phii_of(V, phie_branches()[1][0]), x)),
             ("manufacturedFDABathGradPhiEBath", sp.diff(phie_branches()[0][0], x)))
    for name, gx in grads:
        out[name] = codegen.call("vector", [gx, zero, zero], ind(name))
    return out


HEADER_PATH = ("../../../src/verificationModels/bathBidomainVerification/"
               "manufacturedFDABathBidomainReference.H")
HEADER_PATHS = {"manufacturedFDABathBidomainReference.H": HEADER_PATH}
HEADERS = {"manufacturedFDABathBidomainReference.H": bath_header_regions}


def header_matches(path):
    return codegen.regions_match(__file__, path, HEADERS[os.path.basename(path)])


def main(argv=None):
    return codegen.generator_main(argv, __doc__, __file__, run_self_checks, write_outputs,
                                  HEADERS)


if __name__ == "__main__":
    sys.exit(main())
