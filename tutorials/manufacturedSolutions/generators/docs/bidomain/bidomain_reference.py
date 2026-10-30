#!/usr/bin/env python3
"""Bidomain manufactured solutions: exact fields, LaTeX fragments and checks.

    python bidomain_reference.py --emit
    python bidomain_reference.py --self-check
    python bidomain_reference.py --run-check fda=CASE ... shell=CASE ...

FDA cube.  The field is that of the FDA monodomain problem
(reference/fda_monodomain.py: F_d, G_d, u1, u2, u3) with phiE = -kappa Vm + c
and phiI = Vm + phiE, kappa = Gi_xx/(Gi_xx + Ge_xx).  The NumPy functions below
follow src/verificationModels/bidomainVerification/
manufacturedFDABidomainReference.H and the constants of
src/ionicModels/verificationModels/bidomainFDAManufactured: beta_d is the sum
over the first d axes of -Gm_mm (m pi)^2 with Gm = Gi Ge/(Gi + Ge) entry by
entry, and c makes phiE equal phiEReferenceValue at phiERefPoint.

Sealed shell.  The fields, tensors and sources are those of
generators/shell_bidomain_mms.py (imported); this module adds its LaTeX
fragments, the checks of the bidomain equations on the exact fields and the
comparison with a run.

--emit runs the self-checks and writes latex/*.tex next to this file only if
they pass.  --run-check reads the cell centres (computed from
constant/polyMesh as OpenFOAM does), the fields of the last time directory
and postProcessing/*.dat of finished runs and checks that the exact fields of
this module reproduce the errors the cardiacFoam verifier reported; it
writes nothing.
"""

import argparse
import glob
import importlib.util
import math
import os
import re
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
GENERATORS = os.path.abspath(os.path.join(HERE, "..", ".."))
TUTORIAL = os.path.abspath(os.path.join(GENERATORS, "..", "bidomain"))
sys.path.insert(0, GENERATORS)

import shell_bidomain_mms as sbm  # noqa: E402
import shell_monodomain_mms as smm  # noqa: E402
from mms_common.latex_emit import (emit_latex_fragment,  # noqa: E402
                                   emit_parameter_fragment)


def _load(name):
    path = os.path.join(GENERATORS, "reference", name + ".py")
    spec = importlib.util.spec_from_file_location("ref_" + name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FDA = _load("fda_monodomain")

# ------------------------------------------------------------------ FDA cube ----

# Diagonals of the tutorial dictionary (constant/electroProperties of bidomain/)
GI = (0.111453302, 0.121585420, 0.030396355)
GE = (0.046165469, 0.050362330, 0.012590582)
CHI = 3.0
CM = 2.0
PHIE_REF_POINT = (0.125, 0.05, 0.05)
PHIE_REF_VALUE = 0.0


def kappa():
    """Gi_xx/(Gi_xx + Ge_xx), the k of bidomainFDABenchmarkConstants."""
    return GI[0]/(GI[0] + GE[0])


def beta(d):
    """Reaction coefficient of bidomainFDAManufactured in dimension d."""
    return -sum(GI[m]*GE[m]/(GI[m] + GE[m])*((m + 1)*math.pi)**2 for m in range(d))


def F(d, x, y, z):
    return getattr(FDA, "F_d%d" % d)(x, y, z, 0.0)


def gradF(d, x, y, z):
    x, y, z = np.broadcast_arrays(*(np.asarray(v, float) for v in (x, y, z)))
    pi = math.pi
    c = [np.cos(pi*x), np.cos(2*pi*y), np.cos(3*pi*z)]
    s = [-pi*np.sin(pi*x), -2*pi*np.sin(2*pi*y), -3*pi*np.sin(3*pi*z)]
    g = []
    for a in range(3):
        if a >= d:
            g.append(np.zeros_like(x))
            continue
        term = s[a]
        for b in range(d):
            if b != a:
                term = term*c[b]
        g.append(term)
    return np.stack(g, -1)


def Vm(d, x, y, z, t):
    return np.sqrt(1 + t)*F(d, x, y, z)


def shift(d, t, point=PHIE_REF_POINT, value=PHIE_REF_VALUE):
    """manufacturedFDABidomainShift: the constant c of phiE = -kappa Vm + c."""
    return value + kappa()*math.sqrt(1 + t)*float(F(d, *point))


def phiE(d, x, y, z, t):
    return -kappa()*np.sqrt(1 + t)*F(d, x, y, z) + shift(d, t)


def phiI(d, x, y, z, t):
    return Vm(d, x, y, z, t) + phiE(d, x, y, z, t)


def u1(d, x, y, z, t):
    return getattr(FDA, "u1_d%d" % d)(x, y, z, t)


def u2(d, x, y, z, t):
    return getattr(FDA, "u2_d%d" % d)(x, y, z, t)


def Jion(d, V, a1, a2, a3):
    """ALGEBRAIC[BidomainIion] of bidomainFDAManufactured_2014.H."""
    return -0.5*CM*(a1 + a3 - V)*a2**2*(V - a3) + beta(d)*(V - a3)/CHI


# Symbolic form, for the checks and the printed fragments
x_, y_, z_, t_ = sp.symbols("x y z t", real=True)
XYZ_ = (x_, y_, z_)


def _sym_constants():
    gi = [sp.Rational(repr(v)) for v in GI]
    ge = [sp.Rational(repr(v)) for v in GE]
    return gi, ge


def _sym_F(d):
    return sp.Mul(*[sp.cos((m + 1)*sp.pi*XYZ_[m]) for m in range(d)])


def _sym_G(d):
    return 1 + sp.Mul(*[XYZ_[m]**(m + 1) for m in range(d)])


def fda_residuals(d):
    """Coefficients rho_p, rho_e with residual_p = rho_p V_m and
    residual_e = rho_e V_m on the exact fields (exact rational arithmetic on
    the dictionary decimals); the scales are the coefficients of
    div(Gi grad V_m) and of the time derivative term."""
    gi, ge = _sym_constants()
    kap = gi[0]/(gi[0] + ge[0])
    bet = -sum(gi[m]*ge[m]/(gi[m] + ge[m])*((m + 1)*sp.pi)**2 for m in range(d))
    Fd = _sym_F(d)
    Gd = _sym_G(d)
    S = sp.sqrt(1 + t_)
    V = S*Fd
    pe = -kap*V
    pi_ = V + pe
    a1 = Gd*(1 + t_) + V
    a2 = 1/(sp.sqrt(Gd)*(1 + t_))
    a3 = sp.Integer(0)

    def div(g, u):
        return sum(sp.diff(g[m]*sp.diff(u, XYZ_[m]), XYZ_[m]) for m in range(d))

    chi, cm = sp.Integer(3), sp.Integer(2)
    assert (float(chi), float(cm)) == (CHI, CM)
    J = -sp.Rational(1, 2)*cm*(a1 + a3 - V)*a2**2*(V - a3) + bet*(V - a3)/chi
    I = J/cm
    res_p = chi*cm*sp.diff(V, t_) - (div(gi, pi_) - chi*cm*I)
    gsum = [gi[m] + ge[m] for m in range(d)]
    res_e = div(gi, V) + div(gsum, pe)
    rho_p = sp.simplify(res_p/V)
    rho_e = sp.simplify(res_e/V)
    assert rho_p.free_symbols == set() and rho_e.free_symbols == set()
    rate1 = a2**2*(a1 + a3 - V)**2 + sp.Rational(1, 2)*(a1 + a3 - V)*a2**2*(V - a3)
    rate2 = -(a1 + a3 - V)*a2**3
    rates = (sp.simplify(sp.diff(a1, t_) - rate1), sp.simplify(sp.diff(a2, t_) - rate2))
    scale = sum(gi[m]*((m + 1)*sp.pi)**2 for m in range(d))
    walls = []
    for m in range(d):
        for face in (0, 1):
            for u, g in ((V, gi), (pe, ge), (pi_, gi)):
                walls.append(sp.simplify((g[m]*sp.diff(u, XYZ_[m])).subs(XYZ_[m], face)))
    return {"rho_p": rho_p, "rho_e": rho_e, "scale": scale, "rates": rates,
            "walls": walls, "kappa": kap, "beta": bet}


def read_dictionary_conductivities(path):
    text = open(path).read()
    out = {}
    for name in ("conductivityIntracellular", "conductivityExtracellular"):
        m = re.search(name + r"\s*\[[^\]]*\]\s*\(([^)]*)\)", text)
        out[name] = [float(v) for v in m.group(1).split()]
    m = re.search(r"phiERefPoint\s*\(([^)]*)\)", text)
    out["phiERefPoint"] = tuple(float(v) for v in m.group(1).split())
    return out


def fda_self_checks(verbose=True):
    report = {}
    for d in (1, 2, 3):
        r = fda_residuals(d)
        assert all(sp.simplify(v) == 0 for v in r["rates"]), "ionic rates d=%d" % d
        assert all(v == 0 for v in r["walls"]), "wall flux d=%d" % d
        rel_p = float(abs(r["rho_p"])/r["scale"])
        rel_e = float(abs(r["rho_e"])/r["scale"])
        assert rel_p <= 1e-8 and rel_e <= 1e-8, "residual d=%d" % d
        assert abs(float(r["kappa"]) - kappa()) < 1e-15
        assert abs(float(r["beta"]) - beta(d)) < 1e-12*abs(beta(d))
        report[d] = (rel_p, rel_e)
        # gauge
        assert abs(phiE(d, *PHIE_REF_POINT, 0.37) - PHIE_REF_VALUE) < 1e-15
    gi = np.array(GI)
    ge = np.array(GE)
    assert np.max(np.abs(gi/(gi + ge) - kappa())) <= 1e-6*kappa(), "Ge proportional to Gi"
    case = read_dictionary_conductivities(os.path.join(TUTORIAL, "constant", "electroProperties"))
    gi6, ge6 = case["conductivityIntracellular"], case["conductivityExtracellular"]
    assert gi6 == [GI[0], 0, 0, GI[1], 0, GI[2]], "Gi of the tutorial"
    assert ge6 == [GE[0], 0, 0, GE[1], 0, GE[2]], "Ge of the tutorial"
    assert case["phiERefPoint"] == PHIE_REF_POINT, "phiERefPoint of the tutorial"
    if verbose:
        for d, (p, e) in report.items():
            print("FDA d=%d: |rho_p|/k = %.2e, |rho_e|/k = %.2e" % (d, p, e))
    return report


# ------------------------------------------------------------- sealed shell ----

class Shell:
    """Adapter with the interface of the emitted shell reference modules
    (PARAMS, field(name, x, y, z)) used by docs/_common/figstyle.py."""

    PARAMS = {"default": sbm.default_params()}

    @staticmethod
    def evaluate(name, x, y, z):
        pts = np.stack(np.broadcast_arrays(*(np.asarray(v, float) for v in (x, y, z))), -1)
        return sbm.evaluate(name, pts)

    def field(self, name, x, y, z):
        if name == "sigma":
            return self.tensor("Gi", x, y, z)
        if name == "Vm":
            return self.evaluate("Fi", x, y, z) - self.evaluate("Fe", x, y, z)
        return self.evaluate(name, x, y, z)

    def tensor(self, name, x, y, z):
        """(3, 3, N) array of Gi or Ge."""
        T = sbm._full(self.evaluate(name, x, y, z))
        return np.moveaxis(T.reshape(-1, 3, 3), 0, -1)


SHELL = Shell()


def shell_volume_mean(name, nr=40, nc=40, nphi=80):
    """Volume mean over the shell of a scalar of shell_bidomain_mms."""
    p = sbm.default_params()
    xr, wr = np.polynomial.legendre.leggauss(nr)
    xc, wc = np.polynomial.legendre.leggauss(nc)
    r = 0.5*(p["r2"] - p["r1"])*xr + 0.5*(p["r2"] + p["r1"])
    wr = 0.5*(p["r2"] - p["r1"])*wr*r**2
    ph = 2*math.pi*(np.arange(nphi) + 0.5)/nphi
    R, C, PH = np.meshgrid(r, xc, ph, indexing="ij")
    S = np.sqrt(1 - C**2)
    pts = np.stack([R*S*np.cos(PH), R*S*np.sin(PH), R*C], -1)
    W = np.broadcast_to(wr[:, None, None]*wc[None, :, None]*(2*math.pi/nphi), R.shape)
    volume = 4.0*math.pi*(p["r2"]**3 - p["r1"]**3)/3.0
    assert abs(W.sum() - volume) < 1e-12*volume, "shell quadrature weights"
    return float((W*sbm.evaluate(name, pts)).sum()/W.sum())


def shell_frame_check(n=500):
    """Gi and Ge share the eigenvectors (f, s, n_s); their eigenvalue ratios
    differ, so Gi + Ge is not a multiple of Gi.  Returns the largest commutator
    and the smallest distance of Ge from the span of Gi."""
    pts = sbm._shell_points(n, 0.0, seed=7)
    Gi = sbm._full(sbm.evaluate("Gi", pts))
    Ge = sbm._full(sbm.evaluate("Ge", pts))
    comm = np.abs(Gi @ Ge - Ge @ Gi).max()
    c = np.einsum("pij,pij->p", Ge, Gi)/np.einsum("pij,pij->p", Gi, Gi)
    off = np.linalg.norm(Ge - c[:, None, None]*Gi, axis=(1, 2))/np.linalg.norm(Ge, axis=(1, 2))
    mono = np.abs(Gi - sbm._full(smm.evaluate("G", pts))).max()
    return float(comm), float(off.min()), float(mono)


def shell_vm_equation_residual():
    """Residual of the Vm equation on the exact fields with the ionic closure
    of monodomainFDAManufactured and the source I_s of the verifier, written
    with the spatial factors as symbols."""
    t = sp.Symbol("t", positive=True)
    Fi, Fe, g, b, dGi = sp.symbols("Fi Fe g beta dGi", real=True)
    chi, cm = sp.symbols("chi C_m", positive=True)
    S = sp.sqrt(1 + t)
    V = S*(Fi - Fe)
    a1 = (1 + t)*g + V
    a2 = 1/((1 + t)*sp.sqrt(g))
    J = -sp.Rational(1, 2)*cm*a2**2*V*(a1 - V) + b*V/chi
    Is = S*(b*(Fi - Fe) - dGi)
    res = chi*cm*sp.diff(V, t) - (S*dGi - chi*cm*J/cm + Is)
    rate1 = a2**2*(a1 - V)**2 + sp.Rational(1, 2)*(a1 - V)*a2**2*V
    rate2 = -(a1 - V)*a2**3
    return [sp.simplify(res), sp.simplify(sp.diff(a1, t) - rate1),
            sp.simplify(sp.diff(a2, t) - rate2)]


def shell_self_checks(verbose=True):
    sbm.run_self_checks()
    comm, off, mono = shell_frame_check()
    assert comm < 1e-14, "Gi and Ge share the frame"
    assert off > 0.1, "Ge is not proportional to Gi"
    assert mono < 1e-14, "Gi equals sigma of the monodomain shell"
    assert all(v == 0 for v in shell_vm_equation_residual()), "Vm equation"
    case = read_shell_dictionary()
    p = sbm.default_params()
    rref = float(np.linalg.norm(case["phiERefPoint"]))
    assert p["r1"] < rref < p["r2"], "phiERefPoint inside the shell"
    if verbose:
        mean, absolute = sbm.phiE_source_mean()
        print("shell: [Gi, Ge] max %.1e; min |Ge - c Gi|/|Ge| = %.3f; "
              "phiE source mean %.1e of %.3f" % (comm, off, mean, absolute))


def read_shell_dictionary():
    text = open(os.path.join(TUTORIAL, "sphericalShell", "constant", "electroProperties")).read()
    m = re.search(r"phiERefPoint\s*\(([^)]*)\)", text)
    return {"phiERefPoint": tuple(float(v) for v in m.group(1).split()),
            "ionicModel": re.search(r"ionicModel\s+(\w+)", text).group(1)}


# ------------------------------------------------------------ LaTeX output ----

def _fda_fragments():
    gi, ge = _sym_constants()
    out = {}
    Gi = sp.diag(*[sp.Float(v, 9) for v in GI])
    Ge = sp.diag(*[sp.Float(v, 8) for v in GE])
    out["fda_bidomain_Gi.tex"] = emit_latex_fragment(r"\mathbf G_i", Gi, {})
    out["fda_bidomain_Ge.tex"] = emit_latex_fragment(r"\mathbf G_e", Ge, {})
    gxx, exx, V, c, Fd, t = sp.symbols("gxx exx V c F_d t", real=True)
    names = {gxx: r"G_{i,xx}", exx: r"G_{e,xx}", V: "V_m", c: "c"}
    kap = sp.Symbol("kappa", real=True)
    names[kap] = r"\kappa"
    out["fda_bidomain_kappa.tex"] = emit_latex_fragment(
        r"\kappa", gxx/(gxx + exx), names)
    out["fda_bidomain_phiE.tex"] = emit_latex_fragment(r"\phi_e", -kap*V + c, names)
    out["fda_bidomain_phiI.tex"] = emit_latex_fragment(r"\phi_i", V*(1 - kap) + c, names)
    xr, ve = sp.symbols("xr ve", real=True)
    names.update({xr: r"V_m(\mathbf x_{\mathrm{ref}},t)", ve: r"\phi_{e,\mathrm{ref}}"})
    out["fda_bidomain_c.tex"] = emit_latex_fragment("c", ve + kap*xr, names)
    out["fda_bidomain_beta.tex"] = (
        "$$\n\\beta_d = -\\sum_{m=1}^{d} H_{mm}\\,(m\\pi)^2,\\qquad "
        "H_{mm} = \\frac{G_{i,mm}\\,G_{e,mm}}{G_{i,mm}+G_{e,mm}}\n$$\n")
    out["fda_bidomain_kd.tex"] = (
        "$$\nk_d = \\sum_{m=1}^{d} G_{i,mm}\\,(m\\pi)^2\n$$\n")
    pairs = [(r"C_m", 2), (r"\chi", 3), (r"\kappa", sp.Float(kappa(), 12))]
    pairs += [(r"\beta_%d" % d, sp.Float(beta(d), 10)) for d in (1, 2, 3)]
    out["fda_bidomain_params.tex"] = emit_parameter_fragment(pairs, per_row=3)
    pairs = [(r"\mathbf x_{\mathrm{ref}}", sp.Matrix(PHIE_REF_POINT)),
             (r"\phi_{e,\mathrm{ref}}", 0)]
    out["fda_bidomain_ref.tex"] = emit_parameter_fragment(pairs)
    rows = []
    for d, (p, e) in fda_self_checks(verbose=False).items():
        rows.append(r"d=%d:\ \ \frac{|\rho_p|}{k_d} = %s,\quad \frac{|\rho_e|}{k_d} = %s"
                    % (d, _sci(p), _sci(e)))
    out["fda_bidomain_residual_relative.tex"] = (
        "$$\n\\begin{gathered}\n" + ",\\\\\n".join(rows) + "\n\\end{gathered}\n$$\n")
    return out


def _sci(v):
    if v == 0:
        return "0"
    m, e = ("%.1e" % v).split("e")
    return r"%s\times10^{%d}" % (m, int(e))


SHELL_TEX = dict(smm.SHELL_SYMBOLS)
SHELL_TEX.update({"gil": r"\sigma_{i,f}", "git": r"\sigma_{i,s}", "gin": r"\sigma_{i,n}",
                  "gel": r"\sigma_{e,f}", "get": r"\sigma_{e,s}", "gen": r"\sigma_{e,n}",
                  "fA0": "b_1", "fA1": "b_2", "fA2": "b_3", "fA3": "b_4"})


def _shell_names(extra=()):
    names = {sbm.P[n]: tex for n, tex in SHELL_TEX.items() if n in sbm.P}
    names[sp.Symbol("k", real=True)] = "k"
    names.update(dict(extra))
    return names


def _shell_fragments():
    out = {}
    ff, ss, nn = sp.symbols("ff ss nn", real=True)
    nm = _shell_names([(ff, r"\mathbf f\otimes\mathbf f"), (ss, r"\mathbf s\otimes\mathbf s"),
                       (nn, r"\mathbf n_s\otimes\mathbf n_s")])
    for tag, (l, tr, n) in (("Gi", ("gil", "git", "gin")), ("Ge", ("gel", "get", "gen"))):
        out["shell_bidomain_%s.tex" % tag] = emit_latex_fragment(
            r"\mathbf G_%s" % tag[1], sbm.P[l]*ff + sbm.P[tr]*ss + sbm.P[n]*nn, nm)
    names = _shell_names()
    out["shell_bidomain_hhat_e.tex"] = emit_latex_fragment(
        r"\hat h_e", sbm.expressions()["ehat"], names)
    out["shell_bidomain_hhat_i.tex"] = emit_latex_fragment(r"\hat h_i", smm._hhat(), names)
    for tag, tens in (("i", r"\mathbf G_i"), ("e", r"\mathbf G_e")):
        a, b, dj, lam, m1, m2, hh, om, m = sp.symbols("a b dj lam m1 m2 hh om m", real=True)
        nmj = _shell_names([
            (a, r"\left(\mathbf n_j\cdot %s\nabla\hat h_%s\right)" % (tens, tag)),
            (b, r"\left(\mathbf n_j\cdot %s\,\mathbf n_j\right)" % tens),
            (dj, "d_j"), (lam, r"\lambda"), (m1, "m_{%s,1}" % tag), (m2, "m_{%s,2}" % tag),
            (hh, r"\hat h_%s" % tag), (om, r"\omega"), (m, "m_%s" % tag)])
        out["shell_bidomain_mj_%s.tex" % tag] = emit_latex_fragment(
            "m_{%s,j}" % tag, smm._wall_correction(a, b, dj), nmj)
        out["shell_bidomain_m_%s.tex" % tag] = emit_latex_fragment(
            "m_%s" % tag, smm._blend(lam, m1, m2), nmj)
        out["shell_bidomain_F%s.tex" % tag] = emit_latex_fragment(
            "F_%s" % tag, hh + om*m, nmj)
    t, Fi, Fe, g, b, Vm, dGi, dGe, shift = sp.symbols(
        "t Fi Fe g beta Vm dGi dGe shift", real=True)
    nf = _shell_names([(Fi, "F_i"), (Fe, "F_e"), (g, "g"), (b, r"\beta"), (Vm, "V_m"),
                       (dGi, r"\nabla\!\cdot\!\left(\mathbf G_i\nabla F_i\right)"),
                       (dGe, r"\nabla\!\cdot\!\left(\mathbf G_e\nabla F_e\right)"),
                       (shift, "c")])
    S = sp.sqrt(1 + t)
    out["shell_bidomain_phiI.tex"] = emit_latex_fragment(r"\phi_i", S*Fi + shift, nf)
    xr, ve = sp.symbols("xr ve", real=True)
    nf.update({xr: r"F_e(\mathbf x_{\mathrm{ref}})", ve: r"\phi_{e,\mathrm{ref}}"})
    out["shell_bidomain_c.tex"] = emit_latex_fragment("c", ve - S*xr, nf)
    out["shell_bidomain_phiE.tex"] = emit_latex_fragment(r"\phi_e", S*Fe + shift, nf)
    out["shell_bidomain_Vm.tex"] = emit_latex_fragment("V_m", S*(Fi - Fe), nf)
    out["shell_bidomain_u1.tex"] = emit_latex_fragment("u_1", (1 + t)*g + Vm, nf)
    out["shell_bidomain_u2.tex"] = emit_latex_fragment("u_2", g**sp.Rational(-1, 2)/(1 + t), nf)
    out["shell_bidomain_Is.tex"] = emit_latex_fragment(
        r"I_{\mathrm{app}}", S*(b*(Fi - Fe) - dGi), nf)
    out["shell_bidomain_se.tex"] = emit_latex_fragment(r"s_e", S*(dGi + dGe), nf)
    p = sbm.default_params()
    pairs = []
    for name in sbm.SCALAR_PARAMS:
        v = p[name]
        if name in sbm.ANGLE_PARAMS:
            v = r"%g^\circ" % v
        pairs.append((SHELL_TEX[name], v))
    pairs += [("a_%d" % (i + 1), v) for i, v in enumerate(p["hA"])]
    pairs += [("b_%d" % (i + 1), v) for i, v in enumerate(p["fA"])]
    out["shell_bidomain_params.tex"] = emit_parameter_fragment(pairs, per_row=5)
    ref = read_shell_dictionary()["phiERefPoint"]
    out["shell_bidomain_ref.tex"] = emit_parameter_fragment(
        [(r"\mathbf x_{\mathrm{ref}}", sp.Matrix(ref)), (r"\phi_{e,\mathrm{ref}}", 0),
         (r"\beta", FDA.PARAMS["BETA"][3])])
    return out


def latex_fragments():
    out = _fda_fragments()
    out.update(_shell_fragments())
    return out


def emit():
    run_self_checks()
    outdir = os.path.join(HERE, "latex")
    os.makedirs(outdir, exist_ok=True)
    for name, text in sorted(latex_fragments().items()):
        with open(os.path.join(outdir, name), "w") as fh:
            fh.write(text)
        print("wrote latex/" + name)


def run_self_checks():
    fda_self_checks()
    shell_self_checks()
    print("self-checks passed")


# ------------------------------------------------------- comparison with runs ----

def _foam_lists(path):
    """Every list (count followed by '(') of an OpenFOAM file after its header,
    as raw bytes, in order; binary or ascii."""
    data = open(path, "rb").read()
    fmt = re.search(rb"format\s+(\w+);", data).group(1).decode()
    cls = re.search(rb"class\s+(\w+);", data).group(1).decode()
    body = data[data.index(b"}", data.index(b"FoamFile")) + 1:]
    return fmt, cls, body


_WIDTH = {"scalar": 1, "vector": 3, "symmTensor": 6, "tensor": 9}


def read_internal_field(path):
    fmt, cls, body = _foam_lists(path)
    kind = cls.replace("vol", "").replace("Field", "")
    kind = kind[0].lower() + kind[1:]
    m = re.search(rb"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(", body)
    if m is None:
        u = re.search(rb"internalField\s+uniform\s+([^;]+);", body)
        return float(u.group(1))
    n, w = int(m.group(1)), _WIDTH[kind]
    start = m.end()
    if fmt == "binary":
        a = np.frombuffer(body[start:start + 8*n*w], dtype="<f8")
    else:
        a = np.array(body[start:body.index(b")\n;", start)].replace(b"(", b" ").replace(b")", b" ")
                     .split(), float)
    return a.reshape(n, w) if w > 1 else a


def _read_list(path, dtype, width=1, which=0):
    fmt, cls, body = _foam_lists(path)
    pos, k = 0, 0
    while True:
        m = re.compile(rb"(\d+)\s*\(").search(body, pos)
        n = int(m.group(1))
        start = m.end()
        size = n*width*np.dtype(dtype).itemsize
        if fmt == "binary":
            a = np.frombuffer(body[start:start + size], dtype=dtype)
            end = start + size
        else:
            end = body.index(b"\n)", start)
            text = body[start:end].replace(b"(", b" ").replace(b")", b" ").split()
            a = np.array(text, dtype=dtype)
        if k == which:
            return a.reshape(n, width) if width > 1 else a
        pos, k = end + 1, k + 1


def cell_centres_and_volumes(case):
    """Cell centres and volumes of constant/polyMesh, by the face-triangle and
    face-pyramid decomposition of OpenFOAM's primitiveMeshTools."""
    pm = os.path.join(case, "constant", "polyMesh")
    p = _read_list(os.path.join(pm, "points"), "<f8", 3)
    fmt, cls, _ = _foam_lists(os.path.join(pm, "faces"))
    if cls == "faceCompactList":
        offs = _read_list(os.path.join(pm, "faces"), "<i4", 1, 0)
        labs = _read_list(os.path.join(pm, "faces"), "<i4", 1, 1)
        faces = [labs[offs[i]:offs[i + 1]] for i in range(len(offs) - 1)]
    else:
        raise NotImplementedError("faceList")
    own = _read_list(os.path.join(pm, "owner"), "<i4")
    nei = _read_list(os.path.join(pm, "neighbour"), "<i4")
    nf = len(faces)
    fc = np.zeros((nf, 3))
    fa = np.zeros((nf, 3))
    for i, f in enumerate(faces):
        q = p[f]
        if len(f) == 3:
            fc[i] = q.mean(0)
            fa[i] = 0.5*np.cross(q[1] - q[0], q[2] - q[0])
            continue
        c0 = q.mean(0)
        nxt = np.roll(q, -1, axis=0)
        n = np.cross(nxt - q, c0 - q)
        a = np.linalg.norm(n, axis=1)
        fc[i] = (1.0/3.0)*((a[:, None]*(q + nxt + c0)).sum(0))/a.sum()
        fa[i] = 0.5*n.sum(0)
    ncell = int(max(own.max(), nei.max() if len(nei) else 0)) + 1
    cest = np.zeros((ncell, 3))
    cnt = np.zeros(ncell)
    np.add.at(cest, own, fc)
    np.add.at(cnt, own, 1)
    np.add.at(cest, nei, fc[:len(nei)])
    np.add.at(cnt, nei, 1)
    cest /= cnt[:, None]
    cc = np.zeros((ncell, 3))
    vol = np.zeros(ncell)
    pyr = np.einsum("ij,ij->i", fa, fc - cest[own])
    np.add.at(cc, own, pyr[:, None]*(0.75*fc + 0.25*cest[own]))
    np.add.at(vol, own, pyr)
    ni = np.arange(len(nei))
    pyr = np.einsum("ij,ij->i", fa[ni], cest[nei] - fc[ni])
    np.add.at(cc, nei, pyr[:, None]*(0.75*fc[ni] + 0.25*cest[nei]))
    np.add.at(vol, nei, pyr)
    cc /= vol[:, None]
    return cc, vol/3.0


def _norms(e, V):
    e = np.abs(e)
    return (float((V*e).sum()/V.sum()), float(np.sqrt((V*e*e).sum()/V.sum())), float(e.max()))


def _rows(dat):
    out = {}
    for line in open(dat):
        parts = line.split()
        if len(parts) == 4 and not line.startswith("#"):
            try:
                out[parts[0]] = tuple(float(v) for v in parts[1:])
            except ValueError:
                pass
    return out


def _last_time(case):
    times = [d for d in os.listdir(case) if re.fullmatch(r"[0-9.e+-]+", d) and d != "0"]
    return max(times, key=float)


def _compare(label, mine, reported, rtol):
    worst = max(abs(a - b)/max(abs(b), 1e-300) for a, b in zip(mine, reported))
    ok = worst <= rtol
    print("  %-12s computed %s reported %s  rel. diff %.1e %s" % (
        label, " ".join("%.6g" % v for v in mine), " ".join("%.6g" % v for v in reported),
        worst, "ok" if ok else "FAIL"))
    return ok


def run_check_fda(case):
    """Vm, phiE (with the reference-point constant), phiE_gauge, phiI_gauge
    errors and the gauge offset of an FDA run, recomputed from its fields."""
    T = _last_time(case)
    t = float(T)
    dat = glob.glob(os.path.join(case, "postProcessing", "*_cells.dat"))[0]
    d = int(os.path.basename(dat)[0])
    rows = _rows(dat)
    offset = float(re.search(r"phiE gauge offset\s*=\s*(\S+)", open(dat).read()).group(1))
    C, V = cell_centres_and_volumes(case)
    x, y, z = C.T
    vm = read_internal_field(os.path.join(case, T, "Vm"))
    pe = read_internal_field(os.path.join(case, T, "phiE"))
    u1_0 = read_internal_field(os.path.join(case, T, "preProcess_u1"))
    ex_v, ex_e = Vm(d, x, y, z, t), phiE(d, x, y, z, t)
    off = float((V*(pe - ex_e)).sum()/V.sum())
    print("FDA %s (d=%d, t=%g, %d cells)" % (os.path.basename(case.rstrip("/")), d, t, len(V)))
    ok = _compare("Vm", _norms(vm - ex_v, V), rows["Vm"], 2e-5)
    ok &= _compare("phiE", _norms(pe - ex_e, V), rows["phiE"], 2e-5)
    ok &= _compare("gauge offset", (off,), (offset,), 2e-5)
    ok &= _compare("phiE_gauge", _norms(pe - ex_e - off, V), rows["phiE_gauge"], 2e-5)
    ok &= _compare("phiI_gauge", _norms(vm + pe - (ex_v + ex_e) - off, V), rows["phiI_gauge"], 2e-5)
    e0 = float(np.abs(u1_0 - u1(d, x, y, z, 0.0)).max())
    print("  preProcess_u1 against u1(t=0): max |diff| %.1e" % e0)
    return ok and e0 < 1e-12


def run_check_shell(case):
    """Gi, Ge written by setManufacturedFields, and the Vm, phiE_gauge errors
    of a shell run, recomputed from its fields."""
    T = _last_time(case)
    t = float(T)
    dat = glob.glob(os.path.join(case, "postProcessing", "*_cells.dat"))[0]
    rows = _rows(dat)
    C, V = cell_centres_and_volumes(case)
    print("shell %s (t=%g, %d cells)" % (os.path.basename(case.rstrip("/")), t, len(V)))
    ok = True
    for name, fname in (("Gi", "ConductivityIntracellular"), ("Ge", "ConductivityExtracellular")):
        G = read_internal_field(os.path.join(case, "0", fname))
        ref = sbm.evaluate(name, C)
        rel = float((np.linalg.norm(G - ref, axis=1)/np.linalg.norm(ref, axis=1)).max())
        print("  0/%s against %s at the cell centres: max rel. diff %.1e" % (fname, name, rel))
        ok &= rel < 1e-9
    S = math.sqrt(1 + t)
    Fi, Fe = sbm.evaluate("Fi", C), sbm.evaluate("Fe", C)
    ref = np.array(read_shell_dictionary()["phiERefPoint"])
    shift = -S*float(sbm.evaluate("Fe", ref))
    vm = read_internal_field(os.path.join(case, T, "Vm"))
    pe = read_internal_field(os.path.join(case, T, "phiE"))
    ex_e = S*Fe + shift
    off = float((V*(pe - ex_e)).sum()/V.sum())
    ok &= _compare("Vm", _norms(vm - S*(Fi - Fe), V), rows["Vm"], 1e-8)
    ok &= _compare("phiE_gauge", _norms(pe - ex_e - off, V), rows["phiE_gauge"], 1e-8)
    u1_0 = read_internal_field(os.path.join(case, T, "preProcess_u1"))
    g = sbm.evaluate("g", C)
    e0 = float(np.abs(u1_0 - (g + (Fi - Fe))).max())
    print("  preProcess_u1 against u1(t=0): max |diff| %.1e" % e0)
    return ok and e0 < 1e-12


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--emit", action="store_true")
    grp.add_argument("--self-check", action="store_true")
    grp.add_argument("--run-check", nargs="+", metavar="KIND=CASE",
                     help="KIND is fda or shell; CASE a finished run")
    args = ap.parse_args(argv)
    if args.emit:
        emit()
        return 0
    if args.self_check:
        run_self_checks()
        return 0
    ok = True
    for item in args.run_check:
        kind, case = item.split("=", 1)
        ok &= {"fda": run_check_fda, "shell": run_check_shell}[kind](case)
    print("run checks passed" if ok else "run checks FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
