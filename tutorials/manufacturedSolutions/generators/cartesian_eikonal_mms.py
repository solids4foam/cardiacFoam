#!/usr/bin/env python3
"""Manufactured eikonal field on the unit cube [0, 1]^d.

The activation-time field is

    psi(x) = exp(k_d . x),   k_1 = (1, 0, 0), k_2 = (1, 1, 0), k_3 = (1, 1, 1),

so grad psi = k_d psi.  With the constant tensor M = sigma / (chi C_m) and
K_M = k_d . M k_d the eikonal equation

    -div(M grad psi) + c_0 sqrt(grad psi . M grad psi) - 1 = S_psi

has the source

    S_psi = psi (c_0 sqrt(K_M) - K_M) - 1.

The generator raises for K_M <= 1e-300, the bound of the C++ verifier (the
conductivity metric must be positive along k_d).

fields(d, sigma, chi, cm, c0) returns psi, grad_psi, S_psi, K_M and M for a
3x3 sigma.  The emitted NumPy module holds tau_d<d>, grad_tau_d<d> and
source_d<d> (arguments x, y, z, c0, chi, cm, sxx, sxy, sxz, syy, syz, szz)
and a checked wrapper source(d, point, sigma, c0, chi, cm).

Checks (--self-check): the residual of the eikonal equation vanishes for
d = 1, 2, 3 and c_0 in {1, 5, 66}; K_M <= 0 is rejected.

--write PATH rewrites, and --check PATH compares, the generated regions of the
hand-written manufacturedEikonalReference.H: k_d, M, psi and grad psi for a
runtime k, K_M = k . M k and S_psi.  The positivity guard on K_M stays
hand-written.

    python cartesian_eikonal_mms.py --emit
    python cartesian_eikonal_mms.py --self-check
    python cartesian_eikonal_mms.py --write ../../../src/verificationModels/eikonalVerification/manufacturedEikonalReference.H
    python cartesian_eikonal_mms.py --check ../../../src/verificationModels/eikonalVerification/manufacturedEikonalReference.H
"""

import os
import sys

import sympy as sp

import mms_codegen as codegen
from mms_common.latex_emit import emit_latex_fragment, paper_symbol_names
from mms_common.numpy_emit import emit_numpy_module

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "cartesian_eikonal"
DIMENSIONS = (1, 2, 3)
KM_MIN = 1e-300
CHI = 3
CM = 2

x, y, z = sp.symbols("x y z", real=True)
AXES = (x, y, z)
C0 = sp.Symbol("c0", positive=True)
CHI_SYMBOL = sp.Symbol("chi", positive=True)
CM_SYMBOL = sp.Symbol("cm", positive=True)
SIGMA_SYMBOLS = sp.symbols("sxx sxy sxz syy syz szz", real=True)

# Tutorial tensors (eikonalECG/constant/electroProperties and README).
SIGMA_ROTATED = (0.106875, -0.0084931, -0.022561, 0.116682, -0.0130256, 0.0398782)
SIGMA_AXIS = (0.111453302, 0.0, 0.0, 0.121585420, 0.0, 0.030396355)


def sigma_matrix(entries):
    """Symmetric 3x3 matrix from (xx, xy, xz, yy, yz, zz)."""
    a, b, c, d, e, f = entries
    return sp.Matrix([[a, b, c], [b, d, e], [c, e, f]])


SIGMA_SYMBOLIC = sigma_matrix(SIGMA_SYMBOLS)


def k_vector(d):
    """k_d as a 3-vector."""
    if d not in DIMENSIONS:
        raise ValueError("d must be 1, 2 or 3, got %r" % (d,))
    return sp.Matrix([1 if m < d else 0 for m in range(3)])


def psi_of_k(k):
    """psi = exp(k . x) for a 3-vector k."""
    return sp.exp((sp.Matrix(k).T*sp.Matrix(AXES))[0])


def psi_expr(d):
    return psi_of_k(k_vector(d))


def k_metric(k, M):
    """K_M = k . M k."""
    k = sp.Matrix(k)
    return (k.T*M*k)[0]


def source_of(psi, KM, c0):
    """S_psi = psi (c0 sqrt(K_M) - K_M) - 1."""
    return psi*(c0*sp.sqrt(KM) - KM) - 1


def metric(sigma, chi=CHI, cm=CM):
    """M = sigma / (chi C_m)."""
    return sp.Matrix(sigma)/(chi*cm)


def fields(d, sigma, chi=CHI, cm=CM, c0=1):
    """psi, grad_psi, S_psi, K_M, M of the manufactured eikonal field."""
    k = k_vector(d)
    M = metric(sigma, chi, cm)
    KM = sp.expand(k_metric(k, M))
    if KM.is_number and not float(KM) > KM_MIN:
        raise ValueError("manufactured eikonal field requires K_M > %g, got %s" % (KM_MIN, KM))
    psi = psi_expr(d)
    S = source_of(psi, KM, c0)
    return {"psi": psi, "grad_psi": k*psi, "S_psi": S, "K_M": KM, "M": M}


def residual(d, sigma, chi=CHI, cm=CM, c0=1):
    """-div(M grad psi) + c0 sqrt(grad psi . M grad psi) - 1 - S_psi."""
    f = fields(d, sigma, chi, cm, c0)
    M, g = f["M"], f["grad_psi"]
    flux = M*g
    div = sum(sp.diff(flux[m], AXES[m]) for m in range(3))
    return -div + c0*sp.sqrt((g.T*M*g)[0]) - 1 - f["S_psi"]


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def run_self_checks():
    pt = {x: 0.37, y: 0.58, z: 0.81}
    for d in DIMENSIONS:
        for sigma in (sigma_matrix(SIGMA_ROTATED), sigma_matrix(SIGMA_AXIS)):
            for c0 in (1, 5, 66):
                r = residual(d, sigma, c0=c0).subs(pt)
                _check(abs(float(r)) < 1e-12, "d=%d c0=%s: eikonal residual %s" % (d, c0, r))
    try:
        fields(1, sp.zeros(3, 3))
    except ValueError:
        pass
    else:
        raise AssertionError("K_M = 0 was accepted")


# Source appended to the emitted module.
RUNTIME = '''

KM_MIN = 1e-300


def _sigma_entries(sigma):
    s = numpy.asarray(sigma, dtype=float)
    if s.shape == (3, 3):
        return (s[0, 0], s[0, 1], s[0, 2], s[1, 1], s[1, 2], s[2, 2])
    return tuple(float(v) for v in s)


def k_metric(d, sigma, chi=3.0, cm=2.0):
    """K_M = k_d . M k_d with M = sigma/(chi cm); raises when not positive."""
    sxx, sxy, sxz, syy, syz, szz = _sigma_entries(sigma)
    if d == 1:
        km = sxx
    elif d == 2:
        km = sxx + 2.0*sxy + syy
    else:
        km = sxx + syy + szz + 2.0*(sxy + sxz + syz)
    km = km/(chi*cm)
    if not km > KM_MIN:
        raise ValueError("K_M must exceed %g, got %r" % (KM_MIN, km))
    return km


def source(d, point, sigma, c0, chi=3.0, cm=2.0):
    """S_psi at point (x, y, z); raises when K_M <= 1e-300."""
    k_metric(d, sigma, chi, cm)
    f = globals()["source_d%d" % d]
    return f(*point, c0, chi, cm, *_sigma_entries(sigma))
'''


def numpy_module():
    functions = {}
    args = (x, y, z)
    M = SIGMA_SYMBOLIC/(CHI_SYMBOL*CM_SYMBOL)
    for d in DIMENSIONS:
        psi = psi_expr(d)
        k = k_vector(d)
        KM = sp.expand(k_metric(k, M))
        functions["tau_d%d" % d] = (args, psi)
        functions["grad_tau_d%d" % d] = (args, k*psi)
        functions["source_d%d" % d] = (
            args + (C0, CHI_SYMBOL, CM_SYMBOL) + tuple(SIGMA_SYMBOLS),
            source_of(psi, KM, C0))
    doc = ("Manufactured eikonal field psi = exp(k_d . x) on the unit cube: tau_d<d>, "
           "grad_tau_d<d>(x, y, z) and source_d<d>(x, y, z, c0, chi, cm, sxx, sxy, sxz, "
           "syy, syz, szz); source(d, point, sigma, c0) checks K_M > 0.  Generated by "
           "cartesian_eikonal_mms.py.")
    params = {"SIGMA_ROTATED": list(SIGMA_ROTATED), "SIGMA_AXIS": list(SIGMA_AXIS),
              "CHI": float(CHI), "CM": float(CM), "KM_MIN": KM_MIN}
    return emit_numpy_module(REFERENCE_NAME, functions, params, doc).rstrip("\n") + "\n" + RUNTIME


def latex_fragments():
    psi_s, km_s, kd, msym = sp.symbols("psi K_M k_d M", real=True)
    names = paper_symbol_names([C0, CHI_SYMBOL, CM_SYMBOL, psi_s])
    names.update({km_s: r"K_M", kd: r"\mathbf k_d", msym: r"\mathbf M"})
    out = {}
    for d in DIMENSIONS:
        out["eikonal_psi_d%d.tex" % d] = emit_latex_fragment(r"\psi", psi_expr(d), names)
    out["eikonal_gradpsi.tex"] = emit_latex_fragment(r"\nabla\psi", kd*psi_s, names)
    out["eikonal_source.tex"] = emit_latex_fragment(
        r"S_\psi", psi_s*(C0*sp.sqrt(km_s) - km_s) - 1, names)
    quad = sp.Symbol("kMk", real=True)
    names[quad] = r"\mathbf k_d\cdot\mathbf M\,\mathbf k_d"
    out["eikonal_KM.tex"] = emit_latex_fragment(r"K_M", quad, names)
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


def _sym(name):
    return sp.Symbol(name, real=True)


def _vector_cxx(name):
    return sp.Matrix([_sym("%s.%s()" % (name, c)) for c in "xyz"])


def _tensor_cxx(name):
    return sp.Matrix(3, 3, lambda i, j: _sym("%s.%s%s()" % (name, "xyz"[i], "xyz"[j])))


def eikonal_header_regions(text):
    """{region name: C++ lines} of manufacturedEikonalReference.H."""
    ind = lambda name: codegen.region_indent(text, name)
    out = {}
    for d in DIMENSIONS:
        name = "manufacturedEikonalK_d%d" % d
        out[name] = codegen.call("vector", list(k_vector(d)), ind(name))
    name = "manufacturedEikonalM"
    M = metric(_tensor_cxx("conductivity"), _sym("chi"), _sym("Cm"))
    out[name] = codegen.call("tensor", list(M), ind(name))
    k = _vector_cxx("k")
    point = {x: _sym("point.x()"), y: _sym("point.y()"), z: _sym("point.z()")}
    psi = psi_of_k(k).subs(point, simultaneous=True)
    name = "manufacturedEikonalTau"
    out[name] = codegen.statement("return", psi, ind(name))
    name = "manufacturedEikonalGradTau"
    grad = [sp.diff(psi_of_k(k), a).subs(point, simultaneous=True) for a in AXES]
    lines, red = codegen.cse_lines(grad, ind(name))
    out[name] = lines + codegen.call("vector", red, ind(name))
    name = "manufacturedEikonalSourceTerm_KM"
    out[name] = codegen.assignment(
        "KM", k_metric(k, _tensor_cxx("M")), ind(name), decl="const scalar")
    name = "manufacturedEikonalSourceTerm"
    out[name] = codegen.statement(
        "return", source_of(_sym("tau"), _sym("KM"), _sym("c0")), ind(name))
    return out


HEADER_PATH = ("../../../src/verificationModels/eikonalVerification/"
               "manufacturedEikonalReference.H")
HEADER_PATHS = {"manufacturedEikonalReference.H": HEADER_PATH}
HEADERS = {"manufacturedEikonalReference.H": eikonal_header_regions}


def header_matches(path):
    return codegen.regions_match(__file__, path, HEADERS[os.path.basename(path)])


def main(argv=None):
    return codegen.generator_main(argv, __doc__, __file__, run_self_checks, write_outputs,
                                  HEADERS)


if __name__ == "__main__":
    sys.exit(main())
