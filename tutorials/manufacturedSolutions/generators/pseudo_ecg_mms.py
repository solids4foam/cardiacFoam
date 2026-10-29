#!/usr/bin/env python3
"""Pseudo-ECG reference of the manufactured monodomain fields.

For the fields V_m = sqrt(1 + t) F the pseudo-ECG at an electrode e is

    phi(e, t) = sqrt(1 + t) * integral over [0, 1]^d of g_d dx,
    g_d = (r . sigma grad F) / |r|^3,   r = x - e,

with F = F_d of the FDA fields (fda_monodomain_mms, sigma = diag FDA
conductivity) or F = sin^2(pi x) sin^2(pi y) sin^2(pi z) (rotated_monodomain_mms,
sigma = the rotated tensor).  In d < 3 only the first d coordinates enter r and
grad F.  The integrand is 0 where |r|^2 <= 1e-300, the guard of the C++
reference, so an electrode on a quadrature node contributes nothing and never
produces NaN.  The integral is a tensor-product Gauss-Legendre quadrature.

The emitted NumPy module holds integrand_d1, integrand_d2, integrand_d3 and
integrand_rotated (arguments x, y, z, ex, ey, ez and the six distinct entries
of sigma: xx xy xz yy yz zz), a numpy-only quadrature spatial_reference, and
time_factor.  Signatures, in the order of the generator functions
spatial_reference(d, sigma, electrode, n_quad) and
reference(d, sigma, electrode, t, n_quad):

    spatial_reference(case, sigma, electrode, order)
    reference(case, sigma, electrode, t, order)
    time_factor(t)

with case one of "d1", "d2", "d3", "rotated".

Checks (--self-check): the integrand equals the C++ form (Cartesian sums for
d = 1, 2, 3) and vanishes at the electrode.

--write PATH rewrites, and --check PATH compares, the generated regions of a
hand-written header, chosen by its file name: in manufacturedFDAReference.H
the time factor sqrt(1 + t), grad F_d and the kernel (r . sigma grad F)/|r|^3
for d = 1, 2, 3; in manufacturedAnisotropicMonodomainReference.H grad F and
the kernel of the rotated case.  The guards, the product sigma grad F and the
quadrature stay hand-written.

    python pseudo_ecg_mms.py --emit
    python pseudo_ecg_mms.py --self-check
    python pseudo_ecg_mms.py --write ../../../src/verificationModels/monodomainVerification/manufacturedFDAReference.H
    python pseudo_ecg_mms.py --write ../../../src/verificationModels/monodomainVerification/manufacturedAnisotropicMonodomainReference.H
    python pseudo_ecg_mms.py --check <either header>
"""

import os
import sys

import numpy as np
import sympy as sp

import fda_monodomain_mms as fda
import mms_codegen as codegen
import rotated_monodomain_mms as rot
from mms_common.latex_emit import (emit_latex_fragment, emit_parameter_fragment,
                                   paper_symbol_names)
from mms_common.numpy_emit import emit_numpy_module
from mms_common.quadrature import gauss_legendre_box

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "pseudo_ecg"
DIMENSIONS = (1, 2, 3)
VSMALL = sp.Float("1e-300")
TUTORIAL_QUADRATURE_ORDER = 96

x, y, z, t = fda.x, fda.y, fda.z, fda.t
AXES = (x, y, z)

_SIGMA_NAMES = "sxx sxy sxz syy syz szz".split()
SIGMA_SYMBOLS = sp.symbols(_SIGMA_NAMES, real=True)
ELECTRODE_SYMBOLS = sp.symbols("ex ey ez", real=True)


def sigma_matrix(entries):
    """Symmetric 3x3 matrix from (xx, xy, xz, yy, yz, zz)."""
    a, b, c, d, e, f = entries
    return sp.Matrix([[a, b, c], [b, d, e], [c, e, f]])


SIGMA_SYMBOLIC = sigma_matrix(SIGMA_SYMBOLS)
ELECTRODE_SYMBOLIC = sp.Matrix(ELECTRODE_SYMBOLS)
SIGMA_FDA = sp.diag(*fda.SIGMA_DIAG)


def grad_F(d, F=None):
    """Gradient of F (default: the FDA F_d) as a 3-vector, zero beyond d."""
    F = fda.fields(d)["F"] if F is None else F
    return sp.Matrix([sp.diff(F, AXES[m]) if m < d else 0 for m in range(3)])


def kernel(d, flux, r):
    """(r . flux)/|r|^3 over the first d components, without the guard."""
    dot = sum(r[m]*flux[m] for m in range(d))
    if d == 1:
        return dot/sp.Abs(r[0])**3
    return dot/sum(c**2 for c in r[:d])**sp.Rational(3, 2)


def integrand(d, sigma, electrode, F=None):
    """g_d as a sympy expression in x, y, z, with the VSMALL guard."""
    if d not in DIMENSIONS:
        raise ValueError("d must be 1, 2 or 3, got %r" % (d,))
    sigma = sp.Matrix(sigma)
    electrode = sp.Matrix(electrode)
    flux = sigma*grad_F(d, F)
    r = [AXES[m] - electrode[m] for m in range(d)]
    if d == 1:
        # the C++ guard of d = 1 is on |dx|, not on dx^2
        return sp.Piecewise((sp.Integer(0), sp.Abs(r[0]) <= VSMALL),
                            (kernel(d, flux, r), True))
    rsq = sum(c**2 for c in r)
    return sp.Piecewise((sp.Integer(0), rsq <= VSMALL), (kernel(d, flux, r), True))


def time_factor(tt):
    return np.sqrt(1.0 + tt)


def spatial_reference(d, sigma, electrode, n_quad, F=None):
    """Integral of g_d over [0, 1]^d by n_quad-point Gauss-Legendre per axis."""
    g = sp.lambdify(AXES, integrand(d, sigma, electrode, F), "numpy")

    def f(p):
        full = [p[m] if m < d else np.zeros_like(p[0]) for m in range(3)]
        return g(*full)

    with np.errstate(divide="ignore", invalid="ignore"):
        return gauss_legendre_box(f, np.zeros(d), np.ones(d), n_quad)


def reference(d, sigma, electrode, t, n_quad, F=None):
    """phi(e, t) = sqrt(1 + t) * spatial_reference."""
    return float(time_factor(t)*spatial_reference(d, sigma, electrode, n_quad, F))


def _rotated_F():
    return rot.fields(rot.K_CASE)["F"]


def cxx_integrand(d, sigma, electrode, F=None):
    """The integrand written as the Cartesian component sums of the C++."""
    g = grad_F(d, F)
    s = sp.Matrix(sigma)
    sg = [sum(s[i, j]*g[j] for j in range(3)) for i in range(3)]
    dx, dy, dz = (AXES[m] - electrode[m] for m in range(3))
    if d == 1:
        return sp.Piecewise((sp.Integer(0), sp.Abs(dx) <= VSMALL),
                            (sg[0]*dx/sp.Abs(dx)**3, True))
    if d == 2:
        rsq = dx**2 + dy**2
        num = sg[0]*dx + sg[1]*dy
    else:
        rsq = dx**2 + dy**2 + dz**2
        num = sg[0]*dx + sg[1]*dy + sg[2]*dz
    return sp.Piecewise((sp.Integer(0), rsq <= VSMALL), (num/rsq**sp.Rational(3, 2), True))


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def run_self_checks():
    e = ELECTRODE_SYMBOLIC
    pts = [(0.31, 0.62, 0.17, 0.2), (0.9, 0.05, 0.44, -0.3)]
    for d in DIMENSIONS:
        for sigma in (SIGMA_FDA, rot.K_CASE):
            a = integrand(d, sigma, e)
            b = cxx_integrand(d, sigma, e)
            for p in pts:
                sub = {x: p[0], y: p[1], z: p[2], ELECTRODE_SYMBOLS[0]: p[3],
                       ELECTRODE_SYMBOLS[1]: -0.7, ELECTRODE_SYMBOLS[2]: 1.4}
                _check(abs(float(a.subs(sub)) - float(b.subs(sub))) < 1e-12,
                       "d=%d: integrand differs from the C++ form" % d)
        at = {x: 0.5, y: 0.5, z: 0.5, ELECTRODE_SYMBOLS[0]: 0.5,
              ELECTRODE_SYMBOLS[1]: 0.5, ELECTRODE_SYMBOLS[2]: 0.5}
        _check(integrand(d, SIGMA_FDA, e).subs(at) == 0, "d=%d: integrand not 0 at electrode" % d)


# Source appended to the emitted module: a numpy-only quadrature over [0, 1]^d.
RUNTIME = '''

CASE_DIMENSION = {"d1": 1, "d2": 2, "d3": 3, "rotated": 3}


def _sigma_entries(sigma):
    s = numpy.asarray(sigma, dtype=float)
    if s.shape == (3, 3):
        return (s[0, 0], s[0, 1], s[0, 2], s[1, 1], s[1, 2], s[2, 2])
    return tuple(float(v) for v in s)


def time_factor(t):
    return numpy.sqrt(1.0 + t)


def spatial_reference(case, sigma, electrode, order):
    """Integral of the integrand over [0, 1]^d (d from case: d1, d2, d3,
    rotated) with `order` Gauss-Legendre points per axis.  sigma is a 3x3
    array or (xx, xy, xz, yy, yz, zz); electrode is (ex, ey, ez)."""
    d = CASE_DIMENSION[case]
    nodes, weights = numpy.polynomial.legendre.leggauss(order)
    nodes = 0.5*(nodes + 1.0)
    weights = 0.5*weights
    grid = numpy.meshgrid(*([nodes]*d), indexing="ij")
    w = numpy.ones_like(grid[0])
    for i in range(d):
        shape = [1]*d
        shape[i] = order
        w = w*weights.reshape(shape)
    coords = list(grid) + [numpy.zeros_like(grid[0])]*(3 - d)
    f = globals()["integrand_" + case]
    with numpy.errstate(divide="ignore", invalid="ignore"):
        values = f(*coords, *electrode, *_sigma_entries(sigma))
    return float(numpy.sum(w*values))


def reference(case, sigma, electrode, t, order):
    return time_factor(t)*spatial_reference(case, sigma, electrode, order)
'''


def numpy_module():
    functions = {}
    args = (x, y, z) + tuple(ELECTRODE_SYMBOLS) + tuple(SIGMA_SYMBOLS)
    for d in DIMENSIONS:
        functions["integrand_d%d" % d] = (
            args, integrand(d, SIGMA_SYMBOLIC, ELECTRODE_SYMBOLIC))
    functions["integrand_rotated"] = (
        args, integrand(3, SIGMA_SYMBOLIC, ELECTRODE_SYMBOLIC, _rotated_F()))
    doc = ("Pseudo-ECG integrands of the manufactured monodomain fields on the unit "
           "cube and their Gauss-Legendre reference; integrand_<case>(x, y, z, ex, ey, "
           "ez, sxx, sxy, sxz, syy, syz, szz).  Generated by pseudo_ecg_mms.py.")
    params = {"SIGMA_FDA": [float(v) for v in fda.SIGMA_DIAG],
              "SIGMA_ROTATED": [float(rot.K_CASE[i, j]) for i, j in
                                ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))],
              "QUADRATURE_ORDER": TUTORIAL_QUADRATURE_ORDER, "VSMALL": 1e-300}
    return emit_numpy_module(REFERENCE_NAME, functions, params, doc).rstrip("\n") + "\n" + RUNTIME


def latex_fragments():
    names_syms = {n: sp.Symbol(n, real=True) for n in ("G",)}
    rdot = sp.Symbol("rdot", real=True)
    rnorm = sp.Symbol("rnorm", real=True)
    names = paper_symbol_names(names_syms.values())
    names.update({rdot: r"\mathbf r\cdot\boldsymbol\sigma\nabla F_d",
                  rnorm: r"\lVert\mathbf r\rVert"})
    out = {}
    for d in DIMENSIONS:
        out["pseudo_ecg_gradF_d%d.tex" % d] = emit_latex_fragment(
            r"\nabla F_%d" % d, grad_F(d)[:d, :], names)
    out["pseudo_ecg_gradF_rotated.tex"] = emit_latex_fragment(
        r"\nabla F", grad_F(3, _rotated_F()), names)
    out["pseudo_ecg_integrand.tex"] = emit_latex_fragment(
        "g_d", rdot/rnorm**3, names)
    out["pseudo_ecg_time_factor.tex"] = emit_latex_fragment(
        r"\sqrt{1+t}", sp.sqrt(1 + t), names)
    integral = sp.Symbol("integral", real=True)
    phi0 = sp.Symbol("phi0", real=True)
    out["pseudo_ecg_functional.tex"] = emit_latex_fragment(
        r"\Phi_e(\mathbf r_e,t)", sp.sqrt(1 + t)*phi0,
        {**names, phi0: r"\Phi_e(\mathbf r_e,0)"})
    out["pseudo_ecg_spatial.tex"] = emit_latex_fragment(
        r"\Phi_e(\mathbf r_e,0)", integral,
        {**names, integral: r"\int_{[0,1]^d} g_d\,\mathrm d\mathbf x"})
    out["pseudo_ecg_params.tex"] = emit_parameter_fragment(
        [(r"q_{\mathrm{ref}}", TUTORIAL_QUADRATURE_ORDER),
         (r"g_d\big|_{\lVert\mathbf r\rVert^2\le 10^{-300}}", 0)], per_row=2)
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


def kernel_code(d, flux_name, indent):
    """C++ `return` of the kernel: flux components flux_name.x() ..., r
    components dx, dy, dz, and the guard's |dx| (d = 1) or |r|^2 (d > 1)
    as absDx or rSq."""
    flux = [_sym("%s.%s()" % (flux_name, c)) for c in "xyz"]
    r = [_sym(c) for c in ("dx", "dy", "dz")]
    expr = kernel(d, flux, r)
    if d == 1:
        expr = expr.subs(sp.Abs(r[0]), _sym("absDx"))
    else:
        expr = expr.subs(sum(c**2 for c in r[:d]), _sym("rSq"))
    return codegen.statement("return", expr, indent)


def _grad_code(g, indent):
    return codegen.call("vector", list(g), indent)


def fda_header_regions(text):
    """{region name: C++ lines} of manufacturedFDAReference.H."""
    ind = lambda name: codegen.region_indent(text, name)
    out = {"pseudoECGManufacturedTimeFactor": codegen.statement(
        "return", fda.vm_of(sp.Integer(1)), ind("pseudoECGManufacturedTimeFactor"))}
    for d in DIMENSIONS:
        name = "pseudoECGManufacturedGradF_d%d" % d
        out[name] = _grad_code(grad_F(d), ind(name))
        name = "pseudoECGManufacturedIntegrand_d%d" % d
        out[name] = kernel_code(d, "sigmaGradF", ind(name))
    return out


def anisotropic_header_regions(text):
    """{region name: C++ lines} of manufacturedAnisotropicMonodomainReference.H."""
    ind = lambda name: codegen.region_indent(text, name)
    name = "pseudoECGManufacturedGradFAnisotropic"
    lines, red = codegen.cse_lines(list(grad_F(3, _rotated_F())), ind(name))
    out = {name: lines + _grad_code(red, ind(name))}
    name = "pseudoECGManufacturedIntegrandAnisotropic"
    out[name] = kernel_code(3, "sigmaGradF", ind(name))
    return out


_SRC = "../../../src/verificationModels/monodomainVerification/"
HEADER_PATHS = {
    "manufacturedFDAReference.H": _SRC + "manufacturedFDAReference.H",
    "manufacturedAnisotropicMonodomainReference.H":
        _SRC + "manufacturedAnisotropicMonodomainReference.H",
}
HEADERS = {
    "manufacturedFDAReference.H": fda_header_regions,
    "manufacturedAnisotropicMonodomainReference.H": anisotropic_header_regions,
}


def header_matches(path):
    return codegen.regions_match(__file__, path, HEADERS[os.path.basename(path)])


def main(argv=None):
    return codegen.generator_main(argv, __doc__, __file__, run_self_checks, write_outputs,
                                  HEADERS)


if __name__ == "__main__":
    sys.exit(main())
