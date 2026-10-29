#!/usr/bin/env python3
"""Eikonal ECG reference of the manufactured eikonal field.

For psi = exp(k_d . x) (cartesian_eikonal_mms) and a transmembrane template
U(s) the eikonal ECG at an electrode e is

    phi(e, t) = integral over [0, 1]^d of g_d dx,
    g_d = (r . sigma_g) / |r|^3,   r = x - e,   sigma_g = -sigma grad(psi) U'(t - psi),

where in d < 3 only the first d coordinates enter r.  The integrand is 0 where
|r|^2 <= 1e-300 (|dx| <= 1e-300 for d = 1), the guard of the C++ reference, so
an electrode on a quadrature node contributes nothing and never produces NaN.
The integral is a tensor-product Gauss-Legendre quadrature (tutorial order 48).

Templates (derivative U'):
    singleEvent   w = rise / (2 atanh(0.8)), a = tanh(s / w), b = tanh((s - duration) / w)
                  U' = 0.25/w ((1 - a^2)(1 - b) - (1 + a)(1 - b^2)),
                  duration = 64, rise = 4
    sine          U' = 2 pi cos(2 pi s)

The emitted NumPy module holds integrand_d<d>_<template>(x, y, z, t, ex, ey, ez,
sxx, sxy, sxz, syy, syz, szz, duration, rise), template_derivative_single_event,
template_derivative_sine and a numpy-only quadrature

    spatial_reference(d, template, sigma, electrode, t, order, duration, rise)
    reference(d, template, sigma, electrode, t, order, duration, rise)

(the same value; phi has no separate time factor).

Checks (--self-check): the integrand equals the C++ form and vanishes at the
electrode; the quadrature is converged between orders 48 and 96; the printed
single-event building blocks (H, U = H(s)(1 - H(s - D)), w, a, b) give the
same U' as template_derivative and as the derivative of U.

--write PATH rewrites, and --check PATH compares, the generated regions of the
hand-written manufacturedEikonalReference.H: the kernel (r . sigma_g)/|r|^3 of
manufacturedEikonalECGIntegrand for d = 1, 2, 3 (pseudo_ecg_mms.kernel).  The
guards, sigma_g and the template derivative stay hand-written.

    python eikonal_ecg_mms.py --emit
    python eikonal_ecg_mms.py --self-check
    python eikonal_ecg_mms.py --write ../../../src/verificationModels/eikonalVerification/manufacturedEikonalReference.H
    python eikonal_ecg_mms.py --check ../../../src/verificationModels/eikonalVerification/manufacturedEikonalReference.H
"""

import os
import sys

import numpy as np
import sympy as sp

import cartesian_eikonal_mms as eik
import mms_codegen as codegen
from pseudo_ecg_mms import kernel, kernel_code
from mms_common.latex_emit import (emit_latex_fragment, emit_parameter_fragment,
                                   paper_symbol_names)
from mms_common.numpy_emit import emit_numpy_module
from mms_common.quadrature import gauss_legendre_box

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "eikonal_ecg"
DIMENSIONS = eik.DIMENSIONS
TEMPLATES = ("single_event", "sine")
VSMALL = sp.Float("1e-300")
TUTORIAL_QUADRATURE_ORDER = 48
DEFAULT_DURATION = 64.0
DEFAULT_RISE = 4.0

x, y, z = eik.x, eik.y, eik.z
AXES = eik.AXES
t = sp.Symbol("t", real=True)
s = sp.Symbol("s", real=True)
DURATION = sp.Symbol("duration", positive=True)
RISE = sp.Symbol("rise", positive=True)
SIGMA_SYMBOLS = eik.SIGMA_SYMBOLS
SIGMA_SYMBOLIC = eik.SIGMA_SYMBOLIC
ELECTRODE_SYMBOLS = sp.symbols("ex ey ez", real=True)
ELECTRODE_SYMBOLIC = sp.Matrix(ELECTRODE_SYMBOLS)


def template_derivative(template):
    """U'(s) of a template as a sympy expression in s (duration, rise symbols)."""
    if template == "sine":
        return 2*sp.pi*sp.cos(2*sp.pi*s)
    if template == "single_event":
        w = RISE/(2*sp.atanh(sp.Rational(4, 5)))
        up = sp.tanh(s/w)
        down = sp.tanh((s - DURATION)/w)
        return sp.Rational(1, 4)/w*((1 - up**2)*(1 - down) - (1 + up)*(1 - down**2))
    raise ValueError("template must be one of %r, got %r" % (TEMPLATES, template))


def integrand(d, sigma, electrode, template="single_event"):
    """g_d as a sympy expression in x, y, z, t, with the VSMALL guard."""
    if d not in DIMENSIONS:
        raise ValueError("d must be 1, 2 or 3, got %r" % (d,))
    sigma = sp.Matrix(sigma)
    electrode = sp.Matrix(electrode)
    psi = eik.psi_expr(d)
    grad = eik.k_vector(d)*psi
    dU = template_derivative(template).subs(s, t - psi)
    flux = sigma*(-dU*grad)
    r = [AXES[m] - electrode[m] for m in range(d)]
    if d == 1:
        return sp.Piecewise((sp.Integer(0), sp.Abs(r[0]) <= VSMALL),
                            (kernel(d, flux, r), True))
    rsq = sum(c**2 for c in r)
    return sp.Piecewise((sp.Integer(0), rsq <= VSMALL), (kernel(d, flux, r), True))


def cxx_integrand(d, sigma, electrode, template="single_event"):
    """The integrand as the Cartesian component sums of the C++."""
    sigma = sp.Matrix(sigma)
    psi = eik.psi_expr(d)
    grad = eik.k_vector(d)*psi
    dU = template_derivative(template).subs(s, t - psi)
    sg = sigma*(-dU*grad)
    dx, dy, dz = (AXES[m] - electrode[m] for m in range(3))
    if d == 1:
        return sp.Piecewise((sp.Integer(0), sp.Abs(dx) <= VSMALL),
                            (sg[0]*dx/sp.Abs(dx)**3, True))
    if d == 2:
        rsq, num = dx**2 + dy**2, sg[0]*dx + sg[1]*dy
    else:
        rsq, num = dx**2 + dy**2 + dz**2, sg[0]*dx + sg[1]*dy + sg[2]*dz
    return sp.Piecewise((sp.Integer(0), rsq <= VSMALL), (num/rsq**sp.Rational(3, 2), True))


def reference(d, sigma, electrode, tt, n_quad, template="single_event",
              duration=DEFAULT_DURATION, rise=DEFAULT_RISE):
    """phi(e, tt) by n_quad-point Gauss-Legendre per axis."""
    expr = integrand(d, sigma, electrode, template).subs(
        {t: tt, DURATION: duration, RISE: rise})
    g = sp.lambdify(AXES, expr, "numpy")

    def f(p):
        full = [p[m] if m < d else np.zeros_like(p[0]) for m in range(3)]
        return g(*full)

    with np.errstate(divide="ignore", invalid="ignore"):
        return gauss_legendre_box(f, np.zeros(d), np.ones(d), n_quad)


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def run_self_checks():
    e = ELECTRODE_SYMBOLIC
    sigma = eik.sigma_matrix(eik.SIGMA_ROTATED)
    for d in DIMENSIONS:
        for template in TEMPLATES:
            a = integrand(d, sigma, e, template)
            b = cxx_integrand(d, sigma, e, template)
            sub = {eik.x: 0.31, eik.y: 0.62, eik.z: 0.17, t: 20.0, DURATION: 64, RISE: 4,
                   ELECTRODE_SYMBOLS[0]: 1.2, ELECTRODE_SYMBOLS[1]: -0.7,
                   ELECTRODE_SYMBOLS[2]: 1.4}
            _check(abs(float(a.subs(sub)) - float(b.subs(sub))) < 1e-12,
                   "d=%d %s: integrand differs from the C++ form" % (d, template))
        at = {x: 0.5, y: 0.5, z: 0.5, ELECTRODE_SYMBOLS[0]: 0.5,
              ELECTRODE_SYMBOLS[1]: 0.5, ELECTRODE_SYMBOLS[2]: 0.5}
        _check(integrand(d, sigma, e).subs(at) == 0, "d=%d: integrand not 0 at electrode" % d)
    el = (1.2, 0.23, 0.61)
    lo = reference(3, sigma, el, 30.0, TUTORIAL_QUADRATURE_ORDER)
    hi = reference(3, sigma, el, 30.0, 96)
    _check(abs(lo - hi) <= 1e-10*max(1.0, abs(hi)), "quadrature not converged: %r %r" % (lo, hi))
    _check_single_event_blocks()


# Source appended to the emitted module: a numpy-only quadrature over [0, 1]^d.
RUNTIME = '''

CASE_DIMENSION = {"d1": 1, "d2": 2, "d3": 3}


def _sigma_entries(sigma):
    s = numpy.asarray(sigma, dtype=float)
    if s.shape == (3, 3):
        return (s[0, 0], s[0, 1], s[0, 2], s[1, 1], s[1, 2], s[2, 2])
    return tuple(float(v) for v in s)


def spatial_reference(d, template, sigma, electrode, t, order,
                      duration=64.0, rise=4.0):
    """phi(electrode, t): the integral of integrand_d<d>_<template> over
    [0, 1]^d with `order` Gauss-Legendre points per axis.  template is
    "single_event" or "sine"; sigma is a 3x3 array or (xx, xy, xz, yy, yz, zz);
    electrode is (ex, ey, ez).  An electrode on a quadrature point adds 0."""
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
    f = globals()["integrand_d%d_%s" % (d, template)]
    with numpy.errstate(divide="ignore", invalid="ignore"):
        values = f(*coords, t, *electrode, *_sigma_entries(sigma), duration, rise)
    return float(numpy.sum(w*values))


reference = spatial_reference
'''


def numpy_module():
    functions = {
        "template_derivative_single_event": (
            (s, DURATION, RISE), template_derivative("single_event")),
        "template_derivative_sine": ((s,), template_derivative("sine")),
    }
    args = (x, y, z, t) + tuple(ELECTRODE_SYMBOLS) + tuple(SIGMA_SYMBOLS) + (DURATION, RISE)
    for d in DIMENSIONS:
        for template in TEMPLATES:
            functions["integrand_d%d_%s" % (d, template)] = (
                args, integrand(d, SIGMA_SYMBOLIC, ELECTRODE_SYMBOLIC, template))
    doc = ("Eikonal ECG integrands of the manufactured eikonal field on the unit cube "
           "and their Gauss-Legendre reference; integrand_d<d>_<template>(x, y, z, t, "
           "ex, ey, ez, sxx, sxy, sxz, syy, syz, szz, duration, rise).  Generated by "
           "eikonal_ecg_mms.py.")
    params = {"SIGMA_ROTATED": list(eik.SIGMA_ROTATED), "SIGMA_AXIS": list(eik.SIGMA_AXIS),
              "QUADRATURE_ORDER": TUTORIAL_QUADRATURE_ORDER, "DURATION": DEFAULT_DURATION,
              "RISE": DEFAULT_RISE, "VSMALL": 1e-300}
    return emit_numpy_module(REFERENCE_NAME, functions, params, doc).rstrip("\n") + "\n" + RUNTIME


_W = sp.Symbol("w", positive=True)
_A, _B = sp.symbols("a b", real=True)
_HS, _HSD = sp.symbols("Hs HsD", real=True)
_BLOCK_NAMES = {DURATION: "D", RISE: "r", _W: "w", _A: "a", _B: "b",
                _HS: "H(s)", _HSD: "H(s-D)"}


def _single_event_blocks():
    """The single-event template in the building blocks w, a = tanh(s/w),
    b = tanh((s - D)/w): H, U = H(s)(1 - H(s - D)) and U'."""
    H = (1 + sp.tanh(s/_W))/2
    U = _HS*(1 - _HSD)
    dU = sp.Rational(1, 4)/_W*((1 - _A**2)*(1 - _B) - (1 + _A)*(1 - _B**2))
    return {"H": H, "U": U, "dU": dU}


def _check_single_event_blocks():
    """U' of the building blocks equals template_derivative and the derivative of
    H(s)(1 - H(s - D))."""
    w = RISE/(2*sp.atanh(sp.Rational(4, 5)))
    blocks = _single_event_blocks()
    H = blocks["H"].subs(_W, w)
    dU = blocks["dU"].subs({_A: sp.tanh(s/_W), _B: sp.tanh((s - DURATION)/_W)}).subs(_W, w)
    U = H*(1 - H.subs(s, s - DURATION))
    at = {s: sp.Rational(37, 10), DURATION: 64, RISE: 4}
    for other in (template_derivative("single_event"), sp.diff(U, s)):
        _check(abs(float((dU - other).subs(at))) < 1e-14, "single-event building blocks")


def latex_fragments():
    rdot = sp.Symbol("rdot", real=True)
    rnorm = sp.Symbol("rnorm", real=True)
    u, up = sp.symbols("U Uprime", real=True)
    psi_s = sp.Symbol("psi", real=True)
    names = paper_symbol_names([psi_s])
    names.update({rdot: r"\mathbf r\cdot\boldsymbol\sigma\nabla\psi",
                  rnorm: r"\lVert\mathbf r\rVert", up: r"U'", u: r"U"})
    out = {
        "eikonal_ecg_integrand.tex": emit_latex_fragment(
            "g_d", -up*rdot/rnorm**3, names),
        "eikonal_ecg_template_sine.tex": emit_latex_fragment(
            r"U'(s)", template_derivative("sine"), {}),
        "eikonal_ecg_template_single_event.tex": emit_latex_fragment(
            r"U'(s)", _single_event_blocks()["dU"], _BLOCK_NAMES),
        "eikonal_ecg_template_U.tex": emit_latex_fragment(
            r"U(s)", _single_event_blocks()["U"], _BLOCK_NAMES),
        "eikonal_ecg_template_H.tex": emit_latex_fragment(
            r"H(s)", _single_event_blocks()["H"], _BLOCK_NAMES),
        "eikonal_ecg_template_w.tex": emit_latex_fragment(
            "w", RISE/(2*sp.atanh(sp.Rational(4, 5))), _BLOCK_NAMES),
        "eikonal_ecg_template_ab.tex": emit_parameter_fragment(
            [("a", sp.latex(sp.tanh(s/_W), symbol_names=_BLOCK_NAMES)),
             ("b", sp.latex(sp.tanh((s - DURATION)/_W), symbol_names=_BLOCK_NAMES))]),
        "eikonal_ecg_template_sine_U.tex": emit_latex_fragment(
            r"U(s)", sp.sin(2*sp.pi*s), {}),
        "eikonal_ecg_params.tex": emit_parameter_fragment(
            [("D", DEFAULT_DURATION), ("r", DEFAULT_RISE),
             (r"q_{\mathrm{ref}}", TUTORIAL_QUADRATURE_ORDER)], per_row=3),
    }
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


def eikonal_header_regions(text):
    """{region name: C++ lines} of manufacturedEikonalReference.H."""
    out = {}
    for d in DIMENSIONS:
        name = "manufacturedEikonalECGIntegrand_d%d" % d
        out[name] = kernel_code(d, "sigmaGradVm", codegen.region_indent(text, name))
    return out


HEADER_PATH = eik.HEADER_PATH
HEADER_PATHS = {"manufacturedEikonalReference.H": HEADER_PATH}
HEADERS = {"manufacturedEikonalReference.H": eikonal_header_regions}


def header_matches(path):
    return codegen.regions_match(__file__, path, HEADERS[os.path.basename(path)])


def main(argv=None):
    return codegen.generator_main(argv, __doc__, __file__, run_self_checks, write_outputs,
                                  HEADERS)


if __name__ == "__main__":
    sys.exit(main())
