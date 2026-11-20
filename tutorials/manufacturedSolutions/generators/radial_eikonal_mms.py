#!/usr/bin/env python3
"""Exact source-free eikonal solution in a sealed sphere around a held ball.

With an isotropic M = m I and a field of the radius alone, the eikonal operator
of the solver

    -div(M grad psi) + c0 sqrt(grad psi . M grad psi) = 1

reduces to -m (psi'' + 2 psi'/r) + V psi' = 1, V = c0 sqrt(m), with psi(R) = Tk on
the ball and psi'(Ro) = 0 on the sealed wall. With L = Ro - R, s = r/L, a = R/L,
b = Ro/L, Pe = c0 L/sqrt(m) and w = V psi' this is the linear equation

    (w' + 2 w/s)/Pe = w - 1,   w(b) = 0,

so that, with q(s) = s^2 + 2 s/Pe + 2/Pe^2,

    w(s) = (q(s) - q(b) exp(-Pe (b - s)))/s^2,
    psi(s) = Tk + (L/V) (W(s) - W(a)),
    W(s) = s + 2 log(s)/Pe - 2/(Pe^2 s)
           - q(b) exp(-Pe b) (Pe Ei(Pe s) - exp(Pe s)/s).

The solution extends below R and is singular at r = 0.

A constant anisotropic M = m Q diag(a^2) Q^T keeps the solution exact in the
coordinate xi = diag(1/a) Q^T x, with r = |xi|: the domain and the ball are
ellipsoids, and the sealed condition n.M grad(psi) = 0 maps to psi'(Ro) = 0.

The three seed profiles of the eikonal PVJ coupler (flat, front, curvedFront; the paper
names them flat, front-time and exact-radial) define three problems for this field, which
differ in the Dirichlet data only: Tk at R (flat), Tk + R/V at R (front), and Tk at the
critical radius 2 sqrt(m)/c0, i.e. s = 2/Pe (curvedFront, solved in the field outside it).
reference/radial_eikonal.py evaluates all three, the field in the stretched coordinate, its
gradient, and the seed offsets pvjMapper::setSeedOffsets gives each support cell.

The reference does not use the exponential integral Ei.  As the C++ header
manufacturedRadialEikonalReference.H does, it writes the Ei term of W as the integral of
exp(-Pe (b - t))/t^2 over [s, b] and evaluates it by Gauss-Legendre on panels that follow the
exponential layer at the wall, so it needs NumPy only and does not overflow where the closed
form does (Pe (a + 1) above about 700, where exp(Pe s) leaves a double).  The closed form
here, which uses scipy.special.expi, is the independent check of that evaluation.
test_radial_eikonal_mms.py checks the samples of the C++ reference.

    python radial_eikonal_mms.py --report --PeL 46 --ballRatio 0.25
    python radial_eikonal_mms.py --self-check
    python radial_eikonal_mms.py --emit
"""

import argparse
import os
import sys

import numpy as np
from scipy.special import expi

import mms_codegen as codegen
from mms_common.numpy_emit import emit_numpy_module

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_NAME = "radial_eikonal"

# The closed form needs Pe (a + 1) below this: exp(Pe s) still fits a double.
CLOSED_FORM_LIMIT = 700.0


def q(s, Pe):
    return s*s + 2*s/Pe + 2/Pe**2


def slope(s, Pe, a):
    """w = V dpsi/dr, 1 for a plane front."""
    b = a + 1
    return (q(s, Pe) - q(b, Pe)*np.exp(-Pe*(b - s)))/s**2


def scaled_time(s, Pe, a):
    """V (psi - Tk)/L, zero at s = a (closed form, scipy.special.expi)."""
    b = a + 1

    def W(t):
        return (
            t + 2*np.log(t)/Pe - 2/(Pe*Pe*t)
          - q(b, Pe)*np.exp(-Pe*b)*(Pe*expi(Pe*t) - np.exp(Pe*t)/t)
        )

    return W(s) - W(a)


def rad_arrival(Pe, a):
    return float(scaled_time(np.array(a + 1.0), Pe, a))


# ---- the NumPy reference ------------------------------------------------------
# The module reference/radial_eikonal.py is this text after a numpy-only header.  It
# repeats the C++ evaluation (wallIntegral, scaledTime), so a solver that cannot call
# scipy.special.expi, or whose Pe (a + 1) exceeds CLOSED_FORM_LIMIT, still has the field.

_REFERENCE_FUNCTIONS = '''

CRITICAL_RADIUS = 2.0
PROFILES = ("flat", "front", "exact", "curvedFront")

_GAUSS_X, _GAUSS_W = numpy.polynomial.legendre.leggauss(8)


def parameters(c0, m, L):
    """(V, Pe) of the problem: V = c0 sqrt(m) and Pe = c0 L / sqrt(m), with m the metric
    component of M = m I (sigma/(chi Cm) in the solver) and L the wall thickness Ro - R."""
    root = numpy.sqrt(m)
    return c0*root, c0*L/root


def q(s, Pe):
    return s*s + 2*s/Pe + 2/Pe**2


def slope(s, Pe, a):
    """w = V dpsi/dr at s = r/L, for a ball of radius a = R/L: 1 for a plane front, 0 at
    the sealed wall s = a + 1.  Solves (w' + 2 w/s)/Pe = w - 1."""
    b = a + 1
    s = numpy.asarray(s, dtype=float)
    return (q(s, Pe) - q(b, Pe)*numpy.exp(-Pe*(b - s)))/s**2


def wall_integral(s, Pe, a):
    """Integral from s to b = a + 1 of exp(-Pe (b - t))/t^2 dt: Gauss-Legendre (8 points) on
    panels that follow the exponential layer at the wall and the growth of 1/t^2 at small t,
    as manufacturedRadialEikonalReference.H does."""
    b = a + 1.0
    length = b - numpy.asarray(s, dtype=float)
    tau_low = numpy.zeros_like(length)
    total = numpy.zeros_like(length)
    scale = 2.0/Pe
    with numpy.errstate(divide="ignore", invalid="ignore", over="ignore"):
        while numpy.any(tau_low < length):
            active = tau_low < length
            step = numpy.minimum(scale, 0.5*(b - tau_low))
            tau_high = numpy.minimum(tau_low + step, length)
            half = 0.5*(tau_high - tau_low)
            tau = 0.5*(tau_high + tau_low)[..., None] + half[..., None]*_GAUSS_X
            t = b - tau
            panel = half*numpy.sum(_GAUSS_W*numpy.exp(-Pe*tau)/(t*t), axis=-1)
            total = total + numpy.where(active, panel, 0.0)
            tau_low = numpy.where(active, tau_high, tau_low)
            scale *= 2.0
    return total


def scaled_time(s, Pe, a):
    """V (psi - Tk)/L at s = r/L: zero at s = a, the solution of the radial problem with
    psi(R) = Tk.  Valid for s > 0, also below the ball (s < a)."""
    b = a + 1.0
    s = numpy.asarray(s, dtype=float)

    def regular(t):
        return t + 2*numpy.log(t)/Pe - 2/(Pe*Pe*t)

    return (regular(s) - regular(a)
            - q(b, Pe)*(wall_integral(a, Pe, a) - wall_integral(s, Pe, a)))


def exact_time(s, Pe, a, L, V, terminal_time=0.0, profile="flat"):
    """Exact arrival time at s = r/L of the problem of a seed profile.

    flat (and exact)  Tk + (L/V) W(s)                       psi(R) = Tk
    front             Tk + (L/V) (a + W(s))                 psi(R) = Tk + R/V
    curvedFront       Tk + (L/V) (W(max(s, 2/Pe)) - W(2/Pe))  psi = Tk at r = 2 sqrt(m)/c0

    with W = scaled_time(., Pe, a).  curvedFront needs a >= 2/Pe, as the C++ verifier
    requires: the field holds outside the critical radius, which the ball must contain.
    The maximum only acts below that radius, where the solver holds Tk.
    """
    if profile not in PROFILES:
        raise ValueError("profile %r not in %s" % (profile, ", ".join(PROFILES)))
    s = numpy.asarray(s, dtype=float)
    if profile == "curvedFront":
        if 2.0/Pe > a:
            raise ValueError("profile curvedFront needs ballRatio a >= 2/Pe = %g, got %g"
                             % (2.0/Pe, a))
        core = CRITICAL_RADIUS/Pe
        return terminal_time + (L/V)*(scaled_time(numpy.maximum(s, core), Pe, a)
                                      - scaled_time(core, Pe, a))
    edge = a if profile == "front" else 0.0
    return terminal_time + (L/V)*(edge + scaled_time(s, Pe, a))


def rotation_about_diagonal(angle):
    """Rotation Q by angle about the axis (1, 1, 1)/sqrt(3)."""
    k = 1.0/numpy.sqrt(3.0)
    K = numpy.array([[0.0, -k, k], [k, 0.0, -k], [-k, k, 0.0]])
    return numpy.eye(3) + numpy.sin(angle)*K + (1.0 - numpy.cos(angle))*(K @ K)


def stretched(x, axes=(1.0, 1.0, 1.0), angle=0.0):
    """xi = diag(1/axes) Q^T x for points x of shape (..., 3), the ball centred at the
    origin: the coordinate in which the field depends on r = |xi| alone, for the
    conductivity M = m Q diag(axes^2) Q^T."""
    local = numpy.asarray(x, dtype=float) @ rotation_about_diagonal(angle)
    return local/numpy.asarray(axes, dtype=float)


def exact_field(x, Pe, a, L, V, terminal_time=0.0, profile="flat",
                axes=(1.0, 1.0, 1.0), angle=0.0):
    """Exact arrival time at the points x (..., 3).  Points inside the ball (|xi| < R)
    are given the value at the ball surface, as the verifier's error evaluation does."""
    s = numpy.linalg.norm(stretched(x, axes, angle), axis=-1)/L
    return exact_time(numpy.maximum(s, a), Pe, a, L, V, terminal_time, profile)


def exact_gradient(x, Pe, a, L, V, axes=(1.0, 1.0, 1.0), angle=0.0):
    """Gradient of the exact field in x, shape (..., 3): Q diag(1/axes) (w/V) xi/|xi|;
    zero inside the ball (|xi| <= R).  It does not depend on the seed profile."""
    xi = stretched(x, axes, angle)
    r = numpy.linalg.norm(xi, axis=-1)
    s = r/L
    free = s > a
    safe = numpy.where(free, s, a + 1.0)
    radial = numpy.where(free, slope(safe, Pe, a)/V, 0.0)/numpy.where(free, r, 1.0)
    d_xi = radial[..., None]*xi
    return (d_xi/numpy.asarray(axes, dtype=float)) @ rotation_about_diagonal(angle).T


def metric_radius(d, M):
    """rho = sqrt(max(d . M^-1 d, 0)) for offsets d (..., 3) from the terminal."""
    d = numpy.asarray(d, dtype=float)
    return numpy.sqrt(numpy.maximum(
        numpy.einsum("...i,ij,...j->...", d, numpy.linalg.inv(numpy.asarray(M, dtype=float)), d),
        0.0))


def seed_offset(rho, c0, profile="front"):
    """Seed time of a support cell relative to its terminal, pvjMapper::setSeedOffsets.

    flat 0;  front rho/c0;  curvedFront the radial source-free solution of an unbounded
    medium measured from rho_c = 2/c0 (zero for rho <= rho_c):
    (rho - rho_c)/c0 + 2 log(rho/rho_c)/c0^2 - 2 (1/rho - 1/rho_c)/c0^3.
    rho is the metric distance sqrt(d . M^-1 d) of the cell centre from the terminal."""
    rho = numpy.asarray(rho, dtype=float)
    if profile == "flat":
        return numpy.zeros_like(rho)
    if profile == "front":
        return rho/c0
    if profile != "curvedFront":
        raise ValueError("profile %r not in flat, front, curvedFront" % (profile,))
    rho_c = CRITICAL_RADIUS/c0
    past = rho > rho_c
    r = numpy.where(past, rho, rho_c)
    return numpy.where(
        past, (r - rho_c)/c0 + 2*numpy.log(r/rho_c)/c0**2 - 2*(1/r - 1/rho_c)/c0**3, 0.0)


def seed_times(centres, terminal, terminal_time, M, c0, profile="front"):
    """Seed time of every support cell: terminal_time + seed_offset of its metric distance
    from the terminal (centres (N, 3), M the 3 x 3 activation metric of the tissue)."""
    rho = metric_radius(numpy.asarray(centres, dtype=float) - numpy.asarray(terminal), M)
    return terminal_time + seed_offset(rho, c0, profile)
'''


def numpy_module():
    doc = ("Exact source-free eikonal field around a held ball in a sealed sphere, for the "
           "flat, front and curvedFront seed profiles of the eikonal PVJ coupler: "
           "scaled_time, slope, exact_time, exact_field, exact_gradient, seed_offset, "
           "seed_times.  NumPy only; the wall integral is evaluated by Gauss-Legendre "
           "panels as in manufacturedRadialEikonalReference.H.  Generated by "
           "radial_eikonal_mms.py.")
    base = emit_numpy_module(REFERENCE_NAME, {}, {"CRITICAL_RADIUS": 2.0}, doc)
    return base.rstrip("\n") + "\n" + _REFERENCE_FUNCTIONS


def reference_namespace():
    """The emitted module executed in a fresh namespace."""
    namespace = {"__name__": "radial_eikonal_reference"}
    exec(compile(numpy_module(), "<radial_eikonal reference>", "exec"), namespace)
    return namespace


# ---- self-checks --------------------------------------------------------------

def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def _panel_integral(f, lo, hi, panels, n=8):
    """Composite Gauss-Legendre integral of f over [lo, hi]."""
    x, w = np.polynomial.legendre.leggauss(n)
    edges = np.linspace(lo, hi, panels + 1)
    half = 0.5*np.diff(edges)
    mid = 0.5*(edges[1:] + edges[:-1])
    return float(np.sum(half[:, None]*w[None, :]*f(mid[:, None] + half[:, None]*x[None, :])))


def run_self_checks():
    ref = reference_namespace()
    ref_slope, ref_time = ref["slope"], ref["scaled_time"]

    for Pe in (2.0, 5.0, 21.0, 46.0, 94.0):
        for a in (0.02, 0.1, 0.25, 0.5):
            if Pe*(a + 1) > CLOSED_FORM_LIMIT:
                continue
            s = np.linspace(a, a + 1, 41)[1:]
            # the panel quadrature is accurate to 1.5e-10 at worst (Pe 2, a 0.02), the C++ design
            _check(np.max(np.abs(ref_time(s, Pe, a) - scaled_time(s, Pe, a))) < 3e-10,
                   "reference time differs from the closed form at Pe %g, a %g" % (Pe, a))
            _check(np.max(np.abs(ref_slope(s, Pe, a) - slope(s, Pe, a))) < 1e-13,
                   "reference slope differs from the closed form at Pe %g, a %g" % (Pe, a))
            if a >= 0.1:     # the central difference is too coarse where w ~ 1/s^2 is steep
                ds = 1e-5
                sm = s[:-1]
                w = ref_slope(sm, Pe, a)
                dw = (ref_slope(sm + ds, Pe, a) - ref_slope(sm - ds, Pe, a))/(2*ds)
                _check(np.max(np.abs((dw + 2*w/sm)/Pe - (w - 1))) < 1e-5,
                       "slope does not satisfy (w' + 2 w/s)/Pe = w - 1 at Pe %g, a %g" % (Pe, a))
            _check(abs(float(ref_slope(np.array(a + 1.0), Pe, a))) < 1e-14,
                   "slope does not vanish at the sealed wall")
            _check(abs(float(ref_time(np.array(a), Pe, a))) < 1e-14,
                   "time does not vanish on the ball")

    # Where the closed form overflows the reference still integrates the slope
    for Pe, a in ((1000.0, 0.25), (3000.0, 0.1)):
        _check(Pe*(a + 1) > CLOSED_FORM_LIMIT, "overflow case lies inside the closed form")
        arrival = float(ref_time(np.array(a + 1.0), Pe, a))
        integral = _panel_integral(lambda s: ref_slope(s, Pe, a), a, a + 1, 4000)
        _check(np.isfinite(arrival) and abs(arrival - integral) < 1e-8,
               "reference time at Pe %g differs from the integral of the slope" % Pe)

    # The seed problems: curvedFront is the flat field re-based at the critical radius
    Pe, a, L, V = 46.0, 0.25, 0.1, 0.8
    s = np.array([0.3, 0.7, 1.2])
    flat = ref["exact_time"](s, Pe, a, L, V, 0.5, "flat")
    front = ref["exact_time"](s, Pe, a, L, V, 0.5, "front")
    curved = ref["exact_time"](s, Pe, a, L, V, 0.5, "curvedFront")
    _check(np.allclose(front - flat, L/V*a, rtol=0, atol=1e-14), "front differs from flat by R/V")
    _check(np.allclose(curved - flat, -L/V*ref_time(2/Pe, Pe, a), rtol=0, atol=1e-14),
           "curvedFront differs from flat by the time at the critical radius")
    _check(abs(float(ref["exact_time"](a, Pe, a, L, V, 0.5, "flat")) - 0.5) < 1e-14,
           "flat field is not Tk on the ball")
    try:
        ref["exact_time"](s, 5.0, 0.1, L, V, 0.0, "curvedFront")
    except ValueError:
        pass
    else:
        raise AssertionError("curvedFront with a < 2/Pe was accepted")

    # An anisotropic M = m Q diag(axes^2) Q^T leaves the field exact in the stretched
    # coordinate: the eikonal residual of the field in x vanishes, and so does the
    # conormal flux on the ellipsoid that bounds the domain
    c0, m, axes, angle = 5.0, 0.5, (1.3, 0.8, 1.1), 0.7
    Q = ref["rotation_about_diagonal"](angle)
    M = m*Q @ np.diag(np.square(axes)) @ Q.T
    V, Pe = ref["parameters"](c0, m, L)
    rng = np.random.default_rng(7)
    xi = rng.normal(size=(12, 3))
    xi *= (L*rng.uniform(1.3*a, 0.9*(a + 1), 12)/np.linalg.norm(xi, axis=1))[:, None]
    x = (xi*np.asarray(axes)) @ Q.T
    h = 1e-5
    grad = ref["exact_gradient"](x, Pe, a, L, V, axes, angle)
    for i in range(3):
        e = np.zeros(3)
        e[i] = h
        fd = (ref["exact_field"](x + e, Pe, a, L, V, 0.0, "flat", axes, angle)
              - ref["exact_field"](x - e, Pe, a, L, V, 0.0, "flat", axes, angle))/(2*h)
        _check(np.max(np.abs(fd - grad[:, i])) < 1e-6, "gradient differs from finite differences")
    div = 0.0
    for i in range(3):
        e = np.zeros(3)
        e[i] = h
        flux = lambda p: ref["exact_gradient"](p, Pe, a, L, V, axes, angle) @ M[:, i]
        div = div + (flux(x + e) - flux(x - e))/(2*h)
    residual = -div + c0*np.sqrt(np.einsum("ni,ij,nj->n", grad, M, grad)) - 1.0
    _check(np.max(np.abs(residual)) < 1e-5, "anisotropic field does not solve the eikonal equation")
    directions = np.array([[1.0, 0.0, 0.0], [0.3, 0.8, -0.5], [-0.2, 0.1, 0.9]])
    xi_wall = directions/np.linalg.norm(directions, axis=1)[:, None]*L*(a + 1)
    wall = (xi_wall*np.asarray(axes)) @ Q.T
    normal = (xi_wall/np.asarray(axes)) @ Q.T       # grad_x |xi|^2 up to a factor 2
    flux = ref["exact_gradient"](wall, Pe, a, L, V, axes, angle) @ M
    _check(np.max(np.abs(np.einsum("ni,ni->n", normal, flux))) < 1e-12,
           "conormal flux does not vanish on the wall ellipsoid")

    # The offsets of pvjMapper::setSeedOffsets: the curvedFront offset is the slope s(rho)
    # integrated, with c0 s - s' - 2 s/rho = 1 the unbounded radial equation
    c0 = 30.0
    rho = np.linspace(2/c0*1.05, 5/c0, 9)
    off = ref["seed_offset"]
    eps = 1e-7
    slope_rho = (off(rho + eps, c0, "curvedFront") - off(rho - eps, c0, "curvedFront"))/(2*eps)
    _check(np.max(np.abs(slope_rho - (1/c0 + 2/(c0**2*rho) + 2/(c0**3*rho**2)))) < 1e-6,
           "curvedFront offset is not the integral of the radial slope")
    _check(float(off(2/c0, c0, "curvedFront")) == 0.0 and float(off(0.0, c0, "curvedFront")) == 0.0,
           "curvedFront offset is not zero inside the critical radius")
    _check(np.allclose(off(rho, c0, "front"), rho/c0, rtol=0, atol=1e-15)
           and not np.any(off(rho, c0, "flat")), "front or flat offset")


# ---- emission and command line ------------------------------------------------

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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    grp = ap.add_mutually_exclusive_group()
    grp.add_argument("--emit", action="store_true", help="write reference/radial_eikonal.py")
    grp.add_argument("--self-check", action="store_true", help="run the self-checks")
    grp.add_argument("--report", action="store_true",
                     help="print the arrival time at the wall and the slope on the ball "
                          "(the default)")
    ap.add_argument("--PeL", type=float, default=46.0, help="Pe_L of --report (default 46)")
    ap.add_argument("--ballRatio", type=float, default=0.25,
                    help="R/L of --report (default 0.25)")
    args = ap.parse_args(argv)
    if args.emit or args.self_check:
        try:
            run_self_checks()
        except AssertionError as err:
            print("self-check failed: %s" % err)
            return 1
        if args.self_check:
            print("self-checks passed")
            return 0
        for rel in write_outputs():
            print("wrote %s" % rel)
        return 0
    a = args.ballRatio
    if args.PeL*(a + 1) > CLOSED_FORM_LIMIT:
        ap.error("--PeL times (R/L + 1) above %g overflows the closed form" % CLOSED_FORM_LIMIT)
    print("Pe_L %g, R/L %g: scaled arrival at the wall %.6f, start slope %.4f" % (
        args.PeL, a, rad_arrival(args.PeL, a), slope(np.array(a), args.PeL, a)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
