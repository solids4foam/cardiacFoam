import importlib.util
import os
import re

import numpy as np
import pytest
import sympy as sp

import cartesian_eikonal_mms as eik
import eikonal_ecg_mms as ecg

DIMS = (1, 2, 3)
SIGMAS = {"rotated": eik.SIGMA_ROTATED, "axis": eik.SIGMA_AXIS}
ELECTRODE_FILE = os.path.join(
    ecg.HERE, "..", "eikonalECG", "constant", "electroProperties")
BUILTIN_ELECTRODES = [(-0.5, 0.5, 0.5), (1.5, 0.5, 0.5), (1.2, 0.23, 0.61),
                      (1.35, 0.74, 0.28), (1.55, 0.41, 0.83)]


def tutorial_electrodes():
    """The tutorial electrode set when present, else the built-in set."""
    if not os.path.exists(ELECTRODE_FILE):
        return list(BUILTIN_ELECTRODES)
    pat = re.compile(r"^\s*([ER]\d+)\s+\(([^)]*)\);", re.M)
    with open(ELECTRODE_FILE) as fh:
        found = pat.findall(fh.read())
    return [tuple(float(v) for v in xyz.split()) for _, xyz in found]


def _reference():
    path = os.path.join(ecg.HERE, "reference", "eikonal_ecg.py")
    spec = importlib.util.spec_from_file_location("eikecg_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _sig(name):
    return np.array(eik.sigma_matrix(SIGMAS[name]).tolist(), dtype=float)


@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("template", ecg.TEMPLATES)
def test_integrand_matches_cxx_form(d, template):
    S = eik.sigma_matrix(eik.SIGMA_ROTATED)
    e = ecg.ELECTRODE_SYMBOLIC
    a = ecg.integrand(d, S, e, template)
    b = ecg.cxx_integrand(d, S, e, template)
    rng = np.random.default_rng(1)
    syms = list(ecg.AXES) + list(ecg.ELECTRODE_SYMBOLS) + [ecg.t, ecg.DURATION, ecg.RISE]
    for _ in range(5):
        p = list(rng.uniform(-0.5, 1.5, 6)) + [rng.uniform(0, 40), 64.0, 4.0]
        sub = dict(zip(syms, p))
        assert float(a.subs(sub)) == pytest.approx(float(b.subs(sub)), rel=1e-13)


def test_template_derivative_is_derivative_of_template():
    # U' of the single event is the s-derivative of 0.25 (1 + a)(1 - b)
    w = ecg.RISE/(2*sp.atanh(sp.Rational(4, 5)))
    U = sp.Rational(1, 4)*(1 + sp.tanh(ecg.s/w))*(1 - sp.tanh((ecg.s - ecg.DURATION)/w))
    diff = sp.diff(U, ecg.s) - ecg.template_derivative("single_event")
    for sv in (-3.0, 1.0, 30.0, 66.0):
        assert abs(float(diff.subs({ecg.s: sv, ecg.DURATION: 64, ecg.RISE: 4}))) < 1e-13


def test_default_template_parameters():
    assert ecg.DEFAULT_DURATION == 64 and ecg.DEFAULT_RISE == 4


@pytest.mark.parametrize("template", ecg.TEMPLATES)
def test_eikonal_ecg_quadrature_converged(template):
    S = eik.sigma_matrix(eik.SIGMA_ROTATED)
    electrodes = tutorial_electrodes()
    for e in electrodes[:5] + electrodes[5::40]:
        lo = ecg.reference(3, S, e, 31.0, 48, template)
        hi = ecg.reference(3, S, e, 31.0, 96, template)
        assert abs(lo - hi) <= 1e-9*max(1.0, abs(hi)), (e, lo, hi)


@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("template", ecg.TEMPLATES)
def test_emitted_quadrature_matches_generator(d, template):
    ref = _reference()
    for name in SIGMAS:
        S = eik.sigma_matrix(SIGMAS[name])
        for e in BUILTIN_ELECTRODES[:3]:
            want = ecg.reference(d, S, e, 28.0, 24, template)
            got = ref.spatial_reference(d, template, _sig(name), e, 28.0, 24)
            assert got == pytest.approx(want, rel=1e-13, abs=1e-13), (d, template, e)


@pytest.mark.parametrize("d", DIMS)
def test_coincident_point_returns_zero(d):
    ref = _reference()
    S = _sig("axis")
    for template in ecg.TEMPLATES:
        v = ref.spatial_reference(d, template, S, (0.5, 0.5, 0.5), 10.0, 1)
        assert v == 0.0 and np.isfinite(v)
        with np.errstate(all="ignore"):
            pt = getattr(ref, "integrand_d%d_%s" % (d, template))(
                0.5, 0.5, 0.5, 10.0, 0.5, 0.5, 0.5, *S[[0, 0, 0, 1, 1, 2],
                                                     [0, 1, 2, 1, 2, 2]], 64.0, 4.0)
        assert pt == 0.0


def test_d1_guard_is_on_abs_dx():
    ref = _reference()
    sig = (1.0, 0, 0, 1.0, 0, 1.0)
    with np.errstate(all="ignore"):
        f = ref.integrand_d1_sine
        assert f(1e-301, 0, 0, 0.3, 0, 0, 0, *sig, 64.0, 4.0) == 0.0
        assert f(1e-299, 0, 0, 0.3, 0, 0, 0, *sig, 64.0, 4.0) != 0.0


def test_self_checks_and_cli():
    assert ecg.main(["--check"]) == 0


def test_emitted_files_up_to_date():
    assert ecg.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    ecg.write_outputs(str(tmp_path))
    for rel, text in ecg.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def test_emitted_module_is_numpy_only():
    text = open(os.path.join(ecg.HERE, "reference", "eikonal_ecg.py")).read()
    assert "sympy" not in text and "mms_common" not in text


def _parity():
    try:
        from parity.test_parity import compare, parity_available
    except ImportError:
        pytest.skip("parity harness absent (cardiacFoam-only)")
    ok, why = parity_available()
    if not ok:
        pytest.skip(why)
    return compare


def _parity_points():
    pts = tutorial_electrodes()[:5] + tutorial_electrodes()[5::40]
    rng = np.random.default_rng(5)
    return np.column_stack([np.array(pts), rng.uniform(-6, 92, len(pts))])


def _kernel_scale(d, sigma, electrode, order):
    """Sum of |weight * (r . sigma grad psi) / |r|^3|: the size of the terms of the
    quadrature sum when U' = 1."""
    import cartesian_eikonal_mms as e
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = 0.5*(nodes + 1.0), 0.5*weights
    grid = np.meshgrid(*([nodes]*d), indexing="ij")
    w = np.ones_like(grid[0])
    for i in range(d):
        shape = [1]*d
        shape[i] = order
        w = w*weights.reshape(shape)
    c = list(grid) + [np.zeros_like(grid[0])]*(3 - d)
    k = [1.0 if m < d else 0.0 for m in range(3)]
    psi = np.exp(sum(k[m]*c[m] for m in range(3)))
    sg = sigma @ np.array(k)
    r = [c[m] - electrode[m] for m in range(d)]
    num = sum(r[m]*sg[m]*psi for m in range(d))
    rn = np.sqrt(sum(q**2 for q in r))
    with np.errstate(all="ignore"):
        val = np.where(rn > 1e-150, num/rn**3, 0.0)
    return float(np.sum(np.abs(w*val)))


@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("template", ecg.TEMPLATES)
@pytest.mark.parametrize("name", sorted(SIGMAS))
def test_parity_eikonal_ecg(d, template, name):
    """The quadrature sum cancels and U' is evaluated through tanh/cos in a
    different order than the C++, so a value passes when
    |diff| <= 1e-12 |C++ value| + 1e-14 (sum of |terms| for U' = 1), the rounding
    noise of the terms."""
    _parity()
    from parity.test_parity import run_app
    ref = _reference()
    S = _sig(name)
    pts = _parity_points()
    cpp_template = "sine" if template == "sine" else "singleEvent"
    cpp = run_app("eikonalECG", d, pts, order=ecg.TUTORIAL_QUADRATURE_ORDER,
                  template=cpp_template, conductivity=name, duration=64, rise=4)[:, 0]
    got = np.array([ref.spatial_reference(d, template, S, r[:3], r[3],
                                          ecg.TUTORIAL_QUADRATURE_ORDER) for r in pts])
    floor = 1e-14*np.array([_kernel_scale(d, S, r[:3], ecg.TUTORIAL_QUADRATURE_ORDER)
                            for r in pts])
    assert (np.abs(got - cpp) <= 1e-12*np.abs(cpp) + floor).all(), np.abs(got - cpp).max()
