import importlib.util
import os
import re

import numpy as np
import pytest
import sympy as sp

import fda_monodomain_mms as fda
import pseudo_ecg_mms as ecg
import rotated_monodomain_mms as rot

DIMS = (1, 2, 3)
ELECTRODE_FILE = os.path.join(
    ecg.HERE, "..", "monodomainPseudoECG", "constant", "electroProperties")


BUILTIN_ELECTRODES = [(-0.5, 0.5, 0.5), (1.5, 0.5, 0.5), (1.2, 0.23, 0.61),
                      (0.3, 0.6, 1.4), (0.55, 0.45, 0.5)]


def tutorial_electrodes():
    """The tutorial electrode set when present, else the built-in set."""
    if not os.path.exists(ELECTRODE_FILE):
        return list(BUILTIN_ELECTRODES)
    pat = re.compile(r"^\s*([ER]\d+)\s+\(([^)]*)\);", re.M)
    with open(ELECTRODE_FILE) as fh:
        found = pat.findall(fh.read())
    return [tuple(float(v) for v in xyz.split()) for _, xyz in found]


def _reference():
    path = os.path.join(ecg.HERE, "reference", "pseudo_ecg.py")
    spec = importlib.util.spec_from_file_location("ecg_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("sigma", ["fda", "rotated"])
def test_integrand_matches_cxx_form(d, sigma):
    S = ecg.SIGMA_FDA if sigma == "fda" else rot.K_CASE
    e = ecg.ELECTRODE_SYMBOLIC
    a, b = ecg.integrand(d, S, e), ecg.cxx_integrand(d, S, e)
    rng = np.random.default_rng(1)
    for _ in range(5):
        p = rng.uniform(-0.5, 1.5, 6)
        sub = dict(zip(list(ecg.AXES) + list(ecg.ELECTRODE_SYMBOLS), p))
        assert float(a.subs(sub)) == pytest.approx(float(b.subs(sub)), rel=1e-13)


def test_integrand_uses_the_fda_gradient():
    g = ecg.grad_F(2)
    F = fda.fields(2)["F"]
    assert list(g) == [sp.diff(F, fda.x), sp.diff(F, fda.y), 0]


def test_electrode_arrays_agree_with_sympy_integrand():
    ref = _reference()
    rng = np.random.default_rng(2)
    p = rng.uniform(0, 1, (3, 5))
    e = (1.2, 0.23, 0.61)
    for d in DIMS:
        g = sp.lambdify(ecg.AXES, ecg.integrand(d, ecg.SIGMA_FDA, e), "numpy")
        sig = [float(ecg.SIGMA_FDA[i, j]) for i, j in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))]
        got = getattr(ref, "integrand_d%d" % d)(*p, *e, *sig)
        np.testing.assert_allclose(got, g(*p), rtol=1e-13)


@pytest.mark.parametrize("sigma", ["fda", "rotated"])
def test_quadrature_converged(sigma):
    electrodes = tutorial_electrodes()
    assert len(electrodes) in (5, 161)
    S = ecg.SIGMA_FDA if sigma == "fda" else rot.K_CASE
    F = None if sigma == "fda" else ecg._rotated_F()
    for e in electrodes[:5] + electrodes[5::12]:
        lo = ecg.spatial_reference(3, S, e, 48, F)
        hi = ecg.spatial_reference(3, S, e, 96, F)
        assert abs(lo - hi) <= 1e-10*max(1.0, abs(hi)), (e, lo, hi)


@pytest.mark.parametrize("d", DIMS)
def test_coincident_point_returns_zero(d):
    # order 1 has the single node 1/2: an electrode there sits on it
    S = ecg.SIGMA_FDA
    assert ecg.spatial_reference(d, S, (0.5, 0.5, 0.5), 1) == 0.0
    ref = _reference()
    case = "d%d" % d
    sig = [float(S[i, j]) for i, j in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))]
    v = ref.spatial_reference(case, np.array(S.tolist(), dtype=float), (0.5, 0.5, 0.5), 1)
    assert v == 0.0 and np.isfinite(v)
    with np.errstate(all="ignore"):
        pt = getattr(ref, "integrand_" + case)(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, *sig)
    assert pt == 0.0


@pytest.mark.parametrize("case,d", [("d1", 1), ("d2", 2), ("d3", 3), ("rotated", 3)])
def test_emitted_quadrature_matches_generator(case, d):
    ref = _reference()
    S = ecg.SIGMA_FDA if case != "rotated" else rot.K_CASE
    F = ecg._rotated_F() if case == "rotated" else None
    arr = np.array(S.tolist(), dtype=float)
    for e in BUILTIN_ELECTRODES:
        want = ecg.spatial_reference(d, S, e, 48, F)
        got = ref.spatial_reference(case, arr, e, 48)
        assert got == pytest.approx(want, rel=1e-13, abs=1e-13), (case, e)


def test_d1_guard_is_on_abs_dx():
    # C++ d = 1: |dx| <= 1e-300 returns 0 (dx^2 would test |dx| <= 1e-150)
    ref = _reference()
    sig = (1.0, 0, 0, 1.0, 0, 1.0)
    with np.errstate(all="ignore"):
        assert ref.integrand_d1(1e-301, 0, 0, 0.0, 0, 0, *sig) == 0.0
        assert ref.integrand_d1(1e-299, 0, 0, 0.0, 0, 0, *sig) != 0.0
        assert ref.integrand_d1(1e-200, 0, 0, 0.0, 0, 0, *sig) != 0.0
    g = ecg.integrand(1, ecg.SIGMA_FDA, ecg.ELECTRODE_SYMBOLIC)
    assert g.atoms(sp.Abs)


def test_reference_scales_with_time_factor():
    e = (1.4, 0.3, 0.2)
    a = ecg.reference(3, ecg.SIGMA_FDA, e, 0.0, 12)
    b = ecg.reference(3, ecg.SIGMA_FDA, e, 0.44, 12)
    assert b == pytest.approx(a*np.sqrt(1.44), rel=1e-14)


def test_self_checks_and_cli():
    assert ecg.main(["--check"]) == 0


def test_emitted_files_up_to_date():
    assert ecg.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    ecg.write_outputs(str(tmp_path))
    for rel, text in ecg.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def test_emitted_module_is_numpy_only():
    text = open(os.path.join(ecg.HERE, "reference", "pseudo_ecg.py")).read()
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
    pts = tutorial_electrodes()[:5] + tutorial_electrodes()[5::31]
    rng = np.random.default_rng(5)
    return np.column_stack([np.array(pts), rng.uniform(0, 0.1, len(pts))])


@pytest.mark.parametrize("d", DIMS)
def test_parity_pseudo_ecg(d):
    compare = _parity()
    ref = _reference()
    S = np.array(ecg.SIGMA_FDA.tolist(), dtype=float)

    def numpy_fn(p):
        rows = []
        for row in p:
            s = ref.spatial_reference("d%d" % d, S, row[:3], ecg.TUTORIAL_QUADRATURE_ORDER)
            rows.append((s, ref.time_factor(row[3])*s))
        return np.array(rows)

    compare("pseudoECG", d, numpy_fn, _parity_points(), rtol=1e-12, order=96)


def test_parity_pseudo_ecg_rotated():
    compare = _parity()
    ref = _reference()
    K = np.array(rot.K_CASE.tolist(), dtype=float)

    def numpy_fn(p):
        rows = []
        for row in p:
            s = ref.spatial_reference("rotated", K, row[:3], ecg.TUTORIAL_QUADRATURE_ORDER)
            rows.append((s, ref.time_factor(row[3])*s))
        return np.array(rows)

    compare("pseudoECGRotated", 3, numpy_fn, _parity_points(), rtol=1e-12, order=96)
