import importlib.util
import math
import os

import numpy as np
import pytest

import coupled_1d3d_mms as cp
from mms_common import quadrature


def _reference():
    path = os.path.join(cp.HERE, "reference", "coupled_1d3d.py")
    spec = importlib.util.spec_from_file_location("coupled_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_terminal_current_exact():
    for t in (0.0, 0.07, 0.1):
        a, b = cp.terminal_current_exact(t)
        assert a == pytest.approx(1.5*np.sqrt(1 + t), rel=1e-15)
        assert b == pytest.approx(-1.5*np.sqrt(1 + t), rel=1e-15)
    a, b = cp.terminal_current_exact(0.3, R_pvj=2.0)
    assert a == pytest.approx(0.75*np.sqrt(1.3), rel=1e-15)


def test_terminal_point_values_give_exact_current():
    t = 0.07
    for k in (0, 1):
        vm = float(cp._vm_numpy()(*cp.TERMINALS[k], t))
        assert cp.terminal_vm(k, t) - vm == pytest.approx(cp.terminal_current_exact(t)[k], rel=1e-14)


def test_continuum_sampling_converges():
    for k in (0, 1):
        a = cp.continuum_sampling(k, 0.07, n=32)
        b = cp.continuum_sampling(k, 0.07, n=64)
        assert abs(a - b) < 1e-9


def test_continuum_current_within_second_order_of_point_value():
    # |Vbar - V_m(x_k)| <= 7 pi^2 r^2 sqrt(1+t): dV_m/dx = 0 and the ball is
    # symmetric in y and z, so the first-order term vanishes, and
    # |s^T H s| <= 14 pi^2 |s|^2 with the weighted mean of |s|^2 at most r^2.
    t = 0.07
    errs = {}
    for r in (0.11, 0.055):
        for k in (0, 1):
            point = cp.terminal_current_exact(t)[k]
            err = abs(cp.continuum_current(k, t, radius=r) - point)
            assert err <= 7*np.pi**2*r**2*np.sqrt(1 + t)
            errs[(r, k)] = err
    for k in (0, 1):
        assert 3.0 < errs[(0.11, k)]/errs[(0.055, k)] < 5.0  # second order


def test_kernel_mass_is_half_ball():
    _, mass = quadrature.half_ball_kernel_integral(
        lambda p: np.ones(p.shape[1:]), cp.TERMINALS[1], cp.RADIUS, cp.CUBE_LO, cp.CUBE_HI, 32)
    assert mass == pytest.approx(cp.kernel_mass(), rel=1e-12)


def _deposition_integral(k, t):
    """int of the deposition over the half ball at terminal k, by a tensor
    Gauss rule in spherical coordinates (r, mu, phi) about the terminal with
    the polar axis along x.  The kernel is linear in r, so the rule is exact
    in r and smooth in the angles; the cube cuts the ball only at the face
    through the terminal (y and z margins 1/6 and 1/3 exceed the radius)."""
    n = 16
    u, wu = np.polynomial.legendre.leggauss(n)
    r = 0.5*cp.RADIUS*(u + 1.0)
    wr = 0.5*cp.RADIUS*wu
    lo, hi = (0.0, 1.0) if k == 0 else (-1.0, 0.0)
    mu = 0.5*(hi - lo)*u + 0.5*(hi + lo)
    wm = 0.5*(hi - lo)*wu
    phi = np.pi*(u + 1.0)
    wp = np.pi*wu
    R, MU, PHI = np.meshgrid(r, mu, phi, indexing="ij")
    sin = np.sqrt(1.0 - MU**2)
    pts = np.stack([MU*R, sin*np.cos(PHI)*R, sin*np.sin(PHI)*R], axis=-1) + np.array(cp.TERMINALS[k])
    vals = cp.terminal_deposition(k, pts.reshape(-1, 3), t).reshape(R.shape)
    w = wr[:, None, None]*wm[None, :, None]*wp[None, None, :]*R**2
    return float(np.sum(w*vals))


def test_junction_source_integrates_to_current_per_terminal():
    # Each terminal's deposition integrates to its own current (the sum of the
    # two is trivially zero and would hide a wrong split).
    t = 0.07
    for k in (0, 1):
        assert _deposition_integral(k, t) == pytest.approx(
            cp.continuum_current(k, t), rel=1e-11)
    assert cp.continuum_current(0, t)*cp.continuum_current(1, t) < 0
    # a wrongly normalised kernel (full-ball mass) would halve the integral
    assert _deposition_integral(0, t) != pytest.approx(0.5*cp.continuum_current(0, t), rel=1e-3)


def test_junction_source_is_sum_of_depositions():
    t = 0.05
    rng = np.random.default_rng(3)
    pts = rng.uniform(0, 1, (60, 3))
    pts[:15, 0] *= 0.15          # inside terminal 0's reach
    pts[15:30, 0] = 1 - pts[15:30, 0]*0.15   # inside terminal 1's reach
    pts[:30, 1] = 1/6 + rng.uniform(-0.1, 0.1, 30)
    pts[:30, 2] = 1/3 + rng.uniform(-0.1, 0.1, 30)
    currents = [cp.continuum_current(k, t) for k in (0, 1)]
    mass = math.pi*cp.RADIUS**3/6.0
    want = np.empty(len(pts))
    for i, p in enumerate(pts):
        total = 0.0
        for k in (0, 1):
            c = cp.TERMINALS[k]
            dist = math.sqrt(sum((p[j] - c[j])**2 for j in range(3)))
            total += currents[k]*max(0.0, 1.0 - dist/cp.RADIUS)/mass
        want[i] = total
    assert np.count_nonzero(want) >= 10
    np.testing.assert_allclose(cp.continuum_junction_source(pts, t), want, rtol=1e-14, atol=0)
    assert np.all(cp.continuum_junction_source([[0.5, 0.5, 0.5]], t) == 0)


def test_self_checks_and_cli():
    assert cp.main(["--self-check"]) == 0


def test_emitted_files_up_to_date():
    assert cp.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    cp.write_outputs(str(tmp_path))
    for rel, text in cp.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def test_emitted_reference_is_numpy_only():
    text = open(os.path.join(cp.HERE, "reference", "coupled_1d3d.py")).read()
    assert "sympy" not in text and "mms_common" not in text
    assert "import numpy\n" in text


def test_emitted_quadrature_equals_mms_common():
    ref = _reference()
    for centre in ((0.0, 1/6, 1/3), (1.0, 1/6, 1/3)):
        def f(p):
            return np.cos(np.pi*p[0])*np.cos(2*np.pi*p[1])*np.cos(3*np.pi*p[2])
        got = ref.half_ball_kernel_integral(f, centre, 0.11, np.zeros(3), np.ones(3), 24)
        want = quadrature.half_ball_kernel_integral(f, centre, 0.11, np.zeros(3), np.ones(3), 24)
        np.testing.assert_allclose(got, want, rtol=1e-13, atol=0)


def test_emitted_functions_match_module():
    ref = _reference()
    t = 0.04
    for k in (0, 1):
        assert ref.continuum_sampling(k, t) == pytest.approx(cp.continuum_sampling(k, t), rel=1e-13)
        assert ref.continuum_current(k, t) == pytest.approx(cp.continuum_current(k, t), rel=1e-12)
    pts = np.random.default_rng(1).uniform(0, 1, (20, 3))
    pts[:, 0] *= 0.1
    np.testing.assert_allclose(ref.continuum_junction_source(pts, t),
                               cp.continuum_junction_source(pts, t), rtol=1e-12, atol=1e-14)
    assert ref.terminal_current_exact(t) == cp.terminal_current_exact(t)


def test_parity_coupled_terminal_current():
    try:
        from parity.test_parity import compare, parity_available, sample_points
    except ImportError:
        pytest.skip("parity harness absent (cardiacFoam-only)")
    ok, why = parity_available()
    if not ok:
        pytest.skip(why)
    ref = _reference()

    def numpy_fn(p):
        rows = []
        for t in p[:, 3]:
            row = []
            for k in (0, 1):
                x, y, z = ref.TERMINALS[k]
                v1, v3 = ref.V1D(x, t), ref.Vm3D(x, y, z, t)
                row += [v1, v3, ref.terminal_current_exact(t)[k]]
            rows.append(row)
        return np.array(rows)

    # The graph file stores the terminal y and z to 12 digits (0.166666666667,
    # 0.333333333333) while this reference uses 1/6 and 1/3, so the comparison
    # is at rel 1e-10 rather than 1e-12.
    compare("coupled", 3, numpy_fn, sample_points(), rtol=1e-10)
