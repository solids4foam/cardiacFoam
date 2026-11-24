"""The radial solution satisfies its equation, its wall condition and its scaling."""

import ast
import importlib.util
import os
import re

import numpy as np
import pytest

import radial_eikonal_mms as rad

HEADER = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "src", "verificationModels",
    "eikonalVerification", "manufacturedRadialEikonalReference.H",
)


@pytest.mark.parametrize("Pe", [5.0, 21.0, 46.0, 94.0])
@pytest.mark.parametrize("a", [0.1, 0.25, 0.5])
def test_equation_residual(Pe, a):
    s = np.linspace(a, a + 1, 20001)[1:-1]
    ds = 1e-5
    w = rad.slope(s, Pe, a)
    dw = (rad.slope(s + ds, Pe, a) - rad.slope(s - ds, Pe, a))/(2*ds)
    assert np.abs((dw + 2*w/s)/Pe - (w - 1)).max() < 1e-6


@pytest.mark.parametrize("Pe", [5.0, 46.0, 94.0])
def test_time_integrates_slope(Pe):
    a = 0.25
    s = np.linspace(a, a + 1, 400001)
    w = rad.slope(s, Pe, a)
    quadrature = np.sum(0.5*(w[1:] + w[:-1])*np.diff(s))
    assert abs(quadrature - rad.scaled_time(np.array(a + 1.0), Pe, a)) < 1e-8


def test_boundary_values():
    Pe, a = 46.0, 0.25
    assert abs(rad.slope(np.array(a + 1.0), Pe, a)) < 1e-14
    assert abs(rad.scaled_time(np.array(a), Pe, a)) < 1e-14


def test_front_limit():
    # Away from the ball and the wall, V dpsi/dr tends to 1 + 2/(Pe s), the curvature slowing
    Pe, a, s = 500.0, 0.25, 0.8
    assert abs(rad.slope(np.array(s), Pe, a) - (1 + 2/(Pe*s))) < 1e-4


def test_header_samples():
    # The self-test samples of the C++ reference come from this solution
    with open(HEADER) as f:
        text = f.read()
    rows = re.findall(r"\{([0-9.]+), ([0-9.]+), ([0-9.]+), ([-0-9.e]+), ([-0-9.e]+)\},", text)
    assert len(rows) == 7
    for row in rows:
        Pe, a, s, w, time = map(float, row)
        assert abs(rad.slope(np.array(s), Pe, a) - w) < 1e-13
        assert abs(rad.scaled_time(np.array(s), Pe, a) - time) < 1e-11


def test_seed_offset_slope_solves_unbounded_radial_equation():
    # pvjMapper::setSeedOffsets integrates s = 1/c0 + 2/(c0^2 rho) + 2/(c0^3 rho^2)
    # of c0 s - s' - 2 s/rho = 1 in the metric distance, from rho = 2/c0
    c0 = 30.0
    rho = np.linspace(2/c0, 5/c0, 7)
    s = 1/c0 + 2/(c0**2*rho) + 2/(c0**3*rho**2)
    ds = -2/(c0**2*rho**2) - 4/(c0**3*rho**3)
    assert np.max(np.abs(c0*s - ds - 2*s/rho - 1)) < 1e-12
    rhoC = 2/c0
    offset = lambda r: (r - rhoC)/c0 + 2*np.log(r/rhoC)/c0**2 - 2*(1/r - 1/rhoC)/c0**3
    assert abs(offset(rhoC)) < 1e-15
    eps = 1e-7
    assert np.max(np.abs((offset(rho + eps) - offset(rho - eps))/(2*eps) - s)) < 1e-7


@pytest.mark.parametrize("Pe,a", [(2.0, 0.05), (5.0, 0.02), (21.0, 0.25), (46.0, 0.25),
                                  (94.0, 0.1), (94.0, 0.5)])
def test_parity_radial(Pe, a):
    """The quadrature of the C++ header equals the closed form with the
    exponential integral to 1e-9."""
    try:
        from parity.test_parity import compare, parity_available
    except ImportError:
        pytest.skip("parity harness absent (cardiacFoam-only)")
    ok, why = parity_available()
    if not ok:
        pytest.skip(why)
    s = np.linspace(1.01*a, a + 1.0, 30)
    pts = np.column_stack([s, np.zeros_like(s), np.zeros_like(s), np.zeros_like(s)])

    def numpy_fn(p):
        return np.column_stack([rad.slope(p[:, 0], Pe, a), rad.scaled_time(p[:, 0], Pe, a)])

    compare("radial", 1, numpy_fn, pts, rtol=1e-9, Pe=Pe, a=a)


# ---- the emitted NumPy reference ---------------------------------------------

def _reference():
    path = os.path.join(rad.HERE, "reference", "radial_eikonal.py")
    spec = importlib.util.spec_from_file_location("radial_eikonal_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_self_checks_and_cli(capsys):
    assert rad.main(["--self-check"]) == 0
    assert rad.main([]) == 0
    assert "scaled arrival at the wall 1.049665" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        rad.main(["--PeL", "1000"])


def test_emitted_files_up_to_date():
    assert rad.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    rad.write_outputs(str(tmp_path))
    for rel, text in rad.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def test_emitted_reference_is_numpy_only():
    text = open(os.path.join(rad.HERE, "reference", "radial_eikonal.py")).read()
    for banned in ("sympy", "scipy", "mms_common", "expi"):
        assert banned not in text
    assert "import numpy\n" in text
    imported = [alias.name for node in ast.walk(ast.parse(text))
                for alias in getattr(node, "names", []) if isinstance(node, (ast.Import, ast.ImportFrom))]
    assert imported == ["numpy"] and not any(isinstance(n, ast.ImportFrom) for n in ast.walk(ast.parse(text)))


@pytest.mark.parametrize("Pe", [2.0, 5.0, 21.0, 46.0, 94.0])
@pytest.mark.parametrize("a", [0.02, 0.1, 0.25, 0.5])
def test_reference_matches_closed_form(Pe, a):
    ref = _reference()
    s = np.linspace(a, a + 1, 101)[1:]
    np.testing.assert_allclose(ref.slope(s, Pe, a), rad.slope(s, Pe, a), rtol=0, atol=1e-13)
    # the panel quadrature of the C++ reference is accurate to 1.5e-10 at Pe 2, a 0.02
    np.testing.assert_allclose(ref.scaled_time(s, Pe, a), rad.scaled_time(s, Pe, a),
                               rtol=0, atol=3e-10)


@pytest.mark.parametrize("Pe,a", [(1000.0, 0.25), (3000.0, 0.1)])
def test_reference_beyond_the_closed_form(Pe, a):
    # exp(Pe s) overflows the closed form here; the quadrature does not
    assert Pe*(a + 1) > rad.CLOSED_FORM_LIMIT
    ref = _reference()
    s = np.linspace(a, a + 1, 400001)
    w = ref.slope(s, Pe, a)
    trapezoid = np.sum(0.5*(w[1:] + w[:-1])*np.diff(s))
    assert abs(trapezoid - float(ref.scaled_time(np.array(a + 1.0), Pe, a))) < 1e-6
    t = ref.scaled_time(s[::1000], Pe, a)
    assert np.all(np.isfinite(t)) and np.all(np.diff(t) > 0)


@pytest.mark.skipif(not os.path.exists(HEADER), reason="cardiacFoam source tree absent")
def test_reference_header_samples():
    ref = _reference()
    with open(HEADER) as f:
        text = f.read()
    rows = re.findall(r"\{([0-9.]+), ([0-9.]+), ([0-9.]+), ([-0-9.e]+), ([-0-9.e]+)\},", text)
    assert len(rows) == 7
    for row in rows:
        Pe, a, s, w, time = map(float, row)
        assert abs(ref.slope(np.array(s), Pe, a) - w) < 1e-13
        assert abs(ref.scaled_time(np.array(s), Pe, a) - time) < 1e-10


def test_reference_scalar_and_array_arguments_agree():
    ref = _reference()
    s = np.array([0.3, 0.7, 1.2])
    got = ref.scaled_time(s, 46.0, 0.25)
    for i, v in enumerate(s):
        assert float(ref.scaled_time(v, 46.0, 0.25)) == pytest.approx(got[i], rel=1e-15)
    assert ref.scaled_time(s.reshape(3, 1), 46.0, 0.25).shape == (3, 1)


def test_parameters():
    ref = _reference()
    V, Pe = ref.parameters(c0=5.0, m=0.25, L=0.1)
    assert V == pytest.approx(2.5, rel=1e-15) and Pe == pytest.approx(1.0, rel=1e-15)


def test_seed_profiles_differ_in_the_dirichlet_data():
    ref = _reference()
    Pe, a, L, V, Tk = 46.0, 0.25, 0.1, 0.8, 0.5
    s = np.array([0.3, 0.7, 1.2])
    flat = ref.exact_time(s, Pe, a, L, V, Tk, "flat")
    assert np.all(ref.exact_time(s, Pe, a, L, V, Tk, "exact") == flat)
    assert float(ref.exact_time(a, Pe, a, L, V, Tk, "flat")) == pytest.approx(Tk, abs=1e-15)
    front = ref.exact_time(s, Pe, a, L, V, Tk, "front")
    assert float(ref.exact_time(a, Pe, a, L, V, Tk, "front")) == pytest.approx(Tk + a*L/V, abs=1e-15)
    np.testing.assert_allclose(front - flat, a*L/V, rtol=0, atol=1e-14)
    core = 2/Pe
    curved = ref.exact_time(s, Pe, a, L, V, Tk, "curvedFront")
    assert float(ref.exact_time(core, Pe, a, L, V, Tk, "curvedFront")) == pytest.approx(Tk, abs=1e-15)
    assert float(ref.exact_time(0.5*core, Pe, a, L, V, Tk, "curvedFront")) == pytest.approx(Tk, abs=1e-15)
    np.testing.assert_allclose(curved - flat, -L/V*float(ref.scaled_time(core, Pe, a)),
                               rtol=0, atol=1e-14)
    with pytest.raises(ValueError, match="2/Pe"):
        ref.exact_time(s, 5.0, 0.1, L, V, Tk, "curvedFront")
    with pytest.raises(ValueError, match="profile"):
        ref.exact_time(s, Pe, a, L, V, Tk, "point")


def test_seed_profile_arrival_dependence_on_ball_radius():
    # README of the tutorial coupled1D3DSphere: the time at a fixed point depends on R
    # through R/L for flat and (almost) not at all for curvedFront.  At fixed L and Pe_L,
    # with R/L from 0.1 to 0.5, the exact time at r = L and Pe_L = 46 changes by 0.48
    # (flat), 0.08 (front) and less than 0.01 (curvedFront) L/V.
    ref = _reference()
    Pe, L, V = 46.0, 1.0, 1.0
    change = {}
    for profile in ("flat", "front", "curvedFront"):
        t = [float(ref.exact_time(1.0, Pe, a, L, V, 0.0, profile)) for a in (0.1, 0.5)]
        change[profile] = abs(t[1] - t[0])
    assert change["flat"] == pytest.approx(0.48, abs=0.01)
    assert change["front"] == pytest.approx(0.08, abs=0.01)
    assert change["curvedFront"] < 0.01


def test_anisotropic_field_is_isotropic_in_the_stretched_coordinate():
    ref = _reference()
    Pe, a, L, V = 21.0, 0.25, 0.2, 1.7
    axes, angle = (1.4, 0.7, 1.0), 0.9
    xi = np.array([[0.1, -0.15, 0.2], [0.05, 0.3, -0.2]])
    x = (xi*np.asarray(axes)) @ ref.rotation_about_diagonal(angle).T
    iso = ref.exact_field(xi, Pe, a, L, V, 0.3, "front")
    np.testing.assert_allclose(ref.exact_field(x, Pe, a, L, V, 0.3, "front", axes, angle), iso,
                               rtol=1e-14)
    Q = ref.rotation_about_diagonal(angle)
    np.testing.assert_allclose(Q @ Q.T, np.eye(3), atol=1e-15)
    assert np.linalg.det(Q) == pytest.approx(1.0, rel=1e-14)
    assert np.allclose(Q @ np.ones(3), np.ones(3), atol=1e-14)    # the axis is fixed
    # inside the ball the value is the surface value, and the gradient is zero
    inner = ref.exact_field(np.array([[0.01, 0.0, 0.0]]), Pe, a, L, V, 0.3, "flat")
    assert float(inner[0]) == pytest.approx(0.3, abs=1e-15)
    assert not np.any(ref.exact_gradient(np.array([[0.01, 0.0, 0.0]]), Pe, a, L, V))


def test_seed_offsets_follow_the_mapper():
    # pvjMapper::setSeedOffsets, written out per cell as the C++ does
    ref = _reference()
    c0 = 30.0
    rng = np.random.default_rng(5)
    A = rng.normal(size=(3, 3))
    M = A @ A.T + 3*np.eye(3)
    centres = rng.uniform(-0.2, 0.2, (40, 3))
    terminal = np.array([0.01, -0.02, 0.03])
    for profile in ("flat", "front", "curvedFront"):
        want = []
        for c in centres:
            d = c - terminal
            rho = np.sqrt(max(d @ np.linalg.inv(M) @ d, 0.0))
            rho_c = 2/c0
            if profile == "curvedFront":
                off = 0.0
                if rho > rho_c:
                    off = (rho - rho_c)/c0 + 2*np.log(rho/rho_c)/c0**2 - 2*(1/rho - 1/rho_c)/c0**3
            else:
                off = rho/c0 if profile == "front" else 0.0
            want.append(off)
        got = ref.seed_times(centres, terminal, 4.0, M, c0, profile)
        np.testing.assert_allclose(got, 4.0 + np.array(want), rtol=1e-14, atol=1e-15)
    assert np.any(np.array(want) > 0)
    with pytest.raises(ValueError, match="profile"):
        ref.seed_offset(0.1, c0, "exact")


def test_parity_radial_slope_and_time():
    try:
        from parity.test_parity import compare, parity_available
    except ImportError:
        pytest.skip("parity harness absent (cardiacFoam-only)")
    ok, why = parity_available()
    if not ok:
        pytest.skip(why)
    ref = _reference()
    rng = np.random.default_rng(11)
    # Pe_L (a + 1) up to 1.1e3: the C++ quadrature, like the NumPy one, has no overflow limit
    for Pe, a in ((2.0, 0.05), (21.0, 0.25), (46.0, 0.25), (94.0, 0.5), (900.0, 0.2)):
        s = rng.uniform(a, a + 1, 40)
        pts = np.column_stack([s, np.zeros((40, 3))])

        def numpy_fn(p):
            return np.column_stack([ref.slope(p[:, 0], Pe, a), ref.scaled_time(p[:, 0], Pe, a)])

        compare("radial", 1, numpy_fn, pts, rtol=1e-12, Pe=Pe, a=a)
