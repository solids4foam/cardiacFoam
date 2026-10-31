import importlib.util
import os
import re

import numpy as np
import pytest
import sympy as sp

import fda_bidomain_mms as fbm
import fda_monodomain_mms as fmm

DIMS = (1, 2, 3)
x, y, z, t = fbm.x, fbm.y, fbm.z, fbm.t
TUTORIAL = os.path.join(fbm.HERE, "..", "bidomain", "constant", "electroProperties")


def test_constants():
    assert (fbm.CM, fbm.CHI) == (2, 3)
    assert fbm.GI_DIAG == fmm.SIGMA_DIAG
    assert set(fbm.fields(3)) == {"Vm", "phiE", "phiI", "shift", "u1", "u2", "u3", "Iion",
                                  "F", "G"}


@pytest.mark.parametrize("d", DIMS)
def test_monodomain_field_and_states(d):
    """V_m, u_1, u_2, u_3 are those of the FDA monodomain problem."""
    f, m = fbm.fields(d), fmm.fields(d)
    for name in ("Vm", "u1", "u2", "u3", "F", "G"):
        assert sp.simplify(f[name] - m[name]) == 0


@pytest.mark.parametrize("d", DIMS)
def test_residuals_are_constant_multiples_of_Vm(d):
    """The residuals are rho V_m with rho of order 1e-9 k_d; with Ge exactly
    proportional to Gi both vanish."""
    for rho in fbm.residual_coefficients(d):
        assert rho.free_symbols == set()
    assert max(fbm.relative_residuals(d)) <= fbm.RESIDUAL_RTOL
    if d == 1:
        assert fbm.residual_coefficients(1) == (0, 0)


@pytest.mark.parametrize("d", DIMS)
def test_gauge_and_wall_fluxes(d):
    f = fbm.fields(d)
    at_ref = f["phiE"].subs(dict(zip(fbm.AXES, fbm.PHIE_REF_POINT)))
    assert sp.simplify(at_ref - fbm.PHIE_REF_VALUE) == 0
    for m in range(d):
        for face in (0, 1):
            for u in (f["Vm"], f["phiE"], f["phiI"]):
                assert sp.simplify(sp.diff(u, fbm.AXES[m]).subs(fbm.AXES[m], face)) == 0


def test_beta_is_minus_the_harmonic_eigenvalue():
    for d in DIMS:
        k_h = sum(fbm.harmonic(m)*((m + 1)*sp.pi)**2 for m in range(d))
        assert sp.simplify(fbm.beta(d) + k_h) == 0


def test_self_checks_and_cli():
    assert fbm.main(["--self-check"]) == 0


def test_no_header_options():
    with pytest.raises(SystemExit):
        fbm.main(["--check", "manufacturedFDABidomainReference.H"])


def test_emitted_files_up_to_date():
    assert fbm.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    fbm.write_outputs(str(tmp_path))
    for rel, text in fbm.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def _reference():
    path = os.path.join(fbm.HERE, "reference", "fda_bidomain.py")
    spec = importlib.util.spec_from_file_location("fda_bidomain_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("d", DIMS)
def test_emitted_numpy_matches_sympy(d):
    ref = _reference()
    pts = np.random.default_rng(5).uniform(0, 1, (6, 4))
    for name, expr in fbm.fields(d).items():
        fn = sp.lambdify((x, y, z, t), expr, "numpy")
        want = np.broadcast_to(fn(*pts.T), (6,))
        got = getattr(ref, "%s_d%d" % (name, d))(*pts.T)
        np.testing.assert_allclose(got, want, rtol=1e-13, atol=1e-15)


def test_numpy_constants():
    ref = _reference()
    assert abs(ref.PARAMS["KAPPA"] - float(fbm.kappa())) < 1e-16
    for d in DIMS:
        assert abs(ref.PARAMS["BETA"][d] - float(fbm.beta(d))) <= 1e-15*abs(float(fbm.beta(d)))


@pytest.mark.skipif(not os.path.exists(TUTORIAL), reason="tutorial case absent (cardiacFoam-only)")
def test_tutorial_dictionary_matches():
    text = open(TUTORIAL).read()
    gi = re.search(r"conductivityIntracellular\s*\[[^\]]*\]\s*\(([^)]*)\)", text).group(1)
    ge = re.search(r"conductivityExtracellular\s*\[[^\]]*\]\s*\(([^)]*)\)", text).group(1)
    ref = re.search(r"phiERefPoint\s*\(([^)]*)\)", text).group(1)
    as_diag = lambda v: [sp.Rational(s) for s in v.split()]  # noqa: E731
    assert as_diag(gi) == [fbm.GI_DIAG[0], 0, 0, fbm.GI_DIAG[1], 0, fbm.GI_DIAG[2]]
    assert as_diag(ge) == [fbm.GE_DIAG[0], 0, 0, fbm.GE_DIAG[1], 0, fbm.GE_DIAG[2]]
    assert tuple(sp.Rational(s) for s in ref.split()) == fbm.PHIE_REF_POINT
