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
DOCS_REFERENCE = os.path.join(fbm.HERE, "docs", "bidomain", "bidomain_reference.py")


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


# ---------------------------------------------------- C++ header regions ----

F_S, FREF_S, K_S, C_S, REF_S = (sp.Symbol(n, real=True)
                                for n in ("F", "Fref", "k", "shift", "phiEReferenceValue"))


def _point_subs(d):
    F = fmm._F(d)
    Fref = F.subs(dict(zip(fbm.AXES, fbm.PHIE_REF_POINT)), simultaneous=True)
    return {F_S: F, FREF_S: Fref, K_S: fbm.kappa(), C_S: fbm.shift(d), REF_S: fbm.PHIE_REF_VALUE}


@pytest.mark.parametrize("d", DIMS)
def test_point_expressions_are_the_fields(d):
    """The expressions written into the header (F, k and shift as C++ symbols)
    equal the generator's fields once F = F_d, k = kappa and shift = c."""
    e, f, sub = fbm.point_expressions(), fbm.fields(d), _point_subs(d)
    for name, field in (("VmAt", "Vm"), ("PhiEAt", "phiE"), ("PhiIAt", "phiI"), ("Shift", "shift")):
        assert sp.simplify(e[name].subs(sub) - f[field]) == 0, name
    F = fmm._F(d)
    for name, field in (("GradVm", "Vm"), ("GradPhiE", "phiE"), ("GradPhiI", "phiI")):
        for a in fbm.AXES:
            assert sp.simplify(e[name].subs(sub)*sp.diff(F, a) - sp.diff(f[field], a)) == 0, name


@pytest.mark.parametrize("name", sorted(fbm.HEADER_PATHS))
def test_header_up_to_date(name):
    assert fbm.header_matches(fbm.HEADER_PATHS[name])


@pytest.mark.parametrize("name", sorted(fbm.HEADER_PATHS))
def test_header_check_cli(name):
    assert fbm.main(["--check", fbm.HEADER_PATHS[name]]) == 0


def test_header_write_touches_only_regions(tmp_path):
    import mms_codegen as codegen
    for name, rel in fbm.HEADER_PATHS.items():
        text = open(codegen.resolve(fbm.__file__, rel)).read()
        stale = codegen.splice_regions(text, {n: ["// stale"] for n in fbm.HEADERS[name](text)})
        assert stale != text
        dst = tmp_path/name
        dst.write_text(stale)
        assert fbm.main(["--check", str(dst)]) == 1
        assert fbm.main(["--write", str(dst)]) == 0
        assert dst.read_text() == text


def test_header_without_regions_is_refused(tmp_path):
    for name in fbm.HEADER_PATHS:
        dst = tmp_path/name
        dst.write_text("namespace Foam {}\n")
        assert fbm.main(["--check", str(dst)]) == 1
        assert fbm.main(["--write", str(dst)]) == 1
        assert dst.read_text() == "namespace Foam {}\n"


@pytest.mark.parametrize("d", DIMS)
def test_parity_bidomain(d):
    """The compiled header against reference/fda_bidomain.py and the sympy
    gradients, at k = kappa and the tutorial's gauge."""
    try:
        from parity.test_parity import compare, parity_available, sample_points
    except ImportError:
        pytest.skip("parity harness absent (cardiacFoam-only)")
    ok, why = parity_available()
    if not ok:
        pytest.skip(why)
    ref = _reference()
    f = fbm.fields(d)
    grads = [sp.lambdify((x, y, z, t), sp.diff(f[n], a), "numpy")
             for n in ("Vm", "phiE", "phiI") for a in fbm.AXES]

    def numpy_fn(p):
        a = [p[:, i] for i in range(4)]
        n = p.shape[0]
        cols = [getattr(ref, "%s_d%d" % (m, d))(*a) for m in ("phiE", "phiI", "shift")]
        cols += [getattr(ref, "%s_d%d" % (m, d))(*a) for m in ("phiE", "phiI", "Vm")]
        cols += [np.broadcast_to(g(*a), (n,)) for g in grads]
        return np.column_stack(cols)

    compare("bidomain", d, numpy_fn, sample_points(), rtol=1e-12, k="%.17g" % float(fbm.kappa()))


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


@pytest.mark.skipif(not os.path.exists(TUTORIAL), reason="tutorial case absent (cardiacFoam-only)")
def test_docs_reference_self_check():
    """docs/bidomain/bidomain_reference.py --self-check: the conductivities and
    phiERefPoint of bidomain and bidomain/sphericalShell equal the generators'."""
    spec = importlib.util.spec_from_file_location("bidomain_reference", DOCS_REFERENCE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.main(["--self-check"]) == 0
