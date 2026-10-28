import importlib.util
import os

import numpy as np
import pytest
import sympy as sp

import fda_monodomain_mms as fda

DIMS = (1, 2, 3)
x, y, z, t = fda.x, fda.y, fda.z, fda.t


def test_constants():
    assert fda.BETA == {1: -1.1, 2: -5.9, 3: -8.6}
    assert (fda.CM, fda.CHI) == (2, 3)
    assert set(fda.fields(3)) == {"Vm", "u1", "u2", "u3", "Iion", "F", "G"}


@pytest.mark.parametrize("d", DIMS)
def test_ionic_rates_hold(d):
    f = fda.fields(d)
    rates = fda.ionic_rates(f["Vm"], f["u1"], f["u2"], f["u3"])
    for name in ("u1", "u2", "u3"):
        assert sp.simplify(sp.diff(f[name], t) - rates[name]) == 0


@pytest.mark.parametrize("d", DIMS)
def test_pde_residual_is_k_plus_beta_times_V(d):
    """The residual is exactly (k_d + beta_d) V_m.  beta_d is printed to two
    digits but sigma (nine digits) is tuned to it, so the residual is nonzero
    yet below 1e-8 of |div(sigma grad V_m)| = k_d |V_m|."""
    f = fda.fields(d)
    k, beta = fda.eigenvalue(d), sp.nsimplify(fda.BETA[d], rational=True)
    assert sp.simplify(fda.pde_residual(d) - (k + beta)*f["Vm"]) == 0
    assert abs(float(k) + fda.BETA[d])/float(k) <= fda.RESIDUAL_RTOL
    # with the unrounded beta = -k the residual vanishes identically
    V, u1, u2, u3 = (f[n] for n in ("Vm", "u1", "u2", "u3"))
    exact = fda.ionic_current(V, u1, u2, u3, fda.CM, -k, fda.CHI)
    div = sum(fda.SIGMA_DIAG[m]*sp.diff(V, fda.AXES[m], 2) for m in range(d))
    assert sp.simplify(fda.CHI*fda.CM*sp.diff(V, t) - div + fda.CHI*exact) == 0


@pytest.mark.parametrize("d", DIMS)
def test_zero_flux_on_cube_faces(d):
    V = fda.fields(d)["Vm"]
    for m in range(d):
        for face in (0, 1):
            assert sp.simplify(sp.diff(V, fda.AXES[m]).subs(fda.AXES[m], face)) == 0


def test_self_checks_and_cli():
    assert fda.main(["--self-check"]) == 0


def test_emitted_files_up_to_date():
    assert fda.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    fda.write_outputs(str(tmp_path))
    for rel, text in fda.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def _reference():
    path = os.path.join(fda.HERE, "reference", "fda_monodomain.py")
    spec = importlib.util.spec_from_file_location("fda_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("d", DIMS)
def test_emitted_numpy_matches_sympy(d):
    ref = _reference()
    pts = np.random.default_rng(3).uniform(0, 1, (6, 4))
    for name, expr in fda.fields(d).items():
        fn = sp.lambdify((x, y, z, t), expr, "numpy")
        want = np.broadcast_to(fn(*pts.T), (6,))
        got = getattr(ref, "%s_d%d" % (name, d))(*pts.T)
        np.testing.assert_allclose(got, want, rtol=1e-13, atol=1e-15)


@pytest.mark.parametrize("d", DIMS)
def test_parity_fda(d):
    try:
        from parity.test_parity import compare, parity_available, sample_points
    except ImportError:
        pytest.skip("parity harness absent (cardiacFoam-only)")
    ok, why = parity_available()
    if not ok:
        pytest.skip(why)
    ref = _reference()

    def numpy_fn(p):
        a = [p[:, i] for i in range(4)]
        return np.column_stack([getattr(ref, "%s_d%d" % (n, d))(*a)
                                for n in ("Vm", "u1", "u2", "u3", "Iion")])

    compare("fda", d, numpy_fn, sample_points(), rtol=1e-12)


@pytest.mark.parametrize("name", sorted(fda.HEADER_PATHS))
def test_header_up_to_date(name):
    assert fda.header_matches(fda.HEADER_PATHS[name])


@pytest.mark.parametrize("name", sorted(fda.HEADER_PATHS))
def test_header_check_cli(name):
    assert fda.main(["--check", fda.HEADER_PATHS[name]]) == 0


def test_header_write_touches_only_regions(tmp_path):
    import mms_codegen as codegen
    for name, rel in fda.HEADER_PATHS.items():
        text = open(codegen.resolve(fda.__file__, rel)).read()
        stale = codegen.splice_regions(text, {n: ["// stale"] for n in fda.HEADERS[name](text)})
        assert stale != text
        dst = tmp_path/name
        dst.write_text(stale)
        assert fda.main(["--check", str(dst)]) == 1
        assert fda.main(["--write", str(dst)]) == 0
        assert dst.read_text() == text


def test_header_without_regions_is_refused(tmp_path):
    for name in fda.HEADER_PATHS:
        dst = tmp_path/name
        dst.write_text("namespace Foam {}\n")
        assert fda.main(["--check", str(dst)]) == 1
        assert fda.main(["--write", str(dst)]) == 1
        assert dst.read_text() == "namespace Foam {}\n"
