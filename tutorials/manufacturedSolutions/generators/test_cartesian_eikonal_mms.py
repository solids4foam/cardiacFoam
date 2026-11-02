import importlib.util
import os

import numpy as np
import pytest
import sympy as sp

import cartesian_eikonal_mms as eik

DIMS = (1, 2, 3)
C0S = (1, 5, 66)
SIGMAS = {"rotated": eik.SIGMA_ROTATED, "axis": eik.SIGMA_AXIS}


def _reference():
    path = os.path.join(eik.HERE, "reference", "cartesian_eikonal.py")
    spec = importlib.util.spec_from_file_location("eik_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("c0", C0S)
@pytest.mark.parametrize("name", sorted(SIGMAS))
def test_eikonal_residual_zero(d, c0, name):
    r = eik.residual(d, eik.sigma_matrix(SIGMAS[name]), c0=c0)
    for p in ((0.1, 0.2, 0.3), (0.9, 0.4, 0.75)):
        val = r.subs(dict(zip(eik.AXES, p)))
        assert abs(float(val)) < 1e-12


@pytest.mark.parametrize("d", DIMS)
def test_fields_have_the_documented_form(d):
    f = eik.fields(d, eik.sigma_matrix(SIGMAS["rotated"]), c0=5)
    k = eik.k_vector(d)
    assert sp.simplify(f["psi"] - sp.exp((k.T*sp.Matrix(eik.AXES))[0])) == 0
    assert list(f["grad_psi"]) == list(k*f["psi"])
    assert sp.simplify(f["S_psi"] - (f["psi"]*(5*sp.sqrt(f["K_M"]) - f["K_M"]) - 1)) == 0


@pytest.mark.parametrize("bad", [sp.zeros(3, 3), -sp.eye(3), sp.diag(-1, 5, 5)])
def test_rejects_nonpositive_KM(bad):
    with pytest.raises(ValueError):
        eik.fields(1, bad)
    ref = _reference()
    with pytest.raises(ValueError):
        ref.source(1, (0.1, 0.2, 0.3), np.array(bad.tolist(), dtype=float), 1.0)


@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("name", sorted(SIGMAS))
def test_emitted_matches_generator(d, name):
    ref = _reference()
    sigma = eik.sigma_matrix(SIGMAS[name])
    rng = np.random.default_rng(3)
    p = rng.uniform(0, 1, (3, 6))
    for c0 in C0S:
        f = eik.fields(d, sigma, c0=c0)
        S = sp.lambdify(eik.AXES, f["S_psi"], "numpy")(*p)
        np.testing.assert_allclose(ref.source(d, p, np.array(sigma.tolist(), dtype=float), c0),
                                   S, rtol=1e-13, atol=1e-15)
    psi = sp.lambdify(eik.AXES, eik.psi_expr(d), "numpy")(*p)
    np.testing.assert_allclose(getattr(ref, "tau_d%d" % d)(*p), psi, rtol=1e-14)
    g = getattr(ref, "grad_tau_d%d" % d)(*p)
    np.testing.assert_allclose(g, np.array(eik.k_vector(d).tolist(), dtype=float)*psi, rtol=1e-14)


def test_self_checks_and_cli():
    assert eik.main(["--self-check"]) == 0


def test_emitted_files_up_to_date():
    assert eik.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    eik.write_outputs(str(tmp_path))
    for rel, text in eik.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def test_emitted_module_is_numpy_only():
    text = open(os.path.join(eik.HERE, "reference", "cartesian_eikonal.py")).read()
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


@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("c0", [1, 66])
@pytest.mark.parametrize("name", sorted(SIGMAS))
def test_parity_eikonal(d, c0, name):
    compare = _parity()
    ref = _reference()
    sigma = np.array(eik.sigma_matrix(SIGMAS[name]).tolist(), dtype=float)

    def numpy_fn(p):
        x, y, z = p[:, 0], p[:, 1], p[:, 2]
        tau = getattr(ref, "tau_d%d" % d)(x, y, z)
        g = getattr(ref, "grad_tau_d%d" % d)(x, y, z)
        S = ref.source(d, (x, y, z), sigma, c0)
        return np.column_stack([tau, g[0], g[1], g[2], S])

    rng = np.random.default_rng(11)
    pts = np.column_stack([rng.uniform(0, 1, (20, 3)), rng.uniform(0, 0.1, 20)])
    compare("eikonal", d, numpy_fn, pts, rtol=1e-12, c0=c0, conductivity=name)


@pytest.mark.parametrize("name", sorted(eik.HEADER_PATHS))
def test_header_up_to_date(name):
    assert eik.header_matches(eik.HEADER_PATHS[name])


@pytest.mark.parametrize("name", sorted(eik.HEADER_PATHS))
def test_header_check_cli(name):
    assert eik.main(["--check", eik.HEADER_PATHS[name]]) == 0


def test_header_write_touches_only_regions(tmp_path):
    import mms_codegen as codegen
    for name, rel in eik.HEADER_PATHS.items():
        text = open(codegen.resolve(eik.__file__, rel)).read()
        stale = codegen.splice_regions(text, {n: ["// stale"] for n in eik.HEADERS[name](text)})
        assert stale != text
        dst = tmp_path/name
        dst.write_text(stale)
        assert eik.main(["--check", str(dst)]) == 1
        assert eik.main(["--write", str(dst)]) == 0
        assert dst.read_text() == text


def test_header_without_regions_is_refused(tmp_path):
    for name in eik.HEADER_PATHS:
        dst = tmp_path/name
        dst.write_text("namespace Foam {}\n")
        assert eik.main(["--check", str(dst)]) == 1
        assert eik.main(["--write", str(dst)]) == 1
        assert dst.read_text() == "namespace Foam {}\n"
