import importlib.util
import os

import numpy as np
import pytest
import sympy as sp

import rotated_monodomain_mms as rot

x, y, z, t = rot.x, rot.y, rot.z, rot.t


def test_K_spd_with_offdiagonals():
    K = rot.K_CASE
    assert K == K.T
    assert all(K[i, j] != 0 for i in range(3) for j in range(i))
    assert all(K[:k, :k].det() > 0 for k in (1, 2, 3))
    entries = [float(K[i, j]) for i, j in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))]
    assert entries == [0.106875, -0.0084931, -0.022561, 0.116682, -0.0130256, 0.0398782]


def test_gradient_zero_on_faces():
    F = rot.fields(rot.K_CASE)["F"]
    for a in rot.AXES:
        for face in (0, 1):
            assert sp.simplify(sp.diff(F, a).subs(a, face)) == 0


def test_source_matches_cxx_form():
    S = rot.fields(rot.K_CASE)["S"]
    assert sp.simplify(S - rot.cxx_source(rot.K_CASE)) == 0


def test_source_is_beta_V_minus_div_K_grad_V():
    f = rot.fields(rot.K_CASE)
    div = rot.div_K_grad(f["Vm"], rot.K_CASE)
    assert sp.simplify(f["S"] - (sp.Rational(str(rot.BETA))*f["Vm"] - div)) == 0


def test_rejects_nonsymmetric():
    with pytest.raises(ValueError):
        rot.fields(sp.Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]]))


def test_self_checks_and_cli():
    assert rot.main(["--check"]) == 0


def test_emitted_files_up_to_date():
    assert rot.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    rot.write_outputs(str(tmp_path))
    for rel, text in rot.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def _reference():
    path = os.path.join(rot.HERE, "reference", "rotated_monodomain.py")
    spec = importlib.util.spec_from_file_location("rot_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_emitted_numpy_matches_sympy():
    ref = _reference()
    pts = np.random.default_rng(4).uniform(0, 1, (6, 4))
    for name, expr in rot.fields(rot.K_CASE).items():
        fn = sp.lambdify((x, y, z, t), expr, "numpy")
        want = np.broadcast_to(fn(*pts.T), (6,))
        np.testing.assert_allclose(getattr(ref, name)(*pts.T), want, rtol=1e-13, atol=1e-15)


def test_parity_rotated():
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
        return np.column_stack([getattr(ref, n)(*a) for n in ("Vm", "u1", "u2", "u3", "S")])

    compare("rotated", 3, numpy_fn, sample_points(), rtol=1e-12)
