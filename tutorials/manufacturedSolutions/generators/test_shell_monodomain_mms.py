# test_shell_monodomain_mms.py
"""Self-checks of the heart-only monodomain shell solution and its header."""
import importlib.util
import os

import numpy as np
import pytest

import shell_monodomain_mms as smm

HERE = os.path.dirname(os.path.abspath(__file__))


def test_sealed_on_both_spheres():
    assert smm.max_sealing_residual(n=200) < 1e-25


def test_source_matches_finite_differences():
    assert smm.max_source_fd_error() < 1e-6


def test_positive_denominators_and_g():
    assert smm.min_denominator(band=0.05) > 0.0
    assert smm.min_g(band=0.05) > 0.5


def test_rich_on_the_walls():
    tangential_conormal, tangential_gradient = smm.wall_richness()
    assert tangential_conormal > 1e-2
    assert tangential_gradient > 1e-2


def test_header_up_to_date():
    assert smm.header_matches(smm.HEADER_PATH)


def _reference():
    path = os.path.join(HERE, "reference", "shell_monodomain.py")
    spec = importlib.util.spec_from_file_location("shell_monodomain_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _close(emitted, generator, tol=1e-13):
    """|emitted - generator| <= tol max(1, |generator|) componentwise."""
    emitted, generator = np.asarray(emitted, dtype=float), np.asarray(generator, dtype=float)
    assert emitted.shape == generator.shape
    assert np.all(np.abs(emitted - generator) <= tol*np.maximum(1.0, np.abs(generator)))


@pytest.mark.parametrize("which", ["default", "alternate"])
def test_numpy_reference_matches_lambdify(which):
    ref = _reference()
    params = smm.default_params() if which == "default" else smm.alternate_params()
    pv = smm.param_values(params)
    assert ref.PARAMS[which] == pv
    pts = smm._shell_points(50, 0.0, seed=21)
    X, Y, Z = pts.T
    args = [pv[n] for n in ref.PARAM_ORDER]
    _close(ref.F(X, Y, Z, *args), smm.evaluate("F", pts, params))
    _close(np.moveaxis(ref.gradF(X, Y, Z, *args), 0, -1), smm.evaluate("gradF", pts, params))
    _close(ref.divSigmaGradF(X, Y, Z, *args), smm.evaluate("divGGradF", pts, params))
    _close(ref.g(X, Y, Z, *args), smm.evaluate("g", pts, params))
    S = ref.sigma(X, Y, Z, *args)
    G6 = np.stack([S[0, 0], S[0, 1], S[0, 2], S[1, 1], S[1, 2], S[2, 2]], -1)
    _close(G6, smm.evaluate("G", pts, params))
    _close(ref.field("F", X, Y, Z), smm.evaluate("F", pts, smm.default_params()))


def test_emitted_outputs_up_to_date():
    assert smm.outputs_match()
