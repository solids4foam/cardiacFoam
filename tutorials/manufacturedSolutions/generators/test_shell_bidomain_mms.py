# test_shell_bidomain_mms.py
"""Self-checks of the heart-only bidomain shell solution and its header."""
import os

import numpy as np
import pytest

import shell_bidomain_mms as sbm
import shell_monodomain_mms as smm

HERE = os.path.dirname(os.path.abspath(__file__))


def test_sealed_on_both_spheres():
    assert sbm.max_sealing_residual(n=200) < 1e-25


def test_source_matches_finite_differences():
    assert sbm.max_source_fd_error() < 1e-6


def test_positive_denominators_and_g():
    assert sbm.min_denominator(band=0.05) > 0.0
    assert sbm.min_g(band=0.05) > 0.5


def test_rich_on_the_walls():
    for value in sbm.wall_richness().values():
        assert value > 1e-2


def test_phiE_source_has_zero_mean():
    mean, absolute = sbm.phiE_source_mean()
    assert abs(mean) < sbm.SOURCE_MEAN_TOLERANCE*absolute


def test_Fi_equals_monodomain_F():
    assert sbm.max_difference_to_monodomain_F() < 1e-13


def test_header_up_to_date():
    assert sbm.header_matches(sbm.HEADER_PATH)


def test_other_headers_unchanged():
    assert smm.header_matches(smm.HEADER_PATH)


def test_one_frame_unequal_ratios():
    comm, off, mono = sbm.frame_check()
    assert comm < 1e-14
    assert off > 0.1
    assert mono < 1e-14


def test_vm_equation_and_ionic_rates():
    assert all(v == 0 for v in sbm.vm_equation_residual())


def test_emitted_outputs_up_to_date():
    assert sbm.outputs_match()


def _reference():
    import importlib.util
    path = os.path.join(HERE, "reference", "shell_bidomain.py")
    spec = importlib.util.spec_from_file_location("shell_bidomain_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("which", ["default", "alternate"])
def test_numpy_reference_matches_lambdify(which):
    ref = _reference()
    params = sbm.default_params() if which == "default" else sbm.alternate_params()
    pv = sbm.param_values(params)
    assert ref.PARAMS[which] == pv
    pts = sbm._shell_points(50, 0.0, seed=21)
    X, Y, Z = pts.T
    args = [pv[n] for n in ref.PARAM_ORDER]

    def close(a, b, tol=1e-12):
        a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
        assert a.shape == b.shape
        assert np.all(np.abs(a - b) <= tol*np.maximum(1.0, np.abs(b)))

    for name in ("Fi", "Fe", "divGiGradFi", "divGeGradFe", "g"):
        close(getattr(ref, name)(X, Y, Z, *args), sbm.evaluate(name, pts, params))
    for name in ("gradFi", "gradFe"):
        close(np.moveaxis(getattr(ref, name)(X, Y, Z, *args), 0, -1),
              sbm.evaluate(name, pts, params))
    for name in ("Gi", "Ge"):
        T = getattr(ref, name)(X, Y, Z, *args)
        G6 = np.stack([T[0, 0], T[0, 1], T[0, 2], T[1, 1], T[1, 2], T[2, 2]], -1)
        close(G6, sbm.evaluate(name, pts, params))


SHELL_CASE = os.path.join(HERE, "..", "bidomain", "sphericalShell", "constant",
                          "electroProperties")


@pytest.mark.skipif(not os.path.exists(SHELL_CASE),
                    reason="tutorial case absent (cardiacFoam-only)")
def test_tutorial_reference_point():
    import re
    text = open(SHELL_CASE).read()
    ref = re.search(r"phiERefPoint\s*\(([^)]*)\)", text).group(1)
    assert tuple(float(v) for v in ref.split()) == sbm.PHIE_REF_POINT
