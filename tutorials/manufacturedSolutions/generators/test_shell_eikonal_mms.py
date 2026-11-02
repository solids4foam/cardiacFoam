# test_shell_eikonal_mms.py
"""Self-checks of the heart-only eikonal shell solution and its header."""
import functools
import importlib.util
import os

import numpy as np
import pytest

import shell_bidomain_mms
import shell_eikonal_mms as sem
import shell_monodomain_mms as smm

HERE = os.path.dirname(os.path.abspath(__file__))


@functools.lru_cache(maxsize=None)
def critical_points():
    return sem.critical_points()


def test_sealed_on_both_spheres():
    assert sem.max_sealing_residual(n=200) < 1e-25


def test_stimulus_metric_is_the_inverse_of_M_at_xs():
    for params in (sem.default_params(), sem.alternate_params()):
        G6 = sem.evaluate("G", params["xs"], params, "mpmath")
        M = np.array([[G6[0], G6[1], G6[2]], [G6[1], G6[3], G6[4]], [G6[2], G6[4], G6[5]]],
                     dtype=float)/(params["chi"]*params["cm"])
        m = params["Msi"]
        Msi = np.array([[m[0], m[1], m[2]], [m[1], m[3], m[4]], [m[2], m[4], m[5]]])
        assert np.abs(Msi @ M - np.eye(3)).max() < 1e-12


def test_gradient_matches_finite_differences_of_psi():
    assert sem.max_gradient_fd_error() < 1e-5


def test_source_matches_finite_differences():
    assert sem.max_source_fd_error() < 1e-6


def test_positive_denominators():
    assert sem.min_denominator(band=0.05) > 0.0


def test_source_is_finite():
    stats = sem.source_statistics()
    assert stats["finite"]
    near = [c["x"] + 1e-6*np.array([0.6, -0.5, 0.4]) for c in critical_points()]
    assert np.isfinite(sem.evaluate("source", np.array(near))).all()


def test_critical_points_are_complete_and_nondegenerate():
    cps = critical_points()
    assert sem.morse_count(cps, 1) == 2 and sem.morse_count(cps, 2) == 2
    assert all(c["wall"] or c["inPin"] for c in cps)
    assert sem.min_nondegeneracy(cps) > 1e-3


def test_psi_is_positive_outside_the_pinned_ball():
    cps = critical_points()
    assert sem.smallest_value_outside_pin(cps) > 0
    assert sem.pinned_minimum(cps) > -1e-2


def test_header_up_to_date():
    assert sem.header_matches(sem.HEADER_PATH)


def test_other_headers_unchanged():
    assert smm.header_matches(smm.HEADER_PATH)
    assert shell_bidomain_mms.header_matches(shell_bidomain_mms.HEADER_PATH)


def _with_c0(c0):
    params = sem.default_params()
    params["c0"] = float(c0)
    return params


@pytest.mark.parametrize("c0", [1, 2, 5])
def test_self_checks_pass_at_c0(c0):
    cps = sem.run_self_checks(_with_c0(c0))
    assert sem.morse_count(cps, 1) == 2 and sem.morse_count(cps, 2) == 2
    assert sem.min_nondegeneracy(cps) > 1e-3


def test_c0_scales_psi_and_keeps_the_critical_points():
    pts = sem._shell_points(20, 0.0, seed=5)
    for c0 in (2.0, 5.0):
        assert np.allclose(sem.evaluate("psi", pts, _with_c0(c0)),
                           sem.evaluate("psi", pts)/c0, rtol=1e-13, atol=0)


def test_stimulus_does_not_depend_on_the_summation_order(monkeypatch):
    def left_to_right(values, start=0.0):
        total = start
        for v in values:
            total = total + v
        return total
    before = (sem.default_params(), sem.alternate_params())
    monkeypatch.setattr(sem, "sum", left_to_right, raising=False)
    assert (sem.default_params(), sem.alternate_params()) == before


@pytest.mark.parametrize("argv", [["--report", "--c0", "0"], ["--report", "--c0", "-1"],
                                  ["--report", "--c0", "abc"], ["--report", "--c0", "nan"],
                                  ["--check", sem.HEADER_PATH, "--c0", "2"],
                                  ["--write", "unused.H", "--c0", "2"], ["--emit", "--c0", "2"]])
def test_cli_refuses_bad_or_misplaced_c0(argv):
    with pytest.raises(SystemExit) as err:
        sem.main(argv)
    assert err.value.code == 2


def _reference():
    path = os.path.join(HERE, "reference", "shell_eikonal.py")
    spec = importlib.util.spec_from_file_location("shell_eikonal_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _close(emitted, generator, tol=1e-13):
    """|emitted - generator| <= tol max(1, |generator|) componentwise."""
    emitted, generator = np.asarray(emitted, dtype=float), np.asarray(generator, dtype=float)
    assert emitted.shape == generator.shape
    assert np.all(np.abs(emitted - generator) <= tol*np.maximum(1.0, np.abs(generator)))


@pytest.mark.parametrize("which", ["default", "alternate", "c0=5"])
def test_numpy_reference_matches_lambdify(which):
    ref = _reference()
    params = {"default": sem.default_params(), "alternate": sem.alternate_params(),
              "c0=5": _with_c0(5)}[which]
    pv = sem.param_values(params)
    if which in ref.PARAMS:
        assert ref.PARAMS[which] == pv
    pts = sem._shell_points(50, 0.0, seed=22)
    X, Y, Z = pts.T
    args = [pv[n] for n in ref.PARAM_ORDER]
    _close(ref.psi(X, Y, Z, *args), sem.evaluate("psi", pts, params))
    _close(np.moveaxis(ref.gradPsi(X, Y, Z, *args), 0, -1), sem.evaluate("gradPsi", pts, params))
    _close(ref.source(X, Y, Z, *args), sem.evaluate("source", pts, params))
    S = ref.sigma(X, Y, Z, *args)
    G6 = np.stack([S[0, 0], S[0, 1], S[0, 2], S[1, 1], S[1, 2], S[2, 2]], -1)
    _close(G6, sem.evaluate("G", pts, params))
    _close(ref.field("psi", X, Y, Z, c0=5.0), sem.evaluate("psi", pts, _with_c0(5)))


def test_emitted_outputs_up_to_date():
    assert sem.outputs_match()


def test_latex_names_distinct_without_case():
    names = list(sem.latex_fragments()) + list(smm.latex_fragments())
    assert len({n.lower() for n in names}) == len(names)


if __name__ == "__main__":
    for f in (test_sealed_on_both_spheres, test_stimulus_metric_is_the_inverse_of_M_at_xs,
              test_gradient_matches_finite_differences_of_psi, test_source_matches_finite_differences,
              test_positive_denominators, test_source_is_finite,
              test_critical_points_are_complete_and_nondegenerate,
              test_psi_is_positive_outside_the_pinned_ball, test_header_up_to_date,
              test_other_headers_unchanged):
        f()
    print("ALL GENERATOR TESTS PASS")
