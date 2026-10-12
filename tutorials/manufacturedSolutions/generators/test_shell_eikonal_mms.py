# test_shell_eikonal_mms.py
"""Self-checks of the heart-only eikonal shell solution and its header."""
import functools
import os

import numpy as np

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


if __name__ == "__main__":
    for f in (test_sealed_on_both_spheres, test_stimulus_metric_is_the_inverse_of_M_at_xs,
              test_gradient_matches_finite_differences_of_psi, test_source_matches_finite_differences,
              test_positive_denominators, test_source_is_finite,
              test_critical_points_are_complete_and_nondegenerate,
              test_psi_is_positive_outside_the_pinned_ball, test_header_up_to_date,
              test_other_headers_unchanged):
        f()
    print("ALL GENERATOR TESTS PASS")
