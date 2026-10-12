# test_shell_bidomain_mms.py
"""Self-checks of the heart-only bidomain shell solution and its header."""
import os

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
