# test_shell_monodomain_mms.py
"""Self-checks of the heart-only monodomain shell solution and its header."""
import os

import numpy as np

import shell_mms
import shell_monodomain_mms as smm
import transmission_mms

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


def test_equals_transmission_phiI_at_t0():
    assert smm.max_difference_to_transmission_phiI() < 1e-13


def test_header_up_to_date():
    assert smm.header_matches(smm.HEADER_PATH)


def test_transmission_headers_unchanged():
    assert shell_mms.header_matches("../" + shell_mms.HEADER_NAME)
    assert transmission_mms.header_matches("../" + transmission_mms.HEADER_NAME)
