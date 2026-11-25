import os
import re
import tempfile

import pytest

import bath_shell_bidomain_mms as B


CASE = os.path.join(B.HERE, "..", "bathBidomain", "sphericalShell")


@pytest.mark.skipif(not os.path.isdir(CASE), reason="tutorial case absent (cardiacFoam-only)")
def test_tutorial_matches():
    """The tutorial's radii, ionic constants and defaults are the generator's."""
    p = B.default_params()
    geo = open(os.path.join(CASE, "shell.geo.template")).read()
    for name in ("r1", "r2", "r3"):
        assert float(re.search(r"^%s\s*=\s*([0-9.]+);" % name, geo, re.M).group(1)) == p[name]
    props = open(os.path.join(CASE, "constant", "electroProperties")).read()
    # Cm = 2 and Chi = 3 are the constants of the ionic closure (ionic_residual)
    assert re.search(r"\bchi\s*\[[^\]]*\]\s*(\S+);", props).group(1) == "3"
    assert re.search(r"\bcm\s*\[[^\]]*\]\s*(\S+);", props).group(1) == "2"
    assert re.search(r"\bouterClosure\s+(\S+);", props).group(1) == B.SOLUTIONS[0]
    assert re.search(r"\bparameterSet\s+(\S+);", props).group(1) == "default"
    allrun = open(os.path.join(CASE, "Allrun")).read()
    assert "oneOf -outer \"$outer\" " + " ".join(B.SOLUTIONS) in allrun
    assert "oneOf -parameterSet \"$parameterSet\" default alternate" in allrun


def test_shared_construction():
    B.check_shared_construction()


def test_interfaces():
    B.check_interfaces()


def test_outer_closure():
    B.check_outer_closure()


def test_sources_match_fd():
    B.check_sources_match_fd()


def test_bounded_and_finite():
    B.check_bounded_and_finite()


def test_denominators():
    B.check_denominators()


def test_richness():
    B.check_richness()


def test_compatibility():
    B.check_compatibility()


def test_common_mode_keeps_vm():
    B.check_common_mode_keeps_vm()


def test_ionic_closure():
    B.check_ionic_closure()


def test_reference_point():
    B.check_reference_point()


# Regenerating the header takes about half an hour (sympy cse of the
# sources of both outer closures): run it with pytest -m slow
@pytest.mark.slow
def test_header_check():
    assert B.header_matches(B.HEADER_PATH)


def test_emitted_files_up_to_date():
    assert B.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    B.write_outputs(str(tmp_path))
    for rel, text in B.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def test_write_refuses_after_failed_self_check():
    path = os.path.join(tempfile.mkdtemp(), "h.H")
    saved = B.SELF_CHECKS
    try:
        B.SELF_CHECKS = (lambda: B._require(False, "forced failure"),)
        assert B.main(["--write", path]) == 1
        assert B.main(["--emit"]) == 1
    finally:
        B.SELF_CHECKS = saved
    assert not os.path.exists(path)
