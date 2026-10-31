import os
import tempfile

import shell_mms as S


def test_shared_construction():
    S.check_shared_construction()


def test_interfaces():
    S.check_all_interfaces()


def test_sources_match_fd():
    S.check_sources_match_fd()


def test_bounded_and_finite():
    S.check_bounded_and_finite()


def test_denominators():
    S.check_denominators()


def test_richness():
    S.check_richness()


def test_common_mode_keeps_vm():
    S.check_common_mode_keeps_vm()


def test_ionic_closure():
    S.check_ionic_closure()


def test_header_check():
    assert S.header_matches("../../../src/verificationModels/bathBidomainVerification/manufacturedShellTransmissionBathBidomainReference.H")


def test_emitted_files_up_to_date():
    assert S.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    S.write_outputs(str(tmp_path))
    for rel, text in S.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def test_write_refuses_after_failed_self_check():
    path = os.path.join(tempfile.mkdtemp(), "h.H")
    saved = S.SELF_CHECKS
    try:
        S.SELF_CHECKS = (lambda: S._require(False, "forced failure"),)
        assert S.main(["--write", path]) == 1
        assert S.main(["--emit"]) == 1
    finally:
        S.SELF_CHECKS = saved
    assert not os.path.exists(path)


if __name__ == "__main__":
    for f in (test_shared_construction, test_interfaces, test_sources_match_fd,
              test_bounded_and_finite, test_denominators, test_richness,
              test_common_mode_keeps_vm, test_ionic_closure, test_header_check,
              test_emitted_files_up_to_date, test_write_refuses_after_failed_self_check):
        f()
    print("ALL SHELL GENERATOR TESTS PASS")
