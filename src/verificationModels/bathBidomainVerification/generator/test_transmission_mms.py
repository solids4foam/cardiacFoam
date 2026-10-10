import numpy as np
import transmission_mms as T


def test_interfaces():
    T.check_all_interfaces()

def test_sources_match_fd():
    T.check_sources_match_fd()

def test_bounded_and_finite():
    T.check_bounded_and_finite()

def test_periodic():
    T.check_periodic()

def test_richness():
    T.check_richness()

def test_common_mode_keeps_vm():
    T.check_common_mode_keeps_vm()

def test_affine():
    T.check_affine()

def test_ionic_closure():
    T.check_ionic_closure()

def test_header_check():
    assert T.header_matches("../manufacturedTransmissionBathBidomainReference.H")

def test_write_refuses_after_failed_self_check():
    import os, tempfile
    path = os.path.join(tempfile.mkdtemp(), "h.H")
    saved = T.SELF_CHECKS
    try:
        T.SELF_CHECKS = (lambda: T._require(False, "forced failure"),)
        assert T.main(["--write", path]) == 1
    finally:
        T.SELF_CHECKS = saved
    assert not os.path.exists(path)

if __name__ == "__main__":
    for f in (test_interfaces, test_sources_match_fd, test_bounded_and_finite, test_periodic,
              test_richness, test_common_mode_keeps_vm, test_affine, test_ionic_closure,
              test_header_check, test_write_refuses_after_failed_self_check):
        f()
    print("ALL GENERATOR TESTS PASS")
