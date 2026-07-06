from openfoam_driver.dashboard.run_status import resolve_status


def test_ran_when_numeric_time_dir_present(tutorial_tree):
    case = tutorial_tree / "manufacturedSolutions" / "monodomainPseudoECG"
    st = resolve_status(case)
    assert st.ran is True
    assert st.last_run is not None


def test_not_run_when_no_output(tutorial_tree):
    case = tutorial_tree / "PATHOS" / "RBBB"
    st = resolve_status(case)
    assert st.ran is False


def test_regression_available_detected(tutorial_tree):
    case = tutorial_tree / "PATHOS" / "RBBB"
    (case / "regressionTest.sh").write_text("#!/bin/sh\n")
    (case / "RBBB.reference").write_text("ref\n")
    st = resolve_status(case)
    assert st.regression == "available"


def test_regression_none_without_reference(tutorial_tree):
    case = tutorial_tree / "PATHOS" / "RBBB"
    (case / "regressionTest.sh").write_text("#!/bin/sh\n")
    st = resolve_status(case)
    assert st.regression == "none"
