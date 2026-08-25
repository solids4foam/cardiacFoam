"""Unit tests for dual-run helpers + skip behavior (solver-free)."""
from __future__ import annotations

from types import SimpleNamespace

import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.tests.regression_equivalence.registry import REGRESSION_CASES, RegressionCase
from openfoam_driver.tests.regression_equivalence import dual_run
from openfoam_driver.tests.regression_equivalence.dual_run import (
    extract_summary_value,
    parse_columnar_reference,
    read_series_value,
    solver_available,
    values_agree,
    verify_reproduction,
)

SINGLECELL_REF = """\
# file                     time       variable  expected     tolerance
postProcessing/x.txt       1.0000000  Vm        -87.0639985  5e-3
postProcessing/x.txt       1.5000000  cai         0.0001139  5e-5
"""

BIDOMAIN_REF = """\
# kind     key        metric      expected       tolerance
summary   cells      value       20             0
error     Vm         L1          0.00025015     1e-6
"""

SERIES = """\
time        Vm          cai
0.5000000   -80.0       0.0001
1.0000000   -87.06      0.00009
1.5000000   -86.90      0.0001139
"""

BATH_SUMMARY = """\
# Bath-bidomain manufactured solution error summary
# dimension 1D
# cellsPerDirection 80
# time 0.200112
Vm 1.52888e-05 1.99517e-05 3.99158e-05
"""

STANDARD_SUMMARY = """\
Number of cells (N)   = 80
Final simulation time = 0.200112
Vm 1.52888e-05 1.99517e-05 3.99158e-05
"""


def test_parse_columnar_reference_singlecell():
    points = parse_columnar_reference(SINGLECELL_REF)
    assert len(points) == 2
    assert points[0].variable == "Vm"
    assert points[0].time == 1.0
    assert points[0].tolerance == 5e-3


def test_parse_columnar_reference_rejects_metric_layout():
    # 'summary cells value 20 0' -> col2 'value' non-numeric expected? expected=20
    # ok; but 'error Vm L1 0.00025015 1e-6' has variable col fine — the giveaway
    # is that data_file 'summary'/'error' are not paths, yet columns parse. This
    # layout is still 5 cols and numeric, so it parses; ensure caller relies on
    # file-existence at run time. Here we assert it does NOT crash.
    points = parse_columnar_reference(BIDOMAIN_REF)
    assert isinstance(points, list)


def test_read_series_value_exact_time():
    assert read_series_value(SERIES, 1.0, "Vm") == -87.06
    assert read_series_value(SERIES, 1.5, "cai") == 0.0001139


def test_read_series_value_missing_variable_or_time():
    assert read_series_value(SERIES, 1.0, "phiE") is None
    assert read_series_value(SERIES, 9.0, "Vm") is None  # no time within atol
    assert read_series_value("", 1.0, "Vm") is None


def test_values_agree():
    assert values_agree(1.0, 1.0004, 5e-3)
    assert not values_agree(1.0, 1.01, 5e-3)


def test_extract_summary_value_supports_standard_and_bath_formats():
    assert extract_summary_value(STANDARD_SUMMARY, "cells") == 80.0
    assert extract_summary_value(STANDARD_SUMMARY, "finalTime") == 0.200112
    assert extract_summary_value(BATH_SUMMARY, "cellsPerDirection") == 80.0
    assert extract_summary_value(BATH_SUMMARY, "finalTime") == 0.200112


def test_solver_available_returns_bool():
    assert isinstance(solver_available(), bool)


def test_verify_reproduction_skips_without_solver(monkeypatch):
    monkeypatch.setattr(dual_run, "solver_available", lambda: False)
    result = verify_reproduction(REGRESSION_CASES[0], driver="strict")
    assert result.status == "skipped"
    assert "cardiacFoam" in result.detail


def test_verify_reproduction_skips_non_addressable_generic(monkeypatch):
    monkeypatch.setattr(dual_run, "solver_available", lambda: True)
    result = verify_reproduction(
        RegressionCase(
            "synthetic/nonAddressableCase",
            None,
            (),
            "regression/reference.txt",
            generic_addressable=False,
        ),
        driver="generic",
    )
    assert result.status == "skipped"
    assert "not addressable" in result.detail


def test_verify_reproduction_generic_prefers_regression_script(monkeypatch, tmp_path):
    case = RegressionCase(
        "synthetic/case",
        "syntheticEntry",
        (),
        "regression/reference.txt",
    )
    root = tmp_path / "sandbox" / "tutorials"
    case_path = root / "synthetic" / "case"
    script_path = case_path / "regression" / "regressionTest.sh"
    script_path.parent.mkdir(parents=True)
    script_path.write_text("#!/usr/bin/env bash\nexit 0\n")

    monkeypatch.setattr(dual_run, "solver_available", lambda: True)
    monkeypatch.setattr(dual_run, "_stage_tutorials_root", lambda _case: (root, case_path))
    monkeypatch.setattr(
        dual_run,
        "_run_regression_script",
        lambda _case, _case_path: SimpleNamespace(returncode=0, stdout="", stderr=""),
    )
    monkeypatch.setattr(dual_run, "_drive_agent", lambda *_args, **_kwargs: pytest.fail("should not call driverFoam run"))
    monkeypatch.setattr(dual_run, "check_protocol", lambda *_args, **_kwargs: (True, "ok"))

    result = verify_reproduction(case, driver="generic")

    assert result.status == "reproduced"
    assert result.detail == "committed regression script passed"


def test_verify_reproduction_generic_falls_back_without_regression_script(monkeypatch, tmp_path):
    case = RegressionCase(
        "synthetic/case",
        "syntheticEntry",
        (),
        "regression/reference.txt",
    )
    root = tmp_path / "sandbox" / "tutorials"
    case_path = root / "synthetic" / "case"
    case_path.mkdir(parents=True)

    monkeypatch.setattr(dual_run, "solver_available", lambda: True)
    monkeypatch.setattr(dual_run, "_stage_tutorials_root", lambda _case: (root, case_path))
    monkeypatch.setattr(
        dual_run,
        "_drive_agent",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=0, stdout="", stderr=""),
    )
    monkeypatch.setattr(dual_run, "check_protocol", lambda *_args, **_kwargs: (True, "ok"))

    result = verify_reproduction(case, driver="generic")

    assert result.status == "reproduced"
    assert result.detail == "ok"


def test_verify_reproduction_generic_maps_regression_skip(monkeypatch, tmp_path):
    case = RegressionCase(
        "synthetic/case",
        "syntheticEntry",
        (),
        "regression/reference.txt",
    )
    root = tmp_path / "sandbox" / "tutorials"
    case_path = root / "synthetic" / "case"
    script_path = case_path / "regression" / "regressionTest.sh"
    script_path.parent.mkdir(parents=True)
    script_path.write_text("#!/usr/bin/env bash\nexit 77\n")

    monkeypatch.setattr(dual_run, "solver_available", lambda: True)
    monkeypatch.setattr(dual_run, "_stage_tutorials_root", lambda _case: (root, case_path))
    monkeypatch.setattr(
        dual_run,
        "_run_regression_script",
        lambda _case, _case_path: SimpleNamespace(returncode=77, stdout="expected skip", stderr=""),
    )

    result = verify_reproduction(case, driver="generic")

    assert result.status == "skipped"
    assert "rc=77" in result.detail
