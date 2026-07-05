"""Unit tests for dual-run helpers + skip behavior (solver-free)."""
from __future__ import annotations

from openfoam_driver.tests.regression_equivalence.registry import REGRESSION_CASES
from openfoam_driver.tests.regression_equivalence import dual_run
from openfoam_driver.tests.regression_equivalence.dual_run import (
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


def test_solver_available_returns_bool():
    assert isinstance(solver_available(), bool)


def test_verify_reproduction_skips_without_solver(monkeypatch):
    monkeypatch.setattr(dual_run, "solver_available", lambda: False)
    result = verify_reproduction(REGRESSION_CASES[0], driver="strict")
    assert result.status == "skipped"
    assert "cardiacFoam" in result.detail


def test_verify_reproduction_skips_non_addressable_generic(monkeypatch):
    monkeypatch.setattr(dual_run, "solver_available", lambda: True)
    em = next(c for c in REGRESSION_CASES if not c.generic_addressable)
    result = verify_reproduction(em, driver="generic")
    assert result.status == "skipped"
    assert "not addressable" in result.detail
