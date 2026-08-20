"""Matrix builder tests (solver-free: phase 2 not run)."""
from __future__ import annotations

from types import SimpleNamespace

import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.tests.regression_equivalence import __main__ as cli_matrix
from openfoam_driver.tests.regression_equivalence.__main__ import build_matrix_iter
from openfoam_driver.tests.regression_equivalence.registry import REGRESSION_CASES


def test_matrix_has_row_per_case_driver():
    rows = list(build_matrix_iter(run_phase2=False))
    expected_rows = sum(len(case.drivers) for case in REGRESSION_CASES)
    assert len(rows) == expected_rows
    for row in rows:
        assert row["reproduces"] == "not-run"
        assert row["reproduces_detail"] == ""


def test_strict_rows_resolve_and_are_idempotent():
    rows = list(build_matrix_iter(run_phase2=False))
    strict = [r for r in rows if r["driver"] == "strict"]
    assert len(strict) == sum(1 for case in REGRESSION_CASES if case.mapped)
    assert all(r["resolves"] == "ok" for r in strict)
    assert all(r["idempotent"] == "ok" for r in strict)


def test_electromechanical_generic_row_is_addressable():
    rows = list(build_matrix_iter(run_phase2=False))
    em = [r for r in rows
          if r["case"] == "NiedererEtAl2011/electroMechanicalNiedererEtAl2011"]
    assert len(em) == 1
    assert em[0]["driver"] == "generic"
    assert em[0]["resolves"] == "ok"


def test_phase2_rows_keep_reproduction_detail(monkeypatch):
    def fake_verify_reproduction(case, *, driver):
        return SimpleNamespace(status="run_failed", detail=f"{case.case_dir} via {driver}")

    monkeypatch.setattr(cli_matrix, "verify_reproduction", fake_verify_reproduction)
    rows = list(build_matrix_iter(run_phase2=True))
    strict = next(row for row in rows if row["driver"] == "strict")
    generic = next(row for row in rows if row["driver"] == "generic")
    assert strict["reproduces"] == "skipped"
    assert "committed-case regression bypasses registered-entry defaults" in strict["reproduces_detail"]
    assert generic["reproduces"] == "run_failed"
    assert generic["reproduces_detail"] == f"{generic['case']} via {generic['driver']}"
