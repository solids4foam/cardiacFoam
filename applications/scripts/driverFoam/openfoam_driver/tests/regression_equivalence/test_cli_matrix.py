"""Matrix builder tests (solver-free: phase 2 not run)."""
from __future__ import annotations

from openfoam_driver.tests.regression_equivalence.__main__ import build_matrix


def test_matrix_has_row_per_case_driver():
    rows = build_matrix(run_phase2=False)
    # 6 mapped * 2 drivers + 3 unmapped * 1 driver = 15 rows.
    assert len(rows) == 15
    for row in rows:
        assert row["reproduces"] == "not-run"


def test_strict_rows_resolve_and_are_idempotent():
    rows = build_matrix(run_phase2=False)
    strict = [r for r in rows if r["driver"] == "strict"]
    assert len(strict) == 6
    assert all(r["resolves"] == "ok" for r in strict)
    assert all(r["idempotent"] == "ok" for r in strict)


def test_non_addressable_generic_row_is_unaddressable():
    rows = build_matrix(run_phase2=False)
    em = [r for r in rows
          if r["case"] == "NiedererEtAl2011/electroMechanicalNiedererEtAl2011"]
    assert len(em) == 1
    assert em[0]["driver"] == "generic"
    assert em[0]["resolves"] == "unaddressable"
