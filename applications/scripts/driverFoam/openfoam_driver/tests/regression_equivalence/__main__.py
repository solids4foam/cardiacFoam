"""CLI: print the regression-equivalence matrix.

    python -m openfoam_driver.tests.regression_equivalence [--run-phase2]

Phase 1 (solver-free) always runs: agent dict-layer idempotence + agent
addressability of each case/driver. Phase 2 (--run-phase2) additionally runs the
dual-run numeric comparison (canonical vs agent-driven), which self-skips
without a built cardiacFoam.
"""
from __future__ import annotations

import argparse
from typing import Any

from openfoam_driver.tests.regression_equivalence.dual_run import verify_reproduction
from openfoam_driver.tests.regression_equivalence.registry import REGRESSION_CASES
from openfoam_driver.tests.regression_equivalence.round_trip import (
    electro_build_parse_fixpoint,
)
from openfoam_driver.tests.regression_equivalence.staging import (
    resolve_generic,
    resolve_strict,
)


def _resolves(case, driver: str) -> str:
    try:
        res = resolve_strict(case) if driver == "strict" else resolve_generic(case)
    except KeyError:
        return "unaddressable"
    return "ok" if res.get("is_runnable") else "not-runnable"


def _idempotent(case) -> str:
    if not case.mapped:
        return "n/a"
    once, twice = electro_build_parse_fixpoint(case)
    return "ok" if once == twice else "DRIFT"


def build_matrix(*, run_phase2: bool) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for case in REGRESSION_CASES:
        for driver in case.drivers:
            repro = "not-run"
            if run_phase2:
                repro = verify_reproduction(case, driver=driver).status
            rows.append({
                "case": case.case_dir,
                "driver": driver,
                "resolves": _resolves(case, driver),
                "idempotent": _idempotent(case) if driver == "strict" else "n/a",
                "reproduces": repro,
            })
    return rows


def _print_matrix(rows: list[dict[str, Any]]) -> None:
    width = max(len(r["case"]) for r in rows)
    header = f"{'case':<{width}}  driver   resolves       idempotent  reproduces"
    print(header)
    print("-" * len(header))
    for r in rows:
        print(f"{r['case']:<{width}}  {r['driver']:<7}  "
              f"{r['resolves']:<13}  {r['idempotent']:<10}  {r['reproduces']}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-phase2", action="store_true")
    args = parser.parse_args()
    rows = build_matrix(run_phase2=args.run_phase2)
    _print_matrix(rows)
    failed = [
        r for r in rows
        if r["resolves"] == "not-runnable"
        or r["idempotent"] == "DRIFT"
        or r["reproduces"] in {"mismatch", "run_failed"}
    ]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
