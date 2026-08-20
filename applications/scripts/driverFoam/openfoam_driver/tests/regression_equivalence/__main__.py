"""CLI: print the regression-equivalence matrix.

    python -m openfoam_driver.tests.regression_equivalence [--run-phase2]

Phase 1 (solver-free) always runs: agent dict-layer idempotence + agent
addressability of each case/driver.

Phase 2 (--run-phase2) runs the numeric comparison only through the committed
case-folder path. This is deliberate: numerical regression must execute the
checked-in case state (its committed dictionaries and Allrun), not a registered
tutorial family's default case selection or workflow interpretation. Registered
entry checks therefore remain in Phase 1 (resolution/idempotence), while
Phase 2 reproduces the committed case exactly as authored on disk.
"""
from __future__ import annotations

import argparse
from collections.abc import Iterator
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


def build_matrix_iter(*, run_phase2: bool) -> Iterator[dict[str, Any]]:
    for case in REGRESSION_CASES:
        for driver in case.drivers:
            repro = "not-run"
            repro_detail = ""
            if run_phase2:
                if driver == "generic":
                    result = verify_reproduction(case, driver=driver)
                    repro = result.status
                    repro_detail = result.detail
                else:
                    repro = "skipped"
                    repro_detail = (
                        "Phase 2 committed-case regression bypasses registered-entry "
                        "defaults/workflows and runs only the case_folder path."
                    )
            yield {
                "case": case.case_dir,
                "driver": driver,
                "resolves": _resolves(case, driver),
                "idempotent": _idempotent(case) if driver == "strict" else "n/a",
                "reproduces": repro,
                "reproduces_detail": repro_detail,
            }


def _print_header(width: int) -> None:
    header = f"{'case':<{width}}  driver   resolves       idempotent  reproduces"
    print(header)
    print("-" * len(header))


def _print_row(r: dict[str, Any], width: int) -> None:
    print(f"{r['case']:<{width}}  {r['driver']:<7}  "
          f"{r['resolves']:<13}  {r['idempotent']:<10}  {r['reproduces']}")
    detail = str(r.get("reproduces_detail") or "").strip()
    if detail and r["reproduces"] in {"mismatch", "run_failed", "unsupported_ref"}:
        first_line, *_ = detail.splitlines()
        print(f"{'':<{width}}  {'':<7}  {'':<13}  {'':<10}  detail: {first_line}")



def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-phase2", action="store_true")
    parser.add_argument("--stream", action="store_true", help="Print per-case outcomes live")
    args = parser.parse_args()

    width = max(len(c.case_dir) for c in REGRESSION_CASES) if REGRESSION_CASES else 50
    if args.stream:
        _print_header(width)

    rows = []
    for r in build_matrix_iter(run_phase2=args.run_phase2):
        if args.stream:
            _print_row(r, width)
        rows.append(r)

    if not args.stream:
        _print_header(width)
        for r in rows:
            _print_row(r, width)

    failed = [
        r for r in rows
        if r["resolves"] == "not-runnable"
        or r["idempotent"] == "DRIFT"
        or r["reproduces"] in {"mismatch", "run_failed"}
    ]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
