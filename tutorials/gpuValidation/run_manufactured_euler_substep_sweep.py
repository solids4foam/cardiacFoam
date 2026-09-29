#!/usr/bin/env python3
"""Sweep fixed Euler substeps across the manufactured SBDF2 time ladder."""

from __future__ import annotations

import argparse
import csv
import math
import time
from pathlib import Path

from compare_manufactured_batched_parity import STATE_NAMES, read_states
from compare_slab import field
from run_manufactured_batched_temporal import (
    DT_DEFAULTS,
    FIELDS,
    NORMS,
    prepare_case,
    run_case,
)


SUBSTEP_FIELDS = ("Vm", "ionicCurrent", "activationTime", "u1", "u2", "u3")


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    columns = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in columns}
                         for row in rows)


def make_comparison_rows(root: Path, rows: list[dict[str, str]],
                         end_time: float) -> list[dict[str, str]]:
    completed = {
        (int(row["substeps"]), row["deltaT_s"]): row
        for row in rows if row.get("status") == "completed"
    }
    max_substeps = max(int(row["substeps"]) for row in rows)
    comparisons: list[dict[str, str]] = []

    for dt in sorted({row["deltaT_s"] for row in rows}, key=float, reverse=True):
        reference_row = completed.get((max_substeps, dt))
        if reference_row is None:
            continue
        reference_case = Path(reference_row["case_path"])
        for substeps in sorted({int(row["substeps"]) for row in rows}):
            row = completed.get((substeps, dt))
            if row is None:
                continue
            case = Path(row["case_path"])
            comparison = {
                "deltaT_s": dt,
                "substeps": str(substeps),
                "reference_substeps": str(max_substeps),
                "ionic_dt_s": row["ionic_dt_s"],
                "reference_case": reference_case.name,
            }

            reference_vm = field(reference_case, f"{end_time:g}", "Vm")
            for name in SUBSTEP_FIELDS:
                base = field(reference_case, f"{end_time:g}", name,
                             len(reference_vm))
                current = field(case, f"{end_time:g}", name, len(reference_vm))
                difference = current - base
                comparison[f"{name}_rmse_vs_{max_substeps}"] = (
                    f"{math.sqrt(float((difference*difference).mean())):.9e}"
                )
                comparison[f"{name}_max_vs_{max_substeps}"] = (
                    f"{float(abs(difference).max()):.9e}"
                )

            reference_states = read_states(reference_case, f"{end_time:g}")
            current_states = read_states(case, f"{end_time:g}")
            if current_states.shape != reference_states.shape:
                raise ValueError(f"state array shape mismatch: {case}")
            for state_i, state_name in enumerate(STATE_NAMES):
                difference = (current_states[:, state_i]
                              - reference_states[:, state_i])
                comparison[f"state_{state_name}_rmse_vs_{max_substeps}"] = (
                    f"{math.sqrt(float((difference*difference).mean())):.9e}"
                )
                comparison[f"state_{state_name}_max_vs_{max_substeps}"] = (
                    f"{float(abs(difference).max()):.9e}"
                )
            comparisons.append(comparison)

    return comparisons


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--delta-ts", nargs="+", type=float,
                        default=DT_DEFAULTS)
    parser.add_argument("--substeps", nargs="+", type=int,
                        default=(1, 5, 25, 100, 400, 800))
    parser.add_argument("--cells", type=int, default=640)
    parser.add_argument("--end-time", type=float, default=0.2)
    parser.add_argument("--require-gpu", action="store_true")
    args = parser.parse_args()
    if (not args.delta_ts or min(args.delta_ts) <= 0 or not args.substeps
            or min(args.substeps) < 1 or args.cells < 2 or args.end_time <= 0):
        parser.error("deltaT/endTime/substeps must be positive and cells >= 2")
    if not args.require_gpu:
        parser.error("this GPU substep sweep requires --require-gpu")

    root = args.output_dir.resolve()
    if root.exists():
        parser.error(f"output directory already exists: {root}")
    root.mkdir(parents=True)
    rows: list[dict[str, str]] = []

    for substeps in sorted(set(args.substeps)):
        substep_root = root / f"n{substeps}"
        substep_root.mkdir()
        for dt in sorted(set(args.delta_ts), reverse=True):
            case = substep_root / (
                f"gpu_sbdf2_N{args.cells}_n{substeps}_dt{dt:g}"
            )
            prepare_case(case, dt, args.cells, args.end_time, substeps, "sbdf2")
            started = time.monotonic()
            result = run_case(case, require_gpu=True)
            row = {
                "case": case.name,
                "case_path": str(case),
                "backend_requested": "gpu",
                "model": "monodomainFDAManufacturedBatched",
                "coupling": "sbdf2",
                "ddt_Vm": "backward",
                "deltaT_s": f"{dt:.16g}",
                "substeps": str(substeps),
                "substep_scaling": "fixed",
                "ionic_dt_s": f"{dt / substeps:.16g}",
                "cells_per_direction": str(args.cells),
                "cell_count": str(args.cells * args.cells),
                "end_time_s": f"{args.end_time:.16g}",
                "wall_s": f"{time.monotonic() - started:.6g}",
                **result,
            }
            rows.append(row)

            for field_name in FIELDS:
                for norm in NORMS:
                    key = f"{field_name}_{norm}"
                    row[f"{key}_order"] = ""
            previous = next((old for old in reversed(rows[:-1])
                             if old["substeps"] == str(substeps)
                             and old.get("status") == "completed"), None)
            if previous and row.get("status") == "completed":
                ratio = float(previous["deltaT_s"]) / dt
                for field_name in FIELDS:
                    for norm in NORMS:
                        key = f"{field_name}_{norm}"
                        try:
                            order = math.log(
                                float(previous[key]) / float(row[key])
                            ) / math.log(ratio)
                            row[f"{key}_order"] = f"{order:.6g}"
                        except (KeyError, ValueError, ZeroDivisionError):
                            pass
            write_csv(root / "summary.csv", rows)
            print(
                f"n={substeps} dt={dt:g}: {row.get('status')} "
                f"backend={row.get('backend')} "
                f"Vm_L2={row.get('Vm_L2', '')}",
                flush=True,
            )
            if row.get("log_tail"):
                print(row["log_tail"], flush=True)

    comparison_rows = make_comparison_rows(root, rows, args.end_time)
    write_csv(root / "substep_effect_vs_finest.csv", comparison_rows)
    completed = all(row.get("status") == "completed" for row in rows)
    (root / "run_status.txt").write_text(
        f"completed_cases={sum(row.get('status') == 'completed' for row in rows)}\n"
        f"total_cases={len(rows)}\n"
        f"all_completed={'yes' if completed else 'no'}\n"
        "Each temporal-order gate is descriptive here; fixed-substep Euler "
        "may reduce the observed SBDF2 order.\n"
    )
    return 0 if completed else 1


if __name__ == "__main__":
    raise SystemExit(main())
