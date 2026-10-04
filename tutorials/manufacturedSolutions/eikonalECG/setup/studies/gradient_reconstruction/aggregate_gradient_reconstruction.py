#!/usr/bin/env python3
"""Regenerate the registered "eikonal_gradient_tet" Paper I table.

Reads an omnidriver sweep of sweep_gradient_tet.json: each case's
`workflow_logs/gradientReconstructionOrder.attempt*.stdout.log` is the report
of `gradientReconstructionOrder` (`applications/test/gradientReconstructionOrder/`).
Writes eikonal_gradient_tet.csv into the sweep directory.

Isolated `leastSquares` gradient reconstruction on the tet mesh, N = 10, 20, 40, 80.

Usage:
    python3 aggregate_gradient_reconstruction.py SWEEP_DIR
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import sweep_cases  # noqa: E402

_FIELDS = [
    "scheme", "N", "h", "n_cells",
    "Linf_max", "Linf_mean", "n_cells_Linf_gt_0_05",
    "L2_bulk", "L2_boundary", "L2_total",
]


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: aggregate_gradient_reconstruction.py SWEEP_DIR")
    sweep_dir = Path(sys.argv[1]).resolve()
    rows = sweep_cases.gradient_rows(sweep_dir)
    if not rows:
        raise SystemExit(f"no completed gradientReconstructionOrder case under {sweep_dir}")
    rows.sort(key=lambda r: int(r["N"]))

    dest = sweep_dir / "eikonal_gradient_tet.csv"
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {dest} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
