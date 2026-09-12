#!/usr/bin/env python3
"""Regenerate eikonalECG's gauss-linear-vs-leastSquares gradient comparison CSV.

Reads this study's driverFOAM sweep archive (results/sweepCases, produced by
`driverFoam sweep-run --spec sweep_gradient_tet.json`, see README.md) via
adapters.from_eikonal_gradient_reconstruction, and writes
results/eikonal_gradient_tet.csv -- the full leastSquares-vs-gaussLinear
comparison used in the paper's qualitative discussion.

For the registered "eikonal_gradient_tet" Paper I table (leastSquares only),
see ../gradient_reconstruction/aggregate_gradient_reconstruction.py instead.

Usage (run after the sweep above):
    python3 aggregate_gradient_verification.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

_STUDY_DIR = Path(__file__).resolve().parent               # .../gradientVerification
_REPO_ROOT = _STUDY_DIR.parents[5]                          # repo root
sys.path.insert(0, str(_REPO_ROOT / "applications/scripts/paperI_results"))
import adapters  # noqa: E402

_FIELDS = [
    "scheme", "N", "h", "n_cells",
    "Linf_max", "Linf_mean", "n_cells_Linf_gt_0_05",
    "L2_bulk", "L2_boundary", "L2_total",
]


def main() -> None:
    sweep_cases = _STUDY_DIR / "results" / "sweepCases"
    manifest = _STUDY_DIR / "results" / "sweepRun" / "sweep_manifest.json"
    rows = adapters.from_eikonal_gradient_reconstruction(sweep_cases, manifest)
    rows.sort(key=lambda r: (r["scheme"], int(r["N"])))

    dest = _STUDY_DIR / "results" / "eikonal_gradient_tet.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in _FIELDS})
    print(f"Wrote {dest} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
