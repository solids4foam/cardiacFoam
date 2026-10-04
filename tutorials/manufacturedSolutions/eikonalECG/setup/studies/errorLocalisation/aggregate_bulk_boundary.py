#!/usr/bin/env python3
"""Regenerate eikonalECG's solved-field bulk/boundary CSV.

Reads an omnidriver sweep of sweep_tet_error_localisation.json: each case's
`postProcessing/manufacturedEikonalActivationTime.dat` carries the
`activationTimeSplit` line the verifier writes when `writeErrorField` is on.
Writes eikonal_bulk_boundary_tet.csv into the sweep directory.

This is the SOLVED field's bulk/boundary split, not the gradient-operator-only
split of ../gradientVerification/.

Usage:
    python3 aggregate_bulk_boundary.py SWEEP_DIR
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import sweep_cases  # noqa: E402

_SPLIT = re.compile(r"^activationTimeSplit\s+(\S+)\s+(\S+)\s+(\S+)\s*$", re.MULTILINE)
_FIELDS = [
    "case", "variant", "dim", "N", "h",
    "L2_bulk", "L2_boundary", "L2_total", "boundary_energy_fraction",
]


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: aggregate_bulk_boundary.py SWEEP_DIR")
    sweep_dir = Path(sys.argv[1]).resolve()
    rows = []
    for axes, case_dir in sweep_cases.completed_cases(sweep_dir):
        dat = case_dir / "postProcessing" / "manufacturedEikonalActivationTime.dat"
        if not dat.is_file():
            continue
        text = dat.read_text(errors="ignore")
        split = _SPLIT.search(text)
        if split is None:
            continue
        bulk, boundary, total = (float(value) for value in split.groups())
        n = int(axes["tetNumberCells"])
        rows.append({
            "case": "eikonal_tet_split", "variant": sweep_cases.grad_scheme(case_dir),
            "dim": "3D", "N": str(n), "h": sweep_cases.measured_h(text) or f"{1.0 / n:g}",
            "L2_bulk": f"{bulk:g}", "L2_boundary": f"{boundary:g}", "L2_total": f"{total:g}",
            "boundary_energy_fraction": f"{(boundary / total) ** 2 if total else 0.0:g}",
        })
    if not rows:
        raise SystemExit(f"no completed case with an activationTimeSplit under {sweep_dir}")
    rows.sort(key=lambda r: (r["variant"], -float(r["h"])))

    dest = sweep_dir / "eikonal_bulk_boundary_tet.csv"
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows -> {dest}")


if __name__ == "__main__":
    main()
