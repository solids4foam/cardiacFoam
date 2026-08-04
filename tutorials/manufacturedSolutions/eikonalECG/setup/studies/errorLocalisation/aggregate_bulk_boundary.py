#!/usr/bin/env python3
"""Regenerate eikonalECG's canonical solved-field bulk/boundary CSV.

Reads this study's own driverFOAM sweep archive (results/sweepCases +
results/sweepRun/sweep_manifest.json, produced by run_bulk_boundary_tet.sh)
via adapters.from_eikonal_bulk_boundary, and writes
setup/results/eikonal_bulk_boundary_tet.csv.

Naming follows applications/scripts/paperI_results/README.md's
<physics>_<operator>_<mesh> convention for isolated operator diagnostics
(e.g. eikonal_gradient_tet), which this decomposition is a sibling of: this
one is the SOLVED field's bulk/boundary split (from
manufacturedEikonalVerifier.C, only computed when writeErrorField is
enabled), not the standalone gradient-operator-only split gradientReconstru
ctionOrder / eikonal_gradient_tet.csv reports.

Not a keyset-gated convergence table (no committed reference, no 'field'/
rate columns), so this writes its own bespoke CSV shape directly rather than
going through schema.write_canonical -- see adapters.from_eikonal_bulk_boundary's
docstring.

Usage (normally invoked by run_bulk_boundary_tet.sh after the sweep):
    python3 aggregate_bulk_boundary.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

_STUDY_DIR = Path(__file__).resolve().parent               # .../errorLocalisation
_CASE_DIR = _STUDY_DIR.parents[2]                           # .../eikonalECG
_REPO_ROOT = _STUDY_DIR.parents[5]                          # repo root
sys.path.insert(0, str(_REPO_ROOT / "applications/scripts/paperI_results"))
import adapters  # noqa: E402

_FIELDS = [
    "case", "variant", "dim", "N", "h",
    "L2_bulk", "L2_boundary", "L2_total", "boundary_energy_fraction",
]


def main() -> None:
    sweep_cases = _STUDY_DIR / "results" / "sweepCases"
    manifest = _STUDY_DIR / "results" / "sweepRun" / "sweep_manifest.json"
    rows = adapters.from_eikonal_bulk_boundary(sweep_cases, manifest)
    rows.sort(key=lambda r: (r["variant"], -float(r["h"])))

    dest = _CASE_DIR / "setup" / "results" / "eikonal_bulk_boundary_tet.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in _FIELDS})
    print(f"wrote {len(rows)} rows -> {dest}")


if __name__ == "__main__":
    main()
