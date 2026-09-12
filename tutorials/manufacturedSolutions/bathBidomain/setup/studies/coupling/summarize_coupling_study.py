#!/usr/bin/env python3
"""Summarise the bath predictor-corrector control.

Reads driverFOAM sweep-run's own archive layout for
setup/studies/coupling/sweep_coupling_study.json: each case lands at
<case_root>/<caseId>/<archive_dir_name>/bathBidomainInterfaceMetrics.csv,
where <archive_dir_name> is the spec's own
"setup/studies/coupling/results/sweepCases" and <caseId> is
"<number_cells>_<bath_predictor_corrector>" (e.g. "10_False", "10_True"),
per the spec's case_id_template derive. Verified against a real N=10
baseline/predictor run (2026-08-19).
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
import sys


METRICS = (
    "heartPhiE_L2",
    "bathPhiE_L2",
    "x0FluxJump_L2",
    "x0IntracellularLeak_L2",
    "x0AssembledFlux_L2",
)

ARCHIVE_RELPATH = Path("setup/studies/coupling/results/sweepCases")
CASE_ID_RE = re.compile(r"^(?P<resolution>\d+)_(?P<predictor>True|False)$")


def discover_rows(case_root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for csv_path in sorted(case_root.glob(f"*/{ARCHIVE_RELPATH}/bathBidomainInterfaceMetrics.csv")):
        case_id = csv_path.parents[len(ARCHIVE_RELPATH.parts)].name
        match = CASE_ID_RE.match(case_id)
        if match is None:
            continue
        with csv_path.open() as handle:
            values = next(csv.DictReader(handle))
        row = {
            "resolution": match["resolution"],
            "variant": "predictor" if match["predictor"] == "True" else "baseline",
        }
        row.update({name: values[name] for name in METRICS})
        rows.append(row)
    return rows


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: summarize_coupling_study.py CASE_ROOT")
    case_root = Path(sys.argv[1]).resolve()
    rows = discover_rows(case_root)
    if not rows:
        raise SystemExit(
            f"No sweep-run case output found under {case_root}/*/{ARCHIVE_RELPATH}/ "
            "-- run 'driverFoam sweep-run --spec setup/studies/coupling/sweep_coupling_study.json' first."
        )

    baselines = {
        row["resolution"]: row for row in rows if row["variant"] == "baseline"
    }
    for row in rows:
        baseline = baselines.get(row["resolution"])
        if baseline is None:
            raise SystemExit(f"Missing baseline for N={row['resolution']}")
        for name in METRICS:
            reference = float(baseline[name])
            row[f"{name}_change_percent"] = str(
                100.0 * (float(row[name]) - reference) / reference
            )

    output_dir = case_root / "setup" / "studies" / "coupling" / "results"
    output_dir.mkdir(parents=True, exist_ok=True)

    fields = list(rows[0])
    with (output_dir / "raw_results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Bath coupling-control summary",
        "",
        "Percent changes are relative to `bathPredictorCorrector false` on "
        "the same mesh and with the same shared non-orthogonal count.",
        "",
        "| N | variant | heart phiE L2 | change | bath phiE L2 | change | "
        "assembled current L2 | change |",
        "|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in sorted(rows, key=lambda r: (int(r["resolution"]), r["variant"])):
        values: dict[str, object] = dict(row)
        values.update(
            heart_value=float(row["heartPhiE_L2"]),
            heart_change=float(row["heartPhiE_L2_change_percent"]),
            bath_value=float(row["bathPhiE_L2"]),
            bath_change=float(row["bathPhiE_L2_change_percent"]),
            flux_value=float(row["x0AssembledFlux_L2"]),
            flux_change=float(row["x0AssembledFlux_L2_change_percent"]),
        )
        lines.append(
            "| {resolution} | {variant} | {heart_value:.6e} | "
            "{heart_change:+.2f}% | {bath_value:.6e} | {bath_change:+.2f}% | "
            "{flux_value:.6e} | {flux_change:+.2f}% |".format(**values)
        )
    (output_dir / "summary.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {output_dir / 'raw_results.csv'} and {output_dir / 'summary.md'}")


if __name__ == "__main__":
    main()
