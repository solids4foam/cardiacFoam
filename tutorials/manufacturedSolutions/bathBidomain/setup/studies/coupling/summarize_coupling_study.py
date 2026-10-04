"""Summarise the bath predictor-corrector control.

Reads an omnidriver sweep of sweep_coupling_study.json: each case's
`<sweep>/<caseId>/case_record.json` carries its resolved axis values, and its
metrics are `<sweep>/cases/<caseId>/postProcessing/bathBidomainInterfaceMetrics.csv`.
Writes raw_results.csv and summary.md into the sweep directory.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys


METRICS = (
    "heartPhiE_L2",
    "bathPhiE_L2",
    "x0FluxJump_L2",
    "x0IntracellularLeak_L2",
    "x0AssembledFlux_L2",
)

RESOLUTION_AXIS = "tetNumberCells"
PREDICTOR_AXIS = "constant/electroProperties:bidomainSolverCoeffs.bathPredictorCorrector"


def discover_rows(sweep_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record_path in sorted(sweep_dir.glob("*/case_record.json")):
        case_id = record_path.parent.name
        csv_path = sweep_dir / "cases" / case_id / "postProcessing" / "bathBidomainInterfaceMetrics.csv"
        if not csv_path.is_file():
            continue
        axes = json.loads(record_path.read_text())["resolved_axis_values"]
        with csv_path.open() as handle:
            values = next(csv.DictReader(handle))
        row = {
            "resolution": str(axes[RESOLUTION_AXIS]),
            "variant": "predictor" if axes[PREDICTOR_AXIS] else "baseline",
        }
        row.update({name: values[name] for name in METRICS})
        rows.append(row)
    return rows


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: summarize_coupling_study.py SWEEP_DIR")
    sweep_dir = Path(sys.argv[1]).resolve()
    rows = discover_rows(sweep_dir)
    if not rows:
        raise SystemExit(f"No completed coupling-study case found under {sweep_dir}")

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

    output_dir = sweep_dir

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
