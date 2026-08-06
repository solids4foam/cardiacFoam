#!/usr/bin/env python3
"""Summarise the bath predictor-corrector control."""

from __future__ import annotations

import csv
from pathlib import Path
import sys


METRICS = (
    "heartPhiE_L2",
    "bathPhiE_L2",
    "x0FluxJump_L2",
    "x0IntracellularLeak_L2",
    "x0AssembledFlux_L2",
)


def metadata(path: Path) -> dict[str, str]:
    return dict(line.split("=", 1) for line in path.read_text().splitlines())


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: summarize_coupling_study.py RESULTS_DIR")
    root = Path(sys.argv[1]).resolve()
    rows: list[dict[str, str]] = []
    for meta_path in sorted(root.glob("*/metadata.env")):
        row = metadata(meta_path)
        with (meta_path.parent / "metrics.csv").open() as handle:
            values = next(csv.DictReader(handle))
        row.update({name: values[name] for name in METRICS})
        rows.append(row)
    if not rows:
        raise SystemExit(f"No study results found under {root}")

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

    fields = list(rows[0])
    with (root / "raw_results.csv").open("w", newline="") as handle:
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
    for row in rows:
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
    (root / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
