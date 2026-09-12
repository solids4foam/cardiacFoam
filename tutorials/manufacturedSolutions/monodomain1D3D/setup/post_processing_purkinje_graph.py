#!/usr/bin/env python3
"""Post-process manufactured Purkinje graph convergence sweeps."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


ERROR_FILE_RE = re.compile(r"graph_(?P<dimension>[^_]+)_(?P<nodes>\d+)_nodes\.dat$")
FIELD_RE = re.compile(
    r"^(?P<field>Vm1D|u1|u2)\s+"
    r"(?P<l1>[-+0-9.eE]+)\s+"
    r"(?P<l2>[-+0-9.eE]+)\s+"
    r"(?P<linf>[-+0-9.eE]+)\s*$"
)

SUMMARY_FIELDS = (
    "Graph",
    "Dimension",
    "Nodes",
    "SegmentsPerBranch",
    "h",
    "L1_Vm",
    "L2_Vm",
    "Linf_Vm",
    "L1_u1",
    "L2_u1",
    "Linf_u1",
    "L1_u2",
    "L2_u2",
    "Linf_u2",
)

RATE_FIELDS = (
    "Dimension",
    "Nodes_lower",
    "Nodes_higher",
    "rate_Vm",
    "rate_u1",
    "rate_u2",
)


def _safe_rate(error_lower: float, error_higher: float, h_lower: float, h_higher: float) -> str:
    if error_lower <= 0.0 or error_higher <= 0.0 or h_lower <= 0.0 or h_higher <= 0.0:
        return ""
    return f"{math.log(error_lower / error_higher) / math.log(h_lower / h_higher):.12g}"


def _write_csv(rows: list[dict[str, object]], destination: Path, fieldnames: tuple[str, ...]) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_error_file(path: Path) -> dict[str, object]:
    match = ERROR_FILE_RE.match(path.name)
    if not match:
        raise ValueError(f"Unexpected graph verifier filename: {path}")

    nodes = int(match.group("nodes"))
    if nodes < 3 or (nodes - 1) % 2 != 0:
        raise ValueError(f"Expected Y-graph node count 2*segments+1, got {nodes}")

    segments_per_branch = (nodes - 1) // 2
    row: dict[str, object] = {
        "Graph": path.parent.name,
        "Dimension": match.group("dimension"),
        "Nodes": nodes,
        "SegmentsPerBranch": segments_per_branch,
        "h": f"{1.0 / segments_per_branch:.12g}",
    }

    for line in path.read_text(encoding="utf-8").splitlines():
        field_match = FIELD_RE.match(line.strip())
        if not field_match:
            continue

        field = "Vm" if field_match.group("field") == "Vm1D" else field_match.group("field")
        row[f"L1_{field}"] = float(field_match.group("l1"))
        row[f"L2_{field}"] = float(field_match.group("l2"))
        row[f"Linf_{field}"] = float(field_match.group("linf"))

    missing = [name for name in SUMMARY_FIELDS if name not in row]
    if missing:
        raise ValueError(f"{path} is missing expected columns: {', '.join(missing)}")

    return row


def collect_rows(output_dir: Path, graph_ids: list[str] | None = None) -> list[dict[str, object]]:
    rows = []
    patterns = (
        [f"{graph_id}/graph_*_nodes.dat" for graph_id in graph_ids]
        if graph_ids
        else ["nodes*/graph_*_nodes.dat"]
    )
    for pattern in patterns:
        for path in sorted(output_dir.glob(pattern)):
            rows.append(parse_error_file(path))
    return sorted(rows, key=lambda row: (str(row["Dimension"]), int(row["Nodes"])))


def compute_rates(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["Dimension"]), []).append(row)

    rates = []
    for dimension, group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: int(row["Nodes"]))
        for lower, higher in zip(ordered, ordered[1:]):
            rates.append(
                {
                    "Dimension": dimension,
                    "Nodes_lower": lower["Nodes"],
                    "Nodes_higher": higher["Nodes"],
                    "rate_Vm": _safe_rate(
                        float(lower["Linf_Vm"]),
                        float(higher["Linf_Vm"]),
                        float(lower["h"]),
                        float(higher["h"]),
                    ),
                    "rate_u1": _safe_rate(
                        float(lower["Linf_u1"]),
                        float(higher["Linf_u1"]),
                        float(lower["h"]),
                        float(higher["h"]),
                    ),
                    "rate_u2": _safe_rate(
                        float(lower["Linf_u2"]),
                        float(higher["Linf_u2"]),
                        float(lower["h"]),
                        float(higher["h"]),
                    ),
                }
            )

    return rates


def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    graph_ids: list[str] | tuple[str, ...] | None = None,
    **_: object,
) -> list[dict]:
    del setup_root
    output_dir = Path(output_dir).resolve()
    rows = collect_rows(output_dir, list(graph_ids) if graph_ids else None)
    if not rows:
        print(f"No graph verifier files found under {output_dir}")
        return []

    rates = compute_rates(rows)
    summary_csv = output_dir / "graph_convergence_summary.csv"
    rates_csv = output_dir / "graph_convergence_rates.csv"

    _write_csv(rows, summary_csv, SUMMARY_FIELDS)
    _write_csv(rates, rates_csv, RATE_FIELDS)

    print(f"Wrote {summary_csv}")
    print(f"Wrote {rates_csv}")

    return [
        {
            "path": str(summary_csv),
            "label": "Manufactured Purkinje graph convergence summary",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(rates_csv),
            "label": "Manufactured Purkinje graph convergence rates",
            "kind": "table",
            "format": "csv",
        },
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs" / "1dGraphConvergence",
        help="Directory containing nodes*/graph_*_nodes.dat files",
    )
    args = parser.parse_args()

    run_postprocessing(output_dir=str(args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
