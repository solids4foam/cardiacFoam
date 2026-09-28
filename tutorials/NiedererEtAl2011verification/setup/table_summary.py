"""table_summary.py — Activation time summary table for NiedererEtAl2011.

Reads each swept case's own `postProcessing/Niedererpoints/<time>/activationTime`
(OpenFOAM's raw `probes` functionObject output) and writes
NiedererEtAl2011_summary.csv + .html via TableWriter.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from omnidriver.postprocessing.table_writer import TableMetadata, TableWriter

_PROBE_LINE = re.compile(r"^# Probe (\d+) \(([^)]+)\)")

#: `probes` never activates every point at a short/coarse verification
#: endTime; this repository's sentinel for "never activated" -- never
#: converted s -> ms (CLAUDE.md: "a sentinel is never converted").
_UNACTIVATED_SENTINEL = -1.0


def _read_latest_probe_row(function_object_dir: Path) -> list[float] | None:
    """Return the last data row of a `probes` functionObject's output file.

    `postProcess -latestTime` writes exactly one data row per probe sample
    call; the instance directory it lands under is named for when that
    postProcess invocation itself started (often "0"), not for the
    simulated time in the row, so every time subdirectory must be checked
    and the row's own leading Time column used instead.
    """
    if not function_object_dir.is_dir():
        return None
    candidates = [p for p in function_object_dir.iterdir() if p.is_dir()]
    for instance_dir in sorted(candidates):
        sample_path = instance_dir / "activationTime"
        if not sample_path.is_file():
            continue
        last_row: list[float] | None = None
        for line in sample_path.read_text().splitlines():
            if _PROBE_LINE.match(line) or not line.strip() or line.startswith("#"):
                continue
            last_row = [float(token) for token in line.split()]
        if last_row is not None:
            return last_row[1:]  # drop the leading Time column
    return None


def _iter_swept_cases(output_dir: Path):
    """Yield (case_id, case_dir, resolved_axis_values) for every case in this sweep.

    A tutorial-record case's own case_record.json never carries sweep axis
    values (Core's standalone case record does not scan for them); the
    sweep's own sweep_manifest.json, written once per sweep in output_dir,
    is the authoritative source.
    """
    manifest_path = output_dir / "sweep_manifest.json"
    if not manifest_path.is_file():
        return
    manifest = json.loads(manifest_path.read_text())
    for case in manifest.get("cases", []):
        case_dir = output_dir / "cases" / case["case_id"]
        if case_dir.is_dir():
            yield case["case_id"], case_dir, case.get("resolved_axis_values", {})


def _dx_dt_solver(resolved_axis_values: dict) -> tuple[float, float, str]:
    """Return (DX_mm, DT_ms, solver) from a case's resolved sweep axes.

    `dx` (hex mesh) or `tetDx` (tet mesh) is in metres; `system/controlDict:
    deltaT` is in seconds. Both studies fix `solutionAlgorithm implicit` (the
    only value NiedererEtAl2011verification's studies sweep today).
    """
    dx_m = resolved_axis_values.get("dx", resolved_axis_values.get("tetDx"))
    dt_s = resolved_axis_values.get("system/controlDict:deltaT")
    dx_mm = round(float(dx_m) * 1000.0, 4) if dx_m is not None else float("nan")
    dt_ms = round(float(dt_s) * 1000.0, 5) if dt_s is not None else float("nan")
    return dx_mm, dt_ms, "implicit"


def build_summary_rows(output_dir: Path) -> list[dict]:
    rows: list[dict] = []
    for case_id, case_dir, resolved_axis_values in _iter_swept_cases(output_dir):
        values = _read_latest_probe_row(case_dir / "postProcessing" / "Niedererpoints")
        if values is None:
            continue
        dx_mm, dt_ms, solver = _dx_dt_solver(resolved_axis_values)
        row: dict = {"case_id": case_id, "DX_mm": dx_mm, "DT_ms": dt_ms, "solver": solver}
        for i, value in enumerate(values):
            activation_ms = value if value == _UNACTIVATED_SENTINEL else round(value * 1000.0, 4)
            row[f"point_{i}_activation_ms"] = activation_ms
        rows.append(row)
    if not rows:
        print(f"[NiedererEtAl2011/table_summary] No swept cases with Niedererpoints output in {output_dir}")
    return rows


def run_postprocessing(
    *, output_dir: str, setup_root: str | None = None, **_: object
) -> list[dict]:
    output_path = Path(output_dir)
    rows = build_summary_rows(output_path)
    if not rows:
        return []
    meta = TableMetadata(
        tutorial="NiedererEtAl2011",
        units={"activationTime": "ms", "DX": "mm", "DT": "ms"},
    )
    return TableWriter.write(
        rows,
        output_path,
        "NiedererEtAl2011_summary",
        "Niederer activation time summary",
        meta,
    )


if __name__ == "__main__":
    folder = Path(__file__).resolve().parents[1]
    print(f"[table_summary] Default folder = {folder}")
    run_postprocessing(output_dir=str(folder))
