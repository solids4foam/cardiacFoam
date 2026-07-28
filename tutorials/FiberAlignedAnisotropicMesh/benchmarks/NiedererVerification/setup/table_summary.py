"""table_summary.py — Activation time summary table for NiedererEtAl2011.

Reads *points_DT*_DX*.csv files from output_dir, extracts activation times
per probe point, and writes NiedererEtAl2011_summary.csv + .html via TableWriter.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

TUTORIALS_ROOT = Path(__file__).resolve().parents[4]
if str(TUTORIALS_ROOT) not in sys.path:
    sys.path.insert(0, str(TUTORIALS_ROOT))

from openfoam_driver.postprocessing.plotting_common import extract_dx_dt
from openfoam_driver.postprocessing.table_writer import TableMetadata, TableWriter


def _parse_filename(filename: str) -> tuple[str, float, float, str]:
    """Return (case_id, dx_mm, dt_ms, solver) from a points CSV filename.

    Expected pattern: ``{solver}_{model}_{tissue}_points_DT{tag}_DX{tag}.csv``
    """
    stem = filename.replace(".csv", "")
    dx, dt = extract_dx_dt(stem)
    case_id = stem.split("_points_DT")[0]
    solver = case_id.split("_")[0]
    return case_id, round(dx, 4), round(dt, 5), solver


def build_summary_rows(output_dir: Path) -> list[dict]:
    files = sorted(output_dir.glob("*points_DT*_DX*.csv"))
    if not files:
        print(f"[NiedererEtAl2011/table_summary] No points CSV files in {output_dir}")
        return []

    rows = []
    for fpath in files:
        case_id, dx, dt, solver = _parse_filename(fpath.name)
        df = pd.read_csv(fpath)
        activation_ms = df["activationTime"].values * 1000.0  # s → ms
        row: dict = {
            "case_id": case_id,
            "DX_mm": dx,
            "DT_ms": dt,
            "solver": solver,
        }
        for i, t in enumerate(activation_ms):
            row[f"point_{i}_activation_ms"] = round(float(t), 4)
        rows.append(row)

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
    folder = Path(__file__).resolve().parents[1] / "NiedererFoam"
    print(f"[table_summary] Default folder = {folder}")
    run_postprocessing(output_dir=str(folder))
