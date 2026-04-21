"""table_summary.py — Consolidated restitution table for restitutionCurves tutorial.

Reads per-model ``*_restitution.csv`` files produced by postProcessing_restCurves.py,
merges them into a single table with an ``ionic_model`` column, and writes
restitutionCurves_summary.csv + .html via TableWriter.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

TUTORIALS_ROOT = Path(__file__).resolve().parents[4]
if str(TUTORIALS_ROOT) not in sys.path:
    sys.path.insert(0, str(TUTORIALS_ROOT))

from openfoam_driver.postprocessing.table_writer import TableMetadata, TableWriter


def run_postprocessing(
    *, output_dir: str, setup_root: str | None = None, **_: object
) -> list[dict]:
    output_path = Path(output_dir)
    csv_files = sorted(output_path.glob("*_restitution.csv"))
    if not csv_files:
        print(f"[restitutionCurves/table_summary] No *_restitution.csv in {output_path}")
        return []

    rows: list[dict] = []
    for fpath in csv_files:
        ionic_model = fpath.stem.replace("_restitution", "")
        try:
            df = pd.read_csv(fpath)
        except Exception as exc:
            print(f"[restitutionCurves/table_summary] Could not read {fpath.name}: {exc}")
            continue
        for _, row in df.iterrows():
            rows.append({
                "ionic_model": ionic_model,
                "tissue": str(row.get("tissue", "unknown")),
                "DI_ms": round(float(row.get("DI_ms", 0.0)), 3),
                "APD_ms": round(float(row.get("APD90_ms", 0.0)), 3),
            })

    if not rows:
        return []

    meta = TableMetadata(
        tutorial="restitutionCurves_s1s2Protocol",
        units={"DI_ms": "ms", "APD_ms": "ms"},
    )
    return TableWriter.write(
        rows,
        output_path,
        "restitutionCurves_summary",
        "Restitution curves consolidated summary",
        meta,
    )


if __name__ == "__main__":
    folder = Path(__file__).resolve().parents[2]
    print(f"[table_summary] Default folder = {folder}")
    run_postprocessing(output_dir=str(folder))
