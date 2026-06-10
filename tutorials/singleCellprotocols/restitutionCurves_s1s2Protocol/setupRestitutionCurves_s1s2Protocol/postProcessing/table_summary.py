"""table_summary.py — Restitution curve summary table for restitutionCurves_s1s2Protocol.

Reads per-model *_restitution.csv files from output_dir and consolidates them
into a single restitutionCurves_summary.csv / .html.  The ionic model name is
derived from the filename stem (e.g. TNNP_restitution.csv → ionic_model=TNNP).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

TUTORIALS_ROOT = Path(__file__).resolve().parents[4]
if str(TUTORIALS_ROOT) not in sys.path:
    sys.path.insert(0, str(TUTORIALS_ROOT))

from openfoam_driver.postprocessing.table_writer import TableMetadata, TableWriter

_TUTORIAL_NAME = "restitutionCurves_s1s2Protocol"


def _model_from_stem(stem: str) -> str:
    """Extract ionic model name from a stem like 'TNNP_restitution'."""
    return stem.split("_")[0] if "_" in stem else stem


def run_postprocessing(
    *, output_dir: str, setup_root: str | None = None, **_: object
) -> list[dict]:
    output_path = Path(output_dir)
    csv_files = sorted(output_path.glob("*_restitution.csv"))
    if not csv_files:
        print(f"[{_TUTORIAL_NAME}/table_summary] No *_restitution.csv files found in {output_path}")
        return []

    rows: list[dict] = []
    for fpath in csv_files:
        try:
            df = pd.read_csv(fpath)
        except Exception as exc:
            print(f"[{_TUTORIAL_NAME}/table_summary] Could not read {fpath.name}: {exc}")
            continue
        model = _model_from_stem(fpath.stem)
        df.insert(0, "ionic_model", model)
        rows.extend(df.to_dict(orient="records"))

    if not rows:
        return []

    meta = TableMetadata(
        tutorial=_TUTORIAL_NAME,
        units={"DI_ms": "ms", "APD90_ms": "ms"},
    )
    return TableWriter.write(
        rows,
        output_path,
        "restitutionCurves_summary",
        "Restitution curve APD90 summary",
        meta,
    )


if __name__ == "__main__":
    folder = Path(__file__).resolve().parents[2]
    print(f"[table_summary] Default folder = {folder}")
    run_postprocessing(output_dir=str(folder))
