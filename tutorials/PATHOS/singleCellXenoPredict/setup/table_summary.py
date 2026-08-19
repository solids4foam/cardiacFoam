"""table_summary.py — Voltage and APD summary table for singleCell tutorial.

Reads .txt simulation output files from output_dir.  Each file contains a
space-separated time series with columns: time, Vm, [additional state vars...].
Extracts resting voltage, peak voltage, and APD at 90% repolarisation.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

TUTORIALS_ROOT = Path(__file__).resolve().parents[4]
if str(TUTORIALS_ROOT) not in sys.path:
    sys.path.insert(0, str(TUTORIALS_ROOT))

from openfoam_driver.postprocessing.table_writer import TableMetadata, TableWriter


def _parse_model_and_cell(stem: str) -> tuple[str, str]:
    """Extract ionic model and tissue from a filename stem like TNNP_epicardialCells_run."""
    parts = stem.split("_")
    model = parts[0] if parts else "unknown"
    cell = parts[1] if len(parts) > 1 else "unknown"
    return model, cell

# 10% of excursion above resting remaining = 90% repolarisation
_APD_REPOL_FRACTION = 0.10


def _compute_apd90(time: np.ndarray, vm: np.ndarray) -> float | None:
    """Return APD90 in ms, or None if detection fails."""
    resting = float(vm[0])
    peak_idx = int(np.argmax(vm))
    peak = float(vm[peak_idx])
    if peak <= resting:
        return None
    threshold = resting + _APD_REPOL_FRACTION * (peak - resting)
    for i in range(peak_idx + 1, len(vm)):
        if vm[i] <= threshold:
            # Linear interpolation for sub-sample accuracy
            frac = (threshold - float(vm[i - 1])) / (float(vm[i]) - float(vm[i - 1]))
            t_repol = float(time[i - 1]) + frac * (float(time[i]) - float(time[i - 1]))
            return (t_repol - float(time[peak_idx])) * 1000.0  # s → ms
    return None


def _extract_row(fpath: Path) -> dict | None:
    try:
        df = pd.read_csv(fpath, sep=r"\s+", comment="#", engine="python")
    except Exception as exc:
        print(f"[singleCell/table_summary] Could not read {fpath.name}: {exc}")
        return None
    if df.empty or len(df.columns) < 2:
        return None

    time = df.iloc[:, 0].values.astype(float)
    vm = df.iloc[:, 1].values.astype(float)
    apd = _compute_apd90(time, vm)
    model, cell = _parse_model_and_cell(fpath.stem)

    return {
        "case_id": fpath.stem,
        "ionic_model": model or "unknown",
        "tissue": cell or "unknown",
        "APD_ms": round(apd, 3) if apd is not None else "N/A",
        "peak_voltage_mV": round(float(np.max(vm)), 3),
        "resting_voltage_mV": round(float(vm[0]), 3),
    }


def run_postprocessing(
    *, output_dir: str, setup_root: str | None = None, **_: object
) -> list[dict]:
    output_path = Path(output_dir)
    txt_files = sorted(output_path.glob("*.txt"))
    if not txt_files:
        print(f"[singleCell/table_summary] No .txt files found in {output_path}")
        return []

    rows = [r for f in txt_files if (r := _extract_row(f)) is not None]
    if not rows:
        return []

    meta = TableMetadata(
        tutorial="singleCell",
        units={"APD_ms": "ms", "peak_voltage_mV": "mV", "resting_voltage_mV": "mV"},
    )
    return TableWriter.write(
        rows,
        output_path,
        "singleCell_summary",
        "Single cell voltage and APD summary",
        meta,
    )


if __name__ == "__main__":
    folder = Path(__file__).resolve().parents[2]
    print(f"[table_summary] Default folder = {folder}")
    run_postprocessing(output_dir=str(folder))
