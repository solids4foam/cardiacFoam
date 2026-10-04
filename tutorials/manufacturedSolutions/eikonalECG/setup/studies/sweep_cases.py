"""Read a completed omnidriver sweep of a manufacturedEikonalECG tet study."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterator

_GRAD_SCHEME_LABELS = {"Gauss linear": "GaussLinear", "leastSquares": "leastSquares"}
_GRAD_SCHEMES_DEFAULT = re.compile(r"gradSchemes\s*\{[^}]*?\bdefault\s+([^;]+);")
_NUMBER_OF_CELLS = re.compile(r"Number of cells\s*=\s*(\d+)")
_GRADIENT_REPORT = {
    "n_cells": r"^\s*cells\s*:\s*(\d+)",
    "linf_max": r"^\s*E_inf\s*=\s*([-+0-9.eE]+)",
    "linf_mean": r"^\s*Mean E_inf\s*=\s*([-+0-9.eE]+)",
    "n_above": r"^\s*Cells with error > 0\.05\s*=\s*(\d+)",
    "l2_bulk": r"^\s*L2 Bulk\s*=\s*([-+0-9.eE]+)",
    "l2_boundary": r"^\s*L2 Bound\s*=\s*([-+0-9.eE]+)",
    "l2_total": r"^\s*L2 Total\s*=\s*([-+0-9.eE]+)",
}


def completed_cases(sweep_dir: Path) -> Iterator[tuple[dict, Path]]:
    """Yield (resolved axis values, staged case directory) for each completed case.

    `<sweep_dir>/<caseId>/case_record.json` holds the axis values and
    `<sweep_dir>/cases/<caseId>/` the case as it ran.
    """
    for record_path in sorted(sweep_dir.glob("*/case_record.json")):
        record = json.loads(record_path.read_text())
        if record["status"] == "completed":
            yield record["resolved_axis_values"], sweep_dir / "cases" / record_path.parent.name


def grad_scheme(case_dir: Path) -> str:
    """The case's own `gradSchemes.default`, as the label the tet tables use."""
    match = _GRAD_SCHEMES_DEFAULT.search((case_dir / "system" / "fvSchemes").read_text())
    if match is None:
        raise ValueError(f"no gradSchemes default in {case_dir / 'system' / 'fvSchemes'}")
    return _GRAD_SCHEME_LABELS[match.group(1).strip()]


def measured_h(dat_text: str) -> str | None:
    """Mesh scale from the verifier's `Number of cells`: 1/round(cells^(1/3))."""
    match = _NUMBER_OF_CELLS.search(dat_text)
    if match is None:
        return None
    return f"{1.0 / max(1, int(int(match.group(1)) ** (1.0 / 3.0) + 0.5)):g}"


def gradient_reconstruction_report(case_dir: Path) -> dict[str, str] | None:
    """The values `gradientReconstructionOrder` printed in its last attempt's stdout, or None."""
    logs = sorted((case_dir / "workflow_logs").glob("gradientReconstructionOrder.attempt*.stdout.log"))
    if not logs:
        return None
    text = logs[-1].read_text(errors="ignore")
    report = {}
    for key, pattern in _GRADIENT_REPORT.items():
        match = re.search(pattern, text, re.MULTILINE)
        if match is None:
            return None
        report[key] = match.group(1)
    return report


def gradient_rows(sweep_dir: Path) -> list[dict[str, str]]:
    """One row per completed case that ran `gradientReconstructionOrder`."""
    rows = []
    for axes, case_dir in completed_cases(sweep_dir):
        report = gradient_reconstruction_report(case_dir)
        if report is None:
            continue
        n_cells = int(report["n_cells"])
        rows.append({
            "scheme": grad_scheme(case_dir),
            "N": str(int(axes["tetNumberCells"])),
            "h": f"{1.0 / max(1, round(n_cells ** (1.0 / 3.0))):g}",
            "n_cells": str(n_cells),
            "Linf_max": f"{float(report['linf_max']):g}",
            "Linf_mean": f"{float(report['linf_mean']):g}",
            "n_cells_Linf_gt_0_05": report["n_above"],
            "L2_bulk": f"{float(report['l2_bulk']):g}",
            "L2_boundary": f"{float(report['l2_boundary']):g}",
            "L2_total": f"{float(report['l2_total']):g}",
        })
    return rows
