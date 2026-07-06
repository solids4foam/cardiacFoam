from __future__ import annotations

from pathlib import Path

import numpy as np

_ECG_FILES = (
    "pseudoECG.dat",
    "manufacturedPseudoECG.dat",
    "eikonalECG.dat",
    "manufacturedEikonalECG.dat",
)


def load_dat(path: Path) -> tuple[list[str], np.ndarray]:
    """Parse an OpenFOAM `# col col ...` + numeric-rows .dat file.

    Raises ValueError on ragged/unparseable numeric content.
    """
    cols: list[str] = []
    rows: list[list[float]] = []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            if not cols:
                cols = s.lstrip("#").split()
            continue
        parts = s.split()
        rows.append([float(x) for x in parts])  # ValueError propagates
    if not rows:
        raise ValueError(f"no data rows in {path}")
    width = len(rows[0])
    if any(len(r) != width for r in rows):
        raise ValueError(f"ragged rows in {path}")
    arr = np.asarray(rows, dtype=float)
    if not cols or len(cols) != width:
        cols = [f"c{i}" for i in range(width)]
    return cols, arr


def _ecg_metrics(dat: Path) -> dict[str, str | float]:
    cols, arr = load_dat(dat)
    n = int(arr.shape[0])
    t_end = float(arr[-1][0])
    best_lead, best_pp = None, -1.0
    for j in range(1, arr.shape[1]):
        pp = float(arr[:, j].max() - arr[:, j].min())
        if pp > best_pp:
            best_pp, best_lead = pp, cols[j]
    return {
        "n_samples": n,
        "t_end": t_end,
        "ecg_peak_to_peak": round(best_pp, 6),
        "ecg_lead": best_lead or "n/a",
    }


def extract_metrics(case_id: str, case_dir: Path) -> dict[str, str | float]:
    case_dir = Path(case_dir)
    pp = case_dir / "postProcessing"
    if not pp.is_dir():
        return {}
    for name in _ECG_FILES:
        hits = list(pp.rglob(name))
        if hits:
            try:
                return _ecg_metrics(hits[0])
            except (ValueError, OSError):
                return {}
    return {}
