from __future__ import annotations

from pathlib import Path
from typing import Iterator

from .assets3d import is_final_3d
from .case_scan import scan_cases
from .metrics import extract_metrics
from .models import CaseCard
from .plots import collect_plots
from .run_status import resolve_status

_ARTIFACT_GLOBS = {
    "foam": "*.foam",
    "vtu": "*.vtu",
    "vtk": "*.vtk",
}


def _artifacts(case_dir: Path, root: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for kind, pattern in _ARTIFACT_GLOBS.items():
        hits = sorted(case_dir.glob(pattern))
        if hits:
            out[kind] = hits[0].relative_to(root).as_posix()
    return out


def build_catalog(root: Path) -> Iterator[CaseCard]:
    root = Path(root)
    for card in scan_cases(root):
        case_dir = root / card.case_id
        st = resolve_status(case_dir)
        card.status = "ran" if st.ran else "not_run"
        card.last_run = st.last_run
        card.regression = st.regression
        card.metrics = extract_metrics(card.case_id, case_dir)
        card.plots = [p.relative_to(root).as_posix() for p in collect_plots(case_dir)]
        card.artifacts = _artifacts(case_dir, root)
        card.is_3d = is_final_3d(card.case_id)
        yield card
