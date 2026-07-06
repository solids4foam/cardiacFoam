from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass
class RunStatus:
    ran: bool
    last_run: str | None
    time_dirs: list[str]
    regression: str            # available | none
    decomposed: bool


def _numeric_time_dirs(case_dir: Path) -> list[str]:
    out = []
    for child in case_dir.iterdir():
        if not child.is_dir():
            continue
        try:
            val = float(child.name)
        except ValueError:
            continue
        if val > 0.0:
            out.append(child.name)
    return sorted(out, key=float)


def _iso(ts: float) -> str:
    return datetime.fromtimestamp(ts, tz=timezone.utc).isoformat()


def _latest_mtime(case_dir: Path) -> float | None:
    candidates: list[float] = []
    for pattern in ("log.*",):
        for p in case_dir.glob(pattern):
            candidates.append(p.stat().st_mtime)
    pp = case_dir / "postProcessing"
    if pp.is_dir():
        for p in pp.rglob("*"):
            if p.is_file():
                candidates.append(p.stat().st_mtime)
    for name in _numeric_time_dirs(case_dir):
        candidates.append((case_dir / name).stat().st_mtime)
    return max(candidates) if candidates else None


def _regression(case_dir: Path) -> str:
    has_script = any((case_dir / n).is_file()
                     for n in ("regressionTest.sh", "runRegressionTest.sh"))
    has_ref = any(case_dir.glob("*.reference")) or any(case_dir.glob("*.ref"))
    return "available" if (has_script and has_ref) else "none"


def resolve_status(case_dir: Path) -> RunStatus:
    case_dir = Path(case_dir)
    time_dirs = _numeric_time_dirs(case_dir)
    pp = case_dir / "postProcessing"
    has_pp = pp.is_dir() and any(pp.rglob("*"))
    has_log = any(case_dir.glob("log.*"))
    ran = bool(time_dirs) or has_pp or has_log
    mtime = _latest_mtime(case_dir)
    return RunStatus(
        ran=ran,
        last_run=_iso(mtime) if (ran and mtime) else None,
        time_dirs=time_dirs,
        regression=_regression(case_dir),
        decomposed=(case_dir / "processor0").is_dir(),
    )
