from __future__ import annotations

import re
from pathlib import Path
from typing import Iterator

from .models import CaseCard

FAMILIES = (
    "electrophysiologyProtocols",
    "manufacturedSolutions",
    "NiedererEtAl2011",
    "PATHOS",
    "heartSim3D-1D",
    "template",
)
_EXCLUDE_DIR = re.compile(r"(^|/)(processor\d+|setup|\.venv|build|\.git)(/|$)")


def _is_case(d: Path) -> bool:
    if not (d / "constant").is_dir():
        return False
    if (d / "system" / "controlDict").is_file():
        return True
    return any(d.glob("*.foam"))


def _family_of(rel: str) -> str:
    head = rel.split("/", 1)[0]
    return head if head in FAMILIES else "other"


def _goal_from_readme(case_dir: Path) -> str | None:
    readme = case_dir / "README.md"
    if not readme.is_file():
        return None
    lines = readme.read_text(errors="replace").splitlines()
    h1 = None
    body: list[str] = []
    seen_h1 = False
    for line in lines:
        if line.startswith("# ") and not seen_h1:
            h1 = line[2:].strip()
            seen_h1 = True
            continue
        if seen_h1:
            if line.strip() == "" and not body:
                continue
            if line.startswith("#"):
                break
            if line.strip() == "" and body:
                break
            body.append(line.strip())
    if body:
        return " ".join(body).strip()
    return h1


def _dict_lookup(case_dir: Path, relpath: str, key: str) -> str | None:
    f = case_dir / relpath
    if not f.is_file():
        return None
    try:
        text = f.read_text(errors="replace")
    except OSError:
        return None
    m = re.search(rf"^\s*{re.escape(key)}\s+([A-Za-z0-9_./+-]+)\s*;", text, re.MULTILINE)
    return m.group(1) if m else None


def scan_cases(root: Path) -> Iterator[CaseCard]:
    root = Path(root)
    if not root.is_dir():
        return
    for constant in sorted(root.rglob("constant")):
        case_dir = constant.parent
        rel = case_dir.relative_to(root).as_posix()
        if _EXCLUDE_DIR.search(rel + "/"):
            continue
        if not _is_case(case_dir):
            continue
        yield CaseCard(
            case_id=rel,
            family=_family_of(rel),
            name=case_dir.name,
            goal=_goal_from_readme(case_dir),
            solver=_dict_lookup(case_dir, "constant/electroProperties", "myocardiumSolver"),
            ionic_model=_dict_lookup(case_dir, "constant/electroProperties", "ionicModel"),
        )
