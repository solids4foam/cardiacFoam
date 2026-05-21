"""Post-run artifact reconciliation (plan §10).

`reconcile_artifacts(case_root, predicted)` walks the predicted
`DataArtifact` set and checks which expected files exist on disk under
`case_root`. The resulting `ReconciliationReport` tells the agent which
predictions were realised and which were not — useful for confirming
that a run produced what the predictor advertised, and for catching
predictor drift over time.

Design notes:

* Non-time-indexed artifacts are checked at a single literal path
  derived by stripping the closed-set placeholders. A `{case_id}`
  placeholder is left unresolved at this layer — callers wanting per-case
  expansion should call the reconciler per case_id (the engine wires this
  up once per sweep run).
* Time-indexed artifacts are matched by globbing top-level numeric-looking
  directory names under `case_root` (the OpenFOAM convention) and then
  checking the per-time relative path. This avoids walking
  `constant/`, `system/`, `processor*/`, and any other non-time scope.
* "Extra files" (on-disk content not advertised by any prediction) is out
  of scope for v1 — the case tree is large and the value of enumerating
  every unmodelled file is low. If a real consumer needs it, add it then.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .models import DataArtifact


# OpenFOAM time directories use a decimal-number naming convention.
# `constant`, `system`, `processor*`, and any other named directory is
# excluded by this filter.
_TIME_DIR_RE: re.Pattern[str] = re.compile(r"^-?\d+(\.\d+)?(e[+\-]?\d+)?$")


@dataclass(frozen=True)
class ReconciliationReport:
    """Result of comparing predicted artifacts against on-disk reality."""

    case_root: str
    predicted_count: int
    matched_count: int
    missing_count: int
    artifacts: tuple[dict, ...]
    """One entry per predicted artifact, in input order. Each entry is a
    dict with keys: artifact_id, predicted_path, status (matched|missing),
    matched_files (list of {path, size_bytes}), optional (bool)."""


def _list_time_dirs(case_root: Path) -> list[str]:
    """Return every top-level directory name under case_root that looks
    like an OpenFOAM time directory (numeric)."""
    return sorted(
        child.name for child in case_root.iterdir()
        if child.is_dir() and _TIME_DIR_RE.match(child.name)
    )


def _check_path(case_root: Path, relative: str) -> dict | None:
    """Return a `{path, size_bytes}` dict if the file exists, else None."""
    full = case_root / relative
    if not full.exists() or not full.is_file():
        return None
    return {"path": str(full), "size_bytes": full.stat().st_size}


def _reconcile_artifact(
    case_root: Path,
    artifact: DataArtifact,
) -> dict:
    """Classify one artifact, returning the per-entry report dict."""
    matched_files: list[dict] = []

    if artifact.time_indexed:
        # Glob every time directory and check the post-{time} portion.
        # The pattern is expected to begin with `{time}` for OpenFOAM
        # outputs; everything after that is the per-time relative path.
        for time_name in _list_time_dirs(case_root):
            resolved = artifact.path_pattern.replace("{time}", time_name)
            # Strip any unresolved {case_id} — the per-case fan-out is
            # the engine's responsibility, not the reconciler's.
            resolved = resolved.replace("{case_id}", "")
            hit = _check_path(case_root, resolved)
            if hit is not None:
                matched_files.append(hit)
    else:
        resolved = artifact.path_pattern.replace("{case_id}", "")
        hit = _check_path(case_root, resolved)
        if hit is not None:
            matched_files.append(hit)

    status = "matched" if matched_files else "missing"
    return {
        "artifact_id": artifact.artifact_id,
        "predicted_path": artifact.path_pattern,
        "status": status,
        "matched_files": matched_files,
        "optional": artifact.optional,
    }


def reconcile_artifacts(
    case_root: Path,
    predicted: Iterable[DataArtifact],
) -> ReconciliationReport:
    """Compare `predicted` artifacts against the on-disk state of
    `case_root`. Returns a `ReconciliationReport` enumerating every
    predicted artifact and whether (and how) it was realised."""
    predicted_tuple = tuple(predicted)
    entries: list[dict] = []
    matched_count = 0
    for artifact in predicted_tuple:
        entry = _reconcile_artifact(case_root, artifact)
        if entry["status"] == "matched":
            matched_count += 1
        entries.append(entry)
    return ReconciliationReport(
        case_root=str(case_root),
        predicted_count=len(predicted_tuple),
        matched_count=matched_count,
        missing_count=len(predicted_tuple) - matched_count,
        artifacts=tuple(entries),
    )
