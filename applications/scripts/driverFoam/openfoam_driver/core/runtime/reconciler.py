#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     reconciler
#
# Description
#     Performs post-run reconciliation of expected versus actual generated artifacts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Post-run artifact reconciliation.

`reconcile_artifacts(case_root, predicted, *, case_id=None)` walks the
predicted `DataArtifact` set and reports which expected paths exist on
disk under `case_root`. The resulting `ReconciliationReport` tells the
agent which predictions were realised, which were not, and for matched
entries whether the realisation is a file or a directory.

Design notes (rewritten after the 2026-05-21 audit):

* **Matching uses `pathlib.Path.glob`**, not literal `exists()`. The
  predictor's `path_pattern` may contain placeholders (`{case_id}`,
  `{time}`) AND shell-glob wildcards once substituted. Globbing handles
  both literal paths (degenerate glob) and wildcards uniformly.
* **Directories are valid match targets.** The monodomain/bidomain/eikonal
  predictors emit `path_pattern="{time}"` — the OpenFOAM time directory
  itself, not a file inside it. Earlier code rejected directories via
  `is_file()`; the new code accepts both kinds and tags each match with
  `kind: "file" | "dir"`.
* **`{case_id}` substitution is explicit.** Callers pass `case_id=` when
  they want a literal substitution. When omitted, `{case_id}` is replaced
  with the glob wildcard `*` so the reconciler still finds matches across
  all per-case outputs in a sweep. Earlier code silently stripped
  `{case_id}` to empty string, producing malformed paths like
  `postProcessing/.txt` that never matched anything.
* **Time-indexed iteration** walks top-level directory names matching
  the OpenFOAM time-dir convention (decimal numerals) so `constant/`,
  `system/`, `processor*/` are filtered out by name.
* **"Extra files"** (on-disk content not advertised by any prediction)
  remain out of scope. The case tree is large and the value of
  enumerating every unmodelled file is low.
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .models import DataArtifact


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
    matched_files (list of {path, kind, size_bytes?, entries?}),
    optional (bool)."""

    case_id: str | None = None
    """Case identifier when invoked per-case during a sweep. ``None`` for
    whole-run reconciliations."""

    def to_json(self) -> dict:
        """Plain-dict form for embedding in run payloads and manifests."""
        return {
            "case_root": self.case_root,
            "case_id": self.case_id,
            "predicted_count": self.predicted_count,
            "matched_count": self.matched_count,
            "missing_count": self.missing_count,
            "artifacts": [dict(entry) for entry in self.artifacts],
        }


def _list_time_dirs(case_root: Path) -> list[str]:
    """Return every top-level directory name under case_root that looks
    like an OpenFOAM time directory (numeric)."""
    return sorted(
        child.name for child in case_root.iterdir()
        if child.is_dir() and _TIME_DIR_RE.match(child.name)
    )


def _sha256_of(path: Path) -> str | None:
    """Stream a file's sha256. Returns None if it cannot be read."""
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(65536), b""):
                digest.update(chunk)
    except OSError:
        return None
    return digest.hexdigest()


def _entries_for_match(path: Path) -> dict | None:
    """Build the per-match dict for a single glob hit. Returns None for
    matches that are neither files nor directories (symlinks to nowhere,
    etc.) so the caller can skip them silently."""
    if path.is_file():
        entry = {
            "path": str(path),
            "kind": "file",
            "size_bytes": path.stat().st_size,
        }
        digest = _sha256_of(path)
        if digest is not None:
            entry["sha256"] = digest
        return entry
    if path.is_dir():
        # Count immediate children so an agent inspecting the realized
        # manifest knows whether the matched dir is non-empty.
        try:
            entry_count = sum(1 for _ in path.iterdir())
        except OSError:
            entry_count = 0
        return {
            "path": str(path),
            "kind": "dir",
            "entries": entry_count,
        }
    return None


def _glob_under(case_root: Path, pattern: str) -> list[dict]:
    """Glob `case_root` for `pattern` (relative). Returns one match dict
    per hit. Skips entries that are neither file nor directory."""
    results: list[dict] = []
    try:
        hits = sorted(case_root.glob(pattern))
    except (NotImplementedError, ValueError):
        return results
    for hit in hits:
        entry = _entries_for_match(hit)
        if entry is not None:
            results.append(entry)
    return results


def _substitute_case_id(pattern: str, case_id: str | None) -> str:
    """Replace `{case_id}` with the literal case_id when supplied, else
    with the glob wildcard ``*`` so the matcher still finds per-case
    outputs across a sweep."""
    return pattern.replace("{case_id}", case_id if case_id is not None else "*")


def _reconcile_artifact(
    case_root: Path,
    artifact: DataArtifact,
    *,
    case_id: str | None,
) -> dict:
    """Classify one artifact, returning the per-entry report dict."""
    matched_files: list[dict] = []

    if artifact.time_indexed:
        for time_name in _list_time_dirs(case_root):
            resolved = artifact.path_pattern.replace("{time}", time_name)
            resolved = _substitute_case_id(resolved, case_id)
            matched_files.extend(_glob_under(case_root, resolved))
    else:
        resolved = _substitute_case_id(artifact.path_pattern, case_id)
        # Defensive: a non-time-indexed pattern shouldn't contain {time},
        # but if it does, accept any time dir name via wildcard.
        resolved = resolved.replace("{time}", "*")
        matched_files.extend(_glob_under(case_root, resolved))

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
    *,
    case_id: str | None = None,
) -> ReconciliationReport:
    """Compare `predicted` artifacts against the on-disk state of
    `case_root`.

    Args:
        case_root: the case directory to inspect.
        predicted: the predicted artifacts (typically from
            ``predict_data_artifacts``).
        case_id: when supplied, substituted into `{case_id}` placeholders
            literally. When omitted, `{case_id}` becomes the glob wildcard
            ``*`` so per-case outputs in a sweep still match.

    Returns:
        A `ReconciliationReport` enumerating every predicted artifact and
        whether (and how) it was realised.
    """
    predicted_tuple = tuple(predicted)
    entries: list[dict] = []
    matched_count = 0
    for artifact in predicted_tuple:
        entry = _reconcile_artifact(case_root, artifact, case_id=case_id)
        if entry["status"] == "matched":
            matched_count += 1
        entries.append(entry)
    return ReconciliationReport(
        case_root=str(case_root),
        predicted_count=len(predicted_tuple),
        matched_count=matched_count,
        missing_count=len(predicted_tuple) - matched_count,
        artifacts=tuple(entries),
        case_id=case_id,
    )
