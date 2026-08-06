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
#     fresh
#
# Description
#     Safety-gated deletion for the --fresh flag: deletes a resolved output
#     directory before a run/sweep-run so a stale workflow_state.json /
#     sweep_manifest.json cannot silently resume as "completed". See
#     docs/superpowers/specs/2026-08-06-driverfoam-fresh-flag-design.md for
#     the incident and design rationale.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import shutil
import sys
from pathlib import Path

_DRIVERFOAM_MARKER_NAMES = ("workflow_state.json", "sweep_manifest.json", "run_document.json")


def _has_driverfoam_marker(output_dir: Path) -> bool:
    """True if output_dir contains a recognizable driverFOAM artifact.

    Bounded to the top level and one level of subdirectories -- markers
    always live at a case root or a sweep's per-case root, and an unbounded
    walk would be slow across large mesh trees.
    """
    for name in _DRIVERFOAM_MARKER_NAMES:
        if (output_dir / name).exists():
            return True
    for child in output_dir.iterdir():
        if not child.is_dir():
            continue
        for name in _DRIVERFOAM_MARKER_NAMES:
            if (child / name).exists():
                return True
    return False


def check_fresh_deletion_allowed(output_dir: Path, *, allowed_root: Path | None) -> str | None:
    """Return an error message if deleting output_dir would be unsafe, else None.

    output_dir need not exist -- a nonexistent directory always passes (there
    is nothing to lose).
    """
    resolved = output_dir.resolve()
    if resolved.parent == resolved:
        return f"--fresh refuses to delete the filesystem root ({resolved})."
    if resolved == Path.home().resolve():
        return f"--fresh refuses to delete the home directory ({resolved})."
    depth = len(resolved.parts) - 1
    if depth < 3:
        return (
            f"--fresh refuses to delete {resolved}: fewer than 3 path segments "
            "beneath the filesystem root (too shallow to be a case/sweep output "
            "directory)."
        )
    if allowed_root is not None and not resolved.is_relative_to(allowed_root):
        return (
            f"--fresh refuses to delete {resolved}: outside "
            f"DRIVERFOAM_ALLOWED_RUNS_ROOT ({allowed_root})."
        )
    if resolved.exists() and not _has_driverfoam_marker(resolved):
        return (
            f"--fresh refuses to delete {resolved}: directory exists but contains "
            "no recognizable driverFOAM artifact (workflow_state.json, "
            "sweep_manifest.json, or run_document.json) at its top level or one "
            "level of subdirectories. Check --output-dir for a typo."
        )
    return None


def ensure_fresh_output_dir(
    output_dir: Path, *, fresh: bool, allowed_root: Path | None = None
) -> str | None:
    """Delete output_dir when fresh=True and the safety guards allow it.

    Returns an error message (leaving output_dir untouched) if a guard
    refuses; returns None on success, including the fresh=False no-op and the
    "didn't exist" no-op. Prints the resolved path to stderr before deleting
    it -- the only audit trail, since there is no interactive prompt. Stderr
    (not stdout) so the CLI's stdout stays pure JSON.
    """
    if not fresh:
        return None
    error = check_fresh_deletion_allowed(output_dir, allowed_root=allowed_root)
    if error is not None:
        return error
    resolved = output_dir.resolve()
    if resolved.exists():
        print(f"--fresh: deleting {resolved}", file=sys.stderr)
        shutil.rmtree(resolved)
    return None
