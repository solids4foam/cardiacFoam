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
#     resume_guard
#
# Description
#     Interim, warning-only guard against replaying a completed workflow state
#     that predates a solver rebuild. Superseded by checkpointed provenance in
#     Phase 2; this module is intended to be deleted then.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path
from typing import Any

_CODE = "possibly_stale_resume"


def _finished_epoch(finished_at: str | None) -> float | None:
    if not finished_at:
        return None
    try:
        return datetime.fromisoformat(finished_at).timestamp()
    except ValueError:
        return None


def _executable_mtime(command: str, case_root: Path) -> float | None:
    """Resolve ``command`` the way the executor would and return its mtime.

    Only PATH resolution is attempted. A case-local script is deliberately not
    resolved: its own mtime says nothing about the solver binary it invokes,
    which is the thing this guard is about.
    """
    if not command or "/" in command:
        return None
    resolved = shutil.which(command)
    if resolved is None:
        return None
    try:
        return Path(resolved).stat().st_mtime
    except OSError:
        return None


def stale_resume_warnings(
    state: Any, *, case_root: Path
) -> tuple[dict[str, str], ...]:
    """Warn for each completed step whose executable is newer than the step.

    Warning-only and best-effort by construction: an unresolvable command, an
    unparseable timestamp, or a step that never completed all yield nothing.
    This cannot fail a run and must never be treated as provenance.
    """
    warnings: list[dict[str, str]] = []
    for step in getattr(state, "steps", ()):
        if getattr(step, "status", "") != "completed":
            continue
        finished = _finished_epoch(getattr(step, "finished_at", None))
        if finished is None:
            continue
        mtime = _executable_mtime(getattr(step, "command", ""), case_root)
        if mtime is None or mtime <= finished:
            continue
        warnings.append({
            "code": _CODE,
            "message": (
                f"Step {step.step_id!r} is recorded as completed at "
                f"{step.finished_at}, but its executable "
                f"{step.command!r} was modified more recently. Resuming will "
                "report the previous run's results without re-executing. "
                "Pass --fresh if this run is meant to reflect the rebuild."
            ),
            "field": str(step.step_id),
        })
    return tuple(warnings)
