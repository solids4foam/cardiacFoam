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
#     postprocess_phase
#
# Description
#     Post-DAG hand-off, run once a workflow reaches a terminal state.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""The phase that runs after the deterministic DAG, not part of it.

The DAG (workflow_orchestrator.run_workflow) is deterministic: mesh, solve,
whatever shell steps a tutorial declares. Nothing in that engine, or in the
CLI dispatch that drives it, used to hand off to anything once the DAG
finished -- the previous postprocess mechanism (PostprocessTask,
run_postprocess_tasks, TutorialSpec.postprocess) was never actually called
by the execution engine and was removed 2026-08-18.

This module is the replacement hand-off point. It currently does no real
work -- `run_postprocess_phase` returns a stub outcome -- but it establishes
where the hand-off lives and what shape it returns, so real analysis can
land here later without re-plumbing the call site.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class PostprocessOutcome:
    status: str
    message: str

    def to_json(self) -> dict[str, Any]:
        return {"status": self.status, "message": self.message}


def run_postprocess_phase(*, entry: str, output_dir: Path) -> PostprocessOutcome:
    """Placeholder post-DAG hand-off. Proves the wiring; does no real work yet."""
    return PostprocessOutcome(
        status="stub",
        message=f"postprocess stub: entry={entry} output_dir={output_dir}",
    )
