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
#     execution_context
#
# Description
#     Resolves the case_root/setup_root/output_dir/workflow_state_path a strict
#     plan needs, independent of any legacy CLI action vocabulary.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Neutral path resolver for strict_plan().

strict_plan() used to call describe_launch("sim", entry, ...) purely to get
these four paths -- that re-resolved the entry (resolve_entry + factory) a
second time (strict_plan already has the spec from load_entry_spec) and,
worse, coupled the strict/workflow-DAG execution path (which never runs the
legacy sim/post/all CLI) to describe_launch's VALID_DRIVER_ACTIONS vocabulary.
This module takes the already-built TutorialSpec directly: no re-resolution,
no action string, no dependency on which CLI actions happen to exist.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .models import TutorialSpec


@dataclass(frozen=True)
class ExecutionContext:
    case_root: Path
    setup_root: Path
    output_dir: Path
    workflow_state_path: Path


def resolve_execution_context(spec: TutorialSpec) -> ExecutionContext:
    output_dir = Path(spec.output_dir)
    return ExecutionContext(
        case_root=Path(spec.case_root),
        setup_root=Path(spec.setup_root),
        output_dir=output_dir,
        workflow_state_path=output_dir / "workflow_state.json",
    )
