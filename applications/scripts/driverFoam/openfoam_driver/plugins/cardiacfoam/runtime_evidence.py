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
#     runtime_evidence
#
# Description
#     Where cardiacFoam's runtime evidence lives: which steps actually solve,
#     where their logs land, what extra inputs belong in a provenance
#     snapshot, and how to read values out of solver-specific artifacts.
#     Declaration only -- consumed by later phases, not by Phase 1.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path

# Commands that run the solver. Deliberately excludes
# bathBidomainInterfaceMetrics, which is authorized to run but post-processes
# rather than solves -- expecting solver telemetry from it would be wrong.
_SOLVE_STEP_COMMANDS = frozenset({"cardiacFoam"})

# OpenFOAM's runApplication redirects solver output to log.<app>, so a step
# that runs an Allrun script produces no parseable driver-captured stdout.
# Phase 4's telemetry collector uses these globs to find the real log.
_TELEMETRY_GLOBS: dict[str, tuple[str, ...]] = {
    "Allrun": ("log.cardiacFoam", "log.*"),
    "Allrun.pre": ("log.*",),
    "Allrun.post": ("log.*",),
    "cardiacFoam": (),
}


def solve_step_commands() -> frozenset[str]:
    return _SOLVE_STEP_COMMANDS


def telemetry_source_globs(command: str) -> tuple[str, ...]:
    """Case-relative globs where ``command``'s solver log may land, beyond
    driver-captured stdout. Empty for commands that write no solver log."""
    return _TELEMETRY_GLOBS.get(command, ())


def extra_provenance_paths(case_root: Path) -> tuple[Path, ...]:
    """Additional inputs Phase 2 must digest beyond ``system/`` and
    ``constant/``.

    Empty today: cardiacFoam reads no plugin-owned data file outside the case.
    The hook exists so Phase 2 need not widen the plugin contract."""
    del case_root
    return ()


def artifact_value_reader(artifact_format: str):
    """Reader for a solver-specific artifact format, or ``None``.

    Empty today. Phase 5 registers readers here for cardiac formats such as
    ECG traces and Purkinje time series. Returning ``None`` must make Phase 5
    report ``not_evaluated`` with a reason -- never an implicit pass."""
    del artifact_format
    return None
