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
#     case_provenance
#
# Description
#     cardiacFoam's CaseProvenanceCapability: which case files are inputs the
#     workflow consumes versus generated diagnostics nothing reads.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Solver-declared case classification for the provenance snapshot.

``required_inputs`` intentionally returns ``()`` for now. The full
model-dependent resolution -- which ``0/`` fields are present depends on the
configured solver and ionic model, and their *locations* resolve by a
backward ``Time::findInstance`` search with a ``constant/`` fallback, so a
field's canonical path is not knowable from its name alone -- is Task 2b's
input-enumeration job, not this task's. Returning ``()`` here is still safe:
the resolution precedence an unclassified file falls back to is
``required_input``, so nothing here can under-classify a file that turns out
to matter.

``generated_output_globs`` is the one classification this task does own:
``constant/C``, ``Cx``, ``Cy``, ``Cz`` and ``skewness`` are mesh-diagnostic
byproducts. An exhaustive grep across ``src/`` and ``applications/`` found
zero reads of any of them -- nothing in the solver or its utilities ever
opens these files, so excluding them from the input snapshot (rather than
fingerprinting large diagnostic arrays nothing depends on) does not risk a
silent stale replay.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from openfoam_driver.core.plugin_capabilities import ResolvedInput

# Fixed names, independent of the resolved model -- unlike the 0/ solver
# fields, these mesh-diagnostic byproducts are never renamed by a dictionary
# key.
_GENERATED_OUTPUT_GLOBS: tuple[str, ...] = (
    "constant/C",
    "constant/Cx",
    "constant/Cy",
    "constant/Cz",
    "constant/skewness",
)


def required_inputs(
    case_root: Path,
    resolved_case: dict[str, Any],
    selected_start_time: str,
) -> tuple[ResolvedInput, ...]:
    """Deferred to Task 2b. See the module docstring for why ``()`` is safe."""
    del case_root, resolved_case, selected_start_time
    return ()


def generated_output_globs(
    case_root: Path,
    resolved_case: dict[str, Any],
    selected_start_time: str,
) -> tuple[str, ...]:
    """``constant/C``, ``Cx``, ``Cy``, ``Cz``, ``skewness`` -- read by nothing."""
    del case_root, resolved_case, selected_start_time
    return _GENERATED_OUTPUT_GLOBS
