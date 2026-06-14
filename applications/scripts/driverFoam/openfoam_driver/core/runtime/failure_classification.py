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
#     failure_classification
#
# Description
#     Classifies a failed workflow step as retryable (transient) or fatal
#     (deterministic), from the diagnostic codes the runner emits. Pure: no
#     log scraping, no environment access.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

# workflow_step_timeout is the canonical transient failure (load/contention).
# Everything else is deterministic by default: missing_artifacts (command
# "succeeded" but produced nothing), workflow_step_exec_error (could not
# launch), and generic nonzero exit (FOAM FATAL ERROR, divergence).
RETRYABLE_CODES = frozenset({"workflow_step_timeout"})


def classify_failure(step_state, *, overrides: dict[str, str] | None = None) -> str:
    """Return "retryable" or "fatal" for a failed step.

    overrides maps a diagnostic code to a forced classification — the extension
    seam for problem #3 / agent reclassification. The mechanical retry path
    never populates it.
    """
    overrides = overrides or {}
    for diagnostic in step_state.diagnostics:
        code = diagnostic.get("code")
        if code in overrides:
            return overrides[code]
        if code in RETRYABLE_CODES:
            return "retryable"
    return "fatal"
