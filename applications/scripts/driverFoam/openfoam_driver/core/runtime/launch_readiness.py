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
#     launch_readiness
#
# Description
#     P2.1: the shared predicates every CLI launch/outcome decision must use.
#
#     Before this module, cli.py computed launch/exit decisions from three
#     separate, independently hand-written checks, each over a different
#     diagnostic set:
#       - ``report.status != "ok"``            (structural/semantic plan
#         validity, from ``StrictPlanReport.status`` in strict_planning.py)
#       - ``_environment_errors(diagnostics)``  (execution-readiness, from
#         ``environment_diagnostics``)
#       - ``status == "ok"`` via
#         ``_terminal_status_label(workflow_status)`` (post-execution
#         outcome, from ``WorkflowRunState``/``WorkflowStepState.status``)
#
#     ``is_launchable`` unifies the first two -- both are pre-execution
#     "can this plan be dispatched" questions -- into one predicate, so the
#     classification of a diagnostic as blocking ("error") vs. informational
#     ("warning") lives in exactly one place. The third is a genuinely
#     different question (did an already-dispatched workflow finish
#     successfully), so it gets its own small parallel predicate,
#     ``is_execution_successful``, rather than being forced into
#     ``is_launchable``.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from dataclasses import dataclass

from openfoam_driver.core.planning_types import StrictDiagnostic


@dataclass(frozen=True)
class LaunchReadiness:
    """One consolidated verdict over a plan's pre-execution readiness.

    ``structural_ok`` and ``environment_ok`` are exposed individually, not
    just folded into ``launchable``, because some CLI call sites only ever
    need one half: ``action=plan`` reports readiness without requiring the
    execution environment to be available (planning must work on a machine
    without OpenFOAM installed), while the environment gate only applies
    once a run is actually about to be dispatched -- and by that point
    structural validity has already been established separately. Route both
    kinds of call sites through this same predicate; just read the field
    each site actually needs.
    """

    launchable: bool
    structural_ok: bool
    environment_ok: bool
    has_warnings: bool
    blocking_reason: str | None


def is_launchable(
    *,
    plan_status: str,
    environment_diagnostics: tuple[StrictDiagnostic, ...] = (),
) -> LaunchReadiness:
    """Compute launch readiness from a plan's status and environment diagnostics.

    ``plan_status`` is ``StrictPlanReport.status`` ("ok"/"failed"), already
    computed upstream from structural/semantic ``plan_diagnostics`` only --
    warn-only diagnostic sets (function-object field warnings, environment
    diagnostics) never factor into it by design, so a sampled-field or
    environment warning can never fail a plan.

    ``environment_diagnostics`` is the plan's execution-readiness diagnostics
    (e.g. missing OpenFOAM/MPI/solver binary). Only ``level == "error"``
    entries block launch; ``level == "warning"`` entries are surfaced via
    ``has_warnings`` but never block.
    """
    structural_ok = plan_status == "ok"
    environment_errors = tuple(
        diagnostic for diagnostic in environment_diagnostics if diagnostic.level == "error"
    )
    environment_warnings = tuple(
        diagnostic for diagnostic in environment_diagnostics if diagnostic.level == "warning"
    )
    environment_ok = not environment_errors
    has_warnings = bool(environment_warnings)
    launchable = structural_ok and environment_ok

    blocking_reason: str | None = None
    if not structural_ok:
        blocking_reason = "plan is not structurally/semantically valid"
    elif not environment_ok:
        blocking_reason = "execution environment is not ready"

    return LaunchReadiness(
        launchable=launchable,
        structural_ok=structural_ok,
        environment_ok=environment_ok,
        has_warnings=has_warnings,
        blocking_reason=blocking_reason,
    )


def is_execution_successful(workflow_status: str) -> bool:
    """Whether a workflow status is the strict success terminal state.

    Mirrors the strict success contract shared by step and run: only
    ``"completed"`` counts. Any other status (``pending``, ``running``,
    ``failed``, ``skipped``) is not a success at this boundary -- decisions
    derive from status, never from a raw subprocess exit code.
    """
    return workflow_status == "completed"
