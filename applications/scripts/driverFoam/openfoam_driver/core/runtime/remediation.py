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
#     remediation
#
# Description
#     Suggestion-only remediation ladder for a failed workflow step. Computes a
#     best-effort list of candidate fixes (RemediationHint) from the failure
#     bundle: static diagnostic-code hints first, then signature-gated log
#     interpretation only when no structured diagnostic code is present. Pure:
#     never does I/O, never raises (returns () on any internal problem), never
#     mutates state. The driver only *suggests*; the agent decides and applies.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class RemediationHint:
    diagnostic_code: str
    driver_path: str   # "" => advisory (no mutation)
    change: str        # human-readable transform; descriptive, never executed by the driver
    rationale: str
    source: str        # "static" | "log_signature"
    confidence: str    # "high" | "low"

    def to_json(self) -> dict[str, Any]:
        return {
            "diagnostic_code": self.diagnostic_code,
            "driver_path": self.driver_path,
            "change": self.change,
            "rationale": self.rationale,
            "source": self.source,
            "confidence": self.confidence,
        }


# Static, diagnostic-code-keyed hints. Advisory only (empty driver_path) — these
# failures are not fixed by a catalog mutation.
STATIC_REMEDIATION_HINTS: dict[str, tuple[RemediationHint, ...]] = {
    "missing_artifacts": (
        RemediationHint(
            diagnostic_code="missing_artifacts",
            driver_path="",
            change="",
            rationale=(
                "The command reported success but produced no expected artifact. "
                "Verify the producing utility/pre-solve step actually ran and wrote output."
            ),
            source="static",
            confidence="low",
        ),
    ),
    "workflow_step_exec_error": (
        RemediationHint(
            diagnostic_code="workflow_step_exec_error",
            driver_path="",
            change="",
            rationale=(
                "The step command could not be launched. Check the executable is on PATH "
                "and the OpenFOAM environment is sourced."
            ),
            source="static",
            confidence="low",
        ),
    ),
}


def _static_hints(failure_context: dict[str, Any]) -> tuple[RemediationHint, ...]:
    out: list[RemediationHint] = []
    for diagnostic in failure_context.get("diagnostics", ()):
        code = diagnostic.get("code")
        out.extend(STATIC_REMEDIATION_HINTS.get(code, ()))
    return tuple(out)


# Signatures that point at solver temporal instability. Matched case-insensitively
# against the combined log tails. Multi-word markers use substring match; short markers
# are matched as whole tokens to avoid false positives (e.g. "nan" inside "meaning").
_DIVERGENCE_SUBSTRINGS: tuple[str, ...] = (
    "maximum number of iterations",
    "singularity",
    "floating point exception",
)
_DIVERGENCE_TOKENS: frozenset[str] = frozenset({"nan", "inf"})

_HALVE_DELTAT_HINT = RemediationHint(
    diagnostic_code="",
    driver_path="deltaT",
    change="halve",
    rationale=(
        "Solver divergence / FOAM FATAL with no structured diagnostic code. A smaller "
        "time step is the conservative first remedy; halve controlDict deltaT and rerun."
    ),
    source="log_signature",
    confidence="low",
)


def _has_divergence_signature(failure_context: dict[str, Any]) -> bool:
    blob = (
        f"{failure_context.get('stdout_tail', '') or ''}\n"
        f"{failure_context.get('stderr_tail', '') or ''}"
    ).lower()
    if any(sig in blob for sig in _DIVERGENCE_SUBSTRINGS):
        return True
    return bool(_DIVERGENCE_TOKENS & set(re.findall(r"\w+", blob)))


def interpret_log_signatures(failure_context: dict[str, Any]) -> tuple[RemediationHint, ...]:
    """Last-resort: infer a candidate from the bounded log tails.

    Returns the conservative deltaT hint only when a known divergence signature is present
    in the tails, and () otherwise. By the time this layer runs the deterministic layers
    had no explanation; interpretation is licensed to *recognize* a known failure, never to
    guess at an unrecognized one.
    """
    try:
        if _has_divergence_signature(failure_context):
            return (_HALVE_DELTAT_HINT,)
        return ()
    except Exception:
        return ()


def build_candidate_remediations(failure_context: dict[str, Any]) -> tuple[RemediationHint, ...]:
    """Suggestion-only remediation ladder. Never raises; returns () on any problem."""
    try:
        static = _static_hints(failure_context)
        if static:
            return static
        # Interpretation is licensed only for failures the deterministic layer could not
        # name at all: a nonzero exit with an empty diagnostics tuple. A coded failure with
        # no static hint (e.g. workflow_step_timeout) yields no candidate rather than
        # falling through to a (here, counter-productive) deltaT suggestion.
        if failure_context.get("diagnostics"):
            return ()
        return interpret_log_signatures(failure_context)
    except Exception:
        return ()
