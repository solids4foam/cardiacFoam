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
#     Post-DAG hand-off, run once a workflow or sweep reaches a terminal
#     state. Splits into a "brain" that grounds itself in what the sweep
#     actually did, and an independent postprocessing module that consumes
#     that grounding instead of re-deriving it.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""The phase that runs after the deterministic sweep, not part of it.

The DAG (workflow_orchestrator.run_workflow) and the sweep loop
(sweep_runner.sweep_run) are both deterministic: they decide which cases run
and execute them, nothing more. Nothing in that engine, or in the CLI
dispatch that drives it, used to hand off to anything once the DAG or sweep
finished -- the previous postprocess mechanism (PostprocessTask,
run_postprocess_tasks, TutorialSpec.postprocess) was never actually called
by the execution engine and was removed 2026-08-18.

This module is the replacement hand-off point, split into two independent
pieces on purpose:

* **The brain** (`build_sweep_context`) reads the sweep's own record
  (sweep_manifest.json: when it started, when it finished, case_id ->
  resolved_axis_values, status) and then verifies that record against what
  is actually on disk for each case -- resolving each case's real output
  directory (which may or may not sit under the sweep's own --output-dir;
  entry-mode cases live under the target tutorial's own case_root instead,
  see sweep_runner._relative_or_absolute) and listing the files genuinely
  found there. Its output, `SweepContext`, is the single grounded picture of
  "what ran and where its output went."

* **The postprocessing module** (`run_postprocessing_module`) is a separate,
  independent function that receives that `SweepContext` as input. It never
  re-reads the manifest or re-derives file locations itself -- that is
  exactly the confusion the brain/module split avoids. It currently does no
  real analysis (`run_postprocessing_module` returns a stub outcome), but its
  signature already reflects the real contract: analysis code lands here,
  consuming grounded context, not raw paths it has to re-verify.

The module isn't limited to the flat summary in `SweepContext`, though.
`read_case_workflow_state` and `read_case_output_file` let it (or a
reasoning agent driving it) ask the brain for more -- full workflow_state.json
detail, or a specific output file's content -- without ever bypassing what
the brain already verified: both raise clearly on an unknown case_id, and
`read_case_output_file` only reads a path already present in the brain's own
`CaseRecord.output_files` scan, never an arbitrary path the caller guesses.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .sweep_manifest import read_manifest


@dataclass(frozen=True)
class PostprocessOutcome:
    status: str
    message: str

    def to_json(self) -> dict[str, Any]:
        return {"status": self.status, "message": self.message}


def run_postprocess_phase(*, entry: str | None, output_dir: Path) -> PostprocessOutcome:
    """Placeholder post-DAG hand-off for a single (non-sweep) run.

    A single `run --strict` invocation has no sweep_manifest.json to ground
    itself in -- it's one case, and _execute_run already has that case's
    real workflow_state in hand before calling this. Proves the wiring;
    does no real work yet.
    """
    return PostprocessOutcome(
        status="stub",
        message=f"postprocess stub: entry={entry} output_dir={output_dir}",
    )


@dataclass(frozen=True)
class CaseRecord:
    case_id: str
    resolved_axis_values: dict[str, Any]
    status: str
    outcome: str
    workflow_state_path: str
    case_output_dir: str | None
    output_files: tuple[str, ...]

    def to_json(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "resolved_axis_values": dict(self.resolved_axis_values),
            "status": self.status,
            "outcome": self.outcome,
            "workflow_state_path": self.workflow_state_path,
            "case_output_dir": self.case_output_dir,
            "output_files": list(self.output_files),
        }


@dataclass(frozen=True)
class SweepContext:
    output_dir: str
    sweep_spec_hash: str
    started_at: str
    finished_at: str
    case_count: int
    completed_count: int
    failed_count: int
    cases: tuple[CaseRecord, ...]

    def to_json(self) -> dict[str, Any]:
        return {
            "output_dir": self.output_dir,
            "sweep_spec_hash": self.sweep_spec_hash,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "case_count": self.case_count,
            "completed_count": self.completed_count,
            "failed_count": self.failed_count,
            "cases": [case.to_json() for case in self.cases],
        }


def _resolve_case_output(workflow_state_raw: str, *, output_dir: Path) -> tuple[Path, Path | None, tuple[str, ...]]:
    """Resolve one case's real output directory and list what's on it.

    `workflow_state_path` in sweep_manifest.json is relative to output_dir
    for generic/case-folder sweeps, but already absolute for entry-mode
    sweeps (whose case lives under the target tutorial's own case_root --
    see sweep_runner._relative_or_absolute). Handle both; never assume a
    fixed subpath like "postProcessing/" exists.
    """
    candidate = Path(workflow_state_raw)
    workflow_state_path = candidate if candidate.is_absolute() else (output_dir / candidate)

    if not workflow_state_path.exists():
        return workflow_state_path, None, ()

    case_output_dir = workflow_state_path.parent
    output_files = tuple(sorted(
        str(path.relative_to(case_output_dir))
        for path in case_output_dir.rglob("*")
        if path.is_file()
    ))
    return workflow_state_path, case_output_dir, output_files


def build_sweep_context(output_dir: Path) -> SweepContext:
    """The brain: read sweep_manifest.json, then verify it against disk.

    For every case the manifest records, resolve its real output directory
    (entry-mode and generic-mode sweeps place it differently -- see
    `_resolve_case_output`) and list the files genuinely found there. This
    is the single grounded picture `run_postprocessing_module` consumes; it
    never re-derives any of this itself.
    """
    output_dir = Path(output_dir)
    manifest = read_manifest(output_dir / "sweep_manifest.json")

    cases: list[CaseRecord] = []
    for case_entry in manifest.cases:
        workflow_state_path, case_output_dir, output_files = _resolve_case_output(
            case_entry.workflow_state_path, output_dir=output_dir,
        )
        cases.append(CaseRecord(
            case_id=case_entry.case_id,
            resolved_axis_values=dict(case_entry.resolved_axis_values),
            status=case_entry.status,
            outcome=case_entry.outcome,
            workflow_state_path=str(workflow_state_path),
            case_output_dir=str(case_output_dir) if case_output_dir is not None else None,
            output_files=output_files,
        ))

    return SweepContext(
        output_dir=str(output_dir),
        sweep_spec_hash=manifest.sweep_spec_hash,
        started_at=manifest.created_at,
        finished_at=manifest.updated_at,
        case_count=len(manifest.cases),
        completed_count=sum(1 for case in manifest.cases if case.status == "completed"),
        failed_count=sum(1 for case in manifest.cases if case.status == "failed"),
        cases=tuple(cases),
    )


def run_postprocessing_module(context: SweepContext) -> PostprocessOutcome:
    """Placeholder postprocessing module. Still no real analysis -- but now
    driven entirely by the brain's grounded SweepContext, not by re-reading
    the manifest or the filesystem itself."""
    case_summaries = ", ".join(
        f"{case.case_id}: {len(case.output_files)} file(s) @ {case.case_output_dir}"
        for case in context.cases
    )
    return PostprocessOutcome(
        status="stub",
        message=(
            f"postprocess stub: sweep {context.sweep_spec_hash} "
            f"({context.completed_count}/{context.case_count} completed) -- {case_summaries}"
        ),
    )


def _case_by_id(context: SweepContext, case_id: str) -> CaseRecord:
    for case in context.cases:
        if case.case_id == case_id:
            return case
    raise KeyError(f"{case_id!r} is not a case in this SweepContext (sweep {context.sweep_spec_hash})")


def read_case_workflow_state(context: SweepContext, case_id: str) -> dict[str, Any]:
    """On-demand deeper read: the full workflow_state.json for one case.

    CaseRecord only summarizes status/outcome; this returns the real
    per-step detail (status, attempt, produced_artifacts, diagnostics, ...)
    for callers that need more than the summary.
    """
    case = _case_by_id(context, case_id)
    workflow_state_path = Path(case.workflow_state_path)
    if not workflow_state_path.is_file():
        raise FileNotFoundError(
            f"workflow_state_path for case {case_id!r} no longer exists: {workflow_state_path}"
        )
    return json.loads(workflow_state_path.read_text())


def read_case_output_file(context: SweepContext, case_id: str, relative_path: str) -> str:
    """On-demand deeper read: one case's output file content.

    `relative_path` must be one the brain already found during
    `build_sweep_context` (present in `case.output_files`) -- this can never
    be used to read a file the brain hasn't verified exists.
    """
    case = _case_by_id(context, case_id)
    if relative_path not in case.output_files:
        raise KeyError(
            f"{relative_path!r} is not among the verified output files for case "
            f"{case_id!r}: {list(case.output_files)}"
        )
    if case.case_output_dir is None:
        raise FileNotFoundError(f"case {case_id!r} has no case_output_dir on record")
    return (Path(case.case_output_dir) / relative_path).read_text()
