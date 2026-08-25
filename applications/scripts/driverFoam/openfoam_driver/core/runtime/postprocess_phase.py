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
  independent function that receives that `SweepContext`, plus a `task`
  (what the caller actually wants done), as input. It never re-reads the
  manifest or re-derives file locations itself -- that is exactly the
  confusion the brain/module split avoids. task(sweep) (sweep_run itself)
  is purely deterministic and has no task concept; the task only enters at
  this hand-off, where reasoning actually happens. It currently does no
  real analysis (`run_postprocessing_module` returns a stub outcome), but
  it does list each case's own setup/ postprocessing script catalog via
  `list_postprocess_scripts` -- real dispatch (run whichever cataloged
  script fits the task, or reason freely when none does) lands here next.

Every cataloged script carries a standard description (`PostprocessScriptInfo`
-- path, function_name, description), read statically from its docstring so
a reasoning agent can judge whether a script applies to the task at hand
without running or even fully reading it first.

The module isn't limited to the flat summary in `SweepContext`, though.
`read_case_workflow_state` and `read_case_output_file` let it (or a
reasoning agent driving it) ask the brain for more -- full workflow_state.json
detail, or a specific output file's content -- without ever bypassing what
the brain already verified: both raise clearly on an unknown case_id, and
`read_case_output_file` only reads a path already present in the brain's own
`CaseRecord.output_files` scan, never an arbitrary path the caller guesses.
"""
from __future__ import annotations

import ast
import json
import os
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


def build_standalone_case_record(
    *, entry: str, case_root: Path, setup_root: Path | None, output_dir: Path,
) -> "CaseRecord":
    """Build a standalone (non-sweep) run's case_record.json equivalent.

    Mirrors what build_sweep_context derives per sweep case, but a
    standalone run has no sweep_manifest.json to ground itself in --
    case_root and output_dir are already known directly from the CLI's own
    execution context (_execute_run has resolved them before calling this),
    so there is no manifest indirection to resolve.

    resolved_axis_values is always {} and outcome is always "fresh": a
    standalone run has no sweep axes and no manifest-tracked retry
    bookkeeping to report.
    """
    output_dir = Path(output_dir)
    workflow_state_path = output_dir / "workflow_state.json"
    status = "unknown"
    if workflow_state_path.is_file():
        try:
            status = json.loads(workflow_state_path.read_text()).get("status", "unknown")
        except json.JSONDecodeError:
            status = "unknown"
    output_files = _list_output_files(output_dir) if output_dir.is_dir() else ()
    return CaseRecord(
        case_id=entry,
        resolved_axis_values={},
        status=status,
        outcome="fresh",
        workflow_state_path=str(workflow_state_path),
        case_output_dir=str(output_dir),
        output_files=output_files,
        setup_root=str(setup_root) if setup_root else None,
        case_root=str(case_root),
        run_document_path=None,
        override_hash=None,
        started_at=None,
        updated_at=None,
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
    setup_root: str | None
    case_root: str | None = None
    run_document_path: str | None = None
    override_hash: str | None = None
    started_at: str | None = None
    updated_at: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "resolved_axis_values": dict(self.resolved_axis_values),
            "status": self.status,
            "outcome": self.outcome,
            "workflow_state_path": self.workflow_state_path,
            "case_output_dir": self.case_output_dir,
            "output_files": list(self.output_files),
            "setup_root": self.setup_root,
            "case_root": self.case_root,
            "run_document_path": self.run_document_path,
            "override_hash": self.override_hash,
            "started_at": self.started_at,
            "updated_at": self.updated_at,
        }


def write_case_record(path: Path, record: CaseRecord) -> None:
    """Persist one case's record atomically (os.replace of a .tmp sibling),
    matching sweep_manifest.write_manifest's crash-safety pattern -- a
    reader must never observe a half-written case_record.json.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(path.name + ".tmp")
    tmp_path.write_text(json.dumps(record.to_json(), indent=2))
    os.replace(tmp_path, path)


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


def _list_output_files(case_output_dir: Path) -> tuple[str, ...]:
    return tuple(sorted(
        str(path.relative_to(case_output_dir))
        for path in case_output_dir.rglob("*")
        if path.is_file()
    ))


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
    return workflow_state_path, case_output_dir, _list_output_files(case_output_dir)


def _read_run_document_field(run_document_path_raw: str, field: str, *, output_dir: Path) -> str | None:
    """Read a case's run_document.json and extract launch.<field>.

    run_document_path is always written under the sweep's own output_dir
    (sweep_runner._relative_or_absolute's docstring), unlike
    workflow_state_path -- no relative/absolute branching needed here.
    Returns None if the file is missing, malformed, or the field is absent
    -- never raises, since this is used for record-keeping (setup_root,
    case_root), not for judging whether the case itself succeeded.
    """
    run_document_path = output_dir / run_document_path_raw
    if not run_document_path.is_file():
        return None
    try:
        run_document = json.loads(run_document_path.read_text())
    except json.JSONDecodeError:
        return None
    value = run_document.get("launch", {}).get(field)
    return str(value) if value else None


def _resolve_case_setup_root(run_document_path_raw: str, *, output_dir: Path) -> str | None:
    """Read a case's run_document.json and extract launch.setupRoot."""
    return _read_run_document_field(run_document_path_raw, "setupRoot", output_dir=output_dir)


def _resolve_case_root(run_document_path_raw: str, *, output_dir: Path) -> str | None:
    """Read a case's run_document.json and extract launch.caseRoot -- the
    staged OpenFOAM case root the solver actually ran in (entry mode) or
    the case-folder root (generic mode). Same missing/malformed handling as
    _resolve_case_setup_root.
    """
    return _read_run_document_field(run_document_path_raw, "caseRoot", output_dir=output_dir)


def build_sweep_context(output_dir: Path) -> SweepContext:
    """The brain: read sweep_manifest.json, then verify it against disk.

    For every case the manifest records, resolve its real output directory
    (entry-mode and generic-mode sweeps place it differently -- see
    `_resolve_case_output`) and list the files genuinely found there. This
    is the single grounded picture `run_postprocessing_module` consumes; it
    never re-derives any of this itself.

    As a side effect, persists each case's full CaseRecord to
    `output_dir / case_entry.case_record_path` (skipped for a case whose
    manifest entry predates that field and so has an empty path) -- this is
    the durable case_record.json the design review asked for: one file per
    case that answers "what ran, with which inputs, where's its output"
    without an agent reconstructing paths from the manifest and disk itself.
    """
    output_dir = Path(output_dir)
    manifest = read_manifest(output_dir / "sweep_manifest.json")

    cases: list[CaseRecord] = []
    for case_entry in manifest.cases:
        workflow_state_path, case_output_dir, output_files = _resolve_case_output(
            case_entry.workflow_state_path, output_dir=output_dir,
        )
        setup_root = _resolve_case_setup_root(case_entry.run_document_path, output_dir=output_dir)
        case_root = _resolve_case_root(case_entry.run_document_path, output_dir=output_dir)
        record = CaseRecord(
            case_id=case_entry.case_id,
            resolved_axis_values=dict(case_entry.resolved_axis_values),
            status=case_entry.status,
            outcome=case_entry.outcome,
            workflow_state_path=str(workflow_state_path),
            case_output_dir=str(case_output_dir) if case_output_dir is not None else None,
            output_files=output_files,
            setup_root=setup_root,
            case_root=case_root,
            run_document_path=case_entry.run_document_path,
            override_hash=case_entry.override_hash,
            started_at=case_entry.started_at,
            updated_at=case_entry.updated_at,
        )
        if case_entry.case_record_path:
            write_case_record(output_dir / case_entry.case_record_path, record)
        cases.append(record)

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


@dataclass(frozen=True)
class PostprocessScriptInfo:
    path: str
    function_name: str
    description: str

    def to_json(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "function_name": self.function_name,
            "description": self.description,
        }


def _describe_postprocess_function(tree: ast.Module) -> tuple[str, str] | None:
    """Find the run_postprocessing() function node and its docstring.

    Returns (function_name, description) or None if no such function is
    defined at module level. Falls back to the module's own docstring when
    the function itself has none, since several kept tutorial scripts
    (table_summary.py) document themselves that way instead.
    """
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "run_postprocessing":
            doc = ast.get_docstring(node)
            if not doc:
                doc = ast.get_docstring(tree)
            return "run_postprocessing", (doc.strip() if doc else "(no description provided)")
    return None


def list_postprocess_scripts(setup_root: Path) -> tuple[PostprocessScriptInfo, ...]:
    """The catalog: every setup/ script exposing a run_postprocessing()
    function, each with a standard description. ``TutorialSpec.metadata`` is
    deliberately not a selector here: it is not persisted in RunDocument and
    cannot be the source of truth after a sweep has completed. The setup root
    recorded in that document is the durable discovery boundary. This lets
    the postprocessing
    module (or a reasoning agent) can see what's *available* and judge
    whether it applies, before running or reading anything further.

    Descriptions are read statically via `ast` (module/function docstrings)
    -- never imports or executes a candidate file, since setup/ scripts are
    untrusted tutorial content, not driver code. A script that fails to
    parse (SyntaxError) is skipped rather than breaking the whole catalog;
    that mirrors this codebase's "cataloging is best-effort" convention
    (see registry.py's _registered_tutorial_entry).
    """
    setup_root = Path(setup_root)
    if not setup_root.is_dir():
        return ()

    catalog: list[PostprocessScriptInfo] = []
    for path in sorted(setup_root.glob("*.py")):
        try:
            source = path.read_text()
        except OSError:
            continue
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue
        described = _describe_postprocess_function(tree)
        if described is None:
            continue
        function_name, description = described
        catalog.append(PostprocessScriptInfo(
            path=str(path), function_name=function_name, description=description,
        ))
    return tuple(catalog)


def run_postprocessing_module(context: SweepContext, *, task: str) -> PostprocessOutcome:
    """Placeholder postprocessing module. Still no real analysis -- but now
    driven entirely by the brain's grounded SweepContext plus the task it
    was asked to do. Lists (does not yet invoke) each case's setup/
    postprocessing script catalog -- real dispatch (run whichever script
    fits the task, or reason freely when none does) lands here next."""
    case_summaries = []
    for case in context.cases:
        catalog = list_postprocess_scripts(Path(case.setup_root)) if case.setup_root else ()
        if catalog:
            names = ", ".join(Path(info.path).name for info in catalog)
            found = f"{len(catalog)} script(s) in catalog: {names}"
        else:
            found = "no setup/ scripts in catalog"
        case_summaries.append(f"{case.case_id}: {found}, {len(case.output_files)} file(s)")
    return PostprocessOutcome(
        status="stub",
        message=(
            f"postprocess stub: sweep {context.sweep_spec_hash} "
            f"({context.completed_count}/{context.case_count} completed), task={task!r} -- "
            + "; ".join(case_summaries)
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
