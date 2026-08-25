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
#     run_document_exec
#
# Description
#     Adapts an agent-authored RunDocument v3 into the inputs the strict
#     workflow executor consumes. The producer-side counterpart to
#     strict_planning.strict_plan: instead of deriving the plan from the
#     on-disk case, it executes the plan the agent already authored.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Load and adapt a RunDocument v3 for strict workflow execution.

``load_run_document(path)`` reads and schema-validates a document (migrating
v1 explicitly); raises ``ValueError`` or ``json.JSONDecodeError`` on malformed
input. ``build_execution_inputs(doc)`` turns a document into the same
``(workflow_dag, workflow_state, case_root, output_dir, expected_artifacts)``
tuple ``strict_plan`` produces, so the CLI run/step path is identical for both
producers. It enforces the same command allowlist as ``strict_plan``; anything
that makes the document non-executable is returned as a diagnostic with
``inputs is None``.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from .models import DataArtifact, data_artifact_from_json
# _case_is_runnable is a private helper reused as-is: the plan treats this as
# an accepted pragmatic tradeoff rather than promoting it to a public API
# (out of scope here).
from .registry import _case_is_runnable
from .run_model import RunDocument
from .workflow import normalize_workflow_dag, validate_workflow_commands, workflow_output_artifacts
from .workflow_state import (
    WorkflowRunState,
    initial_workflow_state,
    workflow_state_from_json,
)
from openfoam_driver.core.specs.validation import validate_run

if TYPE_CHECKING:
    from ..plugin_interface import DriverContext


@dataclass(frozen=True)
class RunDocumentExecutionInputs:
    """Everything the strict executor needs, derived from a RunDocument."""

    workflow_dag: dict[str, Any]
    workflow_state: WorkflowRunState
    case_root: Path
    output_dir: Path
    expected_artifacts: tuple[DataArtifact, ...]
    run_document: RunDocument


def load_run_document(path: str | Path) -> RunDocument:
    """Read, schema-validate, and return a RunDocument from ``path``.

    A version-1 document is migrated to v3 via the explicit migration path;
    a version-3 document is validated against ``schemas/run-document.json``.
    Version-2 documents are rejected -- callers must migrate them explicitly
    via ``RunDocument.migrate_v2`` before loading. Raises ``ValueError`` /
    ``json.JSONDecodeError`` on malformed input.
    """
    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError("Run document must be a JSON object")
    if data.get("version") == "1":
        return RunDocument.migrate_v1(data)
    return RunDocument.from_json(data)


def _diag(level: str, code: str, message: str, field: str = "") -> dict[str, Any]:
    return {"level": level, "code": code, "message": message, "field": field}


def _allowed_runs_root(env: dict[str, str] | None = None) -> Path | None:
    """Resolved DRIVERFOAM_ALLOWED_RUNS_ROOT, or None when unset/empty.

    Accepts an explicit ``env`` mapping so callers (including future tests)
    can inject the value without mutating ``os.environ``; the current tests
    use ``mock.patch.dict`` on real ``os.environ`` instead.
    """
    source = env if env is not None else os.environ
    value = source.get("DRIVERFOAM_ALLOWED_RUNS_ROOT")
    if not value:
        return None
    return Path(value).resolve()


def _validate_config_against_plugin_schema(
    run_doc: RunDocument,
    driver_context: "DriverContext | None",
    diagnostics: list[dict[str, Any]],
) -> None:
    """Append a ``plugin_config_schema_violation`` diagnostic per violation.

    Same check, code, and message shape as
    ``run_document_adapter._run_document_from_case`` uses on the emission
    path -- kept symmetric so a config the planner would refuse to emit is
    also a config the executor refuses to ingest.
    """
    import jsonschema

    from ..compatibility import resolve_public_driver_context

    context = resolve_public_driver_context(driver_context)
    config_schema = context.capabilities.run_document_configuration.schema()
    try:
        jsonschema.validate(run_doc.config, config_schema)
    except jsonschema.exceptions.ValidationError as exc:
        diagnostics.append(_diag(
            "error",
            "plugin_config_schema_violation",
            f"Plugin-declared config schema rejected the run document config: "
            f"{exc.message}",
            ".".join(str(part) for part in exc.absolute_path) or "config",
        ))


def build_execution_inputs(
    run_doc: RunDocument,
    *,
    utility_produces: dict[str, tuple[str, ...]] | None = None,
    driver_context: "DriverContext | None" = None,
) -> tuple[RunDocumentExecutionInputs | None, tuple[dict[str, Any], ...]]:
    """Adapt ``run_doc`` into executor inputs.

    Returns ``(inputs, diagnostics)``. ``inputs`` is ``None`` whenever any
    error-level diagnostic is present (config invalid, no/invalid workflow
    DAG, disallowed command, no launch paths, unparseable artifact/state).
    ``diagnostics`` is always the full list, using the shape ``{level, code,
    message, field}``.
    """
    diagnostics: list[dict[str, Any]] = []

    if driver_context is not None and run_doc.plugin is not None:
        planned = run_doc.plugin
        selected = driver_context.identity.to_json()
        mismatched = [
            key for key in ("id", "version", "api_version", "capability_digest")
            if planned.get(key) != selected.get(key)
        ]
        if mismatched:
            diagnostics.append(_diag(
                "error",
                "plugin_identity_mismatch",
                "RunDocument plugin does not match the supplied driver context: "
                + ", ".join(mismatched),
                "plugin",
            ))

    # 1) Config validity against the selected plugin's live dictionary
    # catalog and semantic validators. This must use the same immutable
    # context whose identity was checked by the CLI before execution.
    for err in validate_run(run_doc, driver_context=driver_context):
        message = f"[{err.phase}] {err.message}" if err.phase else err.message
        diagnostics.append(_diag(err.level, "run_validation", message, err.field))

    # 1b) Plugin-declared config schema. run_document_adapter applies this on
    # the *emission* path (config the plugin just built); an ingested,
    # agent-authored document is the untrusted counterpart and must clear the
    # same gate, with the same `plugin_config_schema_violation` code.
    #
    # No generic-case exemption here: unlike run_document_adapter (which reads
    # spec.metadata["generic_case"]), a RunDocument carries no generic-case
    # marker -- resolvedEntry.entryKind/sourceType do not distinguish a core
    # generic case from any other discovered case folder. This mirrors the
    # unconditional validate_run above, and a generic case's all-empty config
    # satisfies the plugin schemas in-tree regardless.
    _validate_config_against_plugin_schema(run_doc, driver_context, diagnostics)

    # 2) Expected artifacts: reconstruct, reporting any malformed entry.
    expected_artifacts: list[DataArtifact] = []
    raw_artifacts = run_doc.expectedArtifacts
    if not isinstance(raw_artifacts, (list, tuple)):
        diagnostics.append(_diag(
            "error", "invalid_expected_artifacts",
            f"expectedArtifacts must be a list, got {type(raw_artifacts).__name__}.",
            "expectedArtifacts",
        ))
        raw_artifacts = []
    for raw in raw_artifacts:
        try:
            expected_artifacts.append(data_artifact_from_json(raw))
        except Exception as exc:  # malformed shape / unknown placeholder
            diagnostics.append(
                _diag("error", "invalid_expected_artifact", str(exc), "expectedArtifacts")
            )

    # 3) Re-normalize the supplied DAG so a hand-authored workflow gets the
    #    same shape guarantees and diagnostics as a strict-planned one.
    dag, wf_diagnostics = normalize_workflow_dag(
        run_doc.workflowDag,
        expected_artifacts=workflow_output_artifacts(expected_artifacts),
        utility_produces=utility_produces,
        driver_context=driver_context,
    )
    for d in wf_diagnostics:
        diagnostics.append(_diag(d.level, d.code, d.message, d.field))

    # 4) Command allowlist — same gate as the --entry path.
    for d in validate_workflow_commands(dag, driver_context=driver_context):
        diagnostics.append(_diag(d.level, d.code, d.message, d.field))

    # 5) Launch paths are mandatory for execution.
    raw_launch = run_doc.launch
    if raw_launch is None:
        launch: dict[str, Any] = {}
    elif isinstance(raw_launch, dict):
        launch = raw_launch
    else:
        launch = {}
        diagnostics.append(_diag(
            "error", "invalid_launch",
            f"Run document launch must be a JSON object, got {type(raw_launch).__name__}.",
            "launch",
        ))
    case_root_raw = launch.get("caseRoot")
    output_dir_raw = launch.get("outputDir")
    if not case_root_raw:
        diagnostics.append(_diag(
            "error", "missing_case_root",
            "Run document launch.caseRoot is required for execution.",
            "launch.caseRoot",
        ))
    if not output_dir_raw:
        diagnostics.append(_diag(
            "error", "missing_output_dir",
            "Run document launch.outputDir is required for execution.",
            "launch.outputDir",
        ))

    # 5b) Validate + canonicalize launch paths against the OpenFOAM layout.
    #     caseRoot is the OpenFOAM/solver output base and must be a runnable
    #     case. Resolve (follow symlinks) so all downstream checks and the
    #     artifact gate operate on one canonical absolute path.
    resolved_case_root: Path | None = None
    if case_root_raw:
        resolved_case_root = Path(case_root_raw).resolve()
        if not resolved_case_root.exists():
            diagnostics.append(_diag(
                "error", "case_root_missing",
                f"Run document launch.caseRoot does not exist: {case_root_raw}.",
                "launch.caseRoot",
            ))
            resolved_case_root = None
        elif not resolved_case_root.is_dir():
            diagnostics.append(_diag(
                "error", "case_root_not_a_directory",
                f"Run document launch.caseRoot is not a directory: {case_root_raw}.",
                "launch.caseRoot",
            ))
            resolved_case_root = None
        elif not _case_is_runnable(
            resolved_case_root, driver_context=driver_context,
        ):
            diagnostics.append(_diag(
                "error", "case_root_not_a_runnable_case",
                f"Run document launch.caseRoot is not a runnable OpenFOAM case "
                f"(missing required constant/ or system/ files): {case_root_raw}.",
                "launch.caseRoot",
            ))
            resolved_case_root = None

    # outputDir is the driver-bookkeeping base. Resolve the same way the driver
    # builds it (resolve_spec_paths): absolute as-is, relative under caseRoot.
    resolved_output_dir: Path | None = None
    if output_dir_raw:
        candidate = Path(output_dir_raw)
        if candidate.is_absolute():
            resolved_output_dir = candidate.resolve()
        elif resolved_case_root is not None:
            resolved_output_dir = (resolved_case_root / candidate).resolve()
        else:
            # resolved_case_root is None here (caseRoot missing/invalid), which
            # already forces `blocked = True` below regardless of outputDir.
            # This CWD-relative resolution is never actually used for
            # execution — computed only so every branch yields a value.
            resolved_output_dir = candidate.resolve()

    # Opt-in hard boundary: when DRIVERFOAM_ALLOWED_RUNS_ROOT is set, both
    # resolved paths must sit under it. Runs on resolved paths, so a symlink or
    # an absolute output_dir cannot escape. Unset -> skipped (no behavior change).
    allowed_root = _allowed_runs_root()
    if allowed_root is not None:
        if resolved_case_root is not None and not resolved_case_root.is_relative_to(allowed_root):
            diagnostics.append(_diag(
                "error", "case_root_outside_allowed_root",
                f"launch.caseRoot resolves outside DRIVERFOAM_ALLOWED_RUNS_ROOT "
                f"({allowed_root}): {resolved_case_root}.",
                "launch.caseRoot",
            ))
        if resolved_output_dir is not None and not resolved_output_dir.is_relative_to(allowed_root):
            diagnostics.append(_diag(
                "error", "output_dir_outside_allowed_root",
                f"launch.outputDir resolves outside DRIVERFOAM_ALLOWED_RUNS_ROOT "
                f"({allowed_root}): {resolved_output_dir}.",
                "launch.outputDir",
            ))

    # 6) Workflow state: prefer the document's snapshot, else derive from DAG.
    workflow_state: WorkflowRunState | None
    if run_doc.workflowState is not None:
        try:
            workflow_state = workflow_state_from_json(run_doc.workflowState)
        except Exception as exc:
            workflow_state = None
            diagnostics.append(
                _diag("error", "invalid_workflow_state", str(exc), "workflowState")
            )
    else:
        # Computed regardless of earlier errors so all diagnostics are gathered
        # before the single blocked check below.
        workflow_state = initial_workflow_state(dag) if dag is not None else None

    blocked = (
        any(d["level"] == "error" for d in diagnostics)
        or dag is None
        or workflow_state is None
        or resolved_case_root is None
        or resolved_output_dir is None
    )
    if blocked:
        return None, tuple(diagnostics)

    inputs = RunDocumentExecutionInputs(
        workflow_dag=dag,
        workflow_state=workflow_state,
        case_root=resolved_case_root,
        output_dir=resolved_output_dir,
        expected_artifacts=tuple(expected_artifacts),
        run_document=run_doc,
    )
    return inputs, tuple(diagnostics)
