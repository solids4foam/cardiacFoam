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
#     Adapts an agent-authored RunDocument v2 into the inputs the strict
#     workflow executor consumes. The producer-side counterpart to
#     strict_planning.strict_plan: instead of deriving the plan from the
#     on-disk case, it executes the plan the agent already authored.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Load and adapt a RunDocument v2 for strict workflow execution.

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
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .models import DataArtifact, data_artifact_from_json
from .run_model import RunDocument
from .workflow import normalize_workflow_dag, validate_workflow_commands
from .workflow_state import (
    WorkflowRunState,
    initial_workflow_state,
    workflow_state_from_json,
)
from ...specs.validation import validate_run


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

    A version-1 document is migrated to v2 via the explicit migration path;
    a version-2 document is validated against ``schemas/run-document.json``.
    Raises ``ValueError`` / ``json.JSONDecodeError`` on malformed input.
    """
    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError("Run document must be a JSON object")
    if data.get("version") == "1":
        return RunDocument.migrate_v1(data)
    return RunDocument.from_json(data)


def _diag(level: str, code: str, message: str, field: str = "") -> dict[str, Any]:
    return {"level": level, "code": code, "message": message, "field": field}


def build_execution_inputs(
    run_doc: RunDocument,
    *,
    utility_produces: dict[str, tuple[str, ...]] | None = None,
) -> tuple[RunDocumentExecutionInputs | None, tuple[dict[str, Any], ...]]:
    """Adapt ``run_doc`` into executor inputs.

    Returns ``(inputs, diagnostics)``. ``inputs`` is ``None`` whenever any
    error-level diagnostic is present (config invalid, no/invalid workflow
    DAG, disallowed command, no launch paths, unparseable artifact/state).
    ``diagnostics`` is always the full list, using the shape ``{level, code,
    message, field}``.
    """
    diagnostics: list[dict[str, Any]] = []

    # 1) Config validity against the live dict-entry catalog.
    for err in validate_run(run_doc):
        message = f"[{err.phase}] {err.message}" if err.phase else err.message
        diagnostics.append(_diag(err.level, "run_validation", message, err.field))

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
        expected_artifacts=tuple(expected_artifacts),
        utility_produces=utility_produces,
    )
    for d in wf_diagnostics:
        diagnostics.append(_diag(d.level, d.code, d.message, d.field))

    # 4) Command allowlist — same gate as the --entry path.
    for d in validate_workflow_commands(dag):
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
        or not case_root_raw
        or not output_dir_raw
    )
    if blocked:
        return None, tuple(diagnostics)

    inputs = RunDocumentExecutionInputs(
        workflow_dag=dag,
        workflow_state=workflow_state,
        case_root=Path(case_root_raw),
        output_dir=Path(output_dir_raw),
        expected_artifacts=tuple(expected_artifacts),
        run_document=run_doc,
    )
    return inputs, tuple(diagnostics)
