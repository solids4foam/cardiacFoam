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
#     run_document_adapter
#
# Description
#     Builds RunDocument payloads from strict-planned case state.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from dataclasses import asdict
from typing import Any

from openfoam_driver.core.planning_types import StrictDiagnostic, artifact_to_json, diagnostic
from openfoam_driver.core.specs.validation import validate_run
from ..plugin_capabilities import RunDocumentConfigurationRequest
from .models import DataArtifact
from .run_model import RunDocument
from .workflow_state import WorkflowRunState


def _run_document_from_case(
    *,
    entry: str,
    spec,
    launch: dict[str, Any],
    workflow_dag: dict[str, Any] | None,
    workflow_state: WorkflowRunState | None,
    expected_artifacts: tuple[DataArtifact, ...],
    driver_context,
) -> tuple[RunDocument, tuple[StrictDiagnostic, ...]]:
    config, configuration_diagnostics = (
        driver_context.capabilities.run_document_configuration.build(
            RunDocumentConfigurationRequest(spec),
        )
    )
    diagnostics: list[StrictDiagnostic] = list(configuration_diagnostics)
    generic_case = bool(spec.metadata.get("generic_case")) if spec.metadata else False

    if not generic_case:
        import jsonschema

        config_schema = driver_context.capabilities.run_document_configuration.schema()
        try:
            jsonschema.validate(config, config_schema)
        except jsonschema.exceptions.ValidationError as exc:
            diagnostics.append(diagnostic(
                "error",
                "plugin_config_schema_violation",
                f"Plugin-declared config schema rejected the built config: {exc.message}",
                field=".".join(str(part) for part in exc.absolute_path) or "",
            ))

    run_doc = RunDocument(
        id=f"plan-{entry}",
        name=entry,
        status="planned" if not diagnostics else "failed",
        intent={"source": "strict_plan"},
        plugin=driver_context.identity.to_json(),
        config=config,
        resolvedEntry={
            "entry": entry,
            "entryKind": spec.metadata.get("entry_kind"),
            "entryPath": spec.metadata.get("entry_path"),
            "resolvedName": spec.metadata.get("entry_name", entry),
            "sourceType": spec.metadata.get("source_type"),
            "workflowFamily": spec.metadata.get("workflow_family"),
            "isRunnable": True,
        },
        workflowDag=workflow_dag,
        workflowState=workflow_state.to_json() if workflow_state else None,
        launch={
            "action": launch.get("action"),
            "command": launch.get("command", []),
            "commandDisplay": launch.get("command_display", ""),
            "workflowStatePath": launch.get("workflow_state_path"),
            "caseRoot": launch.get("case_root"),
            "setupRoot": launch.get("setup_root"),
            "outputDir": launch.get("output_dir"),
        },
        expectedArtifacts=[artifact_to_json(artifact) for artifact in expected_artifacts],
        validation={"status": "not_run", "diagnostics": []},
    )
    # A core generic case has no solver-specific config to validate. Every
    # other case uses the explicit planning context rather than an ambient
    # cardiac compatibility default.
    validator_errors = [] if generic_case else validate_run(
        run_doc, driver_context=driver_context,
    )
    for error in validator_errors:
        diagnostics.append(diagnostic(
            error.level,
            "run_validation",
            error.message,
            field=error.field,
            source=error.phase,
        ))
    run_doc.validation = {
        "status": "ok" if not diagnostics else "failed",
        "diagnostics": [asdict(d) for d in diagnostics],
    }
    run_doc.status = "planned" if not any(d.level == "error" for d in diagnostics) else "failed"
    return run_doc, tuple(diagnostics)
