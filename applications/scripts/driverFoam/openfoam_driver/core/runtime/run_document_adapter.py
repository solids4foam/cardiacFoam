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
from pathlib import Path
from typing import Any

from ...planning_types import StrictDiagnostic, artifact_to_json, diagnostic
from ...specs.dict_builder import (
    build_electro_properties,
    build_physics_properties,
    parse_electro_properties,
    populate_values,
    resolve_context,
    select_applicable_entries,
)
from ...specs.validation import primary_phase, slot_key, validate_run
from .models import DataArtifact
from .run_model import RunDocument
from .workflow_state import WorkflowRunState


def _read_physics_type(path: Path) -> str | None:
    if not path.exists():
        return None
    for line in path.read_text().splitlines():
        stripped = line.split("//", 1)[0].strip()
        if not stripped.startswith("type"):
            continue
        tokens = stripped.rstrip(";").split()
        if len(tokens) >= 2:
            return tokens[1]
    return None


def _run_document_from_case(
    *,
    entry: str,
    spec,
    launch: dict[str, Any],
    workflow_dag: dict[str, Any] | None,
    workflow_state: WorkflowRunState | None,
    expected_artifacts: tuple[DataArtifact, ...],
) -> tuple[RunDocument, tuple[StrictDiagnostic, ...]]:
    diagnostics: list[StrictDiagnostic] = []
    config: dict[str, dict[str, Any]] = {
        "anatomy": {},
        "physics": {},
        "stimulus": {},
        "solver": {},
    }
    case_root = Path(spec.case_root)
    electro_path = case_root / "constant" / "electroProperties"
    physics_path = case_root / "constant" / "physicsProperties"
    physics_type = _read_physics_type(physics_path)
    if physics_type is None:
        diagnostics.append(diagnostic(
            "error",
            "missing_physics_properties",
            f"Could not read physicsProperties type from {physics_path}",
            source=str(physics_path),
            field="type",
        ))
    else:
        config["physics"]["type"] = physics_type
        try:
            build_physics_properties({"type": physics_type})
        except Exception as exc:
            diagnostics.append(diagnostic(
                "error",
                "invalid_physics_properties",
                str(exc),
                source=str(physics_path),
            ))

    if not electro_path.exists():
        diagnostics.append(diagnostic(
            "error",
            "missing_electro_properties",
            f"Missing electroProperties at {electro_path}",
            source=str(electro_path),
        ))
    else:
        try:
            parsed = parse_electro_properties(electro_path)
            selectors = parsed["selectors"]
            overrides = parsed.get("overrides", {})
            try:
                build_electro_properties(selectors, overrides=overrides or None)
            except Exception as exc:
                diagnostics.append(diagnostic(
                    "error",
                    "invalid_electro_properties",
                    str(exc),
                    source=str(electro_path),
                ))
            context = resolve_context(selectors, overrides=overrides or None)
            applicable_entries = select_applicable_entries(context)
            populated = populate_values(applicable_entries, context)
            for entry_obj in applicable_entries:
                key = slot_key(entry_obj.driver_path)
                if entry_obj.dynamic_path and key not in context:
                    continue
                if key not in populated:
                    continue
                phase = primary_phase(entry_obj) or "physics"
                config[phase][key] = populated[key]
        except Exception as exc:
            diagnostics.append(diagnostic(
                "error",
                "unparseable_electro_properties",
                str(exc),
                source=str(electro_path),
            ))

    run_doc = RunDocument(
        id=f"plan-{entry}",
        name=entry,
        status="planned" if not diagnostics else "failed",
        intent={"source": "strict_plan"},
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
            "manifestPath": launch.get("manifest_path"),
            "caseRoot": launch.get("case_root"),
            "setupRoot": launch.get("setup_root"),
            "outputDir": launch.get("output_dir"),
        },
        expectedArtifacts=[artifact_to_json(artifact) for artifact in expected_artifacts],
        validation={"status": "not_run", "diagnostics": []},
    )
    validator_errors = validate_run(run_doc)
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

