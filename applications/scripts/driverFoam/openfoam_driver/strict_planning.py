from __future__ import annotations

import os
import shutil
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .core.runtime.artifacts import predict_data_artifacts
from .core.runtime.models import DataArtifact
from .core.runtime.registry import load_entry_spec
from .core.runtime.run_model import RunDocument
from .core.runtime.workflow import WorkflowDiagnostic, normalize_workflow_dag
from .core.runtime.workflow_state import WorkflowRunState, initial_workflow_state
from .ionic_model_catalog import IONIC_MODEL_CATALOG
from .launch import describe_launch
from .scripts._dict_keys_scanner import strict_dict_key_report
from .specs.common import detect_ionic_model_name, detect_myocardium_solver_name
from .specs.dict_builder import (
    build_electro_properties,
    build_physics_properties,
    parse_electro_properties,
    populate_values,
    resolve_context,
    select_applicable_entries,
)
from .specs.validation import primary_phase, slot_key, validate_run
from .utility_catalog import UTILITY_CATALOG


@dataclass(frozen=True)
class StrictDiagnostic:
    level: str
    code: str
    message: str
    source: str = ""
    field: str = ""


@dataclass(frozen=True)
class StrictPlanReport:
    status: str
    entry: str
    resolved_entry: dict[str, Any]
    validation_diagnostics: tuple[StrictDiagnostic, ...] = ()
    workflow_diagnostics: tuple[StrictDiagnostic, ...] = ()
    catalog_coverage_errors: tuple[StrictDiagnostic, ...] = ()
    artifact_diagnostics: tuple[StrictDiagnostic, ...] = ()
    launch: dict[str, Any] = field(default_factory=dict)
    workflow_dag: dict[str, Any] | None = None
    workflow_state: WorkflowRunState | None = None
    expected_artifacts: tuple[DataArtifact, ...] = ()
    run_document: RunDocument | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "entry": self.entry,
            "resolved_entry": self.resolved_entry,
            "validation_diagnostics": [asdict(d) for d in self.validation_diagnostics],
            "workflow_diagnostics": [asdict(d) for d in self.workflow_diagnostics],
            "catalog_coverage_errors": [asdict(d) for d in self.catalog_coverage_errors],
            "artifact_diagnostics": [asdict(d) for d in self.artifact_diagnostics],
            "launch": self.launch,
            "workflow_dag": self.workflow_dag,
            "workflow_state": self.workflow_state.to_json() if self.workflow_state else None,
            "expected_artifacts": [_artifact_to_json(a) for a in self.expected_artifacts],
            "run_document": self.run_document.to_json() if self.run_document else None,
        }


_OPENFOAM_OR_DRIVER_COMMANDS = frozenset(
    {
        "Allrun",
        "blockMesh",
        "cardiacFoam",
        "decomposePar",
        "postProcess",
        "reconstructPar",
        "setExprFields",
        "topoSet",
    }
)


def _repo_root_from_here() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "src").exists() and (parent / "tutorials").exists():
            return parent
    raise RuntimeError("Could not locate repository root from strict_planning.py")


def _diagnostic(level: str, code: str, message: str, *, source: str = "", field: str = "") -> StrictDiagnostic:
    return StrictDiagnostic(level=level, code=code, message=message, source=source, field=field)


def _artifact_to_json(artifact: DataArtifact) -> dict[str, Any]:
    payload = asdict(artifact)
    payload["variables"] = list(artifact.variables)
    return payload


def _workflow_diagnostic_to_strict(diagnostic: WorkflowDiagnostic) -> StrictDiagnostic:
    return _diagnostic(
        diagnostic.level,
        diagnostic.code,
        diagnostic.message,
        source="workflow_dag",
        field=diagnostic.field,
    )


def _utility_produces_by_command() -> dict[str, tuple[str, ...]]:
    return {
        command: tuple(produce.artifact_id for produce in manifest.produces)
        for command, manifest in UTILITY_CATALOG.items()
        if manifest.produces
    }


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
        diagnostics.append(_diagnostic(
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
            diagnostics.append(_diagnostic(
                "error",
                "invalid_physics_properties",
                str(exc),
                source=str(physics_path),
            ))

    if not electro_path.exists():
        diagnostics.append(_diagnostic(
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
                diagnostics.append(_diagnostic(
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
            diagnostics.append(_diagnostic(
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
        expectedArtifacts=[_artifact_to_json(artifact) for artifact in expected_artifacts],
        validation={"status": "not_run", "diagnostics": []},
    )
    validator_errors = validate_run(run_doc)
    for error in validator_errors:
        diagnostics.append(_diagnostic(
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


def _artifact_diagnostics(
    spec,
    artifacts: tuple[DataArtifact, ...],
    workflow_dag: dict[str, Any] | None,
) -> tuple[StrictDiagnostic, ...]:
    diagnostics: list[StrictDiagnostic] = []
    case_root = Path(spec.case_root)
    electro_path = case_root / "constant" / "electroProperties"

    if not artifacts:
        diagnostics.append(_diagnostic(
            "error",
            "empty_artifact_prediction",
            "Strict planning could not predict any artifacts for this entry.",
            source=str(case_root),
        ))

    if electro_path.exists():
        try:
            solver = detect_myocardium_solver_name(electro_path)
            if solver not in {"singleCellSolver", "monodomainSolver", "bidomainSolver", "eikonalSolver"}:
                diagnostics.append(_diagnostic(
                    "error",
                    "unknown_solver",
                    f"No strict artifact handler is registered for myocardiumSolver {solver!r}.",
                    source=str(electro_path),
                    field="myocardiumSolver",
                ))
        except KeyError as exc:
            diagnostics.append(_diagnostic("error", "missing_solver", str(exc), source=str(electro_path)))

        try:
            ionic_model = detect_ionic_model_name(electro_path)
        except KeyError:
            ionic_model = None
        if ionic_model is not None and ionic_model not in IONIC_MODEL_CATALOG:
            diagnostics.append(_diagnostic(
                "error",
                "unknown_ionic_model",
                f"Ionic model {ionic_model!r} is not in IONIC_MODEL_CATALOG.",
                source=str(electro_path),
                field="ionicModel",
            ))

    for step in (workflow_dag or {}).get("steps", ()):
        command = step.get("command", "")
        if not command:
            diagnostics.append(_diagnostic(
                "error",
                "workflow_step_without_command",
                f"Workflow step {step.get('id', '<unknown>')!r} has no command.",
            ))
            continue
        if command in _OPENFOAM_OR_DRIVER_COMMANDS:
            continue
        manifest = UTILITY_CATALOG.get(command)
        if manifest is None:
            diagnostics.append(_diagnostic(
                "error",
                "unknown_workflow_command",
                f"Workflow command {command!r} is not a known OpenFOAM command or utility manifest.",
                field=str(step.get("id", "")),
            ))
        elif not manifest.produces:
            diagnostics.append(_diagnostic(
                "error",
                "utility_without_produces",
                f"Utility {command!r} has no authoritative produces entries.",
                field=str(step.get("id", "")),
            ))

    return tuple(diagnostics)


def _catalog_diagnostics(repo_root: Path) -> tuple[StrictDiagnostic, ...]:
    report = strict_dict_key_report(repo_root / "src")
    diagnostics: list[StrictDiagnostic] = []
    payload = report.to_json()
    for key in ("absent_keys", "stale_paths", "unmatched_subdicts", "unused_allowlist"):
        for item in payload[key]:
            diagnostics.append(_diagnostic(
                "error",
                f"dict_key_{key}",
                f"Strict dict-key scanner reported {key}: {item}",
                source="scan-dict-keys --strict",
            ))
    return tuple(diagnostics)


def _environment_diagnostics(spec) -> tuple[StrictDiagnostic, ...]:
    if "SKIP_ENV_DIAGNOSTICS" in os.environ:
        return ()
    diagnostics: list[StrictDiagnostic] = []
    if "WM_PROJECT_DIR" not in os.environ:
        diagnostics.append(_diagnostic(
            "error",
            "missing_openfoam_env",
            "WM_PROJECT_DIR is not set. OpenFOAM environment not sourced.",
            source="environment"
        ))
    if not shutil.which("cardiacFoam"):
        diagnostics.append(_diagnostic(
            "error",
            "missing_executable",
            "cardiacFoam not found on PATH.",
            source="environment",
            field="cardiacFoam"
        ))
    if not shutil.which("blockMesh"):
        diagnostics.append(_diagnostic(
            "error",
            "missing_executable",
            "blockMesh not found on PATH.",
            source="environment",
            field="blockMesh"
        ))
    return tuple(diagnostics)


def strict_plan(
    entry: str,
    *,
    entry_kind: str | None = None,
    overrides: dict[str, Any] | None = None,
    config_path: str | Path | None = None,
) -> StrictPlanReport:
    """Build a non-mutating strict simulation plan report."""
    spec = load_entry_spec(entry, entry_kind=entry_kind, overrides=overrides)
    launch = describe_launch(
        "sim",
        entry,
        entry_kind=entry_kind,
        overrides=overrides,
        config_path=config_path,
    )
    artifacts = tuple(predict_data_artifacts(Path(spec.case_root), spec))
    workflow_dag, workflow_diagnostics_raw = normalize_workflow_dag(
        spec.metadata.get("workflow_dag") if spec.metadata else None,
        expected_artifacts=artifacts,
        utility_produces=_utility_produces_by_command(),
    )
    workflow_diagnostics = tuple(
        _workflow_diagnostic_to_strict(diagnostic)
        for diagnostic in workflow_diagnostics_raw
    )
    workflow_state = initial_workflow_state(workflow_dag)
    run_document, validation_diagnostics = _run_document_from_case(
        entry=entry,
        spec=spec,
        launch=launch,
        workflow_dag=workflow_dag,
        workflow_state=workflow_state,
        expected_artifacts=artifacts,
    )
    repo_root = _repo_root_from_here()
    catalog_diagnostics = _catalog_diagnostics(repo_root)
    artifact_diagnostics = _artifact_diagnostics(spec, artifacts, workflow_dag)
    env_diagnostics = _environment_diagnostics(spec)
    all_diagnostics = (
        validation_diagnostics
        + workflow_diagnostics
        + catalog_diagnostics
        + artifact_diagnostics
        + env_diagnostics
    )
    failed = any(
        diagnostic.level == "error"
        for diagnostic in all_diagnostics
    )
    run_document.status = "failed" if failed else "planned"
    run_document.validation = {
        "status": "failed" if failed else "ok",
        "diagnostics": [asdict(diagnostic) for diagnostic in all_diagnostics],
    }
    return StrictPlanReport(
        status="failed" if failed else "ok",
        entry=entry,
        resolved_entry={
            "entry_name": spec.metadata.get("entry_name", entry),
            "entry_kind": spec.metadata.get("entry_kind"),
            "entry_path": spec.metadata.get("entry_path"),
            "source_type": spec.metadata.get("source_type"),
            "workflow_family": spec.metadata.get("workflow_family"),
        },
        validation_diagnostics=validation_diagnostics,
        workflow_diagnostics=workflow_diagnostics,
        catalog_coverage_errors=catalog_diagnostics,
        artifact_diagnostics=artifact_diagnostics,
        launch=launch,
        workflow_dag=workflow_dag,
        workflow_state=workflow_state,
        expected_artifacts=artifacts,
        run_document=run_document,
    )
