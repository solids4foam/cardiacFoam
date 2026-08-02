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
#     strict_planning
#
# Description
#     Evaluates execution plans against strict schema constraints.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import os
import shlex
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

def __get_capabilities():
    from openfoam_driver.core.plugin_interface import get_active_plugin
    return get_active_plugin().get_capabilities()


from .core.runtime.artifacts import predict_data_artifacts
from .core.runtime.environment_preflight import (
    _environment_diagnostics,
    _required_executables,
    _unwrap_mpi_program,
    shutil,
)
from .core.runtime.execution_context import resolve_execution_context
from .core.runtime.models import DataArtifact
from .core.runtime.registry import load_entry_spec
from .core.runtime.run_document_adapter import _run_document_from_case
from .core.runtime.run_model import RunDocument
from .core.runtime.strict_audit import _build_simulation_audit
from .core.runtime.workflow import (
    WorkflowDiagnostic,
    normalize_workflow_dag,
    validate_workflow_commands,
)
from .core.runtime.workflow_state import WorkflowRunState, initial_workflow_state
from .capability_manifest import build_capability_manifest, resolve_case_models
from .planning_types import (
    StrictDiagnostic,
    SimulationAuditItem,
    artifact_to_json as _artifact_to_json,
    diagnostic as _diagnostic,
)
from .scripts._dict_keys_scanner import strict_dict_key_report
from .specs.common import (
    detect_ionic_model_name,
    detect_myocardium_solver_name,
    detect_verification_model_type,
)
from .specs.function_object_fields import function_object_field_diagnostics
from .specs.mesh_geometry import mesh_geometry_diagnostics as _detect_mesh_geometry


@dataclass(frozen=True)
class StrictPlanReport:
    status: str
    entry: str
    resolved_entry: dict[str, Any]
    readiness_score: dict[str, Any] = field(default_factory=dict)
    simulation_audit: tuple[SimulationAuditItem, ...] = ()
    validation_diagnostics: tuple[StrictDiagnostic, ...] = ()
    workflow_diagnostics: tuple[StrictDiagnostic, ...] = ()
    catalog_coverage_errors: tuple[StrictDiagnostic, ...] = ()
    artifact_diagnostics: tuple[StrictDiagnostic, ...] = ()
    environment_diagnostics: tuple[StrictDiagnostic, ...] = ()
    mesh_geometry_diagnostics: tuple[StrictDiagnostic, ...] = ()
    launch: dict[str, Any] = field(default_factory=dict)
    workflow_dag: dict[str, Any] | None = None
    workflow_state: WorkflowRunState | None = None
    expected_artifacts: tuple[DataArtifact, ...] = ()
    run_document: RunDocument | None = None
    capability_manifest: dict[str, Any] = field(default_factory=dict)
    function_object_diagnostics: tuple[StrictDiagnostic, ...] = ()

    def to_json(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "entry": self.entry,
            "resolved_entry": self.resolved_entry,
            "readiness_score": self.readiness_score,
            "simulation_audit": [asdict(item) for item in self.simulation_audit],
            "validation_diagnostics": [asdict(d) for d in self.validation_diagnostics],
            "workflow_diagnostics": [asdict(d) for d in self.workflow_diagnostics],
            "catalog_coverage_errors": [asdict(d) for d in self.catalog_coverage_errors],
            "artifact_diagnostics": [asdict(d) for d in self.artifact_diagnostics],
            "environment_diagnostics": [asdict(d) for d in self.environment_diagnostics],
            "mesh_geometry_diagnostics": [
                asdict(d) for d in self.mesh_geometry_diagnostics
            ],
            "launch": self.launch,
            "workflow_dag": self.workflow_dag,
            "workflow_state": self.workflow_state.to_json() if self.workflow_state else None,
            "expected_artifacts": [_artifact_to_json(a) for a in self.expected_artifacts],
            "run_document": self.run_document.to_json() if self.run_document else None,
            "capability_manifest": self.capability_manifest,
            "function_object_diagnostics": [
                asdict(d) for d in self.function_object_diagnostics
            ],
        }


def _repo_root_from_here() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "src").exists() and (parent / "tutorials").exists():
            return parent
    raise RuntimeError("Could not locate repository root from strict_planning.py")


def _workflow_diagnostic_to_strict(diagnostic: WorkflowDiagnostic) -> StrictDiagnostic:
    return _diagnostic(
        diagnostic.level,
        diagnostic.code,
        diagnostic.message,
        source="workflow_dag",
        field=diagnostic.field,
    )


def _utility_produces_by_command() -> dict[str, tuple[str, ...]]:
    from .utility_catalog import UTILITY_CATALOG

    return {
        command: tuple(produce.artifact_id for produce in manifest.produces)
        for command, manifest in UTILITY_CATALOG.items()
        if manifest.produces
    }


def _artifact_diagnostics(
    spec,
    artifacts: tuple[DataArtifact, ...],
    workflow_dag: dict[str, Any] | None,
) -> tuple[StrictDiagnostic, ...]:
    from .core.runtime.workflow import validate_workflow_commands

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
        if ionic_model is not None and ionic_model not in __get_capabilities().get("ionic_models", {}):
            diagnostics.append(_diagnostic(
                "error",
                "unknown_ionic_model",
                f"Ionic model {ionic_model!r} is not supported by the active plugin..",
                source=str(electro_path),
                field="ionicModel",
            ))

    for diagnostic in validate_workflow_commands(workflow_dag):
        diagnostics.append(_diagnostic(
            diagnostic.level,
            diagnostic.code,
            diagnostic.message,
            field=diagnostic.field,
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


def _is_nondimensional_entry(spec) -> bool:
    """Return True when the SI mesh-scale gate is not meaningful."""
    entry_name = ""
    family = ""
    if spec.metadata:
        entry_name = str(spec.metadata.get("entry_name", "") or "")
        family = str(spec.metadata.get("workflow_family", "") or "")
    haystack = f"{entry_name} {family}".lower()
    if "manufactured" in haystack or "verification" in haystack:
        return True
    electro_path = Path(spec.case_root) / "constant" / "electroProperties"
    if electro_path.exists():
        try:
            if detect_myocardium_solver_name(electro_path) == "singleCellSolver":
                return True
            if detect_verification_model_type(electro_path) is not None:
                return True
        except Exception:
            pass
    return False


def _mesh_geometry_diagnostics(
    case_root: str | Path,
    *,
    exempt: bool = False,
) -> tuple[StrictDiagnostic, ...]:
    """Adapt mesh-scale detection into StrictDiagnostics for the report."""
    if exempt or "SKIP_MESH_DIAGNOSTICS" in os.environ:
        return ()
    return tuple(
        _diagnostic(
            d.level,
            d.code,
            d.message,
            source="mesh_geometry",
            field=d.region,
        )
        for d in _detect_mesh_geometry(Path(case_root))
    )


def _has_error(diagnostics: tuple[StrictDiagnostic, ...]) -> bool:
    return any(diagnostic.level == "error" for diagnostic in diagnostics)


def _run_launch_description(
    entry: str,
    context,
    *,
    entry_kind: str | None,
    config_path: str | Path | None,
) -> dict[str, Any]:
    """Describe the modern `run --strict --entry` invocation for this plan.

    Replaces strict_plan's former reuse of describe_launch("sim", ...):
    that call re-resolved the entry a second time (strict_plan already has
    `spec` from load_entry_spec) purely to read these four paths off it, and
    tied the strict/workflow-DAG path -- which never runs the legacy
    sim/post/all CLI at all -- to describe_launch's action vocabulary.
    `run --strict --entry` is the command that actually executes this exact
    plan today.
    """
    command = [sys.executable, "-m", "openfoam_driver", "run", "--strict", "--entry", entry]
    if entry_kind is not None:
        command.extend(["--entry-kind", entry_kind])
    if config_path is not None:
        command.extend(["--config", str(config_path)])
    return {
        "action": "run",
        "command": command,
        "command_display": shlex.join(command),
        "manifest_path": str(context.manifest_path),
        "case_root": str(context.case_root),
        "setup_root": str(context.setup_root),
        "output_dir": str(context.output_dir),
    }


def strict_plan(
    entry: str,
    *,
    entry_kind: str | None = None,
    overrides: dict[str, Any] | None = None,
    config_path: str | Path | None = None,
    openfoam_bashrc: str | Path | None = None,
) -> StrictPlanReport:
    """Build a non-mutating strict simulation plan report."""
    spec = load_entry_spec(entry, entry_kind=entry_kind, overrides=overrides)
    execution_context = resolve_execution_context(spec)
    launch = _run_launch_description(
        entry, execution_context, entry_kind=entry_kind, config_path=config_path,
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
    env_diagnostics = _environment_diagnostics(
        workflow_dag,
        openfoam_bashrc=str(openfoam_bashrc) if openfoam_bashrc is not None else None,
    )
    mesh_diagnostics = _mesh_geometry_diagnostics(
        spec.case_root, exempt=_is_nondimensional_entry(spec)
    )
    simulation_audit, generation_diagnostics, readiness_score = _build_simulation_audit(
        spec=spec,
        workflow_dag=workflow_dag,
        artifacts=artifacts,
        validation_diagnostics=validation_diagnostics,
        workflow_diagnostics=workflow_diagnostics,
        artifact_diagnostics=artifact_diagnostics,
        environment_diagnostics=env_diagnostics,
        mesh_geometry_diagnostics=mesh_diagnostics,
    )
    plan_diagnostics = (
        generation_diagnostics
        + validation_diagnostics
        + workflow_diagnostics
        + catalog_diagnostics
        + artifact_diagnostics
        + mesh_diagnostics
    )
    resolved_solver, resolved_ionic, resolved_active_tension = resolve_case_models(
        spec.case_root
    )
    capability_manifest = build_capability_manifest(
        resolved_solver=resolved_solver,
        resolved_ionic_model=resolved_ionic,
        resolved_active_tension=resolved_active_tension,
    )
    function_object_diagnostics = function_object_field_diagnostics(
        spec.case_root, samplable=capability_manifest["samplable_fields"]
    )
    # Field diagnostics are warn-only: reported (in all_diagnostics) but never
    # part of plan_diagnostics, so a sampled-field warning cannot fail a plan.
    all_diagnostics = (
        plan_diagnostics + env_diagnostics + function_object_diagnostics
    )
    failed = _has_error(plan_diagnostics)
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
        readiness_score=readiness_score,
        simulation_audit=simulation_audit,
        validation_diagnostics=validation_diagnostics,
        workflow_diagnostics=workflow_diagnostics,
        catalog_coverage_errors=catalog_diagnostics,
        artifact_diagnostics=artifact_diagnostics,
        environment_diagnostics=env_diagnostics,
        mesh_geometry_diagnostics=mesh_diagnostics,
        launch=launch,
        workflow_dag=workflow_dag,
        workflow_state=workflow_state,
        expected_artifacts=artifacts,
        run_document=run_document,
        capability_manifest=capability_manifest,
        function_object_diagnostics=function_object_diagnostics,
    )
