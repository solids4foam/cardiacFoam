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
from dataclasses import asdict, dataclass, field, is_dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .plugin_interface import DriverContext


from .runtime.artifacts import predict_data_artifacts
from .runtime.environment_preflight import (
    _environment_diagnostics,
    _required_executables,
    _unwrap_mpi_program,
    shutil,
)
from .runtime.execution_context import resolve_execution_context
from .runtime.models import DataArtifact
from .runtime.registry import load_entry_spec
from .runtime.run_document_adapter import _run_document_from_case
from .runtime.run_model import RunDocument
from .runtime.strict_audit import _build_simulation_audit
from .runtime.workflow import (
    WorkflowDiagnostic,
    normalize_workflow_dag,
    validate_workflow_commands,
    workflow_output_artifacts,
)
from .runtime.workflow_state import WorkflowRunState, initial_workflow_state
from openfoam_driver.core.planning_types import (
    StrictDiagnostic,
    SimulationAuditItem,
    artifact_to_json as _artifact_to_json,
    diagnostic as _diagnostic,
)
from ..scripts._dict_keys_scanner import strict_dict_key_report
from openfoam_driver.core.specs.function_object_fields import function_object_field_diagnostics
from openfoam_driver.core.specs.mesh_geometry import mesh_geometry_diagnostics as _detect_mesh_geometry


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
    plugin: dict[str, str] = field(default_factory=dict)

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
            "plugin": self.plugin,
        }


def _jsonable(value: Any) -> Any:
    """Convert plugin capability metadata into a report-safe value."""
    if is_dataclass(value):
        return _jsonable(asdict(value))
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, set):
        return sorted(_jsonable(item) for item in value)
    return value


def _repo_root_from_here() -> Path | None:
    """Return the monorepo root, or None when running in a standalone install.

    Mirrors the three-tier logic in ``specs/paths.py:repo_root_default()``.
    Returns ``None`` instead of raising so that callers can gracefully skip
    operations that require the full source tree (e.g. dict-key scanning).
    """
    current = Path(__file__).resolve()
    tier2_candidate: Path | None = None
    for parent in current.parents:
        has_src = (parent / "src").exists()
        has_tutorials = (parent / "tutorials").exists()
        if has_src and has_tutorials:   # Tier 1: full monorepo
            return parent
        if has_tutorials and tier2_candidate is None:  # Tier 2: tutorials-only
            tier2_candidate = parent
    return tier2_candidate  # Tier 3: fully standalone → None


def _workflow_diagnostic_to_strict(diagnostic: WorkflowDiagnostic) -> StrictDiagnostic:
    return _diagnostic(
        diagnostic.level,
        diagnostic.code,
        diagnostic.message,
        source="workflow_dag",
        field=diagnostic.field,
    )


def _utility_produces_by_command(
    driver_context: "DriverContext",
) -> dict[str, tuple[str, ...]]:
    """The active plugin's utilities, keyed to the artifacts they declare."""

    from openfoam_driver.core.capability_manifest import utility_produces

    return utility_produces(
        driver_context.capabilities.command_authorization.utility_manifests()
    )


def _artifact_diagnostics(
    spec,
    artifacts: tuple[DataArtifact, ...],
    workflow_dag: dict[str, Any] | None,
    driver_context: "DriverContext",
) -> tuple[StrictDiagnostic, ...]:
    from .runtime.workflow import validate_workflow_commands

    diagnostics: list[StrictDiagnostic] = []
    case_root = Path(spec.case_root)

    if not artifacts:
        diagnostics.append(_diagnostic(
            "error",
            "empty_artifact_prediction",
            "Strict planning could not predict any artifacts for this entry.",
            source=str(case_root),
        ))

    # Defer domain-specific validation to the selected capability while
    # preserving the public plugin call and diagnostic order.
    from .plugin_capabilities import ConfigurationValidationRequest

    diagnostics.extend(
        driver_context.capabilities.configuration_validator.validate(
            ConfigurationValidationRequest(spec),
        )
    )

    for diagnostic in validate_workflow_commands(
        workflow_dag, driver_context=driver_context,
    ):
        diagnostics.append(_diagnostic(
            diagnostic.level,
            diagnostic.code,
            diagnostic.message,
            field=diagnostic.field,
        ))

    return tuple(diagnostics)


def _catalog_diagnostics(driver_context: "DriverContext") -> tuple[StrictDiagnostic, ...]:
    """Run only the active plugin's reviewed C++↔Python mapping checks."""

    mapping = driver_context.capabilities.cxx_mapping.profile().cxx_mapping
    if mapping is None:
        return ()
    diagnostics: list[StrictDiagnostic] = []
    for source_root in mapping.source_roots:
        if not source_root.is_dir():
            diagnostics.append(_diagnostic(
                "warning",
                "plugin_cxx_source_unavailable",
                f"Plugin C++ source root is unavailable: {source_root}",
                source=driver_context.identity.id,
            ))
            continue
        report = strict_dict_key_report(
            source_root,
            allowlist_path=mapping.allowlist_path,
            entries=driver_context.capabilities.dictionaries.entries(),
        )
        payload = report.to_json()
        for key in ("unmatched_cxx_reads", "stale_paths", "unmatched_subdicts", "unused_allowlist"):
            for item in payload[key]:
                diagnostics.append(_diagnostic(
                    "error",
                    f"plugin_dict_key_{key}",
                    f"Plugin C++/catalog scanner reported {key}: {item}",
                    source=f"{driver_context.identity.id}:{source_root}",
                ))
    return tuple(diagnostics)


def _is_nondimensional_entry(spec, driver_context=None) -> bool:
    """Return True when the SI mesh-scale gate is not meaningful."""
    entry_name = ""
    family = ""
    if spec.metadata:
        entry_name = str(spec.metadata.get("entry_name", "") or "")
        family = str(spec.metadata.get("workflow_family", "") or "")
    haystack = f"{entry_name} {family}".lower()
    if "manufactured" in haystack or "verification" in haystack:
        return True
    from .compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    return driver_context.capabilities.mesh_diagnostic_policy.is_nondimensional(spec)


def _mesh_geometry_diagnostics(
    case_root: str | Path,
    *,
    exempt: bool = False,
    driver_context: "DriverContext | None" = None,
) -> tuple[StrictDiagnostic, ...]:
    """Adapt mesh-scale detection into StrictDiagnostics for the report.

    Core classifies every polyMesh region's scale; the active plugin may add
    checks for point sets that are not mesh regions (cardiacFoam's
    ``constant/purkinjeGraph*``). Both report under the same
    ``mesh_geometry`` source, and both are skipped by the same exemption.
    """
    if exempt or "SKIP_MESH_DIAGNOSTICS" in os.environ:
        return ()
    from .compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    detected = list(_detect_mesh_geometry(Path(case_root)))
    detected.extend(
        driver_context.capabilities.mesh_diagnostic_policy.extra_geometry_diagnostics(
            Path(case_root),
        )
    )
    return tuple(
        _diagnostic(
            d.level,
            d.code,
            d.message,
            source="mesh_geometry",
            field=d.region,
        )
        for d in detected
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
        "workflow_state_path": str(context.workflow_state_path),
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
    driver_context: "DriverContext | None" = None,
) -> StrictPlanReport:
    """Build a non-mutating strict simulation plan report."""
    from .compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    spec = load_entry_spec(
        entry,
        entry_kind=entry_kind,
        overrides=overrides,
        driver_context=driver_context,
    )
    execution_context = resolve_execution_context(spec)
    launch = _run_launch_description(
        entry, execution_context, entry_kind=entry_kind, config_path=config_path,
    )
    artifacts = tuple(
        predict_data_artifacts(
            Path(spec.case_root), spec, driver_context=driver_context,
        )
    )
    workflow_dag, workflow_diagnostics_raw = normalize_workflow_dag(
        spec.metadata.get("workflow_dag") if spec.metadata else None,
        expected_artifacts=workflow_output_artifacts(artifacts),
        utility_produces=_utility_produces_by_command(driver_context),
        driver_context=driver_context,
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
        driver_context=driver_context,
    )
    catalog_diagnostics = _catalog_diagnostics(driver_context)
    artifact_diagnostics = _artifact_diagnostics(
        spec, artifacts, workflow_dag, driver_context,
    )
    env_diagnostics = _environment_diagnostics(
        workflow_dag,
        openfoam_bashrc=str(openfoam_bashrc) if openfoam_bashrc is not None else None,
        driver_context=driver_context,
    )
    mesh_diagnostics = _mesh_geometry_diagnostics(
        spec.case_root,
        exempt=(
            _is_nondimensional_entry(spec, driver_context)
            or bool(spec.metadata.get("generic_case"))
        ),
        driver_context=driver_context,
    )
    simulation_audit, generation_diagnostics, readiness_score = _build_simulation_audit(
        spec=spec,
        driver_context=driver_context,
        workflow_dag=workflow_dag,
        artifacts=artifacts,
        validation_diagnostics=validation_diagnostics,
        workflow_diagnostics=workflow_diagnostics,
        artifact_diagnostics=artifact_diagnostics,
        environment_diagnostics=env_diagnostics,
        mesh_geometry_diagnostics=mesh_diagnostics,
        required_case_files=tuple(
            rule.path
            for rule in driver_context.capabilities.cxx_mapping.profile().case_files
            if rule.required == "always"
        ),
    )
    plan_diagnostics = (
        generation_diagnostics
        + validation_diagnostics
        + workflow_diagnostics
        + catalog_diagnostics
        + artifact_diagnostics
        + mesh_diagnostics
    )
    # The plugin owns solver capabilities and model-specific field exposure.
    # Keep the established payload shape for cardiacFoam compatibility while
    # attaching the immutable identity that supplied it.
    raw_capability_manifest = dict(driver_context.capabilities.manifest.manifest())
    raw_capability_manifest["plugin_identity"] = driver_context.identity.to_json()
    capability_manifest = _jsonable(raw_capability_manifest)
    function_object_diagnostics = function_object_field_diagnostics(
        spec.case_root,
        samplable=raw_capability_manifest.get("samplable_fields", {}),
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
        plugin=driver_context.identity.to_json(),
    )
