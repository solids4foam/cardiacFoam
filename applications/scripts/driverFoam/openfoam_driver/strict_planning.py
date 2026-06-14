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
    environment_diagnostics: tuple[StrictDiagnostic, ...] = ()
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
            "environment_diagnostics": [asdict(d) for d in self.environment_diagnostics],
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

_MPI_LAUNCHERS = frozenset({"mpirun", "mpiexec", "orterun"})
_INTERPRETER_SKIP = frozenset({"python", "python3"})
# Flags that consume their following token (MPI launcher context only).
_MPI_VALUE_FLAGS = frozenset({"-np", "-n", "--np"})


def _unwrap_mpi_program(args: tuple[str, ...]) -> str | None:
    """Return the wrapped program from an MPI launcher's args, or None.

    Skips value-taking launcher flags (``-np 4`` / ``-n 4`` / ``--np 4``) and
    bare flags (``--oversubscribe``), returning the first program token.

    Limitation: only process-count flags are decoded. Other value-taking
    placement flags (``--host``, ``--hostfile``, ...) are not modelled, so the
    token following them would be misidentified as the program. This is
    acceptable for the cardiacFoam workflows we generate, which use plain
    ``mpirun -np N <solver> -parallel``.
    """
    index = 0
    while index < len(args):
        token = args[index]
        if token in _MPI_VALUE_FLAGS:
            index += 2  # skip the flag and its value
            continue
        if token.startswith("-"):
            index += 1  # bare flag
            continue
        return token
    return None


@dataclass(frozen=True)
class _ExecutableRequirements:
    executables: tuple[str, ...]
    is_parallel: bool
    mpi_launcher_in_dag: bool


def _required_executables(workflow_dag: dict[str, Any] | None) -> _ExecutableRequirements:
    """Derive the executables a plan will invoke from its workflow DAG.

    The authoritative source is each step's ``command`` (not ``launch["command"]``,
    which is only the ``python -m openfoam_driver`` re-invocation). MPI launcher
    steps contribute both the launcher and the wrapped program. Parallelism is
    inferred from an MPI launcher command, a ``-parallel`` arg, or a
    ``decomposePar`` step.
    """
    executables: list[str] = []
    is_parallel = False
    mpi_launcher_in_dag = False

    def _add(name: str) -> None:
        if name and name not in _INTERPRETER_SKIP and name not in executables:
            executables.append(name)

    for step in (workflow_dag or {}).get("steps", ()):
        command = str(step.get("command", "")).strip()
        args = tuple(str(arg) for arg in step.get("args", ()))
        if not command:
            continue
        if command in _MPI_LAUNCHERS:
            is_parallel = True
            mpi_launcher_in_dag = True
            _add(command)
            wrapped = _unwrap_mpi_program(args)
            if wrapped is not None:
                _add(wrapped)
            continue
        if command == "decomposePar" or "-parallel" in args:
            is_parallel = True
        _add(command)

    return _ExecutableRequirements(
        executables=tuple(executables),
        is_parallel=is_parallel,
        mpi_launcher_in_dag=mpi_launcher_in_dag,
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


def _environment_diagnostics(
    workflow_dag: dict[str, Any] | None,
) -> tuple[StrictDiagnostic, ...]:
    """Preflight the runtime environment against the plan's actual commands."""
    if "SKIP_ENV_DIAGNOSTICS" in os.environ:
        return ()
    diagnostics: list[StrictDiagnostic] = []

    # OpenFOAM environment.
    if "WM_PROJECT_DIR" not in os.environ:
        diagnostics.append(_diagnostic(
            "error",
            "missing_openfoam_env",
            "WM_PROJECT_DIR is not set. OpenFOAM environment not sourced.",
            source="environment",
        ))
    else:
        # Only meaningful when the base env IS sourced; otherwise the error above
        # already covers a completely unsourced environment.
        for var in ("WM_PROJECT_VERSION", "FOAM_USER_LIBBIN"):
            if var not in os.environ:
                diagnostics.append(_diagnostic(
                    "warning",
                    "partial_openfoam_env",
                    f"{var} is not set. OpenFOAM environment may be partially sourced.",
                    source="environment",
                    field=var,
                ))

    # Command-aware executable resolution.
    requirements = _required_executables(workflow_dag)
    for executable in requirements.executables:
        if not shutil.which(executable):
            diagnostics.append(_diagnostic(
                "error",
                "missing_executable",
                f"{executable} not found on PATH.",
                source="environment",
                field=executable,
            ))

    # MPI launcher when parallel but no launcher command is explicit in the DAG.
    if (
        requirements.is_parallel
        and not requirements.mpi_launcher_in_dag
        and not (shutil.which("mpirun") or shutil.which("mpiexec"))
    ):
        diagnostics.append(_diagnostic(
            "error",
            "missing_mpi",
            "Plan is parallel but no MPI launcher (mpirun/mpiexec) found on PATH.",
            source="environment",
            field="mpirun",
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
    env_diagnostics = _environment_diagnostics(workflow_dag)
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
        environment_diagnostics=env_diagnostics,
        launch=launch,
        workflow_dag=workflow_dag,
        workflow_state=workflow_state,
        expected_artifacts=artifacts,
        run_document=run_document,
    )
