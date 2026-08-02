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
#     artifacts
#
# Description
#     Predicts generated data artifacts based on workflow specifications.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Predict the data artifacts a tutorial run will (or did) produce.

The predictor is the single agent-facing answer to "what raw data does this
case produce?". It is consumed by the engine (to write
``artifacts_manifest.json`` alongside ``run_manifest.json``) and by agents
exploring a case ahead of a real run.

Design discipline (plan v2 section 3):

* **Compose, do not branch.** Solver-aware logic SHOULD live in existing
  catalogs; the predictor reads them rather than reimplementing branching.
  Today the predictor actively consumes
  ``ionic_model_catalog.IONIC_MODEL_CATALOG`` (state + algebraic variables),
  ``specs.common.detect_ionic_export_list`` (user-declared exports),
  ``active_tension_catalog.ACTIVE_TENSION_MODEL_CATALOG`` (AT state variables,
  fired when ``activeTensionModel`` block is present), and
  ``utility_catalog.UTILITY_CATALOG.produces`` (pre/post-solve utility outputs
  declared in ``workflow_dag`` steps).

* **Never raise on shape divergence.** Agents may call the predictor before
  ``apply_case`` has run, or against a partly-mutated case. Missing files,
  unknown solver names, and unknown ionic models all degrade to "return
  what we know" rather than throwing.
* **Static override wins.** A tutorial that knows it produces something the
  predictor cannot derive (e.g. analytic error norms for a manufactured
  solution) declares it via ``spec.metadata['expected_artifacts']``; on
  ``artifact_id`` collision the static entry replaces the derived one.

Adding a new solver means: write a ``_predict_<solver>`` handler and register
it in :data:`_SOLVER_HANDLERS`. Nothing else in this module branches on the
solver name.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable

from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
from ...utility_catalog import UTILITY_CATALOG, ProducesEntry
from .models import DataArtifact, TutorialSpec


SolverHandler = Callable[[Path, TutorialSpec, "str | None"], tuple[DataArtifact, ...]]


def _exported_ionic_variables(
    case_root: Path,
    ionic_model: str | None,
) -> tuple[str, ...]:
    """Return the ionic variables that will actually appear on disk.

    Prefers the declared ``outputVariables.ionic.export`` list
    (this is what the C++ writes); falls back to the catalog's
    ``recommended_exports`` when no declaration is present. Returns ``()``
    only when both the file-side declaration and the catalog entry are
    missing.
    """
    from ...specs.common import detect_ionic_export_list

    properties = case_root / "constant" / "electroProperties"
    if properties.exists():
        declared = detect_ionic_export_list(properties)
        if declared is not None:
            return declared
    if ionic_model is None:
        return ()
    entry = IONIC_MODEL_CATALOG.get(ionic_model)
    if entry is None:
        return ()
    return entry.recommended_exports


def _time_indexed_field_artifact(
    *,
    solver: str,
    field_name: str,
    ionic_model: str | None,
    description: str,
    optional: bool = False,
) -> DataArtifact:
    """Builder for one-file-per-time-dir artifacts (the OpenFOAM AUTO_WRITE
    convention). `solver` is the lowercase solver tag used to prefix the
    artifact_id; `field_name` is the on-disk filename inside `<time>/`."""
    artifact_id = f"{solver}_{field_name.lower()}_series"
    produced_by = {
        "monodomain": "monodomainSolver",
        "bidomain": "bidomainSolver",
        "eikonal": "eikonalSolver",
        "single_cell": "singleCellSolver",
    }[solver]
    return DataArtifact(
        artifact_id=artifact_id,
        path_pattern=f"{{time}}/{field_name}",
        format="openfoam_time_dirs",
        variables=(field_name,),
        description=description,
        produced_by=produced_by,
        time_indexed=True,
        optional=optional,
    )


def _predict_single_cell(
    case_root: Path, spec: TutorialSpec, ionic_model: str | None
) -> tuple[DataArtifact, ...]:
    """Two outputs:

    1. ``postProcessing/<ionicModel>_<tissue>_<protocolSuffix>.txt`` — the
       OFstream-written time-series trace.
    2. ``<time>/Vm`` — optional AUTO_WRITE on the 1-cell mesh.
    """
    if ionic_model is None:
        return ()
    return (
        DataArtifact(
            artifact_id="single_cell_trace",
            path_pattern="postProcessing/*.txt",
            format="csv_sweep",
            variables=_exported_ionic_variables(case_root, ionic_model),
            description=(
                f"Per-case time series produced by singleCellSolver "
                f"(ionicModel={ionic_model})"
            ),
            produced_by="singleCellSolver",
            time_indexed=False,
        ),
        _time_indexed_field_artifact(
            solver="single_cell",
            field_name="Vm",
            ionic_model=ionic_model,
            description=f"Membrane voltage Vm on 1-cell mesh (singleCellSolver, ionicModel={ionic_model})",
            optional=True,
        ),
    )


def _predict_monodomain(
    case_root: Path, spec: TutorialSpec, ionic_model: str | None
) -> tuple[DataArtifact, ...]:
    """Emit one artifact per ``<time>/<field>`` written by the monodomain
    solver: always ``<time>/Vm``, plus one per declared export token.
    """
    if ionic_model is None:
        return ()
    artifacts: list[DataArtifact] = []
    artifacts.append(_time_indexed_field_artifact(
        solver="monodomain",
        field_name="Vm",
        ionic_model=ionic_model,
        description=f"Membrane voltage Vm (monodomainSolver, ionicModel={ionic_model})",
    ))
    for var in _exported_ionic_variables(case_root, ionic_model):
        artifacts.append(_time_indexed_field_artifact(
            solver="monodomain",
            field_name=var,
            ionic_model=ionic_model,
            description=f"Ionic export {var} (monodomainSolver)",
        ))
    return tuple(artifacts)


def _predict_bidomain(
    case_root: Path, spec: TutorialSpec, ionic_model: str | None
) -> tuple[DataArtifact, ...]:
    """Emit one artifact per ``<time>/<field>`` written by the bidomain
    solver: Vm + phiE + phiI plus per-export ionic vars."""
    if ionic_model is None:
        return ()
    artifacts: list[DataArtifact] = []
    for field_name, description in (
        ("Vm", f"Membrane voltage Vm (bidomainSolver, ionicModel={ionic_model})"),
        ("phiE", "Extracellular potential phiE (bidomainSolver)"),
        ("phiI", "Intracellular potential phiI (bidomainSolver)"),
    ):
        artifacts.append(_time_indexed_field_artifact(
            solver="bidomain",
            field_name=field_name,
            ionic_model=ionic_model,
            description=description,
        ))
    for var in _exported_ionic_variables(case_root, ionic_model):
        artifacts.append(_time_indexed_field_artifact(
            solver="bidomain",
            field_name=var,
            ionic_model=ionic_model,
            description=f"Ionic export {var} (bidomainSolver)",
        ))
    return tuple(artifacts)


def _predict_eikonal(
    case_root: Path, spec: TutorialSpec, ionic_model: str | None
) -> tuple[DataArtifact, ...]:
    """Emit one artifact per ``<time>/<field>`` written by eikonalSolver:
    psi (activation time) + Vm (recovered membrane voltage). No ionic
    exports — eikonal does not integrate cell models."""
    return (
        _time_indexed_field_artifact(
            solver="eikonal",
            field_name="activationTime",
            ionic_model=None,
            description="Activation time (eikonalSolver)",
        ),
    )


def _predict_ecg(case_root: Path) -> tuple[DataArtifact, ...]:
    """When ``ecgDomains`` block is declared, predict the ECG time-series
    file written by the chosen ``ecgSolver``."""
    from ...specs.common import electro_properties_has_block

    properties = case_root / "constant" / "electroProperties"
    if not properties.exists():
        return ()
    if not electro_properties_has_block(properties, "ecgDomains"):
        return ()

    text = properties.read_text()
    artifacts: list[DataArtifact] = []
    if "pseudoECG" in text:
        artifacts.append(DataArtifact(
            artifact_id="ecg_pseudo_ecg",
            path_pattern="postProcessing/pseudoECG.dat",
            format="csv_probe",
            description="Pseudo-ECG time series at the declared electrodes",
            produced_by="pseudoECG",
            time_indexed=False,
        ))
    if "torsoECG" in text:
        artifacts.append(DataArtifact(
            artifact_id="ecg_torso_ecg",
            path_pattern="postProcessing/torsoECG.dat",
            format="csv_probe",
            description="Torso-ECG time series at the declared electrodes",
            produced_by="torsoECG",
            time_indexed=False,
        ))
    return tuple(artifacts)


def _predict_purkinje(case_root: Path) -> tuple[DataArtifact, ...]:
    """When ``conductionNetworkDomains`` block is declared, predict the two
    Purkinje outputs: the per-timestep ``.dat`` time series and the
    per-timestep VTK series (6-digit zero-padded timeIndex; globbed with ``*``).
    """
    from ...specs.common import electro_properties_has_block

    properties = case_root / "constant" / "electroProperties"
    if not properties.exists():
        return ()
    if not electro_properties_has_block(properties, "conductionNetworkDomains"):
        return ()

    return (
        DataArtifact(
            artifact_id="purkinje_network_time_series",
            path_pattern="postProcessing/purkinjeNetwork.dat",
            format="csv_probe",
            description=(
                "Purkinje network time-series — node Vm, activation times, "
                "PVJ coupling currents (one row per writeInterval)"
            ),
            produced_by="conductionSystemDomain",
            time_indexed=False,
        ),
        DataArtifact(
            artifact_id="purkinje_network_vtk_series",
            path_pattern="postProcessing/purkinjeNetworkVTK/purkinjeNetwork_*.vtk",
            format="vtk_sequence",
            description=(
                "Per-timestep Purkinje network VTK — one file per write step "
                "named purkinjeNetwork_<6-digit-timeIndex>.vtk"
            ),
            produced_by="conductionSystemDomain",
            time_indexed=False,
        ),
    )


def _predict_verification(case_root: Path) -> tuple[DataArtifact, ...]:
    """When ``verificationModel.type`` is declared, predict the verifier's
    error-summary file. The exact filename encodes (dimension, cells,
    algorithm) chosen at runtime — globbed via ``*`` since those tokens
    aren't recoverable from electroProperties alone.
    """
    from ...specs.common import detect_verification_model_type

    properties = case_root / "constant" / "electroProperties"
    if not properties.exists():
        return ()
    verifier_type = detect_verification_model_type(properties)
    if verifier_type is None:
        return ()
    return (
        DataArtifact(
            artifact_id="verification_error_summary",
            path_pattern="postProcessing/manufactured*Summary*.dat" if "Eikonal" in verifier_type else "postProcessing/*_*_cells_*.dat",
            format="csv_probe",
            description=(
                f"Manufactured-solution L1/L2/Linf error norms emitted by "
                f"{verifier_type}"
            ),
            produced_by=verifier_type,
            time_indexed=False,
        ),
    )


def _predict_active_tension(case_root: Path) -> tuple[DataArtifact, ...]:
    """When an ``activeTensionModel`` block is declared, predict the
    ``<time>/Ta`` field written by the electromechanical solver.

    Variable list is driven by the declared ``outputVariables.activeTension.export``
    block when present; otherwise falls back to the catalog's
    ``recommended_exports`` for the named model.
    """
    from ...specs.common import detect_active_tension_model_name, detect_active_tension_export_list
    from openfoam_driver.plugins.cardiacfoam.active_tension_catalog import ACTIVE_TENSION_MODEL_CATALOG

    properties = case_root / "constant" / "electroProperties"
    if not properties.exists():
        return ()
    at_model = detect_active_tension_model_name(properties)
    if at_model is None:
        return ()

    declared = detect_active_tension_export_list(properties)
    if declared is not None:
        variables = declared
    else:
        entry = ACTIVE_TENSION_MODEL_CATALOG.get(at_model)
        variables = entry.recommended_exports if entry is not None else ("Ta",)

    return tuple(
        DataArtifact(
            artifact_id=f"active_tension_{var}_series",
            path_pattern=f"{{time}}/{var}",
            format="openfoam_time_dirs",
            variables=(var,),
            description=f"Active tension {var} (activeTensionModel={at_model})",
            produced_by="sequentialElectroMechanical",
            time_indexed=True,
        )
        for var in variables
    )


def _produces_entry_to_artifact(
    entry: "ProducesEntry",
    utility_name: str,
) -> DataArtifact:
    """Translate a utility manifest's ProducesEntry into a DataArtifact.

    `produced_by` defaults to the utility name when the manifest leaves
    it blank — agents need to attribute the artifact regardless.
    """
    return DataArtifact(
        artifact_id=entry.artifact_id,
        path_pattern=entry.path_pattern,
        format=entry.format,
        variables=entry.variables,
        description=entry.description,
        produced_by=entry.produced_by or utility_name,
        optional=entry.optional,
        time_indexed=entry.time_indexed,
    )


def _predict_from_workflow_utilities(spec: TutorialSpec) -> tuple[DataArtifact, ...]:
    """Walk spec.metadata['workflow_dag'].steps; for each step whose
    `command` matches a utility in UTILITY_CATALOG, emit its `produces`
    entries as DataArtifacts.

    Returns ``()`` when the spec has no workflow_dag, no steps, or no
    matching utility commands. Unknown command names (e.g. OpenFOAM
    built-ins like ``blockMesh``) are silently skipped.
    """
    dag = spec.metadata.get("workflow_dag") if spec.metadata else None
    if not dag:
        return ()
    steps = dag.get("steps", ())
    if not steps:
        return ()
    derived: list[DataArtifact] = []
    for step in steps:
        command = step.get("command")
        if not command or command not in UTILITY_CATALOG:
            continue
        manifest = UTILITY_CATALOG[command]
        for produce in manifest.produces:
            derived.append(_produces_entry_to_artifact(produce, command))
    return tuple(derived)


_SOLVER_HANDLERS: dict[str, SolverHandler] = {
    "singleCellSolver": _predict_single_cell,
    "monodomainSolver": _predict_monodomain,
    "bidomainSolver": _predict_bidomain,
    "eikonalSolver": _predict_eikonal,
}


def _merge_static_override(
    derived: tuple[DataArtifact, ...],
    static: Iterable[DataArtifact],
) -> tuple[DataArtifact, ...]:
    by_id: dict[str, DataArtifact] = {a.artifact_id: a for a in derived}
    for override in static:
        by_id[override.artifact_id] = override
    return tuple(by_id.values())


def _read_solver_and_ionic(case_root: Path) -> tuple[str, str | None] | None:
    """Return ``(myocardium_solver, ionic_model)`` or ``None`` if the
    electroProperties file is missing entirely."""
    properties = case_root / "constant" / "electroProperties"
    if not properties.exists():
        return None
    from ...specs.common import (
        detect_ionic_model_name,
        detect_myocardium_solver_name,
    )

    try:
        solver = detect_myocardium_solver_name(properties)
    except KeyError:
        return None
    try:
        ionic = detect_ionic_model_name(properties)
    except KeyError:
        ionic = None
    return solver, ionic


def predict_data_artifacts(
    case_root: Path,
    spec: TutorialSpec,
) -> tuple[DataArtifact, ...]:
    """Return the artifacts ``case_root`` will (or does) produce.

    Composes:

    * the static ``spec.metadata['expected_artifacts']`` override (if any),
    * solver-specific derivations driven by ``constant/electroProperties``,
      sourced from the ionic-model and active-tension catalogs.

    Never raises. Returns ``()`` when nothing can be derived and no static
    override is supplied.
    """
    static_override = spec.metadata.get("expected_artifacts", ()) if spec.metadata else ()
    static_tuple = tuple(static_override)

    read = _read_solver_and_ionic(case_root)
    if read is None:
        return static_tuple

    solver, ionic_model = read
    handler = _SOLVER_HANDLERS.get(solver)
    solver_derived: tuple[DataArtifact, ...] = (
        handler(case_root, spec, ionic_model) if handler is not None else ()
    )
    utility_derived = _predict_from_workflow_utilities(spec)
    derived = (
        solver_derived
        + _predict_ecg(case_root)
        + _predict_purkinje(case_root)
        + _predict_verification(case_root)
        + _predict_active_tension(case_root)
        + utility_derived
    )
    return _merge_static_override(derived, static_tuple)
