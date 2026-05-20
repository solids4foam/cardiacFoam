"""Predict the data artifacts a tutorial run will (or did) produce.

The predictor is the single agent-facing answer to "what raw data does this
case produce?". It is consumed by the engine (to write
``artifacts_manifest.json`` alongside ``run_manifest.json``) and by agents
exploring a case ahead of a real run.

Design discipline (plan v2 section 3):

* **Compose, do not branch.** Solver-aware logic lives in the existing
  catalogs (``dict_entries.py``, ``ionic_model_catalog.py``,
  ``active_tension_catalog.py``); the predictor reads them. Reimplementing
  branching here forks the source of truth and is what the plan exists to
  prevent.
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

from ...ionic_model_catalog import IONIC_MODEL_CATALOG
from .models import DataArtifact, TutorialSpec


SolverHandler = Callable[[Path, TutorialSpec, "str | None"], tuple[DataArtifact, ...]]


# Solver-provided PDE field names sourced from the C++ side (plan §3e).
# A C++ rename means updating these constants AND the lock tests that
# assert their presence. There is intentionally no catalog backing these
# yet — see plan §3e for the rationale.
_MONODOMAIN_FIELDS: tuple[str, ...] = ("Vm",)
_BIDOMAIN_FIELDS: tuple[str, ...] = ("Vm", "phiE", "phiI")
_EIKONAL_FIELDS: tuple[str, ...] = ("psi", "Vm")


def _exported_ionic_variables(
    case_root: Path,
    ionic_model: str | None,
) -> tuple[str, ...]:
    """Return the ionic variables that will actually appear on disk.

    Plan §3d-1: prefer the declared ``outputVariables.ionic.export`` list
    (this is what the C++ writes); fall back to the catalog's
    ``recommended_exports`` when no declaration is present. Returns ``()``
    only when both the file-side declaration and the catalog entry are
    missing.

    Aliasing note: declared export tokens are returned verbatim (the user
    chose those names because that is what they want to see in the output
    file). The catalog's ``recommended_exports`` uses C++-internal state
    names, which may differ from the user-facing aliases. Resolving the
    alias table is a future follow-up tracked in plan §3d-1.
    """
    # Late import to keep specs.common off the model-load path.
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


def _predict_single_cell(
    case_root: Path, spec: TutorialSpec, ionic_model: str | None
) -> tuple[DataArtifact, ...]:
    if ionic_model is None:
        return ()
    return (
        DataArtifact(
            artifact_id="single_cell_trace",
            path_pattern="postProcessing/{case_id}.txt",
            format="csv_sweep",
            variables=_exported_ionic_variables(case_root, ionic_model),
            description=(
                f"Per-case time series produced by singleCellSolver "
                f"(ionicModel={ionic_model})"
            ),
            produced_by="singleCellSolver",
            time_indexed=False,
        ),
    )


def _predict_monodomain(
    case_root: Path, spec: TutorialSpec, ionic_model: str | None
) -> tuple[DataArtifact, ...]:
    if ionic_model is None:
        return ()
    return (
        DataArtifact(
            artifact_id="myocardium_time_series",
            path_pattern="{time}",
            format="openfoam_time_dirs",
            variables=_MONODOMAIN_FIELDS + _exported_ionic_variables(case_root, ionic_model),
            description=(
                f"OpenFOAM time directories containing Vm and ionic fields "
                f"(ionicModel={ionic_model})"
            ),
            produced_by="monodomainSolver",
            time_indexed=True,
        ),
    )


def _predict_bidomain(
    case_root: Path, spec: TutorialSpec, ionic_model: str | None
) -> tuple[DataArtifact, ...]:
    if ionic_model is None:
        return ()
    return (
        DataArtifact(
            artifact_id="myocardium_time_series",
            path_pattern="{time}",
            format="openfoam_time_dirs",
            variables=_BIDOMAIN_FIELDS + _exported_ionic_variables(case_root, ionic_model),
            description=(
                f"OpenFOAM time directories containing Vm, phiE, phiI and "
                f"ionic fields (ionicModel={ionic_model})"
            ),
            produced_by="bidomainSolver",
            time_indexed=True,
        ),
    )


def _predict_eikonal(
    case_root: Path, spec: TutorialSpec, ionic_model: str | None
) -> tuple[DataArtifact, ...]:
    """Eikonal cases do not integrate ionic cells; the ``ionic_model``
    argument is ignored (and dict_entries.py forbids it being set when
    myocardiumSolver=eikonalSolver)."""
    return (
        DataArtifact(
            artifact_id="activation_time_field",
            path_pattern="{time}",
            format="openfoam_time_dirs",
            variables=_EIKONAL_FIELDS,
            description="Activation time (psi) and recovered Vm from eikonalSolver",
            produced_by="eikonalSolver",
            time_indexed=True,
        ),
    )


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
    electroProperties file is missing entirely.

    Imported lazily and tolerantly: agents may probe cases that don't yet
    have a constant/electroProperties file. A missing ionicModel inside an
    existing file is returned as ``(solver, None)`` — eikonal cases are the
    canonical example.
    """
    properties = case_root / "constant" / "electroProperties"
    if not properties.exists():
        return None
    # Local import to avoid pulling specs.common at module load time
    # (specs.common imports subprocess and the postprocessing chain).
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
    derived: tuple[DataArtifact, ...] = (
        handler(case_root, spec, ionic_model) if handler is not None else ()
    )
    return _merge_static_override(derived, static_tuple)
