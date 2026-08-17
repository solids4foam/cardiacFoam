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
#     manufactured_fda
#
# Description
#     Defines configuration template for manufactured FDA scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Mapping, Sequence
from functools import partial
from itertools import product
from pathlib import Path

from openfoam_driver.plugins.cardiacfoam.defaults import manufactured_fda as defaults
from openfoam_driver.postprocessing.driver import PostprocessTask, run_postprocess_tasks
from openfoam_driver.plugins.cardiacfoam.overrides import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
    remove_electro_property_dict,
)
from openfoam_driver.specs.common import (
    replace_single_block_mesh_resolution,
    resolve_run_script_path,
    resolve_spec_paths,
    set_delta_t,
)
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec
from openfoam_driver.core.runtime.mutators import update_foam_entry
from openfoam_driver.core.runtime.parallel_execution import solve_steps
from openfoam_driver.specs.tet_mesh_provisioning import render_tet_geo

# Stable driver-side names mapped to the literal OpenFOAM tokens -- "GaussLinear"
# (no space) is a CSV/shell label from the original bash sweep scripts, not a
# valid gradSchemes value; the real tokens are "Gauss linear" (two words) and
# "leastSquares". Never forward a free-form string into the dict file.
_GRAD_SCHEME_TOKENS: dict[str, str] = {
    "gauss_linear": "Gauss linear",
    "least_squares": "leastSquares",
}

# Which system/ overlay files mesh_family="tet" installs from setup/mesh/tet/
# -- not every tet case has the same set. Confirmed: monodomainPseudoECG's
# tet variant additionally installs an fvSolution overlay that bidomain's
# does not have. Explicit per case, never inferred from what happens to
# exist on disk.
_NUMERICS_PROFILES: dict[str, tuple[str, ...]] = {
    "bidomain_tet": ("fvSchemes",),
    "monodomain_tet": ("fvSchemes", "fvSolution"),
}


def _normalize_convergence_axis(convergence_axis: str) -> str:
    axis = str(convergence_axis).strip().lower()
    if axis not in {"spatial", "temporal"}:
        raise ValueError(
            f"Unsupported convergence_axis '{convergence_axis}'. "
            "Expected 'spatial' or 'temporal'."
        )
    return axis


def _format_dt(dt_value: float) -> str:
    token = f"{dt_value:.12g}".replace(".", "p")
    return token.replace("-", "m")


def _case_output_filename(case: CaseConfig, *, convergence_axis: str = "spatial") -> str:
    axis = _normalize_convergence_axis(convergence_axis)
    dimension = str(case.params["dimension"])
    cells = int(case.params["cells"])
    solver = str(case.params["solver"])
    if axis == "temporal":
        dt_value = float(case.params["dt"])
        return f"{dimension}_{cells}_cells_{solver}_DT{_format_dt(dt_value)}.dat"
    return f"{dimension}_{cells}_cells_{solver}.dat"


def _sanitize_archive_tag(value: str) -> str:
    return "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in value)


def _archive_output_dir(case_root: Path, *, archive_tag: str = "default") -> Path:
    return case_root / f"driverPostProcessingArchive_{_sanitize_archive_tag(archive_tag)}"


def _build_cases(
    dt_values: Sequence[float],
    number_cells: Sequence[int],
    dimensions: Sequence[str],
    solver_types: Sequence[str],
    piecewise_sweep: bool,
) -> list[CaseConfig]:
    if piecewise_sweep:
        dt_cells_pairs = list(zip(dt_values, number_cells))
    else:
        dt_cells_pairs = list(product(dt_values, number_cells))

    cases: list[CaseConfig] = []
    for dimension, solver in product(dimensions, solver_types):
        for dt_value, cells in dt_cells_pairs:
            case_id = f"{dimension}_{cells}_cells_{solver}_DT{_format_dt(dt_value)}"
            cases.append(
                CaseConfig(
                    case_id=case_id,
                    params={
                        "dimension": dimension,
                        "solver": solver,
                        "cells": int(cells),
                        "dt": float(dt_value),
                    },
                )
            )
    return cases


def _replace_blockmesh_resolution(block_mesh_dict_path: Path, cells: int, dimension: str) -> None:
    replace_single_block_mesh_resolution(
        block_mesh_dict_path, cells, dimension,
        resolution_by_dimension=defaults.BLOCK_MESH_RESOLUTION_BY_DIMENSION,
    )


def _workflow_dag_for(
    mesh_family: str,
    dimensions_list: list[str],
    *,
    case_root: Path,
    run_in_parallel: bool = False,
) -> dict[str, object]:
    """Build this spec's workflow_dag, branching on mesh_family.

    mesh_family="hex" (default) is unchanged from before this kwarg existed.
    mesh_family="tet" must *not* run blockMesh -- it would silently overwrite
    the tet mesh with a Cartesian block. gmsh/gmshToFoam/checkMesh run as
    real workflow steps (executed only when run --strict/sweep-run actually
    runs), never inside apply_case/materialization.

    run_in_parallel=True wraps the solve step with decomposePar/reconstructPar
    via solve_steps() -- identical mechanism regardless of mesh_family, so
    the mesh-generation steps above are unaffected either way.
    """
    if mesh_family == "tet":
        mesh_steps = [
            {"id": "clean", "command": "Allclean", "depends_on": []},
            {
                "id": "gmsh",
                "command": "gmsh",
                # Bare "gmsh" with no args launches its GUI and hangs
                # forever instead of meshing anything (found via a real,
                # non-mocked sweep-run) -- these are the same args the
                # original bash scripts pass.
                "args": ["-3", "setup/mesh/tet/box.geo", "-o", "box.msh", "-format", "msh2"],
                "depends_on": ["clean"],
            },
            {
                "id": "gmshToFoam",
                "command": "gmshToFoam",
                "args": ["box.msh"],
                "depends_on": ["gmsh"],
            },
            {"id": "checkMesh", "command": "checkMesh", "depends_on": ["gmshToFoam"]},
        ]
        solve_depends_on = ["checkMesh"]
    else:
        mesh_steps = [
            {
                "id": "mesh",
                "command": "blockMesh",
                "args": ["-dict", f"system/blockMeshDict.{dimensions_list[-1]}"],
                "depends_on": [],
            },
        ]
        solve_depends_on = ["mesh"]

    steps, _final_id = solve_steps(
        solve_id="solve",
        solve_command="cardiacFoam",
        depends_on=solve_depends_on,
        run_in_parallel=run_in_parallel,
        case_root=case_root,
    )
    return {"steps": mesh_steps + steps}


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    control_dict_relpath: Path = defaults.CONTROL_DICT_RELPATH,
    electro_properties_relpath: Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    physics_properties_relpath: Path = Path("constant/physicsProperties"),
    electro_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    physics_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    verification_model_type: str = defaults.VERIFICATION_MODEL_TYPE,
    conductivity: str | None = None,
    ecg_enabled: bool = defaults.ECG_ENABLED,
    ecg_reference_quadrature_order: int = defaults.ECG_REFERENCE_QUADRATURE_ORDER,
    ecg_check_quadrature_orders: Sequence[int] = defaults.ECG_CHECK_QUADRATURE_ORDERS,
    ecg_electrodes_by_dimension: Mapping[str, Mapping[str, str]] = defaults.ECG_ELECTRODES_BY_DIMENSION,
    block_mesh_dict_template: str = defaults.BLOCK_MESH_DICT_TEMPLATE,
    mesh_family: str = "hex",
    tet_geo_template_relpath: Path = Path("setup/mesh/tet/box.geo.template"),
    numerics_profile: str | None = None,
    grad_scheme: str | None = None,
    phi_tolerance: float | None = None,
    end_time: float | None = None,
    fv_scheme_overrides: Sequence[Mapping[str, object]] | None = None,
    fv_solution_overrides: Sequence[Mapping[str, object]] | None = None,
) -> None:
    dimension = str(case.params["dimension"])
    solver = str(case.params["solver"])
    cells = int(case.params["cells"])
    dt_value = float(case.params["dt"])

    control_dict = case_root / control_dict_relpath
    electro_properties = case_root / electro_properties_relpath
    physics_properties = case_root / physics_properties_relpath
    block_mesh_dict = case_root / Path(block_mesh_dict_template.format(dimension=dimension))
    case_overrides = {
        f"{electro_properties_scope}.dimension": f'"{dimension}"',
        f"{electro_properties_scope}.solutionAlgorithm": solver,
        f"{electro_properties_scope}.verificationModel.type": verification_model_type,
    }
    if conductivity is not None:
        case_overrides[f"{electro_properties_scope}.conductivity"] = conductivity

    ecg_scope = f"{electro_properties_scope}.ecgDomains.ECG"
    if ecg_enabled:
        try:
            electrodes = ecg_electrodes_by_dimension[dimension]
        except KeyError as exc:
            raise ValueError(f"Missing ECG electrode set for dimension '{dimension}'") from exc

        case_overrides.update(
            {
                f"{ecg_scope}.ecgSolver": "pseudoECG",
                f"{ecg_scope}.manufactured.enabled": True,
                f"{ecg_scope}.manufactured.dimension": f'"{dimension}"',
                f"{ecg_scope}.manufactured.referenceQuadratureOrder": int(
                    ecg_reference_quadrature_order
                ),
                f"{ecg_scope}.manufactured.checkQuadratureOrders": "("
                + " ".join(str(int(value)) for value in ecg_check_quadrature_orders)
                + ")",
            }
        )
        for electrode_name, electrode_position in electrodes.items():
            case_overrides[
                f"{ecg_scope}.electrodePositions.{electrode_name}"
            ] = electrode_position

    if mesh_family == "tet":
        # Render-only: substitutes __LC__ and writes overlay files. gmsh/
        # gmshToFoam/checkMesh are workflow_dag steps, not run here -- see
        # _workflow_dag_for's docstring.
        render_tet_geo(case_root, cells, template_relpath=tet_geo_template_relpath)
        for overlay_name in _NUMERICS_PROFILES.get(numerics_profile or "", ()):
            overlay_source = case_root / "setup" / "mesh" / "tet" / overlay_name
            shutil.copy(overlay_source, case_root / "system" / overlay_name)
    else:
        _replace_blockmesh_resolution(block_mesh_dict, cells, dimension)

    set_delta_t(control_dict, dt_value)
    if end_time is not None:
        update_foam_entry(control_dict, "endTime", end_time)
    if grad_scheme is not None:
        update_foam_entry(
            case_root / "system" / "fvSchemes",
            "default",
            _GRAD_SCHEME_TOKENS[grad_scheme],
            scope=["gradSchemes"],
        )
    if phi_tolerance is not None:
        update_foam_entry(
            case_root / "system" / "fvSolution",
            "tolerance",
            phi_tolerance,
            scope=["solvers", '"phiE|phiEFinal|phiI|phiIFinal"'],
        )
    for entry in fv_scheme_overrides or ():
        update_foam_entry(
            case_root / "system" / "fvSchemes", entry["key"], entry["value"],
            scope=entry.get("scope"),
        )
    for entry in fv_solution_overrides or ():
        update_foam_entry(
            case_root / "system" / "fvSolution", entry["key"], entry["value"],
            scope=entry.get("scope"),
        )

    apply_electro_property_overrides(electro_properties, case_overrides)
    if not ecg_enabled:
        # Entry sweeps reuse one physical tutorial directory. Remove an ECG
        # block left by a preceding ECG-enabled case: the system builder does
        # not route an explicit ecgSolver=none provider.
        remove_electro_property_dict(
            electro_properties,
            "ecgDomains",
            scope=electro_properties_scope,
            missing_ok=True,
        )
    apply_electro_property_overrides(electro_properties, electro_property_overrides)
    apply_physics_property_overrides(physics_properties, physics_property_overrides)


def _run_case(
    case_root: Path,
    setup_root: Path,
    case: CaseConfig,
    *,
    tutorials_root: Path | None = None,
    run_script_relpath: Path = defaults.RUN_SCRIPT_RELPATH,
    run_in_parallel: bool = defaults.RUN_IN_PARALLEL,
    convergence_axis: str = "spatial",
    archive_tag: str = "default",
    verification_model_type: str = defaults.VERIFICATION_MODEL_TYPE,
) -> None:
    del setup_root
    dimension = str(case.params["dimension"])
    run_script = resolve_run_script_path(
        tutorials_root=tutorials_root,
        run_script_relpath=run_script_relpath,
    )
    command = [
        "bash",
        "-l",
        str(run_script),
        "--case-dir",
        str(case_root),
        "--dimension",
        dimension,
    ]
    if run_in_parallel:
        command.append("--parallel")

    try:
        subprocess.run(
            command,
            check=True,
        )
    finally:
        _archive_case_logs(case_root, case)
    _stage_case_output(
        case_root,
        case,
        convergence_axis=convergence_axis,
        archive_tag=archive_tag,
        verification_model_type=verification_model_type,
    )
    _stage_case_ecg_outputs(
        case_root,
        case,
        archive_tag=archive_tag,
    )


def _archive_case_logs(case_root: Path, case: CaseConfig) -> Path | None:
    log_files = sorted(path for path in case_root.glob("log.*") if path.is_file())
    if not log_files:
        return None

    destination_root = case_root / "logs" / case.case_id
    if destination_root.exists():
        shutil.rmtree(destination_root)
    destination_root.mkdir(parents=True, exist_ok=True)

    for source in log_files:
        shutil.copy2(source, destination_root / source.name)

    print(f"Archived {len(log_files)} log file(s) for {case.case_id}: {destination_root}")
    return destination_root


def _stage_case_output(
    case_root: Path,
    case: CaseConfig,
    *,
    convergence_axis: str = "spatial",
    archive_tag: str = "default",
    verification_model_type: str = defaults.VERIFICATION_MODEL_TYPE,
) -> Path:
    # The manufactured verifiers write via Time::globalPath(), so their
    # postProcessing/ output lands in the shared case dir under both serial
    # and parallel (./Allrun parallel) execution. processor0/postProcessing/
    # is kept as a fallback only for output from an unrebuilt/older solver
    # binary that predates that fix.
    filename = _case_output_filename(case, convergence_axis=convergence_axis)
    legacy_filename = _case_output_filename(case, convergence_axis="spatial")
    destination_dir = _archive_output_dir(case_root, archive_tag=archive_tag)
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / filename
    candidates = [
        case_root / "postProcessing" / filename,
        case_root / "processor0" / "postProcessing" / filename,
    ]

    if verification_model_type == "manufacturedAnisotropicMonodomainVerifier":
        rotated_filename = (
            f"rotatedAnisotropy_3D_{int(case.params['cells'])}_cells_"
            f"{case.params['solver']}.dat"
        )
        candidates.extend(
            [
                case_root / "postProcessing" / rotated_filename,
                case_root / "processor0" / "postProcessing" / rotated_filename,
            ]
        )

    if convergence_axis != "spatial" and legacy_filename != filename:
        candidates.extend(
            [
                case_root / "postProcessing" / legacy_filename,
                case_root / "processor0" / "postProcessing" / legacy_filename,
            ]
        )

    for candidate in candidates:
        if not candidate.exists():
            continue
        if candidate == destination:
            return destination
        shutil.copy2(candidate, destination)
        print(f"Archived manufactured output: {candidate} -> {destination}")
        return destination

    checked = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(
        f"Manufactured output '{filename}' not found after run. Checked: {checked}"
    )


def _stage_case_ecg_outputs(
    case_root: Path,
    case: CaseConfig,
    *,
    archive_tag: str = "default",
) -> list[Path]:
    staged_outputs: list[Path] = []
    destination_dir = _archive_output_dir(case_root, archive_tag=archive_tag)
    destination_dir.mkdir(parents=True, exist_ok=True)

    for source_name in (
        "pseudoECG.dat",
        "manufacturedPseudoECG.dat",
        "manufacturedPseudoECGSummary.dat",
    ):
        destination = destination_dir / f"ECG_{case.case_id}_{source_name}"
        candidates = (
            case_root / "postProcessing" / source_name,
            case_root / "processor0" / "postProcessing" / source_name,
        )

        for candidate in candidates:
            if not candidate.exists():
                continue
            if candidate == destination:
                staged_outputs.append(destination)
                break
            if candidate.parent == destination_dir:
                shutil.move(str(candidate), str(destination))
                print(f"Archived ECG output: {candidate} -> {destination} (moved)")
            else:
                shutil.copy2(candidate, destination)
                print(f"Archived ECG output: {candidate} -> {destination}")
            staged_outputs.append(destination)
            break

    return staged_outputs


def _collect_outputs(case_root: Path, output_dir: Path, *, archive_tag: str = "default") -> None:
    archived_dir = _archive_output_dir(case_root, archive_tag=archive_tag)
    archived_outputs = []
    if archived_dir.exists():
        for source in sorted(archived_dir.glob("*.dat")):
            name = source.name
            if name.startswith("ECG_") or (
                name.endswith(".dat") and "_cells_" in name and name[1:2] == "D"
            ):
                archived_outputs.append(source)
    same_output_dir = archived_dir.exists() and archived_dir.resolve() == output_dir.resolve()
    if archived_outputs:
        if same_output_dir:
            print(f"Archived outputs already available in {output_dir}; preserving in place")
        else:
            for stale_output in output_dir.glob("*.dat"):
                stale_output.unlink()
            for source in archived_outputs:
                destination = output_dir / source.name
                shutil.copy2(source, destination)
                print(f"Copied output: {source.name} -> {destination}")
    else:
        print(
            "No archived .dat files found in case root; preserving existing "
            f"outputs in {output_dir}"
        )

    source_logs = case_root / "logs"
    destination_logs = output_dir / "logs"
    if destination_logs.exists():
        shutil.rmtree(destination_logs)
    if source_logs.exists():
        shutil.copytree(source_logs, destination_logs)
        print(f"Copied archived logs -> {destination_logs}")

def _postprocess(
    setup_root: Path,
    output_dir: Path,
    *,
    tutorial_name: str = defaults.TUTORIAL_NAME,
    postprocess_script_relpath: Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    strict_artifacts: bool = False,
) -> None:
    run_postprocess_tasks(
        setup_root=setup_root,
        output_dir=output_dir,
        tutorial_name=tutorial_name,
        strict_artifacts=strict_artifacts,
        tasks=[
            PostprocessTask(
                module_relpath=postprocess_script_relpath,
                function_name=postprocess_function_name,
            )
        ],
    )


def make_spec(
    *,
    tutorials_root: Path | None = None,
    tutorial_name: str = defaults.TUTORIAL_NAME,
    case_dir_name: str = defaults.CASE_DIR_NAME,
    setup_dir_name: str | None = defaults.SETUP_DIR_NAME,
    output_dir_name: str | None = None,
    number_cells: Sequence[int] = defaults.NUMBER_CELLS,
    dt_values: Sequence[float] = defaults.DT_VALUES,
    dimensions: Sequence[str] = defaults.DIMENSIONS,
    solver_types: Sequence[str] = defaults.SOLVER_TYPES,
    piecewise_sweep: bool = defaults.PIECEWISE_SWEEP,
    convergence_axis: str = "spatial",
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    control_dict_relpath: str | Path = defaults.CONTROL_DICT_RELPATH,
    electro_properties_relpath: str | Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    physics_properties_relpath: str | Path = "constant/physicsProperties",
    electro_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    physics_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    verification_model_type: str = defaults.VERIFICATION_MODEL_TYPE,
    conductivity: str | None = None,
    ecg_enabled: bool = defaults.ECG_ENABLED,
    ecg_reference_quadrature_order: int = defaults.ECG_REFERENCE_QUADRATURE_ORDER,
    ecg_check_quadrature_orders: Sequence[int] = defaults.ECG_CHECK_QUADRATURE_ORDERS,
    ecg_electrodes_by_dimension: Mapping[str, Mapping[str, str]] = defaults.ECG_ELECTRODES_BY_DIMENSION,
    block_mesh_dict_template: str = defaults.BLOCK_MESH_DICT_TEMPLATE,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    postprocess_script_relpath: str | Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    run_in_parallel: bool = defaults.RUN_IN_PARALLEL,
    postprocess_strict_artifacts: bool = False,
    mesh_family: str = "hex",
    tet_geo_template_relpath: str | Path = "setup/mesh/tet/box.geo.template",
    numerics_profile: str | None = None,
    grad_scheme: str | None = None,
    phi_tolerance: float | None = None,
    end_time: float | None = None,
    fv_scheme_overrides: Sequence[Mapping[str, object]] | None = None,
    fv_solution_overrides: Sequence[Mapping[str, object]] | None = None,
) -> TutorialSpec:
    # Entry-based sweep-run routes a zip axis's per-case resolved value
    # straight through as a scalar (e.g. dimensions="1D"), not wrapped in a
    # list. Sequence[str] params silently iterate a bare str character-by-
    # character instead of raising, so guard both str-typed sequence kwargs
    # here rather than let that corrupt the case matrix (confirmed via a
    # real, non-mocked sweep-run: dimensions="1D" produced cases for
    # dimension "1" and dimension "D").
    if isinstance(dimensions, str):
        dimensions = [dimensions]
    if isinstance(solver_types, str):
        solver_types = [solver_types]
    if isinstance(number_cells, (int, float)):
        number_cells = [number_cells]
    if isinstance(dt_values, (int, float)):
        dt_values = [dt_values]

    dimensions_list = [str(item) for item in dimensions]
    if not dimensions_list:
        raise ValueError("dimensions cannot be empty")

    if mesh_family not in {"hex", "tet"}:
        raise ValueError(f"mesh_family must be 'hex' or 'tet'; got {mesh_family!r}")
    if mesh_family == "tet" and dimensions_list != ["3D"]:
        raise ValueError(
            f"mesh_family='tet' requires dimensions=['3D']; got {dimensions_list!r} "
            "(the unit-cube tet mesh has no 1D/2D variant)"
        )
    if grad_scheme is not None and grad_scheme not in _GRAD_SCHEME_TOKENS:
        known = ", ".join(sorted(_GRAD_SCHEME_TOKENS))
        raise ValueError(f"grad_scheme must be one of: {known}; got {grad_scheme!r}")
    if phi_tolerance is not None and phi_tolerance <= 0:
        raise ValueError(f"phi_tolerance must be positive; got {phi_tolerance}")
    if numerics_profile is not None and numerics_profile not in _NUMERICS_PROFILES:
        known = ", ".join(sorted(_NUMERICS_PROFILES))
        raise ValueError(f"numerics_profile must be one of: {known}; got {numerics_profile!r}")

    cells_list = [int(item) for item in number_cells]
    dt_values_list = [float(item) for item in dt_values]
    solver_types_list = [str(item) for item in solver_types]
    control_dict_path = Path(control_dict_relpath)
    electro_properties_path = Path(electro_properties_relpath)
    physics_properties_path = Path(physics_properties_relpath)
    tet_geo_template_path = Path(tet_geo_template_relpath)
    run_script_path = Path(run_script_relpath)
    postprocess_script_path = Path(postprocess_script_relpath)
    convergence_axis_normalized = _normalize_convergence_axis(convergence_axis)
    archive_tag = output_dir_name if output_dir_name else defaults.OUTPUT_DIR_NAME

    if piecewise_sweep and len(cells_list) != len(dt_values_list):
        raise ValueError(
            "piecewise_sweep requires number_cells and dt_values to have the same length"
        )

    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        default_output_dir_name=defaults.OUTPUT_DIR_NAME,
    )

    return TutorialSpec(
        name=tutorial_name,
        case_root=case_root,
        setup_root=setup_root,
        output_dir=output_dir,
        build_cases=partial(
            _build_cases,
            dt_values=dt_values_list,
            number_cells=cells_list,
            dimensions=dimensions_list,
            solver_types=solver_types_list,
            piecewise_sweep=piecewise_sweep,
        ),
        apply_case=partial(
            _apply_case,
            electro_properties_scope=electro_properties_scope,
            control_dict_relpath=control_dict_path,
            electro_properties_relpath=electro_properties_path,
            physics_properties_relpath=physics_properties_path,
            electro_property_overrides=electro_property_overrides,
            physics_property_overrides=physics_property_overrides,
            mesh_family=mesh_family,
            tet_geo_template_relpath=tet_geo_template_path,
            numerics_profile=numerics_profile,
            grad_scheme=grad_scheme,
            phi_tolerance=phi_tolerance,
            end_time=end_time,
            fv_scheme_overrides=fv_scheme_overrides,
            fv_solution_overrides=fv_solution_overrides,
            verification_model_type=verification_model_type,
            conductivity=conductivity,
            ecg_enabled=ecg_enabled,
            ecg_reference_quadrature_order=ecg_reference_quadrature_order,
            ecg_check_quadrature_orders=ecg_check_quadrature_orders,
            ecg_electrodes_by_dimension=ecg_electrodes_by_dimension,
            block_mesh_dict_template=block_mesh_dict_template,
        ),
        run_case=partial(
            _run_case,
            tutorials_root=tutorials_root,
            run_script_relpath=run_script_path,
            run_in_parallel=run_in_parallel,
            convergence_axis=convergence_axis_normalized,
            archive_tag=archive_tag,
            verification_model_type=verification_model_type,
        ),
        collect_outputs=partial(
            _collect_outputs,
            archive_tag=archive_tag,
        ),
        postprocess=partial(
            _postprocess,
            tutorial_name=tutorial_name,
            postprocess_script_relpath=postprocess_script_path,
            postprocess_function_name=postprocess_function_name,
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "notes": "Manufactured-solution convergence benchmark",
            "workflow_dag": _workflow_dag_for(
                mesh_family, dimensions_list,
                case_root=case_root, run_in_parallel=run_in_parallel,
            ),
            "dimensions": dimensions_list,
            "solver_types": solver_types_list,
            "piecewise_sweep": piecewise_sweep,
            "convergence_axis": convergence_axis_normalized,
            "control_dict_relpath": str(control_dict_path),
            "electro_properties_relpath": str(electro_properties_path),
            "physics_properties_relpath": str(physics_properties_path),
            "electro_properties_scope": electro_properties_scope,
            "block_mesh_dict_template": block_mesh_dict_template,
            "tet_geo_template_relpath": str(tet_geo_template_path),
            "run_script_relpath": str(run_script_path),
            "run_in_parallel": run_in_parallel,
            "postprocess_script_relpath": str(postprocess_script_path),
            "postprocess_function_name": postprocess_function_name,
            "has_electro_property_overrides": bool(electro_property_overrides),
            "conductivity_override": conductivity,
            "has_physics_property_overrides": bool(physics_property_overrides),
            "ecg_enabled": ecg_enabled,
            "ecg_reference_quadrature_order": ecg_reference_quadrature_order,
            "ecg_check_quadrature_orders": [int(value) for value in ecg_check_quadrature_orders],
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
