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
#     manufactured_eikonal_ecg
#
# Description
#     Defines configuration template for manufactured eikonal ECG scenarios.
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

from openfoam_driver.core.defaults import manufactured_eikonal_ecg as defaults
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec
from openfoam_driver.core.runtime.mutators import update_foam_entry
from openfoam_driver.core.runtime.parallel_execution import solve_steps
from openfoam_driver.postprocessing.driver import PostprocessTask, run_postprocess_tasks
from openfoam_driver.specs.common import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
    replace_single_block_mesh_resolution,
    resolve_run_script_path,
    resolve_spec_paths,
)
from openfoam_driver.specs.tet_mesh_provisioning import render_tet_geo


_GRAD_SCHEME_TOKENS: dict[str, str] = {
    "gauss_linear": "Gauss linear",
    "least_squares": "leastSquares",
}

_TET_NUMERICS_PROFILES: dict[str, tuple[str, ...]] = {
    "eikonal_tet": ("fvSolution",),
}


def _build_cases(
    number_cells: Sequence[int],
    dimensions: Sequence[str],
    solver_types: Sequence[str],
) -> list[CaseConfig]:
    cases: list[CaseConfig] = []
    for dimension, solver, cells in product(dimensions, solver_types, number_cells):
        case_id = f"{dimension}_{int(cells)}_cells_{solver}"
        cases.append(
            CaseConfig(
                case_id=case_id,
                params={
                    "dimension": str(dimension),
                    "solver": str(solver),
                    "cells": int(cells),
                },
            )
        )
    return cases


def _replace_blockmesh_resolution(
    block_mesh_dict_path: Path,
    cells: int,
    dimension: str,
) -> None:
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
    # 1D is cheap enough that the original bash scripts never bothered
    # decomposing it; preserved here rather than in the (now dead) legacy
    # _run_case path, since this is the mechanism sweep-run actually executes.
    effective_run_in_parallel = run_in_parallel and dimensions_list[-1] != "1D"

    if mesh_family == "tet":
        mesh_steps = [
            {"id": "clean", "command": "Allclean", "depends_on": []},
            {
                "id": "gmsh",
                "command": "gmsh",
                "args": [
                    "-3",
                    "setup/mesh/tet/box.geo",
                    "-o",
                    "box.msh",
                    "-format",
                    "msh2",
                ],
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
        run_in_parallel=effective_run_in_parallel,
        case_root=case_root,
    )
    return {"steps": mesh_steps + steps}


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    electro_properties_relpath: Path = Path("constant/electroProperties"),
    physics_properties_relpath: Path = Path("constant/physicsProperties"),
    electro_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    physics_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    verification_model_type: str = defaults.VERIFICATION_MODEL_TYPE,
    ecg_reference_quadrature_order: int = defaults.ECG_REFERENCE_QUADRATURE_ORDER,
    ecg_check_quadrature_orders: Sequence[int] = defaults.ECG_CHECK_QUADRATURE_ORDERS,
    ecg_electrodes_by_dimension: Mapping[str, Mapping[str, str]] = (
        defaults.ECG_ELECTRODES_BY_DIMENSION
    ),
    block_mesh_dict_template: str = defaults.BLOCK_MESH_DICT_TEMPLATE,
    mesh_family: str = "hex",
    numerics_profile: str | None = None,
    grad_scheme: str | None = None,
) -> None:
    dimension = str(case.params["dimension"])
    cells = int(case.params["cells"])

    electro_properties = case_root / electro_properties_relpath
    physics_properties = case_root / physics_properties_relpath
    block_mesh_dict = case_root / Path(block_mesh_dict_template.format(dimension=dimension))
    ecg_scope = f"{electro_properties_scope}.ecgDomains.ECG"

    try:
        electrodes = ecg_electrodes_by_dimension[dimension]
    except KeyError as exc:
        raise ValueError(f"Missing ECG electrode set for dimension '{dimension}'") from exc

    case_overrides = {
        f"{electro_properties_scope}.verificationModel.type": verification_model_type,
        f"{ecg_scope}.ecgSolver": "eikonalECG",
        f"{ecg_scope}.manufacturedEikonalECG.enabled": True,
        f"{ecg_scope}.manufacturedEikonalECG.referenceQuadratureOrder":
            int(ecg_reference_quadrature_order),
        f"{ecg_scope}.manufacturedEikonalECG.checkQuadratureOrders": "("
        + " ".join(str(int(value)) for value in ecg_check_quadrature_orders)
        + ")",
    }

    for electrode_name, electrode_position in electrodes.items():
        case_overrides[f"{ecg_scope}.electrodePositions.{electrode_name}"] = (
            electrode_position
        )

    if mesh_family == "tet":
        render_tet_geo(case_root, cells)
        for overlay_name in _TET_NUMERICS_PROFILES.get(numerics_profile or "", ()):
            overlay_source = case_root / "setup" / "mesh" / "tet" / overlay_name
            shutil.copy(overlay_source, case_root / "system" / overlay_name)
    else:
        _replace_blockmesh_resolution(block_mesh_dict, cells, dimension)

    if grad_scheme is not None:
        update_foam_entry(
            case_root / "system" / "fvSchemes",
            "default",
            _GRAD_SCHEME_TOKENS[grad_scheme],
            scope=["gradSchemes"],
        )
    apply_electro_property_overrides(electro_properties, case_overrides)
    apply_electro_property_overrides(electro_properties, electro_property_overrides)
    apply_physics_property_overrides(physics_properties, physics_property_overrides)


def _archive_output_dir(case_root: Path) -> Path:
    return case_root / "postProcessing"


def _candidate_paths(case_root: Path, source_name: str) -> tuple[Path, ...]:
    # The manufactured verifiers write via Time::globalPath(), so their
    # postProcessing/ output lands in the shared case dir under both serial
    # and parallel (./Allrun parallel) execution. processor0/postProcessing/
    # is kept as a fallback only for output from an unrebuilt/older solver
    # binary that predates that fix.
    return (
        case_root / "postProcessing" / source_name,
        case_root / "processor0" / "postProcessing" / source_name,
    )


def _stage_case_outputs(
    case_root: Path,
    case: CaseConfig,
) -> list[Path]:
    staged_outputs: list[Path] = []
    destination_dir = _archive_output_dir(case_root)
    destination_dir.mkdir(parents=True, exist_ok=True)

    for source_name in (
        "eikonalECG.dat",
        "manufacturedEikonalECG.dat",
        "manufacturedEikonalECGSummary.dat",
        "manufacturedEikonalActivationTime.dat",
    ):
        destination = destination_dir / f"{case.case_id}_{source_name}"
        for candidate in _candidate_paths(case_root, source_name):
            if not candidate.exists():
                continue
            if candidate.parent == destination_dir:
                shutil.move(str(candidate), str(destination))
            else:
                shutil.copy2(candidate, destination)
            print(f"Archived eikonal ECG output: {candidate} -> {destination}")
            staged_outputs.append(destination)
            break

    return staged_outputs


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


def _run_case(
    case_root: Path,
    setup_root: Path,
    case: CaseConfig,
    *,
    tutorials_root: Path | None = None,
    run_script_relpath: Path = defaults.RUN_SCRIPT_RELPATH,
    run_in_parallel: bool = defaults.RUN_IN_PARALLEL,
) -> None:
    del setup_root
    dimension = str(case.params["dimension"])
    run_script = resolve_run_script_path(
        tutorials_root=tutorials_root,
        run_script_relpath=run_script_relpath,
    )
    # 1D manufactured meshes only have 10-80 cells total; decomposing them
    # across the tutorial's 6-way decomposeParDict leaves some ranks with
    # 1-2 cells, which has produced a DILU-preconditioner SIGFPE (degenerate
    # local matrix) at tight PIMPLE tolerance. 6-way decomposition is fine for
    # the 2D/3D cases, so only 1D is forced to run serially.
    case_run_in_parallel = run_in_parallel and dimension != "1D"
    command = [
        "bash",
        "-l",
        str(run_script),
        "--case-dir",
        str(case_root),
        "--dimension",
        dimension,
    ]
    if case_run_in_parallel:
        command.append("--parallel")

    try:
        subprocess.run(command, check=True)
    finally:
        _archive_case_logs(case_root, case)
    _stage_case_outputs(case_root, case)


def _collect_outputs(case_root: Path, output_dir: Path) -> None:
    archived_dir = _archive_output_dir(case_root)
    archived_outputs = []
    if archived_dir.exists():
        archived_outputs = sorted(archived_dir.glob("*.dat"))

    if archived_outputs:
        if archived_dir.resolve() != output_dir.resolve():
            output_dir.mkdir(parents=True, exist_ok=True)
            for stale_output in output_dir.glob("*.dat"):
                stale_output.unlink()
            for source in archived_outputs:
                destination = output_dir / source.name
                shutil.copy2(source, destination)
                print(f"Copied output: {source.name} -> {destination}")

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
    dimensions: Sequence[str] = defaults.DIMENSIONS,
    solver_types: Sequence[str] = defaults.SOLVER_TYPES,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    electro_properties_relpath: str | Path = "constant/electroProperties",
    physics_properties_relpath: str | Path = "constant/physicsProperties",
    electro_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    physics_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    verification_model_type: str = defaults.VERIFICATION_MODEL_TYPE,
    ecg_reference_quadrature_order: int = defaults.ECG_REFERENCE_QUADRATURE_ORDER,
    ecg_check_quadrature_orders: Sequence[int] = defaults.ECG_CHECK_QUADRATURE_ORDERS,
    ecg_electrodes_by_dimension: Mapping[str, Mapping[str, str]] = (
        defaults.ECG_ELECTRODES_BY_DIMENSION
    ),
    block_mesh_dict_template: str = defaults.BLOCK_MESH_DICT_TEMPLATE,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    postprocess_script_relpath: str | Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    run_in_parallel: bool = defaults.RUN_IN_PARALLEL,
    postprocess_strict_artifacts: bool = False,
    mesh_family: str = "hex",
    numerics_profile: str | None = None,
    grad_scheme: str | None = None,
) -> TutorialSpec:
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
    if numerics_profile is not None and numerics_profile not in _TET_NUMERICS_PROFILES:
        known = ", ".join(sorted(_TET_NUMERICS_PROFILES))
        raise ValueError(f"numerics_profile must be one of: {known}; got {numerics_profile!r}")
    if grad_scheme is not None and grad_scheme not in _GRAD_SCHEME_TOKENS:
        known = ", ".join(sorted(_GRAD_SCHEME_TOKENS))
        raise ValueError(f"grad_scheme must be one of: {known}; got {grad_scheme!r}")

    cells_list = [int(item) for item in number_cells]
    solver_types_list = [str(item) for item in solver_types]

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
            number_cells=cells_list,
            dimensions=dimensions_list,
            solver_types=solver_types_list,
        ),
        apply_case=partial(
            _apply_case,
            electro_properties_scope=electro_properties_scope,
            electro_properties_relpath=Path(electro_properties_relpath),
            physics_properties_relpath=Path(physics_properties_relpath),
            electro_property_overrides=electro_property_overrides,
            physics_property_overrides=physics_property_overrides,
            verification_model_type=verification_model_type,
            ecg_reference_quadrature_order=ecg_reference_quadrature_order,
            ecg_check_quadrature_orders=ecg_check_quadrature_orders,
            ecg_electrodes_by_dimension=ecg_electrodes_by_dimension,
            block_mesh_dict_template=block_mesh_dict_template,
            mesh_family=mesh_family,
            numerics_profile=numerics_profile,
            grad_scheme=grad_scheme,
        ),
        run_case=partial(
            _run_case,
            tutorials_root=tutorials_root,
            run_script_relpath=Path(run_script_relpath),
            run_in_parallel=run_in_parallel,
        ),
        collect_outputs=_collect_outputs,
        postprocess=partial(
            _postprocess,
            tutorial_name=tutorial_name,
            postprocess_script_relpath=Path(postprocess_script_relpath),
            postprocess_function_name=postprocess_function_name,
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "notes": "Manufactured eikonal activation and ECG benchmark",
            "workflow_dag": _workflow_dag_for(
                mesh_family, dimensions_list,
                case_root=case_root, run_in_parallel=run_in_parallel,
            ),
            "dimensions": dimensions_list,
            "solver_types": solver_types_list,
            "mesh_family": mesh_family,
            "numerics_profile": numerics_profile,
            "grad_scheme": grad_scheme,
            "electro_properties_scope": electro_properties_scope,
            "block_mesh_dict_template": block_mesh_dict_template,
            "run_script_relpath": str(run_script_relpath),
            "run_in_parallel": run_in_parallel,
            "postprocess_script_relpath": str(postprocess_script_relpath),
            "postprocess_function_name": postprocess_function_name,
            "ecg_reference_quadrature_order": ecg_reference_quadrature_order,
            "ecg_check_quadrature_orders":
                [int(value) for value in ecg_check_quadrature_orders],
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
