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
#     manufactured_electromechanics_bc
#
# Description
#     Defines configuration template for boundary-constrained manufactured
#     electromechanics scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path

from ...core.defaults import manufactured_electromechanics_bc as defaults
from ...core.runtime.models import CaseConfig, TutorialSpec
from ...postprocessing.driver import PostprocessTask, run_postprocess_tasks
from ..common import (
    apply_electro_property_overrides,
    apply_entry_overrides,
    apply_physics_property_overrides,
    resolve_run_script_path,
    resolve_spec_paths,
    set_delta_t,
)
from .manufactured_fda import (
    _archive_case_logs,
    _build_cases,
    _collect_outputs,
    _replace_blockmesh_resolution,
)


def _case_output_filename(case: CaseConfig) -> str:
    dimension = str(case.params["dimension"])
    cells = int(case.params["cells"])
    solver = str(case.params["solver"])
    return f"{dimension}_{cells}_cells_{solver}.dat"


def _archive_output_dir(case_root: Path) -> Path:
    return case_root / "postProcessing"


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    control_dict_relpath: Path = defaults.CONTROL_DICT_RELPATH,
    electro_properties_relpath: Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    electromechanical_properties_relpath: Path = defaults.ELECTROMECHANICAL_PROPERTIES_RELPATH,
    physics_properties_relpath: Path = defaults.PHYSICS_PROPERTIES_RELPATH,
    block_mesh_dict_template: str = defaults.BLOCK_MESH_DICT_TEMPLATE,
    electro_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    electromechanical_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    physics_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    verification_model_type: str = defaults.ELECTROMECHANICAL_VERIFICATION_MODEL_TYPE,
) -> None:
    dimension = str(case.params["dimension"])
    solver = str(case.params["solver"])
    cells = int(case.params["cells"])
    dt_value = float(case.params["dt"])

    control_dict = case_root / control_dict_relpath
    electro_properties = case_root / electro_properties_relpath
    electromechanical_properties = case_root / electromechanical_properties_relpath
    physics_properties = case_root / physics_properties_relpath
    block_mesh_dict = case_root / Path(block_mesh_dict_template.format(dimension=dimension))

    _replace_blockmesh_resolution(block_mesh_dict, cells, dimension)
    set_delta_t(control_dict, dt_value)

    apply_electro_property_overrides(
        electro_properties,
        {
            f"{electro_properties_scope}.dimension": f'"{dimension}"',
            f"{electro_properties_scope}.solutionAlgorithm": solver,
        },
    )
    apply_electro_property_overrides(electro_properties, electro_property_overrides)
    apply_entry_overrides(
        electromechanical_properties,
        {
            (
                "sequentialElectroMechanicalCoeffs."
                "electromechanicalVerificationModel.type"
            ): verification_model_type,
        },
    )
    apply_entry_overrides(
        electromechanical_properties,
        electromechanical_property_overrides,
    )
    apply_physics_property_overrides(physics_properties, physics_property_overrides)


def _stage_case_output(
    case_root: Path,
    case: CaseConfig,
    *,
    run_in_parallel: bool = False,
) -> Path:
    filename = _case_output_filename(case)
    destination_dir = _archive_output_dir(case_root)
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / filename

    source_name = "manufacturedElectromechanicsSummary.dat"
    if run_in_parallel:
        candidates = (
            case_root / "processor0" / "postProcessing" / source_name,
            case_root / "postProcessing" / source_name,
        )
    else:
        candidates = (
            case_root / "postProcessing" / source_name,
            case_root / "processor0" / "postProcessing" / source_name,
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
        f"Manufactured electromechanics output '{source_name}' not found after run. "
        f"Checked: {checked}"
    )


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

    subprocess.run(command, check=True)
    _archive_case_logs(case_root, case)
    _stage_case_output(case_root, case, run_in_parallel=run_in_parallel)


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
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    control_dict_relpath: str | Path = defaults.CONTROL_DICT_RELPATH,
    electro_properties_relpath: str | Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    electromechanical_properties_relpath: str | Path = defaults.ELECTROMECHANICAL_PROPERTIES_RELPATH,
    physics_properties_relpath: str | Path = defaults.PHYSICS_PROPERTIES_RELPATH,
    block_mesh_dict_template: str = defaults.BLOCK_MESH_DICT_TEMPLATE,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    postprocess_script_relpath: str | Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    electro_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    electromechanical_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    physics_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    verification_model_type: str = defaults.ELECTROMECHANICAL_VERIFICATION_MODEL_TYPE,
    run_in_parallel: bool = defaults.RUN_IN_PARALLEL,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    dimensions_list = [str(item) for item in dimensions]
    cells_list = [int(item) for item in number_cells]
    dt_values_list = [float(item) for item in dt_values]
    solver_types_list = [str(item) for item in solver_types]
    control_dict_path = Path(control_dict_relpath)
    electro_properties_path = Path(electro_properties_relpath)
    electromechanical_properties_path = Path(electromechanical_properties_relpath)
    physics_properties_path = Path(physics_properties_relpath)
    run_script_path = Path(run_script_relpath)
    postprocess_script_path = Path(postprocess_script_relpath)

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
            electromechanical_properties_relpath=electromechanical_properties_path,
            physics_properties_relpath=physics_properties_path,
            block_mesh_dict_template=block_mesh_dict_template,
            electro_property_overrides=electro_property_overrides,
            electromechanical_property_overrides=electromechanical_property_overrides,
            physics_property_overrides=physics_property_overrides,
            verification_model_type=verification_model_type,
        ),
        run_case=partial(
            _run_case,
            tutorials_root=tutorials_root,
            run_script_relpath=run_script_path,
            run_in_parallel=run_in_parallel,
        ),
        collect_outputs=_collect_outputs,
        postprocess=partial(
            _postprocess,
            tutorial_name=tutorial_name,
            postprocess_script_relpath=postprocess_script_path,
            postprocess_function_name=postprocess_function_name,
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "notes": (
                "Boundary-constrained manufactured electromechanics benchmark. "
                "Vm is a rigorous MMS target; D and Ta are BC-constrained diagnostics."
            ),
            "workflow_dag": {
                "steps": [
                    {"id": "mesh", "command": "blockMesh", "depends_on": []},
                    {"id": "solve", "command": "cardiacFoam", "depends_on": ["mesh"]},
                ]
            },
            "dimensions": dimensions_list,
            "solver_types": solver_types_list,
            "piecewise_sweep": piecewise_sweep,
            "control_dict_relpath": str(control_dict_path),
            "electro_properties_relpath": str(electro_properties_path),
            "electromechanical_properties_relpath": str(electromechanical_properties_path),
            "physics_properties_relpath": str(physics_properties_path),
            "electro_properties_scope": electro_properties_scope,
            "block_mesh_dict_template": block_mesh_dict_template,
            "run_script_relpath": str(run_script_path),
            "run_in_parallel": run_in_parallel,
            "postprocess_script_relpath": str(postprocess_script_path),
            "postprocess_function_name": postprocess_function_name,
            "has_electro_property_overrides": bool(electro_property_overrides),
            "has_electromechanical_property_overrides": bool(electromechanical_property_overrides),
            "has_physics_property_overrides": bool(physics_property_overrides),
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
