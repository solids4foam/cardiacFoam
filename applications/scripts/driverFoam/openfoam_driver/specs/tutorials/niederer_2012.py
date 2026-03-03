from __future__ import annotations

import os
import re
import subprocess
from collections.abc import Mapping, Sequence
from functools import partial
from itertools import product
from pathlib import Path

from ...core.defaults import niederer_2012 as defaults
from ...postprocessing.driver import PostprocessTask, run_postprocess_tasks
from ..common import (
    resolve_run_script_path,
    resolve_spec_paths,
    set_delta_t,
    set_ionic_model,
    set_solution_algorithm,
    set_tissue,
)
from ...core.runtime.models import CaseConfig, TutorialSpec


def _closest_key(mapping: dict[float, object], value: float) -> float:
    return min(mapping.keys(), key=lambda item: abs(item - value))


def _replace_blockmesh_resolution(
    block_mesh_dict_path: Path,
    dx: float,
    slab_size_mm: Sequence[float] = defaults.SLAB_SIZE_MM,
) -> None:
    if not block_mesh_dict_path.exists():
        raise FileNotFoundError(f"blockMeshDict not found: {block_mesh_dict_path}")

    if dx <= 0:
        raise ValueError(f"dx must be positive; got {dx}")
    if len(slab_size_mm) != 3:
        raise ValueError("slab_size_mm must have exactly 3 entries (x, y, z)")

    axis_cell_counts: list[str] = []
    for axis_length in slab_size_mm:
        raw_cells = float(axis_length) / dx
        rounded_cells = int(round(raw_cells))
        if abs(raw_cells - rounded_cells) > 1e-9:
            raise ValueError(
                f"dx={dx} does not evenly divide slab axis length {axis_length} mm"
            )
        if rounded_cells <= 0:
            raise ValueError(
                f"Computed non-positive cell count for axis length {axis_length} with dx={dx}"
            )
        axis_cell_counts.append(str(rounded_cells))

    replacement_line = (
        f"hex (0 1 2 3 4 5 6 7) ({' '.join(axis_cell_counts)}) simpleGrading (1 1 1)\n"
    )

    lines = block_mesh_dict_path.read_text().splitlines(keepends=True)
    replaced = False
    with block_mesh_dict_path.open("w") as handle:
        for line in lines:
            stripped = line.strip()
            if (
                not replaced
                and stripped.startswith("hex (0 1 2 3 4 5 6 7)")
                and not stripped.startswith("//")
            ):
                handle.write(replacement_line)
                replaced = True
            else:
                handle.write(line)

    if not replaced:
        raise KeyError(f"Did not find target hex line in {block_mesh_dict_path}")


def _update_end_time(control_dict_path: Path, dx: float, end_time_by_dx: Mapping[float, float]) -> None:
    if not control_dict_path.exists():
        raise FileNotFoundError(f"controlDict not found: {control_dict_path}")

    key = _closest_key(dict(end_time_by_dx), dx)
    end_time_value = end_time_by_dx[key]

    lines = control_dict_path.read_text().splitlines(keepends=True)
    pattern = re.compile(r"^\s*endTime\b")
    replaced = False

    with control_dict_path.open("w") as handle:
        for line in lines:
            stripped = line.strip()
            if pattern.match(line) and not stripped.startswith("//"):
                indent = line[: len(line) - len(line.lstrip())]
                handle.write(f"{indent}endTime    {end_time_value};\n")
                replaced = True
            else:
                handle.write(line)

        if not replaced:
            handle.write(f"\nendTime    {end_time_value};\n")


def _build_cases(
    ionic_models: Sequence[str],
    ionic_model_tissue_map: Mapping[str, Sequence[str]],
    dt_values: Sequence[float],
    dx_values: Sequence[float],
    solvers: Sequence[str],
) -> list[CaseConfig]:
    cases = []
    for ionic_model in ionic_models:
        tissues = ionic_model_tissue_map.get(ionic_model)
        if not tissues:
            raise KeyError(f"Missing tissue mapping for ionic model '{ionic_model}'")
        for tissue, dt, dx, solver in product(tissues, dt_values, dx_values, solvers):
            case_id = f"{solver}_{ionic_model}_{tissue}_DT{dt}_DX{dx}"
            cases.append(
                CaseConfig(
                    case_id=case_id,
                    params={
                        "ionicModel": ionic_model,
                        "tissue": tissue,
                        "dt_ms": dt,
                        "dx_mm": dx,
                        "solver": solver,
                    },
                )
            )
    return cases


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    control_dict_relpath: Path = defaults.CONTROL_DICT_RELPATH,
    block_mesh_dict_relpath: Path = defaults.BLOCK_MESH_DICT_RELPATH,
    electro_properties_relpath: Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    slab_size_mm: Sequence[float] = defaults.SLAB_SIZE_MM,
    end_time_by_dx: Mapping[float, float] = defaults.END_TIME_BY_DX,
) -> None:
    control_dict = case_root / control_dict_relpath
    block_mesh_dict = case_root / block_mesh_dict_relpath
    electro_properties = case_root / electro_properties_relpath

    dx_mm = float(case.params["dx_mm"])
    dt_ms = float(case.params["dt_ms"])
    tissue = str(case.params["tissue"])
    ionic_model = str(case.params["ionicModel"])
    solver = str(case.params["solver"])

    _replace_blockmesh_resolution(block_mesh_dict, dx_mm, slab_size_mm=slab_size_mm)
    set_delta_t(control_dict, dt_ms)
    _update_end_time(control_dict, dx_mm, end_time_by_dx=end_time_by_dx)
    set_tissue(electro_properties, tissue, scope=electro_properties_scope)
    set_ionic_model(electro_properties, ionic_model, scope=electro_properties_scope)
    set_solution_algorithm(
        electro_properties,
        solver,
        scope=electro_properties_scope,
    )


def _run_case(
    case_root: Path,
    setup_root: Path,
    case: CaseConfig,
    *,
    tutorials_root: Path | None = None,
    run_script_relpath: Path = defaults.RUN_SCRIPT_RELPATH,
    paraview_script_relpath: Path = defaults.PARAVIEW_SCRIPT_RELPATH,
    points_file_relpath: Path = defaults.POINTS_FILE_RELPATH,
    output_relpath: Path = defaults.OUTPUT_RELPATH,
    case_foam_relpath: Path = defaults.CASE_FOAM_RELPATH,
    pvpython_env_var: str = defaults.PVPYTHON_ENV_VAR,
    pvpython_executable: Path = defaults.PVPYTHON_EXECUTABLE,
) -> None:
    run_script = resolve_run_script_path(
        tutorials_root=tutorials_root,
        run_script_relpath=run_script_relpath,
    )
    paraview_script = setup_root / paraview_script_relpath
    points_file = setup_root / points_file_relpath
    output_dir = case_root / output_relpath
    output_dir.mkdir(parents=True, exist_ok=True)

    pvpython_path = Path(
        os.environ.get(
            pvpython_env_var,
            str(pvpython_executable),
        )
    )
    if not pvpython_path.exists():
        raise FileNotFoundError(
            f"pvpython executable not found at '{pvpython_path}'. "
            f"Set {pvpython_env_var} to your ParaView pvpython binary."
        )

    subprocess.run(
        [
            "bash",
            "-l",
            str(run_script),
            "--case-dir",
            str(case_root),
            "--parallel",
            "--touch-case-foam",
        ],
        check=True,
    )
    case_foam = case_root / case_foam_relpath
    case_foam.touch()

    subprocess.run(
        [
            str(pvpython_path),
            str(paraview_script),
            "--case",
            str(case_foam),
            "--points",
            str(points_file),
            "--outdir",
            str(output_dir),
            "--dx",
            str(case.params["dx_mm"]),
            "--dt",
            str(case.params["dt_ms"]),
            "--tissue",
            str(case.params["tissue"]),
            "--ionicModel",
            str(case.params["ionicModel"]),
            "--solver",
            str(case.params["solver"]),
        ],
        check=True,
    )


def _postprocess(
    setup_root: Path,
    output_dir: Path,
    *,
    line_postprocess_relpath: Path = defaults.LINE_POSTPROCESS_RELPATH,
    points_postprocess_relpath: Path = defaults.POINTS_POSTPROCESS_RELPATH,
    excel_reference_relpath: Path = defaults.EXCEL_REFERENCE_RELPATH,
    line_postprocess_function_name: str = defaults.LINE_POSTPROCESS_FUNCTION,
    points_postprocess_function_name: str = defaults.POINTS_POSTPROCESS_FUNCTION,
    strict_artifacts: bool = False,
) -> None:
    run_postprocess_tasks(
        setup_root=setup_root,
        output_dir=output_dir,
        tutorial_name=defaults.TUTORIAL_NAME,
        strict_artifacts=strict_artifacts,
        tasks=[
            PostprocessTask(
                module_relpath=line_postprocess_relpath,
                function_name=line_postprocess_function_name,
                kwargs={"excel_path": f"$SETUP_ROOT/{excel_reference_relpath.as_posix()}"},
            ),
            PostprocessTask(
                module_relpath=points_postprocess_relpath,
                function_name=points_postprocess_function_name,
            ),
        ],
    )


def make_spec(
    *,
    tutorials_root: Path | None = None,
    case_dir_name: str = defaults.CASE_DIR_NAME,
    setup_dir_name: str | None = None,
    output_dir_name: str = defaults.OUTPUT_DIR_NAME,
    ionic_models: Sequence[str] = defaults.IONIC_MODELS,
    ionic_model_tissue_map: Mapping[str, Sequence[str]] = defaults.IONIC_MODEL_TISSUE_MAP,
    dt_values: Sequence[float] = defaults.DT_VALUES,
    dx_values: Sequence[float] = defaults.DX_VALUES,
    solvers: Sequence[str] = defaults.SOLVERS,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    slab_size_mm: Sequence[float] = defaults.SLAB_SIZE_MM,
    end_time_by_dx: Mapping[float, float] = defaults.END_TIME_BY_DX,
    control_dict_relpath: str | Path = defaults.CONTROL_DICT_RELPATH,
    block_mesh_dict_relpath: str | Path = defaults.BLOCK_MESH_DICT_RELPATH,
    electro_properties_relpath: str | Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    paraview_script_relpath: str | Path = defaults.PARAVIEW_SCRIPT_RELPATH,
    points_file_relpath: str | Path = defaults.POINTS_FILE_RELPATH,
    case_foam_relpath: str | Path = defaults.CASE_FOAM_RELPATH,
    pvpython_env_var: str = defaults.PVPYTHON_ENV_VAR,
    pvpython_executable: str | Path = defaults.PVPYTHON_EXECUTABLE,
    line_postprocess_relpath: str | Path = defaults.LINE_POSTPROCESS_RELPATH,
    points_postprocess_relpath: str | Path = defaults.POINTS_POSTPROCESS_RELPATH,
    excel_reference_relpath: str | Path = defaults.EXCEL_REFERENCE_RELPATH,
    line_postprocess_function_name: str = defaults.LINE_POSTPROCESS_FUNCTION,
    points_postprocess_function_name: str = defaults.POINTS_POSTPROCESS_FUNCTION,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    ionic_models_list = [str(item) for item in ionic_models]
    ionic_model_tissue_map_normalized: dict[str, list[str]] = {}
    for ionic_model in ionic_models_list:
        tissues = ionic_model_tissue_map.get(ionic_model)
        if not tissues:
            raise KeyError(f"Missing tissue mapping for ionic model '{ionic_model}'")
        ionic_model_tissue_map_normalized[ionic_model] = [str(item) for item in tissues]

    dt_values_list = [float(item) for item in dt_values]
    dx_values_list = [float(item) for item in dx_values]
    solvers_list = [str(item) for item in solvers]
    slab_size_mm_list = [float(item) for item in slab_size_mm]
    end_time_by_dx_map = {float(key): float(value) for key, value in end_time_by_dx.items()}

    if not ionic_models_list:
        raise ValueError("ionic_models cannot be empty")
    if not dt_values_list:
        raise ValueError("dt_values cannot be empty")
    if not dx_values_list:
        raise ValueError("dx_values cannot be empty")
    if not solvers_list:
        raise ValueError("solvers cannot be empty")
    if len(slab_size_mm_list) != 3:
        raise ValueError("slab_size_mm must contain exactly 3 entries")
    if any(value <= 0 for value in slab_size_mm_list):
        raise ValueError("slab_size_mm entries must be positive")
    if not end_time_by_dx_map:
        raise ValueError("end_time_by_dx cannot be empty")

    control_dict_path = Path(control_dict_relpath)
    block_mesh_dict_path = Path(block_mesh_dict_relpath)
    electro_properties_path = Path(electro_properties_relpath)
    run_script_path = Path(run_script_relpath)
    paraview_script_path = Path(paraview_script_relpath)
    points_file_path = Path(points_file_relpath)
    case_foam_path = Path(case_foam_relpath)
    pvpython_path = Path(pvpython_executable)
    line_postprocess_path = Path(line_postprocess_relpath)
    points_postprocess_path = Path(points_postprocess_relpath)
    excel_reference_path = Path(excel_reference_relpath)
    output_relpath = Path(output_dir_name)

    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
    )

    return TutorialSpec(
        name=case_dir_name,
        case_root=case_root,
        setup_root=setup_root,
        output_dir=output_dir,
        build_cases=partial(
            _build_cases,
            ionic_models=ionic_models_list,
            ionic_model_tissue_map=ionic_model_tissue_map_normalized,
            dt_values=dt_values_list,
            dx_values=dx_values_list,
            solvers=solvers_list,
        ),
        apply_case=partial(
            _apply_case,
            electro_properties_scope=electro_properties_scope,
            control_dict_relpath=control_dict_path,
            block_mesh_dict_relpath=block_mesh_dict_path,
            electro_properties_relpath=electro_properties_path,
            slab_size_mm=slab_size_mm_list,
            end_time_by_dx=end_time_by_dx_map,
        ),
        run_case=partial(
            _run_case,
            tutorials_root=tutorials_root,
            run_script_relpath=run_script_path,
            paraview_script_relpath=paraview_script_path,
            points_file_relpath=points_file_path,
            output_relpath=output_relpath,
            case_foam_relpath=case_foam_path,
            pvpython_env_var=pvpython_env_var,
            pvpython_executable=pvpython_path,
        ),
        collect_outputs=None,
        postprocess=partial(
            _postprocess,
            line_postprocess_relpath=line_postprocess_path,
            points_postprocess_relpath=points_postprocess_path,
            excel_reference_relpath=excel_reference_path,
            line_postprocess_function_name=line_postprocess_function_name,
            points_postprocess_function_name=points_postprocess_function_name,
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "notes": "Niederer Et Al. 2012 slab benchmark sweep.",
            "pvpython_env": pvpython_env_var,
            "pvpython_executable": str(pvpython_path),
            "dx_values": dx_values_list,
            "dt_values": dt_values_list,
            "slab_size_mm": slab_size_mm_list,
            "ionic_model_tissue_map": ionic_model_tissue_map_normalized,
            "ionic_models": ionic_models_list,
            "solvers": solvers_list,
            "control_dict_relpath": str(control_dict_path),
            "block_mesh_dict_relpath": str(block_mesh_dict_path),
            "electro_properties_relpath": str(electro_properties_path),
            "electro_properties_scope": electro_properties_scope,
            "run_script_relpath": str(run_script_path),
            "paraview_script_relpath": str(paraview_script_path),
            "points_file_relpath": str(points_file_path),
            "output_relpath": str(output_relpath),
            "line_postprocess_relpath": str(line_postprocess_path),
            "points_postprocess_relpath": str(points_postprocess_path),
            "excel_reference_relpath": str(excel_reference_path),
            "line_postprocess_function_name": line_postprocess_function_name,
            "points_postprocess_function_name": points_postprocess_function_name,
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
