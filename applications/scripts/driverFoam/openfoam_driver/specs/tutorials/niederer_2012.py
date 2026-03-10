from __future__ import annotations

import csv
import math
import re
import shutil
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

DEFAULT_POINTS_FUNCTION_OBJECT = getattr(
    defaults, "NIEDERER_POINTS_FUNCTION_OBJECT", "Niedererpoints"
)
DEFAULT_LINE_FUNCTION_OBJECT = getattr(
    defaults, "NIEDERER_LINE_FUNCTION_OBJECT", "Niedererlines"
)
DEFAULT_SAMPLED_FIELD = getattr(defaults, "NIEDERER_SAMPLED_FIELD", "activationTime")
DEFAULT_SAMPLED_POINTS = getattr(
    defaults,
    "NIEDERER_POINTS",
    (
        ("P1", 0.0, 0.0, 0.007),
        ("P2", 0.0, 0.0, 0.0),
        ("P3", 0.019999, 0.0, 0.007),
        ("P4", 0.019999, 0.0, 0.0),
        ("P5", 0.0, 0.003, 0.007),
        ("P6", 0.0, 0.003, 0.0),
        ("P7", 0.019999, 0.003, 0.007),
        ("P8", 0.019999, 0.003, 0.0),
        ("P9", 0.01, 0.0015, 0.0035),
    ),
)
DEFAULT_LINE_START = getattr(defaults, "NIEDERER_LINE_START", (0.0, 0.0, 0.007))
DEFAULT_LINE_END = getattr(defaults, "NIEDERER_LINE_END", (0.02, 0.003, 0.0))
DEFAULT_LINE_N_POINTS = int(getattr(defaults, "NIEDERER_LINE_NUM_POINTS", 101))


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
    # Input dt is provided in milliseconds in the JSON/spec settings.
    set_delta_t(control_dict, dt_ms * 1.0e-3)
    _update_end_time(control_dict, dx_mm, end_time_by_dx=end_time_by_dx)
    set_tissue(electro_properties, tissue, scope=electro_properties_scope)
    set_ionic_model(electro_properties, ionic_model, scope=electro_properties_scope)
    set_solution_algorithm(
        electro_properties,
        solver,
        scope=electro_properties_scope,
    )


def _float_to_case_tag(value: float) -> str:
    tag = f"{value:.7f}".rstrip("0").rstrip(".").replace(".", "")
    return tag or "0"


def _read_latest_probe_values(
    *,
    postprocess_root: Path,
    function_object_name: str,
    sampled_field: str,
) -> list[float]:
    function_root = postprocess_root / function_object_name
    if not function_root.exists():
        raise FileNotFoundError(
            f"Missing functionObject output folder: {function_root}"
        )

    time_dirs: list[tuple[float, Path]] = []
    for child in function_root.iterdir():
        if not child.is_dir():
            continue
        try:
            time_value = float(child.name)
        except ValueError:
            continue
        time_dirs.append((time_value, child))

    if not time_dirs:
        raise FileNotFoundError(
            f"No numeric time folders found in: {function_root}"
        )

    _, latest_time_dir = max(time_dirs, key=lambda item: item[0])
    sampled_file = latest_time_dir / sampled_field
    if not sampled_file.exists():
        raise FileNotFoundError(f"Missing sampled field file: {sampled_file}")

    numeric_rows: list[list[float]] = []
    for line in sampled_file.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            numeric_rows.append([float(token) for token in stripped.split()])
        except ValueError:
            continue

    if not numeric_rows or len(numeric_rows[-1]) < 2:
        raise ValueError(f"No probe values found in {sampled_file}")

    return numeric_rows[-1][1:]


def _write_points_csv(
    *,
    output_dir: Path,
    file_name: str,
    sampled_points: Sequence[tuple[str, float, float, float]],
    activation_values: Sequence[float],
) -> None:
    if len(sampled_points) != len(activation_values):
        raise ValueError(
            f"Points count mismatch: expected {len(sampled_points)}, got {len(activation_values)}"
        )

    output_file = output_dir / file_name
    with output_file.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Label", "Points:0", "Points:1", "Points:2", "activationTime"])
        for (label, x, y, z), activation_time in zip(sampled_points, activation_values):
            writer.writerow([label, x, y, z, activation_time])


def _write_line_csv(
    *,
    output_dir: Path,
    file_name: str,
    activation_values: Sequence[float],
    line_start: Sequence[float],
    line_end: Sequence[float],
    expected_n_points: int,
) -> None:
    if len(activation_values) != expected_n_points:
        raise ValueError(
            f"Line sample count mismatch: expected {expected_n_points}, got {len(activation_values)}"
        )

    x0, y0, z0 = [float(value) for value in line_start]
    x1, y1, z1 = [float(value) for value in line_end]
    line_length = math.dist((x0, y0, z0), (x1, y1, z1))
    n_intervals = max(expected_n_points - 1, 1)

    output_file = output_dir / file_name
    with output_file.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["activationTime", "arc_length", "Points_0", "Points_1", "Points_2"])
        for index, activation_time in enumerate(activation_values):
            t = index / n_intervals
            x = x0 + t * (x1 - x0)
            y = y0 + t * (y1 - y0)
            z = z0 + t * (z1 - z0)
            arc_length = t * line_length
            writer.writerow([activation_time, arc_length, x, y, z])


def _export_openfoam_samples(
    *,
    case_root: Path,
    output_dir: Path,
    case: CaseConfig,
    points_function_object_name: str,
    line_function_object_name: str,
    sampled_field: str,
    sampled_points: Sequence[tuple[str, float, float, float]],
    line_start: Sequence[float],
    line_end: Sequence[float],
    line_n_points: int,
) -> None:
    postprocess_root = case_root / "postProcessing"
    points_values = _read_latest_probe_values(
        postprocess_root=postprocess_root,
        function_object_name=points_function_object_name,
        sampled_field=sampled_field,
    )
    line_values = _read_latest_probe_values(
        postprocess_root=postprocess_root,
        function_object_name=line_function_object_name,
        sampled_field=sampled_field,
    )

    dx_tag = _float_to_case_tag(float(case.params["dx_mm"]))
    dt_tag = _float_to_case_tag(float(case.params["dt_ms"]))
    solver = str(case.params["solver"])
    ionic_model = str(case.params["ionicModel"])
    tissue = str(case.params["tissue"])

    points_file_name = f"{solver}_{ionic_model}_{tissue}_points_DT{dt_tag}_DX{dx_tag}.csv"
    line_file_name = f"{solver}_{ionic_model}_{tissue}_line_DT{dt_tag}_DX{dx_tag}.csv"

    _write_points_csv(
        output_dir=output_dir,
        file_name=points_file_name,
        sampled_points=sampled_points,
        activation_values=points_values,
    )
    _write_line_csv(
        output_dir=output_dir,
        file_name=line_file_name,
        activation_values=line_values,
        line_start=line_start,
        line_end=line_end,
        expected_n_points=line_n_points,
    )


def _run_case(
    case_root: Path,
    setup_root: Path,
    case: CaseConfig,
    *,
    tutorials_root: Path | None = None,
    run_script_relpath: Path = defaults.RUN_SCRIPT_RELPATH,
    output_relpath: Path = defaults.OUTPUT_RELPATH,
    points_function_object_name: str = DEFAULT_POINTS_FUNCTION_OBJECT,
    line_function_object_name: str = DEFAULT_LINE_FUNCTION_OBJECT,
    sampled_field: str = DEFAULT_SAMPLED_FIELD,
    sampled_points: Sequence[tuple[str, float, float, float]] = DEFAULT_SAMPLED_POINTS,
    line_start: Sequence[float] = DEFAULT_LINE_START,
    line_end: Sequence[float] = DEFAULT_LINE_END,
    line_n_points: int = DEFAULT_LINE_N_POINTS,
) -> None:
    run_script = resolve_run_script_path(
        tutorials_root=tutorials_root,
        run_script_relpath=run_script_relpath,
    )
    output_dir = case_root / output_relpath
    output_dir.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "bash",
            "-l",
            str(run_script),
            "--case-dir",
            str(case_root),
            "--parallel",
        ],
        check=True,
    )

    _export_openfoam_samples(
        case_root=case_root,
        output_dir=output_dir,
        case=case,
        points_function_object_name=points_function_object_name,
        line_function_object_name=line_function_object_name,
        sampled_field=sampled_field,
        sampled_points=sampled_points,
        line_start=line_start,
        line_end=line_end,
        line_n_points=line_n_points,
    )

    _cache_case_postprocessing(
        case_root=case_root,
        setup_root=setup_root,
        case_id=case.case_id,
        cache_dir_name=defaults.CASE_POSTPROCESS_CACHE_DIRNAME,
    )


def _cache_case_postprocessing(
    *,
    case_root: Path,
    setup_root: Path,
    case_id: str,
    cache_dir_name: str,
) -> None:
    source = case_root / "postProcessing"
    if not source.exists():
        raise FileNotFoundError(f"Missing postProcessing folder to cache: {source}")

    cache_root = setup_root / cache_dir_name
    cache_root.mkdir(parents=True, exist_ok=True)
    destination = cache_root / case_id

    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)
    print(f"Cached postProcessing for {case_id}: {destination}")


def _postprocess(
    setup_root: Path,
    output_dir: Path,
    *,
    line_postprocess_relpath: Path = defaults.LINE_POSTPROCESS_RELPATH,
    points_postprocess_relpath: Path = defaults.POINTS_POSTPROCESS_RELPATH,
    cache_postprocess_relpath: Path = defaults.CACHE_POSTPROCESS_RELPATH,
    excel_reference_relpath: Path = defaults.EXCEL_REFERENCE_RELPATH,
    cache_postprocess_function_name: str = defaults.CACHE_POSTPROCESS_FUNCTION,
    line_postprocess_function_name: str = defaults.LINE_POSTPROCESS_FUNCTION,
    points_postprocess_function_name: str = defaults.POINTS_POSTPROCESS_FUNCTION,
    case_postprocess_cache_dirname: str = defaults.CASE_POSTPROCESS_CACHE_DIRNAME,
    strict_artifacts: bool = False,
) -> None:
    run_postprocess_tasks(
        setup_root=setup_root,
        output_dir=output_dir,
        tutorial_name=defaults.TUTORIAL_NAME,
        strict_artifacts=strict_artifacts,
        tasks=[
            PostprocessTask(
                module_relpath=cache_postprocess_relpath,
                function_name=cache_postprocess_function_name,
                kwargs={
                    "cache_root": f"$SETUP_ROOT/{case_postprocess_cache_dirname}",
                    "cache_output_subdir": "cachedPostProcessing",
                },
            ),
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
    points_function_object_name: str = DEFAULT_POINTS_FUNCTION_OBJECT,
    line_function_object_name: str = DEFAULT_LINE_FUNCTION_OBJECT,
    sampled_field: str = DEFAULT_SAMPLED_FIELD,
    sampled_points: Sequence[tuple[str, float, float, float]] = DEFAULT_SAMPLED_POINTS,
    line_start: Sequence[float] = DEFAULT_LINE_START,
    line_end: Sequence[float] = DEFAULT_LINE_END,
    line_n_points: int = DEFAULT_LINE_N_POINTS,
    line_postprocess_relpath: str | Path = defaults.LINE_POSTPROCESS_RELPATH,
    points_postprocess_relpath: str | Path = defaults.POINTS_POSTPROCESS_RELPATH,
    cache_postprocess_relpath: str | Path = defaults.CACHE_POSTPROCESS_RELPATH,
    excel_reference_relpath: str | Path = defaults.EXCEL_REFERENCE_RELPATH,
    cache_postprocess_function_name: str = defaults.CACHE_POSTPROCESS_FUNCTION,
    line_postprocess_function_name: str = defaults.LINE_POSTPROCESS_FUNCTION,
    points_postprocess_function_name: str = defaults.POINTS_POSTPROCESS_FUNCTION,
    case_postprocess_cache_dirname: str = defaults.CASE_POSTPROCESS_CACHE_DIRNAME,
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
    line_postprocess_path = Path(line_postprocess_relpath)
    points_postprocess_path = Path(points_postprocess_relpath)
    cache_postprocess_path = Path(cache_postprocess_relpath)
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
            output_relpath=output_relpath,
            points_function_object_name=points_function_object_name,
            line_function_object_name=line_function_object_name,
            sampled_field=sampled_field,
            sampled_points=sampled_points,
            line_start=line_start,
            line_end=line_end,
            line_n_points=line_n_points,
        ),
        collect_outputs=None,
        postprocess=partial(
            _postprocess,
            line_postprocess_relpath=line_postprocess_path,
            points_postprocess_relpath=points_postprocess_path,
            cache_postprocess_relpath=cache_postprocess_path,
            excel_reference_relpath=excel_reference_path,
            cache_postprocess_function_name=cache_postprocess_function_name,
            line_postprocess_function_name=line_postprocess_function_name,
            points_postprocess_function_name=points_postprocess_function_name,
            case_postprocess_cache_dirname=case_postprocess_cache_dirname,
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "notes": (
                "Niederer Et Al. 2012 slab benchmark sweep "
                "with OpenFOAM functionObject sampling."
            ),
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
            "output_relpath": str(output_relpath),
            "points_function_object_name": points_function_object_name,
            "line_function_object_name": line_function_object_name,
            "sampled_field": sampled_field,
            "sampled_points_count": len(sampled_points),
            "line_start": [float(value) for value in line_start],
            "line_end": [float(value) for value in line_end],
            "line_n_points": line_n_points,
            "line_postprocess_relpath": str(line_postprocess_path),
            "points_postprocess_relpath": str(points_postprocess_path),
            "cache_postprocess_relpath": str(cache_postprocess_path),
            "excel_reference_relpath": str(excel_reference_path),
            "cache_postprocess_function_name": cache_postprocess_function_name,
            "line_postprocess_function_name": line_postprocess_function_name,
            "points_postprocess_function_name": points_postprocess_function_name,
            "case_postprocess_cache_dirname": case_postprocess_cache_dirname,
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
