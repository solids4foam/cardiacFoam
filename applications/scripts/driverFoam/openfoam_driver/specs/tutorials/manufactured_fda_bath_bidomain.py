from __future__ import annotations

import shutil
import subprocess
from collections.abc import Sequence
from functools import partial
from pathlib import Path

from ...core.defaults import manufactured_fda_bath_bidomain as defaults
from ...core.runtime.models import CaseConfig, TutorialSpec
from ...postprocessing.driver import PostprocessTask, run_postprocess_tasks
from ..common import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
    remove_electro_property_dict,
    resolve_run_script_path,
    resolve_spec_paths,
    set_delta_t,
)
from .manufactured_fda import _build_cases


def _case_output_filename(case: CaseConfig) -> str:
    dimension = str(case.params["dimension"])
    cells = int(case.params["cells"])
    solver = str(case.params["solver"])
    return f"bathBidomain_{dimension}_{cells}_cells_{solver}.dat"


def _archive_output_dir(case_root: Path) -> Path:
    return case_root / "archivedPostProcessing"


def _replace_blockmesh_resolution(block_mesh_dict_path: Path, cells: int, dimension: str) -> None:
    if not block_mesh_dict_path.exists():
        raise FileNotFoundError(f"Missing mesh dictionary: {block_mesh_dict_path}")

    try:
        cell_counts = defaults.BLOCK_MESH_RESOLUTION_BY_DIMENSION[dimension].format(cells=cells)
    except KeyError as exc:
        raise ValueError(f"Unsupported dimension: {dimension}") from exc

    lines = block_mesh_dict_path.read_text().splitlines(keepends=True)
    replaced = 0

    with block_mesh_dict_path.open("w") as handle:
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("hex (") and not stripped.startswith("//"):
                prefix, _, suffix = line.partition(") (")
                if not suffix:
                    handle.write(line)
                    continue
                _, _, trailing = suffix.partition(") simpleGrading")
                handle.write(f"{prefix}) ({cell_counts}) simpleGrading{trailing}")
                replaced += 1
            else:
                handle.write(line)

    if replaced != 3:
        raise KeyError(
            f"Expected to update 3 hex blocks in {block_mesh_dict_path}, "
            f"updated {replaced}."
        )


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    control_dict_relpath: Path = Path("system/controlDict"),
    electro_properties_relpath: Path = Path("constant/electroProperties"),
    physics_properties_relpath: Path = Path("constant/physicsProperties"),
    electro_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    physics_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    verification_model_type: str = defaults.VERIFICATION_MODEL_TYPE,
    ecg_enabled: bool = False,
    block_mesh_dict_template: str = defaults.BLOCK_MESH_DICT_TEMPLATE,
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
        f"{electro_properties_scope}.verificationModel.groundElectrode": True,
        f"{electro_properties_scope}.manufacturedBidomain.groundElectrode": True,
    }

    if ecg_enabled:
        case_overrides.update(
            {
                f"{electro_properties_scope}.ecgDomains.bodyECG.ecgSolver": "torsoECG",
                f"{electro_properties_scope}.ecgDomains.bodyECG.ecgVerificationModel":
                    "bathECGManufacturedVerifier",
                f"{electro_properties_scope}.ecgDomains.pseudoECGSignals.ecgSolver": "pseudoECG",
            }
        )

    _replace_blockmesh_resolution(block_mesh_dict, cells, dimension)
    set_delta_t(control_dict, dt_value)
    apply_electro_property_overrides(electro_properties, case_overrides)
    if not ecg_enabled:
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
    ecg_enabled: bool = False,
) -> None:
    del setup_root
    dimension = str(case.params["dimension"])
    run_script = resolve_run_script_path(
        tutorials_root=tutorials_root,
        run_script_relpath=run_script_relpath,
    )
    command = ["bash", "-l", str(run_script), "--case-dir", str(case_root), "--dimension", dimension]
    if run_in_parallel:
        command.append("--parallel")

    subprocess.run(command, check=True)
    _archive_case_logs(case_root, case)
    _stage_case_output(case_root, case, run_in_parallel=run_in_parallel)
    if ecg_enabled:
        _stage_case_ecg_outputs(case_root, case, run_in_parallel=run_in_parallel)


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
    run_in_parallel: bool = False,
) -> Path:
    filename = _case_output_filename(case)
    destination_dir = _archive_output_dir(case_root)
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / filename
    if run_in_parallel:
        candidates = (
            case_root / "processor0" / "postProcessing" / filename,
            case_root / "postProcessing" / filename,
        )
    else:
        candidates = (
            case_root / "postProcessing" / filename,
            case_root / "processor0" / "postProcessing" / filename,
        )

    for candidate in candidates:
        if not candidate.exists():
            continue
        if candidate != destination:
            shutil.copy2(candidate, destination)
            print(f"Archived bath manufactured output: {candidate} -> {destination}")
        return destination

    checked = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(
        f"Bath manufactured output '{filename}' not found after run. Checked: {checked}"
    )


def _stage_case_ecg_outputs(
    case_root: Path,
    case: CaseConfig,
    *,
    run_in_parallel: bool = False,
) -> list[Path]:
    staged_outputs: list[Path] = []
    destination_dir = _archive_output_dir(case_root)
    destination_dir.mkdir(parents=True, exist_ok=True)

    ecg_outputs = (
        ("BathECG", "torsoECG.dat"),
        ("BathECG", "manufacturedBathECG.dat"),
        ("BathECG", "manufacturedBathECGSummary.dat"),
        ("PseudoECG", "pseudoECG.dat"),
    )

    for prefix, source_name in ecg_outputs:
        destination = destination_dir / f"{prefix}_{case.case_id}_{source_name}"
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
            if candidate.parent == destination_dir:
                shutil.move(str(candidate), str(destination))
                print(f"Archived bath ECG output: {candidate} -> {destination} (moved)")
            else:
                shutil.copy2(candidate, destination)
                print(f"Archived bath ECG output: {candidate} -> {destination}")
            staged_outputs.append(destination)
            break

    return staged_outputs


def _collect_outputs(case_root: Path, output_dir: Path) -> None:
    archived_dir = _archive_output_dir(case_root)
    archived_outputs = []
    if archived_dir.exists():
        for source in sorted(archived_dir.glob("*.dat")):
            if source.name.startswith(("bathBidomain_", "BathECG_")):
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
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    control_dict_relpath: str | Path = "system/controlDict",
    electro_properties_relpath: str | Path = "constant/electroProperties",
    physics_properties_relpath: str | Path = "constant/physicsProperties",
    electro_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    physics_property_overrides: Sequence[dict[str, object]] | dict[str, object] | None = None,
    verification_model_type: str = defaults.VERIFICATION_MODEL_TYPE,
    block_mesh_dict_template: str = defaults.BLOCK_MESH_DICT_TEMPLATE,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    postprocess_script_relpath: str | Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    run_in_parallel: bool = defaults.RUN_IN_PARALLEL,
    ecg_enabled: bool = False,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    dimensions_list = [str(item) for item in dimensions]
    cells_list = [int(item) for item in number_cells]
    dt_values_list = [float(item) for item in dt_values]
    solver_types_list = [str(item) for item in solver_types]

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
            control_dict_relpath=Path(control_dict_relpath),
            electro_properties_relpath=Path(electro_properties_relpath),
            physics_properties_relpath=Path(physics_properties_relpath),
            electro_property_overrides=electro_property_overrides,
            physics_property_overrides=physics_property_overrides,
            verification_model_type=verification_model_type,
            ecg_enabled=ecg_enabled,
            block_mesh_dict_template=block_mesh_dict_template,
        ),
        run_case=partial(
            _run_case,
            tutorials_root=tutorials_root,
            run_script_relpath=Path(run_script_relpath),
            run_in_parallel=run_in_parallel,
            ecg_enabled=ecg_enabled,
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
            "notes": "FDA bath-bidomain manufactured-solution convergence benchmark",
            "workflow_dag": {
                "steps": [
                    {"id": "mesh", "command": "blockMesh", "depends_on": []},
                    {"id": "topoSet", "command": "topoSet", "depends_on": ["mesh"]},
                    {
                        "id": "setConductivity",
                        "command": "setTorsoOrganConductivityField",
                        "depends_on": ["topoSet"],
                    },
                    {"id": "solve", "command": "cardiacFoam", "depends_on": ["setConductivity"]},
                ]
            },
            "dimensions": dimensions_list,
            "solver_types": solver_types_list,
            "piecewise_sweep": piecewise_sweep,
            "control_dict_relpath": str(control_dict_relpath),
            "electro_properties_relpath": str(electro_properties_relpath),
            "physics_properties_relpath": str(physics_properties_relpath),
            "electro_properties_scope": electro_properties_scope,
            "block_mesh_dict_template": block_mesh_dict_template,
            "run_script_relpath": str(run_script_relpath),
            "run_in_parallel": run_in_parallel,
            "ecg_enabled": ecg_enabled,
            "postprocess_script_relpath": str(postprocess_script_relpath),
            "postprocess_function_name": postprocess_function_name,
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
