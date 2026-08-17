from __future__ import annotations

import subprocess
import sys
from collections.abc import Mapping, Sequence
from functools import partial
from itertools import product
from pathlib import Path

from openfoam_driver.plugins.cardiacfoam.defaults import cable_1d_restitution as defaults
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec, DataArtifact
from openfoam_driver.postprocessing.driver import PostprocessTask, run_postprocess_tasks
from openfoam_driver.plugins.cardiacfoam.overrides import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
)
from openfoam_driver.specs.common import (
    load_python_module,
    replace_block_mesh_resolutions,
    resolve_run_script_path,
    resolve_spec_paths,
    set_delta_t,
    set_end_time,
)
from openfoam_driver.specs.spatial_pacing import generate_spatial_s1_s2_stimulus_lists
from openfoam_driver.specs.mesh_provisioning import cell_counts_from_dx



def _build_cases(
    ionic_models: Sequence[str],
    ionic_model_tissue_map: Mapping[str, Sequence[str]],
    dt_values: Sequence[float],
    dx_values: Sequence[float],
    solvers: Sequence[str],
    conductivity_values: Sequence[str],
    s2_intervals_ms: Sequence[float],
) -> list[CaseConfig]:
    cases: list[CaseConfig] = []
    for conductivity_index, conductivity in enumerate(conductivity_values, start=1):
        for ionic_model in ionic_models:
            tissues = ionic_model_tissue_map.get(ionic_model)
            if not tissues:
                raise KeyError(f"Missing tissue mapping for ionic model '{ionic_model}'")
            for tissue, dt, dx, solver, s2_interval in product(
                tissues, dt_values, dx_values, solvers, s2_intervals_ms
            ):
                case_id = (
                    f"{solver}_{ionic_model}_{tissue}_DT{dt:g}_DX{dx:g}_"
                    f"COND{conductivity_index:02d}_S2{s2_interval:g}"
                )
                cases.append(
                    CaseConfig(
                        case_id=case_id,
                        params={
                            "ionicModel": ionic_model,
                            "tissue": tissue,
                            "dt_ms": float(dt),
                            "dx_mm": float(dx),
                            "solver": solver,
                            "conductivity": conductivity,
                            "conductivity_id": conductivity_index,
                            "s2Interval": float(s2_interval),
                        },
                    )
                )
    return cases


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    s1_interval_ms: float = defaults.S1_INTERVAL_MS,
    n_s1: int = defaults.N_S1,
    n_s2: int = defaults.N_S2,
    end_time_buffer_s: float = defaults.END_TIME_BUFFER_S,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    control_dict_relpath: Path = defaults.CONTROL_DICT_RELPATH,
    block_mesh_dict_relpath: Path = defaults.BLOCK_MESH_DICT_RELPATH,
    electro_properties_relpath: Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    physics_properties_relpath: Path = Path("constant/physicsProperties"),
    electro_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    physics_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    cable_length_mm: float = defaults.CABLE_LENGTH_MM,
    cross_section_cell_counts: Sequence[int] = defaults.CROSS_SECTION_CELL_COUNTS,
) -> None:
    control_dict = case_root / control_dict_relpath
    block_mesh_dict = case_root / block_mesh_dict_relpath
    electro_properties = case_root / electro_properties_relpath
    physics_properties = case_root / physics_properties_relpath

    (cells,) = cell_counts_from_dx(float(case.params["dx_mm"]), (cable_length_mm,))
    cell_counts_str = f"{cells} {int(cross_section_cell_counts[0])} {int(cross_section_cell_counts[1])}"
    replace_block_mesh_resolutions(block_mesh_dict, cell_counts_str)

    set_delta_t(control_dict, float(case.params["dt_ms"]) * 1.0e-3)
    
    end_time = (s1_interval_ms * n_s1 + case.params["s2Interval"] * n_s2) / 1000.0 + end_time_buffer_s
    set_end_time(control_dict, end_time)

    stimulus_arrays = generate_spatial_s1_s2_stimulus_lists(
        s1_interval_ms=s1_interval_ms, n_s1=n_s1,
        s2_interval_ms=case.params["s2Interval"], n_s2=n_s2,
        bounds_min="(0 0 0)", bounds_max="(2e-3 2e-4 2e-4)",
        duration_s="4e-3", intensity="50000"
    )

    case_overrides = {
        f"{electro_properties_scope}.conductivity": str(case.params["conductivity"]),
        f"{electro_properties_scope}.tissue": str(case.params["tissue"]),
        f"{electro_properties_scope}.ionicModel": str(case.params["ionicModel"]),
        f"{electro_properties_scope}.solutionAlgorithm": str(case.params["solver"]),
    }
    for k, v in stimulus_arrays.items():
        case_overrides[f"{electro_properties_scope}.externalStimulus.{k}"] = v

    apply_electro_property_overrides(electro_properties, case_overrides)
    apply_electro_property_overrides(electro_properties, electro_property_overrides)
    apply_physics_property_overrides(physics_properties, physics_property_overrides)

    # Write sentinel so Allrun.post can discover the case_id without template substitution.
    (case_root / ".driverfoam_case_id").write_text(case.case_id)



def _run_case(
    case_root: Path,
    setup_root: Path,
    case: CaseConfig,
    *,
    tutorials_root: Path | None = None,
    run_script_relpath: Path = defaults.RUN_SCRIPT_RELPATH,
    parallel: bool = defaults.PARALLEL,
    output_dir: Path,
    cv_extract_script_relpath: Path = defaults.CV_EXTRACT_SCRIPT_RELPATH,
) -> None:
    run_script = resolve_run_script_path(
        tutorials_root=tutorials_root,
        run_script_relpath=run_script_relpath,
    )
    command = [
        "bash",
        "-l",
        str(case_root / "Allclean"),
    ]
    subprocess.run(command, cwd=case_root, check=True)

    command = [
        "bash",
        "-l",
        str(run_script),
        "--case-dir",
        str(case_root),
    ]
    if parallel:
        command.append("parallel")
    subprocess.run(command, check=True)

    # Extract CV directly in Python (like cable_1d_cv_convergence) so the
    # real case_id is available without any shell-level template substitution.
    output_dir.mkdir(parents=True, exist_ok=True)
    extract_module = load_python_module(
        setup_root / cv_extract_script_relpath,
        module_name=f"{defaults.TUTORIAL_NAME}_extract_cv",
    )
    extract_module.write_summary_files(
        case_dir=case_root,
        output_dir=output_dir,
        case_id=case.case_id,
    )

def _postprocess(
    setup_root: Path,
    output_dir: Path,
    *,
    postprocess_script_relpath: Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    table_summary_relpath: Path = defaults.TABLE_SUMMARY_RELPATH,
    strict_artifacts: bool = False,
) -> None:
    tasks = []
    if (setup_root / postprocess_script_relpath).exists():
        tasks.append(
            PostprocessTask(
                module_relpath=postprocess_script_relpath,
                function_name=postprocess_function_name,
            )
        )
    tasks.append(PostprocessTask(module_relpath=table_summary_relpath))

    run_postprocess_tasks(
        setup_root=setup_root,
        output_dir=output_dir,
        tutorial_name=defaults.TUTORIAL_NAME,
        strict_artifacts=strict_artifacts,
        tasks=tasks,
    )





def make_spec(
    *,
    tutorials_root: Path | None = None,
    case_dir_name: str = defaults.CASE_DIR_NAME,
    setup_dir_name: str | None = defaults.SETUP_DIR_NAME,
    output_dir_name: str | None = defaults.OUTPUT_DIR_NAME,
    ionic_models: Sequence[str] = defaults.IONIC_MODELS,
    ionic_model_tissue_map: Mapping[str, Sequence[str]] = defaults.IONIC_MODEL_TISSUE_MAP,
    dt_values: Sequence[float] = defaults.DT_VALUES,
    dx_values: Sequence[float] = defaults.DX_VALUES,
    solvers: Sequence[str] = defaults.SOLVERS,
    conductivity_values: Sequence[str] = defaults.CONDUCTIVITY_VALUES,
    s1_interval_ms: float = defaults.S1_INTERVAL_MS,
    n_s1: int = defaults.N_S1,
    n_s2: int = defaults.N_S2,
    s2_intervals_ms: Sequence[float] = defaults.S2_INTERVALS_MS,
    end_time_buffer_s: float = defaults.END_TIME_BUFFER_S,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    cable_length_mm: float = defaults.CABLE_LENGTH_MM,
    cross_section_cell_counts: Sequence[int] = defaults.CROSS_SECTION_CELL_COUNTS,
    control_dict_relpath: str | Path = defaults.CONTROL_DICT_RELPATH,
    block_mesh_dict_relpath: str | Path = defaults.BLOCK_MESH_DICT_RELPATH,
    electro_properties_relpath: str | Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    physics_properties_relpath: str | Path = "constant/physicsProperties",
    electro_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    physics_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    cv_extract_script_relpath: str | Path = defaults.CV_EXTRACT_SCRIPT_RELPATH,
    postprocess_script_relpath: str | Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    table_summary_relpath: str | Path = defaults.TABLE_SUMMARY_RELPATH,
    parallel: bool = defaults.PARALLEL,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    ionic_models_list = [str(m) for m in ionic_models]
    if not ionic_models_list:
        raise ValueError("ionic_models cannot be empty")

    for ionic_model in ionic_models_list:
        if not ionic_model_tissue_map.get(ionic_model):
            raise KeyError(f"Missing tissue mapping for ionic model '{ionic_model}'")

    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        default_output_dir_name=defaults.OUTPUT_DIR_NAME,
    )

    return TutorialSpec(
        name=defaults.TUTORIAL_NAME,
        case_root=case_root,
        setup_root=setup_root,
        output_dir=output_dir,
        build_cases=partial(
            _build_cases,
            ionic_models=ionic_models_list,
            ionic_model_tissue_map=ionic_model_tissue_map,
            dt_values=dt_values,
            dx_values=dx_values,
            solvers=solvers,
            conductivity_values=conductivity_values,
            s2_intervals_ms=list(s2_intervals_ms),
        ),
        apply_case=partial(
            _apply_case,
            s1_interval_ms=s1_interval_ms,
            n_s1=n_s1,
            n_s2=n_s2,
            end_time_buffer_s=end_time_buffer_s,
            electro_properties_scope=electro_properties_scope,
            control_dict_relpath=Path(control_dict_relpath),
            block_mesh_dict_relpath=Path(block_mesh_dict_relpath),
            electro_properties_relpath=Path(electro_properties_relpath),
            physics_properties_relpath=Path(physics_properties_relpath),
            electro_property_overrides=electro_property_overrides,
            physics_property_overrides=physics_property_overrides,
            cable_length_mm=cable_length_mm,
            cross_section_cell_counts=cross_section_cell_counts,
        ),
        run_case=partial(
            _run_case,
            tutorials_root=tutorials_root,
            run_script_relpath=Path(run_script_relpath),
            parallel=parallel,
            output_dir=output_dir,
            cv_extract_script_relpath=Path(cv_extract_script_relpath),
        ),
        collect_outputs=None,  # let foamctl collect postProcessing natively
        postprocess=partial(
            _postprocess,
            postprocess_script_relpath=Path(postprocess_script_relpath),
            postprocess_function_name=postprocess_function_name,
            table_summary_relpath=Path(table_summary_relpath),
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "python": sys.executable,
            "expected_artifacts": [
                DataArtifact(
                    artifact_id="cable_probes",
                    path_pattern="postProcessing/cableProbes/*/Vm",
                    format="openfoam_probes",
                    description="Voltage probes along the cable.",
                    produced_by="cardiacFoam"
                )
            ],
            "notes": "1D spatial S1-S2 protocol using monodomain natively inside driverFOAM.",
            "workflow_dag": {
                "steps": [
                    {
                        "id": "clean",
                        "command": "Allclean",
                        "args": [],
                    },
                    {
                        "id": "run",
                        "command": "Allrun",
                        "args": ["parallel"] if parallel else [],
                        "depends_on": ["clean"],
                        "produces": ["monodomain_vm_series", "cable_probes"],
                    },
                    {
                        "id": "extract_cv",
                        "command": "Allrun.post",
                        # --case-id is omitted intentionally: Allrun.post reads
                        # the .driverfoam_case_id sentinel written by apply_case.
                        "args": ["--output-dir", str(output_dir)],
                        "depends_on": ["run"],
                    },
                ]
            },
        },
    )
