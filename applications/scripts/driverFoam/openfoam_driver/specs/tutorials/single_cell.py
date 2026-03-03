from __future__ import annotations

import subprocess
import sys
from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path

from ...core.defaults import single_cell as defaults
from ...postprocessing.driver import PostprocessTask, run_postprocess_tasks
from ..common import (
    collect_outputs_by_pattern,
    resolve_run_script_path,
    resolve_spec_paths,
    set_ionic_model,
    set_stimulus_amplitude,
    set_tissue,
)
from ...core.runtime.models import CaseConfig, TutorialSpec


def _build_cases(
    ionic_models: Sequence[str],
    ionic_model_tissue_map: Mapping[str, Sequence[str]],
) -> list[CaseConfig]:
    cases: list[CaseConfig] = []
    for ionic_model in ionic_models:
        tissues = ionic_model_tissue_map.get(ionic_model)
        if not tissues:
            raise KeyError(f"Missing tissue list for ionic model '{ionic_model}'")
        for tissue in tissues:
            case_id = f"{ionic_model}_{tissue}"
            cases.append(
                CaseConfig(
                    case_id=case_id,
                    params={"ionicModel": ionic_model, "tissue": tissue},
                )
            )
    return cases


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    stimulus_map: Mapping[str, float],
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    electro_properties_relpath: Path = defaults.ELECTRO_PROPERTIES_RELPATH,
) -> None:
    tissue = case.params["tissue"]
    ionic_model = case.params["ionicModel"]

    if ionic_model not in stimulus_map:
        raise KeyError(f"Missing stimulus amplitude for ionic model '{ionic_model}'")

    electro_properties_file = case_root / electro_properties_relpath

    set_tissue(electro_properties_file, tissue, scope=electro_properties_scope)
    set_ionic_model(electro_properties_file, ionic_model, scope=electro_properties_scope)
    set_stimulus_amplitude(
        electro_properties_file,
        stimulus_map[ionic_model],
        scope=electro_properties_scope,
    )


def _run_case(
    case_root: Path,
    setup_root: Path,
    _: CaseConfig,
    *,
    tutorials_root: Path | None = None,
    run_script_relpath: Path = defaults.RUN_SCRIPT_RELPATH,
) -> None:
    run_script = resolve_run_script_path(
        tutorials_root=tutorials_root,
        run_script_relpath=run_script_relpath,
    )
    subprocess.run(
        ["bash", "-l", str(run_script), "--case-dir", str(case_root)],
        check=True,
    )


def _collect_outputs(case_root: Path, output_dir: Path, *, output_glob: str = defaults.OUTPUT_GLOB) -> None:
    collect_outputs_by_pattern(case_root, output_dir, pattern=output_glob)


def _postprocess(
    setup_root: Path,
    output_dir: Path,
    *,
    postprocess_script_relpath: Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    strict_artifacts: bool = False,
) -> None:
    run_postprocess_tasks(
        setup_root=setup_root,
        output_dir=output_dir,
        tutorial_name=defaults.TUTORIAL_NAME,
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
    case_dir_name: str = defaults.TUTORIAL_NAME,
    setup_dir_name: str | None = None,
    output_dir_name: str | None = None,
    ionic_models: Sequence[str] = defaults.IONIC_MODELS,
    ionic_model_tissue_map: Mapping[str, Sequence[str]] = defaults.IONIC_MODEL_TISSUE_MAP,
    stimulus_map: Mapping[str, float] = defaults.STIMULUS_MAP,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    electro_properties_relpath: str | Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    output_glob: str = defaults.OUTPUT_GLOB,
    postprocess_script_relpath: str | Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    ionic_models_list = [str(item) for item in ionic_models]
    if not ionic_models_list:
        raise ValueError("ionic_models cannot be empty")

    for ionic_model in ionic_models_list:
        tissues = ionic_model_tissue_map.get(ionic_model)
        if not tissues:
            raise KeyError(f"Missing tissue mapping for ionic model '{ionic_model}'")
        if ionic_model not in stimulus_map:
            raise KeyError(f"Missing stimulus amplitude for ionic model '{ionic_model}'")

    electro_properties_path = Path(electro_properties_relpath)
    run_script_path = Path(run_script_relpath)
    postprocess_script_path = Path(postprocess_script_relpath)

    default_output_dir_name = defaults.OUTPUT_DIR_NAME
    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        default_output_dir_name=default_output_dir_name,
    )

    return TutorialSpec(
        name=case_dir_name,
        case_root=case_root,
        setup_root=setup_root,
        output_dir=output_dir,
        build_cases=partial(
            _build_cases,
            ionic_models=ionic_models_list,
            ionic_model_tissue_map=ionic_model_tissue_map,
        ),
        apply_case=partial(
            _apply_case,
            stimulus_map=stimulus_map,
            electro_properties_scope=electro_properties_scope,
            electro_properties_relpath=electro_properties_path,
        ),
        run_case=partial(
            _run_case,
            tutorials_root=tutorials_root,
            run_script_relpath=run_script_path,
        ),
        collect_outputs=partial(_collect_outputs, output_glob=output_glob),
        postprocess=partial(
            _postprocess,
            postprocess_script_relpath=postprocess_script_path,
            postprocess_function_name=postprocess_function_name,
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "python": sys.executable,
            "notes": "Single-cell sweep on ionic model and tissue types.",
            "ionic_models": ionic_models_list,
            "electro_properties_relpath": str(electro_properties_path),
            "electro_properties_scope": electro_properties_scope,
            "run_script_relpath": str(run_script_path),
            "output_glob": output_glob,
            "postprocess_script_relpath": str(postprocess_script_path),
            "postprocess_function_name": postprocess_function_name,
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
