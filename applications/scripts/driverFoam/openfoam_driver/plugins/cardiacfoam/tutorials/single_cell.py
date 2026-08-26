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
#     single_cell
#
# Description
#     Defines configuration template for single-cell scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import subprocess
import sys
from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path

from openfoam_driver.plugins.cardiacfoam.tutorials.defaults import single_cell as defaults
from openfoam_driver.plugins.cardiacfoam.overrides import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
)
from openfoam_driver.core.specs.common import (
    resolve_run_script_path,
    resolve_spec_paths,
    set_end_time,
)
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec


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
    control_dict_relpath: Path = Path("system/controlDict"),
    physics_properties_relpath: Path = Path("constant/physicsProperties"),
    stim_start_ms: float | None = None,
    s1_interval_ms: float | None = None,
    n_s1: int | None = None,
    end_time_buffer_s: float = 0.0,
    write_after_time_s: float | None = None,
    ionic_export: Sequence[str] | None = None,
    electro_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    physics_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
) -> None:
    tissue = case.params["tissue"]
    ionic_model = case.params["ionicModel"]

    if ionic_model not in stimulus_map:
        raise KeyError(f"Missing stimulus amplitude for ionic model '{ionic_model}'")

    electro_properties_file = case_root / electro_properties_relpath
    control_dict_file = case_root / control_dict_relpath
    physics_properties_file = case_root / physics_properties_relpath
    case_overrides = {
        f"{electro_properties_scope}.tissue": tissue,
        f"{electro_properties_scope}.ionicModel": ionic_model,
        f"{electro_properties_scope}.singleCellStimulus.stim_amplitude": stimulus_map[ionic_model],
    }

    # Optional pacing controls are deliberately opt-in so the established
    # all-model singleCell sweep retains its authored protocol unchanged.
    stimulus_overrides = {}
    if stim_start_ms is not None:
        stimulus_overrides[
            f"{electro_properties_scope}.singleCellStimulus.stim_start"
        ] = stim_start_ms
    if s1_interval_ms is not None:
        stimulus_overrides[
            f"{electro_properties_scope}.singleCellStimulus.stim_period_S1"
        ] = s1_interval_ms
    if n_s1 is not None:
        stimulus_overrides[
            f"{electro_properties_scope}.singleCellStimulus.nstim1"
        ] = n_s1

    case_overrides.update(stimulus_overrides)

    if s1_interval_ms is not None or n_s1 is not None or stim_start_ms is not None:
        start_ms = 0.0 if stim_start_ms is None else float(stim_start_ms)
        interval_ms = 0.0 if s1_interval_ms is None else float(s1_interval_ms)
        beats = 1 if n_s1 is None else int(n_s1)
        if interval_ms <= 0.0 or beats <= 0:
            raise ValueError("comparison pacing requires a positive CL and n_s1")
        set_end_time(
            control_dict_file,
            (start_ms + interval_ms * beats) / 1000.0 + end_time_buffer_s,
        )

    if write_after_time_s is not None:
        case_overrides[f"{electro_properties_scope}.writeAfterTime"] = write_after_time_s
    if ionic_export is not None:
        case_overrides[
            f"{electro_properties_scope}.outputVariables.ionic.export"
        ] = "(" + " ".join(str(name) for name in ionic_export) + ")"

    apply_electro_property_overrides(electro_properties_file, case_overrides)
    apply_electro_property_overrides(electro_properties_file, electro_property_overrides)
    apply_physics_property_overrides(physics_properties_file, physics_property_overrides)


def make_spec(
    *,
    tutorials_root: Path | None = None,
    case_dir_name: str = defaults.CASE_DIR_NAME,
    setup_dir_name: str | None = defaults.SETUP_DIR_NAME,
    output_dir_name: str | None = None,
    ionic_models: Sequence[str] = defaults.IONIC_MODELS,
    ionic_model_tissue_map: Mapping[str, Sequence[str]] = defaults.IONIC_MODEL_TISSUE_MAP,
    ionic_model: str | None = None,
    tissue: str | None = None,
    stimulus_map: Mapping[str, float] = defaults.STIMULUS_MAP,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    electro_properties_relpath: str | Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    physics_properties_relpath: str | Path = "constant/physicsProperties",
    control_dict_relpath: str | Path = "system/controlDict",
    stim_start_ms: float | None = None,
    s1_interval_ms: float | None = None,
    n_s1: int | None = None,
    end_time_buffer_s: float = 0.0,
    write_after_time_s: float | None = None,
    ionic_export: Sequence[str] | None = None,
    electro_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    physics_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    output_glob: str = defaults.OUTPUT_GLOB,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    if (ionic_model is None) != (tissue is None):
        raise ValueError(
            "ionic_model and tissue must be given together (a single-case "
            "override, for entry-mode sweeps that need build_cases() to "
            "collapse to exactly one case) or not at all (the default "
            "full-catalog ionic_models/ionic_model_tissue_map sweep)"
        )
    if ionic_model is not None:
        ionic_models = [ionic_model]
        ionic_model_tissue_map = {ionic_model: [tissue]}

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
    physics_properties_path = Path(physics_properties_relpath)
    control_dict_path = Path(control_dict_relpath)
    run_script_path = Path(run_script_relpath)

    default_output_dir_name = defaults.OUTPUT_DIR_NAME
    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        default_output_dir_name=default_output_dir_name,
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
        ),
        apply_case=partial(
            _apply_case,
            stimulus_map=stimulus_map,
            electro_properties_scope=electro_properties_scope,
            electro_properties_relpath=electro_properties_path,
            control_dict_relpath=control_dict_path,
            physics_properties_relpath=physics_properties_path,
            stim_start_ms=stim_start_ms,
            s1_interval_ms=s1_interval_ms,
            n_s1=n_s1,
            end_time_buffer_s=end_time_buffer_s,
            write_after_time_s=write_after_time_s,
            ionic_export=ionic_export,
            electro_property_overrides=electro_property_overrides,
            physics_property_overrides=physics_property_overrides,
        ),
        metadata={
            "python": sys.executable,
            "notes": "Single-cell sweep on ionic model and tissue types.",
            "workflow_dag": {
                "steps": [
                    {"id": "mesh", "command": "blockMesh", "depends_on": []},
                    {"id": "solve", "command": "cardiacFoam", "depends_on": ["mesh"]},
                ]
            },
            "ionic_models": ionic_models_list,
            "electro_properties_relpath": str(electro_properties_path),
            "physics_properties_relpath": str(physics_properties_path),
            "control_dict_relpath": str(control_dict_path),
            "stim_start_ms": stim_start_ms,
            "s1_interval_ms": s1_interval_ms,
            "n_s1": n_s1,
            "end_time_buffer_s": end_time_buffer_s,
            "write_after_time_s": write_after_time_s,
            "ionic_export": list(ionic_export) if ionic_export is not None else None,
            "electro_properties_scope": electro_properties_scope,
            "run_script_relpath": str(run_script_path),
            "output_glob": output_glob,
            "has_electro_property_overrides": bool(electro_property_overrides),
            "has_physics_property_overrides": bool(physics_property_overrides),
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
