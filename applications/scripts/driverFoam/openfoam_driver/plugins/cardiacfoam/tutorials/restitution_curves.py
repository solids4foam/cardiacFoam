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
#     restitution_curves
#
# Description
#     Defines configuration template for restitution curve scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import shutil
import subprocess
import sys
from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path

from openfoam_driver.plugins.cardiacfoam.tutorials.defaults import restitution_curves as defaults
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
    s2_intervals_ms: Sequence[int],
) -> list[CaseConfig]:
    cases: list[CaseConfig] = []
    for ionic_model in ionic_models:
        tissues = ionic_model_tissue_map.get(ionic_model)
        if not tissues:
            raise KeyError(f"Missing tissue list for ionic model '{ionic_model}'")
        for tissue in tissues:
            for s2 in s2_intervals_ms:
                case_id = f"{ionic_model}_{tissue}_S2_{s2}"
                cases.append(
                    CaseConfig(
                        case_id=case_id,
                        params={
                            "ionicModel": ionic_model,
                            "tissue": tissue,
                            "s2Interval": s2,
                        },
                    )
                )
    return cases


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    stimulus_map: Mapping[str, float],
    s1_interval_ms: int,
    n_s1: int,
    n_s2: int,
    write_after_time_s: float,
    end_time_buffer_s: float,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    electro_properties_relpath: Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    control_dict_relpath: Path = defaults.CONTROL_DICT_RELPATH,
    physics_properties_relpath: Path = Path("constant/physicsProperties"),
    electro_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    physics_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
) -> None:
    ionic_model = case.params["ionicModel"]
    tissue = case.params["tissue"]
    s2_interval_ms = case.params["s2Interval"]

    if ionic_model not in stimulus_map:
        raise KeyError(f"Missing stimulus amplitude for ionic model '{ionic_model}'")

    electro_properties_file = case_root / electro_properties_relpath
    control_dict_file = case_root / control_dict_relpath
    physics_properties_file = case_root / physics_properties_relpath
    case_overrides = {
        f"{electro_properties_scope}.tissue": tissue,
        f"{electro_properties_scope}.ionicModel": ionic_model,
        f"{electro_properties_scope}.singleCellStimulus.stim_amplitude": stimulus_map[ionic_model],
        f"{electro_properties_scope}.singleCellStimulus.stim_period_S1": s1_interval_ms,
        f"{electro_properties_scope}.singleCellStimulus.nstim1": n_s1,
        f"{electro_properties_scope}.singleCellStimulus.stim_period_S2": s2_interval_ms,
        f"{electro_properties_scope}.singleCellStimulus.nstim2": n_s2,
        f"{electro_properties_scope}.writeAfterTime": write_after_time_s,
    }

    # controlDict: endTime = (S1*n_S1 + S2*n_S2) / 1000 + buffer
    end_time = (s1_interval_ms * n_s1 + s2_interval_ms * n_s2) / 1000.0 + end_time_buffer_s
    set_end_time(control_dict_file, end_time)
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
    stimulus_map: Mapping[str, float] = defaults.STIMULUS_MAP,
    s1_interval_ms: int = defaults.S1_INTERVAL_MS,
    n_s1: int = defaults.N_S1,
    n_s2: int = defaults.N_S2,
    s2_intervals_ms: Sequence[int] = defaults.S2_INTERVALS_MS,
    end_time_buffer_s: float = defaults.END_TIME_BUFFER_S,
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    electro_properties_relpath: str | Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    control_dict_relpath: str | Path = defaults.CONTROL_DICT_RELPATH,
    physics_properties_relpath: str | Path = "constant/physicsProperties",
    electro_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    physics_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    output_glob: str = defaults.OUTPUT_GLOB,
    show_plots: bool = False,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    ionic_models_list = [str(m) for m in ionic_models]
    if not ionic_models_list:
        raise ValueError("ionic_models cannot be empty")

    for ionic_model in ionic_models_list:
        if not ionic_model_tissue_map.get(ionic_model):
            raise KeyError(f"Missing tissue mapping for ionic model '{ionic_model}'")
        if ionic_model not in stimulus_map:
            raise KeyError(f"Missing stimulus amplitude for ionic model '{ionic_model}'")

    if not s2_intervals_ms:
        raise ValueError("s2_intervals_ms cannot be empty")

    electro_properties_path = Path(electro_properties_relpath)
    control_dict_path = Path(control_dict_relpath)
    physics_properties_path = Path(physics_properties_relpath)
    run_script_path = Path(run_script_relpath)

    # write_after_time: start writing 2 s before the end of the S1 phase
    write_after_time_s = (s1_interval_ms * n_s1) / 1000.0 - 2.0

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
            s2_intervals_ms=list(s2_intervals_ms),
        ),
        apply_case=partial(
            _apply_case,
            stimulus_map=stimulus_map,
            s1_interval_ms=s1_interval_ms,
            n_s1=n_s1,
            n_s2=n_s2,
            write_after_time_s=write_after_time_s,
            end_time_buffer_s=end_time_buffer_s,
            electro_properties_scope=electro_properties_scope,
            electro_properties_relpath=electro_properties_path,
            control_dict_relpath=control_dict_path,
            physics_properties_relpath=physics_properties_path,
            electro_property_overrides=electro_property_overrides,
            physics_property_overrides=physics_property_overrides,
        ),
        metadata={
            "python": sys.executable,
            "notes": "S1–S2 restitution protocol sweep on ionic model, tissue, and S2 interval.",
            "workflow_dag": {
                "steps": [
                    {"id": "solve", "command": "cardiacFoam", "depends_on": []},
                ]
            },
            "ionic_models": ionic_models_list,
            "s1_interval_ms": s1_interval_ms,
            "n_s1": n_s1,
            "n_s2": n_s2,
            "s2_intervals_ms": list(s2_intervals_ms),
            "electro_properties_relpath": str(electro_properties_path),
            "electro_properties_scope": electro_properties_scope,
            "control_dict_relpath": str(control_dict_path),
            "physics_properties_relpath": str(physics_properties_path),
            "run_script_relpath": str(run_script_path),
            "output_glob": output_glob,
            "show_plots": show_plots,
            "has_electro_property_overrides": bool(electro_property_overrides),
            "has_physics_property_overrides": bool(physics_property_overrides),
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
