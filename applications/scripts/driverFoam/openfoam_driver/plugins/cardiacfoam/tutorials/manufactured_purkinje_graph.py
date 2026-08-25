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
#     manufactured_purkinje_graph
#
# Description
#     Defines the manufactured Purkinje graph convergence tutorial.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import shutil
import subprocess
from functools import partial
from pathlib import Path
from typing import Sequence

from openfoam_driver.plugins.cardiacfoam.tutorials.defaults import manufactured_purkinje_graph as defaults
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec
from openfoam_driver.specs.common import resolve_spec_paths


def _build_cases(graph_ids: Sequence[str]) -> list[CaseConfig]:
    return [
        CaseConfig(
            case_id=str(graph_id),
            params={"graph_id": str(graph_id)},
        )
        for graph_id in graph_ids
    ]


def _apply_case(case_root: Path, case: CaseConfig) -> None:
    graph_id = str(case.params["graph_id"])
    source = case_root / "constant" / f"purkinjeGraph.{graph_id}"
    destination = case_root / "constant" / "purkinjeGraph"
    if not source.exists():
        raise FileNotFoundError(f"Missing graph file: {source}")
    shutil.copy2(source, destination)


def _ensure_mesh(case_root: Path, block_mesh_dict_relpath: Path) -> None:
    if (case_root / "constant" / "polyMesh").exists():
        return
    with (case_root / "log.blockMesh").open("w") as log:
        subprocess.run(
            [
                "blockMesh",
                "-case",
                str(case_root),
                "-dict",
                str(case_root / block_mesh_dict_relpath),
            ],
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )


def make_spec(
    *,
    tutorials_root: Path | None = None,
    tutorial_name: str = defaults.TUTORIAL_NAME,
    case_dir_name: str = defaults.CASE_DIR_NAME,
    setup_dir_name: str = defaults.SETUP_DIR_NAME,
    output_dir_name: str = defaults.OUTPUT_DIR_NAME,
    graph_ids: Sequence[str] = defaults.GRAPH_IDS,
    n_steps: int = defaults.N_STEPS,
    delta_t: float = defaults.DELTA_T,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        default_output_dir_name=defaults.OUTPUT_DIR_NAME,
    )
    graph_ids_list = [str(item) for item in graph_ids]

    return TutorialSpec(
        name=tutorial_name,
        case_root=case_root,
        setup_root=setup_root,
        output_dir=output_dir,
        build_cases=partial(_build_cases, graph_ids=graph_ids_list),
        apply_case=_apply_case,
        metadata={
            "notes": "Manufactured Purkinje graph convergence benchmark",
            "workflow_dag": {
                "steps": [
                    {
                        "id": "mesh",
                        "command": "blockMesh",
                        "args": ["-dict", "system/blockMeshDict.3D"],
                        "depends_on": [],
                    },
                    {
                        "id": "solve",
                        "command": "runPurkinjeGraph",
                        "depends_on": ["mesh"],
                    },
                ]
            },
            "graph_ids": graph_ids_list,
            "n_steps": int(n_steps),
            "delta_t": float(delta_t),
            "control_dict_relpath": str(defaults.CONTROL_DICT_RELPATH),
            "electro_properties_relpath": str(defaults.ELECTRO_PROPERTIES_RELPATH),
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
