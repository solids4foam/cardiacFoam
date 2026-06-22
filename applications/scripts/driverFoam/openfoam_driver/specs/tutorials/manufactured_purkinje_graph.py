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

from ...core.defaults import manufactured_purkinje_graph as defaults
from ...core.runtime.models import CaseConfig, TutorialSpec
from ...postprocessing.driver import PostprocessTask, run_postprocess_tasks
from ..common import resolve_spec_paths


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


def _run_case(
    case_root: Path,
    setup_root: Path,
    case: CaseConfig,
    *,
    n_steps: int = defaults.N_STEPS,
    delta_t: float = defaults.DELTA_T,
    block_mesh_dict_relpath: Path = defaults.BLOCK_MESH_DICT_RELPATH,
) -> None:
    del setup_root
    graph_id = str(case.params["graph_id"])
    output_dir = case_root / defaults.OUTPUT_DIR_NAME / graph_id
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    _ensure_mesh(case_root, block_mesh_dict_relpath)

    stale_patterns = (
        "postProcessing/graph_*_nodes.dat",
        "postProcessing/purkinjeNetwork.dat",
    )
    for pattern in stale_patterns:
        for path in case_root.glob(pattern):
            path.unlink()

    with (output_dir / "log.runPurkinjeGraph").open("w") as log:
        subprocess.run(
            [
                "runPurkinjeGraph",
                "-case",
                str(case_root),
                "-nSteps",
                str(int(n_steps)),
                "-deltaT",
                f"{float(delta_t):.12g}",
            ],
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )

    purkinje_dat = case_root / "postProcessing" / "purkinjeNetwork.dat"
    if purkinje_dat.exists():
        shutil.copy2(purkinje_dat, output_dir / purkinje_dat.name)

    graph_outputs = sorted((case_root / "postProcessing").glob("graph_*_nodes.dat"))
    if not graph_outputs:
        raise FileNotFoundError(
            f"No manufactured graph verifier output was written for {graph_id}"
        )
    for source in graph_outputs:
        shutil.copy2(source, output_dir / source.name)

    vtk_dir = case_root / "postProcessing" / "purkinjeNetworkVTK"
    if vtk_dir.exists():
        destination_vtk = output_dir / "purkinjeNetworkVTK"
        shutil.copytree(vtk_dir, destination_vtk, dirs_exist_ok=True)


def _collect_outputs(case_root: Path, output_dir: Path) -> None:
    archived_dir = case_root / defaults.OUTPUT_DIR_NAME
    if not archived_dir.exists() or archived_dir.resolve() == output_dir.resolve():
        return
    if output_dir.exists():
        shutil.rmtree(output_dir)
    shutil.copytree(archived_dir, output_dir)


def _postprocess(
    setup_root: Path,
    output_dir: Path,
    *,
    tutorial_name: str = defaults.TUTORIAL_NAME,
    postprocess_script_relpath: Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    graph_ids: Sequence[str] = defaults.GRAPH_IDS,
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
                kwargs={"graph_ids": list(graph_ids)},
            )
        ],
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
    postprocess_script_relpath: str | Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
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
        run_case=partial(_run_case, n_steps=n_steps, delta_t=delta_t),
        collect_outputs=_collect_outputs,
        postprocess=partial(
            _postprocess,
            tutorial_name=tutorial_name,
            postprocess_script_relpath=Path(postprocess_script_relpath),
            postprocess_function_name=postprocess_function_name,
            graph_ids=graph_ids_list,
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "notes": "Manufactured Purkinje graph convergence benchmark",
            "workflow_dag": {
                "steps": [
                    {"id": "mesh", "command": "blockMesh", "depends_on": []},
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
            "postprocess_script_relpath": str(postprocess_script_relpath),
            "postprocess_function_name": postprocess_function_name,
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
