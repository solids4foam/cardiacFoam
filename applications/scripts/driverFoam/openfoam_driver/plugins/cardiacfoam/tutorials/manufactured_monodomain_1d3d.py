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
#     manufactured_monodomain_1d3d
#
# Description
#     Defines the manufactured 1D-3D coupled monodomain tutorial.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Sequence
from functools import partial
from pathlib import Path
from itertools import product

from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec
from openfoam_driver.postprocessing.driver import PostprocessTask, run_postprocess_tasks
from openfoam_driver.plugins.cardiacfoam.overrides import (
    apply_electro_property_overrides,
)
from openfoam_driver.specs.common import (
    replace_block_mesh_resolutions,
    resolve_spec_paths,
    set_delta_t,
)
from openfoam_driver.core.runtime.mutators import update_foam_entry


def _build_cases(
    graph_ids: Sequence[str],
    number_cells: Sequence[int],
    dt_values: Sequence[float],
) -> list[CaseConfig]:
    # We expect these three sequences to be zipped (same length)
    cases: list[CaseConfig] = []
    for graph_id, cells, dt in zip(graph_ids, number_cells, dt_values):
        case_id = f"{cells}"
        cases.append(
            CaseConfig(
                case_id=case_id,
                params={
                    "graph_id": str(graph_id),
                    "cells": int(cells),
                    "dt": float(dt),
                },
            )
        )
    return cases


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    electro_property_overrides: dict[str, object] | None = None,
    end_time: float | None = None,
) -> None:
    graph_id = str(case.params["graph_id"])
    cells = int(case.params["cells"])
    dt_value = float(case.params["dt"])

    # 1. 3D blockMesh resolution
    block_mesh_dict = case_root / "system" / "blockMeshDict.3D"
    block_mesh_active = case_root / "system" / "blockMeshDict.3D.active"
    block_mesh_active.write_text(block_mesh_dict.read_text())
    replace_block_mesh_resolutions(block_mesh_active, f"{cells} {cells} {cells}")

    # 2. 1D graph selection
    source_graph = case_root / "constant" / f"purkinjeGraph.{graph_id}"
    destination_graph = case_root / "constant" / "purkinjeGraph"
    if not source_graph.exists():
        raise FileNotFoundError(f"Missing graph file: {source_graph}")
    shutil.copy2(source_graph, destination_graph)

    # 3. ControlDict deltaT and endTime
    control_dict = case_root / "system" / "controlDict"
    set_delta_t(control_dict, dt_value)
    if end_time is not None:
        update_foam_entry(control_dict, "endTime", end_time)

    # 4. ElectroProperties overrides (rPvj, couplingMode, etc)
    if electro_property_overrides:
        electro_properties = case_root / "constant" / "electroProperties"
        apply_electro_property_overrides(electro_properties, electro_property_overrides)


def _run_case(
    case_root: Path,
    setup_root: Path,
    case: CaseConfig,
) -> None:
    del setup_root

    # Clean previous dynamic outputs before run
    stale_patterns = (
        "postProcessing/graph_*_nodes.dat",
        "postProcessing/3D_*_cells_*.dat",
        "verification/coupled1D3DMonodomain_diagnostics.csv",
    )
    for pattern in stale_patterns:
        for path in case_root.glob(pattern):
            path.unlink()
    # The actual execution is handled by the generic executor running the workflow_dag


def _collect_outputs(case_root: Path, output_dir: Path, case: CaseConfig) -> None:
    cells = int(case.params["cells"])
    # The case output directory
    case_output = output_dir / str(cells)
    if case_output.exists():
        shutil.rmtree(case_output)
    case_output.mkdir(parents=True, exist_ok=True)

    for source in case_root.glob("postProcessing/graph_*_nodes.dat"):
        shutil.copy2(source, case_output / source.name)
    for source in case_root.glob("postProcessing/3D_*_cells_*.dat"):
        shutil.copy2(source, case_output / source.name)
    diag = case_root / "verification" / "coupled1D3DMonodomain_diagnostics.csv"
    if diag.exists():
        shutil.copy2(diag, case_output / "coupling_diagnostics.csv")


def _postprocess(
    setup_root: Path,
    output_dir: Path,
    *,
    strict_artifacts: bool = False,
) -> None:
    # We use a subprocess directly because the original post-processing script
    # expects `--output-dir` and does not export a single entrypoint function
    # that conforms to driverFoam's PostprocessTask kwargs convention easily.
    import subprocess
    post_proc_script = setup_root / "post_processing_coupled_1D3D.py"
    if not post_proc_script.exists():
        print(f"WARN: post-processing script not found at {post_proc_script}")
        return
    
    subprocess.run(
        ["python3", str(post_proc_script), "--output-dir", str(output_dir)],
        check=True
    )


def make_spec(
    *,
    tutorials_root: Path | None = None,
    tutorial_name: str = "manufacturedMonodomain1D3D",
    case_dir_name: str = "manufacturedSolutions/monodomain1D3D",
    setup_dir_name: str = "setup",
    output_dir_name: str = "setup/studies/coupledConvergence/results",
    graph_ids: Sequence[str] = ("nodes011", "nodes021", "nodes041", "nodes081"),
    number_cells: Sequence[int] = (10, 20, 40, 80),
    dt_values: Sequence[float] = (1.40174e-04 * (80/10)**2, 1.40174e-04 * (80/20)**2, 1.40174e-04 * (80/40)**2, 1.40174e-04),
    end_time: float = 0.1,
    electro_property_overrides: dict[str, object] | None = None,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        default_output_dir_name=output_dir_name,
    )
    graph_ids_list = [str(item) for item in graph_ids]
    cells_list = [int(item) for item in number_cells]
    dt_values_list = [float(item) for item in dt_values]

    return TutorialSpec(
        name=tutorial_name,
        case_root=case_root,
        setup_root=setup_root,
        output_dir=output_dir,
        build_cases=partial(
            _build_cases,
            graph_ids=graph_ids_list,
            number_cells=cells_list,
            dt_values=dt_values_list,
        ),
        apply_case=partial(
            _apply_case,
            electro_property_overrides=electro_property_overrides,
            end_time=end_time,
        ),
        run_case=_run_case,
        collect_outputs=lambda c_root, o_dir, case: _collect_outputs(c_root, o_dir, case),
        postprocess=partial(
            _postprocess,
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "notes": "Manufactured coupled 1D-3D monodomain convergence benchmark",
            "workflow_dag": {
                "steps": [
                    {
                        "id": "mesh",
                        "command": "blockMesh",
                        "args": ["-dict", "system/blockMeshDict.3D.active"],
                        "depends_on": [],
                    },
                    {
                        "id": "solve",
                        "command": "cardiacFoam",
                        "depends_on": ["mesh"],
                    },
                ]
            },
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
