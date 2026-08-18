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
#     niederer_2012
#
# Description
#     Defines configuration template for Niederer 2012 benchmarks.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from functools import partial
from itertools import product
from pathlib import Path

from openfoam_driver.plugins.cardiacfoam.tutorials.defaults import niederer_2012 as defaults
from openfoam_driver.plugins.cardiacfoam.overrides import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
)
from openfoam_driver.specs.common import (
    replace_block_mesh_resolutions,
    resolve_spec_paths,
    set_delta_t,
)
from openfoam_driver.specs.mesh_provisioning import cell_counts_from_dx
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec


def _closest_key(mapping: dict[float, object], value: float) -> float:
    return min(mapping.keys(), key=lambda item: abs(item - value))



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


def _workflow_dag_for(mesh_family: str) -> dict[str, object]:
    if mesh_family == "tet":
        steps = [
            {"id": "clean", "command": "Allclean", "depends_on": []},
            {
                "id": "gmsh",
                "command": "gmsh",
                "args": ["-3", "setup/studies/tetConvergence/slab.geo", "-o", "slab.msh", "-format", "msh2"],
                "depends_on": ["clean"],
            },
            {
                "id": "gmshToFoam",
                "command": "gmshToFoam",
                "args": ["slab.msh"],
                "depends_on": ["gmsh"],
            },
            {"id": "checkMesh", "command": "checkMesh", "depends_on": ["gmshToFoam"]},
            {"id": "solve", "command": "cardiacFoam", "depends_on": ["checkMesh"]},
        ]
    else:
        steps = [
            {"id": "mesh", "command": "blockMesh", "depends_on": []},
            {"id": "solve", "command": "cardiacFoam", "depends_on": ["mesh"]},
        ]

    steps.extend([
        {
            "id": "samplePoints",
            "command": "postProcess -func Niedererpoints -latestTime",
            "depends_on": ["solve"],
        },
        {
            "id": "sampleLines",
            "command": "postProcess -func Niedererlines -latestTime",
            "depends_on": ["solve"],
        },
    ])
    return {"steps": steps}


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    mesh_family: str = "hex",
    tet_geo_template_relpath: Path = Path("setup/studies/tetConvergence/slab.geo.template"),
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    control_dict_relpath: Path = defaults.CONTROL_DICT_RELPATH,
    block_mesh_dict_relpath: Path = defaults.BLOCK_MESH_DICT_RELPATH,
    electro_properties_relpath: Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    physics_properties_relpath: Path = Path("constant/physicsProperties"),
    electro_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    physics_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    slab_size_mm: Sequence[float] = defaults.SLAB_SIZE_MM,
    end_time_by_dx: Mapping[float, float] = defaults.END_TIME_BY_DX,
) -> None:
    control_dict = case_root / control_dict_relpath
    block_mesh_dict = case_root / block_mesh_dict_relpath
    electro_properties = case_root / electro_properties_relpath
    physics_properties = case_root / physics_properties_relpath

    dx_mm = float(case.params["dx_mm"])
    dt_ms = float(case.params["dt_ms"])
    tissue = str(case.params["tissue"])
    ionic_model = str(case.params["ionicModel"])
    solver = str(case.params["solver"])
    case_overrides = {
        f"{electro_properties_scope}.tissue": tissue,
        f"{electro_properties_scope}.ionicModel": ionic_model,
        f"{electro_properties_scope}.solutionAlgorithm": solver,
    }

    if mesh_family == "tet":
        template_file = case_root / tet_geo_template_relpath
        if not template_file.exists():
            raise FileNotFoundError(f"Missing tet geo template: {template_file}")
        lc_m = dx_mm * 1e-3
        rendered = template_file.read_text().replace("__LC__", str(lc_m))
        target_file = case_root / "setup" / "studies" / "tetConvergence" / "slab.geo"
        target_file.write_text(rendered)
    else:
        axis_cell_counts = [str(count) for count in cell_counts_from_dx(dx_mm, slab_size_mm)]
        replace_block_mesh_resolutions(block_mesh_dict, " ".join(axis_cell_counts))
    # Input dt is provided in milliseconds in the JSON/spec settings.
    set_delta_t(control_dict, dt_ms * 1.0e-3)
    _update_end_time(control_dict, dx_mm, end_time_by_dx=end_time_by_dx)
    apply_electro_property_overrides(electro_properties, case_overrides)
    apply_electro_property_overrides(electro_properties, electro_property_overrides)
    apply_physics_property_overrides(physics_properties, physics_property_overrides)


def _run_case(case_root: Path, setup_root: Path, case: CaseConfig) -> None:
    """Never called. TutorialSpec.run_case is a required dataclass field, but
    nothing in the current execution engine (run --strict / sweep-run) reads
    it -- both drive a case entirely through workflow_dag's own DAG steps
    (mesh/solve/samplePoints/sampleLines, see _workflow_dag_for), executed by
    core.runtime.workflow_orchestrator.run_workflow. Confirmed: zero
    references to `.run_case` anywhere in core/runtime/*.py or cli.py.

    This used to do real work -- launch the case via a bash wrapper script,
    then convert OpenFOAM's raw probe samples into labeled CSVs with a
    hardcoded probe-label list. All of that duplicated, in dead code, what
    the DAG's own samplePoints/sampleLines steps plus
    setup/convert_raw_samples.py now do for real: the DAG samples the field
    natively; convert_raw_samples.py reads probe labels and coordinates
    directly from the raw sample file's own header instead of a second,
    hand-maintained copy of them. See 9e98e71b.
    """
    del case_root, setup_root, case


def make_spec(
    *,
    tutorials_root: Path | None = None,
    case_dir_name: str = defaults.CASE_DIR_NAME,
    setup_dir_name: str | None = defaults.SETUP_DIR_NAME,
    output_dir_name: str = defaults.OUTPUT_DIR_NAME,
    ionic_models: Sequence[str] = defaults.IONIC_MODELS,
    ionic_model_tissue_map: Mapping[str, Sequence[str]] = defaults.IONIC_MODEL_TISSUE_MAP,
    dt_values: Sequence[float] = defaults.DT_VALUES,
    dx_values: Sequence[float] = defaults.DX_VALUES,
    solvers: Sequence[str] = defaults.SOLVERS,
    mesh_family: str = "hex",
    tet_geo_template_relpath: str | Path = "setup/studies/tetConvergence/slab.geo.template",
    electro_properties_scope: str = defaults.ELECTRO_PROPERTIES_SCOPE,
    slab_size_mm: Sequence[float] = defaults.SLAB_SIZE_MM,
    end_time_by_dx: Mapping[float, float] = defaults.END_TIME_BY_DX,
    control_dict_relpath: str | Path = defaults.CONTROL_DICT_RELPATH,
    block_mesh_dict_relpath: str | Path = defaults.BLOCK_MESH_DICT_RELPATH,
    electro_properties_relpath: str | Path = defaults.ELECTRO_PROPERTIES_RELPATH,
    physics_properties_relpath: str | Path = "constant/physicsProperties",
    electro_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
    physics_property_overrides: Mapping[str, object] | Sequence[Mapping[str, object]] | None = None,
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
    physics_properties_path = Path(physics_properties_relpath)
    tet_geo_template_path = Path(tet_geo_template_relpath)

    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
    )

    return TutorialSpec(
        name=defaults.TUTORIAL_NAME,
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
            mesh_family=mesh_family,
            tet_geo_template_relpath=tet_geo_template_path,
            electro_properties_scope=electro_properties_scope,
            control_dict_relpath=control_dict_path,
            block_mesh_dict_relpath=block_mesh_dict_path,
            electro_properties_relpath=electro_properties_path,
            physics_properties_relpath=physics_properties_path,
            electro_property_overrides=electro_property_overrides,
            physics_property_overrides=physics_property_overrides,
            slab_size_mm=slab_size_mm_list,
            end_time_by_dx=end_time_by_dx_map,
        ),
        run_case=_run_case,
        collect_outputs=None,
        metadata={
            "notes": (
                "Niederer Et Al. 2012 slab benchmark sweep "
                "with OpenFOAM functionObject sampling."
            ),
            "workflow_dag": _workflow_dag_for(mesh_family),
            "mesh_family": mesh_family,
            "tet_geo_template_relpath": str(tet_geo_template_path),
            "dx_values": dx_values_list,
            "dt_values": dt_values_list,
            "slab_size_mm": slab_size_mm_list,
            "ionic_model_tissue_map": ionic_model_tissue_map_normalized,
            "ionic_models": ionic_models_list,
            "solvers": solvers_list,
            "control_dict_relpath": str(control_dict_path),
            "block_mesh_dict_relpath": str(block_mesh_dict_path),
            "electro_properties_relpath": str(electro_properties_path),
            "physics_properties_relpath": str(physics_properties_path),
            "electro_properties_scope": electro_properties_scope,
            "has_electro_property_overrides": bool(electro_property_overrides),
            "has_physics_property_overrides": bool(physics_property_overrides),
        },
    )
