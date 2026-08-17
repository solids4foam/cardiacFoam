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
#     manufactured_bath_bidomain
#
# Description
#     Defines configuration template for manufactured bath bidomain scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path

from openfoam_driver.plugins.cardiacfoam.tutorials.defaults import manufactured_bath_bidomain as defaults
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec
from openfoam_driver.core.runtime.mutators import update_foam_entry
from openfoam_driver.core.runtime.parallel_execution import solve_steps
from openfoam_driver.postprocessing.driver import PostprocessTask, run_postprocess_tasks
from openfoam_driver.plugins.cardiacfoam.overrides import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
    ensure_electro_property_dict,
    remove_electro_property_dict,
)
from openfoam_driver.specs.common import (
    load_python_module,
    replace_block_mesh_resolutions,
    resolve_run_script_path,
    resolve_spec_paths,
    set_delta_t,
    set_end_time,
)
from openfoam_driver.specs.utils import (
    archive_case_logs,
    stage_post_processing_outputs,
)
from openfoam_driver.specs.tet_mesh_provisioning import render_tet_geo
from .manufactured_monodomain_pseudo_ecg import _build_cases


_GRAD_SCHEME_TOKENS: dict[str, str] = {
    "gauss_linear": "Gauss linear",
    "least_squares": "leastSquares",
}




def _phi_e_ref_point_yz(dimension: str, cells: int) -> tuple[float, float]:
    """y/z point for electrodePair's floating-reference cell.

    electrodePair has no ground Dirichlet patch, so phiE floats and
    extracellularPotentialDomain::referenceCell() pins it via whichever
    single mesh cell contains phiERefPoint -- that cell must be owned by
    exactly one processor partition (fatal error otherwise). The exact
    domain midpoint sits exactly on a cell FACE whenever the corresponding
    direction is subdivided by `cells` (every cells value this tutorial
    sweeps is even), which can straddle a processor-decomposition boundary
    for some resolutions (confirmed failing for 2D/40, 3D/20, 3D/40:
    "2 partitions reported a containing cell") while working for others,
    since which decomposition boundaries land where is resolution- and
    method-dependent. Shifting by half a cell width off the midpoint lands
    inside one specific cell's interior instead, which by construction
    belongs to exactly one partition under any decomposition -- robust for
    every resolution, not just the ones tested so far.
    """
    try:
        y_extent, z_extent = defaults.DOMAIN_YZ_EXTENT_BY_DIMENSION[dimension]
        y_subdivided, z_subdivided = defaults.YZ_SUBDIVIDED_BY_DIMENSION[dimension]
    except KeyError as exc:
        raise ValueError(
            f"No phiERefPoint geometry known for dimension {dimension!r}; "
            f"expected one of {sorted(defaults.DOMAIN_YZ_EXTENT_BY_DIMENSION)}."
        ) from exc
    ref_y = y_extent / 2 + (y_extent / (2 * cells) if y_subdivided else 0.0)
    ref_z = z_extent / 2 + (z_extent / (2 * cells) if z_subdivided else 0.0)
    return ref_y, ref_z


_TET_NUMERICS_PROFILES: dict[str, tuple[str, ...]] = {
    "bath_bidomain_tet": ("fvSchemes",),
}


def _case_output_filename(case: CaseConfig) -> str:
    dimension = str(case.params["dimension"])
    cells = int(case.params["cells"])
    solver = str(case.params["solver"])
    return f"bathBidomain_{dimension}_{cells}_cells_{solver}.dat"


_DEFAULT_ECG_DOMAINS_BLOCK = """    ecgDomains
    {
        electrodePositions
        {
            E_left    (-1 0.5 0.5);
            E_right   (2 0.5 0.5);
            E_side    (1.5 1 0.5);
        }

        bodyECG
        {
            ecgSolver    torsoECG;
            ecgVerificationModel    bathECGManufacturedVerifier;
            reportElectrodeLookup    yes;
        }

        pseudoECGSignals
        {
            ecgSolver    pseudoECG;
        }
    }

"""


def _archive_output_dir(case_root: Path) -> Path:
    return case_root / "archivedPostProcessing"



def _workflow_dag_for(
    mesh_family: str,
    dimensions_list: list[str],
    *,
    case_root: Path,
    run_in_parallel: bool = False,
) -> dict[str, object]:
    if mesh_family == "tet":
        mesh_steps = [
            {"id": "clean", "command": "Allclean", "depends_on": []},
            {
                "id": "gmsh",
                "command": "gmsh",
                "args": [
                    "-3",
                    "setup/mesh/tet/three_domain_box.geo",
                    "-o",
                    "three_domain_box.msh",
                    "-format",
                    "msh2",
                ],
                "depends_on": ["clean"],
            },
            {
                "id": "gmshToFoam",
                "command": "gmshToFoam",
                "args": ["three_domain_box.msh"],
                "depends_on": ["gmsh"],
            },
            {"id": "checkMesh", "command": "checkMesh", "depends_on": ["gmshToFoam"]},
            {
                "id": "setConductivity",
                "command": "setTorsoOrganConductivityField",
                "depends_on": ["checkMesh"],
            },
        ]
        solve_depends_on = ["setConductivity"]
        steps, final_id = solve_steps(
            solve_id="solve",
            solve_command="cardiacFoam",
            depends_on=solve_depends_on,
            run_in_parallel=run_in_parallel,
            case_root=case_root,
        )
        # Reads the reconstructed final-time solution -- the live
        # (potentially parallel-decomposed) verifier can't do the
        # heart/bath fvMeshSubset + interface-face analysis itself
        # (applications/utilities/bathBidomainInterfaceMetrics
        # explicitly operates on "a reconstructed serial mesh",
        # which only exists once the solve step(s) have fully exited), so
        # this is its own step rather than folded into the solve step.
        # depends_on final_id (reconstructPar when parallel, solve
        # otherwise), not the literal "solve" id, since the reconstructed
        # mesh only exists after reconstructPar when parallel is enabled.
        interface_metrics_step = {
            "id": "interfaceMetrics",
            "command": "bathBidomainInterfaceMetrics",
            "args": ["-latestTime"],
            "depends_on": [final_id],
        }
        return {"steps": mesh_steps + steps + [interface_metrics_step]}

    mesh_steps = [
        {"id": "clean", "command": "Allclean", "depends_on": []},
        {
            "id": "mesh",
            "command": "blockMesh",
            "args": ["-dict", f"system/blockMeshDict.{dimensions_list[0]}"],
            "depends_on": ["clean"],
        },
        {"id": "topoSet", "command": "topoSet", "depends_on": ["mesh"]},
        {
            "id": "setConductivity",
            "command": "setTorsoOrganConductivityField",
            "depends_on": ["topoSet"],
        },
    ]
    steps, _final_id = solve_steps(
        solve_id="solve",
        solve_command="cardiacFoam",
        depends_on=["setConductivity"],
        run_in_parallel=run_in_parallel,
        case_root=case_root,
    )
    return {"steps": mesh_steps + steps}


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
    bath_predictor_corrector: bool = False,
    fda_bath_variant: str = "groundElectrode",
    mesh_family: str = "hex",
    numerics_profile: str | None = None,
    grad_scheme: str | None = None,
    phi_tolerance: float | None = None,
    end_time: float | None = None,
    fv_scheme_overrides: Sequence[Mapping[str, object]] | None = None,
    fv_solution_overrides: Sequence[Mapping[str, object]] | None = None,
    tet_geo_template_relpath: Path = Path("setup/mesh/tet/three_domain_box.geo.template"),
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
        f"{electro_properties_scope}.bathPredictorCorrector": bool(
            bath_predictor_corrector
        ),
        f"{electro_properties_scope}.verificationModel.type": verification_model_type,
        f"{electro_properties_scope}.verificationModel.fdaBathVariant": fda_bath_variant,
        f"{electro_properties_scope}.manufacturedBidomain.fdaBathVariant": fda_bath_variant,
    }

    if fda_bath_variant not in ("groundElectrode", "electrodePair"):
        raise ValueError(
            "fda_bath_variant must be 'groundElectrode' or 'electrodePair', "
            f"got {fda_bath_variant!r}"
        )

    if mesh_family == "tet":
        render_tet_geo(
            case_root,
            cells,
            template_relpath=tet_geo_template_relpath,
            geo_relpath=Path("setup/mesh/tet/three_domain_box.geo"),
        )
        # Must happen before the variant branch below: this copy resets
        # electro_properties to the tet template's own groundElectrode
        # defaults (its own groundPatches.xMin), which the electrodePair
        # branch then needs to remove. Doing it the other way around lets
        # this copy silently reintroduce the key the removal just cleared.
        shutil.copy(case_root / "setup" / "mesh" / "tet" / "electroProperties", electro_properties)
        for overlay_name in _TET_NUMERICS_PROFILES.get(numerics_profile or "", ()):
            shutil.copy(
                case_root / "setup" / "mesh" / "tet" / overlay_name,
                case_root / "system" / overlay_name,
            )
    else:
        try:
            cell_counts = defaults.BLOCK_MESH_RESOLUTION_BY_DIMENSION[dimension].format(cells=cells)
        except KeyError as exc:
            raise ValueError(f"Unsupported dimension: {dimension}") from exc
        replace_block_mesh_resolutions(block_mesh_dict, cell_counts, expected_blocks=3)

    # The two FDA bidomain-with-bath variants differ in their outer bath
    # boundary conditions, and the dictionary has to follow the verifier or the
    # reported norms describe a different problem than the one solved.
    #   groundElectrode: Dirichlet phiE = 0 at x = -1, I_E = +alpha at x = 2.
    #   electrodePair:   I_E = -alpha at x = -1 and +alpha at x = 2, no ground.
    #                    The integral of I_E over the boundary is zero so the
    #                    problem is solvable, but phiE floats and needs a
    #                    reference point to pin the constant.
    bath_scope = f"{electro_properties_scope}.bathPotentialDomain"
    if fda_bath_variant == "electrodePair":
        # groundPatches and surfaceCurrentPatches are mutually exclusive per
        # patch (extracellularPotentialDomain.C rejects a patch listed in
        # both). The checked-in electroProperties defaults to groundElectrode
        # and so carries groundPatches.xMin; switching variants must remove
        # it, not just add the electrodePair surfaceCurrentPatches.xMin.
        remove_electro_property_dict(
            electro_properties,
            "xMin",
            scope=[electro_properties_scope, "bathPotentialDomain", "groundPatches"],
            missing_ok=True,
        )
        ref_y, ref_z = _phi_e_ref_point_yz(dimension, cells)
        case_overrides.update(
            {
                f"{bath_scope}.surfaceCurrentPatches.xMin": -defaults.FDA_ALPHA,
                f"{bath_scope}.surfaceCurrentPatches.xMax": defaults.FDA_ALPHA,
                f"{bath_scope}.phiERefPoint": f"(-0.9 {ref_y} {ref_z})",
                f"{bath_scope}.phiEReferenceValue": 0.0,
            }
        )
    else:
        # Symmetric cleanup: a prior electrodePair case sharing this
        # case_root may have left surfaceCurrentPatches.xMin behind, which
        # would collide with groundPatches.xMin below the same way.
        remove_electro_property_dict(
            electro_properties,
            "xMin",
            scope=[electro_properties_scope, "bathPotentialDomain", "surfaceCurrentPatches"],
            missing_ok=True,
        )
        case_overrides.update(
            {
                f"{bath_scope}.groundPatches.xMin": 0.0,
                f"{bath_scope}.surfaceCurrentPatches.xMax": defaults.FDA_ALPHA,
            }
        )

    if ecg_enabled:
        ensure_electro_property_dict(
            electro_properties,
            "ecgDomains",
            _DEFAULT_ECG_DOMAINS_BLOCK,
            scope=electro_properties_scope,
        )
        case_overrides.update(
            {
                f"{electro_properties_scope}.ecgDomains.bodyECG.ecgSolver": "torsoECG",
                f"{electro_properties_scope}.ecgDomains.bodyECG.ecgVerificationModel":
                    "bathECGManufacturedVerifier",
                f"{electro_properties_scope}.ecgDomains.pseudoECGSignals.ecgSolver": "pseudoECG",
            }
        )

    set_delta_t(control_dict, dt_value)
    if end_time is not None:
        update_foam_entry(control_dict, "endTime", end_time)
        # writeControl is adjustableRunTime (time-based, not step-count-based)
        # so every case in a temporal-convergence sweep writes a
        # reconstructable time regardless of how few steps its deltaT takes
        # to reach endTime; writeInterval must track an overridden endTime
        # or it stays pinned to the checked-in default and stops matching.
        update_foam_entry(control_dict, "writeInterval", end_time)
    if grad_scheme is not None:
        update_foam_entry(
            case_root / "system" / "fvSchemes",
            "default",
            _GRAD_SCHEME_TOKENS[grad_scheme],
            scope=["gradSchemes"],
        )
    if phi_tolerance is not None:
        update_foam_entry(
            case_root / "system" / "fvSolution",
            "tolerance",
            phi_tolerance,
            scope=["solvers", '"phiE|phiEFinal|phiI|phiIFinal"'],
        )
    for entry in fv_scheme_overrides or ():
        update_foam_entry(
            case_root / "system" / "fvSchemes", entry["key"], entry["value"],
            scope=entry.get("scope"),
        )
    for entry in fv_solution_overrides or ():
        update_foam_entry(
            case_root / "system" / "fvSolution", entry["key"], entry["value"],
            scope=entry.get("scope"),
        )
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

    try:
        subprocess.run(command, check=True)
    finally:
        archive_case_logs(case_root, case.case_id)

    destination_dir = _archive_output_dir(case_root)
    filename = _case_output_filename(case)
    stage_post_processing_outputs(
        case_root, destination_dir, {filename: filename}
    )

    if ecg_enabled:
        ecg_mapping = {
            "torsoECG.dat": f"BathECG_{case.case_id}_torsoECG.dat",
            "manufacturedBathECG.dat": f"BathECG_{case.case_id}_manufacturedBathECG.dat",
            "manufacturedBathECGSummary.dat": f"BathECG_{case.case_id}_manufacturedBathECGSummary.dat",
            "pseudoECG.dat": f"PseudoECG_{case.case_id}_pseudoECG.dat",
        }
        stage_post_processing_outputs(case_root, destination_dir, ecg_mapping)


def _collect_outputs(
    case_root: Path,
    output_dir: Path,
    *,
    ecg_enabled: bool = False,
) -> None:
    archived_dir = _archive_output_dir(case_root)
    archived_outputs = []
    if archived_dir.exists():
        for source in sorted(archived_dir.glob("*.dat")):
            if source.name.startswith("bathBidomain_") or (
                ecg_enabled
                and source.name.startswith(("BathECG_", "PseudoECG_"))
            ):
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
    bath_predictor_corrector: bool = False,
    fda_bath_variant: str = "groundElectrode",
    mesh_family: str = "hex",
    numerics_profile: str | None = None,
    grad_scheme: str | None = None,
    phi_tolerance: float | None = None,
    end_time: float | None = None,
    fv_scheme_overrides: Sequence[Mapping[str, object]] | None = None,
    fv_solution_overrides: Sequence[Mapping[str, object]] | None = None,
    tet_geo_template_relpath: str | Path = "setup/mesh/tet/three_domain_box.geo.template",
) -> TutorialSpec:
    mesh_family = str(mesh_family)
    if mesh_family not in {"hex", "tet"}:
        raise ValueError(f"Unsupported mesh_family '{mesh_family}'. Expected 'hex' or 'tet'.")
    if mesh_family == "tet" and [str(item) for item in dimensions] != ["3D"]:
        raise ValueError('mesh_family="tet" is only supported for dimensions=["3D"].')
    if numerics_profile is not None and numerics_profile not in _TET_NUMERICS_PROFILES:
        raise ValueError(
            f"Unsupported numerics_profile '{numerics_profile}'. "
            f"Expected one of {sorted(_TET_NUMERICS_PROFILES)}."
        )
    if grad_scheme is not None and grad_scheme not in _GRAD_SCHEME_TOKENS:
        raise ValueError(
            f"Unsupported grad_scheme '{grad_scheme}'. "
            f"Expected one of {sorted(_GRAD_SCHEME_TOKENS)}."
        )
    if phi_tolerance is not None and float(phi_tolerance) <= 0.0:
        raise ValueError("phi_tolerance must be positive.")

    dimensions_list = [str(item) for item in dimensions]
    cells_list = [int(item) for item in number_cells]
    dt_values_list = [float(item) for item in dt_values]
    solver_types_list = [str(item) for item in solver_types]

    if piecewise_sweep and len(cells_list) != len(dt_values_list):
        raise ValueError(
            "piecewise_sweep requires number_cells and dt_values to have the same length"
        )

    tet_geo_template_path = Path(tet_geo_template_relpath)

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
            bath_predictor_corrector=bath_predictor_corrector,
            fda_bath_variant=fda_bath_variant,
            mesh_family=mesh_family,
            numerics_profile=numerics_profile,
            grad_scheme=grad_scheme,
            phi_tolerance=phi_tolerance,
            end_time=end_time,
            fv_scheme_overrides=fv_scheme_overrides,
            fv_solution_overrides=fv_solution_overrides,
            tet_geo_template_relpath=tet_geo_template_path,
        ),
        run_case=partial(
            _run_case,
            tutorials_root=tutorials_root,
            run_script_relpath=Path(run_script_relpath),
            run_in_parallel=run_in_parallel,
            ecg_enabled=ecg_enabled,
        ),
        collect_outputs=partial(_collect_outputs, ecg_enabled=ecg_enabled),
        postprocess=partial(
            _postprocess,
            tutorial_name=tutorial_name,
            postprocess_script_relpath=Path(postprocess_script_relpath),
            postprocess_function_name=postprocess_function_name,
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "notes": "FDA bath-bidomain manufactured-solution convergence benchmark",
            "workflow_dag": _workflow_dag_for(
                mesh_family, dimensions_list,
                case_root=case_root, run_in_parallel=run_in_parallel,
            ),
            "dimensions": dimensions_list,
            "solver_types": solver_types_list,
            "piecewise_sweep": piecewise_sweep,
            "mesh_family": mesh_family,
            "numerics_profile": numerics_profile,
            "grad_scheme": grad_scheme,
            "phi_tolerance": phi_tolerance,
            "end_time": end_time,
            "control_dict_relpath": str(control_dict_relpath),
            "electro_properties_relpath": str(electro_properties_relpath),
            "physics_properties_relpath": str(physics_properties_relpath),
            "electro_properties_scope": electro_properties_scope,
            "block_mesh_dict_template": block_mesh_dict_template,
            "run_script_relpath": str(run_script_relpath),
            "run_in_parallel": run_in_parallel,
            "ecg_enabled": ecg_enabled,
            "bath_predictor_corrector": bool(bath_predictor_corrector),
            "fda_bath_variant": str(fda_bath_variant),
            "postprocess_script_relpath": str(postprocess_script_relpath),
            "postprocess_function_name": postprocess_function_name,
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
            "tet_geo_template_relpath": str(tet_geo_template_path),
        },
    )
