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
#     manufactured_bidomain
#
# Description
#     Defines configuration template for manufactured bidomain scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

from openfoam_driver.plugins.cardiacfoam.defaults import manufactured_bidomain as defaults
from .manufactured_monodomain_pseudo_ecg import make_spec as make_base_spec


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
    ecg_enabled: bool = defaults.ECG_ENABLED,
    ecg_reference_quadrature_order: int = defaults.ECG_REFERENCE_QUADRATURE_ORDER,
    ecg_check_quadrature_orders: Sequence[int] = defaults.ECG_CHECK_QUADRATURE_ORDERS,
    ecg_electrodes_by_dimension: Mapping[str, Mapping[str, str]] = defaults.ECG_ELECTRODES_BY_DIMENSION,
    block_mesh_dict_template: str = defaults.BLOCK_MESH_DICT_TEMPLATE,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    postprocess_script_relpath: str | Path = "post_processing_manufactured.py",
    postprocess_function_name: str = "run_postprocessing",
    run_in_parallel: bool = defaults.RUN_IN_PARALLEL,
    postprocess_strict_artifacts: bool = False,
    mesh_family: str = "hex",
    numerics_profile: str | None = None,
    grad_scheme: str | None = None,
    phi_tolerance: float | None = None,
    end_time: float | None = None,
    fv_scheme_overrides: Sequence[Mapping[str, object]] | None = None,
    fv_solution_overrides: Sequence[Mapping[str, object]] | None = None,
):
    return make_base_spec(
        tutorials_root=tutorials_root,
        tutorial_name=tutorial_name,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        number_cells=number_cells,
        dt_values=dt_values,
        dimensions=dimensions,
        solver_types=solver_types,
        piecewise_sweep=piecewise_sweep,
        electro_properties_scope=electro_properties_scope,
        control_dict_relpath=control_dict_relpath,
        electro_properties_relpath=electro_properties_relpath,
        physics_properties_relpath=physics_properties_relpath,
        electro_property_overrides=electro_property_overrides,
        physics_property_overrides=physics_property_overrides,
        verification_model_type=verification_model_type,
        ecg_enabled=ecg_enabled,
        ecg_reference_quadrature_order=ecg_reference_quadrature_order,
        ecg_check_quadrature_orders=ecg_check_quadrature_orders,
        ecg_electrodes_by_dimension=ecg_electrodes_by_dimension,
        block_mesh_dict_template=block_mesh_dict_template,
        run_script_relpath=run_script_relpath,
        postprocess_script_relpath=postprocess_script_relpath,
        postprocess_function_name=postprocess_function_name,
        run_in_parallel=run_in_parallel,
        postprocess_strict_artifacts=postprocess_strict_artifacts,
        mesh_family=mesh_family,
        numerics_profile=numerics_profile,
        grad_scheme=grad_scheme,
        phi_tolerance=phi_tolerance,
        end_time=end_time,
        fv_scheme_overrides=fv_scheme_overrides,
        fv_solution_overrides=fv_solution_overrides,
    )
