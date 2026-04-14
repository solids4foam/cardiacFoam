from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

from ...core.defaults import manufactured_fda_bidomain as defaults
from .manufactured_fda import make_spec as make_base_spec


def make_spec(
    *,
    tutorials_root: Path | None = None,
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
):
    return make_base_spec(
        tutorials_root=tutorials_root,
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
    )
