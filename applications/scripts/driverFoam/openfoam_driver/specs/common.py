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
#     common
#
# Description
#     Provides shared definitions and defaults for specification templates.
#     (Refactored into a facade: imports from detection, overrides, paths, utils)
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from .detection import (
    _AT_EXPORT_RE,
    _BLOCK_DECL_RE,
    _IONIC_EXPORT_RE,
    detect_active_tension_export_list,
    detect_active_tension_model_name,
    detect_electro_coeffs_scope,
    detect_ionic_export_list,
    detect_ionic_model_name,
    detect_myocardium_solver_name,
    detect_verification_model_type,
    electro_properties_has_block,
)
from .overrides import (
    _resolve_scope_tokens,
    apply_electro_property_overrides,
    apply_entry_overrides,
    apply_physics_property_overrides,
    ensure_electro_property_dict,
    normalize_entry_overrides,
    remove_electro_property_dict,
)
from .paths import (
    default_setup_dir_name,
    repo_root_default,
    resolve_run_script_path,
    resolve_spec_paths,
    tutorials_root_default,
)
from .utils import (
    collect_outputs_by_pattern,
    load_python_module,
    replace_single_block_mesh_resolution,
    set_delta_t,
    set_end_time,
)

__all__ = [
    "repo_root_default",
    "tutorials_root_default",
    "default_setup_dir_name",
    "detect_myocardium_solver_name",
    "detect_electro_coeffs_scope",
    "detect_ionic_model_name",
    "detect_ionic_export_list",
    "electro_properties_has_block",
    "detect_verification_model_type",
    "detect_active_tension_model_name",
    "detect_active_tension_export_list",
    "_resolve_scope_tokens",
    "normalize_entry_overrides",
    "apply_entry_overrides",
    "apply_electro_property_overrides",
    "apply_physics_property_overrides",
    "remove_electro_property_dict",
    "ensure_electro_property_dict",
    "resolve_spec_paths",
    "resolve_run_script_path",
    "load_python_module",
    "collect_outputs_by_pattern",
    "set_delta_t",
    "set_end_time",
    "replace_single_block_mesh_resolution",
    "_IONIC_EXPORT_RE",
    "_BLOCK_DECL_RE",
    "_AT_EXPORT_RE",
]
