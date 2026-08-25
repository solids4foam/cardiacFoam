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
#     (Refactored into a facade: imports from paths, utils.)  Solver-specific
#     helpers are NOT re-exported here: cardiac dictionary detection and
#     electro/physics-property overrides live in
#     `plugins/cardiacfoam/{detection,overrides}.py` and must be imported from
#     there directly.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from .paths import (
    default_setup_dir_name,
    repo_root_default,
    resolve_run_script_path,
    resolve_spec_paths,
    tutorials_root_default,
)
from .utils import (
    load_python_module,
    replace_block_mesh_resolutions,
    set_delta_t,
    set_end_time,
)

__all__ = [
    "repo_root_default",
    "tutorials_root_default",
    "default_setup_dir_name",
    "resolve_spec_paths",
    "resolve_run_script_path",
    "load_python_module",
    "set_delta_t",
    "set_end_time",
    "replace_block_mesh_resolutions",
]
