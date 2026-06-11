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
#     shared
#
# Description
#     Provides shared default parameters and constants for case generation.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path


CONTROL_DICT_RELPATH = Path("system/controlDict")
ELECTRO_PROPERTIES_RELPATH = Path("constant/electroProperties")
RUN_CASE_SCRIPT_RELPATH = Path("applications/scripts/driverFoam/openfoam_driver/scripts/run_case.sh")
OUTPUT_DIR_NAME = "postProcessing"
OUTPUT_RELPATH = Path(OUTPUT_DIR_NAME)
