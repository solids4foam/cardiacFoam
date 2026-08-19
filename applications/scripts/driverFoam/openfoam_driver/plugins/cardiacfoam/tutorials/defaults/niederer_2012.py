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

from ..ids import CardiacTutorialID

from pathlib import Path

from .shared import (
    CONTROL_DICT_RELPATH,
    ELECTRO_PROPERTIES_RELPATH as SHARED_ELECTRO_PROPERTIES_RELPATH,
    OUTPUT_DIR_NAME,
    OUTPUT_RELPATH,
)




TUTORIAL_NAME = CardiacTutorialID.NIEDERER_2012.value
CASE_DIR_NAME = "NiedererEtAl2011/NiedererEtAl2011verification"
SETUP_DIR_NAME = "setup"
DX_VALUES = (0.5, 0.2, 0.1)  # in mm
DT_VALUES = (0.01, 0.005, 0.001)  # in ms
IONIC_MODELS = ("TNNP",)
IONIC_MODEL_TISSUE_MAP = {
    "TNNP": ("epicardialCells",),
}
SOLVERS = ("explicit", "implicit")
ELECTRO_PROPERTIES_SCOPE = "monodomainSolverCoeffs"
ELECTRO_PROPERTIES_RELPATH = SHARED_ELECTRO_PROPERTIES_RELPATH
SLAB_SIZE_MM = (20.0, 3.0, 7.0)
END_TIME_BY_DX = {
    0.5: 0.2,
    0.2: 0.08,
    0.1: 0.055,
}
BLOCK_MESH_DICT_RELPATH = Path("system/blockMeshDict")

# Explicit exports used by the Niederer tutorial spec.
#
# Deliberately NOT here: NIEDERER_POINTS/_LINE_START/_LINE_END/_NUM_POINTS,
# NIEDERER_POINTS_FUNCTION_OBJECT/_LINE_FUNCTION_OBJECT, NIEDERER_SAMPLED_FIELD,
# RUN_SCRIPT_RELPATH, and every *_POSTPROCESS_* constant. All of them existed
# solely to feed niederer_2012.py's _run_case/_export_openfoam_samples --
# code that TutorialSpec.run_case never actually called, and that duplicated,
# with hardcoded probe labels, what setup/convert_raw_samples.py now reads
# directly from the raw OpenFOAM probes file's own header. Removed alongside
# that dead code in the niederer_2012.py cleanup that follows 32cc29dd.
# TutorialSpec.run_case itself has since been removed entirely.
__all__ = [
    "TUTORIAL_NAME",
    "CASE_DIR_NAME",
    "SETUP_DIR_NAME",
    "DX_VALUES",
    "DT_VALUES",
    "IONIC_MODELS",
    "IONIC_MODEL_TISSUE_MAP",
    "SOLVERS",
    "ELECTRO_PROPERTIES_SCOPE",
    "ELECTRO_PROPERTIES_RELPATH",
    "CONTROL_DICT_RELPATH",
    "SLAB_SIZE_MM",
    "END_TIME_BY_DX",
    "BLOCK_MESH_DICT_RELPATH",
    "OUTPUT_DIR_NAME",
    "OUTPUT_RELPATH",
]
