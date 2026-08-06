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
#     manufactured_fda_bath_bidomain
#
# Description
#     Defines configuration template for manufactured FDA bath bidomain scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path

from .manufactured_fda import (
    BLOCK_MESH_RESOLUTION_BY_DIMENSION,
    DT_VALUES,
    DIMENSIONS,
    NUMBER_CELLS,
    OUTPUT_DIR_NAME,
    PIECEWISE_SWEEP,
    RUN_CASE_SCRIPT_RELPATH,
)


TUTORIAL_NAME = "manufacturedFDABathBidomain"
CASE_DIR_NAME = "manufacturedSolutions/bathBidomain"
SETUP_DIR_NAME = "setup"
SOLVER_TYPES = ("implicit",)
ELECTRO_PROPERTIES_SCOPE = "bidomainSolverCoeffs"
BLOCK_MESH_DICT_TEMPLATE = "system/blockMeshDict.{dimension}"
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
POSTPROCESS_SCRIPT_RELPATH = Path("post_processing_manufactured_bath.py")
POSTPROCESS_FUNCTION_NAME = "run_postprocessing"
RUN_IN_PARALLEL = True
VERIFICATION_MODEL_TYPE = "manufacturedFDABathBidomainVerifier"

__all__ = [
    "BLOCK_MESH_DICT_TEMPLATE",
    "BLOCK_MESH_RESOLUTION_BY_DIMENSION",
    "CASE_DIR_NAME",
    "DT_VALUES",
    "DIMENSIONS",
    "ELECTRO_PROPERTIES_SCOPE",
    "NUMBER_CELLS",
    "OUTPUT_DIR_NAME",
    "PIECEWISE_SWEEP",
    "POSTPROCESS_FUNCTION_NAME",
    "POSTPROCESS_SCRIPT_RELPATH",
    "RUN_IN_PARALLEL",
    "RUN_SCRIPT_RELPATH",
    "SETUP_DIR_NAME",
    "SOLVER_TYPES",
    "TUTORIAL_NAME",
    "VERIFICATION_MODEL_TYPE",
]
