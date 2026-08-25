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
#     manufactured_monodomain_total_lagrangian_em
#
# Description
#     Defines configuration template for boundary-constrained manufactured
#     electromechanics scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from ..ids import CardiacTutorialID

from pathlib import Path

from .manufactured_monodomain_pseudo_ecg import (
    BLOCK_MESH_RESOLUTION_BY_DIMENSION,
    DT_VALUES,
    DIMENSIONS,
    NUMBER_CELLS,
    OUTPUT_DIR_NAME,
    PIECEWISE_SWEEP,
)
from .shared import CONTROL_DICT_RELPATH, OUTPUT_RELPATH, RUN_CASE_SCRIPT_RELPATH

DIMENSIONS = ["3D"]
NUMBER_CELLS = [10, 20, 40, 80]
DT_VALUES = [0.00892857, 0.00224215, 0.000560538, 0.000140174]
OUTPUT_DIR_NAME = "postProcessing"
OUTPUT_RELPATH = Path("postProcessing")

TUTORIAL_NAME = CardiacTutorialID.MANUFACTURED_MONODOMAIN_TOTAL_LAGRANGIAN_EM.value
CASE_DIR_NAME = "manufacturedSolutions/monodomainTotalLagrangianEM"
SETUP_DIR_NAME = "setup"
SOLVER_TYPES = ("implicit",)
CONTROL_DICT_RELPATH = CONTROL_DICT_RELPATH
ELECTRO_PROPERTIES_RELPATH = Path("constant/electro/electroProperties")
ELECTROMECHANICAL_PROPERTIES_RELPATH = Path("constant/electroMechanicalProperties")
PHYSICS_PROPERTIES_RELPATH = Path("constant/physicsProperties")
BLOCK_MESH_DICT_TEMPLATE = "system/blockMeshDict.{dimension}"
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
RUN_IN_PARALLEL = True
ELECTRO_PROPERTIES_SCOPE = "monodomainSolverCoeffs"
ELECTROMECHANICAL_VERIFICATION_MODEL_TYPE = "manufacturedElectromechanicsVerifier"

__all__ = [
    "BLOCK_MESH_DICT_TEMPLATE",
    "BLOCK_MESH_RESOLUTION_BY_DIMENSION",
    "CASE_DIR_NAME",
    "CONTROL_DICT_RELPATH",
    "DIMENSIONS",
    "DT_VALUES",
    "ELECTRO_PROPERTIES_RELPATH",
    "ELECTRO_PROPERTIES_SCOPE",
    "ELECTROMECHANICAL_PROPERTIES_RELPATH",
    "ELECTROMECHANICAL_VERIFICATION_MODEL_TYPE",
    "NUMBER_CELLS",
    "OUTPUT_DIR_NAME",
    "OUTPUT_RELPATH",
    "PHYSICS_PROPERTIES_RELPATH",
    "PIECEWISE_SWEEP",
    "RUN_IN_PARALLEL",
    "RUN_SCRIPT_RELPATH",
    "SETUP_DIR_NAME",
    "SOLVER_TYPES",
    "TUTORIAL_NAME",
]
