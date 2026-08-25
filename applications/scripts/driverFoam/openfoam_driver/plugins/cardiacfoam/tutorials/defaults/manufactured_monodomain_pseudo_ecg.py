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
#     manufactured_monodomain_pseudo_ecg
#
# Description
#     Defines configuration template for manufactured monodomain with pseudo ECG
#     for FDA verification.
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
    RUN_CASE_SCRIPT_RELPATH,
)


TUTORIAL_NAME = CardiacTutorialID.MANUFACTURED_MONODOMAIN_PSEUDO_ECG.value
CASE_DIR_NAME = "manufacturedSolutions/monodomainPseudoECG"
SETUP_DIR_NAME = "setup"
NUMBER_CELLS = (10, 20, 40, 80)
DT_VALUES = (
    0.00892857,
    0.00224215,
    0.000560538,
    0.000140174,
)
DIMENSIONS = ("1D", "2D", "3D")
SOLVER_TYPES = ("implicit",)
PIECEWISE_SWEEP = True
ELECTRO_PROPERTIES_SCOPE = "monodomainSolverCoeffs"
ELECTRO_PROPERTIES_RELPATH = SHARED_ELECTRO_PROPERTIES_RELPATH
BLOCK_MESH_DICT_TEMPLATE = "system/blockMeshDict.{dimension}"
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
RUN_IN_PARALLEL = True
VERIFICATION_MODEL_TYPE = "manufacturedFDAMonodomainVerifier"
ECG_ENABLED = True
ECG_REFERENCE_QUADRATURE_ORDER = 96
ECG_CHECK_QUADRATURE_ORDERS = (6, 12, 24, 48)
ECG_ELECTRODES_BY_DIMENSION = {
    "1D": {
        "E1": "(-0.5 0 0)",
        "E2": "(1.5 0 0)",
        "E3": "(1.2 0 0)",
        "E4": "(1.35 0 0)",
        "E5": "(1.65 0 0)",
    },
    "2D": {
        "E1": "(-0.5 0.5 0)",
        "E2": "(1.5 0.5 0)",
        "E3": "(1.2 0.23 0)",
        "E4": "(1.35 0.78 0)",
        "E5": "(0.18 1.35 0)",
    },
    "3D": {
        "E1": "(-0.5 0.5 0.5)",
        "E2": "(1.5 0.5 0.5)",
        "E3": "(1.2 0.23 0.61)",
        "E4": "(1.35 0.74 0.28)",
        "E5": "(1.55 0.41 0.83)",
    },
}
BLOCK_MESH_RESOLUTION_BY_DIMENSION = {
    "1D": "{cells} 1 1",
    "2D": "{cells} {cells} 1",
    "3D": "{cells} {cells} {cells}",
}
