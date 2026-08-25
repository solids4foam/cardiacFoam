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
#     monodomain_and_eikonal_1d_cable_cv_convergence
#
# Description
#     Defaults for 1D cable CV convergence studies.
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


TUTORIAL_NAME = CardiacTutorialID.CABLE_1D_CV_CONVERGENCE.value
CASE_DIR_NAME = "electrophysiologyProtocols/cableProtocol/monodomain1DCableCV"
SETUP_DIR_NAME = "setup"
DEFAULT_OUTPUT_DIR_NAME = "outputsCVConvergence"
DX_VALUES = (0.5, 0.2, 0.1)  # mm
DT_VALUES = (0.02, 0.01, 0.005)  # ms
CONDUCTIVITY_VALUES = (
    "[-1 -3 3 0 0 2 0] (0.1334 0 0 0.1334 0 0.1334)",
)
IONIC_MODELS = ("Stewart",)
IONIC_MODEL_TISSUE_MAP = {
    "Stewart": ("myocyte",),
}
SOLVERS = ("implicit",)
PARALLEL = True
END_TIME_S = 0.04
ELECTRO_PROPERTIES_SCOPE = "monodomainSolverCoeffs"
ELECTRO_PROPERTIES_RELPATH = SHARED_ELECTRO_PROPERTIES_RELPATH
BLOCK_MESH_DICT_RELPATH = Path("system/blockMeshDict")
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
CABLE_LENGTH_MM = 20.0
CROSS_SECTION_CELL_COUNTS = (1, 1)
CONTROL_DICT_RELPATH = CONTROL_DICT_RELPATH

__all__ = [
    "TUTORIAL_NAME",
    "CASE_DIR_NAME",
    "SETUP_DIR_NAME",
    "DEFAULT_OUTPUT_DIR_NAME",
    "DX_VALUES",
    "DT_VALUES",
    "CONDUCTIVITY_VALUES",
    "IONIC_MODELS",
    "IONIC_MODEL_TISSUE_MAP",
    "SOLVERS",
    "PARALLEL",
    "END_TIME_S",
    "ELECTRO_PROPERTIES_SCOPE",
    "ELECTRO_PROPERTIES_RELPATH",
    "BLOCK_MESH_DICT_RELPATH",
    "RUN_SCRIPT_RELPATH",
    "CABLE_LENGTH_MM",
    "CROSS_SECTION_CELL_COUNTS",
    "CONTROL_DICT_RELPATH",
    "OUTPUT_DIR_NAME",
]
