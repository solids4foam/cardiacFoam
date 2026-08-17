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
#     Defines configuration template for manufactured bidomain scenarios for FDA verification.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from ..ids import CardiacTutorialID

from .manufactured_monodomain_pseudo_ecg import (
    BLOCK_MESH_DICT_TEMPLATE,
    BLOCK_MESH_RESOLUTION_BY_DIMENSION,
    DT_VALUES,
    DIMENSIONS,
    NUMBER_CELLS,
    OUTPUT_DIR_NAME,
    PIECEWISE_SWEEP,
    RUN_CASE_SCRIPT_RELPATH,
)


TUTORIAL_NAME = CardiacTutorialID.MANUFACTURED_BIDOMAIN.value
CASE_DIR_NAME = "manufacturedSolutions/bidomain"
SETUP_DIR_NAME = "setup"
SOLVER_TYPES = ("implicit",)
ELECTRO_PROPERTIES_SCOPE = "bidomainSolverCoeffs"
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
RUN_IN_PARALLEL = True
VERIFICATION_MODEL_TYPE = "manufacturedFDABidomainVerifier"
ECG_ENABLED = False
ECG_REFERENCE_QUADRATURE_ORDER = 96
ECG_CHECK_QUADRATURE_ORDERS = (6, 12, 24, 48)
ECG_ELECTRODES_BY_DIMENSION = {}

__all__ = [
    "BLOCK_MESH_DICT_TEMPLATE",
    "BLOCK_MESH_RESOLUTION_BY_DIMENSION",
    "CASE_DIR_NAME",
    "DT_VALUES",
    "DIMENSIONS",
    "ECG_CHECK_QUADRATURE_ORDERS",
    "ECG_ELECTRODES_BY_DIMENSION",
    "ECG_ENABLED",
    "ECG_REFERENCE_QUADRATURE_ORDER",
    "ELECTRO_PROPERTIES_SCOPE",
    "NUMBER_CELLS",
    "OUTPUT_DIR_NAME",
    "PIECEWISE_SWEEP",
    "RUN_IN_PARALLEL",
    "RUN_SCRIPT_RELPATH",
    "SETUP_DIR_NAME",
    "SOLVER_TYPES",
    "TUTORIAL_NAME",
    "VERIFICATION_MODEL_TYPE",
]
