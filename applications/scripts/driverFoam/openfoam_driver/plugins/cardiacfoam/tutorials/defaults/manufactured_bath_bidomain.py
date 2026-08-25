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
#     manufactured_bath_bidomain
#
# Description
#     Defines configuration template for manufactured bath bidomain scenarios for FDA verification.
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
    RUN_CASE_SCRIPT_RELPATH,
)

# Surface stimulus magnitude of the FDA bidomain-with-bath problems (A/m^2).
# Section 3.3 of the FDA document: I_E = -alpha at x = -1 and +alpha at x = 2.
FDA_ALPHA = 0.01

# y/z domain extent per dimension, read off each blockMeshDict.<dim>'s own
# vertices block (x always spans -1..2 for all three; only the y/z slab
# thickness differs -- 1D and 2D are thin slivers, not full 1x1 cross
# sections). Paired with which of y/z BLOCK_MESH_RESOLUTION_BY_DIMENSION
# actually subdivides by `cells` (1D: neither: both fixed at 1 cell; 2D:
# y only, z fixed at 1; 3D: both) -- a fixed-at-1-cell direction is always
# safe at its exact center regardless of resolution.
DOMAIN_YZ_EXTENT_BY_DIMENSION: dict[str, tuple[float, float]] = {
    "1D": (0.1, 0.1),
    "2D": (1.0, 0.05),
    "3D": (1.0, 1.0),
}
YZ_SUBDIVIDED_BY_DIMENSION: dict[str, tuple[bool, bool]] = {
    "1D": (False, False),
    "2D": (True, False),
    "3D": (True, True),
}



TUTORIAL_NAME = CardiacTutorialID.MANUFACTURED_BATH_BIDOMAIN.value
CASE_DIR_NAME = "manufacturedSolutions/bathBidomain"
SETUP_DIR_NAME = "setup"
SOLVER_TYPES = ("implicit",)
ELECTRO_PROPERTIES_SCOPE = "bidomainSolverCoeffs"
BLOCK_MESH_DICT_TEMPLATE = "system/blockMeshDict.{dimension}"
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
RUN_IN_PARALLEL = True
VERIFICATION_MODEL_TYPE = "manufacturedFDABathBidomainVerifier"

__all__ = [
    "BLOCK_MESH_DICT_TEMPLATE",
    "BLOCK_MESH_RESOLUTION_BY_DIMENSION",
    "CASE_DIR_NAME",
    "DOMAIN_YZ_EXTENT_BY_DIMENSION",
    "DT_VALUES",
    "DIMENSIONS",
    "ELECTRO_PROPERTIES_SCOPE",
    "FDA_ALPHA",
    "NUMBER_CELLS",
    "OUTPUT_DIR_NAME",
    "PIECEWISE_SWEEP",
    "RUN_IN_PARALLEL",
    "RUN_SCRIPT_RELPATH",
    "SETUP_DIR_NAME",
    "SOLVER_TYPES",
    "TUTORIAL_NAME",
    "VERIFICATION_MODEL_TYPE",
    "YZ_SUBDIVIDED_BY_DIMENSION",
]
