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
#     manufactured_purkinje_graph
#
# Description
#     Defaults for manufactured Purkinje graph convergence scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from ..ids import CardiacTutorialID

from pathlib import Path


TUTORIAL_NAME = CardiacTutorialID.MANUFACTURED_PURKINJE_GRAPH.value
CASE_DIR_NAME = "manufacturedSolutions/monodomain1D3D"
SETUP_DIR_NAME = "setup"
OUTPUT_DIR_NAME = "outputs/1dGraphConvergence"
GRAPH_IDS = ("nodes003", "nodes011", "nodes021", "nodes041", "nodes081", "nodes161")
N_STEPS = 1427
DELTA_T = 0.000140174
CONTROL_DICT_RELPATH = Path("system/controlDict")
ELECTRO_PROPERTIES_RELPATH = Path("constant/electroProperties")
BLOCK_MESH_DICT_RELPATH = Path("system/blockMeshDict.3D")
