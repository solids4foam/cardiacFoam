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
    RUN_CASE_SCRIPT_RELPATH,
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
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
NIEDERER_POINTS_FUNCTION_OBJECT = "Niedererpoints"
NIEDERER_LINE_FUNCTION_OBJECT = "Niedererlines"
NIEDERER_SAMPLED_FIELD = "activationTime"
NIEDERER_POINTS = (
    ("P1", 0.0, 0.0, 0.007),
    ("P2", 0.0, 0.0, 0.0),
    ("P3", 0.019999, 0.0, 0.007),
    ("P4", 0.019999, 0.0, 0.0),
    ("P5", 0.0, 0.003, 0.007),
    ("P6", 0.0, 0.003, 0.0),
    ("P7", 0.019999, 0.003, 0.007),
    ("P8", 0.019999, 0.003, 0.0),
    ("P9", 0.01, 0.0015, 0.0035),
)
NIEDERER_LINE_START = (0.0, 0.0, 0.007)
NIEDERER_LINE_END = (0.02, 0.003, 0.0)
NIEDERER_LINE_NUM_POINTS = 101
LINE_POSTPROCESS_RELPATH = Path("line_postProcessing.py")
POINTS_POSTPROCESS_RELPATH = Path("points_postProcessing.py")
CACHE_POSTPROCESS_RELPATH = Path("convert_raw_samples.py")
LINE_POSTPROCESS_FUNCTION = "run_postprocessing"
POINTS_POSTPROCESS_FUNCTION = "run_postprocessing"
CACHE_POSTPROCESS_FUNCTION = "run_postprocessing"
CASE_POSTPROCESS_CACHE_DIRNAME = "cachedCasePostProcessing"
TABLE_SUMMARY_RELPATH = Path("table_summary.py")
EXCEL_REFERENCE_RELPATH = Path(
    "postProcessing/Niederer_graphs_webplotdigitilizer_points_slab/WebPlotDigitilizerdata.xlsx"
)

# Explicit exports used by the Niederer tutorial spec.
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
    "RUN_SCRIPT_RELPATH",
    "RUN_CASE_SCRIPT_RELPATH",
    "OUTPUT_DIR_NAME",
    "OUTPUT_RELPATH",
    "NIEDERER_POINTS_FUNCTION_OBJECT",
    "NIEDERER_LINE_FUNCTION_OBJECT",
    "NIEDERER_SAMPLED_FIELD",
    "NIEDERER_POINTS",
    "NIEDERER_LINE_START",
    "NIEDERER_LINE_END",
    "NIEDERER_LINE_NUM_POINTS",
    "LINE_POSTPROCESS_RELPATH",
    "POINTS_POSTPROCESS_RELPATH",
    "CACHE_POSTPROCESS_RELPATH",
    "LINE_POSTPROCESS_FUNCTION",
    "POINTS_POSTPROCESS_FUNCTION",
    "CACHE_POSTPROCESS_FUNCTION",
    "CASE_POSTPROCESS_CACHE_DIRNAME",
    "TABLE_SUMMARY_RELPATH",
    "EXCEL_REFERENCE_RELPATH",
]
