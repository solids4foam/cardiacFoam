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
#     single_cell
#
# Description
#     Defines configuration template for single-cell scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path

from .shared import (
    ELECTRO_PROPERTIES_RELPATH as SHARED_ELECTRO_PROPERTIES_RELPATH,
    OUTPUT_DIR_NAME,
    RUN_CASE_SCRIPT_RELPATH,
)
from ...ionic_model_catalog import IONIC_MODEL_CATALOG


TUTORIAL_NAME = "singleCell"
CASE_DIR_NAME = "electrophysiologyProtocols/singleCell"
SETUP_DIR_NAME = "setupSingleCell"
IONIC_MODELS = tuple(
    name
    for name, entry in IONIC_MODEL_CATALOG.items()
    if "manufactured" not in entry.compatible_tissues
)

IONIC_MODEL_TISSUE_MAP = {
    name: entry.compatible_tissues
    for name, entry in IONIC_MODEL_CATALOG.items()
    if "manufactured" not in entry.compatible_tissues
}

STIMULUS_MAP = {}
for name, entry in IONIC_MODEL_CATALOG.items():
    if "manufactured" in entry.compatible_tissues:
        continue
    if name.startswith("Fabbri"):
        STIMULUS_MAP[name] = 0.0
    elif entry.model_type == "phenomenological":
        STIMULUS_MAP[name] = 0.4
    else:
        STIMULUS_MAP[name] = 60.0
ELECTRO_PROPERTIES_SCOPE = "singleCellSolverCoeffs"
ELECTRO_PROPERTIES_RELPATH = SHARED_ELECTRO_PROPERTIES_RELPATH
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
POSTPROCESS_SCRIPT_RELPATH = Path("singleCellinteractivePlots.py")
POSTPROCESS_FUNCTION_NAME = "run_postprocessing"
TABLE_SUMMARY_RELPATH = Path("table_summary.py")
OUTPUT_GLOB = "*.txt"
