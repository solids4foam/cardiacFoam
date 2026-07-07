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
#     restitution_curves
#
# Description
#     Defines configuration template for restitution curve scenarios.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path

from .shared import (
    CONTROL_DICT_RELPATH,
    ELECTRO_PROPERTIES_RELPATH as SHARED_ELECTRO_PROPERTIES_RELPATH,
    OUTPUT_DIR_NAME,
    RUN_CASE_SCRIPT_RELPATH,
)
from ...ionic_model_catalog import IONIC_MODEL_CATALOG, planning_tissues


TUTORIAL_NAME = "restitutionCurves_s1s2Protocol"
CASE_DIR_NAME = "electrophysiologyProtocols/restitutionCurves_s1s2Protocol"
SETUP_DIR_NAME = "setup"

IONIC_MODELS = ("BuenoOrovio",)

IONIC_MODEL_TISSUE_MAP = {
    name: planning_tissues(entry)
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

# S1–S2 protocol parameters (intervals in milliseconds)
S1_INTERVAL_MS = 2000
N_S1 = 10
N_S2 = 2
S2_INTERVALS_MS = (
    1500, 1200, 1000, 800, 600, 400, 390, 380, 370, 360,
    350, 340, 330, 320, 310, 300, 290, 280, 270, 260, 250,
)

# Extra seconds appended to endTime after the last S2 beat
END_TIME_BUFFER_S = 2.0

# File paths
ELECTRO_PROPERTIES_SCOPE = "singleCellSolverCoeffs"
ELECTRO_PROPERTIES_RELPATH = SHARED_ELECTRO_PROPERTIES_RELPATH
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
OUTPUT_GLOB = "*.txt"
POSTPROCESS_SCRIPT_RELPATH = Path("postProcessing_restCurves.py")
POSTPROCESS_FUNCTION_NAME = "run_postprocessing"
TABLE_SUMMARY_RELPATH = Path("table_summary.py")

# Re-export shared paths used by the spec
__all__ = [
    "TUTORIAL_NAME",
    "CASE_DIR_NAME",
    "SETUP_DIR_NAME",
    "IONIC_MODELS",
    "IONIC_MODEL_TISSUE_MAP",
    "STIMULUS_MAP",
    "S1_INTERVAL_MS",
    "N_S1",
    "N_S2",
    "S2_INTERVALS_MS",
    "END_TIME_BUFFER_S",
    "ELECTRO_PROPERTIES_RELPATH",
    "CONTROL_DICT_RELPATH",
    "RUN_SCRIPT_RELPATH",
    "OUTPUT_DIR_NAME",
    "OUTPUT_GLOB",
    "POSTPROCESS_SCRIPT_RELPATH",
    "POSTPROCESS_FUNCTION_NAME",
    "TABLE_SUMMARY_RELPATH",
]
