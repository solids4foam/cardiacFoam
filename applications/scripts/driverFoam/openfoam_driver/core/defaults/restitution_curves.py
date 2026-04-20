from __future__ import annotations

from pathlib import Path

from .shared import (
    CONTROL_DICT_RELPATH,
    ELECTRO_PROPERTIES_RELPATH as SHARED_ELECTRO_PROPERTIES_RELPATH,
    OUTPUT_DIR_NAME,
    RUN_CASE_SCRIPT_RELPATH,
)


TUTORIAL_NAME = "restitutionCurves_s1s2Protocol"
CASE_DIR_NAME = "singleCellprotocols/restitutionCurves_s1s2Protocol"
SETUP_DIR_NAME = "setupRestitutionCurves_s1s2Protocol"

IONIC_MODELS = ("BuenoOrovio",)

IONIC_MODEL_TISSUE_MAP = {
    "TNNP": ("mCells", "endocardialCells", "epicardialCells"),
    "BuenoOrovio": ("mCells", "endocardialCells", "epicardialCells"),
    "Gaur": ("myocyte",),
    "Courtemanche": ("myocyte",),
}

STIMULUS_MAP = {
    "TNNP": 60.0,
    "BuenoOrovio": 0.4,
    "Gaur": 60.0,
    "Courtemanche": 70.0,
}

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
TABLE_SUMMARY_RELPATH = Path("postProcessing/table_summary.py")

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
