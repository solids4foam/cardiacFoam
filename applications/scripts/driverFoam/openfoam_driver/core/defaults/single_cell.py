from __future__ import annotations

from pathlib import Path

from .shared import (
    ELECTRO_PROPERTIES_RELPATH as SHARED_ELECTRO_PROPERTIES_RELPATH,
    OUTPUT_DIR_NAME,
    RUN_CASE_SCRIPT_RELPATH,
)


TUTORIAL_NAME = "singleCell"
CASE_DIR_NAME = "singleCellprotocols/singleCell"
SETUP_DIR_NAME = "setupSingleCell"
IONIC_MODELS = ("TNNP", "Gaur", "Courtemanche", "BuenoOrovio")
IONIC_MODEL_TISSUE_MAP = {
    "TNNP": ("epicardialCells", "mCells", "endocardialCells"),
    "BuenoOrovio": ("epicardialCells", "mCells", "endocardialCells"),
    "Gaur": ("myocyte",),
    "Courtemanche": ("myocyte",),
}
STIMULUS_MAP = {
    "TNNP": 60,
    "BuenoOrovio": 0.4,
    "Gaur": 60,
    "Courtemanche": 60,
}
ELECTRO_PROPERTIES_SCOPE = "singleCellSolverCoeffs"
ELECTRO_PROPERTIES_RELPATH = SHARED_ELECTRO_PROPERTIES_RELPATH
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
POSTPROCESS_SCRIPT_RELPATH = Path("singleCellinteractivePlots.py")
POSTPROCESS_FUNCTION_NAME = "run_postprocessing"
TABLE_SUMMARY_RELPATH = Path("postProcessing/table_summary.py")
OUTPUT_GLOB = "*.txt"
