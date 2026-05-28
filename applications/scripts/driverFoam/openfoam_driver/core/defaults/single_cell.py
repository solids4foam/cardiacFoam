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
IONIC_MODELS = (
    "AlievPanfilov",
    "BuenoOrovio",
    "Courtemanche",
    "Fabbri",
    "Gaur",
    "Grandi",
    "PerisYague",
    "Stewart",
    "TNNP",
    "ToRORd_dynCl",
    "Trovato",
)
IONIC_MODEL_TISSUE_MAP = {
    "AlievPanfilov": ("myocyte",),
    "BuenoOrovio": ("epicardialCells", "mCells", "endocardialCells"),
    "Courtemanche": ("myocyte",),
    "Fabbri": ("myocyte",),
    "Gaur": ("myocyte",),
    "Grandi": ("myocyte",),
    "PerisYague": ("myocyte",),
    "Stewart": ("myocyte",),
    "TNNP": ("epicardialCells", "mCells", "endocardialCells"),
    "ToRORd_dynCl": ("epicardialCells", "mCells", "endocardialCells"),
    "Trovato": ("myocyte",),
}
STIMULUS_MAP = {
    "AlievPanfilov": 0.4,
    "BuenoOrovio": 0.4,
    "Courtemanche": 65,
    "Fabbri": 0,
    "Gaur": 65,
    "Grandi": 65,
    "PerisYague": 65,
    "Stewart": 65,
    "TNNP": 65,
    "ToRORd_dynCl": 65,
    "Trovato": 65,
}
ELECTRO_PROPERTIES_SCOPE = "singleCellSolverCoeffs"
ELECTRO_PROPERTIES_RELPATH = SHARED_ELECTRO_PROPERTIES_RELPATH
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
POSTPROCESS_SCRIPT_RELPATH = Path("singleCellinteractivePlots.py")
POSTPROCESS_FUNCTION_NAME = "run_postprocessing"
TABLE_SUMMARY_RELPATH = Path("postProcessing/table_summary.py")
OUTPUT_GLOB = "*.txt"
