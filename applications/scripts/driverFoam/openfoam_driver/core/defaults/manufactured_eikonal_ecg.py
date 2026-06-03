from __future__ import annotations

from pathlib import Path

from .shared import OUTPUT_DIR_NAME, RUN_CASE_SCRIPT_RELPATH


TUTORIAL_NAME = "manufacturedEikonalECG"
CASE_DIR_NAME = "manufacturedSolutions/eikonalECG"
SETUP_DIR_NAME = "setupManufacturedEikonalECG"
NUMBER_CELLS = (10, 20, 40, 80)
DIMENSIONS = ("1D", "2D", "3D")
SOLVER_TYPES = ("eikonal",)
PIECEWISE_SWEEP = True
ELECTRO_PROPERTIES_SCOPE = "eikonalSolverCoeffs"
BLOCK_MESH_DICT_TEMPLATE = "system/blockMeshDict.{dimension}"
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
POSTPROCESS_SCRIPT_RELPATH = Path("post_processing_manufactured_eikonal_ecg.py")
POSTPROCESS_FUNCTION_NAME = "run_postprocessing"
RUN_IN_PARALLEL = True
VERIFICATION_MODEL_TYPE = "manufacturedEikonalVerifier"
ECG_REFERENCE_QUADRATURE_ORDER = 96
ECG_CHECK_QUADRATURE_ORDERS = (6, 12, 24, 48)
ECG_ELECTRODES_BY_DIMENSION = {
    "1D": {
        "E1": "(-0.5 0 0)",
        "E2": "(1.5 0 0)",
        "E3": "(1.2 0 0)",
        "E4": "(1.35 0 0)",
        "E5": "(1.65 0 0)",
    },
    "2D": {
        "E1": "(-0.5 0.5 0)",
        "E2": "(1.5 0.5 0)",
        "E3": "(1.2 0.23 0)",
        "E4": "(1.35 0.78 0)",
        "E5": "(0.18 1.35 0)",
    },
    "3D": {
        "E1": "(-0.5 0.5 0.5)",
        "E2": "(1.5 0.5 0.5)",
        "E3": "(1.2 0.23 0.61)",
        "E4": "(1.35 0.74 0.28)",
        "E5": "(1.55 0.41 0.83)",
    },
}
BLOCK_MESH_RESOLUTION_BY_DIMENSION = {
    "1D": "{cells} 1 1",
    "2D": "{cells} {cells} 1",
    "3D": "{cells} {cells} {cells}",
}
