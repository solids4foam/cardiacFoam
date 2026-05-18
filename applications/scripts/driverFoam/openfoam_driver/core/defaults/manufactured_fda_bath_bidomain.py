from __future__ import annotations

from pathlib import Path

from .manufactured_fda import (
    BLOCK_MESH_RESOLUTION_BY_DIMENSION,
    DT_VALUES,
    DIMENSIONS,
    NUMBER_CELLS,
    OUTPUT_DIR_NAME,
    PIECEWISE_SWEEP,
    RUN_CASE_SCRIPT_RELPATH,
)


TUTORIAL_NAME = "manufacturedFDABathBidomain"
CASE_DIR_NAME = "manufacturedSolutions/bathBidomain"
SETUP_DIR_NAME = "setupManufacturedFDA"
SOLVER_TYPES = ("implicit",)
ELECTRO_PROPERTIES_SCOPE = "bidomainSolverCoeffs"
BLOCK_MESH_DICT_TEMPLATE = "system/blockMeshDict.{dimension}"
RUN_SCRIPT_RELPATH = RUN_CASE_SCRIPT_RELPATH
POSTPROCESS_SCRIPT_RELPATH = Path("post_processing_manufactured_bath.py")
POSTPROCESS_FUNCTION_NAME = "run_postprocessing"
RUN_IN_PARALLEL = True
VERIFICATION_MODEL_TYPE = "manufacturedFDABathBidomainVerifier"

__all__ = [
    "BLOCK_MESH_DICT_TEMPLATE",
    "BLOCK_MESH_RESOLUTION_BY_DIMENSION",
    "CASE_DIR_NAME",
    "DT_VALUES",
    "DIMENSIONS",
    "ELECTRO_PROPERTIES_SCOPE",
    "NUMBER_CELLS",
    "OUTPUT_DIR_NAME",
    "PIECEWISE_SWEEP",
    "POSTPROCESS_FUNCTION_NAME",
    "POSTPROCESS_SCRIPT_RELPATH",
    "RUN_IN_PARALLEL",
    "RUN_SCRIPT_RELPATH",
    "SETUP_DIR_NAME",
    "SOLVER_TYPES",
    "TUTORIAL_NAME",
    "VERIFICATION_MODEL_TYPE",
]
