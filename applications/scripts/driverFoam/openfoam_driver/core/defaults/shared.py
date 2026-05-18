from __future__ import annotations

from pathlib import Path


CONTROL_DICT_RELPATH = Path("system/controlDict")
ELECTRO_PROPERTIES_RELPATH = Path("constant/electroProperties")
RUN_CASE_SCRIPT_RELPATH = Path("applications/scripts/driverFoam/openfoam_driver/scripts/run_case.sh")
OUTPUT_DIR_NAME = "postProcessing"
OUTPUT_RELPATH = Path(OUTPUT_DIR_NAME)
