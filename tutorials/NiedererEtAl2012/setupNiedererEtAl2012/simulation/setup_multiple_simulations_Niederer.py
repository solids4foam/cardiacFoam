"""
================================================================================
Automated Electrophysiology Simulation Runner (:contentReference[oaicite:0]{index=0} model)
================================================================================

This script provides an automated pipeline to run and test cardiac electrophysiology 
simulations (slab model from :contentReference[oaicite:1]{index=1}) using the 
:contentReference[oaicite:2]{index=2} solver and post-process the results using 
:contentReference[oaicite:3]{index=3}.

Main workflow:
--------------
1. User selects parameter combinations from the CONFIG.
2. For each combination:
   - Update `blockMeshDict` with the corresponding mesh resolution.
   - Update `controlDict` with the chosen Δt.
   - Update `cardiacProperties` with the selected tissue type.
3. Runs the case using a helper shell script (`run_cases.sh`).
4. Exports activation times along a line and specific points via 
   `Export_Paraview_Niederer.py` using `pvpython`.
5. Results are saved in a folder named after the model (default: "Niederer_Activation_time").

Requirements:
-------------
- :contentReference[oaicite:4]{index=4} installed and accessible via `run_cases.sh`
- :contentReference[oaicite:5]{index=5} installed with `pvpython` available
- The following auxiliary files must be in the same directory as this script:
    - `Export_Paraview_Niederer.py`
    - `run_cases.sh`
    - `Niederer_Points.txt`
- Python packages: pandas, tkinter

Usage:
------
- Run this script directly.
- Select your desired Δt, Δx, and tissue type values in the popup GUI.
- Simulations will be executed automatically, and results will be saved.

Note:
-----
- Modify CONFIG paths and parameter lists as needed for different environments.

================================================================================
"""


import subprocess
from itertools import product
import sys
import time
import pandas as pd
from pathlib import Path
import re

# ----------------------------------------------------------
# 🔺 NEW: Import ALL configuration from main.py
# ----------------------------------------------------------
import sys
from pathlib import Path

SETUP_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SETUP_ROOT))

from main_setup_Niederer import CONFIG, DT_VALUES, DX_VALUES, TISSUE_TYPES, IONIC_MODEL, SOLVER, OUTPUT_FOLDER



MESSAGE = (
    "This is an automatic script to run several electrophysiology simulations for "
    "the Niederer et al. 2012 slab benchmark.\n"
    "The user can run the simulations for the spatial and temporal parameters in the paper "
    "and others. (Δt, Δx, tissue type)\n"
    "Make sure you have the auxiliary files in the same directory as this script. Enjoy 😊\n"
    "NOTE: code can be expanded to imported meshes that DO NOT USE blockMeshDict.\n"
)


def type_text(text: str, delay: float = 0.005) -> None:
    """Print text with a typing effect."""
    for char in text:
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)
    print()


def selected_combinations_combo():
    combos = list(product(
        IONIC_MODEL, TISSUE_TYPES, DT_VALUES, DX_VALUES, SOLVER
    ))

    return pd.DataFrame(combos, columns=[
        "ionicModel", "tissueType", "dt", "dx", "solver"
    ])



def update_blockMeshDict_dx(blockMeshDict_path: Path, dx: float) -> None:
    """Replace mesh resolution line according to dx."""
    hex_lines = {
        0.5: "hex (0 1 2 3 4 5 6 7) (40 6 14) simpleGrading (1 1 1)\n",
        0.2: "hex (0 1 2 3 4 5 6 7) (100 15 35) simpleGrading (1 1 1)\n",
        0.1: "hex (0 1 2 3 4 5 6 7) (200 30 70) simpleGrading (1 1 1)\n",
    }

    keys = sorted(hex_lines.keys())
    closest_dx = min(keys, key=lambda k: abs(k - dx))

    lines = blockMeshDict_path.read_text().splitlines(keepends=True)
    replaced = False
    with blockMeshDict_path.open('w') as f:
        for line in lines:
            stripped = line.strip()
            if not replaced and stripped.startswith('hex (0 1 2 3 4 5 6 7)') and not stripped.startswith('//'):
                f.write(hex_lines[closest_dx])
                replaced = True
            else:
                f.write(line)



def update_tissue(electro_tissue_path: Path, tissue_flag: str) -> None:
    """Update the tissue model type in cardiacPropertiesDict."""
    lines = electro_tissue_path.read_text().splitlines(keepends=True)
    with electro_tissue_path.open('w') as f:
        for line in lines:
            if line.strip().startswith("tissue") and not line.strip().startswith("//"):
                f.write(f"tissue  {tissue_flag};\n")
            else:
                f.write(line)



def update_controlDict_dt(controlDict_path: Path, dt_ms: float) -> None:
    """Update deltaT in controlDict."""
    dt_seconds = dt_ms / 1000.0
    lines = controlDict_path.read_text().splitlines(keepends=True)
    with controlDict_path.open('w') as f:
        for line in lines:
            if 'deltaT' in line and not line.strip().startswith('//'):
                f.write(f"deltaT          {dt_seconds};\n")
            else:
                f.write(line)



def update_controlDict_endTime(controlDict_path: Path, dx: float) -> None:
    endTime_map = {
        0.5: 0.2,
        0.2: 0.080,
        0.1: 0.055,
    }

    keys = sorted(endTime_map.keys())
    closest_dx = min(keys, key=lambda k: abs(k - dx))
    endTime_val = endTime_map[closest_dx]

    lines = controlDict_path.read_text().splitlines(keepends=True)
    pattern = re.compile(r"^\s*endTime\s+\S+")

    with controlDict_path.open("w") as f:
        replaced = False
        for line in lines:
            stripped = line.strip()

            if pattern.match(line) and not stripped.startswith("//"):
                indent = line[:len(line) - len(line.lstrip())]
                f.write(f"{indent}endTime    {endTime_val};\n")
                replaced = True
            else:
                f.write(line)

        if not replaced:
            f.write(f"\nendTime    {endTime_val};\n")



def update_ionic_model(electro_tissue_path: Path, model_flag: str) -> None:
    """Update ionic model."""
    lines = electro_tissue_path.read_text().splitlines(keepends=True)
    with electro_tissue_path.open('w') as f:
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("ionicModel") and not stripped.startswith("//"):
                indent = line[:len(line) - len(line.lstrip())]
                f.write(f"{indent}ionicModel  {model_flag};\n")
            else:
                f.write(line)


def update_solver(numericalPropertiesPath: Path, solverType: str) -> None:
    print("SOLVER VAL = ", solverType)
    # Decide flag values
    if solverType == "explicit":
        solveExplicit_flag = "yes"
    else:
        solveExplicit_flag = "no"

    lines = numericalPropertiesPath.read_text().splitlines(keepends=True)

    with numericalPropertiesPath.open('w') as f:
        for line in lines:
            stripped = line.strip()
            # Update solveExplicit
            if stripped.startswith("solveExplicit") and not stripped.startswith("//"):
                indent = line[:len(line) - len(line.lstrip())]
                f.write(f"{indent}solveExplicit    {solveExplicit_flag};\n")
                continue
            else:
                f.write(line)

def run_simulation_cases(selected_combinations: pd.DataFrame) -> None:
    base_dir = Path(CONFIG["case_root"])
    output_dir = base_dir / OUTPUT_FOLDER
    output_dir.mkdir(exist_ok=True)

    controlDict_path = base_dir / "system" / "controlDict"
    blockMeshDict_path = base_dir / "system" / "blockMeshDict"
    tissue_path = base_dir / "constant" / "cardiacProperties"
    solver_path = base_dir / "constant" / "timeIntegrationProperties"
    run_cases_path = CONFIG["run_cases_script"]


    param_list = selected_combinations.values.tolist()
    total = len(param_list)

    for i, (ionic_val, tissue_val, dt_val, dx_val, solver_val) in enumerate(param_list, start=1):

        print(
            f"🔄 Running simulation:{solver_val}, ionicModel={ionic_val}, tissue={tissue_val}, "
            f"Δt={dt_val} ms, Δx={dx_val} mm ({i}/{total})",
            flush=True
        )

        update_blockMeshDict_dx(blockMeshDict_path, dx_val)
        update_controlDict_dt(controlDict_path, dt_val)
        update_controlDict_endTime(controlDict_path, dx_val)
        update_tissue(tissue_path, tissue_val)
        update_ionic_model(tissue_path, ionic_val)
        update_solver(solver_path, solver_val)


        subprocess.run(["bash", "-l", str(run_cases_path), str(base_dir)], check=True)

        (base_dir / "case.foam").touch()

        subprocess.run([
        str(CONFIG["pvpython_path"]),
        str(CONFIG["paraview_script_path"]),
        "--case", str(CONFIG["case_root"] / "case.foam"),
        "--points", str(CONFIG["points_file"]),
        "--outdir", str(OUTPUT_FOLDER),
        "--dx", str(dx_val),
        "--dt", str(dt_val),
        "--tissue", tissue_val,
        "--ionicModel", ionic_val,
        "--solver", solver_val
    ], check=True)


    print("✅ All simulations done.", flush=True)



def main():
    type_text(MESSAGE)
    selected_combinations = selected_combinations_combo()

    if selected_combinations is None or selected_combinations.empty:
        print("No simulation parameters selected. Exiting.")
        return

    run_simulation_cases(selected_combinations)



if __name__ == "__main__":
    main()
