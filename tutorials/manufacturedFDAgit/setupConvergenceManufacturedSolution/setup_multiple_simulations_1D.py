"""
================================================================================
Automated Electrophysiology Simulation Runner (:contentReference[oaicite:0]{index=0} model)
================================================================================

This script provides an automated pipeline to run and test cardiac electrophysiology 
simulations (1Dcable from :contentReference[oaicite:1]{index=1}) using the 
:contentReference[oaicite:2]{index=2} solver and post-process the results using 
:contentReference[oaicite:3]{index=3}.

It offers a GUI to select combinations of simulation parameters:
  - Δt (time step, ms)
  - Δx (spatial resolution, mm)
  - Tissue type (e.g., mCells, endocardialCells, epicardialCells)

Main workflow:
--------------
1. User selects parameter combinations from the GUI.
2. For each combination:
   - Update `blockMeshDict` with the corresponding mesh resolution.
   - Update `controlDict` with the chosen Δt.
   - Update `electroActivationProperties` with the selected tissue type.
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
- The GUI closes automatically after parameter selection.
- Modify CONFIG paths and parameter lists as needed for different environments.

================================================================================
"""

import subprocess
from itertools import product
import sys
import time
import pandas as pd
import tkinter as tk
from tkinter import ttk, messagebox
from pathlib import Path
from typing import List, Tuple, Optional
import tkinter as tk
from tkinter import ttk, messagebox
from itertools import product
import pandas as pd
from typing import Optional
import re
import shutil




def type_text(text: str, delay: float = 0.02) -> None:
    """Print text with a typing effect."""
    for char in text:
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)
    print()



def select_simulation_parameters(
    DX_VALUES, DT_VALUES, piecewise: bool = True) -> pd.DataFrame:
    """
    Generate simulation parameter combinations without GUI.

    Returns a DataFrame with columns:
    ["Δt[ms]", "Cells"]
    """
    if piecewise:
        combos = [
            (dt, dx)
            for (dt, dx) in zip(DT_VALUES, DX_VALUES)
        ]
        mode = "Piecewise (1-to-1 pairing)"
    else:
        combos = list(product(DT_VALUES, DX_VALUES))
        mode = "All combinations"

    df = pd.DataFrame(combos, columns=["Δt[ms]","Cells"])
    print(f"⚙️  Generated simulation combinations ({mode}):")
    print(df)
    return df


def update_blockMeshDict_dx(blockMeshDict_path: Path, dx: float) -> None:
    """
    Replace the hex block line in blockMeshDict according to dx value.
    Automatically generates the hex line based on dx.
    """
    lines = blockMeshDict_path.read_text().splitlines(keepends=True)
    replaced = False

    hex_line_template = "hex (0 1 2 3 4 5 6 7) ({dx} {dx} 1) simpleGrading (1 1 1)\n"

    with blockMeshDict_path.open('w') as f:
        for line in lines:
            stripped = line.strip()
            if not replaced and stripped.startswith('hex (0 1 2 3 4 5 6 7)') and not stripped.startswith('//'):
                f.write(hex_line_template.format(dx=int(dx)))
                replaced = True
            else:
                f.write(line)
    
def update_controlDict_dt(controlDict_path: Path, dt_ms: float) -> None:
    """
    Update the deltaT value in controlDict file with the given dt (in milliseconds).
    """
    dt_seconds = dt_ms / 1000.0  # convert ms to seconds
    lines = controlDict_path.read_text().splitlines(keepends=True)
    with controlDict_path.open('w') as f:
        for line in lines:
            if 'deltaT' in line and not line.strip().startswith('//'):
                f.write(f"deltaT          {dt_seconds};\n")
            else:
                f.write(line)

def run_simulation_cases(selected_combinations: pd.DataFrame, CONFIG) -> None:
    """Run the simulation cases based on selected parameters."""
    
    base_dir = CONFIG["case_file"].parent
    print(f"Base case directory: {base_dir}")

    controlDict_path = base_dir / "system" / "controlDict"
    blockMeshDict_path = base_dir / "system" / "blockMeshDict"
    run_cases_path = CONFIG["run_cases_script"]

    param_list = selected_combinations.values.tolist()
    total = len(param_list)

    for i, (dt_val, dx_val, ) in enumerate(param_list, start=1):
        print(f"🔄 Running simulation"
              f" ({i}/{total})", flush=True)

        update_blockMeshDict_dx(blockMeshDict_path, dx_val)
        update_controlDict_dt(controlDict_path, dt_val)
        # Run the case via shell script
        subprocess.run(
            ["bash", "-l", str(run_cases_path), str(base_dir)],
            check=True)
        
    print("✅ All simulations done.", flush=True)



