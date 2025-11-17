"""
The setup and analysis of convergence tests for simulations using a manufactured solution
for verification of the Electrophysiology solver with the input of dx, dt and dimension.

Main workflow:
--------------
1. User selects parameter combinations from the GUI.
2. For each combination:
   - Update `blockMeshDict."dimension"` with the corresponding mesh resolution.
   - Update `controlDict` with the chosen Δt.
   - Update `electroActivationProperties` with the selected tissue type and conductivity.
3. Runs the case using a helper shell script (`run_cases.sh`).
4. Solves the errors in the solvers and outputs .dat file.

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


def update_blockMeshDict_dx(blockMeshDict_path: Path, dx: float, tissue_type: str) -> None:
    """
    Replace the hex block line in blockMeshDict according to dx value.
    Automatically generates the hex line based on dx.
    """
    lines = blockMeshDict_path.read_text().splitlines(keepends=True)
    replaced = False

    if tissue_type == "1D":
        hex_line_template = "hex (0 1 2 3 4 5 6 7) ({dx} 1 1) simpleGrading (1 1 1)\n"
    elif tissue_type == "2D":
        hex_line_template = "hex (0 1 2 3 4 5 6 7) ({dx} {dx} 1) simpleGrading (1 1 1)\n"
    elif tissue_type == "3D":
        hex_line_template = "hex (0 1 2 3 4 5 6 7) ({dx} {dx} {dx}) simpleGrading (1 1 1)\n"
    else:
        raise ValueError(f"Unsupported tissue type: {tissue_type}")

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
    dt_seconds = dt_ms   # convert ms to seconds
    lines = controlDict_path.read_text().splitlines(keepends=True)
    with controlDict_path.open('w') as f:
        for line in lines:
            if 'deltaT' in line and not line.strip().startswith('//'):
                f.write(f"deltaT          {dt_seconds};\n")
            else:
                f.write(line)


def update_dimension(electro_tissue_path: Path, dimension: str) -> None:
    """
    Update or insert the "dimension" entry in electroActivationProperties.
    Preserves indentation and ignores commented lines.
    """
    lines = electro_tissue_path.read_text().splitlines(keepends=True)
    with electro_tissue_path.open('w') as f:
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("dimension ") and not stripped.startswith("//"):
                indent = line[:len(line) - len(line.lstrip())]
                f.write(f"{indent}dimension \"{dimension}\";\n")
            else:
                f.write(line)



def run_simulation_cases(selected_combinations: pd.DataFrame, CONFIG, TISSUE_TYPES) -> None:
    """Run the simulation cases based on selected parameters."""
    
    base_dir = CONFIG["case_file"].parent
    print(f"Base case directory: {base_dir}")

    tissue_type = TISSUE_TYPES[0]
    print (tissue_type)
    blockMeshDict_path = base_dir / "system" / f"blockMeshDict.{tissue_type}"
    controlDict_path = base_dir / "system" / "controlDict"
    tissue_path = base_dir / "constant" / "electroActivationProperties"
    run_cases_path = CONFIG["run_cases_script"]

    param_list = selected_combinations.values.tolist()
    total = len(param_list)

    for i, (dt_val, dx_val) in enumerate(param_list, start=1):
        print(f"🔄 Running simulation ({i}/{total})", flush=True)

        update_blockMeshDict_dx(blockMeshDict_path, dx_val, tissue_type)
        update_controlDict_dt(controlDict_path, dt_val)
        update_dimension(tissue_path, tissue_type)

        subprocess.run(
            ["bash", "-l", str(run_cases_path), str(base_dir)],
            check=True
        )

    print("✅ All simulations done.", flush=True)
