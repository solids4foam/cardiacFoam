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




def type_text(text: str, delay: float = 0.02) -> None:
    """Print text with a typing effect."""
    for char in text:
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)
    print()



def select_simulation_parameters(
    TISSUE_TYPES, DX_VALUES, DT_VALUES, IONIC_MODELS,
    use_gui: bool = True, piecewise: bool = True
) -> Optional[pd.DataFrame]:
    # --- Skip GUI mode ---
    if not use_gui:
        if piecewise:
            combos = [
                (dt, dx, t, m)
                for (dt, dx) in zip(DT_VALUES, DX_VALUES)
                for t in TISSUE_TYPES
                for m in IONIC_MODELS
            ]
            mode = "Piecewise (1-to-1 pairing)"
        else:
            combos = list(product(DT_VALUES, DX_VALUES, TISSUE_TYPES, IONIC_MODELS))
            mode = "All combinations"

        df = pd.DataFrame(combos, columns=["Δt[ms]", "Cells", "Dimension", "IonicModel"])
        print(f"⚙️  Skipping GUI (use_gui=False) → returning {mode}:")
        print(df)
        return df

    # --- GUI mode ---
    selected_combinations = None

    def select_all(listbox: tk.Listbox) -> None:
        listbox.select_set(0, 'end')

    def run_simulation() -> None:
        nonlocal selected_combinations

        selected_t = [tissue_list.get(i) for i in tissue_list.curselection()]
        selected_dx = [float(dx_list.get(i)) for i in dx_list.curselection()]
        selected_dt = [float(dt_list.get(i)) for i in dt_list.curselection()]
        selected_im = [ionic_list.get(i) for i in ionic_list.curselection()]

        if not selected_t or not selected_dx or not selected_dt or not selected_im:
            messagebox.showwarning("Selection missing", "Please select at least one from each category.")
            return

        if piecewise:
            combos = [
                (dt, dx, t, im)
                for (dt, dx) in zip(selected_dt, selected_dx)
                for t in selected_t
                for im in selected_im
            ]
            mode = "Piecewise (1-to-1 pairing)"
        else:
            combos = list(product(selected_dt, selected_dx, selected_t, selected_im))
            mode = "All combinations"

        selected_combinations = pd.DataFrame(combos, columns=["Δt[ms]", "Cells", "Dimension", "IonicModel"])
        print(f"\nGenerated simulation cases ({mode}):")
        print(selected_combinations)

        messagebox.showinfo(
            "Simulation Combinations",
            f"Generated {len(selected_combinations)} simulations.\nMode: {mode}"
        )
        root.destroy()

    # --- GUI setup ---
    root = tk.Tk()
    root.title("Select Simulation Parameters")

    ttk.Label(root, text="Tissue Types:").grid(row=0, column=0, sticky="w")
    tissue_list = tk.Listbox(root, selectmode="multiple", exportselection=False, height=4)
    for t in TISSUE_TYPES:
        tissue_list.insert("end", t)
    tissue_list.grid(row=1, column=0, sticky="nsew")
    ttk.Button(root, text="Select All", command=lambda: select_all(tissue_list)).grid(row=2, column=0, pady=2)

    ttk.Label(root, text="Number of cells:").grid(row=0, column=1, sticky="w")
    dx_list = tk.Listbox(root, selectmode="multiple", exportselection=False, height=4)
    for x in DX_VALUES:
        dx_list.insert("end", x)
    dx_list.grid(row=1, column=1, sticky="nsew")
    ttk.Button(root, text="Select All", command=lambda: select_all(dx_list)).grid(row=2, column=1, pady=2)

    ttk.Label(root, text="Δt (ms):").grid(row=0, column=2, sticky="w")
    dt_list = tk.Listbox(root, selectmode="multiple", exportselection=False, height=4)
    for t in DT_VALUES:
        dt_list.insert("end", t)
    dt_list.grid(row=1, column=2, sticky="nsew")
    ttk.Button(root, text="Select All", command=lambda: select_all(dt_list)).grid(row=2, column=2, pady=2)

    ttk.Label(root, text="Ionic Model:").grid(row=0, column=3, sticky="w")
    ionic_list = tk.Listbox(root, selectmode="multiple", exportselection=False, height=4)
    for m in IONIC_MODELS:
        ionic_list.insert("end", m)
    ionic_list.grid(row=1, column=3, sticky="nsew")
    ttk.Button(root, text="Select All", command=lambda: select_all(ionic_list)).grid(row=2, column=3, pady=2)

    ttk.Button(root, text="Generate Combinations", command=run_simulation).grid(row=3, column=0, columnspan=4, pady=10)

    root.mainloop()
    return selected_combinations



def update_blockMeshDict_dx(blockMeshDict_path: Path, dx: float) -> None:
    """
    Replace the hex block line in blockMeshDict according to dx value.
    Automatically generates the hex line based on dx.
    """
    lines = blockMeshDict_path.read_text().splitlines(keepends=True)
    replaced = False

    hex_line_template = "hex (0 1 2 3 4 5 6 7) ({dx} 1 1) simpleGrading (1 1 1)\n"

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

def update_electro_model(electro_tissue_path: Path, model_flag: str) -> None:
    """
    Update the 'tissue' model type in electroActivationProperties.
    """
    lines = electro_tissue_path.read_text().splitlines(keepends=True)
    with electro_tissue_path.open('w') as f:
        for line in lines:
            if line.strip().startswith("tissue") and not line.strip().startswith("//"):
                f.write(f"tissue  {model_flag};\n")
            else:
                f.write(line)

def update_ionic_model(electro_tissue_path: Path, model_flag: str) -> None:
    """
    Update or insert the 'ionicModel' entry in electroActivationProperties.
    Preserves indentation and ignores commented lines.
    """
    lines = electro_tissue_path.read_text().splitlines(keepends=True)
    with electro_tissue_path.open('w') as f:
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("ionicModel ") and not stripped.startswith("//"):
                indent = line[:len(line) - len(line.lstrip())]
                f.write(f"{indent}ionicModel  {model_flag};\n")
            else:
                f.write(line)


def run_simulation_cases(selected_combinations: pd.DataFrame, CONFIG) -> None:
    """Run the simulation cases based on selected parameters."""
    
    base_dir = CONFIG["case_file"].parent
    print(f"Base case directory: {base_dir}")
    output_dir = base_dir / CONFIG["model_name"]
    output_dir.mkdir(exist_ok=True)

    controlDict_path = base_dir / "system" / "controlDict"
    blockMeshDict_path = base_dir / "system" / "blockMeshDict"
    tissue_path = base_dir / "constant" / "electroActivationProperties"
    run_cases_path = CONFIG["run_cases_script"]

    param_list = selected_combinations.values.tolist()
    total = len(param_list)

    for i, (dt_val, dx_val, tissue_val, ionicModel_val) in enumerate(param_list, start=1):
        print(f"🔄 Running simulation for the {tissue_val} cell model "
              f"with Δt={dt_val} ms, nCells:{dx_val} ({i}/{total})", flush=True)

        update_blockMeshDict_dx(blockMeshDict_path, dx_val)
        update_controlDict_dt(controlDict_path, dt_val)
        update_electro_model(tissue_path, tissue_val)
        update_ionic_model(tissue_path, ionicModel_val)

        # Run the case via shell script
        subprocess.run(
            ["bash", "-l", str(run_cases_path), str(base_dir)],
            check=True
        )

        # Ensure case.foam exists for ParaView
        case_foam_path = base_dir / "case.foam"
        case_foam_path.touch()

        # Run the ParaView pvpython script
        subprocess.run([
            str(CONFIG["pvpython_path"]),
            str(CONFIG["paraview_script_path"]),
            "--case", str(CONFIG["case_file"]),
            "--outdir", str(output_dir),
            "--dx", str(dx_val),
            "--dt", str(dt_val),
            "--model", tissue_val,
            "--ionicModel", ionicModel_val
        ], check=True)

    print("✅ All simulations done.", flush=True)


