"""
This script sets up and analyzes convergence tests for simulations 
using a manufactured solution for verification of the Eletrophysiology solver.
The user can define a range of spatial resolutions (Δx) and time steps (Δt) and the dimension.
To note that the explicit solver uses a CFL of 0.1 and DT is automatic.

For more information: Pras Pathmanathan, Richard A. Gray et al 2014: 
Verification of omputational models for cardiac electrophysiology,
 
DOI: https://doi.org/10.1002/cnm.2615


Usage:
------
- Run this script directly in the folder.
- Select your desired Δt, Δx, and dimension type in the CONFIG {} below. 
- Simulations will be executed automatically, and results will be saved in .dat files.
- With post processing, results will be organized and copied to a designated folder in the case directory.
"""
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd


#!/usr/bin/env python3
"""
Master controller for Manufactured-Solution Convergence Testing.

Provides a command-line interface similar to main_setup_Niederer:

    python main_setup_convergence.py all
    python main_setup_convergence.py sim
    python main_setup_convergence.py post

Optional flags:

    --case PATH           Override case folder
    --dimension 1D 2D 3D  Limit to specific dimensions

This script orchestrates:
  ✓ Simulation parameter selection
  ✓ Dictionary modification
  ✓ Launching OpenFOAM runs
  ✓ Automatic error extraction & convergence plotting
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from pathlib import Path


from setup_multiple_simulations_manufactured_solution import (
    select_simulation_parameters,
    run_simulation_cases,
)
from post_processing_manufactured import (
    organize_dat_files,
    read_error_dat_files,
    compute_convergence_rates,
    plot_errors_implicit_explicit,
    plot_Vm_across_dimensions,
)

# --------------------- CONFIGURATION ----------------------------

DEFAULT_CONFIG = {
    "case_file": Path.cwd().parent / "case.foam",
    "run_cases_script": Path.cwd() / "run_cases.sh",
}

NUMBER_CELLS = [10, 20, 40, 80, 160, 320]
DT_VALUES = [
    0.00892857, 0.00224215, 0.000560538,
    0.000140174, 3.50471e-05, 8.76201e-06
]

DIMENSION_ALL = ["1D", "2D", "3D"]
SOLVER_TYPES = ["implicit"]


BASE_DIR = Path(__file__).resolve().parent.parent
OUTPUT_FOLDER = BASE_DIR / "error_manufactured_all_dim"


# ================================================================
# --------------------- SIMULATION STAGE --------------------------
# ================================================================

def stage_simulations(case, dimension_override=None):
    print("\n=== Running MANUFACTURED SOLUTION SIMULATIONS ===\n")

    # Determine dimensions to run
    if dimension_override:
        dimensions = [dimension_override]
        print(f"⮞ Restricting to dimension: {dimension_override}")
    else:
        dimensions = DIMENSION_ALL
        print("⮞ Running ALL dimensions:", DIMENSION_ALL)

    CONFIG = DEFAULT_CONFIG.copy()
    CONFIG["case_file"] = Path(case) if case else DEFAULT_CONFIG["case_file"]

    # Select dt/dx combinations
    selected_combinations = select_simulation_parameters(
        NUMBER_CELLS, DT_VALUES, piecewise=True
    )

    if selected_combinations is None or selected_combinations.empty:
        print("❌ No simulation parameters selected. Aborting simulation stage.")
        return

    # Run all cases
    run_simulation_cases(
        selected_combinations,
        CONFIG,
        dimensions,
        SOLVER_TYPES
    )

    print("\n✔ Simulation stage complete.\n")


# --------------------- POSTPROCESSING STAGE ---------------------
def stage_postprocess(dimension_override=None, folder_name=None):
    print("\n=== MANUFACTURED SOLUTION POST-PROCESSING ===\n")

    if folder_name is None:
        print("❌ ERROR: No output folder specified.")
        return

    # Convert to absolute path for debugging clarity
    folder = Path(folder_name).resolve()
    print(f"[DEBUG] Initial folder path (resolved): {folder}")

    # Organize .dat files according to your existing logic
    organize_dat_files(folder)

    # Print full folder path to verify location
    print(f"[DEBUG] After organize_dat_files, using folder: {folder}")

    # Now read the .dat files (also prints its own debug messages)
    print(f"⮞ Reading .dat files from: {folder}")

    df = read_error_dat_files(folder)

    if df is None or df.empty:
        print(f"❌ No error files found in: {folder}")
        return

    # Dimension filtering
    if dimension_override:
        df = df[df["Dimension"] == dimension_override]
        if df.empty:
            print(f"⚠ No data found for dimension {dimension_override} in folder: {folder}")
            return
        print(f"⮞ Filtering only dimension: {dimension_override}")

    # Compute convergence
    print("\n📉 Convergence rates:")
    print(compute_convergence_rates(df))

    # Plot Vm across dimensions
    print("\n📊 Plotting Vm error across dimensions\n")
    plot_Vm_across_dimensions(df)

    # Plot explicit vs implicit
    print("\n📊 Plotting implicit vs explicit per dimension\n")
    plot_errors_implicit_explicit(df)

    print("\n✔ Post-processing complete.\n")



# ------------------------- CLI HANDLER ---------------------------


def parse_args():
    parser = argparse.ArgumentParser(description="Manufactured-solution convergence pipeline")
    parser.add_argument("action", choices=["sim", "post", "all"],
                        help="Choose: sim, post, or all")
    parser.add_argument("--case", type=str, default=None,
                        help="Override OpenFOAM case directory")
    parser.add_argument("--dimension", type=str, choices=["1D", "2D", "3D"],
                        default=None,
                        help="Run only for a specific dimension")

    return parser.parse_args()

def main():
    args = parse_args()

    print("\n--- MANUFACTURED SOLUTION CONTROLLER ---\n")
    print(f"Action    : {args.action}")
    print(f"Case      : {args.case or 'default'}")
    print(f"Dimension : {args.dimension or 'all'}\n")

    if args.action == "sim":
        stage_simulations(args.case, args.dimension)

    elif args.action == "post":
        stage_postprocess(args.dimension, OUTPUT_FOLDER)

    elif args.action == "all":
        stage_simulations(args.case, args.dimension)
        stage_postprocess(args.dimension, OUTPUT_FOLDER)

    print("\n🎉 Completed successfully.\n")


if __name__ == "__main__":
    main()
