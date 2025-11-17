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



from setup_multiple_simulations_manufactured_solution import select_simulation_parameters, run_simulation_cases
from post_processing_manufactured import organize_dat_files, read_error_dat_files, compute_convergence_rates, plot_errors, plot_errors_implicit_explicit, plot_Vm_across_dimensions




# Configuration constants
CONFIG = {
    "case_file": Path.cwd().parent / "case.foam",
    "run_cases_script": Path.cwd() / "run_cases1D.sh",
}
      
    
DX_VALUES = [10,20,40,80,160,320]# number of cells
#DX_VALUES = [10,20,40,80,160]# number of cells
#DT_VALUES = [1,1,1,1,1,1] # Define for implicit solver only. Keep same size list for explicit
#DT_VALUES = [1e-3, 1e-3,5e-4,5e-4,2e-4,1e-4]  # ms
DT_VALUES = [0.00892857,0.00224215,0.000560538, 0.000140174, 3.50471e-05, 8.76201e-06]  # s
#DT_VALUES = [0.00892857,0.00224215,0.000560538, 0.000140174, 3.50471e-05]  # s


TISSUE_TYPES = [ "2D"] 
FOLDER_NAME = [ "error_manufactured_2D_implicit_ODE"] 
"""
Update the Allrun for the specific blockMesh.{dimension} you want for now.
"""


def main():
    selected_combinations = select_simulation_parameters(DX_VALUES, DT_VALUES, piecewise=True)
    if selected_combinations is None or selected_combinations.empty:
        print("No simulation parameters selected. Exiting.")
        return

    run_simulation_cases(selected_combinations, CONFIG, TISSUE_TYPES)

    dfs = []  # store DataFrames for 1D, 2D, 3D.

    for tissue_type, folder_name in zip(TISSUE_TYPES, FOLDER_NAME):
        print(f"\n=== Processing {tissue_type} ({folder_name}) ===")

        organize_dat_files(folder_name)
        df = read_error_dat_files(folder_name)
        dfs.append(df)

        if df is None or df.empty:
            print(f"No data found for {tissue_type}. Skipping.")
            continue

        print(df)

        # Individual plots per tissue type
        #plot_errors(df, solver_type="Explicit")
        #plot_errors(df, solver_type="Implicit")
        plot_errors_implicit_explicit(df, tissue_type)

        # Convergence rates
        rates = compute_convergence_rates(df)
        print(f"\nConvergence rates for {tissue_type}:")
        print(rates)

    # --- Combined Vm plot across tissue types ---
    print("\n=== Plotting Vm comparison across tissues ===")
    plot_Vm_across_dimensions(dfs, TISSUE_TYPES)

        


if __name__ == "__main__":
    main()
