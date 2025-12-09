# main_simulation.py

from pathlib import Path
import re
import numpy as np
import sys
import matplotlib.pyplot as plt


from setup_multiple_simulations_1D import type_text, select_simulation_parameters, run_simulation_cases
# from postProcessing_convergence import (
#     load_csv_files, clean_line_data, get_available_dts_ncells,
#     generate_color_palette, plot_all_by_ncells, plot_by_dt,
#     clean_points_data, organize_points_by_id, plot_activation_2d_with_slopes,
#     base_colors
# )



# Configuration constants
CONFIG = {
    "pvpython_path": Path("/Applications/ParaView-5.13.3.app/Contents/bin/pvpython"),
    "paraview_script_path": Path("Export_Paraview_1D.py"),
    "case_file": Path.cwd().parent / "case.foam",
    "model_name": "Test_1Dconduction",
    "run_cases_script": Path.cwd() / "run_cases1D.sh",
}

         
#DX_VALUES = [20]    # number of cells
#DT_VALUES = [0.01]  # ms     
DX_VALUES = [200 ]# number of cells
DT_VALUES = [0.25]  # ms
TISSUE_TYPES = ["mCells"]  #["mCells", "Endocardial", "Epicardial", "myocyte"] 
IONIC_MODEL = ["BuenoOrovio", "TNNP"]






def main():
    # --- Run simulations ---
    selected_combinations = select_simulation_parameters(TISSUE_TYPES, DX_VALUES, DT_VALUES, IONIC_MODEL, use_gui= False, piecewise=True)
    if selected_combinations is None or selected_combinations.empty:
        print("No simulation parameters selected. Exiting.")
        return

    run_simulation_cases(selected_combinations, CONFIG)

    # --- Load results based on the simulation configuration ---
    print(f"Processing results for model: {CONFIG['model_name']}")

if __name__ == "__main__":
    main()
