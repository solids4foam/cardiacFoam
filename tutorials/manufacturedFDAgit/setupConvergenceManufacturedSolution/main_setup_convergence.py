# main_simulation.py

from pathlib import Path
import re
import numpy as np
import sys
import matplotlib.pyplot as plt
import glob
import os
import shutil
import pandas as pd
from math import sqrt


from setup_multiple_simulations_1D import type_text, select_simulation_parameters, run_simulation_cases
# from postProcessing_convergence import (
#     load_csv_files, clean_line_data, get_available_dts_ncells,
#     generate_color_palette, plot_all_by_ncells, plot_by_dt,
#     clean_points_data, organize_points_by_id, plot_activation_2d_with_slopes,
#     base_colors
# )



# Configuration constants
CONFIG = {
    "case_file": Path.cwd().parent / "case.foam",
    "run_cases_script": Path.cwd() / "run_cases1D.sh",
}
      
    
DX_VALUES = [10,20,40,80,160,320]# number of cells
DT_VALUES = [1,1,1,1,1,1]
#DT_VALUES = [1e-3, 1e-3,5e-4,5e-4,2e-4,1e-4]  # ms
TISSUE_TYPES = ["1D", "2D"]  #["1D", "2D", "3D"] 
FOLDER_NAME = ["error_manufactured", "error_manufactured_2D"] 

def organize_dat_files(folder_name):
    """
    Move all .dat files from the parent directory into a subfolder
    inside the parent directory.

    Parameters:
        folder_name (str): Name of the folder inside the parent directory.
    """
    # Absolute path of parent directory
    parent_dir = os.path.abspath("..")
    dest_dir = os.path.join(parent_dir, folder_name)
    print(f"Looking for .dat files in parent directory: {parent_dir}")

    # List all .dat files in parent directory
    dat_files = [f for f in os.listdir(parent_dir) if f.endswith(".dat")]

    if not dat_files:
        print("No .dat files found in the parent directory.")
        return

    # Create destination folder inside parent directory
    os.makedirs(dest_dir, exist_ok=True)

    # Move each .dat file to the folder
    for f in dat_files:
        src = os.path.join(parent_dir, f)
        dst = os.path.join(dest_dir, f)
        shutil.move(src, dst)
        print(f"Moved {f} -> {dest_dir}/")

    print(f"\nAll {len(dat_files)} .dat files moved to '{dest_dir}'.")





def read_error_dat_files(folder_name):
    """
    Reads all .dat files in parent/error_manufactured and extracts Linf errors for Vm, u1, u2.
    Returns a Pandas DataFrame with N and Solver type.
    """
    parent_dir = os.path.abspath("..")
    folder = os.path.join(parent_dir, folder_name)
    
    if not os.path.exists(folder):
        print("Folder does not exist:", folder)
        return None

    files = [f for f in os.listdir(folder) if f.endswith(".dat")]
    if not files:
        print("No .dat files found in folder:", folder)
        return None

    data = []

    for f in files:
        # Extract N and solver type from filename
        match = re.match(r"(\d+)_cells_(explicit|implicit)\.dat", f)
        if match:
            N = (int(match.group(1)))
            solver = match.group(2).capitalize()
        else:
            continue

        filepath = os.path.join(folder, f)
        with open(filepath, "r") as file:
            content = file.read()

            # Extract Linf errors using regex
            linf_matches = re.findall(
                r"Vm\s+\S+\s+\S+\s+(\S+).*?u1\s+\S+\s+\S+\s+(\S+).*?u2\s+\S+\s+\S+\s+(\S+)", 
                content.replace("\n", " "), re.DOTALL
            )
            if linf_matches:
                Linf_V, Linf_u1, Linf_u2 = map(float, linf_matches[0])
            else:
                Linf_V = Linf_u1 = Linf_u2 = np.nan

        data.append({
            "N": N,
            "Solver": solver,
            "Linf_V": Linf_V,
            "Linf_u1": Linf_u1,
            "Linf_u2": Linf_u2
        })

    df = pd.DataFrame(data)
    df.sort_values(["Solver", "N"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def compute_convergence_rates(df):
    """
    Compute convergence rates for Linf errors of Vm, u1, u2.
    Returns a DataFrame with rates between successive N values for each solver type.
    """
    rates_all = []

    for solver_type in df["Solver"].unique():
        df_solver = df[df["Solver"] == solver_type].sort_values("N").reset_index(drop=True)

        rates = {
            "Solver": [],
            "N_lower": [],
            "N_higher": [],
            "rate_Vm": [],
            "rate_u1": [],
            "rate_u2": []
        }

        for i in range(len(df_solver)-1):
            N1, N2 = df_solver.loc[i, "N"], df_solver.loc[i+1, "N"]
            h1, h2 = 1.0/N1, 1.0/N2

            errV1, errV2 = df_solver.loc[i, "Linf_V"], df_solver.loc[i+1, "Linf_V"]
            errU11, errU12 = df_solver.loc[i, "Linf_u1"], df_solver.loc[i+1, "Linf_u1"]
            errU21, errU22 = df_solver.loc[i, "Linf_u2"], df_solver.loc[i+1, "Linf_u2"]

            rate_V = np.log(errV1 / errV2) / np.log(h1 / h2) if errV1 > 0 and errV2 > 0 else np.nan
            rate_u1 = np.log(errU11 / errU12) / np.log(h1 / h2) if errU11 > 0 and errU12 > 0 else np.nan
            rate_u2 = np.log(errU21 / errU22) / np.log(h1 / h2) if errU21 > 0 and errU22 > 0 else np.nan

            rates["Solver"].append(solver_type)
            rates["N_lower"].append(N1)
            rates["N_higher"].append(N2)
            rates["rate_Vm"].append(rate_V)
            rates["rate_u1"].append(rate_u1)
            rates["rate_u2"].append(rate_u2)

        rates_all.append(pd.DataFrame(rates))

    rate_df = pd.concat(rates_all, ignore_index=True)
    return rate_df


def plot_errors(df, solver_type=None):
    """
    Plot Linf errors for Vm, u1, u2 vs N. Optionally filter by solver type.
    """
    if solver_type:
        df = df[df["Solver"] == solver_type]

    plt.figure(figsize=(8,5))
    plt.loglog(df["N"], df["Linf_V"], marker='o', label="Vm")
    plt.loglog(df["N"], df["Linf_u1"], marker='s', label="u1")
    plt.loglog(df["N"], df["Linf_u2"], marker='^', label="u2")
    plt.xlabel("Number of cells (N)")
    plt.ylabel("Linf Error")
    title = "Manufactured-solution Linf errors"
    if solver_type:
        title += f" ({solver_type})"
    plt.title(title)
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_errors_implicit_explicit(df, TISSUE_TYPE):
    """
    Plot Linf errors for Vm, u1, u2 vs N for both Explicit and Implicit solvers on the same plot.
    Matches column names from read_error_dat_files().
    """
    solvers = df["Solver"].unique()
    markers = {"Explicit": "o", "Implicit": "s"}
    linestyles = {"Explicit": "-", "Implicit": "--"}
    colors = {"Linf_V": "tab:blue", "Linf_u1": "tab:orange", "Linf_u2": "tab:green"}
    labels = {"Linf_V": "Vm", "Linf_u1": "u1", "Linf_u2": "u2"}

    plt.figure(figsize=(8, 6))

    for solver in solvers:
        df_solver = df[df["Solver"] == solver]
        for col, color in colors.items():
            plt.loglog(
                df_solver["N"],
                df_solver[col],
                marker=markers.get(solver, "o"),
                linestyle=linestyles.get(solver, "-"),
                color=color,
                label=f"{labels[col]} ({solver})"
            )

    plt.xlabel("Number of cells (N)")
    plt.ylabel("Linf Error")
    plt.title(f"Manufactured-solution Linf Errors {TISSUE_TYPE} (Explicit vs Implicit)")
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_Vm_across_tissues(dfs, tissue_types):
    """
    Plot Linf_V (Vm error) vs N for multiple tissue types (e.g. 1D, 2D)
    and both solvers (Explicit, Implicit) on the same figure.

    Parameters:
        dfs (list of pd.DataFrame): one DataFrame per tissue type.
        tissue_types (list of str): corresponding tissue type names.
    """
    plt.figure(figsize=(8, 6))

    markers = {"Explicit": "o", "Implicit": "s"}
    linestyles = {"Explicit": "-", "Implicit": "--"}
    colors = {"1D": "tab:blue", "2D": "tab:orange", "3D": "tab:green"}

    for df, tissue in zip(dfs, tissue_types):
        if df is None or df.empty:
            continue

        for solver in df["Solver"].unique():
            df_solver = df[df["Solver"] == solver]
            plt.loglog(
                df_solver["N"],
                df_solver["Linf_V"],
                marker=markers.get(solver, "o"),
                linestyle=linestyles.get(solver, "-"),
                color=colors.get(tissue, "black"),
                label=f"{tissue} ({solver})"
            )

    plt.xlabel("Number of cells (N)")
    plt.ylabel("Linf Error (Vm)")
    plt.title("Linf Error of Vm across tissue types (Explicit vs Implicit)")
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()





def main():
    selected_combinations = select_simulation_parameters(DX_VALUES, DT_VALUES, piecewise=True)
    if selected_combinations is None or selected_combinations.empty:
        print("No simulation parameters selected. Exiting.")
        return

    # run_simulation_cases(selected_combinations, CONFIG)

    dfs = []  # store DataFrames for 1D, 2D, etc.

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
        plot_errors(df, solver_type="Explicit")
        plot_errors(df, solver_type="Implicit")
        plot_errors_implicit_explicit(df, tissue_type)

        # Convergence rates
        rates = compute_convergence_rates(df)
        print(f"\nConvergence rates for {tissue_type}:")
        print(rates)

    # --- Combined Vm plot across tissue types ---
    print("\n=== Plotting Vm comparison across tissues ===")
    plot_Vm_across_tissues(dfs, TISSUE_TYPES)

        


if __name__ == "__main__":
    main()
