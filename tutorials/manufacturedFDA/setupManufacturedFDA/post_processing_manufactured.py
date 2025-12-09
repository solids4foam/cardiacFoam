"""
Post proscessing script for manufactured solution convergence tests. 
Functions for organizing .dat files, reading errors, computing convergence rates, and plotting results.

"""
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import pandas as pd


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
    Reads all .dat files in folder_name and extracts:
        - Dimension  (1D, 2D, 3D)
        - N          (# cells)
        - Solver     (explicit, implicit)
        - Linf errors for Vm, u1, u2

    Returns a Pandas DataFrame: one row per file.
    """

    folder = Path(folder_name)
    if not folder.exists():
        print("Folder does not exist:", folder)
        return None

    files = [f for f in folder.iterdir() if f.suffix == ".dat"]
    if not files:
        print("No .dat files found in folder:", folder)
        return None

    data = []

    for f in files:
        # Expected filename format:
        #   1D_320_cells_explicit.dat
        m = re.match(r"(\dD)_(\d+)_cells_(explicit|implicit)", f.name)
        if not m:
            print("Skipping unrecognized filename:", f.name)
            continue

        dimension = m.group(1)   # "1D"
        N = int(m.group(2))      # 320
        solver = m.group(3)      # "explicit" or "implicit"

        content = f.read_text()

        # Extract Linf errors
        linf_matches = re.findall(
            r"Vm\s+\S+\s+\S+\s+(\S+).*?"
            r"u1\s+\S+\s+\S+\s+(\S+).*?"
            r"u2\s+\S+\s+\S+\s+(\S+)",
            content.replace("\n", " "), re.DOTALL
        )

        if linf_matches:
            Linf_V, Linf_u1, Linf_u2 = map(float, linf_matches[0])
        else:
            Linf_V = Linf_u1 = Linf_u2 = np.nan

        data.append({
            "Dimension": dimension,
            "N": N,
            "Solver": solver,
            "Linf_V": Linf_V,
            "Linf_u1": Linf_u1,
            "Linf_u2": Linf_u2
        })

    df = pd.DataFrame(data)
    df.sort_values(["Dimension", "Solver", "N"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df

def compute_convergence_rates(df):
    """
    Compute convergence rates for Linf errors of Vm, u1, u2.

    - Groups by Dimension (if present) and Solver.
    - Sorts by N.
    - Skips pairs where N_lower == N_higher.
    """

    rates_all = []

    # Decide grouping keys
    group_keys = ["Solver"]
    if "Dimension" in df.columns:
        group_keys = ["Dimension", "Solver"]

    for group_vals, df_group in df.groupby(group_keys):
        df_group = df_group.sort_values("N").reset_index(drop=True)

        # Nice unpack of group labels
        if "Dimension" in df.columns:
            dimension, solver_type = group_vals
        else:
            (solver_type,) = group_vals
            dimension = None

        rates = {
            "Dimension": [],
            "Solver": [],
            "N_lower": [],
            "N_higher": [],
            "rate_Vm": [],
            "rate_u1": [],
            "rate_u2": [],
        }

        # Walk over consecutive pairs
        for i in range(len(df_group) - 1):
            N1 = df_group.loc[i, "N"]
            N2 = df_group.loc[i+1, "N"]

            # Skip pairs with identical N
            if N1 == N2:
                continue

            # In your setup N is "cells in one direction", so dx ~ 1/N
            h1, h2 = 1.0 / N1, 1.0 / N2

            errV1, errV2  = df_group.loc[i,   "Linf_V"],  df_group.loc[i+1, "Linf_V"]
            errU11, errU12 = df_group.loc[i,  "Linf_u1"], df_group.loc[i+1, "Linf_u1"]
            errU21, errU22 = df_group.loc[i,  "Linf_u2"], df_group.loc[i+1, "Linf_u2"]

            def safe_rate(e1, e2):
                if e1 > 0 and e2 > 0 and h1 != h2:
                    return np.log(e1 / e2) / np.log(h1 / h2)
                else:
                    return np.nan

            rate_V  = safe_rate(errV1,  errV2)
            rate_u1 = safe_rate(errU11, errU12)
            rate_u2 = safe_rate(errU21, errU22)

            rates["Dimension"].append(dimension)
            rates["Solver"].append(solver_type)
            rates["N_lower"].append(N1)
            rates["N_higher"].append(N2)
            rates["rate_Vm"].append(rate_V)
            rates["rate_u1"].append(rate_u1)
            rates["rate_u2"].append(rate_u2)

        if rates["N_lower"]:  # only append non-empty groups
            rates_all.append(pd.DataFrame(rates))

    if not rates_all:
        return pd.DataFrame(columns=["Dimension", "Solver", "N_lower", "N_higher",
                                     "rate_Vm", "rate_u1", "rate_u2"])

    return pd.concat(rates_all, ignore_index=True)


def plot_errors(df, solver_type=None):
    """
    Plot Linf errors for Vm, u1, u2 vs N.
    """

    plt.figure(figsize=(8,5))

    for dimension in df["Dimension"].unique():
        df_dim = df[df["Dimension"] == dimension]

        if solver_type:
            df_dim = df_dim[df_dim["Solver"] == solver_type]

        plt.loglog(df_dim["N"], df_dim["Linf_V"], marker='o', label=f"Vm ({dimension})")
        plt.loglog(df_dim["N"], df_dim["Linf_u1"], marker='s', label=f"u1 ({dimension})")
        plt.loglog(df_dim["N"], df_dim["Linf_u2"], marker='^', label=f"u2 ({dimension})")

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

def plot_errors_implicit_explicit(df):
    """
    Plot Linf errors for Vm, u1, u2 vs N for both Explicit and Implicit solvers.
    Grouped per dimension.
    """
    markers = {"explicit": "o", "implicit": "s"}
    linestyles = {"explicit": "-", "implicit": "--"}
    colors = {"Linf_V": "tab:blue", "Linf_u1": "tab:orange", "Linf_u2": "tab:green"}
    labels = {"Linf_V": "Vm", "Linf_u1": "u1", "Linf_u2": "u2"}

    for dimension in df["Dimension"].unique():
        df_dim = df[df["Dimension"] == dimension]

        plt.figure(figsize=(8, 6))
        for solver in df_dim["Solver"].unique():
            df_solver = df_dim[df_dim["Solver"] == solver]

            for col, color in colors.items():
                plt.loglog(
                    df_solver["N"],
                    df_solver[col],
                    marker=markers.get(solver, "o"),
                    linestyle=linestyles.get(solver, "-"),
                    color=color,
                    label=f"{labels[col]} ({solver}, {dimension})"
                )

        plt.xlabel("Number of cells (N)")
        plt.ylabel("Linf Error")
        plt.title(f"Linf Errors per Dimension ({dimension})")
        plt.grid(True, which="both", ls="--", alpha=0.6)
        plt.legend()
        plt.tight_layout()
        plt.show()


def plot_Vm_across_dimensions(df):
    """
    Plot Linf_V (Vm error) vs N across all dimensions and solvers.
    """

    plt.figure(figsize=(8, 6))

    markers = {"explicit": "o", "implicit": "s"}
    linestyles = {"explicit": "-", "implicit": "--"}
    colors = {"1D": "tab:blue", "2D": "tab:orange", "3D": "tab:green"}

    for dimension in df["Dimension"].unique():
        df_dim = df[df["Dimension"] == dimension]

        for solver in df_dim["Solver"].unique():
            df_solver = df_dim[df_dim["Solver"] == solver]

            plt.loglog(
                df_solver["N"],
                df_solver["Linf_V"],
                marker=markers.get(solver, "o"),
                linestyle=linestyles.get(solver, "--"),
                color=colors.get(dimension, "black"),
                label=f"{dimension} ({solver})"
            )

    plt.xlabel("Number of cells (N)")
    plt.ylabel("Linf Error (Vm)")
    plt.title("Linf Error of Vm across dimensions (Explicit vs Implicit)")
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()

