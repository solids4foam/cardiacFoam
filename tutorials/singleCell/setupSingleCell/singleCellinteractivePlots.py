"""
singleCell_Electrophysiology: singleCellinteractivePlots.py

This script plots simulation outputs from multiple model-cell combinations using Plotly.
The input is read from a folder provided to the post-processing entrypoint.

Each input entry in `file_dfs_vars` is a tuple:
    (file_name, dataframe, selected_variables)
where:
    - file_name: The name of the simulation output file (e.g. 'BuenoOrovio_mcells_output.txt')
    - dataframe: A pandas DataFrame containing columns for 'time' and all variable values
    - selected_variables: A dict with keys "States", "Algebraic", "Rates"
                          and values as lists of selected variable names.

The script:
1. Extracts model and cell type from file names:
       - Models are shortened using MODEL_MAP (e.g. 'BuenoOrovio' → 'BO')
       - Cells are shortened using CELL_TYPE_MAP (e.g. 'mcells' → 'Myo')
2. Allows user to optionally rename variable legends before plotting.
3. Creates one Plotly trace per variable per file, with legend name:
       "<Model>-<Cell>: <Variable>"
4. Provides interactive buttons to toggle groups of traces:
       - By category (States, Algebraic, Rates)
       - By model (BO, CO, TNNP, Gaur, etc.)
       - By cell type (Myo, Epi, Endo)
       - A "Show All" button to display everything

Main functions:
    - rename_legends()             → Optional renaming of variable labels
    - collect_variables_for_legend → Collects all variables across files
    - add_traces()                  → Creates and adds traces to the Plotly figure
    - build_buttons()               → Creates button controls for categories/models/cells
    - plot_multiple_files()         → Entry point, builds the full interactive plot

Output:
    - An interactive Plotly figure in the browser with toggleable groups of variables.
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any
from prompt_toolkit import prompt
from prompt_toolkit.shortcuts import checkboxlist_dialog
import plotly.graph_objs as go

TUTORIALS_ROOT = Path(__file__).resolve().parents[2]
if str(TUTORIALS_ROOT) not in sys.path:
    sys.path.insert(0, str(TUTORIALS_ROOT))

from openfoam_driver.postprocessing.plotting_common import (
    build_visibility_mask,
    ordered_unique,
    parse_model_and_cell,
)
from openfoam_driver.postprocessing.style import apply_plotly_layout
from openfoam_driver.postprocessing.style import write_plotly_html


STATE_COUNTS = {
    "Gaur": 29,
    "TNNP": 17,
    "Courtemanche": 21,
    "BuenoOrovio": 4
}
# Model abbreviation mapping
MODEL_MAP = {
    "BuenoOrovio": "BO",
    "Courtemanche": "CO",
    "TNNP": "TNNP",
    "Gaur": "Gaur"
}

# Cell type mapping
CELL_TYPE_MAP = {
    "mCells": "Myo",
    "epicardialCells": "Epi",
    "endocardialCells": "Endo",
    "myocyte": " "
}


# -------------------------------
def list_txt_files(folder):
    """Return all .txt files in the specified folder."""
    if not os.path.isdir(folder):
        print(f"The directory '{folder}' does not exist.")
        return []

    return sorted(f for f in os.listdir(folder) if f.endswith('.txt'))




def select_files(
    files,
    *,
    interactive: bool = True,
    preset_files: list[str] | None = None,
):
    if not interactive:
        if preset_files is None:
            return list(files)
        return [name for name in files if name in set(preset_files)]

    return checkboxlist_dialog(
        title="Select Simulation Output Files",
        text="Choose one or more files to load:",
        values=[(f, f) for f in files]
    ).run() or []

def detect_model_and_states(filename, *, interactive: bool = True):
    model_name = filename.split('_')[0].strip()
    n_states = STATE_COUNTS.get(model_name)
    if n_states is None:
        if not interactive:
            raise ValueError(
                f"Model '{model_name}' is missing from STATE_COUNTS; "
                "add it to STATE_COUNTS for non-interactive post-processing."
            )
        n_states = int(prompt("Model not recognized. Enter number of states: "))
    return model_name, n_states

def load_simulation_data(filename, *, base_folder: str | None = None):
    if not os.path.isabs(filename) and base_folder is not None:
        filename = os.path.join(base_folder, filename)

    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")


    with open(filename, 'r') as f:
        header = f.readline().strip().split()
    
    data = np.loadtxt(filename, skiprows=1)
    df = pd.DataFrame(data, columns=header)
    return df, header

def categorize_columns(header, n_states):

    
    states = header[1:n_states+1]
    rates = header[-n_states:]
    algebraic = header[n_states+1:-n_states]

    return states, rates, algebraic


def select_variables(
    states,
    algebraic,
    rates,
    file_name,
    *,
    interactive: bool = True,
    selected_categories: list[str] | None = None,
):
    categories = selected_categories or ["States", "Algebraic", "Rates"]
    category_set = {name for name in categories}
    if not interactive:
        return {
            "States": list(states) if "States" in category_set else [],
            "Algebraic": list(algebraic) if "Algebraic" in category_set else [],
            "Rates": list(rates) if "Rates" in category_set else [],
        }

    print(f"Selecting variables for file: {file_name}\n")
    selected_states = checkboxlist_dialog(
        title=f"States ({file_name})",
        text="Select State Variables:",
        values=[(s,s) for s in states]
    ).run() or []

    selected_algebraic = checkboxlist_dialog(
        title=f"Algebraic ({file_name})",
        text="Select Algebraic Variables:",
        values=[(a,a) for a in algebraic]
    ).run() or []

    selected_rates = checkboxlist_dialog(
        title=f"Rates ({file_name})",
        text="Select Rate Variables:",
        values=[(r,r) for r in rates]
    ).run() or []

    return {
        "States": selected_states,
        "Algebraic": selected_algebraic,
        "Rates": selected_rates
    }

def rename_legends(variable_list, *, interactive: bool = True):
    """Ask user to rename variables; Enter keeps default."""
    if not interactive:
        return list(variable_list)
    print("Rename variables with same physical meaning. Enter to skip")
    legend_names = []
    for var in variable_list:
        new_name = input(f"{var} (press Enter to keep): ")
        legend_names.append(new_name.strip() if new_name.strip() else var)
    return legend_names

def collect_variables_for_legend(file_dfs_vars, *, interactive: bool = True):
    all_vars = {"States": [], "Algebraic": [], "Rates": []}
    for _, _, selected_variables in file_dfs_vars:
        for category, var_list in selected_variables.items():
            all_vars[category].extend(var_list)
    # Remove duplicates
    for cat in all_vars:
        all_vars[cat] = ordered_unique(all_vars[cat])
    # Ask for legend renaming
    legend_map = {}
    for category, var_list in all_vars.items():
        if var_list:
            legend_map[category] = rename_legends(var_list, interactive=interactive)
        else:
            legend_map[category] = []
    return all_vars, legend_map

def add_traces(fig, file_dfs_vars, all_vars, legend_map):
    categories = ("States", "Algebraic", "Rates")
    category_trace_indices = {category: [] for category in categories}
    model_trace_indices: dict[str, list[int]] = {}
    cell_trace_indices: dict[str, list[int]] = {}

    for file_name, df, selected_variables in file_dfs_vars:
        time = df['time']
        model, cell_type = parse_model_and_cell(
            file_name,
            model_map=MODEL_MAP,
            cell_map=CELL_TYPE_MAP,
        )

        for category, var_list in selected_variables.items():
            if not var_list:
                continue
            for var in var_list:
                idx = all_vars[category].index(var)
                lname = legend_map[category][idx]
                trace = go.Scatter(
                    x=time, y=df[var], mode='lines',
                    name=f"{model}-{cell_type}: {lname}",
                    visible=True
                )
                fig.add_trace(trace)
                trace_index = len(fig.data) - 1
                category_trace_indices.setdefault(category, []).append(trace_index)
                model_trace_indices.setdefault(model, []).append(trace_index)
                cell_trace_indices.setdefault(cell_type, []).append(trace_index)

    return category_trace_indices, model_trace_indices, cell_trace_indices


def build_buttons(fig, category_trace_indices, model_trace_indices, cell_trace_indices):
    categories = ["States", "Algebraic", "Rates"]
    total_traces = len(fig.data)
    all_vis = build_visibility_mask(range(total_traces), total_traces)

    buttons = [dict(label="Show All", method="update", args=[{"visible": all_vis}])]
    buttons += [
        dict(
            label=f"Category: {category}",
            method="update",
            args=[{"visible": build_visibility_mask(category_trace_indices.get(category, []), total_traces)}],
        )
        for category in categories
    ]
    buttons += [
        dict(
            label=f"Model: {model}",
            method="update",
            args=[{"visible": build_visibility_mask(model_trace_indices[model], total_traces)}],
        )
        for model in ordered_unique(model_trace_indices.keys())
    ]
    buttons += [
        dict(
            label=f"Cell: {cell}",
            method="update",
            args=[{"visible": build_visibility_mask(cell_trace_indices[cell], total_traces)}],
        )
        for cell in ordered_unique(cell_trace_indices.keys())
    ]
    return buttons


def plot_multiple_files(
    file_dfs_vars,
    output_html: str | Path | None = None,
    show: bool = True,
    *,
    interactive: bool = True,
):
    fig = go.Figure()
    all_vars, legend_map = collect_variables_for_legend(
        file_dfs_vars,
        interactive=interactive,
    )
    category_traces, model_traces, cell_traces = add_traces(fig, file_dfs_vars, all_vars, legend_map)
    buttons = build_buttons(fig, category_traces, model_traces, cell_traces)

    apply_plotly_layout(
        fig,
        title="Simulation Variables",
        xaxis_title="Time [ms]",
        yaxis_title="Voltage [mV], Currents [pA/pF], Concentrations [mM], Rates [/ms]",
        updatemenus=[dict(type="buttons", x=1.05, y=0.8, buttons=buttons)],
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    if output_html is not None:
        write_plotly_html(fig, output_html)
    if show:
        fig.show()
    return fig


def post_processing_single_cell(
    output_folder,
    *,
    output_html: str | Path | None = None,
    show: bool = True,
    interactive: bool = True,
    selected_files: list[str] | None = None,
    selected_categories: list[str] | None = None,
    rename_legend_labels: bool | None = None,
):
    """Post-process simulation results from the given output folder."""
    rename_legend_labels_resolved = (
        interactive if rename_legend_labels is None else rename_legend_labels
    )
    files = list_txt_files(output_folder)
    if not files:
        print("No .txt files found. Exiting.")
        return

    selected_files_resolved = select_files(
        files,
        interactive=interactive,
        preset_files=selected_files,
    )
    if not selected_files_resolved:
        print("No files selected. Exiting.")
        return

    file_dfs_vars = []
    for file in selected_files_resolved:
        model_name, n_states = detect_model_and_states(file, interactive=interactive)
        print(f"Detected model: {model_name}, Number of states: {n_states}")

        # Construct full file path
        df, header = load_simulation_data(file, base_folder=output_folder)
        states, rates, algebraic = categorize_columns(header, n_states)

        selected_variables = select_variables(
            states,
            algebraic,
            rates,
            file,
            interactive=interactive,
            selected_categories=selected_categories,
        )
        file_dfs_vars.append((file, df, selected_variables))

    plot_multiple_files(
        file_dfs_vars,
        output_html=output_html,
        show=show,
        interactive=rename_legend_labels_resolved,
    )
    artifacts: list[dict[str, Any]] = []
    if output_html is not None:
        artifacts.append(
            {
                "path": str(output_html),
                "label": "singleCell interactive variables",
                "kind": "plot",
                "format": "html",
            }
        )
    return artifacts


def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    files: list[str] | None = None,
    categories: list[str] | None = None,
    rename_legends: bool = False,
    **_: object,
) -> None:
    del setup_root
    output_html = Path(output_dir) / "singleCell_interactive_variables.html"
    return post_processing_single_cell(
        output_dir,
        output_html=output_html,
        show=False,
        interactive=False,
        selected_files=files,
        selected_categories=categories,
        rename_legend_labels=rename_legends,
    )


if __name__ == "__main__":
    OUTPUT_FOLDER = os.path.join(os.path.dirname(os.getcwd()), "singleCellOutputs")
    post_processing_single_cell(OUTPUT_FOLDER)
