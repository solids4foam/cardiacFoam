"""
singleCell_Electrophysiology: singleCellinteractivePlots.py

This script plots simulation outputs from multiple model-cell combinations using Plotly.
The input is read from a folder in the case directory defined by the variable OUTPUT_FOLDER.

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
import numpy as np
import pandas as pd
from prompt_toolkit import prompt
from prompt_toolkit.shortcuts import checkboxlist_dialog
import plotly.graph_objs as go


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

    return [f for f in os.listdir(folder) if f.endswith('.txt')]




def select_files(files):
    return checkboxlist_dialog(
        title="Select Simulation Output Files",
        text="Choose one or more files to load:",
        values=[(f, f) for f in files]
    ).run() or []

def detect_model_and_states(filename):
    model_name = filename.split('_')[0].strip()
    n_states = STATE_COUNTS.get(model_name)
    if n_states is None:
        n_states = int(prompt("Model not recognized. Enter number of states: "))
    return model_name, n_states

def load_simulation_data(filename):
    if not os.path.isabs(filename):
        filename = os.path.join(OUTPUT_FOLDER, filename)

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


def select_variables(states, algebraic, rates, file_name):
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

def rename_legends(variable_list):
    """Ask user to rename variables; Enter keeps default."""
    print("Rename variables with same physical meaning. Enter to skip")
    legend_names = []
    for var in variable_list:
        new_name = input(f"{var} (press Enter to keep): ")
        legend_names.append(new_name.strip() if new_name.strip() else var)
    return legend_names

def collect_variables_for_legend(file_dfs_vars):
    all_vars = {"States": [], "Algebraic": [], "Rates": []}
    for _, _, selected_variables in file_dfs_vars:
        for category, var_list in selected_variables.items():
            all_vars[category].extend(var_list)
    # Remove duplicates
    for cat in all_vars:
        all_vars[cat] = list(dict.fromkeys(all_vars[cat]))
    # Ask for legend renaming
    legend_map = {}
    for category, var_list in all_vars.items():
        if var_list:
            legend_map[category] = rename_legends(var_list)
        else:
            legend_map[category] = []
    return all_vars, legend_map

def add_traces(fig, file_dfs_vars, all_vars, legend_map):
    category_traces = {}
    file_model_map = {}
    file_cell_map = {}

    for file_name, df, selected_variables in file_dfs_vars:
        time = df['time']
        file_base = os.path.splitext(file_name)[0]
        parts = file_base.split('_')
        raw_model = parts[0] if len(parts) > 0 else "UnknownModel"
        raw_cell = parts[1] if len(parts) > 1 else "UnknownCell"

        model = MODEL_MAP.get(raw_model, raw_model)
        cell_type = CELL_TYPE_MAP.get(raw_cell, raw_cell)

        file_model_map[file_name] = model
        file_cell_map[file_name] = cell_type

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
                category_traces.setdefault((category, file_name), []).append(trace)

    return category_traces, file_model_map, file_cell_map

def build_buttons(category_traces, file_dfs_vars, file_model_map, file_cell_map):
    all_traces = [t for traces in category_traces.values() for t in traces]
    categories = ["States", "Algebraic", "Rates"]

    # Category buttons
    category_vis = []
    for cat in categories:
        if any((cat, f) in category_traces for f, _, _ in file_dfs_vars):
            vis = [any(t in category_traces.get((cat, f), []) for f, _, _ in file_dfs_vars)
                   for t in all_traces]
        else:
            vis = [False]*len(all_traces)
        category_vis.append(vis)

    # Model buttons (group all files of the same model)
    unique_models = list(dict.fromkeys(file_model_map.values()))
    model_vis = []
    for model in unique_models:
        vis = [any(t in category_traces.get((cat, f), [])
                   for cat in categories
                   for f in file_model_map if file_model_map[f] == model)
               for t in all_traces]
        model_vis.append(vis)

    # Cell buttons
    unique_cells = list(dict.fromkeys(file_cell_map.values()))
    cell_vis = []
    for cell in unique_cells:
        vis = [any(t in category_traces.get((cat, f), [])
                   for cat in categories
                   for f in file_cell_map if file_cell_map[f] == cell)
               for t in all_traces]
        cell_vis.append(vis)

    # Show All
    all_vis = [True]*len(all_traces)

    # Combine buttons
    buttons = [dict(label="Show All", method="update", args=[{"visible": all_vis}])]
    buttons += [dict(label=f"Category: {cat}", method="update", args=[{"visible": category_vis[i]}])
                for i, cat in enumerate(categories)]
    buttons += [dict(label=f"Model: {model}", method="update", args=[{"visible": model_vis[i]}])
                for i, model in enumerate(unique_models)]
    buttons += [dict(label=f"Cell: {cell}", method="update", args=[{"visible": cell_vis[i]}])
                for i, cell in enumerate(unique_cells)]

    return buttons

# -------------------------------
# Main plotting function
# -------------------------------
def plot_multiple_files(file_dfs_vars):
    fig = go.Figure()

    # Prepare legend names
    all_vars, legend_map = collect_variables_for_legend(file_dfs_vars)

    # Add traces
    category_traces, file_model_map, file_cell_map = add_traces(fig, file_dfs_vars, all_vars, legend_map)

    # Build buttons
    buttons = build_buttons(category_traces, file_dfs_vars, file_model_map, file_cell_map)

    # Update layout
    fig.update_layout(
        title="Simulation Variables",
        xaxis_title="Time [ms]",
        yaxis_title="Voltage [mV], Currents [pA/pF], Concentrations [mM], Rates [/ms]",
        updatemenus=[dict(type="buttons", x=1.05, y=0.8, buttons=buttons)],
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )

    fig.show()



def plot_multiple_files(file_dfs_vars):
    fig = go.Figure()

    all_vars, legend_map = collect_variables_for_legend(file_dfs_vars)
    category_traces, file_cell_map, file_model_map = add_traces(fig, file_dfs_vars, all_vars, legend_map)
    buttons = build_buttons(category_traces, file_dfs_vars, file_cell_map, file_model_map)

    fig.update_layout(
        title="Simulation Variables",
        xaxis_title="Time [ms]",
        yaxis_title="Voltage [mV], Currents [pA/pF], Concentrations [mM], Rates [/ms]",
        updatemenus=[dict(type="buttons", x=1.05, y=0.8, buttons=buttons)],
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )

    fig.show()


def post_processing_single_cell(output_folder):
    """Post-process simulation results from the given output folder."""
    files = list_txt_files(output_folder)
    if not files:
        print("No .txt files found. Exiting.")
        return

    selected_files = select_files(files)
    if not selected_files:
        print("No files selected. Exiting.")
        return

    file_dfs_vars = []
    for file in selected_files:
        model_name, n_states = detect_model_and_states(file)
        print(f"Detected model: {model_name}, Number of states: {n_states}")

        # Construct full file path
        file_path = os.path.join(output_folder, file)

        df, header = load_simulation_data(file_path)
        states, rates, algebraic = categorize_columns(header, n_states, model_name)

        selected_variables = select_variables(states, algebraic, rates, file)
        file_dfs_vars.append((file, df, selected_variables))

    plot_multiple_files(file_dfs_vars)


if __name__ == "__main__":
    OUTPUT_FOLDER = os.path.join(os.path.dirname(os.getcwd()), "singleCellOutputs")
    post_processing_single_cell(OUTPUT_FOLDER)