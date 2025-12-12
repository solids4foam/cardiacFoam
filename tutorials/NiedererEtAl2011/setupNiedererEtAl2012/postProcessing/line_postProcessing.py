import glob
import os
import re
import pandas as pd
import plotly.graph_objects as go
import plotly.colors as pc
from copy import deepcopy


# ----------------------------------------------------------
# --- Helper Functions -------------------------------------
# ----------------------------------------------------------

def lighten_color(color, amount):
    """Lighten a color by amount (0=none, 1=white)."""
    r, g, b = pc.hex_to_rgb(color)
    R = int(r + (255 - r) * amount)
    G = int(g + (255 - g) * amount)
    B = int(b + (255 - b) * amount)
    return f"rgb({R},{G},{B})"


# Base colors per DX group
DX_COLORS = {
    0.1: "#1f77b4",   # blue
    0.2: "#2ca02c",   # green
    0.5: "#eb1616",   # red
}


def load_csv_files(folder):
    """Load line*.csv files into a dict {filename: dataframe}."""
    csv_files = glob.glob(os.path.join(folder, '*line*.csv'))
    line_data = {os.path.basename(fp): pd.read_csv(fp) for fp in csv_files}

    print(f"Found {len(line_data)} CSVs for diagonal Activation time:")
    for name in line_data:
        print(f"  - {name}")

    return line_data


def clean_csv_data(line_data):
    """Drop non-relevant columns from all CSV dataframes."""
    columns_to_drop = [
        'Point ID', 'Points_Magnitude', 'Vm', 'vtkValidPointMask',
        'ionicCurrent', 'Points_0', 'Points_1', 'Points_2'
    ]

    cleaned = {}
    for filename, df in line_data.items():
        to_drop = [c for c in columns_to_drop if c in df.columns]
        cleaned[filename] = df.drop(columns=to_drop)

    return cleaned


def extract_dx_dt(filename):
    """Extract ΔX and ΔT from filename like DX5 / DT005."""
    dx_match = re.search(r'DX(\d+)', filename)
    dt_match = re.search(r'DT(\d+)', filename)

    if dx_match and dt_match:
        dx = int(dx_match.group(1)) / 10.0
        dt_raw = dt_match.group(1)
        dt = int(dt_raw) / (10 ** (len(dt_raw) - 1))
        return dx, dt

    return float("inf"), float("inf")


def load_excel_data(excel_path):
    """Load Niederer Excel columns in pairs."""
    df_excel = pd.read_excel(excel_path)
    num_cols = df_excel.shape[1]

    if num_cols % 2 != 0:
        print("Warning: Excel file has odd number of columns. Ignoring last column.")
        num_cols -= 1

    return df_excel, num_cols


def add_excel_traces(fig, df_excel, num_cols):
    """Add Niederer traces and return their indices."""
    bg_indices = []
    colors = ['red', 'blue', 'green']

    labels = [
        "ΔX=0.1 mm Niederer",
        "ΔX=0.2 mm Niederer",
        "ΔX=0.5 mm Niederer",
    ]

    for i in range(0, num_cols, 2):
        x_col, y_col = df_excel.columns[i], df_excel.columns[i + 1]
        color = colors[(i // 2) % len(colors)]
        label = labels[(i // 2) % len(labels)]

        df_sorted = df_excel[[x_col, y_col]].dropna().sort_values(by=x_col)

        trace = go.Scatter(
            x=df_sorted[x_col], y=df_sorted[y_col],
            mode='lines',
            name=label,
            line=dict(color=color, dash='dash'),
            visible=False
        )

        fig.add_trace(trace)
        bg_indices.append(len(fig.data) - 1)

    return bg_indices


# ----------------------------------------------------------
# --- UPDATED FUNCTION: dynamic DT shading per DX ----------
# ----------------------------------------------------------

def add_csv_traces(fig, sorted_items, extract_dx_dt, dt_target=None):
    """Add CardiacFoam CSV traces, return indices of traces matching dt_target."""
    matching_indices = []

    # --- STEP 1: discover DT values for each DX dynamically ---
    dx_dt_map = {}
    for filename, _ in sorted_items:
        dx, dt = extract_dx_dt(filename)
        dx_dt_map.setdefault(dx, set()).add(dt)

    # Sort DT list for each DX
    for dx in dx_dt_map:
        dx_dt_map[dx] = sorted(dx_dt_map[dx])

    # --- STEP 2: add traces with shading ---
    for filename, df in sorted_items:
        dx, dt = extract_dx_dt(filename)
        full_label = f"ΔX={dx:.1f} mm, ΔT={dt:.3f} ms"

        cols = df.columns.tolist()
        if len(cols) < 2:
            print(f"Warning: Skipping {filename} due to missing columns.")
            continue

        y_col, x_col = cols[0], cols[1]

        # Base color for this DX
        base_color = DX_COLORS.get(dx, "#808080")

        # Determine shade for this DT
        dt_list = dx_dt_map[dx]
        dt_index = dt_list.index(dt)
        count = len(dt_list)
        shade_amount = dt_index / max(count - 1, 1)* 0.4  # 0 → darkest, 1 → lightest

        color = lighten_color(base_color, shade_amount)

        fig.add_trace(go.Scatter(
            x=df[x_col] * 1000,
            y=df[y_col] * 1000,
            mode='lines+markers',
            name=full_label,
            line=dict(color=color),
            marker=dict(color=color)
        ))

        if dt_target is not None and abs(dt - dt_target) < 1e-6:
            matching_indices.append(len(fig.data) - 1)

    return matching_indices


# ----------------------------------------------------------
# --- Toggling system unchanged ----------------------------
# ----------------------------------------------------------

def add_toggle_button(fig, excel_indices, csv_indices, dt_target):
    """Add the Niederer-vs-CSV toggle button with reversible behavior."""

    toggle_set = set(excel_indices + csv_indices)
    visibility_toggle = [(i in toggle_set) for i in range(len(fig.data))]

    visibility_initial = []
    for i in range(len(fig.data)):
        if i in excel_indices:
            visibility_initial.append(False)
        else:
            visibility_initial.append(True)

    cleaned_names = []
    for trace in fig.data:
        name = trace.name
        if ", ΔT=" in name:
            name = name.split(", ΔT=")[0] + " cardiacFoam"
        cleaned_names.append(name)

    original_names = [trace.name for trace in fig.data]

    fig.update_layout(
        updatemenus=[{
            "type": "buttons",
            "direction": "right",
            "x": 0.2,
            "y": 1.0,
            "xanchor": "left",
            "yanchor": "top",
            "showactive": False,
            "buttons": [
                {
                    "label": "All simulations",
                    "method": "update",
                    "args": [
                        {"visible": visibility_initial, "name": original_names}
                    ]
                },
                {
                    "label": f"Niederer VS cardiacFoam (ΔT={dt_target} ms)",
                    "method": "update",
                    "args": [
                        {"visible": visibility_toggle, "name": cleaned_names}
                    ]
                },
            ]
        }]
    )


# ----------------------------------------------------------
# --- Main Function ----------------------------------------
# ----------------------------------------------------------

def plot_line_csvs(folder='.', excel_path=None):
    """Main plotting function (logic unchanged, just organized)."""

    line_data = load_csv_files(folder)
    cleaned_data = clean_csv_data(line_data)

    sorted_items = sorted(cleaned_data.items(), key=lambda item: extract_dx_dt(item[0]))

    fig = go.Figure()
    has_excel = excel_path is not None

    if has_excel:
        df_excel, num_cols = load_excel_data(excel_path)
        excel_indices = add_excel_traces(fig, df_excel, num_cols)

        dt_target = 0.005
        filtered = [name for name, _ in sorted_items
                   if abs(extract_dx_dt(name)[1] - dt_target) < 1e-6]

        if filtered:
            print(f"Found {len(filtered)} CSV files with ΔT={dt_target} ms:")
            for f in filtered:
                print(f"  - {f}")
        else:
            print(f"No CSVs found for ΔT={dt_target} ms.")
    else:
        excel_indices = []
        dt_target = None

    csv_indices = add_csv_traces(fig, sorted_items, extract_dx_dt, dt_target)

    if has_excel:
        add_toggle_button(fig, excel_indices, csv_indices, dt_target)

    fig.update_layout(
        title="Activation time in diagonal line: Niederer N-benchmark",
        xaxis_title="Distance along diagonal line (mm)",
        yaxis_title="Activation Time (ms)",
        legend_title="Resolution",
        template="plotly_white",
        showlegend=True
    )
    

    # --- SAVE INITIAL PLOT (CSV only) ---
    fig_initial = deepcopy(fig)

    for i in range(len(fig_initial.data)):
        if i in excel_indices:   # hide all Niederer
            fig_initial.data[i].visible = False
        else:
            fig_initial.data[i].visible = True

    fig_initial.write_html("cardiacFoam_allSimulations.html")




    # --- SAVE COMPARISON PLOT (Niederer + matching CSV only) ---
    fig_compare = deepcopy(fig)

    toggle_set = set(excel_indices + csv_indices)

    for i in range(len(fig_compare.data)):
        fig_compare.data[i].visible = (i in toggle_set)

    # Remove DT and add "cardiacFoam" in label for CSV traces
    for tr in fig_compare.data:
        if ", ΔT=" in tr.name:
            tr.name = tr.name.split(", ΔT=")[0] + " cardiacFoam"

    fig_compare.write_html("Niederer_vs_cardiacFoam.html")

    fig.show()
    return fig


# ----------------------------------------------------------
# --- Entry Point ------------------------------------------
# ----------------------------------------------------------

if __name__ == "__main__":
    plot_line_csvs(
        folder='.',
        excel_path='Niederer_graphs_webplotdigitilizer_points_slab/WebPlotDigitilizerdata.xlsx'
    )


