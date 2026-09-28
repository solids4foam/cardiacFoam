import json
import math
import re
from copy import deepcopy
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.graph_objects as go

from omnidriver.postprocessing.plotting_common import build_visibility_mask, lighten_hex_color
from omnidriver.postprocessing.style import apply_plotly_layout, write_plotly_html

_PROBE_LINE = re.compile(r"^# Probe (\d+) \(([^)]+)\)")


def rename_cardiacfoam_trace(name: str) -> str:
    """Strip the ΔT suffix and label the trace as cardiacFoam output."""
    if ", ΔT=" in name:
        return name.split(", ΔT=")[0] + " cardiacFoam"
    return name


# ----------------------------------------------------------
# --- Helper Functions -------------------------------------
# ----------------------------------------------------------

# Base colors per DX group
DX_COLORS = {
    0.1: "#eb1616",   # red
    0.2: "#1f77b4",   # blue
    0.5: "#2ca02c",   # green
}


def _read_latest_probe_row(function_object_dir: Path, field: str) -> tuple[list[tuple[float, float, float]], list[float]] | None:
    """Return (probe_xyz, values) from a `probes` functionObject's last row.

    See table_summary.py's identically-named helper for why every instance
    directory must be checked rather than trusting its name.
    """
    if not function_object_dir.is_dir():
        return None
    for instance_dir in sorted(p for p in function_object_dir.iterdir() if p.is_dir()):
        sample_path = instance_dir / field
        if not sample_path.is_file():
            continue
        probes: list[tuple[float, float, float]] = []
        last_row: list[float] | None = None
        for line in sample_path.read_text().splitlines():
            match = _PROBE_LINE.match(line)
            if match:
                _, coords = match.groups()
                x, y, z = (float(v) for v in coords.split())
                probes.append((x, y, z))
                continue
            if not line.strip() or line.startswith("#"):
                continue
            last_row = [float(token) for token in line.split()]
        if last_row is not None:
            return probes, last_row[1:]
    return None


def load_swept_cases(output_dir) -> list[tuple[float, float, pd.DataFrame]]:
    """Build [(DX_mm, DT_ms, DataFrame(arc_length_m, activationTime_s)), ...].

    Replaces the old convention of many `*line*.csv` files (one per dx/dt,
    dropped into one shared folder with dx/dt encoded in the filename) with
    this sweep's own `<output_dir>/cases/case_NNNN/postProcessing/
    Niedererlines/` -- one directory per case, dx/dt read from
    sweep_manifest.json instead of decoded from a name.
    """
    output_dir = Path(output_dir)
    manifest_path = output_dir / "sweep_manifest.json"
    if not manifest_path.is_file():
        print(f"No sweep_manifest.json found in: {output_dir}")
        return []
    manifest = json.loads(manifest_path.read_text())

    cases: list[tuple[float, float, pd.DataFrame]] = []
    for case in manifest.get("cases", []):
        case_dir = output_dir / "cases" / case["case_id"]
        axes = case.get("resolved_axis_values", {})
        dx_m = axes.get("dx", axes.get("tetDx"))
        dt_s = axes.get("system/controlDict:deltaT")
        if dx_m is None or dt_s is None:
            continue
        sample = _read_latest_probe_row(case_dir / "postProcessing" / "Niedererlines", "activationTime")
        if sample is None:
            continue
        probes, values = sample
        origin = probes[0] if probes else (0.0, 0.0, 0.0)
        rows = [
            {"activationTime": value, "arc_length": math.dist(origin, xyz)}
            for xyz, value in zip(probes, values)
            if value != -1.0  # this repository's "never activated" sentinel; never plotted as a distance
        ]
        if not rows:
            continue
        cases.append((float(dx_m) * 1000.0, float(dt_s) * 1000.0, pd.DataFrame(rows)))

    print(f"Found {len(cases)} swept cases with diagonal Activation time:")
    for dx, dt, _df in cases:
        print(f"  - DX={dx:.3f} mm, DT={dt:.4f} ms")
    return cases


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
# --- dynamic DT shading per DX ----------
# ----------------------------------------------------------

def add_case_traces(fig, sorted_items, dt_target=None):
    """Add CardiacFoam case traces, return indices of traces matching dt_target."""
    matching_indices = []

    # --- STEP 1: discover DT values for each DX dynamically ---
    dx_dt_map: dict[float, set[float]] = {}
    for dx, dt, _df in sorted_items:
        dx_dt_map.setdefault(dx, set()).add(dt)

    # Sort DT list for each DX
    for dx in dx_dt_map:
        dx_dt_map[dx] = sorted(dx_dt_map[dx])

    # --- STEP 2: add traces with shading ---
    for dx, dt, df in sorted_items:
        full_label = f"ΔX={dx:.1f} mm, ΔT={dt:.3f} ms"

        # activationTime, arc_length columns; both still in SI units (s, m)
        y_col, x_col = "activationTime", "arc_length"

        # Base color for this DX
        base_color = DX_COLORS.get(dx, "#808080")

        # Determine shade for this DT
        dt_list = dx_dt_map[dx]
        dt_index = dt_list.index(dt)
        count = len(dt_list)
        shade_amount = dt_index / max(count - 1, 1) * 0.4  # 0 → darkest, 1 → lightest

        color = lighten_hex_color(base_color, shade_amount)
        dash_style = 'dashdot'

        fig.add_trace(go.Scatter(
            x=df[x_col] * 1000,
            y=df[y_col] * 1000,
            mode='lines+markers',
            name=full_label,
            line=dict(color=color, dash=dash_style),
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
    visibility_toggle = build_visibility_mask(toggle_set, len(fig.data))
    visibility_initial = build_visibility_mask(
        set(range(len(fig.data))).difference(excel_indices),
        len(fig.data),
    )

    cleaned_names = [rename_cardiacfoam_trace(trace.name) for trace in fig.data]

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

def plot_line_csvs(folder='.', excel_path=None, show: bool = True):
    """Main plotting function (logic unchanged, just organized)."""

    output_folder = Path(folder)
    sorted_items = sorted(load_swept_cases(output_folder), key=lambda item: (item[0], item[1]))

    fig = go.Figure()
    has_excel = excel_path is not None

    if has_excel:
        df_excel, num_cols = load_excel_data(excel_path)
        excel_indices = add_excel_traces(fig, df_excel, num_cols)

        dt_target = 0.005
        filtered = [dx for dx, dt, _df in sorted_items if abs(dt - dt_target) < 1e-6]

        if filtered:
            print(f"Found {len(filtered)} cases with ΔT={dt_target} ms:")
        else:
            print(f"No cases found for ΔT={dt_target} ms.")
    else:
        excel_indices = []
        dt_target = None

    csv_indices = add_case_traces(fig, sorted_items, dt_target=dt_target)

    if has_excel:
        add_toggle_button(fig, excel_indices, csv_indices, dt_target)

    apply_plotly_layout(
        fig,
        title="Activation time in diagonal line: Niederer N-benchmark",
        xaxis_title="Distance along diagonal line (mm)",
        yaxis_title="Activation Time (ms)",
        legend_title="Resolution",
        template="plotly_white",
        showlegend=True
    )

    # Move legend outside to avoid overlap
    fig.update_layout(
        legend=dict(
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02
        )
    )

    # --- SAVE INITIAL PLOT (CSV only) ---
    fig_initial = deepcopy(fig)

    for i in range(len(fig_initial.data)):
        if i in excel_indices:   # hide all Niederer
            fig_initial.data[i].visible = False
        else:
            fig_initial.data[i].visible = True

    write_plotly_html(fig_initial, output_folder / "cardiacFoam_allSimulations.html")
    with open(output_folder / "cardiacFoam_allSimulations.json", "w") as f:
        f.write(fig_initial.to_json())

    # --- SAVE COMPARISON PLOT (Niederer + matching CSV only) ---
    fig_compare = deepcopy(fig)

    toggle_set = set(excel_indices + csv_indices)

    for i in range(len(fig_compare.data)):
        fig_compare.data[i].visible = (i in toggle_set)

    # Remove DT and add "cardiacFoam" in label for CSV traces
    for tr in fig_compare.data:
        tr.name = rename_cardiacfoam_trace(tr.name)

    write_plotly_html(fig_compare, output_folder / "Niederer_vs_cardiacFoam.html")
    with open(output_folder / "Niederer_vs_cardiacFoam.json", "w") as f:
        f.write(fig_compare.to_json())

    if show:
        fig.show()
    return fig


def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    excel_path: str | None = None,
    **_: object,
):
    """Plot activation time along the benchmark's diagonal probe line.

    Reads every case in this sweep's own `postProcessing/Niedererlines/`
    output, shades traces by dx/dt, and -- when excel_path points at the
    digitized Niederer et al. 2012 reference curves -- overlays them for
    comparison. Writes cardiacFoam_allSimulations.html (every case) and
    Niederer_vs_cardiacFoam.html (cases matching the reference's dt only).
    """
    del setup_root
    plot_line_csvs(folder=output_dir, excel_path=excel_path, show=False)
    artifacts: list[dict[str, Any]] = [
        {
            "path": str(Path(output_dir) / "cardiacFoam_allSimulations.html"),
            "label": "Niederer line: all simulations",
            "kind": "plot",
            "format": "html",
        },
        {
            "path": str(Path(output_dir) / "Niederer_vs_cardiacFoam.html"),
            "label": "Niederer line: benchmark comparison",
            "kind": "plot",
            "format": "html",
        },
    ]
    return artifacts


# ----------------------------------------------------------
# --- Entry Point ------------------------------------------
# ----------------------------------------------------------

if __name__ == "__main__":
    plot_line_csvs(
        folder='.',
        excel_path='Niederer_graphs_webplotdigitilizer_points_slab/WebPlotDigitilizerdata.xlsx'
    )
