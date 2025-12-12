from pathlib import Path
import re
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# -----------------------------------------------------------
# Extract DX and DT
# -----------------------------------------------------------
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


# -----------------------------------------------------------
# Load all point CSVs
# -----------------------------------------------------------
def load_all_point_files(folder):
    folder = Path(folder)
    files = sorted(folder.glob("*points_DT*_DX*.csv"))

    if not files:
        print(f"No point CSV files found in: {folder}")
        return []

    all_data = []

    for fpath in files:
        filename = fpath.name
        dx, dt = extract_dx_dt(filename)
        df = pd.read_csv(fpath)

        activation = df["activationTime"].values * 1000.0  # to ms

        all_data.append({
            "DX": dx,
            "DT": dt,
            "activation": activation,
            "df": df
        })

    print(f"Loaded {len(all_data)} point files from {folder}")
    for e in all_data:
        print(f"  - DX={e['DX']:.3f} mm, DT={e['DT']:.4f} ms")

    return all_data


# -----------------------------------------------------------
# Find earliest activated point
# -----------------------------------------------------------
def find_earliest_point(all_data):
    if not all_data:
        raise ValueError("Cannot find earliest point: no data.")

    num_points = len(all_data[0]["activation"])
    avg_times = []

    for p in range(num_points):
        vals = [e["activation"][p] for e in all_data]
        avg_times.append(np.mean(vals))

    return int(np.argmin(avg_times))


# -----------------------------------------------------------
# Load geometry (shared)
# -----------------------------------------------------------
def load_point_geometry_from_sample(folder):
    folder = Path(folder)
    sample = sorted(folder.glob("points_DT*_DX*.csv"))
    if not sample:
        raise FileNotFoundError("No CSV files for geometry.")

    df = pd.read_csv(sample[0])

    xs = df["Points:0"].values
    ys = df["Points:1"].values
    zs = df["Points:2"].values

    return xs, ys, zs


def plot_3d_points_and_grid(folder="."):
    all_data = load_all_point_files(folder)
    if not all_data:
        return

    earliest = find_earliest_point(all_data)
    num_points = len(all_data[0]["activation"])
    points_to_plot = [p for p in range(num_points) if p != earliest]

    palette = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
        "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"
    ]
    point_colors = {p: palette[p] for p in range(num_points)}

    DX_all = np.array([e["DX"] for e in all_data])
    DT_all = np.array([e["DT"] for e in all_data])

    # ------------------------------------------
    # USE make_subplots WITH 3D SUBPLOTS
    # ------------------------------------------
    fig = make_subplots(
        rows=2, cols=4,
        specs=[[{'type': 'scene'}] * 4,
               [{'type': 'scene'}] * 4],
        horizontal_spacing=0.04,   # spacing between columns
        vertical_spacing=0.08      # spacing between rows
    )

    scene_id = 1

    for idx, p in enumerate(points_to_plot):
        row = idx // 4 + 1
        col = idx % 4 + 1

        Z = np.array([e["activation"][p] for e in all_data])
        color_here = point_colors[p]

        # scatter
        fig.add_trace(
            go.Scatter3d(
                x=DX_all, y=DT_all, z=Z,
                mode="markers",
                marker=dict(size=4, color=color_here),
                showlegend=False
            ),
            row=row, col=col
        )

        # DX lines
        for dx_val in np.unique(DX_all):
            mask = (DX_all == dx_val)
            fig.add_trace(
                go.Scatter3d(
                    x=DX_all[mask], y=DT_all[mask], z=Z[mask],
                    mode="lines",
                    line=dict(color=color_here, width=3),
                    showlegend=False
                ),
                row=row, col=col
            )

        # DT lines
        for dt_val in np.unique(DT_all):
            mask = (DT_all == dt_val)
            fig.add_trace(
                go.Scatter3d(
                    x=DX_all[mask], y=DT_all[mask], z=Z[mask],
                    mode="lines",
                    line=dict(color=color_here, width=3, dash="dash"),
                    showlegend=False
                ),
                row=row, col=col
            )

        scene_name = f"scene{scene_id}"
        fig.update_layout({
            scene_name: dict(
                xaxis_title="ΔX (mm)",
                yaxis_title="ΔT (ms)",
                zaxis_title="Activation (ms)",

                # Shrink box inside its subplot
                aspectmode="manual",
                aspectratio=dict(x=1, y=1, z=1),

                xaxis=dict(range=[DX_all.min() - 0.05, DX_all.max() + 0.05]),
                yaxis=dict(range=[DT_all.min() - 0.005, DT_all.max() + 0.005]),
                zaxis=dict(range=[Z.min() - 10, Z.max() + 10]),

                camera=dict(
                    eye=dict(x=-2.2, y=1.75, z=0.9)
                )
            )
        })

        scene_id += 1

    fig.update_layout(
        height=900,
        width=1600,
        title="3D Activation Surfaces"
    )

    fig.show()


# -----------------------------------------------------------
# Manual run
# -----------------------------------------------------------
if __name__ == "__main__":
    folder = Path(__file__).resolve().parent.parent / "NiedererFoam"
    print(f"[points_postProcessing] Default folder = {folder}")
    plot_3d_points_and_grid(folder)
