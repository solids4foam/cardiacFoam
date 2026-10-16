import json
import re
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from omnidriver.postprocessing.style import apply_plotly_layout, write_plotly_html

_PROBE_LINE = re.compile(r"^# Probe (\d+) \(([^)]+)\)")

#: Every study here fixes solutionAlgorithm to implicit.
_SOLVER = "implicit"


def _read_latest_probe_row(function_object_dir: Path, field: str):
    """Return (probe_xyz, values) from a `probes` functionObject's last row, checking every instance directory."""
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


# -----------------------------------------------------------
# Load every swept case's benchmark point samples
# -----------------------------------------------------------
def load_all_point_cases(output_dir):
    output_dir = Path(output_dir)
    manifest_path = output_dir / "sweep_manifest.json"
    if not manifest_path.is_file():
        print(f"No sweep_manifest.json found in: {output_dir}")
        return []
    manifest = json.loads(manifest_path.read_text())

    all_data = []
    for case in manifest.get("cases", []):
        case_dir = output_dir / "cases" / case["case_id"]
        axes = case.get("resolved_axis_values", {})
        dx_m = axes.get("dx", axes.get("tetDx"))
        dt_s = axes.get("system/controlDict:deltaT")
        if dx_m is None or dt_s is None:
            continue
        sample = _read_latest_probe_row(case_dir / "postProcessing" / "Niedererpoints", "activationTime")
        if sample is None:
            continue
        _probes, values = sample
        # s -> ms, except the -1 never-activated sentinel, which is never converted.
        activation = np.array([v if v == -1.0 else v * 1000.0 for v in values])

        all_data.append({
            "DX": round(float(dx_m) * 1000.0, 4),
            "DT": round(float(dt_s) * 1000.0, 5),
            "solver": _SOLVER,
            "activation": activation,
        })

    print(f"Loaded {len(all_data)} swept cases from {output_dir}")
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


def plot_3d_points_and_grid(folder=".", show: bool = True):
    all_data = load_all_point_cases(folder)
    if not all_data:
        return None

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
    Solver_all = np.array([e["solver"] for e in all_data])

    fig = make_subplots(
        rows=2, cols=4,
        specs=[[{'type': 'scene'}] * 4,
               [{'type': 'scene'}] * 4],
        horizontal_spacing=0.04,
        vertical_spacing=0.08
    )

    scene_id = 1

    for idx, p in enumerate(points_to_plot):
        row = idx // 4 + 1
        col = idx % 4 + 1

        Z = np.array([e["activation"][p] for e in all_data])
        color_here = point_colors[p]

        for solver in [_SOLVER]:
            mask = (Solver_all == solver)
            if not np.any(mask): continue
            symbol = 'square'
            fig.add_trace(
                go.Scatter3d(
                    x=DX_all[mask], y=DT_all[mask], z=Z[mask],
                    mode="markers",
                    marker=dict(size=4, color=color_here, symbol=symbol),
                    showlegend=False
                ),
                row=row, col=col
            )

        for dx_val in np.unique(DX_all):
            for solver in [_SOLVER]:
                mask = (DX_all == dx_val) & (Solver_all == solver)
                if not np.any(mask): continue
                idx_sorted = np.argsort(DT_all[mask])
                fig.add_trace(
                    go.Scatter3d(
                        x=DX_all[mask][idx_sorted], y=DT_all[mask][idx_sorted], z=Z[mask][idx_sorted],
                        mode="lines",
                        line=dict(color=color_here, width=3, dash='dash'),
                        showlegend=False
                    ),
                    row=row, col=col
                )

        for dt_val in np.unique(DT_all):
            for solver in [_SOLVER]:
                mask = (DT_all == dt_val) & (Solver_all == solver)
                if not np.any(mask): continue
                idx_sorted = np.argsort(DX_all[mask])
                fig.add_trace(
                    go.Scatter3d(
                        x=DX_all[mask][idx_sorted], y=DT_all[mask][idx_sorted], z=Z[mask][idx_sorted],
                        mode="lines",
                        line=dict(color=color_here, width=3, dash='dot'),
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

    apply_plotly_layout(
        fig,
        height=900,
        width=1600,
        title="3D Activation Surfaces",
        xaxis_title="",
        yaxis_title="",
        showlegend=False,
    )

    output_html = Path(folder) / "activation_surfaces_3d.html"
    write_plotly_html(fig, output_html)
    with open(Path(folder) / "activation_surfaces_3d.json", "w") as f:
        f.write(fig.to_json())
    if show:

        fig.show()
    return output_html


def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **_: object) -> None:
    """Plot 3D activation-time surfaces over the swept dx/dt grid for the 9 benchmark probe points.

    The earliest-activated point (the stimulus reference) is excluded; writes activation_surfaces_3d.html.
    """
    del setup_root
    output_html = plot_3d_points_and_grid(output_dir, show=False)
    if output_html is None:
        return []
    return [
        {
            "path": str(output_html),
            "label": "Niederer point activation surfaces",
            "kind": "plot",
            "format": "html",
        }
    ]


# -----------------------------------------------------------
# Manual run
# -----------------------------------------------------------
if __name__ == "__main__":
    folder = Path(__file__).resolve().parent.parent
    print(f"[points_postProcessing] Default folder = {folder}")
    plot_3d_points_and_grid(folder)
