#!/usr/bin/env python3
"""Plot ionicHeterogeneityProbe Vm traces.

Creates:
- 2D transition-band overlays for endo-M and M-epi bands.
- 3D Vm(time, transmural distance) surface.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

ENDO_COLOR = "#1565c0"
M_CELL_COLOR = "#c2185b"
EPI_COLOR = "#f57c00"


def read_vm_traces(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    traces: dict[float, list[tuple[float, float]]] = defaultdict(list)

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            traces[float(row["t"])].append((float(row["time"]), float(row["Vm"])))

    t_values = np.array(sorted(traces), dtype=float)
    first_t = float(t_values[0])
    time_values = np.array([time for time, _ in traces[first_t]], dtype=float)
    vm = np.zeros((t_values.size, time_values.size), dtype=float)

    for i, t in enumerate(t_values):
        rows = traces[float(t)]
        vm[i, :] = np.array([value for _, value in rows], dtype=float)

    return time_values, t_values, vm


def nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def smoothstep(x: float) -> float:
    x = min(1.0, max(0.0, x))
    return x*x*(3.0 - 2.0*x)


def blend_rgba(lower_color: str, upper_color: str, weight: float) -> np.ndarray:
    lower = np.array(to_rgba(lower_color), dtype=float)
    upper = np.array(to_rgba(upper_color), dtype=float)
    return (1.0 - weight)*lower + weight*upper


def transmural_rgba(
    t: float,
    endo_m: float,
    m_epi: float,
    width: float,
) -> np.ndarray:
    endo_m_max = endo_m + width
    m_epi_max = m_epi + width

    if t <= endo_m:
        return np.array(to_rgba(ENDO_COLOR), dtype=float)

    if t < endo_m_max:
        x = (t - endo_m)/width
        return blend_rgba(ENDO_COLOR, M_CELL_COLOR, smoothstep(x))

    if t <= m_epi:
        return np.array(to_rgba(M_CELL_COLOR), dtype=float)

    if t < m_epi_max:
        x = (t - m_epi)/width
        return blend_rgba(M_CELL_COLOR, EPI_COLOR, smoothstep(x))

    return np.array(to_rgba(EPI_COLOR), dtype=float)


def transmural_facecolors(
    t_values: np.ndarray,
    n_time: int,
    endo_m: float,
    m_epi: float,
    width: float,
) -> np.ndarray:
    colors = np.array(
        [transmural_rgba(t, endo_m, m_epi, width) for t in t_values],
        dtype=float,
    )
    colors[:, 3] = 0.9
    return np.repeat(colors[:, np.newaxis, :], n_time, axis=1)


def plot_transition_panel(
    ax,
    time: np.ndarray,
    t_values: np.ndarray,
    vm: np.ndarray,
    lower_plateau_t: float,
    upper_plateau_t: float,
    interface: float,
    width: float,
    lower_label: str,
    upper_label: str,
    lower_color: str,
    upper_color: str,
    faint_color: str,
) -> None:
    band_min = interface
    band_max = interface + width

    for i, t in enumerate(t_values):
        if band_min <= t <= band_max:
            ax.plot(time, vm[i], color=faint_color, alpha=0.28, linewidth=0.9)

    lower_i = nearest_index(t_values, lower_plateau_t)
    upper_i = nearest_index(t_values, upper_plateau_t)
    interface_i = nearest_index(t_values, interface)

    ax.plot(
        time,
        vm[lower_i],
        color=lower_color,
        linewidth=2.4,
        label=f"{lower_label} t={t_values[lower_i]:.3g}",
    )
    ax.plot(
        time,
        vm[upper_i],
        color=upper_color,
        linewidth=2.4,
        label=f"{upper_label} t={t_values[upper_i]:.3g}",
    )
    ax.plot(
        time,
        vm[interface_i],
        color="black",
        linewidth=1.6,
        linestyle="--",
        label=f"interface t={t_values[interface_i]:.3g}",
    )

    ax.set_xlabel("time [ms]")
    ax.set_ylabel("Vm [mV]")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)


def make_2d_plot(
    output: Path,
    time: np.ndarray,
    t_values: np.ndarray,
    vm: np.ndarray,
    endo_m: float,
    m_epi: float,
    width: float,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

    plot_transition_panel(
        axes[0],
        time,
        t_values,
        vm,
        lower_plateau_t=endo_m,
        upper_plateau_t=min(1.0, endo_m + width),
        interface=endo_m,
        width=width,
        lower_label="Endo plateau",
        upper_label="M-cell plateau",
        lower_color=ENDO_COLOR,
        upper_color=M_CELL_COLOR,
        faint_color="#7e57c2",
    )
    axes[0].set_title("Endo to M-cell transition")

    plot_transition_panel(
        axes[1],
        time,
        t_values,
        vm,
        lower_plateau_t=m_epi,
        upper_plateau_t=min(1.0, m_epi + width),
        interface=m_epi,
        width=width,
        lower_label="M-cell plateau",
        upper_label="Epi plateau",
        lower_color=M_CELL_COLOR,
        upper_color=EPI_COLOR,
        faint_color="#ffb74d",
    )
    axes[1].set_title("M-cell to Epi transition")

    fig.suptitle("Vm traces across transmural transition bands")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def make_all_traces_plot(
    endo_m_output: Path,
    m_epi_output: Path,
    time: np.ndarray,
    t_values: np.ndarray,
    vm: np.ndarray,
    endo_m: float,
    m_epi: float,
    width: float,
) -> None:
    m_cell_t = 0.5*(endo_m + m_epi)

    plot_zone_traces(
        endo_m_output,
        time,
        t_values,
        vm,
        zone_min=0.0,
        zone_max=m_cell_t,
        reference_targets=(0.0, m_cell_t),
        reference_labels=("Endo reference", "M-cell reference"),
        reference_colors=(ENDO_COLOR, M_CELL_COLOR),
        interface_targets=(endo_m,),
        title="All computed Vm traces: Endo to M-cell area",
        endo_m=endo_m,
        m_epi=m_epi,
        width=width,
    )
    plot_zone_traces(
        m_epi_output,
        time,
        t_values,
        vm,
        zone_min=m_cell_t,
        zone_max=1.0,
        reference_targets=(m_cell_t, 1.0),
        reference_labels=("M-cell reference", "Epi reference"),
        reference_colors=(M_CELL_COLOR, EPI_COLOR),
        interface_targets=(m_epi,),
        title="All computed Vm traces: M-cell to Epi area",
        endo_m=endo_m,
        m_epi=m_epi,
        width=width,
    )


def plot_zone_traces(
    output: Path,
    time: np.ndarray,
    t_values: np.ndarray,
    vm: np.ndarray,
    zone_min: float,
    zone_max: float,
    reference_targets: tuple[float, float],
    reference_labels: tuple[str, str],
    reference_colors: tuple[str, str],
    interface_targets: tuple[float, ...],
    title: str,
    endo_m: float,
    m_epi: float,
    width: float,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))

    for i, t in enumerate(t_values):
        if not zone_min <= t <= zone_max:
            continue
        color = transmural_rgba(float(t), endo_m, m_epi, width)
        color[3] = 0.22
        ax.plot(time, vm[i], color=color, linewidth=0.8)

    for target, label, color in zip(
        reference_targets,
        reference_labels,
        reference_colors,
    ):
        i = nearest_index(t_values, target)
        ax.plot(
            time,
            vm[i],
            color=color,
            linewidth=2.8,
            label=f"{label} t={t_values[i]:.3g}",
        )

    for interface in interface_targets:
        i = nearest_index(t_values, interface)
        ax.plot(
            time,
            vm[i],
            color="black",
            linestyle="--",
            linewidth=1.3,
            label=f"interface t={t_values[i]:.3g}",
        )

    ax.set_xlabel("time [ms]")
    ax.set_ylabel("Vm [mV]")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def make_3d_plot(
    output: Path,
    time: np.ndarray,
    t_values: np.ndarray,
    vm: np.ndarray,
    endo_m: float,
    m_epi: float,
    width: float,
) -> None:
    max_time_points = 300
    max_t_points = 101
    time_idx = np.linspace(
        0,
        time.size - 1,
        min(time.size, max_time_points),
        dtype=int,
    )
    t_idx = np.linspace(
        0,
        t_values.size - 1,
        min(t_values.size, max_t_points),
        dtype=int,
    )
    time_plot = time[time_idx]
    t_plot = t_values[t_idx]
    vm_plot = vm[np.ix_(t_idx, time_idx)]
    t_grid, time_grid = np.meshgrid(t_plot, time_plot, indexing="ij")
    facecolors = transmural_facecolors(
        t_plot,
        time_plot.size,
        endo_m,
        m_epi,
        width,
    )

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(
        time_grid,
        t_grid,
        vm_plot,
        facecolors=facecolors,
        shade=False,
        linewidth=0,
        antialiased=True,
    )

    for interface, label in ((endo_m, "Endo-M interface"), (m_epi, "M-Epi interface")):
        i = nearest_index(t_plot, interface)
        ax.plot(
            time_plot,
            np.full_like(time_plot, t_plot[i]),
            vm_plot[i, :],
            color="black",
            linestyle="--",
            linewidth=1.4,
            label=label,
        )

    ax.set_xlabel("time [ms]")
    ax.set_ylabel("transmural distance t")
    ax.set_zlabel("Vm [mV]")
    ax.set_title("Vm(time, transmural distance), colored by tissue band")
    ax.legend(
        handles=[
            Patch(facecolor=ENDO_COLOR, label="Endo"),
            Patch(facecolor=M_CELL_COLOR, label="M-cell"),
            Patch(facecolor=EPI_COLOR, label="Epi"),
        ],
        loc="upper left",
        frameon=False,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def make_apd_plot(output: Path, metrics_path: Path) -> None:
    if not metrics_path.exists():
        return

    t_values: list[float] = []
    apd30: list[float] = []
    apd50: list[float] = []
    apd70: list[float] = []
    apd90: list[float] = []

    with metrics_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            t_values.append(float(row["t"]))
            apd30.append(float(row["APD30"]))
            apd50.append(float(row["APD50"]))
            apd70.append(float(row["APD70"]))
            apd90.append(float(row["APD90"]))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t_values, apd30, label="APD30", color="#43a047")
    ax.plot(t_values, apd50, label="APD50", color="#1e88e5")
    ax.plot(t_values, apd70, label="APD70", color="#8e24aa")
    ax.plot(t_values, apd90, label="APD90", color="#e53935", linewidth=2)
    ax.set_xlabel("transmural distance t")
    ax.set_ylabel("APD [ms]")
    ax.set_title("APD metrics across transmural distance")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot ionicHeterogeneityProbe CSV outputs."
    )
    parser.add_argument(
        "probe_dir",
        nargs="?",
        default="postProcessing/ionicHeterogeneityProbe",
        type=Path,
        help="Directory containing Vm_traces.csv.",
    )
    parser.add_argument("--endo-m-interface", type=float, default=0.3)
    parser.add_argument("--m-epi-interface", type=float, default=0.7)
    parser.add_argument("--transition-width", type=float, default=0.1)
    parser.add_argument(
        "--all-traces",
        action="store_true",
        help="Also write separate all-trace plots for the Endo-M and M-Epi areas.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for plots. Defaults to probe_dir.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    probe_dir = args.probe_dir
    output_dir = args.output_dir or probe_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    traces_path = probe_dir / "Vm_traces.csv"
    metrics_path = probe_dir / "AP_metrics.csv"

    time, t_values, vm = read_vm_traces(traces_path)

    make_2d_plot(
        output_dir / "Vm_transition_bands_2D.png",
        time,
        t_values,
        vm,
        args.endo_m_interface,
        args.m_epi_interface,
        args.transition_width,
    )
    if args.all_traces:
        make_all_traces_plot(
            output_dir / "Vm_endo_m_all_traces_2D.png",
            output_dir / "Vm_m_epi_all_traces_2D.png",
            time,
            t_values,
            vm,
            args.endo_m_interface,
            args.m_epi_interface,
            args.transition_width,
        )
    make_3d_plot(
        output_dir / "Vm_transmural_surface_3D.png",
        time,
        t_values,
        vm,
        args.endo_m_interface,
        args.m_epi_interface,
        args.transition_width,
    )
    make_apd_plot(output_dir / "APD_transmural_metrics.png", metrics_path)

    print(f"Wrote {output_dir / 'Vm_transition_bands_2D.png'}")
    if args.all_traces:
        print(f"Wrote {output_dir / 'Vm_endo_m_all_traces_2D.png'}")
        print(f"Wrote {output_dir / 'Vm_m_epi_all_traces_2D.png'}")
    print(f"Wrote {output_dir / 'Vm_transmural_surface_3D.png'}")
    print(f"Wrote {output_dir / 'APD_transmural_metrics.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
