#!/usr/bin/env python3
"""plot_niederer_bueno_orovio_batched.py

Produce four matplotlib plots comparing activation times across modes for the
NiedererEtAl2011 Niederer benchmark:

1. line_activation_profile.png
   Arc-length vs activation time for CPU and all batched modes.

2. line_abs_error_vs_arclength.png
   |mode − cpu| vs arc-length for each batched mode.

3. scatter_points_wall_time_vs_max_error.png
   Mean wall time (x) vs max point activation error (y), one dot per mode.

4. bar_point_errors.png
   Per-benchmark-point absolute error, grouped by probe index, one bar per mode.

Usage
-----
    python3 plot_niederer_bueno_orovio_batched.py \
        <run-root> <output-dir> <mode> [<mode> ...]
"""

import math
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="matplotlib-"))

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# ---------------------------------------------------------------------------
# Re-use the probe reader from the comparison script (same directory)
# ---------------------------------------------------------------------------

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from compare_niederer_bueno_orovio_batched import (  # noqa: E402
    read_activation_times,
    read_wall_time,
    run_dirs,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Colour palette consistent with the singleCell plotter
MODE_COLORS = {
    "cpu":           "#333333",
    "batched_euler": "#1f77b4",
    "batched_heun":  "#ff7f0e",
    "batched_rl":    "#2ca02c",
    "batched_soa":   "#d62728",
}

# Niederer benchmark point labels (order matches Niedererpoints probeLocations)
POINT_LABELS = [
    "P0 (0,0,7)",
    "P1 (0,0,0)",
    "P2 (20,0,7)",
    "P3 (20,0,0)",
    "P4 (0,3,7)",
    "P5 (0,3,0)",
    "P6 (20,3,7)",
    "P7 (20,3,0)",
    "P8 (10,1.5,3.5)",
]


def _color(mode: str) -> str:
    return MODE_COLORS.get(mode, "#888888")


def _latest_activation(run_root: Path, mode: str, probe_name: str) -> list[float]:
    """Return activation times from the last run of *mode* for *probe_name*."""
    mode_dir = run_root / mode
    runs = run_dirs(mode_dir)
    last_run = runs[-1]
    return read_activation_times(last_run / "postProcessing" / probe_name)


def _mean_wall_time(run_root: Path, mode: str) -> float:
    mode_dir = run_root / mode
    runs = run_dirs(mode_dir)
    wts = [read_wall_time(r) for r in runs]
    return sum(wts) / len(wts)


def _arc_lengths(n_points: int) -> list[float]:
    """Compute cumulative arc-length (mm) for the Niedererlines diagonal.

    The diagonal runs from (0, 0, 7) to (20, 3, 0) in mm coordinates, with
    n_points equally-spaced probe locations.
    """
    # Direction vector
    dx, dy, dz = 20.0, 3.0, -7.0
    total = math.sqrt(dx**2 + dy**2 + dz**2)
    step = total / (n_points - 1) if n_points > 1 else 0.0
    return [i * step for i in range(n_points)]


# ---------------------------------------------------------------------------
# Plot 1 – line activation profile
# ---------------------------------------------------------------------------

def plot_line_activation_profile(
    run_root: Path, modes: list[str], output_dir: Path
) -> None:
    cpu_line = _latest_activation(run_root, "cpu", "Niedererlines")
    n = len(cpu_line)
    arclens = _arc_lengths(n)
    arclens_mm = arclens  # already in mm

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(
        arclens_mm,
        [v * 1000 for v in cpu_line],
        color=_color("cpu"),
        linewidth=2.0,
        label="cpu",
        zorder=5,
    )

    for mode in modes:
        if mode == "cpu":
            continue
        vals = _latest_activation(run_root, mode, "Niedererlines")
        ax.plot(
            arclens_mm,
            [v * 1000 for v in vals],
            color=_color(mode),
            linewidth=1.5,
            linestyle="--",
            label=mode,
        )

    ax.set_xlabel("Arc-length along diagonal (mm)")
    ax.set_ylabel("Activation time (ms)")
    ax.set_title("Niederer benchmark – line activation profile (cpu vs batched)")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    out = output_dir / "line_activation_profile.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"  Saved {out}")


# ---------------------------------------------------------------------------
# Plot 2 – line absolute error vs arc-length
# ---------------------------------------------------------------------------

def plot_line_abs_error(
    run_root: Path, modes: list[str], output_dir: Path
) -> None:
    cpu_line = _latest_activation(run_root, "cpu", "Niedererlines")
    n = len(cpu_line)
    arclens_mm = _arc_lengths(n)

    fig, ax = plt.subplots(figsize=(10, 5))

    for mode in modes:
        if mode == "cpu":
            continue
        vals = _latest_activation(run_root, mode, "Niedererlines")
        errs_ms = [abs(v - r) * 1000 for r, v in zip(cpu_line, vals)]
        ax.plot(
            arclens_mm,
            errs_ms,
            color=_color(mode),
            linewidth=1.5,
            label=mode,
        )
        # Mark peak
        max_i = max(range(n), key=errs_ms.__getitem__)
        ax.scatter([arclens_mm[max_i]], [errs_ms[max_i]], s=30, color=_color(mode), zorder=4)

    ax.set_xlabel("Arc-length along diagonal (mm)")
    ax.set_ylabel("|mode − cpu|  activation time (ms)")
    ax.set_title("Niederer benchmark – line absolute error vs CPU reference")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    out = output_dir / "line_abs_error_vs_arclength.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"  Saved {out}")


# ---------------------------------------------------------------------------
# Plot 3 – scatter: wall time vs max point error
# ---------------------------------------------------------------------------

def plot_scatter_wall_vs_error(
    run_root: Path, modes: list[str], output_dir: Path
) -> None:
    cpu_points = _latest_activation(run_root, "cpu", "Niedererpoints")

    scatter_x: list[float] = []
    scatter_y: list[float] = []
    labels: list[str] = []
    colors: list[str] = []

    for mode in modes:
        if mode == "cpu":
            continue
        vals = _latest_activation(run_root, mode, "Niedererpoints")
        errs_ms = [abs(v - r) * 1000 for r, v in zip(cpu_points, vals)]
        max_err = max(errs_ms)
        scatter_x.append(_mean_wall_time(run_root, mode))
        scatter_y.append(max_err)
        labels.append(mode)
        colors.append(_color(mode))

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.scatter(scatter_x, scatter_y, s=80, c=colors, zorder=3)

    for x_val, y_val, label in zip(scatter_x, scatter_y, labels):
        ax.annotate(
            label,
            (x_val, y_val),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xlabel("Mean cardiacFoam wall time (s)")
    ax.set_ylabel("Max |mode − cpu| activation time (ms)  [benchmark points]")
    ax.set_title("Niederer benchmark – speed vs accuracy tradeoff")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    out = output_dir / "scatter_points_wall_time_vs_max_error.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"  Saved {out}")


# ---------------------------------------------------------------------------
# Plot 4 – per-point error bar chart
# ---------------------------------------------------------------------------

def plot_bar_point_errors(
    run_root: Path, modes: list[str], output_dir: Path
) -> None:
    cpu_points = _latest_activation(run_root, "cpu", "Niedererpoints")
    n_points = len(cpu_points)

    batched_modes = [m for m in modes if m != "cpu"]
    n_modes = len(batched_modes)

    if n_modes == 0:
        print("  No batched modes to plot – skipping bar chart.")
        return

    x = np.arange(n_points)
    width = 0.8 / n_modes

    fig, ax = plt.subplots(figsize=(max(10, n_points * 1.2), 5.5))

    for i, mode in enumerate(batched_modes):
        vals = _latest_activation(run_root, mode, "Niedererpoints")
        errs_ms = [abs(v - r) * 1000 for r, v in zip(cpu_points, vals)]
        offset = (i - n_modes / 2.0 + 0.5) * width
        ax.bar(
            x + offset,
            errs_ms,
            width,
            label=mode,
            color=_color(mode),
            alpha=0.85,
        )

    labels = POINT_LABELS[:n_points] if n_points <= len(POINT_LABELS) else [
        f"P{i}" for i in range(n_points)
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel("|mode − cpu|  activation time (ms)")
    ax.set_title("Niederer benchmark – per-point error vs CPU reference")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.25, axis="y")
    fig.tight_layout()
    out = output_dir / "bar_point_errors.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"  Saved {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    if len(sys.argv) < 4:
        print(
            "Usage: plot_niederer_bueno_orovio_batched.py "
            "<run-root> <output-dir> <mode> [<mode> ...]",
            file=sys.stderr,
        )
        return 1

    run_root   = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    modes      = sys.argv[3:]

    if "cpu" not in modes:
        raise RuntimeError("Mode list must include 'cpu' as the reference")

    output_dir.mkdir(parents=True, exist_ok=True)

    print("Generating line activation profile ...")
    plot_line_activation_profile(run_root, modes, output_dir)

    print("Generating line absolute error plot ...")
    plot_line_abs_error(run_root, modes, output_dir)

    batched_modes = [m for m in modes if m != "cpu"]
    if batched_modes:
        print("Generating wall-time vs error scatter ...")
        plot_scatter_wall_vs_error(run_root, modes, output_dir)

        print("Generating per-point error bar chart ...")
        plot_bar_point_errors(run_root, modes, output_dir)
    else:
        print("Only 'cpu' mode present – skipping scatter and bar charts.")

    print(f"\nAll plots written to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
