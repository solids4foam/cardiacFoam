#!/usr/bin/env python3
"""plot_substep_sweep.py

Speed-vs-accuracy tradeoff plots for the euler and soa substep sweep.

Reads two JSON files:
  - substep_sweep_metrics.json   (from run_substep_sweep.sh, substeps 50/20/10/5)
  - comparison_metrics.json      (from run_...batched_comparison.sh, has sub100
                                   stored as 'batched_euler' and 'batched_soa')

The sub100 data points are pulled from the first comparison so we never need
to re-run them. The CPU reference wall time is also taken from there.

Produces four plots:

1. substep_scatter_wall_vs_point_error.png
   X = mean cardiacFoam wall time (s)
   Y = max benchmark-point activation error vs CPU (ms)
   One connected line per model (euler, soa). Each point labelled N=<substeps>.
   Dotted vertical line = CPU reference wall time.

2. substep_scatter_wall_vs_line_error.png
   Same but Y = max line-probe error.

3. substep_convergence_point_error.png
   X = substep count (log scale)
   Y = max point error. Shows convergence as substeps increase.

4. substep_convergence_line_error.png
   Same but Y = max line error.

Usage
-----
    python3 plot_substep_sweep.py \\
        <substep_sweep_metrics.json> \\
        <comparison_metrics.json> \\
        <output-dir> \\
        <substep> [<substep> ...]

    The substep list should be the ones in the sweep (e.g. 50 20 10 5).
    Sub100 is pulled automatically from comparison_metrics.json.
"""

import json
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="matplotlib-"))

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
MODEL_STYLE = {
    "euler": {"color": "#1f77b4", "marker": "o", "label": "Batched Euler"},
    "soa":   {"color": "#d62728", "marker": "s", "label": "Batched SoA (Euler)"},
}

# Names used for sub100 in the first comparison JSON
FIRST_COMPARISON_NAME = {
    "euler": "batched_euler",
    "soa":   "batched_soa",
}


def _mode_name(family: str, substeps: int) -> str:
    """Return the mode key used in the sweep JSON, e.g. 'euler_sub050'."""
    return f"{family}_sub{substeps:03d}"


# ---------------------------------------------------------------------------
# JSON loading
# ---------------------------------------------------------------------------

def _load(path: Path) -> dict:
    return json.loads(path.read_text())


def _safe_get(metrics: dict, mode: str, key: str, subkey: str):
    return (
        metrics.get("cases", {}).get(mode, {})
               .get(key, {}).get(subkey)
    )


def _accuracy(metrics: dict, mode: str, probe: str, metric: str):
    return (
        metrics.get("accuracy_vs_cpu", {}).get(mode, {})
               .get(probe, {}).get(metric)
    )


# ---------------------------------------------------------------------------
# Merge data from both JSON files into per-family series
# ---------------------------------------------------------------------------

def collect_series(
    sweep_metrics: dict,
    first_metrics: dict,
    substep_list: list[int],
) -> dict:
    """Return a dict: family → {substeps, wall_times, pt_errors_ms, ln_errors_ms}.

    sub100 is pulled from first_metrics (batched_euler / batched_soa).
    All other substeps are pulled from sweep_metrics (euler_subXXX / soa_subXXX).
    Points with missing data are skipped with a warning.
    """
    series: dict[str, dict] = {
        f: {"substeps": [], "wall_times": [], "pt_errors_ms": [], "ln_errors_ms": []}
        for f in MODEL_STYLE
    }

    all_substeps = sorted(set(substep_list) | {100})

    for substeps in all_substeps:
        for family in MODEL_STYLE:
            if substeps == 100:
                # Pull from first comparison JSON
                mode = FIRST_COMPARISON_NAME[family]
                source = first_metrics
            else:
                mode = _mode_name(family, substeps)
                source = sweep_metrics

            wall_time = _safe_get(source, mode, "external_wall_time_s", "mean")
            pt_err    = _accuracy(source, mode, "points", "max_abs_diff_ms")
            ln_err    = _accuracy(source, mode, "line",   "max_abs_diff_ms")

            if wall_time is None or pt_err is None or ln_err is None:
                print(f"  Warning: incomplete data for {family} sub{substeps:03d} "
                      f"(mode='{mode}') – skipping")
                continue

            series[family]["substeps"].append(substeps)
            series[family]["wall_times"].append(wall_time)
            series[family]["pt_errors_ms"].append(pt_err)
            series[family]["ln_errors_ms"].append(ln_err)

    return series


def cpu_wall_time(first_metrics: dict) -> float | None:
    return _safe_get(first_metrics, "cpu", "external_wall_time_s", "mean")


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _annotate(ax, x_vals, y_vals, substep_vals, color):
    for x, y, s in zip(x_vals, y_vals, substep_vals):
        ax.annotate(
            f"N={s}",
            (x, y),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=7,
            color=color,
        )


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"  Saved {path}")


def _legend_handles(families: list[str]) -> list:
    handles = []
    for family in families:
        style = MODEL_STYLE[family]
        handles.append(
            mlines.Line2D([], [],
                          color=style["color"],
                          marker=style["marker"],
                          linewidth=1.5,
                          label=style["label"])
        )
    return handles


# ---------------------------------------------------------------------------
# Plot 1 + 2: wall time vs error (tradeoff scatter with connected lines)
# ---------------------------------------------------------------------------

def plot_wall_vs_error(
    series: dict,
    output_dir: Path,
    error_key: str,
    ylabel: str,
    filename: str,
    ref_wall_time: float | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))

    active_families = []
    for family, style in MODEL_STYLE.items():
        s = series[family]
        if not s["substeps"]:
            continue

        # Sort by substep count so the line goes low → high
        order = sorted(range(len(s["substeps"])), key=lambda i: s["substeps"][i])
        xs = [s["wall_times"][i]  for i in order]
        ys = [s[error_key][i]     for i in order]
        ns = [s["substeps"][i]    for i in order]

        ax.plot(xs, ys, color=style["color"], linewidth=1.5, zorder=3)
        ax.scatter(xs, ys, marker=style["marker"], s=55,
                   color=style["color"], zorder=4)
        _annotate(ax, xs, ys, ns, style["color"])
        active_families.append(family)

    if ref_wall_time is not None:
        ax.axvline(ref_wall_time, color="#333333", linewidth=1.0,
                   linestyle=":", alpha=0.8)
        ax.annotate(
            "CPU ref",
            (ref_wall_time, ax.get_ylim()[1] if ax.get_ylim()[1] != 1.0 else 1.0),
            xytext=(4, -12),
            textcoords="offset points",
            fontsize=7,
            color="#333333",
        )
        cpu_h = mlines.Line2D([], [], color="#333333", linewidth=1.0,
                              linestyle=":", label="CPU reference")
        handles = _legend_handles(active_families) + [cpu_h]
    else:
        handles = _legend_handles(active_families)

    ax.set_xlabel("Mean cardiacFoam wall time (s)")
    ax.set_ylabel(ylabel)
    ax.set_title("Niederer benchmark – substep tradeoff: speed vs accuracy")
    ax.legend(handles=handles, loc="best")
    ax.grid(True, alpha=0.25)
    _save(fig, output_dir / filename)


# ---------------------------------------------------------------------------
# Plot 3 + 4: substep count vs error (convergence)
# ---------------------------------------------------------------------------

def plot_convergence(
    series: dict,
    output_dir: Path,
    error_key: str,
    ylabel: str,
    filename: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))

    active_families = []
    for family, style in MODEL_STYLE.items():
        s = series[family]
        if not s["substeps"]:
            continue

        order = sorted(range(len(s["substeps"])), key=lambda i: s["substeps"][i])
        xs = [s["substeps"][i]  for i in order]
        ys = [s[error_key][i]   for i in order]

        ax.plot(xs, ys, color=style["color"], linewidth=1.5, zorder=3)
        ax.scatter(xs, ys, marker=style["marker"], s=55,
                   color=style["color"], zorder=4)
        active_families.append(family)

    ax.set_xscale("log")
    ax.set_xlabel("Number of substeps  (log scale)")
    ax.set_ylabel(ylabel)
    ax.set_title("Niederer benchmark – activation error vs substep count")
    ax.legend(handles=_legend_handles(active_families), loc="best")
    ax.grid(True, alpha=0.25, which="both")
    _save(fig, output_dir / filename)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    if len(sys.argv) < 5:
        print(
            "Usage: plot_substep_sweep.py "
            "<substep_sweep_metrics.json> "
            "<comparison_metrics.json> "
            "<output-dir> "
            "<substep> [<substep> ...]",
            file=sys.stderr,
        )
        return 1

    sweep_json      = Path(sys.argv[1])
    comparison_json = Path(sys.argv[2])
    output_dir      = Path(sys.argv[3])
    substep_list    = [int(s) for s in sys.argv[4:]]

    sweep_metrics  = _load(sweep_json)
    first_metrics  = _load(comparison_json)

    print(f"Sweep metrics    : {sweep_json}")
    print(f"Comparison source: {comparison_json}  (sub100 + CPU reference)")

    series   = collect_series(sweep_metrics, first_metrics, substep_list)
    ref_wall = cpu_wall_time(first_metrics)
    if ref_wall:
        print(f"CPU reference wall time: {ref_wall:.2f} s")

    output_dir.mkdir(parents=True, exist_ok=True)

    print("Generating wall-time vs point error tradeoff ...")
    plot_wall_vs_error(
        series, output_dir,
        error_key="pt_errors_ms",
        ylabel="Max benchmark-point error vs CPU (ms)",
        filename="substep_scatter_wall_vs_point_error.png",
        ref_wall_time=ref_wall,
    )

    print("Generating wall-time vs line error tradeoff ...")
    plot_wall_vs_error(
        series, output_dir,
        error_key="ln_errors_ms",
        ylabel="Max line-probe error vs CPU (ms)",
        filename="substep_scatter_wall_vs_line_error.png",
        ref_wall_time=ref_wall,
    )

    print("Generating substep convergence (point error) ...")
    plot_convergence(
        series, output_dir,
        error_key="pt_errors_ms",
        ylabel="Max benchmark-point error vs CPU (ms)",
        filename="substep_convergence_point_error.png",
    )

    print("Generating substep convergence (line error) ...")
    plot_convergence(
        series, output_dir,
        error_key="ln_errors_ms",
        ylabel="Max line-probe error vs CPU (ms)",
        filename="substep_convergence_line_error.png",
    )

    print(f"\nAll substep plots written to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
