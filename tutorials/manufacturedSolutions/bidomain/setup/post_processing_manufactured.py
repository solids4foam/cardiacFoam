"""Post-processing helpers for manufactured-solution convergence tests."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import re
import sys

try:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

TUTORIALS_ROOT = Path(__file__).resolve().parents[2]
if str(TUTORIALS_ROOT) not in sys.path:
    sys.path.insert(0, str(TUTORIALS_ROOT))

from openfoam_driver.postprocessing.style import (
    configure_matplotlib_defaults,
    finalize_matplotlib_figure,
    style_matplotlib_axes,
)

RATE_FIELDS = (
    "Dimension",
    "N_lower",
    "N_higher",
    "rate_Vm",
    "rate_phiE",
    "rate_u1",
    "rate_u2",
)
FILENAME_PATTERN = re.compile(r"(\dD)_(\d+)_cells(?:_DT[^_]+)?\.dat$")
FIELD_COLORS = {
    "Linf_V": "tab:blue",
    "Linf_phiE": "tab:orange",
    "Linf_u1": "tab:green",
    "Linf_u2": "tab:red",
}
FIELD_LABELS = {
    "Linf_V": "Vm",
    "Linf_phiE": "phiE (gauge)",
    "Linf_u1": "u1",
    "Linf_u2": "u2",
}
DIMENSION_COLORS = {"1D": "tab:blue", "2D": "tab:orange", "3D": "tab:green"}


def _has_matplotlib() -> bool:
    return plt is not None


def _unique_values(rows, key):
    seen = set()
    ordered = []
    for row in rows:
        value = row[key]
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def _filter_rows(rows, **criteria):
    return [row for row in rows if all(row[field] == value for field, value in criteria.items())]


def _safe_rate(e1: float, e2: float, h1: float, h2: float) -> float:
    if any(math.isnan(value) for value in (e1, e2)):
        return float("nan")
    if e1 <= 0 or e2 <= 0 or h1 == h2:
        return float("nan")
    return math.log(e1 / e2) / math.log(h1 / h2)


def _write_csv(rows, destination: Path, fieldnames) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _field_linf(content: str, field_name: str) -> float:
    match = re.search(
        rf"^{re.escape(field_name)}\s+\S+\s+\S+\s+(\S+)",
        content,
        re.MULTILINE,
    )
    if not match:
        return float("nan")
    return float(match.group(1))


def read_error_dat_files(folder_name):
    """
    Reads all .dat files in folder_name and extracts:
        - Dimension  (1D, 2D, 3D)
        - N          (# cells)
        - Linf errors for Vm, gauge-corrected phiE, u1, u2

    Returns one row per file.
    """

    folder = Path(folder_name)
    if not folder.exists():
        print("Folder does not exist:", folder)
        return []

    files = [f for f in folder.iterdir() if f.suffix == ".dat"]
    if not files:
        print("No .dat files found in folder:", folder)
        return []

    data = []

    for f in files:
        # Expected filename format:
        #   1D_320_cells.dat
        m = FILENAME_PATTERN.match(f.name)
        if not m:
            print("Skipping unrecognized filename:", f.name)
            continue

        dimension = m.group(1)   # "1D"
        N = int(m.group(2))      # 320

        content = f.read_text()

        Linf_V = _field_linf(content, "Vm")
        Linf_phiE_raw = _field_linf(content, "phiE")
        Linf_phiE_gauge = _field_linf(content, "phiE_gauge")
        Linf_phiE = Linf_phiE_gauge
        if math.isnan(Linf_phiE):
            Linf_phiE = Linf_phiE_raw
        Linf_u1 = _field_linf(content, "u1")
        Linf_u2 = _field_linf(content, "u2")

        data.append({
            "Dimension": dimension,
            "N": N,
            "Linf_V": Linf_V,
            "Linf_phiE_raw": Linf_phiE_raw,
            "Linf_phiE": Linf_phiE,
            "Linf_u1": Linf_u1,
            "Linf_u2": Linf_u2
        })

    return sorted(data, key=lambda row: (row["Dimension"], row["N"]))


def compute_convergence_rates(rows):
    """
    Compute convergence rates for Linf errors of Vm, phiE, u1, u2.

    - Groups by Dimension (if present).
    - Sorts by N.
    - Skips pairs where N_lower == N_higher.
    """

    grouped_rows = {}
    for row in rows:
        key = row["Dimension"]
        grouped_rows.setdefault(key, []).append(row)

    convergence_rows = []
    for dimension, group_rows in sorted(grouped_rows.items()):
        ordered = sorted(group_rows, key=lambda row: row["N"])
        for lower, higher in zip(ordered, ordered[1:]):
            N1 = int(lower["N"])
            N2 = int(higher["N"])
            if N1 == N2:
                continue

            h1, h2 = 1.0 / N1, 1.0 / N2
            convergence_rows.append(
                {
                    "Dimension": dimension,
                    "N_lower": N1,
                    "N_higher": N2,
                    "rate_Vm": _safe_rate(lower["Linf_V"], higher["Linf_V"], h1, h2),
                    "rate_phiE": _safe_rate(lower["Linf_phiE"], higher["Linf_phiE"], h1, h2),
                    "rate_u1": _safe_rate(lower["Linf_u1"], higher["Linf_u1"], h1, h2),
                    "rate_u2": _safe_rate(lower["Linf_u2"], higher["Linf_u2"], h1, h2),
                }
            )

    return convergence_rows


def plot_convergence_rates(
    rows,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_matplotlib():
        print("matplotlib is not available; skipping manufactured convergence-rate plot.")
        return None

    if not rows:
        return None

    configure_matplotlib_defaults()
    labels = [
        f"{row['Dimension']} {row['N_lower']}-{row['N_higher']}"
        for row in rows
    ]
    x_positions = range(len(labels))
    rate_fields = ("rate_Vm", "rate_phiE", "rate_u1", "rate_u2")
    width = 0.16
    offsets = [
        width*(index - (len(rate_fields) - 1)/2)
        for index in range(len(rate_fields))
    ]

    fig, ax = plt.subplots(figsize=(max(9, len(labels)*1.35), 5.5))
    for offset, field in zip(offsets, rate_fields):
        values = [row.get(field, float("nan")) for row in rows]
        label = {
            "rate_Vm": "Vm",
            "rate_phiE": "phiE (gauge)",
            "rate_u1": "u1",
            "rate_u2": "u2",
        }[field]
        ax.bar([x + offset for x in x_positions], values, width=width, label=label)

    ax.axhline(1.0, color="0.35", linewidth=0.8, linestyle="--")
    ax.axhline(2.0, color="0.35", linewidth=0.8, linestyle=":")
    ax.set_xticks(list(x_positions))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    style_matplotlib_axes(
        ax,
        title="Manufactured-solution observed convergence rates",
        xlabel="Refinement pair",
        ylabel="Observed rate",
        grid_kwargs={"axis": "y", "alpha": 0.25},
    )
    ax.legend(ncols=4)
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def _finalize_axis_legend(ax) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend()



def _plot_vm_across_dimensions_on_axis(ax, rows) -> bool:
    plotted = False
    for dimension in _unique_values(rows, "Dimension"):
        dimension_rows = _filter_rows(rows, Dimension=dimension)

        plotted_rows = dimension_rows
        if not plotted_rows:
            continue

        ax.loglog(
            [row["N"] for row in plotted_rows],
            [row["Linf_V"] for row in plotted_rows],
            marker="s",
            linestyle="--",
            color=DIMENSION_COLORS.get(dimension, "black"),
            label=f"{dimension}",
        )
        plotted = True

    if not plotted:
        ax.text(
            0.5,
            0.5,
            "No Vm data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color="0.4",
        )

    style_matplotlib_axes(
        ax,
        title="Linf Error of Vm across dimensions",
        xlabel="Number of cells (N)",
        ylabel="Linf Error (Vm)",
        legend=False,
        grid_kwargs={"which": "both", "ls": "--", "alpha": 0.6},
    )
    _finalize_axis_legend(ax)
    return plotted


def _plot_dimension_errors_on_axis(ax, rows, dimension: str) -> bool:
    """Plot the Linf errors of Vm, phiE, u1 and u2 vs N for a single dimension."""
    dimension_rows = _filter_rows(rows, Dimension=dimension)

    plotted = False
    for field, marker in (
        ("Linf_V", "o"),
        ("Linf_phiE", "s"),
        ("Linf_u1", "^"),
        ("Linf_u2", "d"),
    ):
        points = [
            (row["N"], row[field])
            for row in dimension_rows
            if not math.isnan(row.get(field, float("nan")))
        ]
        if not points:
            continue

        ax.loglog(
            [N for N, _ in points],
            [value for _, value in points],
            marker=marker,
            color=FIELD_COLORS.get(field, "black"),
            label=FIELD_LABELS.get(field, field),
        )
        plotted = True

    if not plotted:
        ax.text(
            0.5,
            0.5,
            f"No {dimension} data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color="0.4",
        )

    style_matplotlib_axes(
        ax,
        title=f"Linf errors ({dimension})",
        xlabel="Number of cells (N)",
        ylabel="Linf Error",
        legend=False,
        grid_kwargs={"which": "both", "ls": "--", "alpha": 0.6},
    )
    _finalize_axis_legend(ax)
    return plotted


def plot_Vm_across_dimensions(
    rows,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    """
    Plot Linf_V (Vm error) vs N across all dimensions.
    """
    if not _has_matplotlib():
        print("matplotlib is not available; skipping manufactured Vm plot.")
        return None

    configure_matplotlib_defaults()

    fig, ax = plt.subplots(figsize=(8, 6))

    _plot_vm_across_dimensions_on_axis(ax, rows)
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def plot_summary_dashboard(
    rows,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    """Create a 2x2 manufactured-summary dashboard."""
    if not _has_matplotlib():
        print("matplotlib is not available; skipping manufactured summary dashboard.")
        return None

    configure_matplotlib_defaults()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Manufactured-solution convergence summary", fontsize=14)

    _plot_vm_across_dimensions_on_axis(axes[0, 0], rows)
    _plot_dimension_errors_on_axis(axes[0, 1], rows, "1D")
    _plot_dimension_errors_on_axis(axes[1, 0], rows, "2D")
    _plot_dimension_errors_on_axis(axes[1, 1], rows, "3D")

    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **_: object) -> list[dict]:
    del setup_root
    output_path = Path(output_dir)

    error_rows = read_error_dat_files(output_dir)
    if not error_rows:
        print(f"No .dat files found to post-process in: {output_dir}")
        return []

    print("\nConvergence rates:")
    convergence_rates = compute_convergence_rates(error_rows)
    if convergence_rates:
        for row in convergence_rates:
            print(row)
    else:
        print("No convergence-rate pairs were found.")

    rates_csv = output_path / "manufactured_convergence_rates.csv"
    _write_csv(convergence_rates, rates_csv, RATE_FIELDS)

    artifacts = [
        {
            "path": str(rates_csv),
            "label": "Manufactured convergence rates table",
            "kind": "table",
            "format": "csv",
        }
    ]

    if not _has_matplotlib():
        print("matplotlib is not available; generated CSV only.")
        return artifacts

    vm_plot = plot_Vm_across_dimensions(
        error_rows,
        save_path=output_path / "manufactured_vm_across_dimensions.png",
        show=False,
    )
    summary_plot = plot_summary_dashboard(
        error_rows,
        save_path=output_path / "manufactured_summary_dashboard.png",
        show=False,
    )
    rate_plot = plot_convergence_rates(
        convergence_rates,
        save_path=output_path / "manufactured_convergence_rates.png",
        show=False,
    )
    if vm_plot is not None:
        artifacts.append(
            {
                "path": str(vm_plot),
                "label": "Manufactured Vm across dimensions",
                "kind": "plot",
                "format": "png",
            }
        )
    if summary_plot is not None:
        artifacts.append(
            {
                "path": str(summary_plot),
                "label": "Manufactured summary dashboard",
                "kind": "plot",
                "format": "png",
            }
        )
    if rate_plot is not None:
        artifacts.append(
            {
                "path": str(rate_plot),
                "label": "Manufactured convergence rates",
                "kind": "plot",
                "format": "png",
            }
        )

    return artifacts
