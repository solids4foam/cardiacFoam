#!/usr/bin/env python3
#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     post_processing_manufactured_monodomain_total_lagrangian_em
#----------------------------------------------------------------------------#
"""Post-processing helpers for manufactured electromechanics convergence tests."""

from __future__ import annotations

import csv
import json
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
    "Solver",
    "N_lower",
    "N_higher",
    "rate_Vm",
    "rate_D",
    "rate_lambda",
    "rate_Ta",
)
SUMMARY_FIELDS = (
    "Dimension",
    "N",
    "Solver",
    "L1_Vm",
    "L2_Vm",
    "Linf_Vm",
    "L1_D",
    "L2_D",
    "Linf_D",
    "L1_lambda",
    "L2_lambda",
    "Linf_lambda",
    "L1_Ta",
    "L2_Ta",
    "Linf_Ta",
)
FILENAME_PATTERN = re.compile(r"(\dD)_(\d+)_cells_(explicit|implicit)")
FIELD_LINE_PATTERN = re.compile(
    r"^(Vm|D|lambda|Ta)\s+(\S+)\s+(\S+)\s+(\S+)\s*$",
    re.MULTILINE,
)
FIELD_COLORS = {
    "Linf_Vm": "tab:blue",
    "Linf_D": "tab:orange",
    "Linf_lambda": "tab:green",
    "Linf_Ta": "tab:red",
}
FIELD_LABELS = {
    "Linf_Vm": "Vm",
    "Linf_D": "D (diagnostic)",
    "Linf_lambda": "lambda (diagnostic)",
    "Linf_Ta": "Ta (diagnostic)",
}
DIMENSION_COLORS = {"1D": "tab:blue", "2D": "tab:orange", "3D": "tab:green"}


def _has_matplotlib() -> bool:
    return plt is not None


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


def _load_expected_filenames(output_dir):
    return None


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


def read_error_dat_files(folder_name, *, expected_filenames: set[str] | None = None):
    folder = Path(folder_name)
    if not folder.exists():
        print("Folder does not exist:", folder)
        return []

    files = [f for f in folder.iterdir() if f.suffix == ".dat"]
    if expected_filenames is not None:
        files = [f for f in files if f.name in expected_filenames]

    data = []
    for f in files:
        match = FILENAME_PATTERN.match(f.name)
        if not match:
            continue

        values = {}
        for field, l1, l2, linf in FIELD_LINE_PATTERN.findall(f.read_text()):
            values[field] = (float(l1), float(l2), float(linf))

        if {"Vm", "D", "lambda", "Ta"} - values.keys():
            print("Skipping incomplete electromechanics summary:", f.name)
            continue

        data.append(
            {
                "Dimension": match.group(1),
                "N": int(match.group(2)),
                "Solver": match.group(3),
                "L1_Vm": values["Vm"][0],
                "L2_Vm": values["Vm"][1],
                "Linf_Vm": values["Vm"][2],
                "L1_D": values["D"][0],
                "L2_D": values["D"][1],
                "Linf_D": values["D"][2],
                "L1_lambda": values["lambda"][0],
                "L2_lambda": values["lambda"][1],
                "Linf_lambda": values["lambda"][2],
                "L1_Ta": values["Ta"][0],
                "L2_Ta": values["Ta"][1],
                "Linf_Ta": values["Ta"][2],
            }
        )

    return sorted(data, key=lambda row: (row["Dimension"], row["Solver"], row["N"]))


def compute_convergence_rates(rows):
    grouped_rows = {}
    for row in rows:
        key = (row["Dimension"], row["Solver"])
        grouped_rows.setdefault(key, []).append(row)

    convergence_rows = []
    for (dimension, solver_type), group_rows in sorted(grouped_rows.items()):
        ordered = sorted(group_rows, key=lambda row: row["N"])
        for lower, higher in zip(ordered, ordered[1:]):
            n1 = int(lower["N"])
            n2 = int(higher["N"])
            if n1 == n2:
                continue
            h1, h2 = 1.0 / n1, 1.0 / n2
            convergence_rows.append(
                {
                    "Dimension": dimension,
                    "Solver": solver_type,
                    "N_lower": n1,
                    "N_higher": n2,
                    "rate_Vm": _safe_rate(lower["Linf_Vm"], higher["Linf_Vm"], h1, h2),
                    "rate_D": _safe_rate(lower["Linf_D"], higher["Linf_D"], h1, h2),
                    "rate_lambda": _safe_rate(
                        lower["Linf_lambda"], higher["Linf_lambda"], h1, h2
                    ),
                    "rate_Ta": _safe_rate(lower["Linf_Ta"], higher["Linf_Ta"], h1, h2),
                }
            )
    return convergence_rows


def plot_errors_by_dimension(rows, *, save_dir: Path, show: bool = False) -> list[Path]:
    configure_matplotlib_defaults()
    saved_paths: list[Path] = []

    for dimension in _unique_values(rows, "Dimension"):
        dimension_rows = _filter_rows(rows, Dimension=dimension)
        if not dimension_rows:
            continue

        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        for field_name in ("Linf_Vm", "Linf_D", "Linf_lambda", "Linf_Ta"):
            for solver_type in _unique_values(dimension_rows, "Solver"):
                solver_rows = _filter_rows(dimension_rows, Solver=solver_type)
                ax.loglog(
                    [row["N"] for row in solver_rows],
                    [row[field_name] for row in solver_rows],
                    marker="o",
                    color=FIELD_COLORS[field_name],
                    linestyle="-" if solver_type == "implicit" else "--",
                    label=f"{FIELD_LABELS[field_name]} ({solver_type})",
                )

        style_matplotlib_axes(
            ax,
            title=f"Manufactured Electromechanics Errors ({dimension})",
            xlabel="Number of cells (N)",
            ylabel="Linf error",
        )
        ax.legend()

        output_path = save_dir / f"manufactured_electromechanics_errors_{dimension.lower()}.png"
        finalize_matplotlib_figure(fig, save_path=output_path, show=show)
        saved_paths.append(output_path)

    return saved_paths


def plot_vm_across_dimensions(rows, *, save_path: Path, show: bool = False) -> Path | None:
    configure_matplotlib_defaults()
    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    plotted = False
    for dimension in _unique_values(rows, "Dimension"):
        dimension_rows = _filter_rows(rows, Dimension=dimension, Solver="implicit")
        if not dimension_rows:
            continue

        ax.loglog(
            [row["N"] for row in dimension_rows],
            [row["Linf_Vm"] for row in dimension_rows],
            marker="o",
            color=DIMENSION_COLORS.get(dimension, "tab:blue"),
            linestyle="-",
            label=dimension,
        )
        plotted = True

    if not plotted:
        plt.close(fig)
        return None

    style_matplotlib_axes(
        ax,
        title="Manufactured Electromechanics Vm Across Dimensions",
        xlabel="Number of cells (N)",
        ylabel="Vm Linf error",
    )
    ax.legend()
    finalize_matplotlib_figure(fig, save_path=save_path, show=show)
    return save_path


def plot_summary_dashboard(rows, *, save_path: Path, show: bool = False) -> Path | None:
    configure_matplotlib_defaults()
    dimensions = _unique_values(rows, "Dimension")
    if not dimensions:
        return None

    fig, axes = plt.subplots(1, len(dimensions), figsize=(5.5 * len(dimensions), 5.0))
    if len(dimensions) == 1:
        axes = [axes]

    for ax, dimension in zip(axes, dimensions):
        dimension_rows = _filter_rows(rows, Dimension=dimension, Solver="implicit")
        for field_name in ("Linf_Vm", "Linf_D", "Linf_lambda", "Linf_Ta"):
            ax.loglog(
                [row["N"] for row in dimension_rows],
                [row[field_name] for row in dimension_rows],
                marker="o",
                color=FIELD_COLORS[field_name],
                linestyle="-",
                label=FIELD_LABELS[field_name],
            )
        style_matplotlib_axes(
            ax,
            title=dimension,
            xlabel="Number of cells (N)",
            ylabel="Linf error",
        )
        ax.legend()

    fig.suptitle(
        "Manufactured Electromechanics Summary\n"
        "Vm, D, lambda, and Ta are rigorous MMS targets.",
        fontsize=12,
    )
    finalize_matplotlib_figure(fig, save_path=save_path, show=show)
    return save_path


def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **_: object) -> list[dict]:
    del setup_root
    output_path = Path(output_dir)
    expected_filenames = _load_expected_filenames(output_path)
    if expected_filenames is not None:
        available_filenames = {
            path.name
            for path in output_path.glob("*.dat")
            if FILENAME_PATTERN.match(path.name)
        }
        unexpected = sorted(available_filenames - expected_filenames)
        missing = sorted(expected_filenames - available_filenames)
        if unexpected:
            print("Ignoring stale manufactured outputs:", ", ".join(unexpected))
        if missing:
            print("Missing expected manufactured outputs:", ", ".join(missing))

    error_rows = read_error_dat_files(output_dir, expected_filenames=expected_filenames)
    if not error_rows:
        print(f"No .dat files found to post-process in: {output_dir}")
        return []

    print(
        "\nElectromechanics convergence note: Vm, D, lambda, and Ta are rigorous MMS targets."
    )

    convergence_rates = compute_convergence_rates(error_rows)
    print("\nConvergence rates:")
    if convergence_rates:
        for row in convergence_rates:
            print(row)
    else:
        print("No convergence-rate pairs were found.")

    summary_csv = output_path / "manufactured_electromechanics_summary.csv"
    rates_csv = output_path / "manufactured_electromechanics_convergence_rates.csv"
    _write_csv(error_rows, summary_csv, SUMMARY_FIELDS)
    _write_csv(convergence_rates, rates_csv, RATE_FIELDS)

    artifacts = [
        {
            "path": str(summary_csv),
            "label": "Manufactured electromechanics summary table",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(rates_csv),
            "label": "Manufactured electromechanics convergence rates table",
            "kind": "table",
            "format": "csv",
        },
    ]

    if not _has_matplotlib():
        print("matplotlib is not available; generated CSV only.")
        return artifacts

    vm_plot = plot_vm_across_dimensions(
        error_rows,
        save_path=output_path / "manufactured_electromechanics_vm_across_dimensions.png",
        show=False,
    )
    summary_plot = plot_summary_dashboard(
        error_rows,
        save_path=output_path / "manufactured_electromechanics_summary_dashboard.png",
        show=False,
    )
    error_plots = plot_errors_by_dimension(
        error_rows,
        save_dir=output_path,
        show=False,
    )

    if vm_plot is not None:
        artifacts.append(
            {
                "path": str(vm_plot),
                "label": "Manufactured electromechanics Vm across dimensions",
                "kind": "plot",
                "format": "png",
            }
        )
    if summary_plot is not None:
        artifacts.append(
            {
                "path": str(summary_plot),
                "label": "Manufactured electromechanics summary dashboard",
                "kind": "plot",
                "format": "png",
            }
        )
    artifacts.extend(
        {
            "path": str(path),
            "label": f"Manufactured electromechanics errors {path.stem}",
            "kind": "plot",
            "format": "png",
        }
        for path in error_plots
    )

    return artifacts
