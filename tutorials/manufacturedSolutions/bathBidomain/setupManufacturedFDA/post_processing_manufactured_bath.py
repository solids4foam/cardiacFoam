"""Post-processing helpers for bath-bidomain manufactured convergence tests."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import re
import sys

REPO_ROOT = Path(__file__).resolve().parents[4]
DRIVER_ROOT = REPO_ROOT / "applications" / "scripts" / "driverFoam"
if str(DRIVER_ROOT) not in sys.path:
    sys.path.insert(0, str(DRIVER_ROOT))

from openfoam_driver.postprocessing.style import (
    configure_matplotlib_defaults,
    finalize_matplotlib_figure,
    style_matplotlib_axes,
)

FIELD_NAMES = ("Vm", "phiE", "phiI", "u1", "u2", "u3")
RATE_FIELDS = (
    "Dimension",
    "Solver",
    "N_lower",
    "N_higher",
    *(f"rate_{name}" for name in FIELD_NAMES),
    *(f"rate_{norm}_{name}" for name in FIELD_NAMES for norm in ("L1", "L2", "Linf")),
)
ERROR_FIELDS = (
    "Dimension",
    "N",
    "Solver",
    *(f"{norm}_{name}" for name in FIELD_NAMES for norm in ("L1", "L2", "Linf")),
)
FILENAME_PATTERN = re.compile(r"bathBidomain_(\dD)_(\d+)_cells_(explicit|implicit)\.dat$")
BATH_ECG_SUMMARY_PATTERN = re.compile(
    r"BathECG_(?P<dimension>\dD)_(?P<cells>\d+)_cells_"
    r"(?P<solver>explicit|implicit)_DT[^_]+_manufacturedBathECGSummary\.dat$"
)
DISABLE_PLOT_ENV_VAR = "BATH_BIDOMAIN_DISABLE_PLOTS"
PLOT_DISABLED = os.environ.get(DISABLE_PLOT_ENV_VAR, "").strip().lower() in {
    "1",
    "true",
    "yes",
    "on",
}
plt = None

if not PLOT_DISABLED:
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        plt = None

SOLVER_MARKERS = {"explicit": "o", "implicit": "s"}
SOLVER_LINESTYLES = {"explicit": "-", "implicit": "--"}
FIELD_COLORS = {
    "Vm": "tab:blue",
    "phiE": "tab:orange",
    "phiI": "tab:green",
    "u1": "tab:red",
    "u2": "tab:purple",
    "u3": "tab:brown",
}
DIMENSION_COLORS = {"1D": "tab:blue", "2D": "tab:orange", "3D": "tab:green"}


def _write_csv(rows, destination: Path, fieldnames) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _has_matplotlib() -> bool:
    return not PLOT_DISABLED and plt is not None


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


def _finalize_axis_legend(axis) -> None:
    handles, labels = axis.get_legend_handles_labels()
    if handles:
        axis.legend(handles, labels)


def _positive_xy(x_values, y_values) -> tuple[list[float], list[float]]:
    pairs = [
        (float(x_value), float(y_value))
        for x_value, y_value in zip(x_values, y_values)
        if float(x_value) > 0 and float(y_value) > 0 and math.isfinite(float(y_value))
    ]
    if not pairs:
        return [], []
    x_positive, y_positive = zip(*pairs)
    return list(x_positive), list(y_positive)


def _safe_rate(e1: float, e2: float, h1: float, h2: float) -> float:
    if any(math.isnan(value) for value in (e1, e2)):
        return float("nan")
    if e1 <= 0 or e2 <= 0 or h1 == h2:
        return float("nan")
    return math.log(e1/e2)/math.log(h1/h2)


def _load_expected_filenames(output_dir: Path) -> set[str] | None:
    manifest_path = output_dir / "run_manifest.json"
    if not manifest_path.exists():
        return None

    manifest = json.loads(manifest_path.read_text())
    expected = set()
    for result in manifest.get("results", []):
        status = result.get("status")
        if status not in {"ok", "planned"}:
            continue
        params = result.get("params", {})
        dimension = params.get("dimension")
        cells = params.get("cells")
        solver = params.get("solver")
        if dimension is None or cells is None or solver is None:
            continue
        expected.add(f"bathBidomain_{dimension}_{int(cells)}_cells_{solver}.dat")
    return expected or None


def _load_expected_bath_ecg_summary_filenames(output_dir: Path) -> set[str] | None:
    manifest_path = output_dir / "run_manifest.json"
    if not manifest_path.exists():
        return None

    manifest = json.loads(manifest_path.read_text())
    expected = set()
    for result in manifest.get("results", []):
        status = result.get("status")
        if status not in {"ok", "planned"}:
            continue
        case_id = result.get("case_id")
        if not case_id:
            continue
        expected.add(f"BathECG_{case_id}_manufacturedBathECGSummary.dat")
    return expected or None


def read_error_dat_files(folder_name, expected_filenames: set[str] | None = None):
    folder = Path(folder_name)
    if not folder.exists():
        print("Folder does not exist:", folder)
        return []

    files = sorted(path for path in folder.glob("*.dat") if FILENAME_PATTERN.match(path.name))
    if expected_filenames is not None:
        files = [path for path in files if path.name in expected_filenames]

    rows = []
    for path in files:
        match = FILENAME_PATTERN.match(path.name)
        if match is None:
            continue

        content = path.read_text(encoding="utf-8", errors="ignore")
        row = {
            "Dimension": match.group(1),
            "N": int(match.group(2)),
            "Solver": match.group(3),
        }

        compact = content.replace("\n", " ")
        for field_name in FIELD_NAMES:
            field_match = re.search(rf"{field_name}\s+(\S+)\s+(\S+)\s+(\S+)", compact)
            row[f"L1_{field_name}"] = (
                float(field_match.group(1)) if field_match else float("nan")
            )
            row[f"L2_{field_name}"] = (
                float(field_match.group(2)) if field_match else float("nan")
            )
            row[f"Linf_{field_name}"] = (
                float(field_match.group(3)) if field_match else float("nan")
            )

        rows.append(row)

    return sorted(rows, key=lambda row: (row["Dimension"], row["Solver"], row["N"]))


def read_bath_ecg_summary_files(folder_name, expected_filenames: set[str] | None = None):
    folder = Path(folder_name)
    rows = []

    files = sorted(folder.glob("BathECG_*_manufacturedBathECGSummary.dat"))
    if expected_filenames is not None:
        files = [path for path in files if path.name in expected_filenames]

    for path in files:
        match = BATH_ECG_SUMMARY_PATTERN.match(path.name)
        if match is None:
            print("Skipping unrecognized bath ECG summary filename:", path.name)
            continue

        metadata = {}
        with path.open(encoding="utf-8", errors="ignore") as handle:
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("Manufactured bath ECG summary"):
                    continue
                if line.startswith("Electrode "):
                    break
                key, _, value = line.partition(" ")
                if key and value:
                    metadata[key] = value.strip()

        rows.append(
            {
                "Dimension": match.group("dimension"),
                "N": int(match.group("cells")),
                "Solver": match.group("solver"),
                "samples": int(metadata.get("samples", "0")),
                "field_L1": float(metadata.get("field_L1", "nan")),
                "field_L2": float(metadata.get("field_L2", "nan")),
                "field_Linf": float(metadata.get("field_Linf", "nan")),
            }
        )

    return sorted(rows, key=lambda row: (row["Dimension"], row["Solver"], row["N"]))


def compute_convergence_rates(rows):
    grouped = {}
    for row in rows:
        grouped.setdefault((row["Dimension"], row["Solver"]), []).append(row)

    rate_rows = []
    for (dimension, solver), group_rows in sorted(grouped.items()):
        ordered = sorted(group_rows, key=lambda row: row["N"])
        for lower, higher in zip(ordered, ordered[1:]):
            n1 = int(lower["N"])
            n2 = int(higher["N"])
            if n1 == n2:
                continue

            h1 = 1.0/n1
            h2 = 1.0/n2
            row = {
                "Dimension": dimension,
                "Solver": solver,
                "N_lower": n1,
                "N_higher": n2,
            }
            for field_name in FIELD_NAMES:
                for norm_name in ("L1", "L2", "Linf"):
                    row[f"rate_{norm_name}_{field_name}"] = _safe_rate(
                        lower[f"{norm_name}_{field_name}"],
                        higher[f"{norm_name}_{field_name}"],
                        h1,
                        h2,
                    )
                row[f"rate_{field_name}"] = row[f"rate_Linf_{field_name}"]
            rate_rows.append(row)

    return rate_rows


def compute_bath_ecg_convergence_rates(rows):
    grouped = {}
    for row in rows:
        grouped.setdefault((row["Dimension"], row["Solver"]), []).append(row)

    rate_rows = []
    for (dimension, solver), group_rows in sorted(grouped.items()):
        ordered = sorted(group_rows, key=lambda row: row["N"])
        for lower, higher in zip(ordered, ordered[1:]):
            n1 = int(lower["N"])
            n2 = int(higher["N"])
            if n1 == n2:
                continue
            h1 = 1.0/n1
            h2 = 1.0/n2
            rate_rows.append(
                {
                    "Dimension": dimension,
                    "Solver": solver,
                    "N_lower": n1,
                    "N_higher": n2,
                    "rate_field_L1": _safe_rate(lower["field_L1"], higher["field_L1"], h1, h2),
                    "rate_field_L2": _safe_rate(lower["field_L2"], higher["field_L2"], h1, h2),
                    "rate_field_Linf": _safe_rate(
                        lower["field_Linf"],
                        higher["field_Linf"],
                        h1,
                        h2,
                    ),
                }
            )

    return rate_rows


def _plot_dimension_errors_on_axis(axis, rows, dimension: str) -> bool:
    dimension_rows = _filter_rows(rows, Dimension=dimension)
    if not dimension_rows:
        axis.text(
            0.5,
            0.5,
            f"No {dimension} data",
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color="0.4",
        )
        style_matplotlib_axes(
            axis,
            title=f"Myocardium subdomain Linf errors ({dimension})",
            xlabel="Number of cells (N)",
            ylabel="Linf error",
            legend=False,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.6},
        )
        return False

    plotted = False
    for solver in _unique_values(dimension_rows, "Solver"):
        solver_rows = _filter_rows(dimension_rows, Solver=solver)
        for field_name in FIELD_NAMES:
            ns, errors = _positive_xy(
                [row["N"] for row in solver_rows],
                [row[f"Linf_{field_name}"] for row in solver_rows],
            )
            if not ns:
                continue
            axis.loglog(
                ns,
                errors,
                marker=SOLVER_MARKERS.get(solver, "o"),
                linestyle=SOLVER_LINESTYLES.get(solver, "-"),
                color=FIELD_COLORS.get(field_name),
                label=f"{field_name} myocardium ({solver}, {dimension})",
            )
            plotted = True

    if not plotted:
        axis.text(
            0.5,
            0.5,
            f"No positive {dimension} error data",
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color="0.4",
        )

    style_matplotlib_axes(
        axis,
        title=f"Myocardium subdomain Linf errors ({dimension})",
        xlabel="Number of cells (N)",
        ylabel="Linf error",
        legend=False,
        grid_kwargs={"which": "both", "ls": "--", "alpha": 0.6},
    )
    _finalize_axis_legend(axis)
    return True


def _plot_vm_across_dimensions_on_axis(axis, rows) -> bool:
    plotted = False
    for dimension in _unique_values(rows, "Dimension"):
        dimension_rows = _filter_rows(rows, Dimension=dimension)
        for solver in _unique_values(dimension_rows, "Solver"):
            solver_rows = _filter_rows(dimension_rows, Solver=solver)
            if not solver_rows:
                continue
            ns, errors = _positive_xy(
                [row["N"] for row in solver_rows],
                [row["Linf_Vm"] for row in solver_rows],
            )
            if not ns:
                continue
            axis.loglog(
                ns,
                errors,
                marker=SOLVER_MARKERS.get(solver, "o"),
                linestyle=SOLVER_LINESTYLES.get(solver, "--"),
                color=DIMENSION_COLORS.get(dimension, "black"),
                label=f"{dimension} ({solver})",
            )
            plotted = True

    if not plotted:
        axis.text(
            0.5,
            0.5,
            "No Vm data",
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color="0.4",
        )

    style_matplotlib_axes(
        axis,
        title="Myocardium subdomain Linf error of Vm across dimensions",
        xlabel="Number of cells (N)",
        ylabel="Linf error (Vm)",
        legend=False,
        grid_kwargs={"which": "both", "ls": "--", "alpha": 0.6},
    )
    _finalize_axis_legend(axis)
    return plotted


def _plot_field_errors(rows, destination: Path) -> Path | None:
    if not _has_matplotlib():
        return None

    configure_matplotlib_defaults()
    grouped = {}
    for row in rows:
        grouped.setdefault((row["Dimension"], row["Solver"]), []).append(row)

    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    fig.suptitle("Myocardium subdomain manufactured errors", fontsize=13)
    flat_axes = axes.ravel()
    for axis, field_name in zip(flat_axes, FIELD_NAMES):
        plotted = False
        for (dimension, solver), group_rows in sorted(grouped.items()):
            ordered = sorted(group_rows, key=lambda item: item["N"])
            ns, errors = _positive_xy(
                [int(item["N"]) for item in ordered],
                [item[f"Linf_{field_name}"] for item in ordered],
            )
            if not ns:
                continue
            axis.loglog(ns, errors, marker="o", label=f"{dimension} {solver}")
            plotted = True
        if not plotted:
            axis.text(
                0.5,
                0.5,
                "No positive error data",
                transform=axis.transAxes,
                ha="center",
                va="center",
                fontsize=10,
                color="0.4",
            )
        axis.set_title(field_name)
        axis.set_xlabel("N")
        axis.set_ylabel("Linf error")
        axis.grid(True, which="both", alpha=0.25)

    handles, labels = flat_axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncols=min(3, len(handles)))
    destination.parent.mkdir(parents=True, exist_ok=True)
    finalize_matplotlib_figure(fig, save_path=destination, show=False, close=True)
    return destination


def _plot_convergence_rates(rows, destination: Path) -> Path | None:
    if not _has_matplotlib() or not rows:
        return None

    configure_matplotlib_defaults()
    labels = [f"{row['Dimension']} {row['Solver']} {row['N_lower']}-{row['N_higher']}" for row in rows]
    x_positions = range(len(labels))
    fig, axis = plt.subplots(figsize=(max(8, len(labels) * 1.4), 5))
    width = 0.12
    offsets = [width * (index - (len(FIELD_NAMES) - 1) / 2) for index in range(len(FIELD_NAMES))]
    for offset, field_name in zip(offsets, FIELD_NAMES):
        values = [row.get(f"rate_{field_name}", float("nan")) for row in rows]
        axis.bar([x + offset for x in x_positions], values, width=width, label=field_name)
    axis.axhline(1.0, color="0.35", linewidth=0.8, linestyle="--")
    axis.axhline(2.0, color="0.35", linewidth=0.8, linestyle=":")
    axis.set_xticks(list(x_positions))
    axis.set_xticklabels(labels, rotation=30, ha="right")
    axis.set_ylabel("Observed rate")
    axis.set_title("Myocardium subdomain observed convergence rates")
    axis.grid(True, axis="y", alpha=0.25)
    axis.legend(ncols=3)
    destination.parent.mkdir(parents=True, exist_ok=True)
    finalize_matplotlib_figure(fig, save_path=destination, show=False, close=True)
    return destination


def _plot_bath_ecg_errors(rows, destination: Path) -> Path | None:
    if not _has_matplotlib() or not rows:
        return None

    configure_matplotlib_defaults()
    grouped = {}
    for row in rows:
        grouped.setdefault((row["Dimension"], row["Solver"]), []).append(row)

    fig, axis = plt.subplots(figsize=(8, 5))
    for (dimension, solver), group_rows in sorted(grouped.items()):
        ordered = sorted(group_rows, key=lambda item: item["N"])
        ns, errors = _positive_xy(
            [int(item["N"]) for item in ordered],
            [item["field_Linf"] for item in ordered],
        )
        if not ns:
            continue
        axis.loglog(
            ns,
            errors,
            marker="o",
            label=f"{dimension} {solver}",
        )
    axis.set_xlabel("N")
    axis.set_ylabel("Bath-domain phiE Linf error")
    axis.set_title("Bath-domain phiE manufactured error")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend()
    destination.parent.mkdir(parents=True, exist_ok=True)
    finalize_matplotlib_figure(fig, save_path=destination, show=False, close=True)
    return destination


def plot_vm_across_dimensions(
    rows,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_matplotlib():
        print("matplotlib is not available; skipping bath-bidomain Vm plot.")
        return None

    configure_matplotlib_defaults()
    fig, axis = plt.subplots(figsize=(8, 6))
    _plot_vm_across_dimensions_on_axis(axis, rows)
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def plot_errors_by_dimension(
    rows,
    *,
    save_dir: str | Path | None = None,
    show: bool = True,
) -> list[Path]:
    if not _has_matplotlib():
        print("matplotlib is not available; skipping bath-bidomain dimension plots.")
        return []

    configure_matplotlib_defaults()
    output_paths: list[Path] = []
    for dimension in _unique_values(rows, "Dimension"):
        fig, axis = plt.subplots(figsize=(8, 6))
        _plot_dimension_errors_on_axis(axis, rows, dimension)
        save_path = None
        if save_dir is not None:
            save_path = Path(save_dir) / f"bath_bidomain_errors_{dimension.lower()}_implicit_explicit.png"
            output_paths.append(save_path)
        finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return output_paths


def plot_summary_dashboard(
    rows,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_matplotlib():
        print("matplotlib is not available; skipping bath-bidomain summary dashboard.")
        return None

    configure_matplotlib_defaults()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Bath-bidomain manufactured-solution convergence summary", fontsize=14)

    _plot_vm_across_dimensions_on_axis(axes[0, 0], rows)
    _plot_dimension_errors_on_axis(axes[0, 1], rows, "1D")
    _plot_dimension_errors_on_axis(axes[1, 0], rows, "2D")
    _plot_dimension_errors_on_axis(axes[1, 1], rows, "3D")

    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **_: object) -> list[dict]:
    del setup_root
    output_path = Path(output_dir)
    expected_filenames = _load_expected_filenames(output_path)
    error_rows = read_error_dat_files(output_path, expected_filenames=expected_filenames)
    if not error_rows:
        print(f"No bath-bidomain .dat files found to post-process in: {output_dir}")
        return []

    errors_csv = output_path / "bath_bidomain_errors.csv"
    _write_csv(error_rows, errors_csv, ERROR_FIELDS)

    convergence_rates = compute_convergence_rates(error_rows)
    rates_csv = output_path / "bath_bidomain_convergence_rates.csv"
    _write_csv(convergence_rates, rates_csv, RATE_FIELDS)

    artifacts = [
        {
            "path": str(errors_csv),
            "label": "Myocardium subdomain Linf error table",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(rates_csv),
            "label": "Myocardium subdomain convergence rates table",
            "kind": "table",
            "format": "csv",
        },
    ]

    if _has_matplotlib():
        field_plot = _plot_field_errors(
            error_rows,
            output_path / "bath_bidomain_field_errors.png",
        )
        vm_plot = plot_vm_across_dimensions(
            error_rows,
            save_path=output_path / "bath_bidomain_vm_across_dimensions.png",
            show=False,
        )
        summary_plot = plot_summary_dashboard(
            error_rows,
            save_path=output_path / "bath_bidomain_summary_dashboard.png",
            show=False,
        )
        rate_plot = _plot_convergence_rates(
            convergence_rates,
            output_path / "bath_bidomain_convergence_rates.png",
        )
        dimension_plots = plot_errors_by_dimension(
            error_rows,
            save_dir=output_path,
            show=False,
        )
        for path, label in (
            (field_plot, "Myocardium subdomain field errors"),
            (vm_plot, "Myocardium subdomain Vm across dimensions"),
            (summary_plot, "Myocardium subdomain summary dashboard"),
            (rate_plot, "Myocardium subdomain convergence rates"),
        ):
            if path is None:
                continue
            artifacts.append(
                {
                    "path": str(path),
                    "label": label,
                    "kind": "plot",
                    "format": "png",
                }
            )
        artifacts.extend(
            {
                "path": str(path),
                "label": f"Myocardium subdomain errors {path.stem}",
                "kind": "plot",
                "format": "png",
            }
            for path in dimension_plots
        )
    else:
        if not PLOT_DISABLED:
            print("matplotlib is not available; generated bath-bidomain CSV files only.")
        else:
            print(
                "Bath-bidomain plots disabled; unset "
                f"{DISABLE_PLOT_ENV_VAR} to export PNG plots."
            )

    bath_ecg_rows = read_bath_ecg_summary_files(
        output_path,
        expected_filenames=_load_expected_bath_ecg_summary_filenames(output_path),
    )
    if bath_ecg_rows:
        bath_ecg_csv = output_path / "bath_ecg_errors.csv"
        _write_csv(
            bath_ecg_rows,
            bath_ecg_csv,
            ("Dimension", "N", "Solver", "samples", "field_L1", "field_L2", "field_Linf"),
        )
        bath_ecg_rate_rows = compute_bath_ecg_convergence_rates(bath_ecg_rows)
        bath_ecg_rates_csv = output_path / "bath_ecg_convergence_rates.csv"
        _write_csv(
            bath_ecg_rate_rows,
            bath_ecg_rates_csv,
            (
                "Dimension",
                "Solver",
                "N_lower",
                "N_higher",
                "rate_field_L1",
                "rate_field_L2",
                "rate_field_Linf",
            ),
        )
        artifacts.extend(
            [
                {
                    "path": str(bath_ecg_csv),
                    "label": "Bath-domain phiE error table",
                    "kind": "table",
                    "format": "csv",
                },
                {
                    "path": str(bath_ecg_rates_csv),
                    "label": "Bath-domain phiE convergence rates table",
                    "kind": "table",
                    "format": "csv",
                },
            ]
        )
        bath_ecg_plot = _plot_bath_ecg_errors(
            bath_ecg_rows,
            output_path / "bath_ecg_phiE_errors.png",
        )
        if bath_ecg_plot is not None:
            artifacts.append(
                {
                    "path": str(bath_ecg_plot),
                    "label": "Bath-domain phiE errors",
                    "kind": "plot",
                    "format": "png",
                }
            )

    print("\nBath-bidomain convergence rates:")
    for row in convergence_rates:
        print(row)

    return artifacts
