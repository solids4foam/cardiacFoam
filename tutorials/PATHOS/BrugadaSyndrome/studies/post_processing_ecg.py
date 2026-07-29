"""Post-processing helpers for ECG and manufactured pseudo-ECG outputs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import re
import sys

try:
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

REF_PATTERN = re.compile(r"refQ(?P<q>\d+)_(?P<electrode>.+)")
ERR_PATTERN = re.compile(r"errQ(?P<q>\d+)_(?P<electrode>.+)")


def _has_matplotlib() -> bool:
    return plt is not None


def _discover_file(output_dir: Path, filename: str) -> Path | None:
    candidates = [output_dir / filename]
    case_root = output_dir.parent
    candidates.extend(sorted(case_root.glob(f"processor*/{output_dir.name}/{filename}")))

    existing = [candidate for candidate in candidates if candidate.exists()]
    if not existing:
        return None

    non_empty = [candidate for candidate in existing if candidate.stat().st_size > 0]
    if non_empty:
        return non_empty[0]
    return existing[0]


def _parse_timeseries(path: Path) -> tuple[list[str], list[list[float]]]:
    header: list[str] | None = None
    rows: list[list[float]] = []

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("#"):
                tokens = line.lstrip("#").strip().split()
                if tokens:
                    header = tokens
                continue

            tokens = line.split()
            if header is None:
                header = ["time"] + [f"signal_{idx}" for idx in range(1, len(tokens))]

            if len(tokens) != len(header):
                continue

            try:
                rows.append([float(token) for token in tokens])
            except ValueError:
                continue

    if header is None:
        raise ValueError(f"No header found in {path}")
    if not rows:
        raise ValueError(f"No numeric rows found in {path}")
    return header, rows


def _columns_from_rows(header: list[str], rows: list[list[float]]) -> dict[str, list[float]]:
    columns = {name: [] for name in header}
    for row in rows:
        for name, value in zip(header, row):
            columns[name].append(value)
    return columns


def _group_manufactured_columns(columns: dict[str, list[float]]) -> dict[str, dict[str, object]]:
    electrodes = sorted(name.removeprefix("numeric_") for name in columns if name.startswith("numeric_"))
    grouped: dict[str, dict[str, object]] = {}

    for electrode in electrodes:
        refs: list[tuple[int, str, list[float]]] = []
        errs: list[tuple[int, str, list[float]]] = []

        for name, values in columns.items():
            ref_match = REF_PATTERN.fullmatch(name)
            if ref_match and ref_match.group("electrode") == electrode:
                refs.append((int(ref_match.group("q")), name, values))
                continue

            err_match = ERR_PATTERN.fullmatch(name)
            if err_match and err_match.group("electrode") == electrode:
                errs.append((int(err_match.group("q")), name, values))

        grouped[electrode] = {
            "numeric": columns.get(f"numeric_{electrode}", []),
            "refs": sorted(refs),
            "errs": sorted(errs),
            "delta": columns.get(f"deltaQuadrature_{electrode}", []),
        }

    return grouped


def _parse_summary(path: Path) -> tuple[dict[str, str], list[dict[str, float | str]]]:
    metadata: dict[str, str] = {}
    rows: list[dict[str, float | str]] = []
    table_columns: list[str] | None = None

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("Manufactured pseudo-ECG summary"):
                continue
            if line.startswith("Electrode "):
                table_columns = line.split()
                continue

            if table_columns is None:
                key, _, value = line.partition(" ")
                if key and value:
                    metadata[key] = value.strip()
                continue

            tokens = line.split()
            if len(tokens) != len(table_columns):
                continue

            row: dict[str, float | str] = {"Electrode": tokens[0]}
            for column, token in zip(table_columns[1:], tokens[1:]):
                row[column] = float(token)
            rows.append(row)

    return metadata, rows


def _subplot_shape(count: int, max_cols: int = 3) -> tuple[int, int]:
    if count <= 0:
        return 1, 1
    cols = min(max_cols, count)
    rows = int(math.ceil(count / cols))
    return rows, cols


def _common_y_limits(series_list: list[list[float]]) -> tuple[float, float] | None:
    values = [value for series in series_list for value in series]
    if not values:
        return None
    ymin = min(values)
    ymax = max(values)
    if math.isclose(ymin, ymax, rel_tol=0.0, abs_tol=1e-14):
        pad = max(1e-9, abs(ymin) * 0.05 + 1e-9)
        return ymin - pad, ymax + pad
    pad = 0.05 * (ymax - ymin)
    return ymin - pad, ymax + pad


def _positive_floor(series_list: list[list[float]]) -> float:
    positives = [value for series in series_list for value in series if value > 0]
    if not positives:
        return 1e-16
    return min(positives) * 0.5


def _plot_trace_grid(
    times: list[float],
    grouped_series: dict[str, list[tuple[str, list[float]]]],
    *,
    title: str,
    ylabel: str,
    save_path: Path,
    log_y: bool = False,
) -> Path | None:
    if not grouped_series or not _has_matplotlib():
        return None

    configure_matplotlib_defaults()
    names = list(grouped_series.keys())
    nrows, ncols = _subplot_shape(len(names))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(4.2 * ncols, 3.1 * nrows),
        sharex=True,
        squeeze=False,
    )
    flat_axes = axes.flatten()

    if log_y:
        floor = _positive_floor(
            [series for series_group in grouped_series.values() for _, series in series_group]
        )
    else:
        y_limits = _common_y_limits(
            [series for series_group in grouped_series.values() for _, series in series_group]
        )

    for idx, electrode in enumerate(names):
        ax = flat_axes[idx]
        for label, series in grouped_series[electrode]:
            values = [max(value, floor) for value in series] if log_y else series
            if log_y:
                ax.semilogy(times, values, linewidth=1.2, label=label)
            else:
                ax.plot(times, values, linewidth=1.2, label=label)
        if not log_y and y_limits is not None:
            ax.set_ylim(*y_limits)
        style_matplotlib_axes(
            ax,
            title=electrode,
            xlabel="time (s)",
            ylabel=ylabel,
            legend=True,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
        )
        ax.tick_params(labelsize=8)

    for idx in range(len(names), len(flat_axes)):
        flat_axes[idx].axis("off")

    fig.suptitle(title, fontsize=12)
    finalize_matplotlib_figure(fig, save_path=save_path, show=False, close=True)
    return save_path


def _plot_summary_bars(
    metadata: dict[str, str],
    rows: list[dict[str, float | str]],
    *,
    save_path: Path,
) -> Path | None:
    if not rows or not _has_matplotlib():
        return None

    configure_matplotlib_defaults()
    electrodes = [str(row["Electrode"]) for row in rows]
    x = list(range(len(electrodes)))

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

    ref_vals = [float(row["Linf_err_ref"]) for row in rows]
    check_vals = [float(row["Linf_err_check"]) for row in rows]
    delta_vals = [float(row["Linf_delta_ref"]) for row in rows]

    width = 0.38
    axes[0].bar([xi - width / 2 for xi in x], ref_vals, width=width, label="Linf err ref")
    axes[0].bar([xi + width / 2 for xi in x], check_vals, width=width, label="Linf err check")
    axes[0].set_yscale("log")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(electrodes, rotation=30, ha="right")
    style_matplotlib_axes(
        axes[0],
        title="Per-electrode Linf errors",
        xlabel="Electrode",
        ylabel="Linf error",
        legend=True,
        grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
    )

    axes[1].bar(x, delta_vals, color="tab:purple", label="Linf quadrature delta")
    axes[1].set_yscale("log")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(electrodes, rotation=30, ha="right")
    style_matplotlib_axes(
        axes[1],
        title="Per-electrode Linf quadrature delta",
        xlabel="Electrode",
        ylabel="Linf delta",
        legend=True,
        grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
    )

    dim_txt = metadata.get("dimension", "?")
    q_check = metadata.get("qCheck", "?")
    q_ref = metadata.get("qReference", "?")
    fig.suptitle(
        f"Manufactured pseudo-ECG summary ({dim_txt}, qCheck={q_check}, qRef={q_ref})",
        fontsize=12,
    )
    finalize_matplotlib_figure(fig, save_path=save_path, show=False, close=True)
    return save_path


def _plot_dashboard(
    raw_columns: dict[str, list[float]],
    manufactured_groups: dict[str, dict[str, object]] | None,
    summary_rows: list[dict[str, float | str]] | None,
    *,
    save_path: Path,
) -> Path | None:
    if not _has_matplotlib():
        return None

    configure_matplotlib_defaults()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    times = raw_columns["time"]
    for electrode, series in raw_columns.items():
        if electrode == "time":
            continue
        axes[0, 0].plot(times, series, linewidth=1.2, label=electrode)
    style_matplotlib_axes(
        axes[0, 0],
        title="pseudoECG traces",
        xlabel="time (s)",
        ylabel="pseudoECG",
        legend=True,
        grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
    )

    if manufactured_groups:
        error_floor = _positive_floor(
            [
                series
                for group in manufactured_groups.values()
                for _, _, series in group["errs"]  # type: ignore[index]
            ]
        )
        manufactured_time = next(iter(manufactured_groups.values()))["numeric"]  # type: ignore[index]
        del manufactured_time
        manufactured_times = raw_columns["time"]
        for electrode, group in manufactured_groups.items():
            ref_errs = group["errs"]  # type: ignore[assignment]
            if not ref_errs:
                continue
            highest_q = ref_errs[-1][2]
            axes[0, 1].semilogy(
                manufactured_times,
                [max(value, error_floor) for value in highest_q],
                linewidth=1.2,
                label=electrode,
            )
        style_matplotlib_axes(
            axes[0, 1],
            title="Reference error vs time",
            xlabel="time (s)",
            ylabel="abs. error",
            legend=True,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
        )
    else:
        axes[0, 1].text(0.5, 0.5, "No manufactured data", ha="center", va="center")
        style_matplotlib_axes(
            axes[0, 1],
            title="Reference error vs time",
            xlabel="time (s)",
            ylabel="abs. error",
            legend=False,
        )

    if summary_rows:
        electrodes = [str(row["Electrode"]) for row in summary_rows]
        x = list(range(len(electrodes)))
        ref_vals = [float(row["Linf_err_ref"]) for row in summary_rows]
        delta_vals = [float(row["Linf_delta_ref"]) for row in summary_rows]

        axes[1, 0].bar(x, ref_vals, color="tab:red")
        axes[1, 0].set_yscale("log")
        axes[1, 0].set_xticks(x)
        axes[1, 0].set_xticklabels(electrodes, rotation=30, ha="right")
        style_matplotlib_axes(
            axes[1, 0],
            title="Linf reference error",
            xlabel="Electrode",
            ylabel="Linf error",
            legend=False,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
        )

        axes[1, 1].bar(x, delta_vals, color="tab:purple")
        axes[1, 1].set_yscale("log")
        axes[1, 1].set_xticks(x)
        axes[1, 1].set_xticklabels(electrodes, rotation=30, ha="right")
        style_matplotlib_axes(
            axes[1, 1],
            title="Linf quadrature delta",
            xlabel="Electrode",
            ylabel="Linf delta",
            legend=False,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
        )
    else:
        axes[1, 0].text(0.5, 0.5, "No summary file", ha="center", va="center")
        axes[1, 1].text(0.5, 0.5, "No summary file", ha="center", va="center")
        style_matplotlib_axes(
            axes[1, 0],
            title="Linf reference error",
            xlabel="Electrode",
            ylabel="Linf error",
            legend=False,
        )
        style_matplotlib_axes(
            axes[1, 1],
            title="Linf quadrature delta",
            xlabel="Electrode",
            ylabel="Linf delta",
            legend=False,
        )

    fig.suptitle("ECG verification dashboard", fontsize=13)
    finalize_matplotlib_figure(fig, save_path=save_path, show=False, close=True)
    return save_path


def _write_summary_csv(rows: list[dict[str, float | str]], destination: Path) -> Path:
    fieldnames = [
        "Electrode",
        "L1_err_ref",
        "L2_err_ref",
        "Linf_err_ref",
        "L1_err_check",
        "L2_err_check",
        "Linf_err_check",
        "L1_delta_ref",
        "L2_delta_ref",
        "Linf_delta_ref",
    ]
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    return destination


def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **_: object):
    del setup_root
    output_path = Path(output_dir)
    artifacts: list[dict[str, str]] = []

    pseudo_path = _discover_file(output_path, "pseudoECG.dat")
    if pseudo_path is None:
        print(f"pseudoECG.dat not found in: {output_dir}")
        return artifacts

    raw_header, raw_rows = _parse_timeseries(pseudo_path)
    raw_columns = _columns_from_rows(raw_header, raw_rows)
    raw_grouped = {
        name: [(name, values)]
        for name, values in raw_columns.items()
        if name != "time"
    }

    raw_plot = _plot_trace_grid(
        raw_columns["time"],
        raw_grouped,
        title="pseudoECG traces by electrode",
        ylabel="pseudoECG",
        save_path=output_path / "pseudoECG_traces_grid.png",
    )
    if raw_plot is not None:
        artifacts.append(
            {
                "path": str(raw_plot),
                "label": "pseudoECG traces grid",
                "kind": "plot",
                "format": "png",
            }
        )

    manufactured_groups: dict[str, dict[str, object]] | None = None
    manufactured_summary_rows: list[dict[str, float | str]] | None = None

    manufactured_path = _discover_file(output_path, "manufacturedPseudoECG.dat")
    if manufactured_path is not None:
        manufactured_header, manufactured_rows = _parse_timeseries(manufactured_path)
        manufactured_columns = _columns_from_rows(manufactured_header, manufactured_rows)
        manufactured_groups = _group_manufactured_columns(manufactured_columns)
        manufactured_times = manufactured_columns["time"]

        overlay_series = {
            electrode: [
                ("numeric", group["numeric"]),  # type: ignore[list-item]
                *[
                    (name, values)
                    for _, name, values in group["refs"]  # type: ignore[index]
                ],
            ]
            for electrode, group in manufactured_groups.items()
        }
        overlay_plot = _plot_trace_grid(
            manufactured_times,
            overlay_series,
            title="Manufactured pseudoECG: numeric vs reference",
            ylabel="pseudoECG",
            save_path=output_path / "manufacturedPseudoECG_reference_overlay.png",
        )
        if overlay_plot is not None:
            artifacts.append(
                {
                    "path": str(overlay_plot),
                    "label": "Manufactured pseudoECG numeric vs reference",
                    "kind": "plot",
                    "format": "png",
                }
            )

        error_series = {
            electrode: [
                *[
                    (name, values)
                    for _, name, values in group["errs"]  # type: ignore[index]
                ],
                ("deltaQuadrature", group["delta"]),  # type: ignore[list-item]
            ]
            for electrode, group in manufactured_groups.items()
        }
        error_plot = _plot_trace_grid(
            manufactured_times,
            error_series,
            title="Manufactured pseudoECG: error diagnostics",
            ylabel="abs. error",
            save_path=output_path / "manufacturedPseudoECG_error_diagnostics.png",
            log_y=True,
        )
        if error_plot is not None:
            artifacts.append(
                {
                    "path": str(error_plot),
                    "label": "Manufactured pseudoECG error diagnostics",
                    "kind": "plot",
                    "format": "png",
                }
            )

    summary_path = _discover_file(output_path, "manufacturedPseudoECGSummary.dat")
    if summary_path is not None:
        summary_metadata, manufactured_summary_rows = _parse_summary(summary_path)
        summary_csv = _write_summary_csv(
            manufactured_summary_rows,
            output_path / "manufacturedPseudoECGSummary.csv",
        )
        artifacts.append(
            {
                "path": str(summary_csv),
                "label": "Manufactured pseudoECG summary table",
                "kind": "table",
                "format": "csv",
            }
        )

        summary_plot = _plot_summary_bars(
            summary_metadata,
            manufactured_summary_rows,
            save_path=output_path / "manufacturedPseudoECG_summary_bars.png",
        )
        if summary_plot is not None:
            artifacts.append(
                {
                    "path": str(summary_plot),
                    "label": "Manufactured pseudoECG summary bars",
                    "kind": "plot",
                    "format": "png",
                }
            )

    dashboard_plot = _plot_dashboard(
        raw_columns,
        manufactured_groups,
        manufactured_summary_rows,
        save_path=output_path / "ecg_verification_dashboard.png",
    )
    if dashboard_plot is not None:
        artifacts.append(
            {
                "path": str(dashboard_plot),
                "label": "ECG verification dashboard",
                "kind": "plot",
                "format": "png",
            }
        )

    return artifacts


def main() -> int:
    parser = argparse.ArgumentParser(description="Post-process ECG and manufactured pseudo-ECG outputs.")
    parser.add_argument(
        "--output-dir",
        default="postProcessing",
        help="Directory containing pseudoECG.dat and optional manufactured outputs.",
    )
    args = parser.parse_args()

    artifacts = run_postprocessing(output_dir=args.output_dir)
    print(f"Generated {len(artifacts)} artifact(s).")
    for artifact in artifacts:
        print(f"- {artifact['label']}: {artifact['path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
