"""Post-processing helpers for manufactured-solution convergence tests."""

from __future__ import annotations

import csv
import bisect
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None
try:
    from matplotlib.patches import Patch
except ModuleNotFoundError:
    Patch = None
try:
    from matplotlib.colors import LogNorm
except ModuleNotFoundError:
    LogNorm = None

try:
    import numpy as np
except ModuleNotFoundError:
    np = None

TUTORIALS_ROOT = Path(__file__).resolve().parents[2]
if str(TUTORIALS_ROOT) not in sys.path:
    sys.path.insert(0, str(TUTORIALS_ROOT))

from openfoam_driver.postprocessing.style import (
    configure_matplotlib_defaults,
    finalize_matplotlib_figure,
    style_matplotlib_axes,
)
from openfoam_driver.core.defaults import manufactured_fda as driver_defaults

RATE_FIELDS = ("Dimension", "Solver", "N_lower", "N_higher", "rate_Vm", "rate_u1", "rate_u2")
FILENAME_PATTERN = re.compile(r"(\dD)_(\d+)_cells_(explicit|implicit)")
ECG_SUMMARY_PATTERN = re.compile(
    r"ECG_(?P<dimension>\dD)_(?P<cells>\d+)_cells_(?P<solver>explicit|implicit)_DT[^_]+_"
    r"manufacturedPseudoECGSummary\.dat$"
)
ECG_TIMESERIES_PATTERN = re.compile(
    r"ECG_(?P<dimension>\dD)_(?P<cells>\d+)_cells_(?P<solver>explicit|implicit)_DT[^_]+_"
    r"manufacturedPseudoECG\.dat$"
)
REF_PATTERN = re.compile(r"refQ(?P<q>\d+)_(?P<electrode>.+)")
ERR_PATTERN = re.compile(r"errQ(?P<q>\d+)_(?P<electrode>.+)")
DELTA_PATTERN = re.compile(
    r"deltaQuadratureQ(?P<qcheck>\d+)_Q(?P<qreference>\d+)_(?P<electrode>.+)"
)
SOLVER_MARKERS = {"explicit": "o", "implicit": "s"}
SOLVER_LINESTYLES = {"explicit": "-", "implicit": "--"}
FIELD_COLORS = {"Linf_V": "tab:blue", "Linf_u1": "tab:orange", "Linf_u2": "tab:green"}
FIELD_LABELS = {"Linf_V": "Vm", "Linf_u1": "u1", "Linf_u2": "u2"}
DIMENSION_COLORS = {"1D": "tab:blue", "2D": "tab:orange", "3D": "tab:green"}
SUPPORTED_ECG_POSTPROCESS_DIMENSIONS = ("3D",)
ECG_RATE_FIELDS = (
    "Dimension",
    "Solver",
    "N_lower",
    "N_higher",
    "rate_max_Linf_err_ref",
    "rate_mean_Linf_err_ref",
    "rate_max_Linf_delta_ref",
    "rate_mean_Linf_delta_ref",
)
ECG_SUMMARY_FIELDS = (
    "Dimension",
    "N",
    "Solver",
    "samples",
    "electrodes",
    "qCheck",
    "qChecks",
    "qReference",
    "max_Linf_err_ref",
    "mean_Linf_err_ref",
    "max_Linf_err_check",
    "mean_Linf_err_check",
    "max_Linf_delta_ref",
    "mean_Linf_delta_ref",
)
SUMMARY_ERR_Q_PATTERN = re.compile(r"Linf_err_q(?P<q>\d+)$")
SUMMARY_DELTA_Q_PATTERN = re.compile(r"Linf_delta_q(?P<q>\d+)_ref$")


def _has_matplotlib() -> bool:
    return plt is not None


def _has_surface_plotting() -> bool:
    return plt is not None and np is not None


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
        expected.add(f"{dimension}_{int(cells)}_cells_{solver}.dat")
    return expected or None


def _parse_ecg_summary_file(path: Path) -> tuple[dict[str, str], list[dict[str, float | str]]]:
    metadata: dict[str, str] = {}
    rows: list[dict[str, float | str]] = []
    table_columns: list[str] | None = None

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("Manufactured pseudo-ECG summary"):
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
        deltas: list[tuple[int, int, str, list[float]]] = []

        for name, values in columns.items():
            ref_match = REF_PATTERN.fullmatch(name)
            if ref_match and ref_match.group("electrode") == electrode:
                refs.append((int(ref_match.group("q")), name, values))
                continue

            err_match = ERR_PATTERN.fullmatch(name)
            if err_match and err_match.group("electrode") == electrode:
                errs.append((int(err_match.group("q")), name, values))
                continue

            delta_match = DELTA_PATTERN.fullmatch(name)
            if delta_match and delta_match.group("electrode") == electrode:
                deltas.append(
                    (
                        int(delta_match.group("qcheck")),
                        int(delta_match.group("qreference")),
                        name,
                        values,
                    )
                )

        grouped[electrode] = {
            "numeric": columns.get(f"numeric_{electrode}", []),
            "refs": sorted(refs),
            "errs": sorted(errs),
            "delta": deltas[0][3]
            if deltas
            else columns.get(f"deltaQuadrature_{electrode}", []),
            "deltas": sorted(deltas),
        }

    return grouped


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


def _reference_label(q: int) -> str:
    return f"reference(q={q})"


def _reference_error_label(q: int) -> str:
    return f"abs error to reference(q={q})"


def _quadrature_orders_from_group(group: dict[str, object]) -> tuple[int | None, int | None]:
    refs = group.get("refs", [])
    q_values = sorted({int(q) for q, _, _ in refs})
    if not q_values:
        errs = group.get("errs", [])
        q_values = sorted({int(q) for q, _, _ in errs})
    if not q_values:
        return None, None
    return q_values[0], q_values[-1]


def _check_quadrature_orders_from_group(group: dict[str, object]) -> list[int]:
    refs = sorted({int(q) for q, _, _ in group.get("refs", [])})
    deltas = sorted({int(q_check) for q_check, _, _, _ in group.get("deltas", [])})
    if deltas:
        return deltas
    if len(refs) > 1:
        return refs[:-1]
    errs = sorted({int(q) for q, _, _ in group.get("errs", [])})
    if len(errs) > 1:
        return errs[:-1]
    return errs


def _primary_delta_series(group: dict[str, object]) -> tuple[int | None, int | None, list[float]]:
    deltas = group.get("deltas", [])
    if deltas:
        q_check, q_reference, _, values = deltas[0]
        return int(q_check), int(q_reference), list(values)
    q_check, q_reference = _quadrature_orders_from_group(group)
    return q_check, q_reference, list(group.get("delta", []))


def _delta_series_for_q(
    group: dict[str, object],
    preferred_q: int | None,
) -> tuple[int | None, int | None, list[float]]:
    deltas = group.get("deltas", [])
    if preferred_q is not None:
        for q_check, q_reference, _name, values in deltas:
            if int(q_check) == int(preferred_q):
                return int(q_check), int(q_reference), list(values)
    return _primary_delta_series(group)


def _quadrature_difference_label(q_check: int | None, q_reference: int | None) -> str:
    if q_check is None or q_reference is None:
        return "quadrature difference"
    return f"quadrature difference |q={q_check} - q={q_reference}|"


def _column_label(name: str) -> str:
    ref_match = REF_PATTERN.fullmatch(name)
    if ref_match:
        return _reference_label(int(ref_match.group("q")))

    err_match = ERR_PATTERN.fullmatch(name)
    if err_match:
        return _reference_error_label(int(err_match.group("q")))

    if name == "deltaQuadrature":
        return "quadrature difference"

    return name


def _parse_vector_literal(literal: str) -> tuple[float, float, float]:
    values = tuple(float(token) for token in literal.strip().strip("()").split())
    if len(values) != 3:
        raise ValueError(f"Expected 3-vector literal, got: {literal}")
    return values


def _interpolate_series(
    source_times: list[float],
    source_values: list[float],
    target_times: list[float],
) -> list[float]:
    if not source_times or not source_values or len(source_times) != len(source_values):
        return [0.0 for _ in target_times]
    if len(source_times) == 1:
        return [source_values[0] for _ in target_times]

    result: list[float] = []
    for target in target_times:
        if target <= source_times[0]:
            result.append(source_values[0])
            continue
        if target >= source_times[-1]:
            result.append(source_values[-1])
            continue

        upper = bisect.bisect_left(source_times, target)
        lower = upper - 1
        t0 = source_times[lower]
        t1 = source_times[upper]
        v0 = source_values[lower]
        v1 = source_values[upper]
        weight = (target - t0) / (t1 - t0)
        result.append(v0 + weight * (v1 - v0))

    return result


def _surface_mesh(
    x_values: list[float],
    y_values: list[float],
    z_matrix: list[list[float]],
):
    x_grid, y_grid = np.meshgrid(np.asarray(x_values, dtype=float), np.asarray(y_values, dtype=float))
    z_grid = np.asarray(z_matrix, dtype=float)
    return x_grid, y_grid, z_grid


def _style_surface_axis(
    ax,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    zlabel: str,
    ytick_values: list[float] | None = None,
    ytick_labels: list[str] | None = None,
) -> None:
    ax.set_title(title)
    ax.set_xlabel(xlabel, labelpad=10)
    ax.set_ylabel(ylabel, labelpad=10)
    ax.set_zlabel(zlabel, labelpad=12)
    if ytick_values is not None:
        ax.set_yticks(ytick_values)
    if ytick_labels is not None:
        ax.set_yticklabels(ytick_labels)
    ax.tick_params(axis="x", pad=2, labelsize=8)
    ax.tick_params(axis="y", pad=2, labelsize=8)
    ax.tick_params(axis="z", pad=4, labelsize=8)
    ax.grid(True, linestyle="--", alpha=0.35)
    if hasattr(ax, "xaxis") and hasattr(ax.xaxis, "pane"):
        ax.xaxis.pane.set_facecolor((0.96, 0.96, 0.96, 0.65))
        ax.yaxis.pane.set_facecolor((0.97, 0.97, 0.97, 0.65))
        ax.zaxis.pane.set_facecolor((0.98, 0.98, 0.98, 0.45))
    ax.set_box_aspect((1.5, 1.0, 0.85))
    ax.view_init(elev=26, azim=-132)


def _plot_error_surface(
    ax,
    *,
    x_values: list[float],
    y_values: list[float],
    z_matrix: list[list[float]],
    title: str,
    ylabel: str,
    zlabel: str,
    cmap: str,
    ytick_labels: list[str] | None = None,
) -> None:
    x_grid, y_grid, z_grid = _surface_mesh(x_values, y_values, z_matrix)
    surface = ax.plot_surface(
        x_grid,
        y_grid,
        z_grid,
        cmap=cmap,
        linewidth=0.0,
        antialiased=True,
        alpha=0.92,
    )
    _style_surface_axis(
        ax,
        title=title,
        xlabel="time (s)",
        ylabel=ylabel,
        zlabel=zlabel,
        ytick_values=y_values,
        ytick_labels=ytick_labels,
    )
    return surface


def _prepare_log_surface_matrix(z_matrix: list[list[float]], floor: float):
    return [[max(value, floor) for value in row] for row in z_matrix]


def _plot_dual_error_surfaces(
    ax,
    *,
    x_values: list[float],
    y_values: list[float],
    primary_matrix: list[list[float]],
    secondary_matrix: list[list[float]],
    title: str,
    ylabel: str,
    ytick_labels: list[str] | None = None,
    primary_label: str,
    secondary_label: str,
) -> tuple[object, object]:
    floor = _positive_floor(primary_matrix + secondary_matrix)
    primary_safe = _prepare_log_surface_matrix(primary_matrix, floor)
    secondary_safe = _prepare_log_surface_matrix(secondary_matrix, floor)

    x_grid, y_grid, primary_grid = _surface_mesh(x_values, y_values, primary_safe)
    _, _, secondary_grid = _surface_mesh(x_values, y_values, secondary_safe)
    zmax = max(max(row) for row in (primary_safe + secondary_safe))
    if math.isclose(zmax, floor, rel_tol=0.0, abs_tol=1e-30):
        zmax = floor * 10.0
    norm = LogNorm(vmin=floor, vmax=zmax) if LogNorm is not None else None

    primary_surface = ax.plot_surface(
        x_grid,
        y_grid,
        primary_grid,
        cmap="viridis",
        norm=norm,
        linewidth=0.35,
        edgecolor="black",
        antialiased=True,
        alpha=0.86,
    )
    secondary_surface = ax.plot_surface(
        x_grid,
        y_grid,
        secondary_grid,
        cmap="viridis",
        norm=norm,
        linewidth=0.35,
        edgecolor="white",
        antialiased=True,
        alpha=0.52,
    )
    _style_surface_axis(
        ax,
        title=title,
        xlabel="time (s)",
        ylabel=ylabel,
        zlabel="error magnitude",
        ytick_values=y_values,
        ytick_labels=ytick_labels,
    )
    ax.set_zscale("log")
    ax.set_zlim(floor, zmax)

    if Patch is not None:
        ax.legend(
            handles=[
                Patch(
                    facecolor=plt.get_cmap("viridis")(0.75),
                    edgecolor="black",
                    linewidth=1.0,
                    alpha=0.86,
                    label=primary_label,
                ),
                Patch(
                    facecolor=plt.get_cmap("viridis")(0.45),
                    edgecolor="white",
                    linewidth=1.0,
                    alpha=0.52,
                    label=secondary_label,
                ),
            ],
            loc="upper left",
            fontsize=8,
        )
    return primary_surface, secondary_surface


def _slugify_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", value.strip())
    return token.strip("_").lower()


def _write_surface_polydata_vtp(
    destination: Path,
    *,
    x_values: list[float],
    y_values: list[float],
    z_matrix: list[list[float]],
    title: str,
    scalar_name: str,
) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)

    ny = len(y_values)
    nx = len(x_values)
    if ny == 0 or nx == 0:
        raise ValueError("Surface export requires non-empty x and y coordinates")
    if len(z_matrix) != ny or any(len(row) != nx for row in z_matrix):
        raise ValueError("Surface export requires z_matrix to match y/x dimensions")

    points: list[tuple[float, float, float]] = []
    values: list[float] = []
    for row_index, y_value in enumerate(y_values):
        row = z_matrix[row_index]
        for col_index, x_value in enumerate(x_values):
            z_value = float(row[col_index])
            points.append((float(x_value), float(y_value), z_value))
            values.append(z_value)

    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    zs = [point[2] for point in points]
    bounds = (
        min(xs),
        max(xs),
        min(ys),
        max(ys),
        min(zs),
        max(zs),
    )

    triangles: list[tuple[int, int, int]] = []
    if ny > 1 and nx > 1:
        for row_index in range(ny - 1):
            for col_index in range(nx - 1):
                p00 = row_index * nx + col_index
                p01 = p00 + 1
                p10 = (row_index + 1) * nx + col_index
                p11 = p10 + 1
                triangles.append((p00, p01, p11))
                triangles.append((p00, p11, p10))

    fill_scalar = min(values) if values else 0.0
    line_segments = _append_line_geometry(
        points,
        values,
        _box_edges_from_bounds(bounds) + _axis_lines_from_bounds(bounds),
        fill_scalar=fill_scalar,
    )

    connectivity = " ".join(f"{i0} {i1} {i2}" for i0, i1, i2 in triangles)
    offsets = " ".join(str(3 * (index + 1)) for index in range(len(triangles)))
    line_connectivity = " ".join(f"{i0} {i1}" for i0, i1 in line_segments)
    line_offsets = " ".join(str(2 * (index + 1)) for index in range(len(line_segments)))
    point_data = " ".join(f"{value:.9e}" for value in values)
    raw_x = " ".join(f"{x_value:.9e}" for x_value, _, _ in points)
    raw_y = " ".join(f"{y_value:.9e}" for _, y_value, _ in points)
    raw_z = " ".join(f"{z_value:.9e}" for _, _, z_value in points)
    normalized_points = _normalize_points(points, bounds)
    point_coords = " ".join(
        f"{x_value:.9e} {y_value:.9e} {z_value:.9e}" for x_value, y_value, z_value in points
    )
    normalized_point_coords = " ".join(
        f"{x_value:.9e} {y_value:.9e} {z_value:.9e}"
        for x_value, y_value, z_value in normalized_points
    )

    with destination.open("w", encoding="utf-8") as handle:
        handle.write("<?xml version=\"1.0\"?>\n")
        handle.write(
            "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        )
        handle.write("  <PolyData>\n")
        handle.write(
            f"    <Piece NumberOfPoints=\"{len(points)}\" NumberOfVerts=\"0\" "
            f"NumberOfLines=\"{len(line_segments)}\" NumberOfStrips=\"0\" NumberOfPolys=\"{len(triangles)}\">\n"
        )
        handle.write("      <PointData Scalars=\"{0}\">\n".format(scalar_name))
        handle.write(
            f"        <DataArray type=\"Float32\" Name=\"{scalar_name}\" format=\"ascii\">\n"
        )
        handle.write(f"          {point_data}\n")
        handle.write("        </DataArray>\n")
        handle.write("        <DataArray type=\"Float32\" Name=\"raw_x\" format=\"ascii\">\n")
        handle.write(f"          {raw_x}\n")
        handle.write("        </DataArray>\n")
        handle.write("        <DataArray type=\"Float32\" Name=\"raw_y\" format=\"ascii\">\n")
        handle.write(f"          {raw_y}\n")
        handle.write("        </DataArray>\n")
        handle.write("        <DataArray type=\"Float32\" Name=\"raw_z\" format=\"ascii\">\n")
        handle.write(f"          {raw_z}\n")
        handle.write("        </DataArray>\n")
        handle.write("      </PointData>\n")
        handle.write("      <Points>\n")
        handle.write("        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n")
        handle.write(f"          {normalized_point_coords}\n")
        handle.write("        </DataArray>\n")
        handle.write("      </Points>\n")
        handle.write("      <Lines>\n")
        handle.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
        handle.write(f"          {line_connectivity}\n")
        handle.write("        </DataArray>\n")
        handle.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
        handle.write(f"          {line_offsets}\n")
        handle.write("        </DataArray>\n")
        handle.write("      </Lines>\n")
        handle.write("      <Polys>\n")
        handle.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
        handle.write(f"          {connectivity}\n")
        handle.write("        </DataArray>\n")
        handle.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
        handle.write(f"          {offsets}\n")
        handle.write("        </DataArray>\n")
        handle.write("      </Polys>\n")
        handle.write("    </Piece>\n")
        handle.write("  </PolyData>\n")
        handle.write("</VTKFile>\n")

    return destination


def _export_dual_surface_polydata(
    *,
    base_path: Path,
    x_values: list[float],
    y_values: list[float],
    primary_matrix: list[list[float]],
    secondary_matrix: list[list[float]],
    title_prefix: str,
    primary_slug: str,
    secondary_slug: str,
) -> list[Path]:
    primary_path = base_path.with_name(f"{base_path.stem}_{primary_slug}.vtp")
    secondary_path = base_path.with_name(f"{base_path.stem}_{secondary_slug}.vtp")
    return [
        _write_surface_polydata_vtp(
            primary_path,
            x_values=x_values,
            y_values=y_values,
            z_matrix=primary_matrix,
            title=f"{title_prefix} - {primary_slug}",
            scalar_name=primary_slug,
        ),
        _write_surface_polydata_vtp(
            secondary_path,
            x_values=x_values,
            y_values=y_values,
            z_matrix=secondary_matrix,
            title=f"{title_prefix} - {secondary_slug}",
            scalar_name=secondary_slug,
        ),
    ]


def _unit_cube_edges() -> list[tuple[tuple[float, float, float], tuple[float, float, float]]]:
    return [
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
        ((0.0, 1.0, 0.0), (1.0, 1.0, 0.0)),
        ((0.0, 0.0, 1.0), (1.0, 0.0, 1.0)),
        ((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
        ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        ((1.0, 0.0, 0.0), (1.0, 1.0, 0.0)),
        ((0.0, 0.0, 1.0), (0.0, 1.0, 1.0)),
        ((1.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
        ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
        ((1.0, 0.0, 0.0), (1.0, 0.0, 1.0)),
        ((0.0, 1.0, 0.0), (0.0, 1.0, 1.0)),
        ((1.0, 1.0, 0.0), (1.0, 1.0, 1.0)),
    ]


def _box_edges_from_bounds(
    bounds: tuple[float, float, float, float, float, float]
) -> list[tuple[tuple[float, float, float], tuple[float, float, float]]]:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    p000 = (xmin, ymin, zmin)
    p100 = (xmax, ymin, zmin)
    p010 = (xmin, ymax, zmin)
    p110 = (xmax, ymax, zmin)
    p001 = (xmin, ymin, zmax)
    p101 = (xmax, ymin, zmax)
    p011 = (xmin, ymax, zmax)
    p111 = (xmax, ymax, zmax)
    return [
        (p000, p100),
        (p010, p110),
        (p001, p101),
        (p011, p111),
        (p000, p010),
        (p100, p110),
        (p001, p011),
        (p101, p111),
        (p000, p001),
        (p100, p101),
        (p010, p011),
        (p110, p111),
    ]


def _axis_lines_from_bounds(
    bounds: tuple[float, float, float, float, float, float]
) -> list[tuple[tuple[float, float, float], tuple[float, float, float]]]:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    origin = (xmin, ymin, zmin)
    return [
        (origin, (xmax, ymin, zmin)),
        (origin, (xmin, ymax, zmin)),
        (origin, (xmin, ymin, zmax)),
    ]


def _append_line_geometry(
    points: list[tuple[float, float, float]],
    scalar_values: list[float],
    lines: list[tuple[tuple[float, float, float], tuple[float, float, float]]],
    *,
    fill_scalar: float,
) -> list[tuple[int, int]]:
    index_map = {
        (float(x_value), float(y_value), float(z_value)): index
        for index, (x_value, y_value, z_value) in enumerate(points)
    }
    connectivity: list[tuple[int, int]] = []

    for start_point, end_point in lines:
        start_key = tuple(float(value) for value in start_point)
        end_key = tuple(float(value) for value in end_point)
        if start_key not in index_map:
            index_map[start_key] = len(points)
            points.append(start_key)
            scalar_values.append(fill_scalar)
        if end_key not in index_map:
            index_map[end_key] = len(points)
            points.append(end_key)
            scalar_values.append(fill_scalar)
        connectivity.append((index_map[start_key], index_map[end_key]))

    return connectivity


def _normalize_points(
    points: list[tuple[float, float, float]],
    bounds: tuple[float, float, float, float, float, float],
) -> list[tuple[float, float, float]]:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    xspan = xmax - xmin
    yspan = ymax - ymin
    zspan = zmax - zmin

    def _normalize(value: float, lower: float, span: float) -> float:
        if math.isclose(span, 0.0, rel_tol=0.0, abs_tol=1e-30):
            return 0.0
        return (value - lower) / span

    return [
        (
            _normalize(x_value, xmin, xspan),
            _normalize(y_value, ymin, yspan),
            _normalize(z_value, zmin, zspan),
        )
        for x_value, y_value, z_value in points
    ]


def export_ecg_electrode_geometry_vtp(
    *,
    save_path: str | Path | None = None,
):
    if save_path is None:
        return None

    destination = Path(save_path)
    destination.parent.mkdir(parents=True, exist_ok=True)

    cube_vertices = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (0.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
    ]
    vertex_index = {vertex: index for index, vertex in enumerate(cube_vertices)}
    cube_lines = [(vertex_index[start], vertex_index[end]) for start, end in _unit_cube_edges()]

    electrode_items = sorted(driver_defaults.ECG_ELECTRODES_BY_DIMENSION["3D"].items())
    electrode_points = [_parse_vector_literal(literal) for _, literal in electrode_items]

    points = cube_vertices + electrode_points
    electrode_offset = len(cube_vertices)
    electrode_indices = [electrode_offset + index for index in range(len(electrode_points))]
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    zs = [point[2] for point in points]
    bounds = (
        min(xs),
        max(xs),
        min(ys),
        max(ys),
        min(zs),
        max(zs),
    )

    line_points = list(points)
    scalar_values = [0.0 for _ in cube_vertices] + [float(index + 1) for index in range(len(electrode_points))]
    cube_and_axis_lines = cube_lines + _append_line_geometry(
        line_points,
        scalar_values,
        _axis_lines_from_bounds(bounds),
        fill_scalar=0.0,
    )
    raw_x = " ".join(f"{x_value:.9e}" for x_value, _, _ in line_points)
    raw_y = " ".join(f"{y_value:.9e}" for _, y_value, _ in line_points)
    raw_z = " ".join(f"{z_value:.9e}" for _, _, z_value in line_points)
    normalized_points = _normalize_points(line_points, bounds)
    point_coords = " ".join(
        f"{x_value:.9e} {y_value:.9e} {z_value:.9e}"
        for x_value, y_value, z_value in normalized_points
    )
    line_connectivity = " ".join(f"{start} {end}" for start, end in cube_and_axis_lines)
    line_offsets = " ".join(str(2 * (index + 1)) for index in range(len(cube_and_axis_lines)))
    vert_connectivity = " ".join(str(index) for index in electrode_indices)
    vert_offsets = " ".join(str(index + 1) for index in range(len(electrode_indices)))
    entity_ids = " ".join(f"{value:.9e}" for value in scalar_values)

    with destination.open("w", encoding="utf-8") as handle:
        handle.write("<?xml version=\"1.0\"?>\n")
        handle.write(
            "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        )
        handle.write("  <PolyData>\n")
        handle.write(
            f"    <Piece NumberOfPoints=\"{len(line_points)}\" NumberOfVerts=\"{len(electrode_indices)}\" "
            f"NumberOfLines=\"{len(cube_and_axis_lines)}\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
        )
        handle.write("      <PointData Scalars=\"entity_id\">\n")
        handle.write("        <DataArray type=\"Float32\" Name=\"entity_id\" format=\"ascii\">\n")
        handle.write(f"          {entity_ids}\n")
        handle.write("        </DataArray>\n")
        handle.write("        <DataArray type=\"Float32\" Name=\"raw_x\" format=\"ascii\">\n")
        handle.write(f"          {raw_x}\n")
        handle.write("        </DataArray>\n")
        handle.write("        <DataArray type=\"Float32\" Name=\"raw_y\" format=\"ascii\">\n")
        handle.write(f"          {raw_y}\n")
        handle.write("        </DataArray>\n")
        handle.write("        <DataArray type=\"Float32\" Name=\"raw_z\" format=\"ascii\">\n")
        handle.write(f"          {raw_z}\n")
        handle.write("        </DataArray>\n")
        handle.write("      </PointData>\n")
        handle.write("      <Points>\n")
        handle.write("        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n")
        handle.write(f"          {point_coords}\n")
        handle.write("        </DataArray>\n")
        handle.write("      </Points>\n")
        handle.write("      <Verts>\n")
        handle.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
        handle.write(f"          {vert_connectivity}\n")
        handle.write("        </DataArray>\n")
        handle.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
        handle.write(f"          {vert_offsets}\n")
        handle.write("        </DataArray>\n")
        handle.write("      </Verts>\n")
        handle.write("      <Lines>\n")
        handle.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
        handle.write(f"          {line_connectivity}\n")
        handle.write("        </DataArray>\n")
        handle.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
        handle.write(f"          {line_offsets}\n")
        handle.write("        </DataArray>\n")
        handle.write("      </Lines>\n")
        handle.write("    </Piece>\n")
        handle.write("  </PolyData>\n")
        handle.write("</VTKFile>\n")

    return destination


def organize_dat_files(folder_name):
    """
    Move all .dat files from the parent directory into a subfolder
    inside the parent directory.

    Parameters:
        folder_name (str): Name of the folder inside the parent directory.
    """
    # Absolute path of parent directory
    parent_dir = os.path.abspath("..")
    dest_dir = os.path.join(parent_dir, folder_name)
    print(f"Looking for .dat files in parent directory: {parent_dir}")

    # List all .dat files in parent directory
    dat_files = [f for f in os.listdir(parent_dir) if f.endswith(".dat")]

    if not dat_files:
        print("No .dat files found in the parent directory.")
        return

    # Create destination folder inside parent directory
    os.makedirs(dest_dir, exist_ok=True)

    # Move each .dat file to the folder
    for f in dat_files:
        src = os.path.join(parent_dir, f)
        dst = os.path.join(dest_dir, f)
        shutil.move(src, dst)
        print(f"Moved {f} -> {dest_dir}/")

    print(f"\nAll {len(dat_files)} .dat files moved to '{dest_dir}'.")




def read_error_dat_files(folder_name, *, expected_filenames: set[str] | None = None):
    """
    Reads all .dat files in folder_name and extracts:
        - Dimension  (1D, 2D, 3D)
        - N          (# cells)
        - Solver     (explicit, implicit)
        - Linf errors for Vm, u1, u2

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

    if expected_filenames is not None:
        files = [f for f in files if f.name in expected_filenames]
        if not files:
            print("No expected .dat files found in folder:", folder)
            return []

    data = []

    for f in files:
        # Expected filename format:
        #   1D_320_cells_explicit.dat
        m = FILENAME_PATTERN.match(f.name)
        if not m:
            print("Skipping unrecognized filename:", f.name)
            continue

        dimension = m.group(1)   # "1D"
        N = int(m.group(2))      # 320
        solver = m.group(3)      # "explicit" or "implicit"

        content = f.read_text()

        # Extract Linf errors
        linf_matches = re.findall(
            r"Vm\s+\S+\s+\S+\s+(\S+).*?"
            r"u1\s+\S+\s+\S+\s+(\S+).*?"
            r"u2\s+\S+\s+\S+\s+(\S+)",
            content.replace("\n", " "), re.DOTALL
        )

        if linf_matches:
            Linf_V, Linf_u1, Linf_u2 = map(float, linf_matches[0])
        else:
            Linf_V = Linf_u1 = Linf_u2 = float("nan")

        data.append({
            "Dimension": dimension,
            "N": N,
            "Solver": solver,
            "Linf_V": Linf_V,
            "Linf_u1": Linf_u1,
            "Linf_u2": Linf_u2
        })

    return sorted(data, key=lambda row: (row["Dimension"], row["Solver"], row["N"]))


def read_ecg_summary_dat_files(folder_name):
    """Read archived ECG manufactured summaries and aggregate per-case metrics."""

    folder = Path(folder_name)
    if not folder.exists():
        print("Folder does not exist:", folder)
        return []

    rows = []
    for path in sorted(folder.glob("ECG_*_manufacturedPseudoECGSummary.dat")):
        match = ECG_SUMMARY_PATTERN.match(path.name)
        if not match:
            print("Skipping unrecognized ECG summary filename:", path.name)
            continue

        metadata, electrode_rows = _parse_ecg_summary_file(path)
        if not electrode_rows:
            continue

        q_checks = [
            int(token)
            for token in metadata.get("qChecks", metadata.get("qCheck", "")).split()
            if token.strip()
        ]
        q_check = min(q_checks) if q_checks else 0

        def _metric_values(prefix: str) -> list[float]:
            values: list[float] = []
            for row in electrode_rows:
                for key, value in row.items():
                    if isinstance(key, str) and key.startswith(prefix):
                        values.append(float(value))
            return values

        delta_by_q: dict[int, dict[str, float]] = {}
        err_by_q: dict[int, dict[str, float]] = {}
        for q_value in q_checks:
            delta_values = [
                float(row[key])
                for row in electrode_rows
                for key in row
                if isinstance(key, str)
                and SUMMARY_DELTA_Q_PATTERN.fullmatch(key)
                and int(SUMMARY_DELTA_Q_PATTERN.fullmatch(key).group("q")) == q_value
            ]
            if delta_values:
                delta_by_q[q_value] = {
                    "max": max(delta_values),
                    "mean": sum(delta_values) / len(delta_values),
                }

            err_values = [
                float(row[key])
                for row in electrode_rows
                for key in row
                if isinstance(key, str)
                and SUMMARY_ERR_Q_PATTERN.fullmatch(key)
                and int(SUMMARY_ERR_Q_PATTERN.fullmatch(key).group("q")) == q_value
            ]
            if err_values:
                err_by_q[q_value] = {
                    "max": max(err_values),
                    "mean": sum(err_values) / len(err_values),
                }

        def _mean(key: str) -> float:
            return sum(float(row[key]) for row in electrode_rows) / len(electrode_rows)

        def _max(key: str) -> float:
            return max(float(row[key]) for row in electrode_rows)

        linf_err_check_values = _metric_values("Linf_err_q")
        linf_delta_values = _metric_values("Linf_delta_q")

        rows.append(
            {
                "Dimension": match.group("dimension"),
                "N": int(match.group("cells")),
                "Solver": match.group("solver"),
                "samples": int(metadata.get("samples", "0")),
                "electrodes": len(electrode_rows),
                "qCheck": q_check,
                "qChecks": " ".join(str(value) for value in q_checks),
                "qReference": int(metadata.get("qReference", "0")),
                "max_Linf_err_ref": _max("Linf_err_ref"),
                "mean_Linf_err_ref": _mean("Linf_err_ref"),
                "max_Linf_err_check": max(linf_err_check_values, default=float("nan")),
                "mean_Linf_err_check": (
                    sum(linf_err_check_values) / len(linf_err_check_values)
                    if linf_err_check_values
                    else float("nan")
                ),
                "max_Linf_delta_ref": max(linf_delta_values, default=float("nan")),
                "mean_Linf_delta_ref": (
                    sum(linf_delta_values) / len(linf_delta_values)
                    if linf_delta_values
                    else float("nan")
                ),
                "delta_by_q": delta_by_q,
                "err_by_q": err_by_q,
            }
        )

    return sorted(rows, key=lambda row: (row["Dimension"], row["Solver"], row["N"]))


def read_ecg_timeseries_dat_files(folder_name):
    folder = Path(folder_name)
    if not folder.exists():
        print("Folder does not exist:", folder)
        return []

    cases = []
    for path in sorted(folder.glob("ECG_*_manufacturedPseudoECG.dat")):
        match = ECG_TIMESERIES_PATTERN.match(path.name)
        if not match:
            print("Skipping unrecognized ECG timeseries filename:", path.name)
            continue

        header, rows = _parse_timeseries(path)
        columns = _columns_from_rows(header, rows)
        cases.append(
            {
                "path": path,
                "Dimension": match.group("dimension"),
                "N": int(match.group("cells")),
                "Solver": match.group("solver"),
                "times": columns["time"],
                "groups": _group_manufactured_columns(columns),
            }
        )

    return sorted(cases, key=lambda row: (row["Dimension"], row["Solver"], row["N"]))


def _filter_supported_ecg_rows(rows, *, source_label: str):
    supported = []
    skipped_dimensions = set()
    for row in rows:
        dimension = str(row["Dimension"])
        if dimension in SUPPORTED_ECG_POSTPROCESS_DIMENSIONS:
            supported.append(row)
        else:
            skipped_dimensions.add(dimension)

    if skipped_dimensions:
        supported_list = ", ".join(SUPPORTED_ECG_POSTPROCESS_DIMENSIONS)
        skipped_list = ", ".join(sorted(skipped_dimensions))
        print(
            "Warning: skipping ECG post-processing for dimensions "
            f"{skipped_list}. Supported ECG post-processing dimensions: {supported_list}. "
            "Reason: the numerical pseudoECG is accumulated as a 3D cell-volume sum, "
            "while the current 1D/2D manufactured ECG references are lower-dimensional "
            "integrals and are not a consistent verification target."
        )

    if not supported:
        print(f"No supported ECG {source_label} data found after filtering.")

    return supported


def _cleanup_unsupported_ecg_archives(output_dir: Path) -> list[Path]:
    removed: list[Path] = []
    for path in sorted(output_dir.glob("ECG_*.dat")):
        match = re.match(r"ECG_(?P<dimension>\dD)_", path.name)
        if match is None:
            continue
        dimension = match.group("dimension")
        if dimension in SUPPORTED_ECG_POSTPROCESS_DIMENSIONS:
            continue
        path.unlink()
        removed.append(path)

    if removed:
        removed_names = ", ".join(path.name for path in removed)
        supported_list = ", ".join(SUPPORTED_ECG_POSTPROCESS_DIMENSIONS)
        print(
            "Removed unsupported ECG archived cases from postProcessing: "
            f"{removed_names}. Supported ECG post-processing dimensions: {supported_list}. "
            "Reason: ECG verification is currently defined only for 3D because the "
            "numerical pseudoECG is accumulated as a 3D cell-volume sum, while the "
            "current 1D/2D manufactured ECG references are lower-dimensional integrals."
        )

    return removed


def _cleanup_stale_ecg_plot_artifacts(output_dir: Path) -> list[Path]:
    stale_patterns = (
        "manufactured_ecg_*_1d_*.png",
        "manufactured_ecg_*_2d_*.png",
        "manufactured_ecg_*_1d_*.vtk",
        "manufactured_ecg_*_2d_*.vtk",
        "manufactured_ecg_*_1d_*.vtp",
        "manufactured_ecg_*_2d_*.vtp",
        "manufactured_ecg_max_reference_error_across_dimensions.png",
        "manufactured_ecg_max_reference_quadrature_gap_across_dimensions.png",
        "manufactured_ecg_summary_dashboard.png",
        "manufactured_ecg_error_vs_gap_dashboard.png",
    )

    removed: list[Path] = []
    for pattern in stale_patterns:
        for path in sorted(output_dir.glob(pattern)):
            if not path.is_file():
                continue
            path.unlink()
            removed.append(path)

    if removed:
        print(
            "Removed stale ECG plot artifacts from postProcessing: "
            + ", ".join(path.name for path in removed)
        )

    return removed


def select_representative_ecg_cases(cases):
    grouped = {}
    for case in cases:
        key = (case["Dimension"], case["Solver"])
        current = grouped.get(key)
        if current is None or int(case["N"]) > int(current["N"]):
            grouped[key] = case
    return [grouped[key] for key in sorted(grouped)]


def group_ecg_timeseries_cases(cases):
    grouped = {}
    for case in cases:
        key = (case["Dimension"], case["Solver"])
        grouped.setdefault(key, []).append(case)
    for key in grouped:
        grouped[key] = sorted(grouped[key], key=lambda row: int(row["N"]))
    return grouped


def compute_convergence_rates(rows):
    """
    Compute convergence rates for Linf errors of Vm, u1, u2.

    - Groups by Dimension (if present) and Solver.
    - Sorts by N.
    - Skips pairs where N_lower == N_higher.
    """

    grouped_rows = {}
    for row in rows:
        key = (row["Dimension"], row["Solver"])
        grouped_rows.setdefault(key, []).append(row)

    convergence_rows = []
    for (dimension, solver_type), group_rows in sorted(grouped_rows.items()):
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
                    "Solver": solver_type,
                    "N_lower": N1,
                    "N_higher": N2,
                    "rate_Vm": _safe_rate(lower["Linf_V"], higher["Linf_V"], h1, h2),
                    "rate_u1": _safe_rate(lower["Linf_u1"], higher["Linf_u1"], h1, h2),
                    "rate_u2": _safe_rate(lower["Linf_u2"], higher["Linf_u2"], h1, h2),
                }
            )

    return convergence_rows


def compute_ecg_convergence_rates(rows):
    grouped_rows = {}
    for row in rows:
        key = (row["Dimension"], row["Solver"])
        grouped_rows.setdefault(key, []).append(row)

    convergence_rows = []
    for (dimension, solver_type), group_rows in sorted(grouped_rows.items()):
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
                    "Solver": solver_type,
                    "N_lower": N1,
                    "N_higher": N2,
                    "rate_max_Linf_err_ref": _safe_rate(
                        lower["max_Linf_err_ref"], higher["max_Linf_err_ref"], h1, h2
                    ),
                    "rate_mean_Linf_err_ref": _safe_rate(
                        lower["mean_Linf_err_ref"], higher["mean_Linf_err_ref"], h1, h2
                    ),
                    "rate_max_Linf_delta_ref": _safe_rate(
                        lower["max_Linf_delta_ref"], higher["max_Linf_delta_ref"], h1, h2
                    ),
                    "rate_mean_Linf_delta_ref": _safe_rate(
                        lower["mean_Linf_delta_ref"], higher["mean_Linf_delta_ref"], h1, h2
                    ),
                }
            )

    return convergence_rows


def _finalize_axis_legend(ax) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend()


def _plot_dimension_errors_on_axis(ax, rows, dimension: str) -> bool:
    dimension_rows = _filter_rows(rows, Dimension=dimension)
    if not dimension_rows:
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
            title=f"Linf Errors per Dimension ({dimension})",
            xlabel="Number of cells (N)",
            ylabel="Linf Error",
            legend=False,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.6},
        )
        return False

    for solver in _unique_values(dimension_rows, "Solver"):
        solver_rows = _filter_rows(dimension_rows, Solver=solver)
        if not solver_rows:
            continue

        for col, color in FIELD_COLORS.items():
            ax.loglog(
                [row["N"] for row in solver_rows],
                [row[col] for row in solver_rows],
                marker=SOLVER_MARKERS.get(solver, "o"),
                linestyle=SOLVER_LINESTYLES.get(solver, "-"),
                color=color,
                label=f"{FIELD_LABELS[col]} ({solver}, {dimension})",
            )

    style_matplotlib_axes(
        ax,
        title=f"Linf Errors per Dimension ({dimension})",
        xlabel="Number of cells (N)",
        ylabel="Linf Error",
        legend=False,
        grid_kwargs={"which": "both", "ls": "--", "alpha": 0.6},
    )
    _finalize_axis_legend(ax)
    return True


def _plot_vm_across_dimensions_on_axis(ax, rows) -> bool:
    plotted = False
    for dimension in _unique_values(rows, "Dimension"):
        dimension_rows = _filter_rows(rows, Dimension=dimension)

        for solver in _unique_values(dimension_rows, "Solver"):
            solver_rows = _filter_rows(dimension_rows, Solver=solver)
            if not solver_rows:
                continue

            ax.loglog(
                [row["N"] for row in solver_rows],
                [row["Linf_V"] for row in solver_rows],
                marker=SOLVER_MARKERS.get(solver, "o"),
                linestyle=SOLVER_LINESTYLES.get(solver, "--"),
                color=DIMENSION_COLORS.get(dimension, "black"),
                label=f"{dimension} ({solver})",
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


def _plot_ecg_metric_on_axis(ax, rows, *, value_key: str, title: str, ylabel: str) -> bool:
    plotted = False
    for dimension in _unique_values(rows, "Dimension"):
        dimension_rows = _filter_rows(rows, Dimension=dimension)

        for solver in _unique_values(dimension_rows, "Solver"):
            solver_rows = _filter_rows(dimension_rows, Solver=solver)
            if not solver_rows:
                continue

            ax.loglog(
                [row["N"] for row in solver_rows],
                [row[value_key] for row in solver_rows],
                marker=SOLVER_MARKERS.get(solver, "o"),
                linestyle=SOLVER_LINESTYLES.get(solver, "--"),
                color=DIMENSION_COLORS.get(dimension, "black"),
                label=f"{dimension} ({solver})",
            )
            plotted = True

    if not plotted:
        ax.text(
            0.5,
            0.5,
            "No ECG data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color="0.4",
        )

    style_matplotlib_axes(
        ax,
        title=title,
        xlabel="Number of cells (N)",
        ylabel=ylabel,
        legend=False,
        grid_kwargs={"which": "both", "ls": "--", "alpha": 0.6},
    )
    _finalize_axis_legend(ax)
    return plotted


def plot_ecg_quadrature_summary(
    rows,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_matplotlib():
        print("matplotlib is not available; skipping ECG quadrature-only summary.")
        return None
    if not rows:
        return None

    configure_matplotlib_defaults()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.0), squeeze=False)
    flat_axes = axes.flatten()
    q_reference = None
    q_checks = sorted(
        {
            q_value
            for row in rows
            for q_value in row.get("delta_by_q", {}).keys()
        }
    )
    if rows:
        q_references = sorted({int(row["qReference"]) for row in rows if int(row["qReference"]) > 0})
        q_reference = q_references[0] if q_references else None

    for axis, reducer_key, title_prefix in (
        (flat_axes[0], "max", "Max"),
        (flat_axes[1], "mean", "Mean"),
    ):
        plotted = False
        for solver in _unique_values(rows, "Solver"):
            solver_rows = _filter_rows(rows, Solver=solver)
            for q_check in q_checks:
                selected_rows = [
                    row for row in solver_rows
                    if q_check in row.get("delta_by_q", {})
                ]
                if not selected_rows:
                    continue
                axis.loglog(
                    [row["N"] for row in selected_rows],
                    [row["delta_by_q"][q_check][reducer_key] for row in selected_rows],
                    marker=SOLVER_MARKERS.get(solver, "o"),
                    linestyle=SOLVER_LINESTYLES.get(solver, "-"),
                    label=f"q={q_check} vs qRef={q_reference} ({solver})",
                )
                plotted = True

        if not plotted:
            axis.text(
                0.5,
                0.5,
                "No quadrature summary data",
                transform=axis.transAxes,
                ha="center",
                va="center",
                fontsize=11,
                color="0.4",
            )

        style_matplotlib_axes(
            axis,
            title=f"{title_prefix} quadrature difference vs N",
            xlabel="Number of cells (N)",
            ylabel="quadrature difference",
            legend=True,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
        )

    fig.suptitle(
        "Manufactured pseudoECG: quadrature-only convergence summary",
        fontsize=12,
    )
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def plot_ecg_quadrature_reference_overlay_case(
    case,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_matplotlib():
        print("matplotlib is not available; skipping ECG quadrature reference overlay plot.")
        return None

    groups = case["groups"]
    if not groups:
        return None

    configure_matplotlib_defaults()
    names = list(groups.keys())
    nrows, ncols = _subplot_shape(len(names))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(4.6 * ncols, 3.2 * nrows),
        sharex=True,
        squeeze=False,
    )
    flat_axes = axes.flatten()
    y_limits = _common_y_limits(
        [
            values
            for group in groups.values()
            for _, _, values in group.get("refs", [])
        ]
    )

    for idx, electrode in enumerate(names):
        ax = flat_axes[idx]
        group = groups[electrode]
        for q_value, _name, values in group.get("refs", []):
            line_kwargs = {"linewidth": 1.15}
            if q_value == group.get("refs", [])[-1][0]:
                line_kwargs["linestyle"] = "--"
                line_kwargs["linewidth"] = 1.35
            ax.plot(
                case["times"],
                values,
                label=_reference_label(q_value),
                **line_kwargs,
            )
        if y_limits is not None:
            ax.set_ylim(*y_limits)
        style_matplotlib_axes(
            ax,
            title=electrode,
            xlabel="time (s)",
            ylabel="pseudoECG reference",
            legend=True,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
        )
        ax.tick_params(labelsize=8)

    for idx in range(len(names), len(flat_axes)):
        flat_axes[idx].axis("off")

    fig.suptitle(
        f"Manufactured pseudoECG: quadrature-only reference overlays ({case['Dimension']}, {case['Solver']}, N={case['N']})",
        fontsize=12,
    )
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def plot_ecg_quadrature_final_time_heatmap(
    cases,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_surface_plotting():
        print("matplotlib/numpy is not available; skipping ECG quadrature final-time heatmap.")
        return None
    if not cases:
        return None

    configure_matplotlib_defaults()
    grouped_cases = group_ecg_timeseries_cases(cases)
    ncols = max(1, len(grouped_cases))
    fig, axes = plt.subplots(1, ncols, figsize=(6.2 * ncols, 5.2), squeeze=False)
    flat_axes = axes.flatten()

    for axis, ((dimension, solver), group_cases) in zip(flat_axes, sorted(grouped_cases.items())):
        group_cases = sorted(group_cases, key=lambda row: int(row["N"]))
        q_checks = sorted(
            {
                int(q_check)
                for case in group_cases
                for group in case["groups"].values()
                for q_check, _q_ref, _name, _values in group.get("deltas", [])
            }
        )
        if not q_checks:
            axis.text(
                0.5,
                0.5,
                "No quadrature-difference data",
                transform=axis.transAxes,
                ha="center",
                va="center",
                fontsize=11,
                color="0.4",
            )
            style_matplotlib_axes(
                axis,
                title=f"Final-time max quadrature difference ({dimension}, {solver})",
                xlabel="N",
                ylabel="qCheck",
                legend=False,
                grid_kwargs={"which": "both", "ls": "--", "alpha": 0.30},
            )
            continue

        n_values = [int(case["N"]) for case in group_cases]
        matrix: list[list[float]] = []
        for q_check in q_checks:
            row_values: list[float] = []
            for case in group_cases:
                final_values: list[float] = []
                for group in case["groups"].values():
                    for current_q, _q_ref, _name, values in group.get("deltas", []):
                        if int(current_q) == q_check and values:
                            final_values.append(float(values[-1]))
                row_values.append(max(final_values) if final_values else 0.0)
            matrix.append(row_values)

        floor = _positive_floor(matrix)
        safe_matrix = np.asarray(_prepare_log_surface_matrix(matrix, floor), dtype=float)
        norm = LogNorm(vmin=floor, vmax=float(np.max(safe_matrix))) if LogNorm is not None else None
        image = axis.imshow(
            safe_matrix,
            aspect="auto",
            origin="lower",
            interpolation="nearest",
            cmap="viridis",
            norm=norm,
        )
        axis.set_xticks(range(len(n_values)))
        axis.set_xticklabels([str(value) for value in n_values])
        axis.set_yticks(range(len(q_checks)))
        axis.set_yticklabels([str(value) for value in q_checks])
        style_matplotlib_axes(
            axis,
            title=f"Final-time max quadrature difference ({dimension}, {solver})",
            xlabel="N",
            ylabel="qCheck",
            legend=False,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.25},
        )
        fig.colorbar(
            image,
            ax=axis,
            pad=0.02,
            shrink=0.85,
            label="max over electrodes at final time",
        )

    fig.suptitle(
        "Manufactured pseudoECG: final-time max quadrature-difference heatmap",
        fontsize=12,
    )
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def plot_ecg_reference_overlay_case(
    case,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_matplotlib():
        print("matplotlib is not available; skipping ECG reference overlay plot.")
        return None

    groups = case["groups"]
    if not groups:
        return None

    configure_matplotlib_defaults()
    names = list(groups.keys())
    nrows, ncols = _subplot_shape(len(names))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(4.4 * ncols, 3.2 * nrows),
        sharex=True,
        squeeze=False,
    )
    flat_axes = axes.flatten()
    y_limits = _common_y_limits(
        [
            series
            for group in groups.values()
            for series in (
                [group["numeric"]]
                + ([group["refs"][-1][2]] if group["refs"] else [])
            )
        ]
    )

    for idx, electrode in enumerate(names):
        ax = flat_axes[idx]
        group = groups[electrode]
        ax.plot(case["times"], group["numeric"], linewidth=1.3, label="numeric")
        if group["refs"]:
            q, name, values = group["refs"][-1]
            del name
            ax.plot(
                case["times"],
                values,
                linewidth=1.1,
                linestyle="--",
                label=_reference_label(q),
            )
        if y_limits is not None:
            ax.set_ylim(*y_limits)
        style_matplotlib_axes(
            ax,
            title=electrode,
            xlabel="time (s)",
            ylabel="pseudoECG",
            legend=True,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
        )
        ax.tick_params(labelsize=8)

    for idx in range(len(names), len(flat_axes)):
        flat_axes[idx].axis("off")

    fig.suptitle(
        f"Manufactured pseudoECG: numeric vs reference(qReference) ({case['Dimension']}, {case['Solver']}, N={case['N']})",
        fontsize=12,
    )
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def plot_ecg_error_timeseries_case(
    case,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_matplotlib():
        print("matplotlib is not available; skipping ECG error-timeseries plot.")
        return None

    groups = case["groups"]
    if not groups:
        return None

    configure_matplotlib_defaults()
    names = list(groups.keys())
    nrows, ncols = _subplot_shape(len(names))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(4.4 * ncols, 3.2 * nrows),
        sharex=True,
        squeeze=False,
    )
    flat_axes = axes.flatten()
    floor = _positive_floor(
        [
            series
            for group in groups.values()
            for series in ([values for _, _, values in group["errs"]] + [group["delta"]])
        ]
    )

    for idx, electrode in enumerate(names):
        ax = flat_axes[idx]
        group = groups[electrode]
        if group["errs"]:
            q, name, values = group["errs"][-1]
            del name
            ax.semilogy(
                case["times"],
                [max(value, floor) for value in values],
                linewidth=1.2,
                label=_reference_error_label(q),
            )
        deltas = group.get("deltas", [])
        if deltas:
          for q_check, q_reference, _name, values in deltas:
              ax.semilogy(
                  case["times"],
                  [max(value, floor) for value in values],
                  linewidth=1.1,
                  linestyle=":",
                  label=_quadrature_difference_label(q_check, q_reference),
              )
        else:
          q_check, q_reference = _quadrature_orders_from_group(group)
          ax.semilogy(
              case["times"],
              [max(value, floor) for value in group["delta"]],
              linewidth=1.2,
              linestyle=":",
              label=_quadrature_difference_label(q_check, q_reference),
          )
        style_matplotlib_axes(
            ax,
            title=electrode,
            xlabel="time (s)",
            ylabel="abs. error",
            legend=True,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
        )
        ax.tick_params(labelsize=8)

    for idx in range(len(names), len(flat_axes)):
        flat_axes[idx].axis("off")

    fig.suptitle(
        f"Manufactured pseudoECG: abs. error to reference(qReference) and quadrature difference vs time ({case['Dimension']}, {case['Solver']}, N={case['N']})",
        fontsize=12,
    )
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def plot_ecg_error_heatmap_case(
    case,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_surface_plotting():
        print("matplotlib/numpy is not available; skipping ECG error surface plot.")
        return None

    groups = case["groups"]
    if not groups:
        return None

    configure_matplotlib_defaults()
    names = list(groups.keys())
    err_matrix = []
    delta_matrix = []
    q_check = None
    q_reference = None
    for electrode in names:
        group = groups[electrode]
        if group["errs"]:
            err_matrix.append(group["errs"][-1][2])
        else:
            err_matrix.append([0.0 for _ in case["times"]])
        current_q_check, current_q_reference, current_delta = _primary_delta_series(group)
        delta_matrix.append(current_delta)
        if q_check is None or q_reference is None:
            q_check, q_reference = current_q_check, current_q_reference

    fig = plt.figure(figsize=(10.5, 6.2))
    axis = fig.add_subplot(1, 1, 1, projection="3d")
    electrode_indices = list(range(len(names)))

    primary_surface, _secondary_surface = _plot_dual_error_surfaces(
        axis,
        x_values=list(case["times"]),
        y_values=electrode_indices,
        primary_matrix=err_matrix,
        secondary_matrix=delta_matrix,
        title="Numerical/reference(qReference) error and quadrature-difference surfaces",
        ylabel="Electrode",
        ytick_labels=names,
        primary_label="numerical/reference(qReference) error",
        secondary_label=_quadrature_difference_label(q_check, q_reference),
    )
    fig.colorbar(
        primary_surface,
        ax=axis,
        pad=0.02,
        shrink=0.78,
        label="shared error heatmap",
    )

    fig.suptitle(
        f"Manufactured pseudoECG time-electrode error surfaces ({case['Dimension']}, {case['Solver']}, N={case['N']}, log z-scale)",
        fontsize=12,
    )
    fig.subplots_adjust(left=0.04, right=0.88, bottom=0.08, top=0.90)
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def export_ecg_error_surface_case_vtp(
    case,
    *,
    save_dir: str | Path | None = None,
):
    if save_dir is None:
        return []

    groups = case["groups"]
    if not groups:
        return []

    names = list(groups.keys())
    err_matrix = []
    delta_matrix = []
    q_check = None
    q_reference = None
    for electrode in names:
        group = groups[electrode]
        if group["errs"]:
            err_matrix.append(group["errs"][-1][2])
        else:
            err_matrix.append([0.0 for _ in case["times"]])
        delta_matrix.append(group["delta"])
        if q_check is None or q_reference is None:
            q_check, q_reference = _quadrature_orders_from_group(group)

    base_path = (
        Path(save_dir)
        / f"manufactured_ecg_reference_error_surface_{str(case['Dimension']).lower()}_{str(case['Solver']).lower()}"
    )
    return _export_dual_surface_polydata(
        base_path=base_path,
        x_values=list(case["times"]),
        y_values=list(range(len(names))),
        primary_matrix=err_matrix,
        secondary_matrix=delta_matrix,
        title_prefix=(
            f"Manufactured pseudoECG representative surface ({case['Dimension']}, {case['Solver']}, N={case['N']})"
        ),
        primary_slug="numerical_reference_error",
        secondary_slug=_slugify_token(_quadrature_difference_label(q_check, q_reference)),
    )


def plot_ecg_error_heatmap_sweep(
    cases,
    *,
    save_dir: str | Path | None = None,
    show: bool = True,
):
    if not _has_surface_plotting():
        print("matplotlib/numpy is not available; skipping ECG sweep surfaces.")
        return []
    if not cases:
        return []

    configure_matplotlib_defaults()
    output_paths: list[Path] = []
    grouped_cases = group_ecg_timeseries_cases(cases)

    for (dimension, solver), group_cases in grouped_cases.items():
        group_cases = sorted(group_cases, key=lambda row: int(row["N"]))
        reference_case = max(group_cases, key=lambda row: len(row["times"]))
        target_times = list(reference_case["times"])
        n_values = [int(case["N"]) for case in group_cases]
        electrode_names = list(reference_case["groups"].keys())
        target_q_checks = [None, 24]

        for electrode in electrode_names:
            fig = plt.figure(figsize=(15.5, 6.2))
            axes = [
                fig.add_subplot(1, 2, panel_index + 1, projection="3d")
                for panel_index in range(len(target_q_checks))
            ]

            for panel_index, preferred_q in enumerate(target_q_checks):
                ref_error_matrix = []
                ref_gap_matrix = []
                q_check = None
                q_reference = None

                for case in group_cases:
                    group = case["groups"].get(electrode)
                    if group is None:
                        ref_error_matrix.append([0.0 for _ in target_times])
                        ref_gap_matrix.append([0.0 for _ in target_times])
                        continue

                    if group["errs"]:
                        ref_error_values = group["errs"][-1][2]
                    else:
                        ref_error_values = [0.0 for _ in case["times"]]

                    current_q_check, current_q_reference, current_delta = _delta_series_for_q(
                        group, preferred_q
                    )
                    ref_gap_values = current_delta if current_delta else [0.0 for _ in case["times"]]
                    ref_error_matrix.append(
                        _interpolate_series(list(case["times"]), list(ref_error_values), target_times)
                    )
                    ref_gap_matrix.append(
                        _interpolate_series(list(case["times"]), list(ref_gap_values), target_times)
                    )
                    if q_check is None or q_reference is None:
                        q_check, q_reference = current_q_check, current_q_reference

                axis = axes[panel_index]
                primary_surface, _secondary_surface = _plot_dual_error_surfaces(
                    axis,
                    x_values=target_times,
                    y_values=n_values,
                    primary_matrix=ref_error_matrix,
                    secondary_matrix=ref_gap_matrix,
                    title=(
                        f"{electrode}: numerical/reference(qReference) error "
                        f"with {_quadrature_difference_label(q_check, q_reference)}"
                    ),
                    ylabel="N",
                    primary_label="numerical/reference(qReference) error",
                    secondary_label=_quadrature_difference_label(q_check, q_reference),
                )
                fig.colorbar(
                    primary_surface,
                    ax=axis,
                    pad=0.02,
                    shrink=0.70,
                    label="shared error heatmap",
                )

            fig.suptitle(
                f"Manufactured pseudoECG sweep surfaces ({dimension}, {solver}, {electrode}, log z-scale)",
                fontsize=12,
            )
            fig.subplots_adjust(left=0.03, right=0.94, bottom=0.08, top=0.90, wspace=0.10)

            save_path = None
            if save_dir is not None:
                save_path = (
                    Path(save_dir)
                    / f"manufactured_ecg_sweep_surface_{dimension.lower()}_{solver.lower()}_{electrode}.png"
                )
                output_paths.append(save_path)
            finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)

    return output_paths


def export_ecg_error_sweep_vtp(
    cases,
    *,
    save_dir: str | Path | None = None,
):
    if save_dir is None or not cases:
        return []

    output_paths: list[Path] = []
    grouped_cases = group_ecg_timeseries_cases(cases)

    for (dimension, solver), group_cases in grouped_cases.items():
        group_cases = sorted(group_cases, key=lambda row: int(row["N"]))
        reference_case = max(group_cases, key=lambda row: len(row["times"]))
        target_times = list(reference_case["times"])
        n_values = [int(case["N"]) for case in group_cases]
        electrode_names = list(reference_case["groups"].keys())

        for electrode in electrode_names:
            ref_error_matrix = []
            ref_gap_matrix = []
            q_check = None
            q_reference = None

            for case in group_cases:
                group = case["groups"].get(electrode)
                if group is None:
                    ref_error_matrix.append([0.0 for _ in target_times])
                    ref_gap_matrix.append([0.0 for _ in target_times])
                    continue

                if group["errs"]:
                    ref_error_values = group["errs"][-1][2]
                else:
                    ref_error_values = [0.0 for _ in case["times"]]

                current_q_check, current_q_reference, current_delta = _primary_delta_series(group)
                ref_gap_values = current_delta if current_delta else [0.0 for _ in case["times"]]
                ref_error_matrix.append(
                    _interpolate_series(list(case["times"]), list(ref_error_values), target_times)
                )
                ref_gap_matrix.append(
                    _interpolate_series(list(case["times"]), list(ref_gap_values), target_times)
                )
                if q_check is None or q_reference is None:
                    q_check, q_reference = current_q_check, current_q_reference

            base_path = (
                Path(save_dir)
                / f"manufactured_ecg_sweep_surface_{dimension.lower()}_{solver.lower()}_{electrode}"
            )
            output_paths.extend(
                _export_dual_surface_polydata(
                    base_path=base_path,
                    x_values=target_times,
                    y_values=n_values,
                    primary_matrix=ref_error_matrix,
                    secondary_matrix=ref_gap_matrix,
                    title_prefix=f"Manufactured pseudoECG sweep surface ({dimension}, {solver}, {electrode})",
                    primary_slug="numerical_reference_error",
                    secondary_slug=_slugify_token(_quadrature_difference_label(q_check, q_reference)),
                )
            )

    return output_paths


def plot_ecg_error_heatmap_sweep_overview(
    cases,
    *,
    save_dir: str | Path | None = None,
    show: bool = True,
):
    if not _has_surface_plotting():
        print("matplotlib/numpy is not available; skipping ECG sweep surface overviews.")
        return []
    if not cases:
        return []

    configure_matplotlib_defaults()
    output_paths: list[Path] = []
    grouped_cases = group_ecg_timeseries_cases(cases)

    for (dimension, solver), group_cases in grouped_cases.items():
        group_cases = sorted(group_cases, key=lambda row: int(row["N"]))
        reference_case = max(group_cases, key=lambda row: len(row["times"]))
        target_times = list(reference_case["times"])
        n_values = [int(case["N"]) for case in group_cases]
        electrode_names = list(reference_case["groups"].keys())
        target_q_checks = [None, 24]

        nrows = len(electrode_names)
        fig = plt.figure(figsize=(20.0, max(4.2 * nrows, 8.0)))
        axes: list[object] = []
        for row_index in range(nrows):
            for panel_index in range(len(target_q_checks)):
                axes.append(
                    fig.add_subplot(
                        nrows,
                        len(target_q_checks),
                        row_index * len(target_q_checks) + panel_index + 1,
                        projection="3d",
                    )
                )

        for row_index, electrode in enumerate(electrode_names):
            for panel_index, preferred_q in enumerate(target_q_checks):
                ref_error_matrix = []
                ref_gap_matrix = []
                q_check = None
                q_reference = None

                for case in group_cases:
                    group = case["groups"].get(electrode)
                    if group is None:
                        ref_error_matrix.append([0.0 for _ in target_times])
                        ref_gap_matrix.append([0.0 for _ in target_times])
                        continue

                    if group["errs"]:
                        ref_error_values = group["errs"][-1][2]
                    else:
                        ref_error_values = [0.0 for _ in case["times"]]

                    current_q_check, current_q_reference, current_delta = _delta_series_for_q(
                        group, preferred_q
                    )
                    ref_gap_values = current_delta if current_delta else [0.0 for _ in case["times"]]
                    ref_error_matrix.append(
                        _interpolate_series(list(case["times"]), list(ref_error_values), target_times)
                    )
                    ref_gap_matrix.append(
                        _interpolate_series(list(case["times"]), list(ref_gap_values), target_times)
                    )
                    if q_check is None or q_reference is None:
                        q_check, q_reference = current_q_check, current_q_reference

                axis = axes[row_index * len(target_q_checks) + panel_index]
                primary_surface, _secondary_surface = _plot_dual_error_surfaces(
                    axis,
                    x_values=target_times,
                    y_values=n_values,
                    primary_matrix=ref_error_matrix,
                    secondary_matrix=ref_gap_matrix,
                    title=(
                        f"{electrode}: numerical/reference(qReference) error "
                        f"with {_quadrature_difference_label(q_check, q_reference)}"
                    ),
                    ylabel="N",
                    primary_label="numerical/reference(qReference) error",
                    secondary_label=_quadrature_difference_label(q_check, q_reference),
                )
                fig.colorbar(
                    primary_surface,
                    ax=axis,
                    pad=0.02,
                    shrink=0.62,
                    label="shared error heatmap",
                )

        fig.suptitle(
            f"Manufactured pseudoECG sweep surfaces for all electrodes ({dimension}, {solver}, log z-scale)",
            fontsize=12,
        )
        fig.subplots_adjust(left=0.03, right=0.95, bottom=0.04, top=0.96, hspace=0.30, wspace=0.08)

        save_path = None
        if save_dir is not None:
            save_path = (
                Path(save_dir)
                / f"manufactured_ecg_sweep_surface_{dimension.lower()}_{solver.lower()}_all_electrodes.png"
            )
            output_paths.append(save_path)
        finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)

    return output_paths


def plot_ecg_electrode_geometry(
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_matplotlib():
        print("matplotlib is not available; skipping ECG electrode geometry plot.")
        return None

    configure_matplotlib_defaults()
    fig = plt.figure(figsize=(7.8, 6.2))
    ax3 = fig.add_subplot(1, 1, 1, projection="3d")

    three_d = driver_defaults.ECG_ELECTRODES_BY_DIMENSION["3D"]
    for start, end in _unit_cube_edges():
        ax3.plot(
            [start[0], end[0]],
            [start[1], end[1]],
            [start[2], end[2]],
            color="black",
            linewidth=1.2,
        )
    for name, literal in three_d.items():
        x, y, z = _parse_vector_literal(literal)
        ax3.scatter([x], [y], [z], s=45)
        ax3.text(x, y, z + 0.05, name, fontsize=8)
    ax3.set_title("3D electrodes vs unit cube")
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.set_zlabel("z")
    ax3.grid(True, linestyle="--", alpha=0.35)
    ax3.set_box_aspect((1.6, 1.4, 1.4))

    fig.suptitle("Manufactured pseudoECG 3D electrode locations", fontsize=13)
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.05, top=0.90)
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def plot_errors(rows, solver_type=None, *, save_path: str | Path | None = None, show: bool = True):
    """
    Plot Linf errors for Vm, u1, u2 vs N.
    """
    if not _has_matplotlib():
        print("matplotlib is not available; skipping manufactured error plot.")
        return None

    configure_matplotlib_defaults()
    fig, ax = plt.subplots(figsize=(8, 5))

    for dimension in _unique_values(rows, "Dimension"):
        dimension_rows = _filter_rows(rows, Dimension=dimension)

        if solver_type:
            dimension_rows = _filter_rows(dimension_rows, Solver=solver_type)
        if not dimension_rows:
            continue

        ax.loglog(
            [row["N"] for row in dimension_rows],
            [row["Linf_V"] for row in dimension_rows],
            marker="o",
            label=f"Vm ({dimension})",
        )
        ax.loglog(
            [row["N"] for row in dimension_rows],
            [row["Linf_u1"] for row in dimension_rows],
            marker="s",
            label=f"u1 ({dimension})",
        )
        ax.loglog(
            [row["N"] for row in dimension_rows],
            [row["Linf_u2"] for row in dimension_rows],
            marker="^",
            label=f"u2 ({dimension})",
        )

    title = "Manufactured-solution Linf errors"
    if solver_type:
        title += f" ({solver_type})"
    style_matplotlib_axes(
        ax,
        title=title,
        xlabel="Number of cells (N)",
        ylabel="Linf Error",
        grid_kwargs={"which": "both", "ls": "--"},
    )
    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def plot_errors_implicit_explicit(
    rows,
    *,
    save_dir: str | Path | None = None,
    show: bool = True,
):
    """
    Plot Linf errors for Vm, u1, u2 vs N for both Explicit and Implicit solvers.
    Grouped per dimension.
    """
    if not _has_matplotlib():
        print("matplotlib is not available; skipping manufactured implicit/explicit plots.")
        return []

    configure_matplotlib_defaults()

    output_paths: list[Path] = []
    for dimension in _unique_values(rows, "Dimension"):
        fig, ax = plt.subplots(figsize=(8, 6))
        _plot_dimension_errors_on_axis(ax, rows, dimension)
        save_path = None
        if save_dir is not None:
            save_path = Path(save_dir) / f"manufactured_errors_{dimension.lower()}_implicit_explicit.png"
            output_paths.append(save_path)
        finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return output_paths


def plot_Vm_across_dimensions(
    rows,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    """
    Plot Linf_V (Vm error) vs N across all dimensions and solvers.
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


def plot_ecg_error_vs_gap_summary(
    rows,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
):
    if not _has_matplotlib():
        print("matplotlib is not available; skipping ECG error-vs-gap summary.")
        return None

    configure_matplotlib_defaults()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.0), squeeze=False)
    flat_axes = axes.flatten()
    q_check = None
    q_reference = None
    if rows:
        q_checks = sorted({int(row["qCheck"]) for row in rows if int(row["qCheck"]) > 0})
        q_references = sorted({int(row["qReference"]) for row in rows if int(row["qReference"]) > 0})
        q_check = q_checks[0] if q_checks else None
        q_reference = q_references[0] if q_references else None
    reference_error_label = (
        f"|numeric - reference(q={q_reference})|"
        if q_reference is not None
        else "|numeric - reference(qReference)|"
    )
    quadrature_label = _quadrature_difference_label(q_check, q_reference)
    fig.suptitle(
        "Manufactured pseudoECG (3D): numerical reference error vs quadrature difference",
        fontsize=13,
    )

    metric_specs = (
        (
            "max_Linf_err_ref",
            "max_Linf_delta_ref",
            "3D max metrics",
            "max time-Linf error magnitude",
            f"max {reference_error_label}",
            f"max {quadrature_label}",
        ),
        (
            "mean_Linf_err_ref",
            "mean_Linf_delta_ref",
            "3D mean metrics",
            "mean time-Linf error magnitude",
            f"mean {reference_error_label}",
            f"mean {quadrature_label}",
        ),
    )

    for idx, (error_key, gap_key, title, ylabel, error_label, gap_label) in enumerate(metric_specs):
        ax = flat_axes[idx]
        if not rows:
            ax.text(
                0.5,
                0.5,
                "No 3D ECG data",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=11,
                color="0.4",
            )
            style_matplotlib_axes(
                ax,
                title=title,
                xlabel="Number of cells (N)",
                ylabel=ylabel,
                legend=False,
                grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
            )
            continue

        for solver in _unique_values(rows, "Solver"):
            solver_rows = _filter_rows(rows, Solver=solver)
            ax.loglog(
                [row["N"] for row in solver_rows],
                [max(row[error_key], 1e-30) for row in solver_rows],
                marker="o",
                linestyle="-",
                label=f"{error_label} ({solver})",
            )
            ax.loglog(
                [row["N"] for row in solver_rows],
                [max(row[gap_key], 1e-30) for row in solver_rows],
                marker="s",
                linestyle="--",
                label=f"{gap_label} ({solver})",
            )

        style_matplotlib_axes(
            ax,
            title=title,
            xlabel="Number of cells (N)",
            ylabel=ylabel,
            legend=True,
            grid_kwargs={"which": "both", "ls": "--", "alpha": 0.45},
        )

    finalize_matplotlib_figure(fig, save_path=save_path, show=show, close=not show)
    return Path(save_path) if save_path is not None else None


def write_ecg_plot_guide(
    output_dir: Path,
    *,
    q_checks: str,
    q_reference: int | None,
) -> Path:
    content = f"""# Manufactured pseudoECG plot guide

This folder contains both monodomain manufactured-solution plots and pseudoECG manufactured-verification plots.

## ECG dimension scope

ECG post-processing is currently generated only for `3D`.

If `1D` or `2D` ECG artifacts are present, they are intentionally ignored with a warning because:
- the numerical pseudoECG in OpenFOAM is accumulated as a 3D cell-volume sum
- the current `1D` and `2D` manufactured ECG references are lower-dimensional integrals

So `1D/2D` remain valid for the manufactured monodomain field verification, but they are not used for ECG verification plots.
During post-processing, archived ECG case files for unsupported dimensions are removed from this folder to keep the ECG artifact set consistent with the 3D-only verification path.

## ECG sweep summary plots

- `manufactured_ecg_electrode_geometry.png`
  Shows the configured 3D electrode locations relative to the unit cube only.

- `manufactured_ecg_electrode_geometry.vtp`
  VTK XML PolyData export of the 3D unit cube and ECG electrode points.
  The file also includes axis lines built from the geometry bounds.
  The point coordinates are normalized to a unit viewing box, and the original coordinates are stored as point-data arrays `raw_x`, `raw_y`, and `raw_z`.
  Open this in ParaView if you want to choose your own camera for the geometry view.

- `manufactured_ecg_error_vs_gap_summary.png`
  Two sweep-level 3D-only summary panels versus `N`:
  1. max over electrodes of time-`Linf` numerical/reference error and quadrature difference
  2. mean over electrodes of time-`Linf` numerical/reference error and quadrature difference

  This is the main sweep-level ECG summary plot.
  The numerical/reference error should decrease with refinement, while the quadrature difference should stay much smaller.

- `manufactured_ecg_quadrature_only_summary.png`
  Two quadrature-only sweep panels versus `N`:
  1. max over electrodes of time-`Linf` quadrature difference for each `qCheck`
  2. mean over electrodes of time-`Linf` quadrature difference for each `qCheck`

  This isolates the reference-integral convergence without the numerical pseudoECG error.

- `manufactured_ecg_quadrature_final_time_heatmap.png`
  2D heatmap of the final-time maximum quadrature difference over electrodes:
  - x-axis: `N`
  - y-axis: `qCheck`
  - value: `max_electrodes |reference(qCheck) - reference(qReference)|` at the final time

  This is the direct final-time `N x qCheck` view of the quadrature-only error.

## ECG representative time-series plots

For each dimension/solver pair, the postprocess chooses the highest available `N` case and produces:

- `manufactured_ecg_overlay_<dimension>_<solver>.png`
  Per-electrode time traces of:
  - `numeric`
  - `reference(q=qReference)`

- `manufactured_ecg_quadrature_overlay_<dimension>_<solver>.png`
  Per-electrode time traces of the manufactured reference for all available quadrature orders:
  - `reference(q=6)`
  - `reference(q=12)`
  - `reference(q=24)`
  - `reference(q=48)`
  - `reference(q=qReference)`

  This is the direct “quadratures only” comparison at the electrode points.

- `manufactured_ecg_reference_error_timeseries_<dimension>_<solver>.png`
  Per-electrode semilogy time traces of:
  - `abs error to reference(q=qReference) = |numeric - reference(q=qReference)|`
  - `quadrature difference = |reference(q=qCheck) - reference(q=qReference)|`

  In the raw OpenFOAM columns:
  - `errQ96_<electrode>` means the absolute error to the manufactured reference evaluated with quadrature order 96
  - `deltaQuadratureQ6_Q96_<electrode> = |refQ6_<electrode> - refQ96_<electrode>|`
  - more generally, `errQk_<electrode> = |numeric_<electrode> - refQk_<electrode>|`

- `manufactured_ecg_reference_error_surface_<dimension>_<solver>.png`
  One 3D plot for the representative case with two overlaid colormapped surfaces:
  - viridis surface with black mesh outline: `|numeric - refQreference|`
  - viridis surface with white mesh outline: `|refQcheck_primary - refQreference|`

  Surface axes:
  - x-axis: time
  - y-axis: electrode index
  - z-axis: error magnitude
  - z-axis uses a log scale so both surfaces remain visible when their magnitudes differ
  - both surfaces use one shared heatmap/colorbar, and the legend distinguishes them by outline

- `manufactured_ecg_reference_error_surface_<dimension>_<solver>_*.vtp`
  VTK XML PolyData exports of the representative 3D surfaces.
  These files include the surface triangles plus bounding-box and axis line geometry built from the surface bounds.
  The exported point coordinates are normalized to a unit viewing box, while the original coordinates are stored as point-data arrays `raw_x`, `raw_y`, and `raw_z`.
  Open these in ParaView (or another VTK viewer) if you want to choose your own camera and render the figure there.

- `manufactured_ecg_sweep_surface_<dimension>_<solver>_<electrode>.png`
  Two side-by-side 3D sweep plots for one electrode:
  - left subplot uses the primary check quadrature difference
  - right subplot uses the `q=24` quadrature difference

  In each subplot:
  - x-axis: time
  - y-axis: mesh size `N`
  - z-axis: error magnitude
  - viridis surface with black mesh outline: `|numeric - refQreference|`
  - viridis surface with white mesh outline: `|refQcheck_selected - refQreference|`

  The z-axis uses a log scale so the numerical/reference error and the quadrature gap can be read on the same plot without the smaller surface collapsing to zero visually.

  These are the plots to inspect if you want to see how the numerical reference error evolves jointly in time and mesh resolution for a fixed electrode.

- `manufactured_ecg_sweep_surface_<dimension>_<solver>_<electrode>_*.vtp`
  VTK XML PolyData exports of the per-electrode sweep surfaces.
  These files include the surface triangles plus bounding-box and axis line geometry built from the sweep bounds.
  The exported point coordinates are normalized to a unit viewing box, while the original coordinates are stored as point-data arrays `raw_x`, `raw_y`, and `raw_z`.
  These use the raw surface values, not the Matplotlib camera or the PNG view.

- `manufactured_ecg_sweep_surface_<dimension>_<solver>_all_electrodes.png`
  Multi-panel overview with all electrodes in one figure.
  Each row is one electrode and contains two 3D plots:
  - left column: primary check quadrature difference
  - right column: `q=24` quadrature difference

  In each subplot:
  - viridis surface with black mesh outline: `|numeric - refQreference|`
  - viridis surface with white mesh outline: `|refQcheck_selected - refQreference|`

## Configured quadrature orders

- `qChecks = {q_checks if q_checks else "unknown"}`
- `qReference = {q_reference if q_reference is not None else "unknown"}`

All ECG plots in this folder are derived only from OpenFOAM outputs built from those configured quadrature orders.
"""
    destination = output_dir / "manufactured_ecg_plot_guide.md"
    destination.write_text(content, encoding="utf-8")
    return destination


def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **_: object) -> None:
    del setup_root
    output_path = Path(output_dir)
    _cleanup_unsupported_ecg_archives(output_path)
    _cleanup_stale_ecg_plot_artifacts(output_path)
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
    error_plots = plot_errors_implicit_explicit(
        error_rows,
        save_dir=output_path,
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
    artifacts.extend(
        {
            "path": str(path),
            "label": f"Manufactured errors {path.stem}",
            "kind": "plot",
            "format": "png",
        }
        for path in error_plots
    )

    ecg_rows = _filter_supported_ecg_rows(
        read_ecg_summary_dat_files(output_dir),
        source_label="summary",
    )
    if not ecg_rows:
        print(f"No ECG manufactured summary files found to post-process in: {output_dir}")
        return artifacts

    print("\nECG convergence rates:")
    ecg_convergence_rates = compute_ecg_convergence_rates(ecg_rows)
    if ecg_convergence_rates:
        for row in ecg_convergence_rates:
            print(row)
    else:
        print("No ECG convergence-rate pairs were found.")

    ecg_summary_csv = output_path / "manufactured_ecg_convergence_summary.csv"
    _write_csv(ecg_rows, ecg_summary_csv, ECG_SUMMARY_FIELDS)
    artifacts.append(
        {
            "path": str(ecg_summary_csv),
            "label": "Manufactured ECG convergence summary",
            "kind": "table",
            "format": "csv",
        }
    )

    ecg_rates_csv = output_path / "manufactured_ecg_convergence_rates.csv"
    _write_csv(ecg_convergence_rates, ecg_rates_csv, ECG_RATE_FIELDS)
    artifacts.append(
        {
            "path": str(ecg_rates_csv),
            "label": "Manufactured ECG convergence rates",
            "kind": "table",
            "format": "csv",
        }
    )

    if not _has_matplotlib():
        print("matplotlib is not available; generated ECG CSV files only.")
        return artifacts

    ecg_geometry_plot = plot_ecg_electrode_geometry(
        save_path=output_path / "manufactured_ecg_electrode_geometry.png",
        show=False,
    )
    ecg_geometry_vtp = export_ecg_electrode_geometry_vtp(
        save_path=output_path / "manufactured_ecg_electrode_geometry.vtp",
    )
    ecg_error_vs_gap = plot_ecg_error_vs_gap_summary(
        ecg_rows,
        save_path=output_path / "manufactured_ecg_error_vs_gap_summary.png",
        show=False,
    )
    ecg_quadrature_only_summary = plot_ecg_quadrature_summary(
        ecg_rows,
        save_path=output_path / "manufactured_ecg_quadrature_only_summary.png",
        show=False,
    )
    q_check = None
    q_checks_display = ""
    q_reference = None
    if ecg_rows:
        q_checks = sorted({int(row["qCheck"]) for row in ecg_rows if int(row["qCheck"]) > 0})
        q_checks_tokens = sorted(
            {
                token
                for row in ecg_rows
                for token in str(row.get("qChecks", "")).split()
                if token.strip()
            },
            key=int,
        )
        q_references = sorted({int(row["qReference"]) for row in ecg_rows if int(row["qReference"]) > 0})
        q_check = q_checks[0] if q_checks else None
        q_checks_display = " ".join(q_checks_tokens)
        q_reference = q_references[0] if q_references else None
    for path, label in (
        (ecg_geometry_plot, "Manufactured ECG electrode geometry"),
        (ecg_error_vs_gap, "Manufactured ECG 3D error vs gap summary"),
        (
            ecg_quadrature_only_summary,
            "Manufactured ECG quadrature-only summary",
        ),
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
    if ecg_geometry_vtp is not None:
        artifacts.append(
            {
                "path": str(ecg_geometry_vtp),
                "label": "Manufactured ECG electrode geometry VTP",
                "kind": "data",
                "format": "vtp",
            }
        )

    ecg_timeseries_cases = _filter_supported_ecg_rows(
        read_ecg_timeseries_dat_files(output_dir),
        source_label="time-series",
    )
    ecg_quadrature_final_time_heatmap = plot_ecg_quadrature_final_time_heatmap(
        ecg_timeseries_cases,
        save_path=output_path / "manufactured_ecg_quadrature_final_time_heatmap.png",
        show=False,
    )
    representative_cases = select_representative_ecg_cases(ecg_timeseries_cases)
    sweep_surface_vtp = export_ecg_error_sweep_vtp(
        ecg_timeseries_cases,
        save_dir=output_path,
    )
    sweep_surfaces = plot_ecg_error_heatmap_sweep(
        ecg_timeseries_cases,
        save_dir=output_path,
        show=False,
    )
    sweep_surface_overviews = plot_ecg_error_heatmap_sweep_overview(
        ecg_timeseries_cases,
        save_dir=output_path,
        show=False,
    )
    for case in representative_cases:
        dimension_token = str(case["Dimension"]).lower()
        solver_token = str(case["Solver"]).lower()
        surface_vtp_paths = export_ecg_error_surface_case_vtp(
            case,
            save_dir=output_path,
        )
        overlay_path = plot_ecg_reference_overlay_case(
            case,
            save_path=output_path / f"manufactured_ecg_overlay_{dimension_token}_{solver_token}.png",
            show=False,
        )
        quadrature_overlay_path = plot_ecg_quadrature_reference_overlay_case(
            case,
            save_path=output_path / f"manufactured_ecg_quadrature_overlay_{dimension_token}_{solver_token}.png",
            show=False,
        )
        error_time_path = plot_ecg_error_timeseries_case(
            case,
            save_path=output_path / f"manufactured_ecg_reference_error_timeseries_{dimension_token}_{solver_token}.png",
            show=False,
        )
        surface_path = plot_ecg_error_heatmap_case(
            case,
            save_path=output_path / f"manufactured_ecg_reference_error_surface_{dimension_token}_{solver_token}.png",
            show=False,
        )

        for path, label in (
            (overlay_path, f"Manufactured ECG overlay {case['Dimension']} {case['Solver']}"),
            (
                quadrature_overlay_path,
                f"Manufactured ECG quadrature overlay {case['Dimension']} {case['Solver']}",
            ),
            (
                error_time_path,
                f"Manufactured ECG reference errors vs time {case['Dimension']} {case['Solver']}",
            ),
            (
                surface_path,
                f"Manufactured ECG reference error surface {case['Dimension']} {case['Solver']}",
            ),
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
                "label": f"Manufactured ECG VTP surface {path.stem}",
                "kind": "data",
                "format": "vtp",
            }
            for path in surface_vtp_paths
        )
    artifacts.extend(
        {
            "path": str(path),
            "label": f"Manufactured ECG sweep surface {path.stem}",
            "kind": "plot",
            "format": "png",
        }
        for path in sweep_surfaces
    )
    artifacts.extend(
        {
            "path": str(path),
            "label": f"Manufactured ECG sweep surface overview {path.stem}",
            "kind": "plot",
            "format": "png",
        }
        for path in sweep_surface_overviews
    )
    if ecg_quadrature_final_time_heatmap is not None:
        artifacts.append(
            {
                "path": str(ecg_quadrature_final_time_heatmap),
                "label": "Manufactured ECG final-time quadrature heatmap",
                "kind": "plot",
                "format": "png",
            }
        )
    artifacts.extend(
        {
            "path": str(path),
            "label": f"Manufactured ECG VTP sweep surface {path.stem}",
            "kind": "data",
            "format": "vtp",
        }
        for path in sweep_surface_vtp
    )
    guide_path = write_ecg_plot_guide(
        output_path,
        q_checks=q_checks_display,
        q_reference=q_reference,
    )
    artifacts.append(
        {
            "path": str(guide_path),
            "label": "Manufactured ECG plot guide",
            "kind": "report",
            "format": "md",
        }
    )
    return artifacts
