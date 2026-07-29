"""table_summary.py — CV convergence tables and plots for monodomain1DCableCV."""
from __future__ import annotations

import csv
import io
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

TUTORIALS_ROOT = Path(__file__).resolve().parents[5]
DRIVER_ROOT = TUTORIALS_ROOT / "applications" / "scripts" / "driverFoam"
if str(DRIVER_ROOT) not in sys.path:
    sys.path.insert(0, str(DRIVER_ROOT))

try:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

from openfoam_driver.postprocessing.style import (
    configure_matplotlib_defaults,
    finalize_matplotlib_figure,
    style_matplotlib_axes,
)
from openfoam_driver.postprocessing.table_writer import TableMetadata


def _parse_filename(stem: str) -> dict[str, object]:
    match = re.match(
        r"(?P<solver>[^_]+)_(?P<model>[^_]+)_(?P<tissue>.+)_DT(?P<dt>[0-9.]+)_DX(?P<dx>[0-9.]+)_COND(?P<cond>\d+)_cv_summary",
        stem,
    )
    if not match:
        raise ValueError(f"Unexpected CV summary filename: {stem}")
    return {
        "solver": match.group("solver"),
        "ionic_model": match.group("model"),
        "tissue": match.group("tissue"),
        "DT_ms": float(match.group("dt")),
        "DX_mm": float(match.group("dx")),
        "conductivity_id": int(match.group("cond")),
    }


def build_summary_rows(output_dir: Path) -> list[dict]:
    files = sorted(output_dir.glob("*_cv_summary.json"))
    if not files:
        print(f"[monodomainAndEikonal1DCableCVConvergence/table_summary] No *_cv_summary.json files in {output_dir}")
        return []

    rows: list[dict] = []
    grouped: dict[tuple[str, str, str, int], list[dict]] = defaultdict(list)

    for fpath in files:
        payload = json.loads(fpath.read_text())
        name_info = _parse_filename(fpath.stem)
        central = payload["central_cv"]
        row = {
            "case_id": payload["case_id"],
            **name_info,
            "central_dx_mm": round(1e3 * float(central["dx_m"]), 6),
            "central_dt_ms": round(1e3 * float(central["dt_s"]), 6),
            "central_cv_m_per_s": round(float(central["cv_m_per_s"]), 8),
        }
        rows.append(row)
        grouped[
            (
                str(row["solver"]),
                str(row["ionic_model"]),
                str(row["tissue"]),
                int(row["conductivity_id"]),
            )
        ].append(row)

    for group_rows in grouped.values():
        reference = min(
            group_rows,
            key=lambda row: (float(row["DX_mm"]), float(row["DT_ms"])),
        )
        reference_cv = float(reference["central_cv_m_per_s"])
        for row in group_rows:
            row["reference_case_id"] = reference["case_id"]
            row["reference_cv_m_per_s"] = round(reference_cv, 8)
            row["abs_error_m_per_s"] = round(
                abs(float(row["central_cv_m_per_s"]) - reference_cv),
                8,
            )
            row["rel_error_percent"] = (
                round(
                    100.0 * abs(float(row["central_cv_m_per_s"]) - reference_cv) / reference_cv,
                    6,
                )
                if reference_cv > 0
                else 0.0
            )

    rows.sort(
        key=lambda row: (
            int(row["conductivity_id"]),
            str(row["solver"]),
            str(row["ionic_model"]),
            str(row["tissue"]),
            float(row["DX_mm"]),
            float(row["DT_ms"]),
        )
    )
    return rows


def _discover_summary_dirs(output_root: Path) -> list[tuple[Path, Path]]:
    if not output_root.exists():
        print(f"[monodomainAndEikonal1DCableCVConvergence/table_summary] Output directory does not exist: {output_root}")
        return []

    files = sorted(output_root.rglob("*_cv_summary.json"))
    if not files:
        print(f"[monodomainAndEikonal1DCableCVConvergence/table_summary] No *_cv_summary.json files in {output_root}")
        return []

    parents = sorted({file_path.parent for file_path in files})
    return [(parent, parent.relative_to(output_root)) for parent in parents]


def _write_summary_csv(
    rows: list[dict],
    *,
    output_dir: Path,
    filename_stem: str,
    metadata: TableMetadata,
    label: str,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{filename_stem}.csv"
    html_path = output_dir / f"{filename_stem}.html"

    fieldnames = list(rows[0].keys()) if rows else []
    lines = [
        f"# entry: {metadata.entry}",
        f"# generated_at: {metadata.generated_at}",
        f"# units: {json.dumps(metadata.units)}",
    ]
    if fieldnames:
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        lines.extend(buffer.getvalue().strip().splitlines())
    csv_path.write_text("\n".join(lines) + "\n")

    if html_path.exists():
        html_path.unlink()

    return {
        "path": csv_path.name,
        "label": label,
        "kind": "table",
        "format": "csv",
    }


def _rebase_artifacts(
    artifacts: list[dict[str, object]],
    *,
    relative_dir: Path,
) -> list[dict[str, object]]:
    if relative_dir == Path("."):
        return artifacts

    rebased: list[dict[str, object]] = []
    for artifact in artifacts:
        updated = dict(artifact)
        updated["path"] = str(relative_dir / str(artifact["path"]))
        rebased.append(updated)
    return rebased


def _group_key(row: dict[str, object]) -> tuple[str, str, str, int]:
    return (
        str(row["solver"]),
        str(row["ionic_model"]),
        str(row["tissue"]),
        int(row["conductivity_id"]),
    )


def _group_tag(group: tuple[str, str, str, int]) -> str:
    solver, ionic_model, tissue, conductivity_id = group
    tissue_slug = re.sub(r"[^A-Za-z0-9]+", "-", tissue).strip("-")
    return f"{solver}_{ionic_model}_{tissue_slug}_COND{conductivity_id:02d}"


def _group_title(group: tuple[str, str, str, int]) -> str:
    solver, ionic_model, tissue, conductivity_id = group
    return (
        f"{ionic_model} | {tissue} | {solver} | conductivity #{conductivity_id:02d}"
    )


def _has_matplotlib() -> bool:
    return plt is not None


def _positive_errors(values: list[tuple[float, float]]) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    for x_value, y_value in values:
        if y_value > 0.0:
            xs.append(x_value)
            ys.append(y_value)
    return xs, ys


def _plot_group_convergence(
    group: tuple[str, str, str, int],
    group_rows: list[dict],
    *,
    output_dir: Path,
) -> list[dict[str, object]]:
    if not _has_matplotlib():
        print("[monodomainAndEikonal1DCableCVConvergence/table_summary] matplotlib unavailable; writing CSV only.")
        return []

    configure_matplotlib_defaults()
    tag = _group_tag(group)
    title = _group_title(group)
    artifacts: list[dict[str, object]] = []

    reference_row = min(group_rows, key=lambda row: (float(row["DX_mm"]), float(row["DT_ms"])))
    reference_cv = float(reference_row["reference_cv_m_per_s"])

    dt_groups: dict[float, list[dict]] = defaultdict(list)
    dx_groups: dict[float, list[dict]] = defaultdict(list)
    for row in group_rows:
        dt_groups[float(row["DT_ms"])].append(row)
        dx_groups[float(row["DX_mm"])].append(row)

    fig_dx, axes_dx = plt.subplots(1, 2, figsize=(13.5, 4.8))
    for dt_value in sorted(dt_groups):
        series = sorted(dt_groups[dt_value], key=lambda row: float(row["DX_mm"]))
        dx_values = [float(row["DX_mm"]) for row in series]
        cv_values = [float(row["central_cv_m_per_s"]) for row in series]
        err_pairs = [
            (float(row["DX_mm"]), float(row["abs_error_m_per_s"]))
            for row in series
        ]
        axes_dx[0].plot(dx_values, cv_values, marker="o", linewidth=1.5, label=f"dt={dt_value:g} ms")
        err_x, err_y = _positive_errors(err_pairs)
        if err_x:
            axes_dx[1].plot(err_x, err_y, marker="o", linewidth=1.5, label=f"dt={dt_value:g} ms")

    axes_dx[0].axhline(reference_cv, color="black", linestyle="--", linewidth=1.0, label="reference CV")
    axes_dx[1].set_xscale("log")
    axes_dx[1].set_yscale("log")
    style_matplotlib_axes(
        axes_dx[0],
        title="Central CV vs DX",
        xlabel="DX (mm)",
        ylabel="CV (m/s)",
        legend=True,
    )
    style_matplotlib_axes(
        axes_dx[1],
        title="Absolute CV error vs DX",
        xlabel="DX (mm)",
        ylabel="|CV - CV_ref| (m/s)",
        legend=True,
    )
    fig_dx.suptitle(f"1D cable CV mesh convergence\n{title}", fontsize=12)
    dx_path = output_dir / f"monodomainAndEikonal1DCableCVConvergence_dx_{tag}.png"
    finalize_matplotlib_figure(fig_dx, save_path=dx_path, show=False, close=True)
    artifacts.append(
        {
            "path": dx_path.name,
            "label": f"1D cable CV mesh convergence ({title})",
            "kind": "plot",
            "format": "png",
        }
    )

    fig_dt, axes_dt = plt.subplots(1, 2, figsize=(13.5, 4.8))
    for dx_value in sorted(dx_groups):
        series = sorted(dx_groups[dx_value], key=lambda row: float(row["DT_ms"]))
        dt_values = [float(row["DT_ms"]) for row in series]
        cv_values = [float(row["central_cv_m_per_s"]) for row in series]
        err_pairs = [
            (float(row["DT_ms"]), float(row["abs_error_m_per_s"]))
            for row in series
        ]
        axes_dt[0].plot(dt_values, cv_values, marker="o", linewidth=1.5, label=f"dx={dx_value:g} mm")
        err_x, err_y = _positive_errors(err_pairs)
        if err_x:
            axes_dt[1].plot(err_x, err_y, marker="o", linewidth=1.5, label=f"dx={dx_value:g} mm")

    axes_dt[0].axhline(reference_cv, color="black", linestyle="--", linewidth=1.0, label="reference CV")
    axes_dt[1].set_xscale("log")
    axes_dt[1].set_yscale("log")
    style_matplotlib_axes(
        axes_dt[0],
        title="Central CV vs DT",
        xlabel="DT (ms)",
        ylabel="CV (m/s)",
        legend=True,
    )
    style_matplotlib_axes(
        axes_dt[1],
        title="Absolute CV error vs DT",
        xlabel="DT (ms)",
        ylabel="|CV - CV_ref| (m/s)",
        legend=True,
    )
    fig_dt.suptitle(f"1D cable CV time-step convergence\n{title}", fontsize=12)
    dt_path = output_dir / f"monodomainAndEikonal1DCableCVConvergence_dt_{tag}.png"
    finalize_matplotlib_figure(fig_dt, save_path=dt_path, show=False, close=True)
    artifacts.append(
        {
            "path": dt_path.name,
            "label": f"1D cable CV time-step convergence ({title})",
            "kind": "plot",
            "format": "png",
        }
    )

    return artifacts


def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **_: object) -> list[dict]:
    del setup_root
    output_path = Path(output_dir)
    meta = TableMetadata(
        tutorial="monodomainAndEikonal1DCableCVConvergence",
        units={
            "DX_mm": "mm",
            "DT_ms": "ms",
            "central_dx_mm": "mm",
            "central_dt_ms": "ms",
            "central_cv_m_per_s": "m/s",
            "reference_cv_m_per_s": "m/s",
            "abs_error_m_per_s": "m/s",
            "rel_error_percent": "%",
        },
    )
    artifacts: list[dict[str, object]] = []

    for summary_dir, relative_dir in _discover_summary_dirs(output_path):
        rows = build_summary_rows(summary_dir)
        if not rows:
            continue

        dir_artifacts: list[dict[str, object]] = [
            _write_summary_csv(
                rows,
                output_dir=summary_dir,
                filename_stem="monodomainAndEikonal1DCableCVConvergence_summary",
                label="1D cable CV convergence summary",
                metadata=meta,
            )
        ]

        grouped_rows: dict[tuple[str, str, str, int], list[dict]] = defaultdict(list)
        for row in rows:
            grouped_rows[_group_key(row)].append(row)

        for group, rows_for_group in sorted(grouped_rows.items()):
            dir_artifacts.extend(
                _plot_group_convergence(
                    group,
                    rows_for_group,
                    output_dir=summary_dir,
                )
            )

        artifacts.extend(_rebase_artifacts(dir_artifacts, relative_dir=relative_dir))

    return artifacts


if __name__ == "__main__":
    folder = Path(__file__).resolve().parents[1] / "outputsCVConvergence"
    print(f"[table_summary] Default folder = {folder}")
    run_postprocessing(output_dir=str(folder))
