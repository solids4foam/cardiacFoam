"""table_summary.py — CV convergence tables and plots for eikonal1DCableCV."""
from __future__ import annotations

import csv
import io
import json
import re
from collections import defaultdict
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

from omnidriver.postprocessing.style import (
    configure_matplotlib_defaults,
    finalize_matplotlib_figure,
    style_matplotlib_axes,
)
from omnidriver.postprocessing.table_writer import TableMetadata

#: Every study this tutorial sweeps today fixes solver/ionicModel/tissue/
#: conductivity across the whole sweep -- only dx/dt vary -- so there is
#: exactly one convergence group per sweep. conductivity_id stays a constant
#: instead of decoding a COND<n> filename token (nothing generates one any
#: more); extend this if a future study actually sweeps conductivity.
_CONDUCTIVITY_ID = 1


def _dict_scalar(text: str, keyword: str) -> str | None:
    match = re.search(rf"^\s*{re.escape(keyword)}\s+([^;]+);", text, re.MULTILINE)
    return match.group(1).strip() if match else None


def _case_labels(case_dir: Path) -> tuple[str, str, str]:
    """Read (solver, ionic_model, tissue) from this case's own committed dict.

    These used to be decoded from a `solver_model_tissue_..._cv_summary.json`
    filename that driverFoam's sweep wrapper produced by dumping every case's
    output into one shared folder. A staged sweep case owns its own
    electroProperties, so reading it directly is correct whether a value is
    swept or -- as in every study today -- fixed across the whole sweep.
    """
    text = (case_dir / "constant" / "electroProperties").read_text()
    solver = (_dict_scalar(text, "myocardiumSolver") or "unknown").removesuffix("Solver")
    ionic_model = _dict_scalar(text, "ionicModel") or "unknown"
    tissue = _dict_scalar(text, "tissue") or "unknown"
    return solver, ionic_model, tissue


def _iter_cv_summaries(output_dir: Path):
    """Yield (case_dir, resolved_axis_values, payload) per case in this sweep.

    `<output_dir>/cases/case_NNNN/postProcessing/case_NNNN_cv_summary.json`
    replaces the old shared-folder convention; a case's own
    case_record.json carries no axis values for a tutorial-record sweep, so
    dx/dt come from the sweep's own sweep_manifest.json instead.
    """
    axis_values_by_case: dict[str, dict] = {}
    manifest_path = output_dir / "sweep_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        axis_values_by_case = {
            case["case_id"]: case.get("resolved_axis_values", {})
            for case in manifest.get("cases", [])
        }
    for summary_path in sorted(output_dir.glob("cases/*/postProcessing/*_cv_summary.json")):
        case_dir = summary_path.parent.parent
        payload = json.loads(summary_path.read_text())
        yield case_dir, axis_values_by_case.get(case_dir.name, {}), payload


def build_summary_rows(output_dir: Path) -> list[dict]:
    cases = list(_iter_cv_summaries(output_dir))
    if not cases:
        print(f"[cable1DCVConvergence/table_summary] No *_cv_summary.json files in {output_dir}")
        return []

    rows: list[dict] = []
    grouped: dict[tuple[str, str, str, int], list[dict]] = defaultdict(list)

    for case_dir, axis_values, payload in cases:
        solver, ionic_model, tissue = _case_labels(case_dir)
        central = payload["central_cv"]
        dx_m = axis_values.get("dx", axis_values.get("tetDx", central["dx_m"]))
        dt_s = axis_values.get("system/controlDict:deltaT", central["dt_s"])
        row = {
            "case_id": payload["case_id"],
            "solver": solver,
            "ionic_model": ionic_model,
            "tissue": tissue,
            "DT_ms": round(float(dt_s) * 1e3, 6),
            "DX_mm": round(float(dx_m) * 1e3, 6),
            "conductivity_id": _CONDUCTIVITY_ID,
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
        print("[cable1DCVConvergence/table_summary] matplotlib unavailable; writing CSV only.")
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
    dx_path = output_dir / f"cable1DCVConvergence_dx_{tag}.png"
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
    dt_path = output_dir / f"cable1DCVConvergence_dt_{tag}.png"
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
        tutorial="cable1DCVConvergence",
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
    rows = build_summary_rows(output_path)
    if not rows:
        return []

    artifacts: list[dict[str, object]] = [
        _write_summary_csv(
            rows,
            output_dir=output_path,
            filename_stem="cable1DCVConvergence_summary",
            label="1D cable CV convergence summary",
            metadata=meta,
        )
    ]

    grouped_rows: dict[tuple[str, str, str, int], list[dict]] = defaultdict(list)
    for row in rows:
        grouped_rows[_group_key(row)].append(row)

    for group, rows_for_group in sorted(grouped_rows.items()):
        artifacts.extend(
            _plot_group_convergence(
                group,
                rows_for_group,
                output_dir=output_path,
            )
        )

    return artifacts


if __name__ == "__main__":
    folder = Path(__file__).resolve().parents[1] / "outputsCVConvergence"
    print(f"[table_summary] Default folder = {folder}")
    run_postprocessing(output_dir=str(folder))
