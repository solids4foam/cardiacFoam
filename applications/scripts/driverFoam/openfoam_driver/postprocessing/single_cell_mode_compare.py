from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence


os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="matplotlib-"))


@dataclass(frozen=True)
class TraceKey:
    model: str
    tissue: str


@dataclass
class TraceData:
    mode: str
    key: TraceKey
    raw_model: str
    case_id: str
    path: Path
    header: list[str]
    rows: list[list[float]]


def _strip_batched_suffix(model: str) -> str:
    if model.endswith("compactBatched"):
        return model.removesuffix("Batched")
    return model.removesuffix("Batched")


def _reference_key_for(key: TraceKey) -> TraceKey:
    if key.model.endswith("compact"):
        return TraceKey(key.model.removesuffix("compact"), key.tissue)
    return key


def _mode_sort_key(mode: str) -> tuple[int, int, str]:
    preferred = {
        "cpu": 0,
        "batched": 1,
        "batched_euler": 1,
        "euler": 1,
        "rl": 2,
        "batched_rl": 2,
        "rush_larsen": 2,
        "rushLarsen": 2,
        "soa": 3,
        "batched_soa": 3,
    }
    step_match = re.match(r"^(euler|rl|soa)_(?:step|sub)(\d+)$", mode)
    if step_match:
        family_order = {"euler": 1, "rl": 2, "soa": 3}[step_match.group(1)]
        # Larger step counts are more stable/reference-like, so show them first.
        return family_order, -int(step_match.group(2)), mode
    return preferred.get(mode, 100), 0, mode


def _safe_stem(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "trace"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected JSON object in {path}")
    return payload


def _read_trace(path: Path) -> tuple[list[str], list[list[float]]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"Trace file is empty: {path}")

    header = lines[0].lstrip("#").strip().split()
    if not header:
        raise ValueError(f"Trace file has no header: {path}")

    rows: list[list[float]] = []
    for line_no, line in enumerate(lines[1:], start=2):
        values = line.split()
        if len(values) != len(header):
            raise ValueError(
                f"{path}:{line_no}: expected {len(header)} columns, found {len(values)}"
            )
        rows.append([float(value) for value in values])

    if not rows:
        raise ValueError(f"Trace file has header but no data rows: {path}")
    return header, rows


def _known_tissue_names() -> tuple[str, ...]:
    return ("epicardialCells", "endocardialCells", "mCells", "myocyte")


def _infer_key_from_filename(path: Path) -> TraceKey | None:
    stem = path.stem
    for tissue in _known_tissue_names():
        marker = f"_{tissue}_"
        if marker in f"{stem}_":
            model = stem.split(marker, 1)[0]
            return TraceKey(_strip_batched_suffix(model), tissue)
    return None


def _matching_trace_for_result(
    output_dir: Path,
    *,
    case_id: str,
    ionic_model: str,
    tissue: str,
) -> Path | None:
    exact_patterns = [
        f"{ionic_model}_{tissue}_*.txt",
        f"{case_id}_*.txt",
        f"{_strip_batched_suffix(ionic_model)}_{tissue}_*.txt",
    ]
    for pattern in exact_patterns:
        matches = sorted(output_dir.glob(pattern))
        if matches:
            return matches[-1]

    normalized_model = _strip_batched_suffix(ionic_model)
    for trace_path in sorted(output_dir.glob("*.txt")):
        inferred = _infer_key_from_filename(trace_path)
        if inferred == TraceKey(normalized_model, tissue):
            return trace_path
    return None


def _traces_from_manifest(mode: str, manifest_path: Path) -> tuple[list[TraceData], list[str]]:
    manifest = _load_json(manifest_path)
    output_dir = Path(str(manifest.get("output_dir", manifest_path.parent)))
    if not output_dir.is_absolute():
        cwd_relative = output_dir
        if cwd_relative.exists():
            output_dir = cwd_relative
        else:
            output_dir = manifest_path.parent / output_dir

    warnings: list[str] = []
    traces: list[TraceData] = []
    manifest_results = manifest.get("results", [])
    has_manifest_results = isinstance(manifest_results, list) and bool(manifest_results)
    for result in manifest_results if isinstance(manifest_results, list) else []:
        if not isinstance(result, dict):
            continue
        status = str(result.get("status", ""))
        if status not in {"ok", "completed", ""}:
            warnings.append(
                f"{mode}: skipping case {result.get('case_id')} with status {status}"
            )
            continue

        params = result.get("params", {})
        if not isinstance(params, dict):
            params = {}
        ionic_model = str(params.get("ionicModel", ""))
        tissue = str(params.get("tissue", ""))
        case_id = str(result.get("case_id", f"{ionic_model}_{tissue}"))
        if not ionic_model or not tissue:
            warnings.append(f"{mode}: result {case_id} has no ionicModel/tissue params")
            continue

        trace_path = _matching_trace_for_result(
            output_dir,
            case_id=case_id,
            ionic_model=ionic_model,
            tissue=tissue,
        )
        if trace_path is None:
            warnings.append(f"{mode}: no trace file found for case {case_id}")
            continue

        header, rows = _read_trace(trace_path)
        traces.append(
            TraceData(
                mode=mode,
                key=TraceKey(_strip_batched_suffix(ionic_model), tissue),
                raw_model=ionic_model,
                case_id=case_id,
                path=trace_path,
                header=header,
                rows=rows,
            )
        )

    if traces:
        return traces, warnings

    if has_manifest_results:
        if not traces:
            warnings.append(f"{mode}: no successful trace files found in {output_dir}")
        return traces, warnings

    for trace_path in sorted(output_dir.glob("*.txt")):
        inferred = _infer_key_from_filename(trace_path)
        if inferred is None:
            warnings.append(f"{mode}: cannot infer model/tissue from {trace_path.name}")
            continue
        header, rows = _read_trace(trace_path)
        traces.append(
            TraceData(
                mode=mode,
                key=inferred,
                raw_model=inferred.model,
                case_id=f"{inferred.model}_{inferred.tissue}",
                path=trace_path,
                header=header,
                rows=rows,
            )
        )

    if not traces:
        warnings.append(f"{mode}: no trace files found in {output_dir}")
    return traces, warnings


def _time_column(header: Sequence[str]) -> int:
    try:
        return list(header).index("time")
    except ValueError as exc:
        raise ValueError(f"Trace header has no 'time' column: {header}") from exc


def _column(header: Sequence[str], name: str) -> int | None:
    try:
        return list(header).index(name)
    except ValueError:
        return None


def _series(trace: TraceData, variable: str) -> tuple[list[float], list[float]]:
    t_col = _time_column(trace.header)
    v_col = _column(trace.header, variable)
    if v_col is None:
        raise KeyError(f"{trace.path.name} has no variable {variable!r}")
    return [row[t_col] for row in trace.rows], [row[v_col] for row in trace.rows]


def _interp_at(times: Sequence[float], values: Sequence[float], target: float) -> float:
    if target <= times[0]:
        return values[0]
    if target >= times[-1]:
        return values[-1]

    lo = 0
    hi = len(times) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if times[mid] <= target:
            lo = mid
        else:
            hi = mid

    t0 = times[lo]
    t1 = times[hi]
    if t1 == t0:
        return values[lo]
    weight = (target - t0) / (t1 - t0)
    return values[lo] + weight * (values[hi] - values[lo])


def _compare_variable(
    reference: TraceData,
    candidate: TraceData,
    variable: str,
    *,
    reference_mode: str,
) -> dict[str, Any]:
    if _column(reference.header, variable) is None:
        return {
            "status": "missing_reference_variable",
            "note": f"{reference_mode} trace has no {variable}",
        }
    if _column(candidate.header, variable) is None:
        return {
            "status": "missing_candidate_variable",
            "note": f"{candidate.mode} trace has no {variable}",
        }

    ref_times, ref_values = _series(reference, variable)
    cand_times, cand_values = _series(candidate, variable)

    start = max(ref_times[0], cand_times[0])
    end = min(ref_times[-1], cand_times[-1])
    if end < start:
        return {
            "status": "no_time_overlap",
            "note": f"No overlapping time range for {reference.path.name} and {candidate.path.name}",
        }

    diffs: list[float] = []
    abs_diffs: list[float] = []
    sample_times: list[float] = []
    nonfinite_count = 0

    for time_value, ref_value in zip(ref_times, ref_values):
        if time_value < start or time_value > end:
            continue
        cand_value = _interp_at(cand_times, cand_values, time_value)
        diff = cand_value - ref_value
        if not (
            math.isfinite(time_value)
            and math.isfinite(ref_value)
            and math.isfinite(cand_value)
            and math.isfinite(diff)
        ):
            nonfinite_count += 1
        diffs.append(diff)
        abs_diffs.append(abs(diff))
        sample_times.append(time_value)

    if not diffs:
        return {
            "status": "no_samples",
            "note": f"No reference samples inside overlap for {candidate.path.name}",
        }

    if nonfinite_count:
        mean_abs = math.nan
        rmse = math.nan
        max_abs = math.nan
        max_index = len(abs_diffs) - 1
    else:
        mean_abs = sum(abs_diffs) / len(abs_diffs)
        rmse = math.sqrt(sum(value * value for value in diffs) / len(diffs))
        max_index = max(range(len(abs_diffs)), key=abs_diffs.__getitem__)
        max_abs = abs_diffs[max_index]

    return {
        "status": "ok",
        "note": "",
        "n_samples": len(diffs),
        "time_start": sample_times[0],
        "time_end": sample_times[-1],
        "mean_abs_diff": mean_abs,
        "max_abs_diff": max_abs,
        "rmse": rmse,
        "final_abs_diff": abs_diffs[-1] if not nonfinite_count else math.nan,
        "max_abs_diff_time": sample_times[max_index],
        "reference_final": ref_values[-1],
        "candidate_final": cand_values[-1],
        "nonfinite_count": nonfinite_count,
    }


def _all_variables(reference: TraceData, candidate: TraceData) -> list[str]:
    ref_names = set(reference.header)
    return [
        name
        for name in candidate.header
        if name != "time" and name in ref_names
    ]


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _format_cell(value: Any) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.9g}"
    return str(value)


def _write_html_table(path: Path, rows: list[dict[str, Any]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    parts = [
        "<!doctype html>",
        "<meta charset=\"utf-8\">",
        "<title>Single-cell mode comparison</title>",
        "<style>",
        "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:24px;}",
        "table{border-collapse:collapse;font-size:13px;}",
        "th,td{border:1px solid #ddd;padding:6px 8px;text-align:right;}",
        "th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left;}",
        "th{background:#f5f5f5;position:sticky;top:0;}",
        "</style>",
        "<h1>Single-cell mode comparison</h1>",
        "<table>",
        "<thead><tr>",
    ]
    parts.extend(f"<th>{html.escape(name)}</th>" for name in fieldnames)
    parts.append("</tr></thead><tbody>")
    for row in rows:
        parts.append("<tr>")
        for name in fieldnames:
            parts.append(f"<td>{html.escape(_format_cell(row.get(name, '')))}</td>")
        parts.append("</tr>")
    parts.extend(["</tbody></table>"])
    path.write_text("\n".join(parts))


def _import_pyplot():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    return plt


def _plot_model(
    output_dir: Path,
    model: str,
    traces: list[TraceData],
    *,
    modes: Sequence[str],
    variable: str,
) -> Path | None:
    tissues = sorted({trace.key.tissue for trace in traces})
    if not tissues:
        return None

    plt = _import_pyplot()
    fig, axes = plt.subplots(
        len(tissues),
        1,
        figsize=(10, max(3.5, 3.0 * len(tissues))),
        squeeze=False,
        sharex=False,
    )

    plotted_any = False
    for row, tissue in enumerate(tissues):
        ax = axes[row][0]
        tissue_traces = [trace for trace in traces if trace.key.tissue == tissue]
        by_mode = {trace.mode: trace for trace in tissue_traces}

        for mode in modes:
            trace = by_mode.get(mode)
            if trace is None or _column(trace.header, variable) is None:
                continue
            times, values = _series(trace, variable)
            ax.plot(times, values, linewidth=1.4, label=mode)
            plotted_any = True

        ax.set_title(f"{model} / {tissue}")
        ax.set_ylabel(variable)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best", fontsize=8)

    axes[-1][0].set_xlabel("time")
    fig.tight_layout()

    if not plotted_any:
        plt.close(fig)
        return None

    path = output_dir / f"single_cell_{_safe_stem(model)}_{_safe_stem(variable)}.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def _manifest_mode_from_path(path: Path) -> str:
    parent = path.parent
    if parent.name in {"postProcessing", "setupSingleCell"} and parent.parent.name:
        return parent.parent.name
    return parent.name


def discover_manifests(run_root: Path) -> dict[str, Path]:
    manifests: dict[str, Path] = {}
    for path in sorted(run_root.rglob("run_manifest.json")):
        mode = _manifest_mode_from_path(path)
        if mode in manifests:
            raise ValueError(
                f"Multiple manifests inferred for mode {mode!r}: {manifests[mode]} and {path}"
            )
        manifests[mode] = path
    return manifests


def compare_single_cell_modes(
    manifests: dict[str, Path],
    output_dir: Path,
    *,
    reference_mode: str = "cpu",
    plot_variable: str = "Vm",
    make_plots: bool = True,
    strict: bool = False,
) -> dict[str, Any]:
    if reference_mode not in manifests:
        raise KeyError(f"Reference mode {reference_mode!r} is not in manifest set")

    output_dir.mkdir(parents=True, exist_ok=True)

    traces: list[TraceData] = []
    warnings: list[str] = []
    for mode, manifest_path in manifests.items():
        mode_traces, mode_warnings = _traces_from_manifest(mode, Path(manifest_path))
        traces.extend(mode_traces)
        warnings.extend(mode_warnings)

    if not traces:
        raise RuntimeError("No single-cell traces found in the supplied manifests")

    modes = sorted({trace.mode for trace in traces}, key=_mode_sort_key)
    by_key_mode: dict[tuple[TraceKey, str], TraceData] = {}
    for trace in traces:
        by_key_mode[(trace.key, trace.mode)] = trace

    keys = sorted({trace.key for trace in traces}, key=lambda item: (item.model, item.tissue))
    metrics_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    for key in keys:
        reference_key = _reference_key_for(key)
        reference = by_key_mode.get((reference_key, reference_mode))
        if reference is None:
            warnings.append(
                f"No {reference_mode} reference trace for {key.model}/{key.tissue}"
            )
            continue

        for mode in modes:
            if mode == reference_mode and key == reference_key:
                continue
            candidate = by_key_mode.get((key, mode))
            if candidate is None:
                row = {
                    "model": key.model,
                    "tissue": key.tissue,
                    "mode": mode,
                    "reference_mode": reference_mode,
                    "variable": plot_variable,
                    "status": "missing_mode_trace",
                    "note": f"No {mode} trace for {key.model}/{key.tissue}",
                    "reference_trace_file": str(reference.path),
                    "candidate_trace_file": "",
                }
                metrics_rows.append(row)
                summary_rows.append(row)
                continue

            variables = _all_variables(reference, candidate)
            has_plot_variable = plot_variable in variables
            for variable in variables:
                metric = _compare_variable(
                    reference,
                    candidate,
                    variable,
                    reference_mode=reference_mode,
                )
                row = {
                    "model": key.model,
                    "tissue": key.tissue,
                    "mode": mode,
                    "reference_mode": reference_mode,
                    "variable": variable,
                    "status": metric.get("status", "ok"),
                    "note": metric.get("note", ""),
                    "n_samples": metric.get("n_samples", ""),
                    "time_start": metric.get("time_start", ""),
                    "time_end": metric.get("time_end", ""),
                    "mean_abs_diff": metric.get("mean_abs_diff", ""),
                    "max_abs_diff": metric.get("max_abs_diff", ""),
                    "rmse": metric.get("rmse", ""),
                    "final_abs_diff": metric.get("final_abs_diff", ""),
                    "max_abs_diff_time": metric.get("max_abs_diff_time", ""),
                    "nonfinite_count": metric.get("nonfinite_count", ""),
                    "reference_final": metric.get("reference_final", ""),
                    "candidate_final": metric.get("candidate_final", ""),
                    "reference_trace_file": str(reference.path),
                    "candidate_trace_file": str(candidate.path),
                }
                metrics_rows.append(row)
                if variable == plot_variable:
                    summary_rows.append(row)

            if not has_plot_variable:
                metric = _compare_variable(
                    reference,
                    candidate,
                    plot_variable,
                    reference_mode=reference_mode,
                )
                summary_rows.append(
                    {
                        "model": key.model,
                        "tissue": key.tissue,
                        "mode": mode,
                        "reference_mode": reference_mode,
                        "variable": plot_variable,
                        "status": metric.get("status", "ok"),
                        "note": metric.get("note", ""),
                        "reference_trace_file": str(reference.path),
                        "candidate_trace_file": str(candidate.path),
                    }
                )

    fieldnames = [
        "model",
        "tissue",
        "mode",
        "reference_mode",
        "variable",
        "status",
        "n_samples",
        "time_start",
        "time_end",
        "mean_abs_diff",
        "max_abs_diff",
        "rmse",
        "final_abs_diff",
        "max_abs_diff_time",
        "nonfinite_count",
        "reference_final",
        "candidate_final",
        "note",
        "reference_trace_file",
        "candidate_trace_file",
    ]
    summary_fieldnames = [
        "model",
        "tissue",
        "mode",
        "reference_mode",
        "variable",
        "status",
        "mean_abs_diff",
        "max_abs_diff",
        "rmse",
        "final_abs_diff",
        "nonfinite_count",
        "note",
    ]

    metrics_csv = output_dir / "single_cell_variable_metrics.csv"
    summary_csv = output_dir / "single_cell_average_differences.csv"
    summary_html = output_dir / "single_cell_average_differences.html"
    summary_json = output_dir / "single_cell_mode_comparison.json"
    _write_csv(metrics_csv, metrics_rows, fieldnames)
    _write_csv(summary_csv, summary_rows, summary_fieldnames)
    _write_html_table(summary_html, summary_rows, summary_fieldnames)

    plot_paths: list[Path] = []
    if make_plots:
        plots_dir = output_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        for model in sorted({trace.key.model for trace in traces}):
            model_traces = [trace for trace in traces if trace.key.model == model]
            plot_path = _plot_model(
                plots_dir,
                model,
                model_traces,
                modes=modes,
                variable=plot_variable,
            )
            if plot_path is not None:
                plot_paths.append(plot_path)

    payload = {
        "schema_version": "1.0",
        "reference_mode": reference_mode,
        "plot_variable": plot_variable,
        "modes": modes,
        "trace_count": len(traces),
        "warnings": warnings,
        "artifacts": {
            "summary_csv": str(summary_csv),
            "summary_html": str(summary_html),
            "metrics_csv": str(metrics_csv),
            "plots": [str(path) for path in plot_paths],
        },
        "summary_rows": summary_rows,
    }
    summary_json.write_text(json.dumps(payload, indent=2, allow_nan=True))

    bad_rows = [
        row
        for row in summary_rows
        if row.get("status") != "ok"
        or row.get("nonfinite_count") not in {"", 0, "0"}
    ]
    if strict and (warnings or bad_rows):
        details = "\n".join(
            [*warnings, *(f"{row['model']}/{row['tissue']} {row['mode']}: {row.get('status')}" for row in bad_rows)]
        )
        raise RuntimeError(f"Single-cell comparison failed strict validation:\n{details}")

    return {
        "artifacts": [
            {"path": str(summary_csv), "kind": "table", "format": "csv", "label": "Average differences"},
            {"path": str(summary_html), "kind": "table", "format": "html", "label": "Average differences"},
            {"path": str(metrics_csv), "kind": "table", "format": "csv", "label": "Variable metrics"},
            {"path": str(summary_json), "kind": "summary", "format": "json", "label": "Comparison summary"},
            *(
                {"path": str(path), "kind": "plot", "format": "png", "label": path.stem}
                for path in plot_paths
            ),
        ],
        "warnings": warnings,
    }


def _parse_manifest_arg(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            "Manifest must be passed as MODE=/path/to/run_manifest.json"
        )
    mode, path_text = value.split("=", 1)
    mode = mode.strip()
    if not mode:
        raise argparse.ArgumentTypeError("Manifest mode cannot be empty")
    return mode, Path(path_text).expanduser()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare singleCell driverFOAM runs across CPU/batched/RL modes "
            "using run_manifest.json files."
        )
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--manifest",
        action="append",
        type=_parse_manifest_arg,
        metavar="MODE=PATH",
        help="Mode-labelled driverFOAM run manifest. Repeat for each mode.",
    )
    source.add_argument(
        "--run-root",
        type=Path,
        help="Directory to recursively scan for run_manifest.json files.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--reference-mode", default="cpu")
    parser.add_argument("--plot-variable", default="Vm")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--strict", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.manifest:
        manifests = dict(args.manifest)
    else:
        manifests = discover_manifests(args.run_root)

    result = compare_single_cell_modes(
        manifests,
        args.output_dir,
        reference_mode=args.reference_mode,
        plot_variable=args.plot_variable,
        make_plots=not args.no_plots,
        strict=args.strict,
    )
    artifacts = result.get("artifacts", [])
    print(f"Wrote {len(artifacts)} comparison artifact(s) to {args.output_dir}")
    return 0


def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    manifests: dict[str, str] | None = None,
    run_root: str | None = None,
    reference_mode: str = "cpu",
    plot_variable: str = "Vm",
    make_plots: bool = True,
    strict: bool = False,
) -> dict[str, Any]:
    if manifests is None:
        if run_root is None:
            raise ValueError("run_postprocessing requires manifests or run_root")
        manifest_paths = discover_manifests(Path(run_root))
    else:
        manifest_paths = {mode: Path(path) for mode, path in manifests.items()}

    return compare_single_cell_modes(
        manifest_paths,
        Path(output_dir),
        reference_mode=reference_mode,
        plot_variable=plot_variable,
        make_plots=make_plots,
        strict=strict,
    )


if __name__ == "__main__":
    raise SystemExit(main())
