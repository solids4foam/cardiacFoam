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
#     run-single-cell-substep-sweep
#
# Description
#     Executes parameter sweeps across individual cellular models.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


ROOT_DIR = Path(__file__).resolve().parents[4]
DRIVER_DIR = ROOT_DIR / "applications" / "scripts" / "driverFoam"
sys.path.insert(0, str(DRIVER_DIR))

from openfoam_driver.core.runtime.mutators import update_control_dict, update_foam_entry
from openfoam_driver.postprocessing.single_cell_mode_compare import compare_single_cell_modes


DEFAULT_MODELS = [
    "AlievPanfilovcompactBatched",
    "BuenoOroviocompactBatched",
    "CourtemanchecompactBatched",
    "GaurcompactBatched",
    "GrandicompactBatched",
    "PerisYaguecompactBatched",
    "StewartcompactBatched",
    "TNNPcompactBatched",
    "ToRORd_dynClcompactBatched",
    "TrovatocompactBatched",
    "TWorldcompactBatched",
]

EXCLUDED_CPU_MODELS = {"Fabbri"}
EXCLUDED_BATCHED_MODELS = {"FabbricompactBatched"}

CPU_STIMULUS_MAP = {
    "AlievPanfilov": 0.5,
    "BuenoOrovio": 0.4,
    "Courtemanche": 60,
    "Gaur": 60,
    "Grandi": 60,
    "PerisYague": 60,
    "Stewart": 60,
    "TNNP": 60,
    "ToRORd_dynCl": 60,
    "Trovato": 60,
}

BATCHED_STIMULUS_MAP = {
    "AlievPanfilovcompactBatched": 0.5,
    "BuenoOroviocompactBatched": 0.4,
    "CourtemanchecompactBatched": 60,
    "GaurcompactBatched": 60,
    "GrandicompactBatched": 60,
    "PerisYaguecompactBatched": 60,
    "StewartcompactBatched": 60,
    "TNNPcompactBatched": 60,
    "ToRORd_dynClcompactBatched": 60,
    "TrovatocompactBatched": 60,
    "TWorldcompactBatched": 80000,
}


def _remove_excluded_models(
    single_cell_config: dict[str, Any],
    excluded_models: set[str],
) -> None:
    models = [
        str(model)
        for model in single_cell_config.get("ionic_models", [])
        if str(model) not in excluded_models
    ]
    single_cell_config["ionic_models"] = models

    tissue_map = single_cell_config.get("ionic_model_tissue_map")
    if isinstance(tissue_map, dict):
        single_cell_config["ionic_model_tissue_map"] = {
            key: value
            for key, value in tissue_map.items()
            if str(key) not in excluded_models
        }

    stimulus_map = single_cell_config.get("stimulus_map")
    if isinstance(stimulus_map, dict):
        single_cell_config["stimulus_map"] = {
            key: value
            for key, value in stimulus_map.items()
            if str(key) not in excluded_models
        }


def _parse_steps(value: str) -> list[int]:
    steps = [int(part) for part in re.split(r"[,\s]+", value.strip()) if part]
    if not steps or any(step < 1 for step in steps):
        raise argparse.ArgumentTypeError("steps must contain positive integers")
    return sorted(set(steps), reverse=True)


def _float_text(value: float) -> str:
    return f"{value:.12g}"


def _mode_name(family: str, steps: int) -> str:
    return f"{family}_step{steps:03d}"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected object in {path}")
    return payload


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def _copy_case_template(source_case: Path, target_case: Path) -> None:
    if target_case.exists():
        shutil.rmtree(target_case)
    target_case.mkdir(parents=True)
    shutil.copytree(source_case / "constant", target_case / "constant")
    shutil.copytree(source_case / "system", target_case / "system")


def _set_single_cell_voltage_export(electro_properties: Path) -> None:
    update_foam_entry(
        electro_properties,
        "export",
        "(Vm)",
        scope=("singleCellSolverCoeffs", "outputVariables", "ionic"),
    )


def _set_initial_ode_step(electro_properties: Path, value: str) -> None:
    update_foam_entry(
        electro_properties,
        "initialODEStep",
        value,
        scope="singleCellSolverCoeffs",
    )


def _prepare_cpu_case(
    *,
    source_case: Path,
    target_case: Path,
    delta_t: str,
    end_time: str,
    cpu_initial_step: str,
) -> None:
    _copy_case_template(source_case, target_case)
    electro_properties = target_case / "constant" / "electroProperties"
    update_control_dict(
        target_case / "system" / "controlDict",
        delta_t=delta_t,
        end_time=end_time,
    )
    _set_initial_ode_step(electro_properties, cpu_initial_step)
    _set_single_cell_voltage_export(electro_properties)


def _prepare_batched_case(
    *,
    family: str,
    steps: int,
    source_case: Path,
    target_case: Path,
    delta_t: str,
    end_time: str,
    initial_step: str,
) -> None:
    _copy_case_template(source_case, target_case)

    constant_dir = target_case / "constant"
    if family == "euler":
        template_name = "electroProperties.batched_euler"
    elif family == "rl":
        template_name = "electroProperties.batched_rl"
    elif family == "soa":
        template_name = "electroProperties.batched_soa"
    else:
        raise ValueError(f"Unknown family: {family}")

    template_path = constant_dir / template_name
    if template_path.exists():
        shutil.copyfile(template_path, constant_dir / "electroProperties")

    electro_properties = constant_dir / "electroProperties"
    update_control_dict(
        target_case / "system" / "controlDict",
        delta_t=delta_t,
        end_time=end_time,
    )
    _set_initial_ode_step(electro_properties, initial_step)
    _set_single_cell_voltage_export(electro_properties)
    update_foam_entry(
        electro_properties,
        "batchedSubsteps",
        steps,
        scope="singleCellSolverCoeffs",
    )
    if family == "soa":
        update_foam_entry(
            electro_properties,
            "useSoAEvaluator",
            "true",
            scope="singleCellSolverCoeffs",
        )
        update_foam_entry(
            electro_properties,
            "batchedIntegrator",
            "euler",
            scope="singleCellSolverCoeffs",
        )
    elif family == "rl":
        update_foam_entry(
            electro_properties,
            "batchedIntegrator",
            "rushLarsen",
            scope="singleCellSolverCoeffs",
        )
    else:
        update_foam_entry(
            electro_properties,
            "batchedIntegrator",
            "euler",
            scope="singleCellSolverCoeffs",
        )


def _load_config_template(path: Path) -> dict[str, Any]:
    payload = _load_json(path)
    if "singleCell" not in payload or not isinstance(payload["singleCell"], dict):
        raise KeyError(f"{path} does not contain a singleCell config object")
    return payload


def _write_mode_config(
    *,
    template_config: dict[str, Any],
    config_path: Path,
    case_dir_name: str,
    stimulus_map: dict[str, float],
    excluded_models: set[str],
) -> None:
    payload = json.loads(json.dumps(template_config))
    payload["singleCell"]["case_dir_name"] = case_dir_name
    payload["singleCell"]["stimulus_map"] = stimulus_map
    _remove_excluded_models(payload["singleCell"], excluded_models)
    _write_json(config_path, payload)


def _run_driver(config_path: Path, tutorials_root: Path) -> None:
    command = [
        str(ROOT_DIR / "applications" / "scripts" / "driverFoam" / "bin" / "driverFoam"),
        "sim",
        "--entry",
        "singleCell",
        "--config",
        str(config_path),
        "--tutorials-root",
        str(tutorials_root),
        "--continue-on-error",
    ]
    subprocess.run(command, cwd=ROOT_DIR, check=True)


def _manifest_path(case_root: Path) -> Path:
    return case_root / "driverOutput" / "run_manifest.json"


def _manifest_status(manifest_path: Path) -> tuple[int, int, list[str]]:
    payload = _load_json(manifest_path)
    completed = int(payload.get("completed_cases", 0))
    failed = int(payload.get("failed_cases", 0))
    failed_cases = [
        str(result.get("case_id"))
        for result in payload.get("results", [])
        if isinstance(result, dict) and result.get("status") != "ok"
    ]
    return completed, failed, failed_cases


def _aggregate_rows(
    *,
    summary_csv: Path,
    manifests: dict[str, Path],
    delta_t: float,
    steps: list[int],
) -> list[dict[str, Any]]:
    rows = list(csv.DictReader(summary_csv.open()))
    by_mode: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_mode.setdefault(row["mode"], []).append(row)

    aggregate_rows: list[dict[str, Any]] = []
    for family in ("euler", "rl", "soa"):
        for step in steps:
            mode = _mode_name(family, step)
            mode_rows = by_mode.get(mode, [])
            ok_rows = [
                row
                for row in mode_rows
                if row.get("status") == "ok" and row.get("mean_abs_diff")
            ]
            values = [float(row["mean_abs_diff"]) for row in ok_rows]
            completed, failed, failed_cases = _manifest_status(manifests[mode])
            aggregate_rows.append(
                {
                    "mode": mode,
                    "family": family,
                    "steps_per_pde_step": step,
                    "deltaT": delta_t,
                    "effective_dt": delta_t / step,
                    "completed_cases": completed,
                    "failed_cases": failed,
                    "failed_case_ids": ";".join(failed_cases),
                    "ok_comparisons": len(values),
                    "mean_mean_abs_diff": statistics.fmean(values) if values else "",
                    "median_mean_abs_diff": statistics.median(values) if values else "",
                    "max_mean_abs_diff": max(values) if values else "",
                }
            )
    return aggregate_rows


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "mode",
        "family",
        "steps_per_pde_step",
        "deltaT",
        "effective_dt",
        "completed_cases",
        "failed_cases",
        "failed_case_ids",
        "ok_comparisons",
        "mean_mean_abs_diff",
        "median_mean_abs_diff",
        "max_mean_abs_diff",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _plot_sweep(aggregate_rows: list[dict[str, Any]], output_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    output_dir.mkdir(parents=True, exist_ok=True)
    styles = {
        "euler": {"label": "Euler cell-major", "marker": "o"},
        "rl": {"label": "Rush-Larsen + Euler", "marker": "s"},
        "soa": {"label": "SoA Euler", "marker": "^"},
    }
    plot_paths: list[Path] = []

    for metric, ylabel, filename in [
        ("mean_mean_abs_diff", "Mean Vm mean-abs-diff vs CPU", "step_mean_error.png"),
        ("max_mean_abs_diff", "Max per-case Vm mean-abs-diff vs CPU", "step_max_error.png"),
    ]:
        fig, ax = plt.subplots(figsize=(8.5, 5.2))
        for family, style in styles.items():
            rows = [
                row
                for row in aggregate_rows
                if row["family"] == family and row[metric] != ""
            ]
            rows.sort(key=lambda row: int(row["steps_per_pde_step"]))
            if not rows:
                continue
            xs = [int(row["steps_per_pde_step"]) for row in rows]
            ys = [float(row[metric]) for row in rows]
            ax.plot(xs, ys, marker=style["marker"], linewidth=1.6, label=style["label"])
            for x, y, row in zip(xs, ys, rows):
                failed = int(row["failed_cases"])
                if failed:
                    ax.annotate(
                        f"{failed} fail",
                        (x, y),
                        xytext=(5, 5),
                        textcoords="offset points",
                        fontsize=8,
                    )
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_xlabel("ODE steps per PDE step")
        ax.set_ylabel(ylabel)
        ax.set_title("Single-cell fixed-step sweep")
        ax.grid(True, alpha=0.25, which="both")
        ax.legend(loc="best")
        path = output_dir / filename
        fig.tight_layout()
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot_paths.append(path)

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    for family, style in styles.items():
        rows = [row for row in aggregate_rows if row["family"] == family]
        rows.sort(key=lambda row: int(row["steps_per_pde_step"]))
        xs = [int(row["steps_per_pde_step"]) for row in rows]
        ys = [int(row["failed_cases"]) for row in rows]
        ax.plot(xs, ys, marker=style["marker"], linewidth=1.6, label=style["label"])
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("ODE steps per PDE step")
    ax.set_ylabel("Failed single-cell cases")
    ax.set_title("Single-cell fixed-step sweep stability")
    ax.grid(True, alpha=0.25, which="both")
    ax.legend(loc="best")
    path = output_dir / "step_failures.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    plot_paths.append(path)

    return plot_paths


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a single-cell Euler/Rush-Larsen/SoA fixed-step sweep."
    )
    parser.add_argument(
        "--steps",
        "--substeps",
        dest="steps",
        default="50 40 30 20 10",
        type=_parse_steps,
        help="ODE steps per PDE step. --substeps is kept as a compatibility alias.",
    )
    parser.add_argument("--delta-t", default="1e-4")
    parser.add_argument("--end-time", default="1")
    parser.add_argument("--cpu-initial-step", default="1e-6")
    parser.add_argument("--batched-initial-step", default="1e-6")
    parser.add_argument(
        "--tutorials-root",
        default="tutorials",
        type=Path,
    )
    parser.add_argument(
        "--source-case-root",
        default=Path("tutorials/comparisonResults/single_cell_driver_cases_euler_rl"),
        type=Path,
    )
    parser.add_argument(
        "--source-config-root",
        default=Path("tutorials/comparisonResults/single_cell_driver_configs_euler_rl"),
        type=Path,
    )
    parser.add_argument(
        "--output-root",
        default=Path("tutorials/comparisonResults/single_cell_step_sweep"),
        type=Path,
    )
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--skip-run", action="store_true")
    return parser


def main() -> int:
    parser = _build_arg_parser()
    args = parser.parse_args()

    os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="matplotlib-"))

    tutorials_root = (ROOT_DIR / args.tutorials_root).resolve()
    source_case_root = (ROOT_DIR / args.source_case_root).resolve()
    source_config_root = (ROOT_DIR / args.source_config_root).resolve()
    output_root = (ROOT_DIR / args.output_root).resolve()
    cases_root = output_root / "cases"
    configs_root = output_root / "configs"
    comparison_root = output_root / "comparison"
    plots_root = output_root / "plots"

    delta_t_text = str(args.delta_t)
    end_time_text = str(args.end_time)
    delta_t = float(delta_t_text)

    cpu_config_template = _load_config_template(source_config_root / "cpu.json")
    batched_config_template = _load_config_template(source_config_root / "batched.json")

    if not (source_case_root / "cpu").is_dir():
        raise FileNotFoundError(f"Missing source CPU case: {source_case_root / 'cpu'}")
    if not (source_case_root / "batched").is_dir():
        raise FileNotFoundError(f"Missing source batched case: {source_case_root / 'batched'}")

    print(f"Preparing single-cell fixed-step sweep in {output_root}")
    print(f"deltaT={delta_t_text}, endTime={end_time_text}, steps={args.steps}")

    cases_root.mkdir(parents=True, exist_ok=True)
    configs_root.mkdir(parents=True, exist_ok=True)

    cpu_case = cases_root / "cpu"
    _prepare_cpu_case(
        source_case=source_case_root / "cpu",
        target_case=cpu_case,
        delta_t=delta_t_text,
        end_time=end_time_text,
        cpu_initial_step=str(args.cpu_initial_step),
    )
    _write_mode_config(
        template_config=cpu_config_template,
        config_path=configs_root / "cpu.json",
        case_dir_name=str(cpu_case.relative_to(tutorials_root)),
        stimulus_map=CPU_STIMULUS_MAP,
        excluded_models=EXCLUDED_CPU_MODELS,
    )

    manifests: dict[str, Path] = {"cpu": _manifest_path(cpu_case)}
    mode_order = ["cpu"]

    for family in ("euler", "rl", "soa"):
        for steps in args.steps:
            mode = _mode_name(family, steps)
            target_case = cases_root / mode
            _prepare_batched_case(
                family=family,
                steps=steps,
                source_case=source_case_root / "batched",
                target_case=target_case,
                delta_t=delta_t_text,
                end_time=end_time_text,
                initial_step=str(args.batched_initial_step),
            )
            _write_mode_config(
                template_config=batched_config_template,
                config_path=configs_root / f"{mode}.json",
                case_dir_name=str(target_case.relative_to(tutorials_root)),
                stimulus_map=BATCHED_STIMULUS_MAP,
                excluded_models=EXCLUDED_BATCHED_MODELS,
            )
            manifests[mode] = _manifest_path(target_case)
            mode_order.append(mode)

    if args.prepare_only:
        print(f"Prepared configs in {configs_root}")
        return 0

    if not args.skip_run:
        for mode in mode_order:
            print(f"\n=== Running {mode} ===")
            _run_driver(configs_root / f"{mode}.json", tutorials_root)

    missing = [mode for mode, path in manifests.items() if not path.exists()]
    if missing:
        raise RuntimeError(f"Missing run manifests for modes: {', '.join(missing)}")

    print("\n=== Comparing modes ===")
    compare_single_cell_modes(
        manifests,
        comparison_root,
        reference_mode="cpu",
        plot_variable="Vm",
        make_plots=True,
        strict=False,
    )

    aggregate_rows = _aggregate_rows(
        summary_csv=comparison_root / "single_cell_average_differences.csv",
        manifests=manifests,
        delta_t=delta_t,
        steps=args.steps,
    )
    aggregate_csv = output_root / "single_cell_step_summary.csv"
    _write_csv(aggregate_csv, aggregate_rows)
    plot_paths = _plot_sweep(aggregate_rows, plots_root)

    print(f"\nComparison artifacts: {comparison_root}")
    print(f"Step summary       : {aggregate_csv}")
    print(f"Step plots         : {plots_root}")
    for path in plot_paths:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
