#!/usr/bin/env python3

import json
import math
import re
import statistics
import sys
from pathlib import Path


EXEC_RE = re.compile(
    r"ExecutionTime\s*=\s*([0-9eE+\-.]+)\s*s\s+ClockTime\s*=\s*([0-9eE+\-.]+)\s*s"
)


def read_trace(path):
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    header = lines[0].split()
    rows = [[float(value) for value in line.split()] for line in lines[1:]]
    return header, rows


def trace_file(run_dir):
    traces = sorted((run_dir / "postProcessing").glob("*.txt"))
    if len(traces) != 1:
        raise RuntimeError(f"Expected one trace in {run_dir}, found {len(traces)}")
    return traces[0]


def run_dirs(mode_dir):
    runs = sorted(
        item for item in mode_dir.iterdir()
        if item.is_dir() and item.name.startswith("run")
    )
    if not runs:
        raise RuntimeError(f"No run directories found in {mode_dir}")
    return runs


def read_wall_time(run_dir):
    for line in (run_dir / "time.txt").read_text().splitlines():
        if line.startswith("real "):
            return float(line.split()[1])
    raise RuntimeError(f"Could not find wall time in {run_dir / 'time.txt'}")


def read_openfoam_time(run_dir):
    execution_time = None
    clock_time = None

    for line in (run_dir / "log.cardiacFoam").read_text().splitlines():
        match = EXEC_RE.search(line)
        if match:
            execution_time = float(match.group(1))
            clock_time = float(match.group(2))

    if execution_time is None or clock_time is None:
        raise RuntimeError(f"Could not read OpenFOAM time in {run_dir}")

    return execution_time, clock_time


def stats(values):
    return {
        "mean": statistics.fmean(values),
        "min": min(values),
        "max": max(values),
        "stddev": statistics.stdev(values) if len(values) > 1 else 0.0,
    }


def case_info(mode_dir):
    runs = run_dirs(mode_dir)
    execution_times = []
    clock_times = []
    wall_times = []

    for run_dir in runs:
        execution_time, clock_time = read_openfoam_time(run_dir)
        execution_times.append(execution_time)
        clock_times.append(clock_time)
        wall_times.append(read_wall_time(run_dir))

    last_run = runs[-1]
    return {
        "n_runs": len(runs),
        "execution_time_s": stats(execution_times),
        "clock_time_s": stats(clock_times),
        "external_wall_time_s": stats(wall_times),
        "trace_file": str(trace_file(last_run)),
        "log_file": str(last_run / "log.cardiacFoam"),
    }


def compare_traces(cpu_header, cpu_rows, mode_header, mode_rows):
    if cpu_header != mode_header:
        raise RuntimeError(f"Trace headers do not match: {cpu_header} != {mode_header}")
    if len(cpu_rows) != len(mode_rows):
        raise RuntimeError("Trace lengths do not match")

    variables = {}

    for col, name in enumerate(cpu_header):
        if name == "time":
            continue

        diffs = [mode[col] - cpu[col] for cpu, mode in zip(cpu_rows, mode_rows)]
        abs_diffs = [abs(value) for value in diffs]
        max_i = max(range(len(abs_diffs)), key=abs_diffs.__getitem__)
        rmse = math.sqrt(sum(value*value for value in diffs)/len(diffs))

        variables[name] = {
            "max_abs_diff": abs_diffs[max_i],
            "rmse": rmse,
            "final_abs_diff": abs_diffs[-1],
            "reference_final": cpu_rows[-1][col],
            "test_final": mode_rows[-1][col],
            "max_abs_diff_time": cpu_rows[max_i][0],
            "reference_at_max_abs_diff": cpu_rows[max_i][col],
            "test_at_max_abs_diff": mode_rows[max_i][col],
        }

    return {
        "n_samples": len(cpu_rows),
        "time_end_s": cpu_rows[-1][0],
        "variables": variables,
    }


def number(value):
    if value is None:
        return "n/a"
    return f"{value:.6g}"


def speed_text(speedup):
    if speedup is None:
        return "n/a"
    if speedup >= 1.0:
        return f"{speedup:.3f}x faster"
    return f"{1.0/speedup:.3f}x slower"


def summary_row(mode, case, accuracy=None):
    if mode == "cpu":
        return {
            "mode": "cpu",
            "role": "reference",
            "wall_time_mean_s": case["external_wall_time_s"]["mean"],
            "wall_time_vs_cpu": "reference",
            "max_abs_vm_diff_mV": 0.0,
            "final_abs_vm_diff_mV": 0.0,
            "max_abs_s_diff": 0.0,
            "final_abs_s_diff": 0.0,
            "plain_english": "CPU adaptive BuenoOrovio reference.",
        }

    vm = accuracy["variables"].get("Vm", {})
    s_var = accuracy["variables"].get("s", {})
    speedup = case["wall_time_speedup_vs_cpu"]
    max_vm = vm.get("max_abs_diff")
    final_vm = vm.get("final_abs_diff")

    return {
        "mode": mode,
        "role": "comparison",
        "wall_time_mean_s": case["external_wall_time_s"]["mean"],
        "wall_time_speedup_vs_cpu": speedup,
        "wall_time_vs_cpu": speed_text(speedup),
        "max_abs_vm_diff_mV": max_vm,
        "final_abs_vm_diff_mV": final_vm,
        "max_abs_s_diff": s_var.get("max_abs_diff"),
        "final_abs_s_diff": s_var.get("final_abs_diff"),
        "plain_english": (
            f"{mode}: {speed_text(speedup)} than CPU; "
            f"max Vm difference {number(max_vm)} mV, "
            f"final Vm difference {number(final_vm)} mV."
        ),
    }


def print_summary(summary):
    print()
    print("Comparison summary")
    print("------------------")
    print(
        "mode              wall_mean_s  vs_cpu        maxVm_mV   finalVm_mV  "
        "max_s      final_s"
    )

    for row in summary:
        print(
            f"{row['mode']:17s} "
            f"{number(row['wall_time_mean_s']):>11s}  "
            f"{row['wall_time_vs_cpu']:>12s}  "
            f"{number(row['max_abs_vm_diff_mV']):>9s}  "
            f"{number(row['final_abs_vm_diff_mV']):>10s}  "
            f"{number(row['max_abs_s_diff']):>9s}  "
            f"{number(row['final_abs_s_diff']):>9s}"
        )

    print()
    for row in summary:
        print(f"- {row['plain_english']}")


def main():
    if len(sys.argv) < 4:
        print(
            "Usage: compare_bueno_orovio_batched.py "
            "<run-root> <output-json> <mode> [<mode> ...]",
            file=sys.stderr,
        )
        return 1

    run_root = Path(sys.argv[1])
    output_json = Path(sys.argv[2])
    modes = sys.argv[3:]

    if "cpu" not in modes:
        raise RuntimeError("Mode list must include cpu as the reference")

    cases = {mode: case_info(run_root / mode) for mode in modes}
    cpu_header, cpu_rows = read_trace(Path(cases["cpu"]["trace_file"]))
    cpu_wall = cases["cpu"]["external_wall_time_s"]["mean"]

    accuracy = {}
    summary = [summary_row("cpu", cases["cpu"])]

    for mode in modes:
        if mode == "cpu":
            continue

        header, rows = read_trace(Path(cases[mode]["trace_file"]))
        cases[mode]["wall_time_speedup_vs_cpu"] = (
            cpu_wall / cases[mode]["external_wall_time_s"]["mean"]
        )
        accuracy[mode] = compare_traces(cpu_header, cpu_rows, header, rows)
        summary.append(summary_row(mode, cases[mode], accuracy[mode]))

    result = {
        "reference_mode": "cpu",
        "summary": summary,
        "cases": cases,
        "accuracy_vs_cpu": accuracy,
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2))
    print_summary(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
