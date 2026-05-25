#!/usr/bin/env python3
"""compare_niederer_bueno_orovio_batched.py

Compare activation times sampled at Niederer benchmark points and line across
cpu / batched_euler / batched_heun / batched_rl / batched_soa modes.

Usage
-----
    python3 compare_niederer_bueno_orovio_batched.py \
        <run-root> <output-json> <mode> [<mode> ...]

The run-root must contain one directory per mode, each with sub-directories
run01/, run02/, … produced by run_niederer_bueno_orovio_batched_comparison.sh.

Outputs
-------
    <output-json>   – JSON file with summary table and per-point differences
"""

import json
import math
import re
import statistics
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# OpenFOAM probe file reader
# ---------------------------------------------------------------------------

def _parse_foam_probe_file(path: Path) -> tuple[list[float], list[list[float]]]:
    """Parse an OpenFOAM probes output file.

    The file format written by ``postProcess -func <probes>`` is::

        # Probe 0 at (x y z)
        # Probe 1 at (x y z)
        ...
        # Time   probe0   probe1  ...
        0.04     val0     val1    ...

    Returns
    -------
    times : list[float]
        Simulation times present in the file (usually just the one latestTime).
    values : list[list[float]]
        One inner list per time step; each inner list has one value per probe.
    """
    times = []
    values = []
    in_data = False

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#"):
            # The last comment line is the column header; mark start of data
            in_data = True
            continue
        if in_data:
            parts = line.split()
            if not parts:
                continue
            try:
                row = [float(v) for v in parts]
            except ValueError:
                continue
            times.append(row[0])
            values.append(row[1:])

    return times, values


def read_activation_times(probe_dir: Path) -> list[float]:
    """Return the activation-time vector at the *latest* written time step.

    ``probe_dir`` is the postProcessing/<funcName>/ directory.  OpenFOAM
    writes one sub-directory per written time step; we pick the largest one.
    Inside it, we read the file named ``activationTime``.

    Returns a flat list with one float per probe location.
    """
    time_dirs = sorted(
        (d for d in probe_dir.iterdir() if d.is_dir()),
        key=lambda d: float(d.name),
    )
    if not time_dirs:
        raise RuntimeError(f"No time directories in {probe_dir}")

    latest = time_dirs[-1]
    probe_file = latest / "activationTime"
    if not probe_file.exists():
        raise RuntimeError(f"activationTime not found in {latest}")

    times, values = _parse_foam_probe_file(probe_file)
    if not values:
        raise RuntimeError(f"No data rows parsed from {probe_file}")

    # Return the last time-step row (typically the only one for -latestTime)
    return values[-1]


# ---------------------------------------------------------------------------
# Wall-time and OpenFOAM execution-time readers
# ---------------------------------------------------------------------------

EXEC_RE = re.compile(
    r"ExecutionTime\s*=\s*([0-9eE+\-.]+)\s*s\s+ClockTime\s*=\s*([0-9eE+\-.]+)\s*s"
)


def read_wall_time(run_dir: Path) -> float:
    for line in (run_dir / "time.txt").read_text().splitlines():
        if line.startswith("real "):
            return float(line.split()[1])
    raise RuntimeError(f"Could not find wall time in {run_dir / 'time.txt'}")


def read_openfoam_time(run_dir: Path) -> tuple[float, float]:
    execution_time = clock_time = None
    for line in (run_dir / "log.cardiacFoam").read_text().splitlines():
        m = EXEC_RE.search(line)
        if m:
            execution_time = float(m.group(1))
            clock_time = float(m.group(2))
    if execution_time is None:
        raise RuntimeError(f"Could not read OpenFOAM timing from {run_dir}")
    return execution_time, clock_time


# ---------------------------------------------------------------------------
# Run-directory discovery
# ---------------------------------------------------------------------------

def run_dirs(mode_dir: Path) -> list[Path]:
    runs = sorted(
        d for d in mode_dir.iterdir()
        if d.is_dir() and d.name.startswith("run")
    )
    if not runs:
        raise RuntimeError(f"No run directories found in {mode_dir}")
    return runs


# ---------------------------------------------------------------------------
# Per-mode aggregation
# ---------------------------------------------------------------------------

def _stats(values: list[float]) -> dict:
    return {
        "mean":   statistics.fmean(values),
        "min":    min(values),
        "max":    max(values),
        "stddev": statistics.stdev(values) if len(values) > 1 else 0.0,
    }


def case_info(mode_dir: Path) -> dict:
    runs = run_dirs(mode_dir)
    wall_times = []
    exec_times = []
    clock_times = []

    for run_dir in runs:
        wall_times.append(read_wall_time(run_dir))
        et, ct = read_openfoam_time(run_dir)
        exec_times.append(et)
        clock_times.append(ct)

    last_run = runs[-1]
    return {
        "n_runs": len(runs),
        "external_wall_time_s": _stats(wall_times),
        "openfoam_execution_time_s": _stats(exec_times),
        "openfoam_clock_time_s": _stats(clock_times),
        "last_run_dir": str(last_run),
        "log_file": str(last_run / "log.cardiacFoam"),
    }


# ---------------------------------------------------------------------------
# Accuracy comparison helpers
# ---------------------------------------------------------------------------

def _abs_diff_stats(
    ref: list[float],
    test: list[float],
    scale: float = 1000.0,
) -> dict:
    """Compute per-point absolute differences (converted to ms by *scale*).

    Parameters
    ----------
    ref, test : list[float]
        Activation time vectors in seconds.
    scale : float
        Multiply differences by this factor before storing (1000 → ms).
    """
    if len(ref) != len(test):
        raise RuntimeError(
            f"Probe vector length mismatch: ref={len(ref)}, test={len(test)}"
        )

    diffs_ms = [abs(t - r) * scale for r, t in zip(ref, test)]

    max_i = max(range(len(diffs_ms)), key=diffs_ms.__getitem__)
    rmse_ms = math.sqrt(sum(d * d for d in diffs_ms) / len(diffs_ms))

    return {
        "max_abs_diff_ms":  diffs_ms[max_i],
        "mean_abs_diff_ms": statistics.fmean(diffs_ms),
        "rmse_ms":          rmse_ms,
        "max_abs_diff_at_probe": max_i,
        "per_probe_diff_ms": diffs_ms,
    }


# ---------------------------------------------------------------------------
# Summary row builders
# ---------------------------------------------------------------------------

def _number(v) -> str:
    return "n/a" if v is None else f"{v:.6g}"


def _speed_text(speedup) -> str:
    if speedup is None:
        return "n/a"
    if speedup >= 1.0:
        return f"{speedup:.3f}x faster"
    return f"{1.0 / speedup:.3f}x slower"


def summary_row(mode: str, case: dict, accuracy: dict | None = None) -> dict:
    if mode == "cpu":
        return {
            "mode":                  "cpu",
            "role":                  "reference",
            "wall_time_mean_s":      case["external_wall_time_s"]["mean"],
            "wall_time_vs_cpu":      "reference",
            "max_pt_err_ms":         0.0,
            "mean_pt_err_ms":        0.0,
            "max_ln_err_ms":         0.0,
            "mean_ln_err_ms":        0.0,
            "plain_english":         "CPU adaptive BuenoOrovio reference.",
        }

    pts = accuracy["points"]
    lns = accuracy["line"]
    speedup = case.get("wall_time_speedup_vs_cpu")

    return {
        "mode":              mode,
        "role":              "comparison",
        "wall_time_mean_s":  case["external_wall_time_s"]["mean"],
        "wall_time_speedup_vs_cpu": speedup,
        "wall_time_vs_cpu":  _speed_text(speedup),
        "max_pt_err_ms":     pts["max_abs_diff_ms"],
        "mean_pt_err_ms":    pts["mean_abs_diff_ms"],
        "max_ln_err_ms":     lns["max_abs_diff_ms"],
        "mean_ln_err_ms":    lns["mean_abs_diff_ms"],
        "plain_english": (
            f"{mode}: {_speed_text(speedup)} than CPU; "
            f"max point error {_number(pts['max_abs_diff_ms'])} ms, "
            f"max line error {_number(lns['max_abs_diff_ms'])} ms."
        ),
    }


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

def print_summary(summary: list[dict]) -> None:
    print()
    print("Niederer benchmark comparison summary")
    print("--------------------------------------")
    print(
        f"{'mode':17s}  {'wall_s':>10s}  {'vs_cpu':>14s}  "
        f"{'max_pt_ms':>10s}  {'mean_pt_ms':>10s}  "
        f"{'max_ln_ms':>10s}  {'mean_ln_ms':>10s}"
    )
    for row in summary:
        print(
            f"{row['mode']:17s}  "
            f"{_number(row['wall_time_mean_s']):>10s}  "
            f"{row['wall_time_vs_cpu']:>14s}  "
            f"{_number(row['max_pt_err_ms']):>10s}  "
            f"{_number(row['mean_pt_err_ms']):>10s}  "
            f"{_number(row['max_ln_err_ms']):>10s}  "
            f"{_number(row['mean_ln_err_ms']):>10s}"
        )
    print()
    for row in summary:
        print(f"- {row['plain_english']}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    if len(sys.argv) < 4:
        print(
            "Usage: compare_niederer_bueno_orovio_batched.py "
            "<run-root> <output-json> <mode> [<mode> ...]",
            file=sys.stderr,
        )
        return 1

    run_root    = Path(sys.argv[1])
    output_json = Path(sys.argv[2])
    modes       = sys.argv[3:]

    if "cpu" not in modes:
        raise RuntimeError("Mode list must include 'cpu' as the reference")

    # --- Collect per-mode case info ---
    cases: dict[str, dict] = {}
    for mode in modes:
        mode_dir = run_root / mode
        if not mode_dir.is_dir():
            raise RuntimeError(f"Mode directory not found: {mode_dir}")
        print(f"Reading mode: {mode} ...")
        cases[mode] = case_info(mode_dir)

    # --- Read CPU reference activation times (last run) ---
    cpu_last_run = Path(cases["cpu"]["last_run_dir"])
    cpu_points = read_activation_times(cpu_last_run / "postProcessing" / "Niedererpoints")
    cpu_line   = read_activation_times(cpu_last_run / "postProcessing" / "Niedererlines")
    cpu_wall   = cases["cpu"]["external_wall_time_s"]["mean"]

    print(f"CPU reference: {len(cpu_points)} benchmark points, {len(cpu_line)} line probes")

    # --- Compare each non-cpu mode ---
    accuracy: dict[str, dict] = {}
    for mode in modes:
        if mode == "cpu":
            continue

        last_run = Path(cases[mode]["last_run_dir"])
        mode_points = read_activation_times(last_run / "postProcessing" / "Niedererpoints")
        mode_line   = read_activation_times(last_run / "postProcessing" / "Niedererlines")

        cases[mode]["wall_time_speedup_vs_cpu"] = (
            cpu_wall / cases[mode]["external_wall_time_s"]["mean"]
        )

        accuracy[mode] = {
            "points": _abs_diff_stats(cpu_points, mode_points),
            "line":   _abs_diff_stats(cpu_line,   mode_line),
        }

    # --- Build summary ---
    summary = [summary_row("cpu", cases["cpu"])]
    for mode in modes:
        if mode == "cpu":
            continue
        summary.append(summary_row(mode, cases[mode], accuracy[mode]))

    result = {
        "reference_mode": "cpu",
        "mesh_dx_mm":     0.2,
        "delta_t_s":      5e-5,
        "end_time_s":     0.05,
        "n_subdomains":   6,
        "summary":        summary,
        "cases":          cases,
        "accuracy_vs_cpu": accuracy,
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2))
    print(f"\nMetrics written to {output_json}")
    print_summary(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
