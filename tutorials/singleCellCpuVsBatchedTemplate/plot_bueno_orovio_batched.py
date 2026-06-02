#!/usr/bin/env python3

import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="matplotlib-"))

import matplotlib.pyplot as plt


def read_trace(path):
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    header = lines[0].split()
    rows = [[float(value) for value in line.split()] for line in lines[1:]]
    return header, rows


def latest_trace(run_root, mode):
    mode_dir = run_root / mode
    runs = sorted(
        item for item in mode_dir.iterdir()
        if item.is_dir() and item.name.startswith("run")
    )
    if not runs:
        raise RuntimeError(f"No run directories found for mode {mode}")

    traces = sorted((runs[-1] / "postProcessing").glob("*.txt"))
    if len(traces) != 1:
        raise RuntimeError(f"Expected one trace for mode {mode}, found {len(traces)}")
    return traces[0]


def mean_wall_time(run_root, mode):
    mode_dir = run_root / mode
    runs = sorted(
        item for item in mode_dir.iterdir()
        if item.is_dir() and item.name.startswith("run")
    )
    wall_times = []

    for run_dir in runs:
        for line in (run_dir / "time.txt").read_text().splitlines():
            if line.startswith("real "):
                wall_times.append(float(line.split()[1]))
                break

    if not wall_times:
        raise RuntimeError(f"No wall times found for mode {mode}")
    return sum(wall_times) / len(wall_times)


def check_trace(cpu_header, cpu_rows, mode, mode_header, mode_rows):
    if cpu_header != mode_header:
        raise RuntimeError(f"Trace headers do not match for {mode}")
    if len(cpu_rows) != len(mode_rows):
        raise RuntimeError(f"Trace lengths do not match for {mode}")


def main():
    if len(sys.argv) < 4:
        print(
            "Usage: plot_bueno_orovio_batched.py "
            "<run-root> <output-dir> <mode> [<mode> ...]",
            file=sys.stderr,
        )
        return 1

    run_root = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    modes = sys.argv[3:]

    if "cpu" not in modes:
        raise RuntimeError("Mode list must include cpu as the reference")

    output_dir.mkdir(parents=True, exist_ok=True)

    cpu_header, cpu_rows = read_trace(latest_trace(run_root, "cpu"))
    times = [row[0] for row in cpu_rows]

    traces = {}
    wall_times = {}

    for mode in modes:
        if mode == "cpu":
            continue
        header, rows = read_trace(latest_trace(run_root, mode))
        check_trace(cpu_header, cpu_rows, mode, header, rows)
        traces[mode] = rows
        wall_times[mode] = mean_wall_time(run_root, mode)

    for col, variable in enumerate(cpu_header):
        if variable == "time":
            continue

        cpu_values = [row[col] for row in cpu_rows]

        fig, ax = plt.subplots(figsize=(10, 5.5))
        scatter_x = []
        scatter_y = []
        scatter_labels = []

        for mode, rows in traces.items():
            mode_values = [row[col] for row in rows]
            abs_diff = [
                abs(mode_value - cpu_value)
                for cpu_value, mode_value in zip(cpu_values, mode_values)
            ]
            max_i = max(range(len(abs_diff)), key=abs_diff.__getitem__)

            ax.plot(times, abs_diff, label=mode, linewidth=1.5)
            ax.scatter([times[max_i]], [abs_diff[max_i]], s=22, zorder=3)

            scatter_x.append(wall_times[mode])
            scatter_y.append(abs_diff[max_i])
            scatter_labels.append(mode)

        ax.set_xlabel("time [s]")
        ax.set_ylabel(f"|mode - cpu| for {variable}")
        ax.set_title(f"{variable} accuracy vs CPU reference")
        ax.legend(loc="best")
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(output_dir / f"combined_{variable}_accuracy_vs_cpu.png", dpi=160)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5.5))
        ax.scatter(scatter_x, scatter_y, s=55)
        for x_value, y_value, label in zip(scatter_x, scatter_y, scatter_labels):
            ax.annotate(
                label,
                (x_value, y_value),
                xytext=(5, 4),
                textcoords="offset points",
                fontsize=8,
            )
        ax.set_xlabel("mean wall time [s]")
        ax.set_ylabel(f"max |mode - cpu| for {variable}")
        ax.set_title(f"{variable}: speed vs accuracy tradeoff")
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            output_dir / f"scatter_{variable}_wall_time_vs_accuracy.png",
            dpi=160,
        )
        plt.close(fig)

    print(f"Comparison plots written to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
