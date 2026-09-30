#!/usr/bin/env python3
"""Summarize a slab ionic-integrator sweep and select provisional settings."""

import argparse
import csv
import math
import re
from pathlib import Path


METRIC = re.compile(
    r"(?P<time>[0-9.eE+-]+)s cells=(?P<cells>\d+) "
    r"activated_scalar=(?P<scalar>\d+) activated_batched=(?P<batched>\d+) "
    r"activation_mask_mismatch=(?P<mask_mismatch>\d+) "
    r"Vm_RMSE_mV=(?P<rmse>[0-9.eE+-]+) Vm_max_mV=(?P<maximum>[0-9.eE+-]+) "
    r"activation_p95_ms=(?P<p95>[0-9.eE+-]+|nan)"
)


def metrics(text):
    return [
        {
            "cells": int(match["cells"]),
            "count_diff": abs(int(match["scalar"]) - int(match["batched"])),
            "mask_mismatch": int(match["mask_mismatch"]),
            "rmse": float(match["rmse"]),
            "maximum": float(match["maximum"]),
            "p95": float(match["p95"]),
        }
        for match in METRIC.finditer(text)
    ]


def passes(metrics_rows):
    # Require three configured output times, limited activation-mask mismatch,
    # and limited p95 timing difference. Do not threshold pointwise Vm error,
    # which can be dominated by a traveling upstroke.
    if len(metrics_rows) != 3:
        return False
    for metric in metrics_rows:
        mismatch_limit = max(5, math.ceil(0.01 * metric["cells"]))
        if metric["mask_mismatch"] > mismatch_limit:
            return False
        if not math.isfinite(metric["p95"]) or metric["p95"] > 0.1:
            return False
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    args = parser.parse_args()
    root = args.output_root.resolve()
    with (root / "summary.csv").open(newline="") as source:
        rows = list(csv.DictReader(source))

    models = sorted({row["model"] for row in rows})
    report = [
        "# Slab ionic-integrator sweep: provisional fastest settings",
        "",
        "This report is a short slab screen. A provisional pass requires "
        "a completed CUDA run, three valid configured output times, activation "
        "mask mismatch no greater than max(5 cells, 1% of slab cells), and "
        "p95 activation-time difference at most 0.1 ms against both matched-step "
        "and fine-step scalar RKF45. Pointwise Vm RMSE is reported but is not a "
        "screen gate because the traveling upstroke can dominate it. Passing this "
        "screen does not establish APD or full action-potential accuracy.",
        "",
        "| Model | Fastest passing GPU setting | Ionic substep | Runtime (s) | "
        "Speedup vs matched RKF45 | Worst activation p95 (ms) | Worst Vm RMSE (mV) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]

    for model in models:
        candidates = []
        for row in rows:
            if row["model"] != model or row.get("integrator") == "scalar_RKF45":
                continue
            if row.get("status") != "completed-GPU":
                continue
            matched = metrics(row.get("matched_RKF45", ""))
            fine = metrics(row.get("fine_RKF45", ""))
            if not passes(matched) or not passes(fine):
                continue
            candidates.append((row, matched, fine))

        if not candidates:
            report.append(
                f"| {model} | No provisional pass | — | — | — | — | — |"
            )
            continue

        row, matched, fine = min(
            candidates, key=lambda item: float(item[0]["runtime_s"])
        )
        dt = float(row["tissue_dt_s"])
        steps = int(row["substeps"])
        method = "Euler" if row["integrator"] == "euler" else "RL + Euler"
        label = f"{method}, dt={dt:.3g} s, {steps} substeps"
        scalar_runtime = float(row.get("matched_scalar_runtime_s") or "nan")
        runtime = float(row["runtime_s"])
        speedup = scalar_runtime/runtime if runtime > 0 else float("nan")
        p95_values = [
            metric["p95"] for metric in matched + fine
            if math.isfinite(metric["p95"])
        ]
        p95_text = f"{max(p95_values):.4g}" if p95_values else "n/a"
        report.append(
            f"| {model} | {label} | {dt/steps:.3g} s | {runtime:.2f} | "
            f"{speedup:.2f}x | {p95_text} | "
            f"{max(m['rmse'] for m in matched + fine):.4g} |"
        )

    report.extend([
        "",
        "The full per-case records, including failures and comparisons with the "
        "fine GPU Rush–Larsen baseline, are in `summary.csv`. Runtime speedup is "
        "relative to matched-step scalar RKF45 and includes mesh generation and "
        "the CPU tissue solve. Confirm selected "
        "settings with longer single-cell APD/state/current checks and a full "
        "propagation run before using them as validated production settings.",
    ])
    (root / "best_settings.md").write_text("\n".join(report) + "\n")
    print(root / "best_settings.md")


if __name__ == "__main__":
    main()
