#!/usr/bin/env python3
"""Separate outer-deltaT and CUDA Euler-substep convergence on a fixed mesh."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import time
from pathlib import Path

from run_manufactured_cuda_smoke import prepare_case


FIELDS = ("Vm", "u1", "u2")
NORMS = ("L1", "L2", "Linf")


def run(case: Path, label: str, require_gpu: bool) -> dict[str, object]:
    for command in ("blockMesh", "cardiacFoam"):
        with (case / f"log.{command}").open("w") as log:
            result = subprocess.run([command, "-case", str(case)], cwd=case,
                                    stdout=log, stderr=subprocess.STDOUT,
                                    check=False)
        if result.returncode:
            return {"status": f"{command} exit {result.returncode}",
                    "log_tail": (case / f"log.{command}").read_text(
                        errors="replace")[-3000:]}

    log = (case / "log.cardiacFoam").read_text(errors="replace")
    backend = "GPU" if "using CUDA device" in log else "CPU"
    if require_gpu and backend != "GPU":
        return {"status": "GPU not selected", "backend": backend,
                "log_tail": log[-3000:]}
    errors = {}
    for line in log.splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[0] in FIELDS:
            errors.update({f"{parts[0]}_{norm}": value
                           for norm, value in zip(NORMS, parts[1:4])})
    if len(errors) != len(FIELDS) * len(NORMS):
        return {"status": "missing manufactured errors", "backend": backend,
                "log_tail": log[-3000:]}
    execution = re.findall(r"ExecutionTime\s*=\s*([0-9.eE+-]+)", log)
    return {"status": "completed", "backend": backend,
            "runtime_s": execution[-1] if execution else "", **errors,
            "log_tail": ""}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--delta-ts", nargs="+", type=float,
                        default=[2.5e-6, 5e-6, 1e-5])
    parser.add_argument("--substeps", nargs="+", type=int, default=[1, 5, 25])
    parser.add_argument("--end-time", type=float, default=5e-5)
    parser.add_argument("--cells", type=int, default=8)
    parser.add_argument("--require-gpu", action="store_true")
    args = parser.parse_args()
    if (not args.delta_ts or min(args.delta_ts) <= 0 or not args.substeps
            or min(args.substeps) < 1 or args.end_time <= 0 or args.cells < 2):
        parser.error("deltaTs/endTime must be positive, substeps >= 1, cells >= 2")

    root = args.output_dir.resolve()
    if root.exists():
        parser.error(f"output directory already exists: {root}")
    root.mkdir(parents=True)
    rows = []
    for dt in sorted(set(args.delta_ts)):
        dt_tag = f"{dt:.0e}"
        scalar = root / f"scalar_dt{dt_tag}"
        prepare_case(scalar, "monodomainFDAManufactured", dt,
                     args.end_time, args.cells, 1)
        started = time.monotonic()
        result = run(scalar, "scalar", False)
        result.update({"backend": "scalar", "runtime_wall_s":
                       f"{time.monotonic() - started:.6g}"})
        rows.append({"case": scalar.name, "model": "scalar_RKF45",
                     "deltaT_s": dt, "substeps": "", **result})
        print(f"{scalar.name}: {result['status']}", flush=True)

        for steps in sorted(set(args.substeps)):
            case = root / f"gpu_dt{dt_tag}_{steps}substeps"
            prepare_case(case, "monodomainFDAManufacturedBatched", dt,
                         args.end_time, args.cells, steps)
            started = time.monotonic()
            result = run(case, "gpu", args.require_gpu)
            result["runtime_wall_s"] = f"{time.monotonic() - started:.6g}"
            rows.append({"case": case.name, "model": "batched_Euler",
                         "deltaT_s": dt, "substeps": steps, **result})
            print(f"{case.name}: {result['status']} ({result.get('backend', '')})",
                  flush=True)
        with (root / "summary.csv").open("w", newline="") as stream:
            fields = sorted({key for row in rows for key in row})
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    fields = sorted({key for row in rows for key in row})
    with (root / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return 0 if all(row.get("status") == "completed" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
