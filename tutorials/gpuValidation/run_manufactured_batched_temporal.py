#!/usr/bin/env python3
"""Run the fixed-mesh manufactured temporal ladder with the batched model."""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import time
from pathlib import Path

from run_sbdf2_manufactured_temporal import (
    DT_DEFAULTS,
    ERROR_RE,
    prepare_case as prepare_scalar_case,
)


FIELDS = ("Vm", "u1", "u2")
NORMS = ("L1", "L2", "Linf")


def prepare_case(case: Path, dt: float, cells: int, end_time: float,
                 substeps: int, scheme: str) -> None:
    prepare_scalar_case(case, dt, cells, end_time)

    electro_path = case / "constant/electroProperties"
    electro = electro_path.read_text()
    electro, count = re.subn(
        r"(?m)^(\s*ionicModel\s+)monodomainFDAManufactured\s*;",
        rf"\g<1>monodomainFDAManufacturedBatched;\n"
        rf"    batchedIntegrator euler;\n    batchedSubsteps {substeps};",
        electro,
        count=1,
    )
    if count != 1:
        raise RuntimeError(f"Could not set batched manufactured model in {electro_path}")
    if scheme == "godunov":
        electro, count = re.subn(
            r"(?m)^(\s*timeCouplingScheme\s+)sbdf2\s*;",
            r"\g<1>godunov;",
            electro,
            count=1,
        )
        if count != 1:
            raise RuntimeError(f"Could not set Godunov coupling in {electro_path}")
    electro_path.write_text(electro)

    schemes_path = case / "system/fvSchemes"
    schemes = schemes_path.read_text()
    if scheme == "godunov":
        schemes, count = re.subn(
            r"(?m)^(\s*ddt\(Vm\)\s+)backward\s*;",
            r"\g<1>Euler;",
            schemes,
            count=1,
        )
        if count != 1:
            raise RuntimeError(f"Could not set Euler ddt(Vm) in {schemes_path}")
    schemes_path.write_text(schemes)


def run_case(case: Path, require_gpu: bool) -> dict[str, str]:
    for command in ("blockMesh", "cardiacFoam"):
        with (case / f"log.{command}").open("w") as log:
            result = subprocess.run(
                [command, "-case", str(case)], cwd=case,
                stdout=log, stderr=subprocess.STDOUT, check=False,
            )
        if result.returncode:
            return {
                "status": f"{command} exit {result.returncode}",
                "backend": "",
                "log_tail": (case / f"log.{command}").read_text(
                    errors="replace"
                )[-4000:],
            }

    log = (case / "log.cardiacFoam").read_text(errors="replace")
    backend = "GPU" if "using CUDA device" in log else "CPU"
    if require_gpu and backend != "GPU":
        return {"status": "GPU not selected", "backend": backend,
                "log_tail": log[-4000:]}

    errors = {field: values for field, *values in ERROR_RE.findall(log)}
    if set(errors) != set(FIELDS):
        return {"status": "missing manufactured errors", "backend": backend,
                "log_tail": log[-4000:]}

    row = {"status": "completed", "backend": backend, "log_tail": ""}
    for field, values in errors.items():
        for norm, value in zip(NORMS, values):
            row[f"{field}_{norm}"] = value
    runtime = re.findall(r"ExecutionTime\s*=\s*([0-9.eE+-]+)", log)
    steps = re.search(r"Number of steps\s*=\s*(\d+)", log)
    row["runtime_s"] = runtime[-1] if runtime else ""
    row["steps"] = steps.group(1) if steps else ""
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--backend", choices=("cpu", "gpu"), required=True)
    parser.add_argument("--scheme", choices=("godunov", "sbdf2"),
                        default="godunov")
    parser.add_argument("--delta-ts", nargs="+", type=float,
                        default=DT_DEFAULTS)
    parser.add_argument("--substeps", type=int, default=25)
    parser.add_argument("--substep-scaling", choices=("fixed", "inverse-deltaT"),
                        default="fixed")
    parser.add_argument("--cells", type=int, default=640)
    parser.add_argument("--end-time", type=float, default=0.2)
    parser.add_argument("--require-gpu", action="store_true")
    args = parser.parse_args()
    if (not args.delta_ts or min(args.delta_ts) <= 0 or args.substeps < 1
            or args.cells < 2 or args.end_time <= 0):
        parser.error("deltaT/endTime/substeps must be positive and cells >= 2")
    if args.backend == "gpu" and not args.require_gpu:
        parser.error("GPU phase must pass --require-gpu")

    root = args.output_dir.resolve()
    if root.exists():
        parser.error(f"output directory already exists: {root}")
    root.mkdir(parents=True)
    rows: list[dict[str, str]] = []

    for dt in sorted(set(args.delta_ts), reverse=True):
        substeps = args.substeps
        if args.substep_scaling == "inverse-deltaT":
            substeps = round(args.substeps * max(args.delta_ts) / dt)
        case = root / (
            f"{args.backend}_{args.scheme}_N{args.cells}_n{substeps}_dt{dt:g}"
        )
        prepare_case(case, dt, args.cells, args.end_time, substeps, args.scheme)
        started = time.monotonic()
        result = run_case(case, args.require_gpu)
        row = {
            "case": case.name,
            "backend_requested": args.backend,
            "model": "monodomainFDAManufacturedBatched",
            "coupling": args.scheme,
            "ddt_Vm": "Euler" if args.scheme == "godunov" else "backward",
            "deltaT_s": f"{dt:.16g}",
            "substeps": str(substeps),
            "substep_scaling": args.substep_scaling,
            "ionic_dt_s": f"{dt / substeps:.16g}",
            "cells_per_direction": str(args.cells),
            "cell_count": str(args.cells * args.cells),
            "end_time_s": f"{args.end_time:.16g}",
            "wall_s": f"{time.monotonic() - started:.6g}",
            **result,
        }
        rows.append(row)

        if len(rows) > 1 and all(
            rows[-i].get("status") == "completed" for i in (1, 2)
        ):
            coarse = rows[-2]
            fine = rows[-1]
            ratio = float(coarse["deltaT_s"]) / float(fine["deltaT_s"])
            for field in FIELDS:
                for norm in NORMS:
                    error_key = f"{field}_{norm}"
                    try:
                        error_ratio = (
                            float(coarse[error_key]) / float(fine[error_key])
                        )
                        fine[f"{field}_{norm}_order"] = (
                            f"{math.log(error_ratio) / math.log(ratio):.6g}"
                        )
                    except (KeyError, ValueError, ZeroDivisionError):
                        fine[f"{field}_{norm}_order"] = ""

        with (root / "summary.csv").open("w", newline="") as stream:
            fields = sorted({key for item in rows for key in item
                             if key != "log_tail"})
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows({key: item.get(key, "") for key in fields}
                             for item in rows)
        print(
            f"{case.name}: {row['status']} ({row.get('backend', '')}) "
            f"Vm_L2={row.get('Vm_L2', '')}",
            flush=True,
        )
        if row.get("log_tail"):
            print(row["log_tail"], flush=True)

    coarse_window = rows[:4]
    try:
        vm_l2 = [float(row["Vm_L2"]) for row in coarse_window]
        vm_l2_orders = [float(row["Vm_L2_order"]) for row in coarse_window[1:]]
        errors_present = all(row.get("status") == "completed"
                             for row in coarse_window)
    except (KeyError, ValueError):
        vm_l2 = []
        vm_l2_orders = []
        errors_present = False
    order_bounds = (0.8, 1.2) if args.scheme == "godunov" else (1.8, 2.3)
    checked_fields = ("Vm",) if args.scheme == "godunov" else FIELDS
    field_gates = {}
    for field in checked_fields:
        try:
            errors = [float(row[f"{field}_L2"]) for row in coarse_window]
            orders = [float(row[f"{field}_L2_order"])
                      for row in coarse_window[1:]]
            field_gates[field] = (
                all(errors[i] < errors[i - 1] for i in range(1, 4))
                and len(orders) == 3
                and all(order_bounds[0] <= order <= order_bounds[1]
                        for order in orders)
            )
        except (KeyError, ValueError):
            field_gates[field] = False
    accepted = errors_present and all(field_gates.values())
    (root / "acceptance.txt").write_text(
        f"criterion: on the four coarsest N=640 levels, {args.scheme} "
        f"L2 errors for {', '.join(checked_fields)} decrease monotonically "
        f"and each of the three halving orders is in "
        f"[{order_bounds[0]}, {order_bounds[1]}].\n"
        "The dt=0.0015625 point is reported but excluded from the order gate "
        "because it may approach the fixed-mesh spatial floor.\n"
        + "\n".join(
            f"{field}_L2_orders=" + ",".join(
                row.get(f"{field}_L2_order", "") for row in coarse_window[1:]
            )
            for field in checked_fields
        ) + "\n"
        f"result={'PASS' if accepted else 'NOT PASS'}\n"
    )
    return 0 if accepted and all(
        row.get("status") == "completed" for row in rows
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
