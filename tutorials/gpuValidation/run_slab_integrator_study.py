#!/usr/bin/env python3
"""Sweep tissue dt and ionic substeps for GPU ionic integrators."""

import argparse
import csv
import subprocess
import time
from pathlib import Path

from prepare_slab_case import DEFAULTS
from run_slab_matrix import MODELS


INTEGRATORS = ("euler", "rushLarsen")


def call(command, cwd=None):
    return subprocess.run(command, cwd=cwd, capture_output=True, text=True)


def prepare(generator, model, backend, case, end_time, delta_t, substeps=None,
            integrator="rushLarsen"):
    command = ["python3", str(generator), model, backend, str(case),
               "--end-time", str(end_time), "--delta-t", str(delta_t)]
    if backend == "batched":
        command.extend(("--substeps", str(substeps), "--integrator", integrator))
    return call(command)


def run_case(case, require_gpu):
    started = time.monotonic()
    result = call([str(case / "Allrun")], cwd=case)
    elapsed = time.monotonic() - started
    log_path = case / "log.cardiacFoam"
    log = log_path.read_text(errors="replace") if log_path.exists() else ""
    if result.returncode or "End" not in log:
        return "failed", elapsed, (result.stderr + log[-3000:])[-5000:]
    selected_gpu = "using CUDA device" in log
    if require_gpu and not selected_gpu:
        return "GPU unavailable", elapsed, log[-2000:]
    return ("completed-GPU" if selected_gpu else "completed-CPU"), elapsed, ""


def compare(comparator, first, second, end_time):
    result = call([
        "python3", str(comparator), str(first), str(second), "--times",
        "0.005", "0.01", str(end_time),
    ])
    return result.stdout.strip() if result.returncode == 0 else result.stderr.strip()


def save(rows, root):
    if not rows:
        return
    with (root / "summary.csv").open("w", newline="") as output:
        fieldnames = sorted({field for row in rows for field in row})
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS))
    parser.add_argument("--steps", nargs="+", type=int, default=[1, 5, 10, 25])
    parser.add_argument("--dt-multipliers", nargs="+", type=float, default=[1, 2, 5])
    parser.add_argument(
        "--integrators", nargs="+", choices=INTEGRATORS,
        default=["rushLarsen"],
        help="default focuses on Rush-Larsen gates plus Euler for remaining states",
    )
    parser.add_argument("--end-time", type=float, default=0.015)
    parser.add_argument("--require-gpu", action="store_true")
    args = parser.parse_args()
    if (args.end_time <= 0 or not args.steps or min(args.steps) < 1
            or not args.dt_multipliers or min(args.dt_multipliers) <= 0):
        parser.error("end-time, substeps, and dt multipliers must be positive")
    root = args.output_root.resolve()
    if root.exists():
        parser.error(f"output root already exists: {root}")
    root.mkdir(parents=True)
    here = Path(__file__).resolve().parent
    generator = here / "prepare_slab_case.py"
    comparator = here / "compare_slab.py"
    rows = []

    for model in args.models:
        model_root = root / model
        model_root.mkdir()
        _, base_dt, base_steps = DEFAULTS[model]
        fine_scalar = model_root / "scalar_fine_dt"
        made = prepare(generator, model, "scalar", fine_scalar,
                       args.end_time, base_dt)
        if made.returncode:
            raise RuntimeError(f"could not prepare fine scalar case: {made.stderr}")
        fine_status, fine_runtime, fine_detail = run_case(fine_scalar, False)
        if fine_status == "failed":
            rows.append({"model": model, "tissue_dt_s": base_dt,
                         "integrator": "scalar_RKF45", "substeps": "",
                         "status": fine_status, "runtime_s": fine_runtime,
                         "matched_RKF45": fine_detail, "fine_RKF45": fine_detail,
                         "fine_GPU_RL": "not run"})
            save(rows, root)
            continue

        # Matched scalar cases isolate ionic integration at each tissue dt.
        scalar_by_dt = {}
        for multiplier in sorted(set(args.dt_multipliers)):
            dt = base_dt*multiplier
            if multiplier == 1:
                scalar_by_dt[multiplier] = (fine_scalar, fine_runtime)
                continue
            case = model_root / f"scalar_RKF45_dt{dt:.0e}"
            made = prepare(generator, model, "scalar", case, args.end_time, dt)
            if made.returncode:
                scalar_by_dt[multiplier] = (None, "")
                rows.append({"model": model, "tissue_dt_s": dt,
                             "integrator": "scalar_RKF45", "substeps": "",
                             "status": "generation failed", "runtime_s": "",
                             "matched_RKF45": made.stderr.strip(),
                             "fine_RKF45": "not run", "fine_GPU_RL": "not run"})
                continue
            status, runtime, detail = run_case(case, False)
            scalar_by_dt[multiplier] = (case, runtime) if status != "failed" else (None, runtime)
            rows.append({"model": model, "tissue_dt_s": dt,
                         "integrator": "scalar_RKF45", "substeps": "",
                         "status": status, "runtime_s": runtime,
                         "matched_RKF45": "reference", "fine_RKF45":
                         compare(comparator, fine_scalar, case, args.end_time)
                         if status != "failed" else detail,
                         "fine_GPU_RL": "not run"})
            save(rows, root)

        # Fine batched GPU baseline is RL at five nominal ionic substeps.
        fine_gpu = model_root / f"gpu_rushLarsen_{base_steps}substeps_dt{base_dt:.0e}"
        made = prepare(generator, model, "batched", fine_gpu,
                       args.end_time, base_dt, base_steps, "rushLarsen")
        if made.returncode:
            fine_gpu_status, fine_gpu_runtime, fine_gpu_detail = (
                "generation failed", "", made.stderr.strip()
            )
        else:
            fine_gpu_status, fine_gpu_runtime, fine_gpu_detail = run_case(
                fine_gpu, args.require_gpu
            )
        fine_gpu_ok = fine_gpu_status == "completed-GPU"

        for multiplier in sorted(set(args.dt_multipliers)):
            dt = base_dt*multiplier
            matched_scalar, scalar_runtime = scalar_by_dt[multiplier]
            for integrator in args.integrators:
                for substeps in args.steps:
                    if multiplier == 1 and integrator == "rushLarsen" and substeps == base_steps:
                        case = fine_gpu
                        status, runtime, detail = (
                            fine_gpu_status, fine_gpu_runtime, fine_gpu_detail
                        )
                    else:
                        case = model_root / (
                            f"gpu_{integrator}_{substeps}substeps_"
                            f"dt{dt:.0e}"
                        )
                        made = prepare(generator, model, "batched", case,
                                       args.end_time, dt, substeps, integrator)
                        if made.returncode:
                            status, runtime, detail = (
                                "generation failed", "", made.stderr.strip()
                            )
                        else:
                            status, runtime, detail = run_case(
                                case, args.require_gpu
                            )

                    good = status == "completed-GPU"
                    if not good:
                        (model_root / f"{case.name}.failure.txt").write_text(
                            detail + "\n"
                        )
                    matched = (
                        compare(comparator, matched_scalar, case, args.end_time)
                        if good and matched_scalar is not None else "matched scalar unavailable"
                    )
                    fine_scalar_cmp = (
                        compare(comparator, fine_scalar, case, args.end_time)
                        if good else "not run: GPU case did not complete"
                    )
                    fine_gpu_cmp = (
                        compare(comparator, fine_gpu, case, args.end_time)
                        if good and fine_gpu_ok else "fine GPU baseline unavailable"
                    )
                    rows.append({
                        "model": model, "tissue_dt_s": dt,
                        "integrator": integrator, "substeps": substeps,
                        "status": status, "runtime_s": runtime,
                        "matched_RKF45": matched,
                        "fine_RKF45": fine_scalar_cmp,
                        "fine_GPU_RL": fine_gpu_cmp,
                        "matched_scalar_runtime_s": scalar_runtime,
                    })
                    save(rows, root)
                    print(f"{model}/dt={dt:g}/{integrator}/{substeps}: "
                          f"{status} ({runtime}s)", flush=True)
    save(rows, root)
    print(f"Wrote {root / 'summary.csv'}")


if __name__ == "__main__":
    main()
