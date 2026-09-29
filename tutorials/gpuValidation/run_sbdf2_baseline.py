#!/usr/bin/env python3
"""Compare Godunov and SBDF2 coupling for matched scalar and CUDA slab cases."""

import argparse
import csv
import subprocess
import time
from pathlib import Path

from run_slab_matrix import MODELS


def command(args, cwd=None):
    return subprocess.run(args, cwd=cwd, capture_output=True, text=True)


def run_case(generator, model, backend, scheme, case, delta_t, substeps,
             end_time, require_gpu, write_precision):
    prepare_args = [
        "python3", str(generator), model, backend, str(case),
        "--delta-t", str(delta_t), "--end-time", str(end_time),
        "--time-coupling-scheme", scheme,
    ]
    if write_precision is not None:
        prepare_args.extend(("--write-precision", str(write_precision)))
    if backend == "batched":
        prepare_args.extend(("--integrator", "rushLarsen",
                             "--substeps", str(substeps)))
    generated = command(prepare_args)
    if generated.returncode:
        return "generation failed", "", generated.stderr[-3000:]

    started = time.monotonic()
    result = command([str(case / "Allrun")], cwd=case)
    elapsed = time.monotonic() - started
    log_path = case / "log.cardiacFoam"
    log = log_path.read_text(errors="replace") if log_path.exists() else ""
    if result.returncode or "End" not in log:
        return "failed", elapsed, (result.stderr + log[-3000:])[-5000:]
    cuda = "using CUDA device" in log
    if backend == "batched" and require_gpu and not cuda:
        return "GPU unavailable", elapsed, log[-2500:]
    return ("completed-GPU" if cuda else "completed-CPU"), elapsed, ""


def compare(comparator, first, second, end_time):
    result = command([
        "python3", str(comparator), str(first), str(second),
        "--times", "0.005", "0.01", str(end_time),
    ])
    return (result.stdout.strip() if result.returncode == 0
            else "COMPARISON_FAILED: " + result.stderr.strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--models", nargs="+", choices=MODELS,
                        default=["TNNP", "TWorld"])
    parser.add_argument("--delta-t", type=float, default=2e-6)
    parser.add_argument("--substeps", type=int, default=5)
    parser.add_argument("--end-time", type=float, default=0.015)
    parser.add_argument("--require-gpu", action="store_true")
    parser.add_argument("--write-precision", type=int)
    args = parser.parse_args()
    if args.delta_t <= 0 or args.substeps < 1 or args.end_time < 0.01:
        parser.error("delta-t and substeps must be positive; end-time must be >= 0.01")
    root = args.output_root.resolve()
    if root.exists():
        parser.error(f"output root already exists: {root}")
    root.mkdir(parents=True)
    here = Path(__file__).resolve().parent
    generator = here / "prepare_slab_case.py"
    comparator = here / "compare_slab.py"
    rows = []

    for model in args.models:
        (root / model).mkdir()
        cases = {
            (scheme, backend): root / model / f"{scheme}_{backend}"
            for scheme in ("godunov", "sbdf2")
            for backend in ("scalar", "batched")
        }
        statuses = {}
        for (scheme, backend), case in cases.items():
            status, runtime, detail = run_case(
                generator, model, backend, scheme, case, args.delta_t,
                args.substeps, args.end_time, args.require_gpu,
                args.write_precision,
            )
            statuses[scheme, backend] = status
            rows.append({"model": model, "scheme": scheme, "backend": backend,
                         "status": status, "runtime_s": runtime})
            if detail:
                (root / model / f"{scheme}_{backend}.failure.txt").write_text(
                    detail + "\n"
                )
            print(f"{model}/{scheme}/{backend}: {status} ({runtime}s)", flush=True)
            with (root / "summary.csv").open("w", newline="") as output:
                writer = csv.DictWriter(output, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

        pairs = (
            ("scalar_godunov_vs_scalar_sbdf2", "godunov", "scalar", "sbdf2", "scalar"),
            ("gpu_godunov_vs_gpu_sbdf2", "godunov", "batched", "sbdf2", "batched"),
            ("scalar_godunov_vs_gpu_godunov", "godunov", "scalar", "godunov", "batched"),
            ("scalar_sbdf2_vs_gpu_sbdf2", "sbdf2", "scalar", "sbdf2", "batched"),
        )
        for name, scheme_a, backend_a, scheme_b, backend_b in pairs:
            if (statuses[scheme_a, backend_a].startswith("completed")
                    and statuses[scheme_b, backend_b].startswith("completed")):
                output = compare(comparator, cases[scheme_a, backend_a],
                                 cases[scheme_b, backend_b], args.end_time)
            else:
                output = "not run: a required case did not complete"
            (root / model / f"{name}.txt").write_text(output + "\n")
            print(f"{model}/{name}:\n{output}", flush=True)

    print(f"Wrote {root / 'summary.csv'}")


if __name__ == "__main__":
    main()
