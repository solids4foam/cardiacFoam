#!/usr/bin/env python3
"""Compare scalar RKF45 with CUDA batched models on Niederer geometry."""

import argparse
import csv
import subprocess
import time
from pathlib import Path


MODELS = (
    "AlievPanfilov", "BuenoOrovio", "Courtemanche", "Gaur", "Grandi",
    "PerisYague", "Stewart", "TNNP", "TWorld", "ToRORd_dynCl", "Trovato",
)


def run(command, cwd=None):
    return subprocess.run(command, cwd=cwd, capture_output=True, text=True)


def run_case(generator, model, backend, case, substeps, delta_t, end_time):
    make_args = [
        "python3", str(generator), model, backend, str(case),
    ]
    if backend == "gpu":
        make_args.extend(("--substeps", str(substeps)))
    if delta_t is not None:
        make_args.extend(("--delta-t", str(delta_t)))
    if end_time is not None:
        make_args.extend(("--end-time", str(end_time)))
    made = run(make_args)
    if made.returncode:
        return "generation failed", "", made.stderr[-2000:]
    started = time.monotonic()
    result = run([str(case / "Allrun")], cwd=case)
    elapsed = time.monotonic() - started
    log_path = case / "log.cardiacFoam"
    log = log_path.read_text(errors="replace") if log_path.exists() else ""
    if result.returncode or "End" not in log:
        return f"failed (exit {result.returncode})", elapsed, (result.stderr + log[-2000:])[-4000:]
    if backend == "gpu" and "using CUDA device" not in log:
        return "GPU unavailable (CPU fallback)", elapsed, log[-2000:]
    return "completed", elapsed, ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS))
    parser.add_argument("--substeps", type=int, default=25)
    parser.add_argument("--delta-t", type=float, default=1e-5)
    parser.add_argument("--end-time", type=float, default=0.015)
    args = parser.parse_args()
    if (args.substeps < 1 or (args.delta_t is not None and args.delta_t <= 0)
            or (args.end_time is not None and args.end_time <= 0)):
        parser.error("substeps, delta-t, and end-time must be positive")
    root = args.output_root.resolve()
    if root.exists():
        parser.error(f"output root already exists: {root}")
    root.mkdir(parents=True)
    here = Path(__file__).resolve().parent
    generator = here / "prepare_niederer_case.py"
    comparator = here / "compare_slab.py"
    rows = []

    for model in args.models:
        cases = {backend: root / model / backend for backend in ("scalar", "gpu")}
        result_row = {"model": model, "deltaT_s": args.delta_t,
                      "endTime_s": args.end_time, "substeps": args.substeps}
        statuses = {}
        for backend, case in cases.items():
            status, elapsed, detail = run_case(
                generator, model, backend, case, args.substeps,
                args.delta_t, args.end_time,
            )
            statuses[backend] = status
            result_row[f"{backend}_status"] = status
            result_row[f"{backend}_runtime_s"] = elapsed
            if detail:
                (root / f"{model}_{backend}.failure.txt").write_text(detail + "\n")
            print(f"{model}/{backend}: {status} ({elapsed}s)", flush=True)

        if all(statuses[b] == "completed" for b in ("scalar", "gpu")):
            comparison = run([
                "python3", str(comparator), str(cases["scalar"]),
                str(cases["gpu"]), "--times",
                *[str(t) for t in sorted(set((
                    min(0.005, args.end_time or 0.015),
                    min(0.01, args.end_time or 0.015),
                    args.end_time or 0.015,
                )))],
            ])
            result_row["comparison"] = (
                comparison.stdout.strip() if comparison.returncode == 0
                else comparison.stderr.strip()
            )
        else:
            result_row["comparison"] = "not run: backend did not complete"
        rows.append(result_row)
        with (root / "summary.csv").open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"{model}: {result_row['comparison']}", flush=True)

    print(f"Wrote {root / 'summary.csv'}")


if __name__ == "__main__":
    main()
