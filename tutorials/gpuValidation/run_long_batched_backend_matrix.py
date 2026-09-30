#!/usr/bin/env python3
"""Run long CPU-batched or CUDA-batched cases for every in-scope model."""

import argparse
import csv
import re
import subprocess
import time
from pathlib import Path


MODELS = (
    "AlievPanfilov", "BuenoOrovio", "Courtemanche", "Gaur", "Grandi",
    "PerisYague", "Stewart", "TNNP", "TWorld", "ToRORd_dynCl", "Trovato",
)


def set_entry(path: Path, name: str, value: str) -> None:
    text = path.read_text()
    pattern = rf"(?m)^(\s*{re.escape(name)}\s+).+?;\s*$"
    changed, count = re.subn(pattern, rf"\g<1>{value};", text, count=1)
    if count != 1:
        raise ValueError(f"could not set {name} in {path}")
    path.write_text(changed)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geometry", choices=("slab", "niederer"))
    parser.add_argument("backend", choices=("cpu", "cuda"))
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS))
    parser.add_argument("--delta-t", type=float, default=1e-5)
    parser.add_argument("--substeps", type=int, default=25)
    parser.add_argument("--end-time", type=float, default=0.5)
    parser.add_argument("--scheme", choices=("godunov", "sbdf2"), default="sbdf2")
    args = parser.parse_args()
    if min(args.delta_t, args.end_time) <= 0 or args.substeps < 1:
        parser.error("deltaT/endTime must be positive and substeps >= 1")

    root = args.output_root.resolve()
    if root.exists():
        parser.error(f"output root already exists: {root}")
    root.mkdir(parents=True)
    here = Path(__file__).resolve().parent
    rows = []
    overall = 0
    for model in args.models:
        case = root / model
        if args.geometry == "slab":
            generator = here / "prepare_slab_case.py"
            command = [
                "python3", str(generator), model, "batched", str(case),
                "--delta-t", str(args.delta_t), "--substeps", str(args.substeps),
                "--integrator", "rushLarsen", "--time-coupling-scheme", args.scheme,
                "--end-time", str(args.end_time), "--write-precision", "15",
            ]
            output_interval = 1e-3
        else:
            generator = here / "prepare_niederer_case.py"
            case_backend = "gpu" if args.backend == "cuda" else "batched"
            command = [
                "python3", str(generator), model, case_backend, str(case),
                "--delta-t", str(args.delta_t), "--substeps", str(args.substeps),
                "--end-time", str(args.end_time),
                "--time-coupling-scheme", args.scheme,
            ]
            output_interval = 1e-3

        made = subprocess.run(command, capture_output=True, text=True)
        row = {
            "model": model,
            "geometry": args.geometry,
            "backend": args.backend,
            "scheme": args.scheme,
            "deltaT_s": args.delta_t,
            "substeps": args.substeps,
            "ionic_dt_s": args.delta_t / args.substeps,
            "endTime_s": args.end_time,
            "writeInterval_s": output_interval,
        }
        if made.returncode:
            row.update(status="case-generation-failed", runtime_s="")
            (root / f"{model}.generation-failure.txt").write_text(made.stderr)
            overall = 1
        else:
            control = case / "system/controlDict"
            set_entry(control, "writeInterval", str(output_interval))
            started = time.monotonic()
            with (case / "allrun.stdout").open("w") as log_file:
                result = subprocess.run([str(case / "Allrun")], cwd=case,
                                        stdout=log_file, stderr=subprocess.STDOUT)
            elapsed = time.monotonic() - started
            log_path = case / "log.cardiacFoam"
            log = log_path.read_text(errors="replace") if log_path.exists() else ""
            if result.returncode or "End" not in log:
                status = f"failed-exit-{result.returncode}"
                (root / f"{model}.failure.txt").write_text(log[-10000:])
                overall = 1
            elif args.backend == "cuda" and "using CUDA device" not in log:
                status = "failed-no-CUDA-device"
                overall = 1
            elif args.backend == "cpu" and "using CUDA device" in log:
                status = "failed-unexpected-CUDA"
                overall = 1
            else:
                status = "completed"
            row.update(status=status, runtime_s=f"{elapsed:.6f}")
        rows.append(row)
        with (root / "summary.csv").open("w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"{model}/{args.geometry}/{args.backend}: {row['status']} "
              f"runtime={row['runtime_s']} s", flush=True)
    print(f"Wrote {root / 'summary.csv'}", flush=True)
    return overall


if __name__ == "__main__":
    raise SystemExit(main())
