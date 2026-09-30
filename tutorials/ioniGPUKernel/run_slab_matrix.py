#!/usr/bin/env python3
"""Run matched scalar/batched CPU slab screens for selected ionic models."""

import argparse
import concurrent.futures
import csv
import subprocess
import time
from pathlib import Path


MODELS = (
    "AlievPanfilov", "BuenoOrovio", "Courtemanche", "Fabbri", "Gaur",
    "Grandi", "PerisYague", "Stewart", "TNNP", "TWorld",
    "ToRORd_dynCl", "Trovato",
)
TISSUE = {
    "AlievPanfilov": "myocyte", "BuenoOrovio": "epicardialCells",
    "Courtemanche": "myocyte", "Fabbri": "myocyte", "Gaur": "myocyte",
    "Grandi": "myocyte", "PerisYague": "myocyte", "Stewart": "myocyte",
    "TNNP": "epicardialCells", "TWorld": "epicardialCells",
    "ToRORd_dynCl": "epicardialCells", "Trovato": "myocyte",
}


def run_pair(model, root, end_time, require_gpu):
    here = Path(__file__).resolve().parent
    generator = here / "prepare_slab_case.py"
    compare = here / "slab2D/setup/post_processing_slab_field_parity.py"
    cases = {backend: root / model / backend for backend in ("scalar", "batched")}
    results = {}
    for backend, case in cases.items():
        created = subprocess.run(
            ["python3", str(generator), model, backend, str(case),
             "--end-time", str(end_time)],
            capture_output=True, text=True,
        )
        if created.returncode:
            results[backend] = {"status": "case generation failed", "detail": created.stderr.strip()}
            continue
        start = time.monotonic()
        completed = subprocess.run(
            [str(case / "Allrun")], cwd=case, capture_output=True, text=True,
        )
        elapsed = time.monotonic() - start
        log_path = case / "log.cardiacFoam"
        log = log_path.read_text(errors="replace") if log_path.exists() else ""
        results[backend] = {
            "status": "completed" if completed.returncode == 0 else f"exit {completed.returncode}",
            "backend": "GPU" if "using CUDA device" in log else ("CPU" if backend == "batched" else "scalar"),
            "runtime_s": elapsed,
            "detail": "",
        }
        if (
            backend == "batched"
            and require_gpu
            and results[backend]["status"] == "completed"
            and results[backend]["backend"] != "GPU"
        ):
            results[backend]["status"] = "GPU unavailable"
    comparison = "not run"
    if all(results[b]["status"] == "completed" for b in cases):
        result = subprocess.run(
            ["python3", str(compare), str(cases["scalar"]), str(cases["batched"]),
             "--times", str(end_time)],
            capture_output=True, text=True,
        )
        comparison = result.stdout.strip() if result.returncode == 0 else result.stderr.strip()
    print(f"{model}: scalar={results['scalar']['status']} "
          f"batched={results['batched']['status']} comparison={comparison}", flush=True)
    return {
        "model": model,
        "tissue": TISSUE[model],
        "scalar_status": results["scalar"]["status"],
        "batched_status": results["batched"]["status"],
        "batched_backend": results["batched"].get("backend", ""),
        "scalar_runtime_s": results["scalar"].get("runtime_s", ""),
        "batched_runtime_s": results["batched"].get("runtime_s", ""),
        "comparison": comparison,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS))
    parser.add_argument("--end-time", type=float, default=0.005)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--require-gpu", action="store_true")
    args = parser.parse_args()
    root = args.output_root.resolve()
    if root.exists():
        parser.error(f"output root already exists: {root}")
    if args.end_time <= 0 or args.workers < 1:
        parser.error("end-time and workers must be positive")
    root.mkdir(parents=True)
    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [
            pool.submit(run_pair, model, root, args.end_time, args.require_gpu)
            for model in args.models
        ]
        for future in concurrent.futures.as_completed(futures):
            rows.append(future.result())
    rows.sort(key=lambda row: MODELS.index(row["model"]))
    with (root / "summary.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {root / 'summary.csv'}")


if __name__ == "__main__":
    main()
