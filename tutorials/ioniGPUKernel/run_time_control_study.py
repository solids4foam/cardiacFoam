#!/usr/bin/env python3
"""Run a small matched CPU/GPU slab study of tissue dt and ionic subcycling."""

import argparse
import csv
import subprocess
import time
from pathlib import Path


CONTROLS = (
    ("baseline", 2e-6, 5),
    ("same_tissue_dt_one_substep", 2e-6, 1),
    ("double_tissue_dt_same_ionic_dt", 4e-6, 10),
)


def run_case(generator, model, backend, case, delta_t, substeps, end_time):
    command = [
        "python3", str(generator), model, backend, str(case),
        "--delta-t", str(delta_t), "--end-time", str(end_time),
    ]
    if backend == "batched":
        command.extend(("--substeps", str(substeps)))
    made = subprocess.run(command, capture_output=True, text=True)
    if made.returncode:
        raise RuntimeError(f"case generation failed for {case}: {made.stderr}")

    started = time.monotonic()
    result = subprocess.run(
        [str(case / "Allrun")], cwd=case, capture_output=True, text=True
    )
    elapsed = time.monotonic() - started
    if result.returncode:
        raise RuntimeError(
            f"case failed for {case} ({result.returncode}); "
            f"see {case / 'log.cardiacFoam'}"
        )
    log = (case / "log.cardiacFoam").read_text(errors="replace")
    if backend == "batched" and "using CUDA device" not in log:
        raise RuntimeError(f"GPU required but not selected for {case}")
    return elapsed


def compare(comparator, reference, candidate, end_time):
    result = subprocess.run(
        ["python3", str(comparator), str(reference), str(candidate),
         "--times", str(end_time)],
        capture_output=True, text=True,
    )
    if result.returncode:
        raise RuntimeError(
            f"comparison failed for {reference} and {candidate}: {result.stderr}"
        )
    return result.stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--model", default="TNNP")
    parser.add_argument("--end-time", type=float, default=0.015)
    args = parser.parse_args()
    if args.end_time <= 0:
        parser.error("end-time must be positive")

    here = Path(__file__).resolve().parent
    generator = here / "prepare_slab_case.py"
    comparator = here / "slab2D/setup/post_processing_slab_field_parity.py"
    output = args.output_dir.resolve()
    if output.exists():
        parser.error(f"output directory already exists: {output}")
    output.mkdir(parents=True)

    rows = []
    cases = {}
    for label, delta_t, substeps in CONTROLS:
        cases[label] = {}
        for backend in ("scalar", "batched"):
            case = output / label / backend
            case.parent.mkdir(parents=True, exist_ok=True)
            elapsed = run_case(
                generator, args.model, backend, case, delta_t, substeps,
                args.end_time,
            )
            cases[label][backend] = case
            rows.append({
                "control": label,
                "backend": backend,
                "tissue_deltaT_s": delta_t,
                "batched_substeps": substeps if backend == "batched" else "",
                "runtime_s": elapsed,
                "gpu_backend": backend == "batched",
            })
            print(f"{label}/{backend}: completed in {elapsed:.3f} s", flush=True)

    comparisons = []
    baseline = cases["baseline"]["batched"]
    for label, _, _ in CONTROLS:
        matched = compare(
            comparator, cases[label]["scalar"], cases[label]["batched"],
            args.end_time,
        )
        comparisons.append((label + "_scalar_vs_gpu", matched))
        if label != "baseline":
            control_effect = compare(comparator, baseline, cases[label]["batched"],
                                     args.end_time)
            comparisons.append((label + "_gpu_vs_baseline", control_effect))

    with (output / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "comparisons.txt").write_text(
        "\n".join(f"{label}: {result}" for label, result in comparisons) + "\n"
    )
    print(f"Wrote {output / 'summary.csv'} and {output / 'comparisons.txt'}")


if __name__ == "__main__":
    main()
