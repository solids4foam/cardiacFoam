#!/usr/bin/env python3
"""Run matching scalar and CUDA manufactured-monodomain smoke cases."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "tutorials/manufacturedSolutions/monodomainPseudoECG"


def prepare_case(destination: Path, model: str, delta_t: float, end_time: float,
                 cells: int, substeps: int) -> None:
    (destination / "constant").mkdir(parents=True)
    (destination / "system").mkdir()
    for name in ("electroProperties", "physicsProperties"):
        shutil.copy2(TEMPLATE / "constant" / name, destination / "constant")
    for name in ("controlDict", "fvSchemes", "fvSolution"):
        shutil.copy2(TEMPLATE / "system" / name, destination / "system")
    shutil.copy2(TEMPLATE / "system/blockMeshDict.3D",
                 destination / "system/blockMeshDict")

    electro_path = destination / "constant/electroProperties"
    electro = electro_path.read_text()
    electro = re.sub(r"(?m)^(\s*ionicModel\s+)\w+(\s*;)",
                     rf"\g<1>{model}\g<2>", electro, count=1)
    if model.endswith("Batched"):
        electro = re.sub(
            r"(?m)^(\s*ionicModel\s+\w+;)",
            rf"\1\n    batchedIntegrator euler;\n    batchedSubsteps {substeps};",
            electro,
            count=1,
        )
    electro_path.write_text(electro)

    control_path = destination / "system/controlDict"
    control = control_path.read_text()
    control = re.sub(r"(?m)^endTime\s+[^;]+;", f"endTime {end_time};", control)
    control = re.sub(r"(?m)^deltaT\s+[^;]+;", f"deltaT {delta_t};", control)
    control = re.sub(r"(?m)^writeInterval\s+[^;]+;",
                     f"writeInterval {end_time};", control)
    control = re.sub(r"(?m)^writeFormat\s+[^;]+;", "writeFormat ascii;", control)
    control_path.write_text(control)

    mesh_path = destination / "system/blockMeshDict"
    mesh = mesh_path.read_text()
    mesh, substitutions = re.subn(
        r"(?m)^(\s*hex\s*\([^)]*\)\s*\()\d+\s+\d+\s+\d+(\)\s+simpleGrading)",
        rf"\g<1>{cells} {cells} {cells}\g<2>",
        mesh,
        count=1,
    )
    if substitutions != 1:
        raise RuntimeError(f"Could not set 3D mesh resolution in {mesh_path}")
    mesh_path.write_text(mesh)


def run_case(case: Path, label: str, require_gpu: bool) -> dict[str, str]:
    for command in ("blockMesh", "cardiacFoam"):
        with (case / f"log.{command}").open("w") as log:
            result = subprocess.run(
                [command, "-case", str(case)],
                cwd=case,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if result.returncode:
            tail = (case / f"log.{command}").read_text(errors="replace")[-3000:]
            return {"case": label, "status": f"{command} exit {result.returncode}",
                    "backend": "", "runtime_s": "", "log_tail": tail}

    run_log = (case / "log.cardiacFoam").read_text(errors="replace")
    backend = "GPU" if "using CUDA device" in run_log else "scalar/CPU"
    if require_gpu and label == "gpu" and backend != "GPU":
        return {"case": label, "status": "GPU not selected", "backend": backend,
                "runtime_s": "", "log_tail": run_log[-3000:]}

    error_lines = {line.split()[0]: line.split()[1:4]
                   for line in run_log.splitlines()
                   if line.split() and line.split()[0] in {"Vm", "u1", "u2"}
                   and len(line.split()) >= 4}
    if set(error_lines) != {"Vm", "u1", "u2"}:
        return {"case": label, "status": "missing manufactured Vm error",
                "backend": backend, "runtime_s": "", "log_tail": run_log[-3000:]}
    errors = {f"{field}_{norm}": value
              for field, values in error_lines.items()
              for norm, value in zip(("L1", "L2", "Linf"), values)}
    runtime_lines = [line for line in run_log.splitlines()
                     if line.startswith("ExecutionTime =")]
    runtime = runtime_lines[-1].split("=", 1)[1].split()[0] if runtime_lines else ""
    return {"case": label, "status": "completed", "backend": backend,
            "runtime_s": runtime, **errors, "log_tail": ""}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--delta-t", type=float, default=1.0e-5)
    parser.add_argument("--end-time", type=float, default=5.0e-5)
    parser.add_argument("--cells", type=int, default=8)
    parser.add_argument("--substeps", type=int, default=5)
    parser.add_argument("--require-gpu", action="store_true")
    args = parser.parse_args()
    if args.delta_t <= 0 or args.end_time <= 0 or args.cells < 2 or args.substeps < 1:
        parser.error("deltaT/endTime/substeps must be positive and cells >= 2")

    output = args.output_dir.resolve()
    if output.exists():
        parser.error(f"output directory already exists: {output}")
    output.mkdir(parents=True)
    rows = []
    for label, model in (("scalar", "monodomainFDAManufactured"),
                         ("gpu", "monodomainFDAManufacturedBatched")):
        case = output / label
        prepare_case(case, model, args.delta_t, args.end_time,
                     args.cells, args.substeps)
        rows.append(run_case(case, label, args.require_gpu))

    with (output / "summary.csv").open("w", newline="") as stream:
        fieldnames = ("case", "status", "backend", "runtime_s",
                      *(f"{field}_{norm}" for field in ("Vm", "u1", "u2")
                        for norm in ("L1", "L2", "Linf")))
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fieldnames}
                         for row in rows)
    for row in rows:
        print(f"{row['case']}: {row['status']} ({row['backend']}); Vm L2={row['Vm_L2']}")
        if row.get("log_tail"):
            print(row["log_tail"])
    return 0 if all(row["status"] == "completed" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
