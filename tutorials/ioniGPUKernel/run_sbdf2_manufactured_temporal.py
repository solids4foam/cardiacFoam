#!/usr/bin/env python3
"""Run the backup's coarse 2D manufactured SBDF2 temporal ladder."""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "tutorials/manufacturedSolutions/monodomainPseudoECG"
DT_DEFAULTS = (0.025, 0.0125, 0.00625, 0.003125, 0.0015625)
ERROR_RE = re.compile(
    r"^\s*(Vm|u1|u2)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)",
    re.MULTILINE,
)


def electro_properties() -> str:
    return r'''FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    location "constant";
    object electroProperties;
}

myocardiumSolver monodomainSolver;

monodomainSolverCoeffs
{
    conductivity [-1 -3 3 0 0 2 0] (0.111453302 0 0 0.121585420 0 0.030396355);
    chi [0 -1 0 0 0 0 0] 3;
    cm [-1 -4 4 0 0 2 0] 2;

    ionicModel monodomainFDAManufactured;
    dimension "2D";
    solutionAlgorithm implicit;
    timeCouplingScheme sbdf2;
    solver RKF45;
    initialODEStep 1e-5;
    maxSteps 1000000000;

    outputVariables
    {
        ionic
        {
            export (u1 u2 u3);
            debug (V u1 u2 u3 Iion);
        }
    }

    verificationModel
    {
        type manufacturedFDAMonodomainVerifier;
    }

    externalStimulus
    {
        stimulusLocationMin (0 0 0);
        stimulusLocationMax (0 0 0);
        stimulusDuration [0 0 1 0 0 0 0] 0;
        stimulusIntensity [0 -3 0 0 0 1 0] 0;
        stimulusStartTime 0;
    }
}
'''


def prepare_case(case: Path, dt: float, cells: int, end_time: float) -> None:
    (case / "constant").mkdir(parents=True)
    (case / "system").mkdir()
    shutil.copy2(TEMPLATE / "constant/physicsProperties", case / "constant")
    shutil.copy2(TEMPLATE / "system/fvSolution", case / "system")
    shutil.copy2(TEMPLATE / "system/fvSchemes", case / "system")
    shutil.copy2(TEMPLATE / "system/blockMeshDict.2D", case / "system/blockMeshDict")
    (case / "constant/electroProperties").write_text(electro_properties())

    mesh_path = case / "system/blockMeshDict"
    mesh, count = re.subn(
        r"(?m)^(\s*hex\s*\([^)]*\)\s*\()\d+\s+\d+\s+\d+(\)\s+simpleGrading)",
        rf"\g<1>{cells} {cells} 1\g<2>", mesh_path.read_text(), count=1,
    )
    if count != 1:
        raise RuntimeError(f"Could not set 2D mesh size in {mesh_path}")
    mesh_path.write_text(mesh)

    schemes_path = case / "system/fvSchemes"
    schemes = schemes_path.read_text()
    schemes, count = re.subn(r"(?m)^\s*ddt\(Vm\)\s+\w+\s*;",
                             "    ddt(Vm) backward;", schemes, count=1)
    if count != 1:
        raise RuntimeError(f"Could not set backward Vm derivative in {schemes_path}")
    schemes_path.write_text(schemes)

    control = (TEMPLATE / "system/controlDict").read_text()
    replacements = {
        "endTime": end_time,
        "deltaT": dt,
        "writeInterval": end_time,
    }
    for key, value in replacements.items():
        control, count = re.subn(rf"(?m)^{key}\s+[^;]+;", f"{key} {value:.16g};",
                                 control, count=1)
        if count != 1:
            raise RuntimeError(f"Could not set {key} in controlDict")
    control = re.sub(r"(?m)^writeFormat\s+\w+;", "writeFormat ascii;", control)
    control = re.sub(r"(?m)^writePrecision\s+\d+;", "writePrecision 15;", control)
    (case / "system/controlDict").write_text(control)


def run_case(case: Path, dt: float, cells: int, end_time: float) -> dict[str, str]:
    for command in ("blockMesh", "cardiacFoam"):
        with (case / f"log.{command}").open("w") as log:
            result = subprocess.run([command, "-case", str(case)], cwd=case,
                                    stdout=log, stderr=subprocess.STDOUT,
                                    check=False)
        if result.returncode:
            return {"status": f"{command} exit {result.returncode}",
                    "log_tail": (case / f"log.{command}").read_text(
                        errors="replace")[-4000:]}

    log = (case / "log.cardiacFoam").read_text(errors="replace")
    errors = {field: values for field, *values in ERROR_RE.findall(log)}
    if set(errors) != {"Vm", "u1", "u2"}:
        return {"status": "missing manufactured errors",
                "log_tail": log[-4000:]}
    steps_match = re.search(r"Number of steps\s*=\s*(\d+)", log)
    final_match = re.search(r"Final simulation time\s*=\s*([0-9.eE+-]+)", log)
    runtime_matches = re.findall(r"ExecutionTime\s*=\s*([0-9.eE+-]+)", log)
    cpu_matches = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", log)
    row: dict[str, str] = {
        "status": "completed",
        "steps": steps_match.group(1) if steps_match else "",
        "final_time": final_match.group(1) if final_match else "",
        "runtime_s": runtime_matches[-1] if runtime_matches else "",
        "clock_s": cpu_matches[-1] if cpu_matches else "",
        "log_tail": "",
    }
    for field, values in errors.items():
        for norm, value in zip(("L1", "L2", "Linf"), values):
            row[f"{field}_{norm}"] = value
    row.update({"deltaT_s": f"{dt:.16g}", "cells_per_direction": str(cells),
                "cell_count": str(cells * cells), "end_time": f"{end_time:.16g}"})
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--delta-ts", nargs="+", type=float, default=DT_DEFAULTS)
    parser.add_argument("--cells", type=int, default=640)
    parser.add_argument("--end-time", type=float, default=0.2)
    args = parser.parse_args()
    if args.cells < 2 or args.end_time <= 0 or not args.delta_ts or min(args.delta_ts) <= 0:
        parser.error("cells/end-time/deltaT must be positive; cells >= 2")
    if any(not math.isclose(args.end_time / dt, round(args.end_time / dt),
                            rel_tol=0.0, abs_tol=1e-10) for dt in args.delta_ts):
        parser.error("end-time must be an integer multiple of every deltaT")

    root = args.output_dir.resolve()
    if root.exists():
        parser.error(f"output directory already exists: {root}")
    root.mkdir(parents=True)
    rows = []
    for dt in sorted(set(args.delta_ts), reverse=True):
        case = root / f"sbdf2_backward_N{args.cells}_dt{dt:g}"
        prepare_case(case, dt, args.cells, args.end_time)
        started = time.monotonic()
        row = run_case(case, dt, args.cells, args.end_time)
        row.update({"case": case.name, "wall_s": f"{time.monotonic() - started:.6g}",
                    "method": "implicit SBDF2 + backward ddt(Vm) + scalar RKF45"})
        rows.append(row)
        with (root / "summary.csv").open("w", newline="") as stream:
            fields = sorted({key for item in rows for key in item if key != "log_tail"})
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows({key: item.get(key, "") for key in fields} for item in rows)
        print(f"{case.name}: {row['status']} Vm_L2={row.get('Vm_L2', '')}", flush=True)
        if row.get("log_tail"):
            print(row["log_tail"], flush=True)

    previous = None
    for row in rows:
        if (previous is not None
                and previous.get("status") == "completed"
                and row.get("status") == "completed"
                and previous.get("Vm_L2") and row.get("Vm_L2")):
            error = float(row["Vm_L2"])
            previous_error = float(previous["Vm_L2"])
            previous_dt = float(previous["deltaT_s"])
            current_dt = float(row["deltaT_s"])
            observed_order = math.log(previous_error / error) / math.log(
                previous_dt / current_dt
            )
            row["Vm_L2_order_from_previous_dt"] = f"{observed_order:.6g}"
        else:
            row["Vm_L2_order_from_previous_dt"] = ""
        previous = row
    fields = sorted({key for item in rows for key in item if key != "log_tail"})
    with (root / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: item.get(key, "") for key in fields} for item in rows)

    clean_window = rows[:4]
    orders = []
    for row in clean_window[1:]:
        try:
            orders.append(float(row["Vm_L2_order_from_previous_dt"]))
        except (KeyError, ValueError):
            pass
    has_clean_errors = (len(clean_window) == 4 and all(
        row.get("status") == "completed" and row.get("Vm_L2")
        for row in clean_window
    ))
    accepted = (len(clean_window) == 4
                and has_clean_errors
                and all(float(clean_window[i]["Vm_L2"]) <
                        float(clean_window[i - 1]["Vm_L2"])
                        for i in range(1, len(clean_window)))
                and len(orders) == 3
                and all(1.8 <= order <= 2.3 for order in orders))
    (root / "acceptance.txt").write_text(
        "criterion: Vm L2 decreases monotonically over dt=0.025, 0.0125, "
        "0.00625, 0.003125 and each halving order is within [1.8, 2.3].\n"
        f"observed_orders={','.join(f'{order:.6g}' for order in orders)}\n"
        f"result={'PASS' if accepted else 'NOT PASS'}\n"
        "The dt=0.0015625 point is reported but excluded from the order gate; "
        "the backup reference shows that level approaching the fixed-mesh "
        "spatial floor.\n"
    )
    return 0 if all(row.get("status") == "completed" for row in rows) and accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
