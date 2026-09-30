#!/usr/bin/env python3
"""Create a model-specific scalar or batched case on the Niederer slab."""

import argparse
import shutil
from pathlib import Path


TISSUE = {
    "AlievPanfilov": "myocyte",
    "BuenoOrovio": "epicardialCells",
    "Courtemanche": "myocyte",
    "Gaur": "myocyte",
    "Grandi": "myocyte",
    "PerisYague": "myocyte",
    "Stewart": "myocyte",
    "TNNP": "epicardialCells",
    "TWorld": "epicardialCells",
    "ToRORd_dynCl": "epicardialCells",
    "Trovato": "myocyte",
}


def replace_entry(text, name, value):
    import re

    pattern = rf"(?m)^(\s*{re.escape(name)}\s+).+?;\s*$"
    updated, count = re.subn(pattern, rf"\g<1>{value};", text, count=1)
    if count != 1:
        raise ValueError(f"could not find dictionary entry {name}")
    return updated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", choices=tuple(TISSUE))
    parser.add_argument("backend", choices=("scalar", "batched", "gpu"))
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--substeps", type=int, default=25)
    parser.add_argument("--delta-t", type=float)
    parser.add_argument("--end-time", type=float)
    parser.add_argument(
        "--time-coupling-scheme", choices=("godunov", "sbdf2"),
        default="godunov",
    )
    args = parser.parse_args()
    if args.substeps < 1 or (args.delta_t is not None and args.delta_t <= 0) \
            or (args.end_time is not None and args.end_time <= 0):
        parser.error("substeps, delta-t, and end-time must be positive")
    if args.backend == "scalar" and args.substeps != 25:
        parser.error("substeps apply only to batched cases")

    source = Path(__file__).resolve().parents[1] / "NiedererEtAl2011" / "NiedererEtAl2011verification"
    target = args.output_dir.resolve()
    if target.exists():
        parser.error(f"output directory already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, target)

    properties = target / "constant/electroProperties"
    data = properties.read_text()
    batched = args.backend in ("batched", "gpu")
    model_name = args.model + ("compactBatched" if batched else "")
    data = replace_entry(data, "ionicModel", model_name)
    data = replace_entry(data, "tissue", TISSUE[args.model])
    if args.time_coupling_scheme == "sbdf2":
        if "timeCouplingScheme" in data:
            data = replace_entry(data, "timeCouplingScheme", "sbdf2")
        else:
            data = data.replace(
                "    solutionAlgorithm    implicit;",
                "    timeCouplingScheme sbdf2;\n"
                "    solutionAlgorithm    implicit;",
                1,
            )
    data = replace_entry(data, "export", "(Vm)")
    if batched:
        if "batchedIntegrator" in data:
            data = replace_entry(data, "batchedIntegrator", "rushLarsen")
        else:
            data = data.replace(
                "    solver RKF45;",
                "    solver RKF45;\n"
                "    batchedIntegrator rushLarsen;",
                1,
            )
        if "batchedSubsteps" in data:
            data = replace_entry(data, "batchedSubsteps", str(args.substeps))
        else:
            data = data.replace(
                "    batchedIntegrator rushLarsen;",
                "    batchedIntegrator rushLarsen;\n"
                f"    batchedSubsteps {args.substeps};",
                1,
            )
    else:
        import re

        data = re.sub(r"(?m)^\s*batchedIntegrator\s+.*?;\s*\n", "", data)
        data = re.sub(r"(?m)^\s*batchedSubsteps\s+.*?;\s*\n", "", data)
    properties.write_text(data)

    control = target / "system/controlDict"
    control_data = control.read_text()
    if args.delta_t is not None:
        control_data = replace_entry(control_data, "deltaT", str(args.delta_t))
    if args.end_time is not None:
        control_data = replace_entry(control_data, "endTime", str(args.end_time))
    control.write_text(control_data)
    print(target)


if __name__ == "__main__":
    main()
