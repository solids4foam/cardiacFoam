#!/usr/bin/env python3
"""Create a fresh, matched 2D slab case from the supplied template."""

import argparse
import shutil
from pathlib import Path


DEFAULTS = {
    "AlievPanfilov": ("myocyte", 2e-6, 5),
    "BuenoOrovio": ("epicardialCells", 2e-5, 5),
    "Courtemanche": ("myocyte", 2e-6, 5),
    "Fabbri": ("myocyte", 2e-6, 5),
    "Gaur": ("myocyte", 2e-6, 5),
    "Grandi": ("myocyte", 2e-6, 5),
    "PerisYague": ("myocyte", 2e-6, 5),
    "Stewart": ("myocyte", 2e-6, 5),
    "TNNP": ("epicardialCells", 2e-6, 5),
    "TWorld": ("epicardialCells", 2e-6, 5),
    "ToRORd_dynCl": ("epicardialCells", 2e-6, 5),
    "Trovato": ("myocyte", 2e-6, 5),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", choices=DEFAULTS)
    parser.add_argument("backend", choices=("scalar", "batched"))
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--delta-t", type=float)
    parser.add_argument("--substeps", type=int)
    parser.add_argument(
        "--integrator", choices=("euler", "rushLarsen"),
        default="rushLarsen",
    )
    parser.add_argument(
        "--time-coupling-scheme", choices=("godunov", "sbdf2"),
        default="godunov",
    )
    parser.add_argument("--end-time", type=float, default=0.015)
    parser.add_argument("--write-precision", type=int)
    args = parser.parse_args()
    tissue, default_dt, default_substeps = DEFAULTS[args.model]
    delta_t = args.delta_t if args.delta_t is not None else default_dt
    substeps = args.substeps if args.substeps is not None else default_substeps
    if delta_t <= 0 or substeps < 1 or args.end_time <= 0:
        parser.error("delta-t, end-time, and substeps must be positive")
    if args.backend == "scalar" and args.substeps is not None:
        parser.error("substeps apply only to the batched backend")

    template = Path(__file__).resolve().parent / "slab2D"
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        parser.error(f"output directory already exists: {output_dir}")
    shutil.copytree(template, output_dir)
    model_name = args.model + ("compactBatched" if args.backend == "batched" else "")
    properties = output_dir / "constant/electroProperties"
    data = properties.read_text().replace("TNNPcompactBatched", model_name)
    data = data.replace("tissue    epicardialCells;", f"tissue    {tissue};")
    if args.backend == "batched":
        data = data.replace(
            "batchedIntegrator rushLarsen;",
            f"batchedIntegrator {args.integrator};",
        )
        data = data.replace("batchedSubsteps 5;", f"batchedSubsteps {substeps};")
    else:
        data = data.replace("batchedIntegrator rushLarsen;", "")
        data = data.replace("batchedSubsteps 5;", "")
    properties.write_text(data)
    if args.time_coupling_scheme != "godunov":
        properties.write_text(
            properties.read_text().replace(
                "solutionAlgorithm    implicit;",
                "timeCouplingScheme sbdf2;\n    solutionAlgorithm    implicit;",
            )
        )
    control = output_dir / "system/controlDict"
    control_data = (
        control.read_text()
        .replace("deltaT    2e-6;", f"deltaT    {delta_t};")
        .replace("endTime    0.015;", f"endTime    {args.end_time};")
    )
    if args.write_precision is not None:
        if args.write_precision < 1:
            parser.error("write-precision must be positive")
        import re
        control_data, count = re.subn(
            r"(?m)^(\s*writePrecision\s+)[^;]+;",
            rf"\g<1>{args.write_precision};",
            control_data,
            count=1,
        )
        if count != 1:
            parser.error("could not find writePrecision in controlDict")
    control.write_text(control_data)
    print(output_dir)


if __name__ == "__main__":
    main()
