#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
SRC_DIR = SCRIPT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pipeline import (
    detect_start_stage,
    run_pipeline,
)
from sort_folder import sort_folder


UTILITY_NAME = "cellML2foam"


def foam_banner(args, start, stages):
    print(f"============================================================")
    print(f" cardiacFoam : {UTILITY_NAME} v2")
    print(f"============================================================")
    print(f"  Input     : {args.input}")
    print(f"  Model     : {args.model}")
    print(f"  Stages    : {' -> '.join(stages)}")
    print(f"------------------------------------------------------------")


def print_manual_steps_notice(model_id):
    print(f"""
============================================================
POST-GENERATION CHECKLIST FOR: {model_id}
============================================================

The model has been generated, but manual verification of 
physiological semantics is required:

1) Verify the ALGEBRAIC[Iion_cm] mapping:
   Ensure Iion_cm in {model_id}.H correctly represents 
   the sum of all ionic currents:
   
   ALGEBRAIC[Iion_cm] = ALGEBRAIC[sum_of_currents];

2) Clean up redundant symbols:
   - Remove unused Myokit-generated RATES[0] (Vm derivative)
   - Remove stimulus protocol symbols if desired (e.g. pace, 
     t_end, t_amplitude) from Names.H and .H files.

3) Configure Tissue Flags (if applicable):
   If you use local tissue-specific data, map them to:
   
   tissueFlag == 1 : endo
   tissueFlag == 2 : mid
   tissueFlag == 3 : epi

4) Implement Voltage-Dependent Parameters:
   If parameters depend on Vm, update the lookup in:
   {model_id}DependencyMap().

Model generation complete. Proceed to compile in cardiacFoam.
============================================================
""")


def parse_args():
    p = argparse.ArgumentParser(
        prog=UTILITY_NAME,
        description="Convert CellML models to OpenFOAM-ready ionic model code (Modernized v2)",
    )

    p.add_argument(
        "input",
        type=Path,
        help="Input file (.cellml, .mmt, .c)",
    )

    p.add_argument(
        "--from",
        dest="start",
        help="Pipeline stage to start from (auto-detected if omitted)",
    )

    p.add_argument(
        "--to",
        dest="to",
        default="openfoam",
        help="Pipeline stage to end at (default: openfoam)",
    )

    p.add_argument(
        "--model",
        help="Model name (used for output files)",
    )

    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("."),
        help="Output directory (default: current directory)",
    )

    p.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose output",
    )

    return p.parse_args()


def main():
    args = parse_args()

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}")
        sys.exit(1)

    try:
        start = args.start or detect_start_stage(args.input)

        # Execute pipeline
        result = run_pipeline(
            args.input,
            start=start,
            end=args.to,
            model=args.model,
            verbose=args.verbose,
        )

        if result is not None:
            # Reconstruct stages for the banner if it returned a dict
            stages_list = ["openfoam"] # Minimalist hint
            foam_banner(args, start, stages_list)
            
            # Only create ionic folder if we reached openfoam
            if args.to == "openfoam":
                if args.model is None:
                    raise ValueError("--model <ModelName_Year> is required")

                if "_" not in args.model:
                    raise ValueError("--model must be ModelName_Year")

                model, year = args.model.rsplit("_", 1)
                year = int(year)

                sort_result = sort_folder(
                    model=model,
                    year=year,
                    outdir=args.outdir,
                    verbose=args.verbose,
                )
                if args.verbose:
                    print(f"    Created: {len(sort_result['created'])} files")
                    print(f"    Moved:   {len(sort_result['moved'])} files")
                print_manual_steps_notice(args.model)

    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
