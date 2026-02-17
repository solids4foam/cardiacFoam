#!/usr/bin/env python3
"""Diffusivity tensor tagging.

Inputs:
- Legacy ASCII VTK with CELL_DATA fiber/sheet (optional purkinjeLayer).

Outputs:
- VTK with CELL_DATA Diffusivity (default outputs/Diffusion_purkinjeLayer.vtk).

Main:
- Wraps diffusionTensor_vtk.py using config defaults.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(ROOT / "configs"))

from Diffusivity_config import INPUT, OUTPUT  # noqa: E402
from conduction_preproc.diffusivity.diffusionTensor_vtk import main as diff_main  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Add diffusivity tensor to a VTK mesh.")
    parser.add_argument("--input", default=INPUT, help="Input VTK.")
    parser.add_argument("--output", default=OUTPUT, help="Output VTK.")
    args, unknown = parser.parse_known_args()

    argv = [sys.argv[0], "--input", args.input, "--output", args.output]
    argv.extend(unknown)
    sys.argv = argv
    diff_main()


if __name__ == "__main__":
    main()
