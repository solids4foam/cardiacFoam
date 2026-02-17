#!/usr/bin/env python3
"""Purkinje slab tagging (layer).

Inputs:
- Legacy ASCII VTK with POINT_DATA uvc_transmural + uvc_intraventricular
  and CELL_DATA Diffusivity (from the diffusivity step).

Outputs:
- VTK with CELL_DATA purkinjeLayer (default outputs/purkinjeLayer.vtk).
  If Diffusivity exists in the input, it is scaled in purkinjeLayer cells.

Main:
- Wraps purkinje_slab.py with --mode layer using config defaults.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(ROOT / "configs"))

from purkinjeSlab_config import INPUT, OUTPUT  # noqa: E402
from conduction_preproc.purkinje_network.purkinje_slab.purkinje_slab import (  # noqa: E402
    main as slab_main,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Purkinje slab tagging.")
    parser.add_argument("--input", default=INPUT, help="Input VTK.")
    parser.add_argument("--output", default=OUTPUT, help="Output VTK.")
    args, unknown = parser.parse_known_args()

    sys.argv = [sys.argv[0], "--input", args.input, "--output", args.output, *unknown]
    slab_main()


if __name__ == "__main__":
    main()
