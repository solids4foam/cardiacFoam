#!/usr/bin/env python3
"""Scar tagging from selection mesh.

Inputs:
- Full mesh VTK with CELL_DATA GlobalCellIds.
- Selection VTU with matching GlobalCellIds that have the scar area.

Outputs:
- VTK with CELL_DATA Scar tag (Diffusivity scaled by DIFFUSIVITY_SCAR_MULTIPLIER).

Main:
- Wraps conduction_preproc.scar.scar_vtk with config defaults.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(ROOT / "configs"))

from scar_config import (  # noqa: E402
    FULL_MESH,
    SELECTION,
    OUTPUT,
    DIFFUSIVITY_MODE,
    DIFFUSIVITY_SCAR_MULTIPLIER,
)
from conduction_preproc.scar.scar_vtk import main as scar_main  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Run scar tagging.")
    parser.add_argument("--full-mesh", default=FULL_MESH, help="Full mesh VTK.")
    parser.add_argument("--selection", default=SELECTION, help="Selection VTU.")
    parser.add_argument("--output", default=OUTPUT, help="Output VTK.")
    args, unknown = parser.parse_known_args()

    sys.argv = [
        sys.argv[0],
        "--full-mesh",
        args.full_mesh,
        "--selection",
        args.selection,
        "--output",
        args.output,
        *unknown,
    ]
    if DIFFUSIVITY_MODE == "constant":
        sys.argv.extend(["--diffusivity-scale", str(DIFFUSIVITY_SCAR_MULTIPLIER)])
    scar_main()


if __name__ == "__main__":
    main()
