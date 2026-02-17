#!/usr/bin/env python3
"""Purkinje fractal tree generation (Costabal 2015).

Inputs:
- Triangulated surface mesh with ENDO_LV/ENDO_RV tags (ASCII legacy VTK or gmsh).
- Surface tag mapping in configs/purkinjeFractalTree_config.py.

Outputs:
- Fractal tree node/line files and geometry VTK under outputs/.

Main:
- Runs Costabal2015 fractalTree_purkinjeNetwork.py with config defaults.
  The fractal script will open an interactive viewer unless --no-view-final is set.
"""
import argparse
import sys
from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs"

sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(CONFIGS))

from purkinjeFractalTree_config import MESH  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Purkinje fractal tree generator.")
    parser.add_argument("--input", default=MESH, help="Input mesh.")
    parser.add_argument(
        "--config",
        default=str(ROOT / "configs" / "purkinjeFractalTree_config.py"),
        help="Config file.",
    )
    args, unknown = parser.parse_known_args()

    sys.argv = [
        sys.argv[0],
        "--input",
        args.input,
        "--config",
        args.config,
        *unknown,
    ]
    runpy.run_path(
        str(
            ROOT
            / "src"
            / "conduction_preproc"
            / "purkinje_network"
            / "Costabal2015_purkinjeNetwork"
            / "fractalTree_purkinjeNetwork.py"
        ),
        run_name="__main__",
    )


if __name__ == "__main__":
    main()
