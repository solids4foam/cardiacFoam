#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run_step(cmd):
    print("\n>>", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run Purkinje slab -> diffusivity -> scar -> conversion pipeline."
    )
    parser.add_argument(
        "--input",
        default="ASCIIlegacy.vtk",
        help="Input VTK for Purkinje slab tagging.",
    )
    parser.add_argument(
        "--purkinje-output",
        default="output/purkinjeLayer.vtk",
        help="Purkinje slab output VTK.",
    )
    parser.add_argument(
        "--diffusivity-output",
        default="output/Diffusion_purkinjeLayer.vtk",
        help="Diffusivity output VTK.",
    )
    parser.add_argument(
        "--scar-full-mesh",
        default="output/purkinjeLayer_Diffusivity_IDsGlobal.vtk",
        help="Full mesh VTK for scar tagging (must contain GlobalCellIds).",
    )
    parser.add_argument(
        "--scar-selection",
        default="output/scar_tissue_region.vtu",
        help="Selection VTU file for scar tagging.",
    )
    parser.add_argument(
        "--scar-output",
        default="output/purkinjeLayer_Diffusivity_scar.vtk",
        help="Scar-tagged output VTK.",
    )
    parser.add_argument(
        "--converted-output",
        default="output/purkinjeLayer_Diffusivity_scar_VtkUnstructured.vtk",
        help="Converted VTK output.",
    )
    parser.add_argument(
        "--skip-purkinje",
        action="store_true",
        help="Skip Purkinje slab tagging step.",
    )
    parser.add_argument(
        "--skip-diffusivity",
        action="store_true",
        help="Skip diffusivity tensor step.",
    )
    parser.add_argument(
        "--skip-scar",
        action="store_true",
        help="Skip scar tagging step.",
    )
    parser.add_argument(
        "--skip-convert",
        action="store_true",
        help="Skip conversion step.",
    )
    args = parser.parse_args()

    py = sys.executable

    if not args.skip_purkinje:
        run_step([
            py,
            str(ROOT / "purkinje_network" / "purkinje_slab" / "purkinje_vtk.py"),
            "--mode", "layer",
            "--input", args.input,
            "--output", args.purkinje_output,
        ])

    if not args.skip_diffusivity:
        run_step([
            py,
            str(ROOT / "diffusivity" / "diffusionTensor_vtk.py"),
            "--input", args.purkinje_output,
            "--output", args.diffusivity_output,
        ])

    if not args.skip_scar:
        run_step([
            py,
            str(ROOT / "scar_fibrosis" / "scar_vtk.py"),
            "--full-mesh", args.scar_full_mesh,
            "--selection", args.scar_selection,
            "--output", args.scar_output,
        ])

    if not args.skip_convert:
        run_step([
            py,
            str(ROOT / "fileConversion" / "ASCIIlegacyToVtkUnstructured.py"),
            "--input", args.scar_output,
            "--output", args.converted_output,
        ])


if __name__ == "__main__":
    main()
