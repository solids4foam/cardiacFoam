#!/usr/bin/env python3
"""Diffusivity tensor tagging.

Inputs:
- Legacy ASCII VTK with CELL_DATA fiber/sheet (optional purkinjeLayer).

Outputs:
- VTK with CELL_DATA Diffusivity (default outputs/Diffusion_purkinjeLayer.vtk).
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(ROOT / "configs"))

from Diffusivity_config import INPUT, OUTPUT, Ventricular_DF, Ventricular_DS, Ventricular_DN  # noqa: E402
from conduction_preproc.steps.diffusivity import DiffusivityOptions, run_diffusivity  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Add diffusivity tensor to a VTK mesh.")
    parser.add_argument("--input", default=INPUT, metavar="INPUT_VTK", help="Input VTK file.")
    parser.add_argument("--output", default=OUTPUT, metavar="OUTPUT_VTK", help="Output VTK file.")
    parser.add_argument("--df", type=float, default=Ventricular_DF, metavar="DF_ventricle", help="Diffusivity along fiber.")
    parser.add_argument("--ds", type=float, default=Ventricular_DS, metavar="DS_ventricle", help="Diffusivity along sheet.")
    parser.add_argument("--dn", type=float, default=Ventricular_DN, metavar="DN_ventricle", help="Diffusivity along normal.")
    parser.add_argument(
        "--purkinje-mult",
        type=float,
        default=None,
        metavar="Diffusion_multiplier_purkinje",
        help="Optional diffusivity multiplier for purkinjeLayer cells.",
    )
    parser.add_argument("--no-convert-fields", action="store_true", help="Skip converting FIELD arrays.")
    parser.add_argument("--no-remove-blank-lines", action="store_true", help="Skip removing blank lines.")
    parser.add_argument("--no-inspect", action="store_true", help="Skip inspecting fields before/after.")
    args = parser.parse_args()

    output = args.output
    if output == OUTPUT and output is None:
        import pyvista as pv
        mesh = pv.read(args.input)
        has_purkinje = "purkinjeLayer" in mesh.cell_data
        output = "outputs/Diffusion_purkinjeLayer.vtk" if has_purkinje else "outputs/Diffusion.vtk"

    result = run_diffusivity(
        DiffusivityOptions(
            input_path=args.input,
            output_path=output,
            df=args.df,
            ds=args.ds,
            dn=args.dn,
            purkinje_mult=args.purkinje_mult,
            inspect=not args.no_inspect,
            convert_fields=not args.no_convert_fields,
            remove_blank_lines=not args.no_remove_blank_lines,
        )
    )
    print(f"Mesh saved to {result.output_mesh} and is ready for OpenFOAM.")


if __name__ == "__main__":
    main()
