"""
diffusivity_cli.py

Standalone CLI to add a diffusivity tensor and optionally convert FIELD arrays.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pyvista as pv

from config_Diffusivity import (
    Ventricular_DF,
    Ventricular_DS,
    Ventricular_DN,
    PURKINJE_DIFFUSION_MULTIPLIER,
)
from lib.diffusivity_tensor import add_diffusivity_tensor_ventricles
from utils.vtk_convert_arrays_to_fields import convert_vtk_file
from utils.vtk_utils import remove_blank_lines, inspect_fields


def main() -> None:
    parser = argparse.ArgumentParser(description="Add diffusivity tensor to a VTK mesh.")
    parser.add_argument("--input", required=True, metavar="INPUT_VTK", help="Input VTK file.")
    parser.add_argument("--output", default=None, metavar="OUTPUT_VTK", help="Output VTK file.")
    parser.add_argument("--df", type=float, default=Ventricular_DF, metavar="DF_ventricle", help="Diffusivity along fiber.")
    parser.add_argument("--ds", type=float, default=Ventricular_DS, metavar="DS_ventricle", help="Diffusivity along sheet.")
    parser.add_argument("--dn", type=float, default=Ventricular_DN, metavar="DN_ventricle", help="Diffusivity along normal.")
    parser.add_argument("--purkinje-mult", type=float, default=PURKINJE_DIFFUSION_MULTIPLIER, metavar="Diffusion_multiplier_purkinje",
                        help="Diffusivity multiplier for purkinjeLayer cells.")
    parser.add_argument("--no-convert-fields", action="store_true", help="Skip converting FIELD arrays.")
    parser.add_argument("--no-remove-blank-lines", action="store_true", help="Skip removing blank lines.")
    parser.add_argument("--no-inspect", action="store_true", help="Skip inspecting fields before/after.")
    args = parser.parse_args()

    if not args.no_inspect:
        print("\n-------------------------\nInspecting Input Fields\n-------------------------")
        inspect_fields(args.input)

    mesh = pv.read(args.input)
    has_purkinje = "purkinjeLayer" in mesh.cell_data
    if args.output is None:
        args.output = (
            "output/Diffusion_purkinjeLayer.vtk" if has_purkinje else "output/Diffusion.vtk"
        )
    mesh = add_diffusivity_tensor_ventricles(
        mesh, args.df, args.ds, args.dn, purkinje_mult=args.purkinje_mult
    )
    mesh.save(args.output, binary=False)
    if has_purkinje:
        print(
            "Added diffusivity tensor with PURKINJE layer scaled by "
            f"{args.purkinje_mult}x ventricular diffusion."
        )
    else:
        print("Added diffusivity tensor (no purkinjeLayer present).")
    print(f"Mesh saved to {args.output} and is ready for OpenFOAM.")

    if not args.no_remove_blank_lines or not args.no_convert_fields:
        print("\n-------------------------\nPost-processing VTK\n-------------------------")

    if not args.no_remove_blank_lines:
        remove_blank_lines(args.output, args.output)

    if not args.no_convert_fields:
        convert_vtk_file(args.output, args.output)
        print("Converted FIELD arrays to SCALARS/VECTORS/TENSORS.")

    if not args.no_inspect:
        print("\n-------------------------\nInspecting Output Fields\n-------------------------")
        inspect_fields(args.output)


if __name__ == "__main__":
    main()
