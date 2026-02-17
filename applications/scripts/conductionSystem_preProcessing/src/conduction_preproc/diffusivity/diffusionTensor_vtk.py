"""
diffusivity_cli.py

Standalone CLI to add a diffusivity tensor and optionally convert FIELD arrays.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "configs"))
import pyvista as pv

from Diffusivity_config import (  # noqa: E402
    Ventricular_DF,
    Ventricular_DS,
    Ventricular_DN,
)
from conduction_preproc.lib.diffusivity_tensor import (  # noqa: E402
    add_diffusivity_tensor_ventricles,
)
from conduction_preproc.utils.vtk_convert_arrays_to_fields import (  # noqa: E402
    convert_vtk_file,
)
from conduction_preproc.utils.vtk_utils import (  # noqa: E402
    remove_blank_lines,
    inspect_fields,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Add diffusivity tensor to a VTK mesh.")
    parser.add_argument("--input", required=True, metavar="INPUT_VTK", help="Input VTK file.")
    parser.add_argument("--output", default=None, metavar="OUTPUT_VTK", help="Output VTK file.")
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

    if not args.no_inspect:
        print("\n-------------------------\nInspecting Input Fields\n-------------------------")
        inspect_fields(args.input)

    mesh = pv.read(args.input)
    has_purkinje = "purkinjeLayer" in mesh.cell_data
    if has_purkinje:
        print(
            "Warning: purkinjeLayer found in input. "
            "Preferred order is diffusivity -> purkinje_slab."
        )
    if args.output is None:
        args.output = (
            "outputs/Diffusion_purkinjeLayer.vtk"
            if has_purkinje
            else "outputs/Diffusion.vtk"
        )
    purkinje_mult = 1.0 if args.purkinje_mult is None else args.purkinje_mult
    mesh = add_diffusivity_tensor_ventricles(
        mesh, args.df, args.ds, args.dn, purkinje_mult=purkinje_mult
    )
    mesh.save(args.output, binary=False)
    if has_purkinje and args.purkinje_mult is not None:
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
