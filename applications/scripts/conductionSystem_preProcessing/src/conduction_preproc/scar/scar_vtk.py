"""
scar_vtk.py

Standalone CLI to tag scar cells from a selection mesh and optionally zero diffusivity.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from conduction_preproc.steps.scar import ScarOptions, run_scar


def main() -> None:
    parser = argparse.ArgumentParser(description="Add scar tags from a selection mesh.")
    parser.add_argument(
        "--full-mesh",
        required=True,
        metavar="FULL_MESH_VTK",
        help="Full mesh VTK file (legacy ASCII recommended).",
    )
    parser.add_argument(
        "--selection",
        required=True,
        metavar="SELECTION_VTU",
        help="Selection mesh VTU file (subset of cells).",
    )
    parser.add_argument(
        "--output",
        default="outputs/purkinjeLayer_Diffusivity_scar.vtk",
        metavar="OUTPUT_VTK",
        help="Output VTK file.",
    )
    parser.add_argument(
        "--id-array",
        default="GlobalCellIds",
        metavar="ID_ARRAY",
        help="Cell data array name for global cell IDs.",
    )
    parser.add_argument(
        "--scar-name",
        default="Scar",
        metavar="SCAR_NAME",
        help="Cell data array name for the scar tag.",
    )
    parser.add_argument(
        "--scar-value",
        type=float,
        default=1.0,
        metavar="SCAR_VALUE",
        help="Value assigned to scar cells (typically 0.0 to 1.0).",
    )
    parser.add_argument(
        "--diffusivity-name",
        default="Diffusivity",
        metavar="DIFFUSIVITY_NAME",
        help="Cell data array name for diffusivity tensor.",
    )
    parser.add_argument(
        "--diffusivity-scale",
        type=float,
        default=0.0,
        metavar="SCAR_DIFF_SCALE",
        help="Scale diffusivity in scar cells (0.0 = zero, 1.0 = unchanged).",
    )
    parser.add_argument(
        "--print-arrays",
        action="store_true",
        help="Print available cell data arrays for both meshes.",
    )
    parser.add_argument(
        "--no-convert-fields",
        action="store_true",
        help="Skip converting FIELD arrays.",
    )
    parser.add_argument(
        "--no-remove-blank-lines",
        action="store_true",
        help="Skip removing blank lines.",
    )
    parser.add_argument(
        "--no-inspect",
        action="store_true",
        help="Skip inspecting fields before/after.",
    )
    args = parser.parse_args()

    if args.print_arrays:
        import pyvista as pv

        full = pv.read(args.full_mesh)
        sel = pv.read(args.selection)
        print("Full mesh cell data arrays:", list(full.cell_data.keys()))
        print("Selection cell data arrays:", list(sel.cell_data.keys()))

    run_scar(
        ScarOptions(
            full_mesh_path=args.full_mesh,
            selection_path=args.selection,
            output_path=args.output,
            id_array=args.id_array,
            scar_name=args.scar_name,
            scar_value=args.scar_value,
            diffusivity_name=args.diffusivity_name,
            diffusivity_scale=args.diffusivity_scale,
            inspect=not args.no_inspect,
            convert_fields=not args.no_convert_fields,
            remove_blank_lines=not args.no_remove_blank_lines,
        )
    )


if __name__ == "__main__":
    main()
