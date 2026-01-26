"""
scar_vtk.py

Standalone CLI to tag scar cells from a selection mesh and optionally zero diffusivity.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pyvista as pv

from lib.scar import add_scar_from_selection


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
        default="output/purkinjeLayer_Diffusivity_scar.vtk",
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
    args = parser.parse_args()

    full = pv.read(args.full_mesh)
    sel = pv.read(args.selection)

    if args.print_arrays:
        print("Full mesh cell data arrays:", list(full.cell_data.keys()))
        print("Selection cell data arrays:", list(sel.cell_data.keys()))

    full, scar_mask = add_scar_from_selection(
        full,
        sel,
        id_array=args.id_array,
        scar_name=args.scar_name,
        diffusivity_name=args.diffusivity_name,
        scar_value=args.scar_value,
        diffusivity_scale=args.diffusivity_scale,
        return_mask=True,
    )

    print("Scar cells:", int(scar_mask.sum()), "out of", full.n_cells)

    full.save(args.output, binary=False)
    print(f"Scar-tagged mesh written to {args.output}")


if __name__ == "__main__":
    main()
