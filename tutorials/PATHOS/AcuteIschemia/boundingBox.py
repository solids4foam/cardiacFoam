import argparse
from pathlib import Path

import numpy as np
import pyvista as pv


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Read a VTK file and compute the bounding box of electrode_endo_rv points.",
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the VTK file.",
    )
    args = parser.parse_args()

    vtk_path = Path(args.input)
    mesh = pv.read(str(vtk_path))

    tag = "electrode_endo_rv"

    if tag not in mesh.point_data:
        available = list(mesh.point_data.keys())
        raise KeyError(f"Tag not found. Available point-data tags: {available}")

    data = mesh.point_data[tag]

    indices = np.where(data == 1)[0]

    print(f"Found {len(indices)} electrode points")

    if len(indices) == 0:
        return 0

    coords = mesh.points[indices]

    min_corner = coords.min(axis=0)
    max_corner = coords.max(axis=0)

    print("Bounding box min:", min_corner)
    print("Bounding box max:", max_corner)

    bbox_path = vtk_path.with_name(f"{vtk_path.stem}_bounding_box.txt")

    bbox_path.write_text(
        f"Bounding box min: {min_corner[0]} {min_corner[1]} {min_corner[2]}\n"
        f"Bounding box max: {max_corner[0]} {max_corner[1]} {max_corner[2]}\n"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())