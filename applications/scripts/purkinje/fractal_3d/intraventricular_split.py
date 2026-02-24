"""Split a VTK mesh into uvc_intraventricular == 1 and == -1 outputs."""

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Write two VTKs split by uvc_intraventricular value (1 and -1)."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input legacy ASCII VTK.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional output directory (defaults to input file folder).",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    mesh = pv.read(str(input_path))

    field_name = "uvc_intraventricular"
    if field_name not in mesh.point_data:
        raise KeyError(f"Missing point data field '{field_name}' in {input_path}")

    values = np.asarray(mesh.point_data[field_name]).reshape(-1)
    out_dir = Path(args.output_dir) if args.output_dir else input_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    for label, target in (("plus1", 1.0), ("minus1", -1.0)):
        mask = np.isclose(values, target, atol=1e-6)
        count = int(np.count_nonzero(mask))
        print(f"Points with {field_name} == {target}: {count} / {values.size}")
        if count == 0:
            print(f"Skipping {label} output (no matching points).")
            continue
        filtered = mesh.extract_points(mask, adjacent_cells=False)
        output_path = out_dir / f"{input_path.stem}_{label}.vtk"
        filtered.save(str(output_path), binary=False)
        print(f"Wrote {label} mesh to: {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
