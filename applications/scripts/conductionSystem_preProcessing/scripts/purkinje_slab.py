#!/usr/bin/env python3
"""Purkinje slab tagging (layer).

Inputs:
- Legacy ASCII VTK with POINT_DATA uvc_transmural + uvc_intraventricular
  and CELL_DATA Diffusivity (from the diffusivity step).

Outputs:
- VTK with CELL_DATA purkinjeLayer (default outputs/purkinjeLayer.vtk).
  If Diffusivity exists in the input, it is scaled in purkinjeLayer cells.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(ROOT / "configs"))

import pyvista as pv  # noqa: E402

from purkinjeSlab_config import (  # noqa: E402
    INPUT,
    OUTPUT,
    PURKINJE_TRANSMURAL_MIN,
    PURKINJE_TRANSMURAL_MAX,
    LV_INTRAVENTRICULAR_TAG,
    RV_INTRAVENTRICULAR_TAG,
    PURKINJE_FIELD_NAME,
    PURKINJE_DIFFUSION_MULTIPLIER,
    WALL_TAG_NAME,
    WALL_TAG_VALUE,
    INFLATE_SEED_POINT,
    INFLATE_INSIDE_TAG_NAME,
)
from conduction_preproc.purkinje_network.purkinje_fractal.slab import add_purkinje_layer  # noqa: E402
from conduction_preproc.utils.vtk_utils import inspect_fields, remove_blank_lines  # noqa: E402
from conduction_preproc.utils.vtk_convert_arrays_to_fields import convert_vtk_file  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Purkinje slab tagging CLI.")

    parser.add_argument("--input", default=INPUT, metavar="INPUT_VTK", help="Input VTK file.")
    parser.add_argument("--output", default=OUTPUT, metavar="OUTPUT_VTK", help="Output VTK file.")
    parser.add_argument("--transmural-min", type=float, default=None, metavar="MIN_DEPTH", help="Min transmural.")
    parser.add_argument("--transmural-max", type=float, default=None, metavar="MAX_DEPTH", help="Max transmural.")
    parser.add_argument("--lv-value", type=int, default=None, metavar="LV_TAG", help="LV tag in uvc_intraventricular.")
    parser.add_argument("--rv-value", type=int, default=None, metavar="RV_TAG", help="RV tag in uvc_intraventricular.")
    parser.add_argument("--field-name", default=None, metavar="FIELD_NAME", help="Cell data field name.")
    parser.add_argument("--wall-tag-name", default=WALL_TAG_NAME, metavar="WALL_TAG_NAME", help="Optional wall/boundary tag name.")
    parser.add_argument("--wall-tag-value", type=float, default=WALL_TAG_VALUE, metavar="WALL_TAG_VALUE", help="Optional numeric value for wall tag filtering.")
    parser.add_argument(
        "--inflate-seed-point",
        default=INFLATE_SEED_POINT,
        metavar="SEED_POINT",
        help="Seed point 'x,y,z' for inflate-from-point wall detection.",
    )
    parser.add_argument("--inflate-inside-tag-name", default=INFLATE_INSIDE_TAG_NAME, metavar="INSIDE_TAG", help="Internal tag name for the inflate-from-point region.")
    parser.add_argument("--purkinje-mult", type=float, default=PURKINJE_DIFFUSION_MULTIPLIER, metavar="PURKINJE_MULT", help="Multiplier for Diffusivity in purkinjeLayer cells.")
    parser.add_argument("--no-inspect", action="store_true", help="Skip inspecting fields before/after.")
    parser.add_argument("--no-convert-fields", action="store_true", help="Skip converting FIELD arrays.")
    parser.add_argument("--no-remove-blank-lines", action="store_true", help="Skip removing blank lines.")

    args = parser.parse_args()

    if not args.no_inspect:
        print("\n-------------------------\nInspecting Input Fields\n-------------------------")
        inspect_fields(args.input)

    print("\n-------------------------\nAdding Purkinje Layer\n-------------------------")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    t_min = PURKINJE_TRANSMURAL_MIN if args.transmural_min is None else args.transmural_min
    t_max = PURKINJE_TRANSMURAL_MAX if args.transmural_max is None else args.transmural_max
    print(f"Transmural range: {t_min} to {t_max} ({(t_max - t_min) * 100.0:g}% depth)")
    if args.wall_tag_name:
        print(f"Wall restriction: tag='{args.wall_tag_name}' value={args.wall_tag_value}")
    if args.inflate_seed_point:
        print(f"Inflate-from-point seed: {args.inflate_seed_point}")

    mesh = pv.read(args.input)
    mesh = add_purkinje_layer(
        mesh,
        transmural_min=t_min,
        transmural_max=t_max,
        lv_value=LV_INTRAVENTRICULAR_TAG if args.lv_value is None else args.lv_value,
        rv_value=RV_INTRAVENTRICULAR_TAG if args.rv_value is None else args.rv_value,
        field_name=PURKINJE_FIELD_NAME if args.field_name is None else args.field_name,
        purkinje_mult=args.purkinje_mult,
        wall_tag_name=args.wall_tag_name,
        wall_tag_value=args.wall_tag_value,
        inflate_seed_point=args.inflate_seed_point,
        inflate_inside_tag_name=args.inflate_inside_tag_name,
    )
    mesh.save(args.output, binary=False)
    print(f"Purkinje layer written to {args.output}")

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
