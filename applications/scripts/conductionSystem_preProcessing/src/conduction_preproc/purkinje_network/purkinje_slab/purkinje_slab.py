"""
purkinje_slab.py

Standalone CLI for Purkinje slab tagging.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "configs"))
import pyvista as pv

from purkinjeSlab_config import (  # noqa: E402
    PURKINJE_TRANSMURAL_MIN,
    PURKINJE_TRANSMURAL_MAX,
    LV_INTRAVENTRICULAR_TAG ,
    RV_INTRAVENTRICULAR_TAG ,
    PURKINJE_FIELD_NAME,
    PURKINJE_DIFFUSION_MULTIPLIER,
)
from conduction_preproc.lib.purkinje_slab import add_purkinje_layer  # noqa: E402
from conduction_preproc.utils.vtk_utils import inspect_fields, remove_blank_lines  # noqa: E402
from conduction_preproc.utils.vtk_convert_arrays_to_fields import convert_vtk_file  # noqa: E402



def main() -> None:
    parser = argparse.ArgumentParser(description="Purkinje slab tagging CLI.")

    # Layer options
    parser.add_argument("--input", required=True, metavar="INPUT_VTK", help="Input VTK file.")
    parser.add_argument("--output", default="outputs/purkinjeLayer.vtk", metavar="OUTPUT_VTK", help="Output VTK file.")
    # Layer parameters (optional overrides).
    parser.add_argument("--transmural-min", type=float, default=None, metavar="MIN_DEPTH", help="Min transmural.")
    parser.add_argument("--transmural-max", type=float, default=None, metavar="MAX_DEPTH", help="Max transmural.")
    parser.add_argument("--lv-value", type=int, default=None, metavar="LV_TAG", help="LV tag in uvc_intraventricular.")
    parser.add_argument("--rv-value", type=int, default=None, metavar="RV_TAG", help="RV tag in uvc_intraventricular.")
    parser.add_argument("--field-name", default=None, metavar="FIELD_NAME", help="Cell data field name.")
    parser.add_argument(
        "--purkinje-mult",
        type=float,
        default=PURKINJE_DIFFUSION_MULTIPLIER,
        metavar="PURKINJE_MULT",
        help="Multiplier for Diffusivity in purkinjeLayer cells (if present).",
    )
    parser.add_argument(
        "--no-inspect",
        action="store_true",
        help="Skip inspecting fields before/after.",
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

    args = parser.parse_args()

    if not args.no_inspect:
        print("\n-------------------------\nInspecting Input Fields\n-------------------------")
        inspect_fields(args.input)

    print("\n-------------------------\nAdding Purkinje Layer\n-------------------------")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    t_min = PURKINJE_TRANSMURAL_MIN if args.transmural_min is None else args.transmural_min
    t_max = PURKINJE_TRANSMURAL_MAX if args.transmural_max is None else args.transmural_max
    print(
        "Transmural range: "
        f"{t_min} to {t_max} ({(t_max - t_min) * 100.0:g}% depth)"
    )

    mesh = pv.read(args.input)
    mesh = add_purkinje_layer(
        mesh,
        transmural_min=PURKINJE_TRANSMURAL_MIN if args.transmural_min is None else args.transmural_min,
        transmural_max=PURKINJE_TRANSMURAL_MAX if args.transmural_max is None else args.transmural_max,
        lv_value=LV_INTRAVENTRICULAR_TAG if args.lv_value is None else args.lv_value,
        rv_value=RV_INTRAVENTRICULAR_TAG if args.rv_value is None else args.rv_value,
        field_name=PURKINJE_FIELD_NAME if args.field_name is None else args.field_name,
        purkinje_mult=args.purkinje_mult,
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
