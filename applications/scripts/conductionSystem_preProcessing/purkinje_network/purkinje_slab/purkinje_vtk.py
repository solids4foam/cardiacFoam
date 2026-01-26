"""
purkinje_cli.py

Standalone CLI for Purkinje processing (layer or fractal).
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import pyvista as pv

from config_purkinjeSlab import (
    PURKINJE_TRANSMURAL_MIN,
    PURKINJE_TRANSMURAL_MAX,
    LV_INTRAVENTRICULAR_TAG ,
    RV_INTRAVENTRICULAR_TAG ,
    PURKINJE_FIELD_NAME,
)
from lib.purkinje_slab import add_purkinje_layer
from lib.fractal_purkinje import generate_fractal_purkinje



def main() -> None:
    parser = argparse.ArgumentParser(description="Purkinje processing CLI.")
    parser.add_argument("--mode", choices=["layer", "fractal"], required=True)

    # Layer options
    parser.add_argument("--input", required=True, metavar="INPUT_VTK", help="Input VTK file.")
    parser.add_argument("--output", default="output/purkinjeLayer.vtk", metavar="OUTPUT_VTK", help="Output VTK file.")
    # Layer parameters (optional overrides).
    parser.add_argument("--transmural-min", type=float, default=None, metavar="MIN_DEPTH", help="Min transmural.")
    parser.add_argument("--transmural-max", type=float, default=None, metavar="MAX_DEPTH", help="Max transmural.")
    parser.add_argument("--lv-value", type=int, default=None, metavar="LV_TAG", help="LV tag in uvc_intraventricular.")
    parser.add_argument("--rv-value", type=int, default=None, metavar="RV_TAG", help="RV tag in uvc_intraventricular.")
    parser.add_argument("--field-name", default=None, metavar="FIELD_NAME", help="Cell data field name.")

    args = parser.parse_args()

    if args.mode == "layer":
        mesh = pv.read(args.input)
        mesh = add_purkinje_layer(
            mesh,
            transmural_min=PURKINJE_TRANSMURAL_MIN if args.transmural_min is None else args.transmural_min,
            transmural_max=PURKINJE_TRANSMURAL_MAX if args.transmural_max is None else args.transmural_max,
            lv_value=LV_INTRAVENTRICULAR_TAG if args.lv_value is None else args.lv_value,
            rv_value=RV_INTRAVENTRICULAR_TAG if args.rv_value is None else args.rv_value,
            field_name=PURKINJE_FIELD_NAME if args.field_name is None else args.field_name,
        )
        mesh.save(args.output, binary=False)
        print(f"Purkinje layer written to {args.output}")

    else:
        cfg = config_from_module().fractal
        generate_fractal_purkinje(cfg)


if __name__ == "__main__":
    main()
