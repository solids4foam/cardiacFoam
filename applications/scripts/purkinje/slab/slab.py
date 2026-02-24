#!/usr/bin/env python3
"""Purkinje slab product CLI."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "cardiac_preproc" / "src"
DEFAULT_CONFIG = Path(__file__).resolve().with_name("purkinjeSlab_config.py")

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from cardiac_preproc.steps.purkinje_slab import PurkinjeSlabOptions, run_purkinje_slab  # noqa: E402


def _load_config(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load config from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG),
        help="Path to Purkinje slab config file.",
    )
    pre_args, _ = pre_parser.parse_known_args()
    cfg = _load_config(Path(pre_args.config))

    parser = argparse.ArgumentParser(
        description="Purkinje slab tagging CLI.",
        parents=[pre_parser],
    )
    parser.add_argument("--input", default=cfg.INPUT, metavar="INPUT_VTK", help="Input VTK file.")
    parser.add_argument("--output", default=cfg.OUTPUT, metavar="OUTPUT_VTK", help="Output VTK file.")
    parser.add_argument("--transmural-min", type=float, default=None, metavar="MIN_DEPTH", help="Min transmural.")
    parser.add_argument("--transmural-max", type=float, default=None, metavar="MAX_DEPTH", help="Max transmural.")
    parser.add_argument("--lv-value", type=int, default=None, metavar="LV_TAG", help="LV tag in uvc_intraventricular.")
    parser.add_argument("--rv-value", type=int, default=None, metavar="RV_TAG", help="RV tag in uvc_intraventricular.")
    parser.add_argument("--field-name", default=None, metavar="FIELD_NAME", help="Cell data field name.")
    parser.add_argument("--wall-tag-name", default=cfg.WALL_TAG_NAME, metavar="WALL_TAG_NAME", help="Optional wall/boundary tag name.")
    parser.add_argument("--wall-tag-value", type=float, default=cfg.WALL_TAG_VALUE, metavar="WALL_TAG_VALUE", help="Optional numeric value for wall tag filtering.")
    parser.add_argument(
        "--inflate-seed-point",
        default=cfg.INFLATE_SEED_POINT,
        metavar="SEED_POINT",
        help="Seed point 'x,y,z' for inflate-from-point wall detection.",
    )
    parser.add_argument("--inflate-inside-tag-name", default=cfg.INFLATE_INSIDE_TAG_NAME, metavar="INSIDE_TAG", help="Internal tag name for the inflate-from-point region.")
    parser.add_argument("--purkinje-mult", type=float, default=cfg.PURKINJE_DIFFUSION_MULTIPLIER, metavar="PURKINJE_MULT", help="Multiplier for Diffusivity in purkinjeLayer cells.")
    parser.add_argument("--no-inspect", action="store_true", help="Skip inspecting fields before/after.")
    parser.add_argument("--no-convert-fields", action="store_true", help="Skip converting FIELD arrays.")
    parser.add_argument("--no-remove-blank-lines", action="store_true", help="Skip removing blank lines.")
    args = parser.parse_args()

    print("\n-------------------------\nAdding Purkinje Layer\n-------------------------")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    t_min = cfg.PURKINJE_TRANSMURAL_MIN if args.transmural_min is None else args.transmural_min
    t_max = cfg.PURKINJE_TRANSMURAL_MAX if args.transmural_max is None else args.transmural_max
    print(f"Transmural range: {t_min} to {t_max} ({(t_max - t_min) * 100.0:g}% depth)")
    if args.wall_tag_name:
        print(f"Wall restriction: tag='{args.wall_tag_name}' value={args.wall_tag_value}")
    if args.inflate_seed_point:
        print(f"Inflate-from-point seed: {args.inflate_seed_point}")

    run_purkinje_slab(
        PurkinjeSlabOptions(
            input_path=args.input,
            output_path=args.output,
            transmural_min=t_min,
            transmural_max=t_max,
            lv_value=cfg.LV_INTRAVENTRICULAR_TAG if args.lv_value is None else args.lv_value,
            rv_value=cfg.RV_INTRAVENTRICULAR_TAG if args.rv_value is None else args.rv_value,
            field_name=cfg.PURKINJE_FIELD_NAME if args.field_name is None else args.field_name,
            purkinje_mult=args.purkinje_mult,
            wall_tag_name=args.wall_tag_name,
            wall_tag_value=args.wall_tag_value,
            inflate_seed_point=args.inflate_seed_point,
            inflate_inside_tag_name=args.inflate_inside_tag_name,
            inspect=not args.no_inspect,
            convert_fields=not args.no_convert_fields,
            remove_blank_lines=not args.no_remove_blank_lines,
        )
    )


if __name__ == "__main__":
    main()
