#!/usr/bin/env python3
"""Scar product CLI."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT
DEFAULT_CONFIG = Path(__file__).resolve().with_name("config_scar.py")

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from cardiac_core.steps.scar import ScarOptions, run_scar  # noqa: E402


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
        help="Path to scar config file.",
    )
    pre_args, _ = pre_parser.parse_known_args()
    cfg = _load_config(Path(pre_args.config))
    default_diff_scale = (
        cfg.DIFFUSIVITY_SCAR_MULTIPLIER
        if cfg.DIFFUSIVITY_MODE == "constant"
        else 0.0
    )

    parser = argparse.ArgumentParser(
        description="Add scar tags from a selection mesh.",
        parents=[pre_parser],
    )
    parser.add_argument("--full-mesh", default=cfg.FULL_MESH, metavar="FULL_MESH_VTK", help="Full mesh VTK file.")
    parser.add_argument("--selection", default=cfg.SELECTION, metavar="SELECTION_VTU", help="Selection mesh VTU file.")
    parser.add_argument("--output", default=cfg.OUTPUT, metavar="OUTPUT_VTK", help="Output VTK file.")
    parser.add_argument("--id-array", default="GlobalCellIds", metavar="ID_ARRAY", help="Cell data array name for global cell IDs.")
    parser.add_argument("--scar-name", default="Scar", metavar="SCAR_NAME", help="Cell data array name for the scar tag.")
    parser.add_argument("--scar-value", type=float, default=1.0, metavar="SCAR_VALUE", help="Value assigned to scar cells.")
    parser.add_argument("--diffusivity-name", default="Diffusivity", metavar="DIFFUSIVITY_NAME", help="Cell data array name for diffusivity tensor.")
    parser.add_argument(
        "--diffusivity-scale",
        type=float,
        default=default_diff_scale,
        metavar="SCAR_DIFF_SCALE",
        help="Scale diffusivity in scar cells (0.0 = zero, 1.0 = unchanged).",
    )
    parser.add_argument("--print-arrays", action="store_true", help="Print available cell data arrays for both meshes.")
    parser.add_argument("--no-convert-fields", action="store_true", help="Skip converting FIELD arrays.")
    parser.add_argument("--no-remove-blank-lines", action="store_true", help="Skip removing blank lines.")
    parser.add_argument("--no-inspect", action="store_true", help="Skip inspecting fields before/after.")
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
