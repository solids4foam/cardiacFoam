#!/usr/bin/env python3
"""Diffusivity product CLI."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT
DEFAULT_CONFIG = Path(__file__).resolve().with_name("config_diffusivity.py")

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from cardiac_core.steps.diffusivity import DiffusivityOptions, run_diffusivity  # noqa: E402


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
        help="Path to Diffusivity config file.",
    )
    pre_args, _ = pre_parser.parse_known_args()
    cfg = _load_config(Path(pre_args.config))

    parser = argparse.ArgumentParser(
        description="Add diffusivity tensor to a VTK mesh.",
        parents=[pre_parser],
    )
    parser.add_argument("--input", default=cfg.INPUT, metavar="INPUT_VTK", help="Input VTK file.")
    parser.add_argument("--output", default=cfg.OUTPUT, metavar="OUTPUT_VTK", help="Output VTK file.")
    parser.add_argument(
        "--df",
        type=float,
        default=getattr(cfg, "Ventricular_DF", 0.1143),
        metavar="DF_ventricle",
        help="Diffusivity along fiber.",
    )
    parser.add_argument(
        "--ds",
        type=float,
        default=getattr(cfg, "Ventricular_DS", 0.052),
        metavar="DS_ventricle",
        help="Diffusivity along sheet.",
    )
    parser.add_argument(
        "--dn",
        type=float,
        default=getattr(cfg, "Ventricular_DN", 0.016),
        metavar="DN_ventricle",
        help="Diffusivity along normal.",
    )
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

    output = args.output
    if output == cfg.OUTPUT and output is None:
        import pyvista as pv

        mesh = pv.read(args.input)
        has_purkinje = "purkinjeLayer" in mesh.cell_data
        output = "outputs/Diffusion_purkinjeLayer.vtk" if has_purkinje else "outputs/Diffusion.vtk"

    result = run_diffusivity(
        DiffusivityOptions(
            input_path=args.input,
            output_path=output,
            df=args.df,
            ds=args.ds,
            dn=args.dn,
            purkinje_mult=args.purkinje_mult,
            inspect=not args.no_inspect,
            convert_fields=not args.no_convert_fields,
            remove_blank_lines=not args.no_remove_blank_lines,
        )
    )
    print(f"Mesh saved to {result.output_mesh} and is ready for OpenFOAM.")


if __name__ == "__main__":
    main()
