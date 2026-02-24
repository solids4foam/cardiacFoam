"""CLI orchestration for conduction preprocessing pipelines."""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from types import ModuleType

from cardiac_engine.paths import resolve_paths
from cardiac_engine.runtime import CardiacPreprocEngine


def _load_config(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load config from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _run_cmd(cmd: list[str]) -> None:
    print("\n>>", " ".join(cmd))
    subprocess.run(cmd, check=True)


def _normalize_steps(steps: list[str]) -> list[str]:
    aliases = {
        "conversion": "convert",
        "convert": "convert",
        "file_conversion": "convert",
    }
    return [aliases.get(step, step) for step in steps]


def _validate_step_order(steps: list[str]) -> None:
    order = ["diffusivity", "purkinje_slab", "scar", "convert"]
    last_index = -1
    for step in order:
        if step in steps:
            idx = steps.index(step)
            if idx < last_index:
                raise RuntimeError(
                    "Invalid step order. Use: diffusivity -> purkinje_slab -> scar -> convert."
                )
            last_index = idx


def _print_requirements(steps: list[str], paths: dict) -> None:
    print("\nPlanned steps and requirements:")
    if "purkinje_slab" in steps:
        req_in = paths.get("purkinje_slab", {}).get("input", "VTK input")
        req_out = paths.get("purkinje_slab", {}).get("output", "VTK output")
        print("- purkinje_slab:")
        print(f"  input: {req_in} (needs POINT_DATA uvc_transmural + uvc_intraventricular)")
        print(f"  output: {req_out} (writes CELL_DATA purkinjeLayer)")
        print("  note: run after diffusivity to scale Diffusivity on purkinjeLayer cells.")
    if "diffusivity" in steps:
        req_in = paths.get("diffusivity", {}).get("input", "VTK input")
        req_out = paths.get("diffusivity", {}).get("output", "VTK output")
        print("- diffusivity:")
        print(f"  input: {req_in} (needs CELL_DATA fiber/sheet)")
        print(f"  output: {req_out} (writes CELL_DATA Diffusivity)")
    if "scar" in steps:
        full_mesh = paths.get("scar", {}).get("full_mesh", "full mesh VTK")
        selection = paths.get("scar", {}).get("selection", "selection VTU")
        req_out = paths.get("scar", {}).get("output", "VTK output")
        print("- scar:")
        print(f"  input: {full_mesh} (needs CELL_DATA GlobalCellIds)")
        print(f"  selection: {selection} (needs CELL_DATA GlobalCellIds)")
        print(f"  output: {req_out} (writes CELL_DATA Scar)")
    if "convert" in steps:
        req_in = paths.get("convert", {}).get("input", "VTK input")
        req_out = paths.get("convert", {}).get("output", "VTK output")
        print("- convert:")
        print(f"  input: {req_in} (VTK to convert/inspect)")
        print(f"  output: {req_out} (FIELD arrays converted)")
    if "purkinje_fractal" in steps:
        mesh = paths.get("purkinje_fractal", {}).get("mesh", "surface mesh")
        cfg = paths.get("purkinje_fractal", {}).get("config", "fractal config")
        video = paths.get("purkinje_fractal", {}).get("video", {})
        print("- purkinje_fractal:")
        print(f"  mesh: {mesh} (triangulated surface, ENDO_LV/ENDO_RV tags)")
        print(f"  config: {cfg}")
        if video.get("enabled"):
            print("  video: enabled (videoEditorPurkinje/render_purkinje_frames.py via pvpython)")


def _default_output_dir(outputs_root: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return outputs_root / f"run_{stamp}"


def build_parser(*, project_root: Path) -> argparse.ArgumentParser:
    paths = resolve_paths(project_root=project_root)
    parser = argparse.ArgumentParser(description="Run conduction preprocessing pipeline.")
    parser.add_argument(
        "--input",
        default=str(paths.default_input_mesh),
        help="Input VTK for Purkinje slab tagging.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (defaults to outputs/run_YYYYMMDD_HHMMSS).",
    )
    parser.add_argument("--purkinje-output", default=None, help="Purkinje slab output VTK.")
    parser.add_argument("--diffusivity-output", default=None, help="Diffusivity output VTK.")
    parser.add_argument(
        "--scar-full-mesh",
        default=None,
        help="Full mesh VTK for scar tagging (must contain GlobalCellIds).",
    )
    parser.add_argument("--scar-selection", default=None, help="Selection VTU file for scar tagging.")
    parser.add_argument("--scar-output", default=None, help="Scar-tagged output VTK.")
    parser.add_argument("--converted-output", default=None, help="Converted VTK output.")
    parser.add_argument(
        "--steps",
        nargs="+",
        required=True,
        help=(
            "Steps to run (e.g. diffusivity purkinje_slab scar convert). "
            "For conversion-only runs, use file_conversion."
        ),
    )
    parser.add_argument(
        "--fractal-mesh",
        default=str(paths.default_fractal_mesh),
        help="Mesh for fractal-tree generation.",
    )
    parser.add_argument(
        "--fractal-config",
        default=str(paths.purkinje_fractal_config),
        help="Config for fractal-tree generation.",
    )
    return parser


def run_pipeline_from_args(*, project_root: Path, args: argparse.Namespace) -> None:
    paths = resolve_paths(project_root=project_root)
    py = sys.executable

    steps = _normalize_steps([s.lower() for s in args.steps])
    if "purkinje_fractal" in steps and "purkinje_slab" in steps:
        steps = [s for s in steps if s != "purkinje_slab"]
    _validate_step_order(steps)

    cfg_purkinje = _load_config(paths.purkinje_slab_config) if "purkinje_slab" in steps else None
    cfg_diff = _load_config(paths.diffusivity_config) if "diffusivity" in steps else None
    cfg_scar = _load_config(paths.scar_config) if "scar" in steps else None
    cfg_convert = _load_config(paths.convert_config) if "convert" in steps else None
    cfg_fractal = _load_config(paths.purkinje_fractal_config) if "purkinje_fractal" in steps else None

    requirements: dict[str, dict[str, str]] = {}
    if cfg_purkinje is not None:
        requirements["purkinje_slab"] = {"input": cfg_purkinje.INPUT, "output": cfg_purkinje.OUTPUT}
    if cfg_diff is not None:
        requirements["diffusivity"] = {"input": cfg_diff.INPUT, "output": cfg_diff.OUTPUT}
    if cfg_scar is not None:
        requirements["scar"] = {
            "full_mesh": cfg_scar.FULL_MESH,
            "selection": cfg_scar.SELECTION,
            "output": cfg_scar.OUTPUT,
        }
    if cfg_convert is not None:
        requirements["convert"] = {"input": cfg_convert.INPUT, "output": cfg_convert.OUTPUT}
    if cfg_fractal is not None:
        requirements["purkinje_fractal"] = {"mesh": cfg_fractal.MESH, "config": str(paths.purkinje_fractal_config)}
    _print_requirements(steps, requirements)

    output_dir = Path(args.output_dir) if args.output_dir else _default_output_dir(paths.default_outputs_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    engine = CardiacPreprocEngine(project_root=project_root, output_dir=output_dir)

    if args.purkinje_output is None and cfg_purkinje is not None:
        args.purkinje_output = cfg_purkinje.OUTPUT
    if args.diffusivity_output is None and cfg_diff is not None:
        args.diffusivity_output = cfg_diff.OUTPUT
    if args.scar_selection is None and cfg_scar is not None:
        args.scar_selection = cfg_scar.SELECTION
    if args.scar_output is None and cfg_scar is not None:
        args.scar_output = cfg_scar.OUTPUT
    if args.converted_output is None and cfg_convert is not None:
        args.converted_output = cfg_convert.OUTPUT

    current_mesh = args.input
    for step in steps:
        if step == "purkinje_slab":
            result = engine.run_step(
                "purkinje_slab",
                current_mesh=current_mesh,
                options={
                    "input_path": current_mesh,
                    "output_path": args.purkinje_output,
                    "transmural_min": float(getattr(cfg_purkinje, "PURKINJE_TRANSMURAL_MIN", 0.0)),
                    "transmural_max": float(getattr(cfg_purkinje, "PURKINJE_TRANSMURAL_MAX", 0.1)),
                    "lv_value": int(getattr(cfg_purkinje, "LV_INTRAVENTRICULAR_TAG", -1)),
                    "rv_value": int(getattr(cfg_purkinje, "RV_INTRAVENTRICULAR_TAG", 1)),
                    "field_name": str(getattr(cfg_purkinje, "PURKINJE_FIELD_NAME", "purkinjeLayer")),
                    "purkinje_mult": getattr(cfg_purkinje, "PURKINJE_DIFFUSION_MULTIPLIER", None),
                    "wall_tag_name": getattr(cfg_purkinje, "WALL_TAG_NAME", None),
                    "wall_tag_value": getattr(cfg_purkinje, "WALL_TAG_VALUE", None),
                    "inflate_seed_point": getattr(cfg_purkinje, "INFLATE_SEED_POINT", None),
                    "inflate_inside_tag_name": str(
                        getattr(cfg_purkinje, "INFLATE_INSIDE_TAG_NAME", "inside_shared_boundary")
                    ),
                    "inspect": False,
                    "convert_fields": True,
                    "remove_blank_lines": True,
                },
            )
            current_mesh = result.output_mesh or args.purkinje_output
            continue

        if step == "purkinje_fractal":
            _run_cmd([py, str(paths.purkinje_fractal_runner), "--input", args.fractal_mesh, "--config", args.fractal_config])
            video = getattr(cfg_fractal, "VIDEO", {})
            if video.get("enabled"):
                pvpython = video.get("pvpython", "pvpython")
                script = video.get("script")
                if script and shutil.which(pvpython):
                    _run_cmd([pvpython, script])
                else:
                    print(f"Warning: pvpython or render script not found (pvpython={pvpython}, script={script}).")
            continue

        if step == "diffusivity":
            result = engine.run_step(
                "diffusivity",
                current_mesh=current_mesh,
                options={
                    "input_path": current_mesh,
                    "output_path": args.diffusivity_output,
                    "df": cfg_diff.Ventricular_DF,
                    "ds": cfg_diff.Ventricular_DS,
                    "dn": cfg_diff.Ventricular_DN,
                    "purkinje_mult": None,
                    "inspect": False,
                    "convert_fields": True,
                    "remove_blank_lines": True,
                },
            )
            current_mesh = result.output_mesh or args.diffusivity_output
            continue

        if step == "scar":
            full_mesh = args.scar_full_mesh or current_mesh
            diffusivity_scale = 1.0
            if getattr(cfg_scar, "DIFFUSIVITY_MODE", "").lower() == "constant":
                diffusivity_scale = float(getattr(cfg_scar, "DIFFUSIVITY_SCAR_MULTIPLIER", 0.0))
            result = engine.run_step(
                "scar",
                current_mesh=full_mesh,
                options={
                    "full_mesh_path": full_mesh,
                    "selection_path": args.scar_selection,
                    "output_path": args.scar_output,
                    "id_array": "GlobalCellIds",
                    "scar_name": "Scar",
                    "scar_value": 1.0,
                    "diffusivity_name": "Diffusivity",
                    "diffusivity_scale": diffusivity_scale,
                    "inspect": False,
                    "convert_fields": True,
                    "remove_blank_lines": True,
                },
            )
            current_mesh = result.output_mesh or args.scar_output
            continue

        if step == "convert":
            _run_cmd([py, str(paths.convert_runner), "--input", current_mesh, "--output", args.converted_output])
            current_mesh = args.converted_output


def run_pipeline_cli(*, project_root: Path, argv: list[str] | None = None) -> None:
    parser = build_parser(project_root=project_root)
    args = parser.parse_args(argv)
    run_pipeline_from_args(project_root=project_root, args=args)
