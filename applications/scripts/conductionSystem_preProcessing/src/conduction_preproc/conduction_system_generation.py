#!/usr/bin/env python3

import argparse
import importlib.util
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = ROOT.parents[1]


def default_output_dir() -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return PROJECT_ROOT / "outputs" / f"run_{stamp}"


def load_config(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load config from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_step(cmd):
    print("\n>>", " ".join(cmd))
    subprocess.run(cmd, check=True)


def normalize_steps(steps: list[str]) -> list[str]:
    aliases = {
        "conversion": "convert",
        "convert": "convert",
        "file_conversion": "convert",
    }
    return [aliases.get(step, step) for step in steps]


def print_requirements(steps: list[str], paths: dict) -> None:
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


def validate_step_order(steps: list[str]) -> None:
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Run conduction preprocessing pipeline.")
    parser.add_argument(
        "--input",
        default=str(PROJECT_ROOT / "inputs" / "meshes" / "ASCIIlegacy.vtk"),
        help="Input VTK for Purkinje slab tagging.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (defaults to outputs/run_YYYYMMDD_HHMMSS).",
    )
    parser.add_argument(
        "--purkinje-output",
        default=None,
        help="Purkinje slab output VTK.",
    )
    parser.add_argument(
        "--diffusivity-output",
        default=None,
        help="Diffusivity output VTK.",
    )
    parser.add_argument(
        "--scar-full-mesh",
        default=None,
        help="Full mesh VTK for scar tagging (must contain GlobalCellIds).",
    )
    parser.add_argument(
        "--scar-selection",
        default=None,
        help="Selection VTU file for scar tagging.",
    )
    parser.add_argument(
        "--scar-output",
        default=None,
        help="Scar-tagged output VTK.",
    )
    parser.add_argument(
        "--converted-output",
        default=None,
        help="Converted VTK output.",
    )
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
        default=str(
            PROJECT_ROOT
            / "inputs"
            / "meshes"
            / "biv_ellipsoid.msh"
        ),
        help="Mesh for fractal-tree generation.",
    )
    parser.add_argument(
        "--fractal-config",
        default=str(
            PROJECT_ROOT
            / "configs"
            / "purkinjeFractalTree_config.py"
        ),
        help="Config for fractal-tree generation.",
    )
    args = parser.parse_args()

    py = sys.executable
    steps = normalize_steps([s.lower() for s in args.steps])
    if "purkinje_fractal" in steps and "purkinje_slab" in steps:
        steps = [s for s in steps if s != "purkinje_slab"]
    validate_step_order(steps)

    cfg_purkinje = (
        load_config(PROJECT_ROOT / "configs" / "purkinjeSlab_config.py")
        if "purkinje_slab" in steps
        else None
    )
    cfg_diff = (
        load_config(PROJECT_ROOT / "configs" / "Diffusivity_config.py")
        if "diffusivity" in steps
        else None
    )
    cfg_scar = (
        load_config(PROJECT_ROOT / "configs" / "scar_config.py")
        if "scar" in steps
        else None
    )
    cfg_convert = (
        load_config(PROJECT_ROOT / "configs" / "convert_config.py")
        if "convert" in steps
        else None
    )
    cfg_fractal = (
        load_config(PROJECT_ROOT / "configs" / "purkinjeFractalTree_config.py")
        if "purkinje_fractal" in steps
        else None
    )

    requirements = {}
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
        requirements["purkinje_fractal"] = {
            "mesh": cfg_fractal.MESH,
            "config": str(PROJECT_ROOT / "configs" / "purkinjeFractalTree_config.py"),
        }
    print_requirements(steps, requirements)
    output_dir = Path(args.output_dir) if args.output_dir else default_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.purkinje_output is None and cfg_purkinje is not None:
        args.purkinje_output = cfg_purkinje.OUTPUT
    if args.diffusivity_output is None and cfg_diff is not None:
        args.diffusivity_output = cfg_diff.OUTPUT
    if args.scar_full_mesh is None:
        args.scar_full_mesh = None
    if args.scar_selection is None and cfg_scar is not None:
        args.scar_selection = cfg_scar.SELECTION
    if args.scar_output is None and cfg_scar is not None:
        args.scar_output = cfg_scar.OUTPUT
    if args.converted_output is None and cfg_convert is not None:
        args.converted_output = cfg_convert.OUTPUT
    current_mesh = args.input

    for step in steps:
        if step == "purkinje_slab":
            run_step([
                py,
                str(ROOT / "purkinje_network" / "purkinje_slab" / "purkinje_slab.py"),
                "--input", current_mesh,
                "--output", args.purkinje_output,
            ])
            current_mesh = args.purkinje_output
            continue

        if step == "purkinje_fractal":
            run_step([
                py,
                str(
                    PROJECT_ROOT
                    / "src"
                    / "conduction_preproc"
                    / "purkinje_network"
                    / "Costabal2015_purkinjeNetwork"
                    / "fractalTree_purkinjeNetwork.py"
                ),
                "--input", args.fractal_mesh,
                "--config", args.fractal_config,
            ])
            video = getattr(cfg_fractal, "VIDEO", {})
            if video.get("enabled"):
                pvpython = video.get("pvpython", "pvpython")
                script = video.get("script")
                if script and shutil.which(pvpython):
                    run_step([pvpython, script])
                else:
                    print(
                        f"Warning: pvpython or render script not found (pvpython={pvpython}, script={script})."
                    )
            continue

        if step == "diffusivity":
            run_step([
                py,
                str(ROOT / "diffusivity" / "diffusionTensor_vtk.py"),
                "--input", current_mesh,
                "--output", args.diffusivity_output,
            ])
            current_mesh = args.diffusivity_output
            continue

        if step == "scar":
            full_mesh = args.scar_full_mesh or current_mesh
            run_step([
                py,
                str(PROJECT_ROOT / "scripts" / "scar_vtk.py"),
                "--full-mesh", full_mesh,
                "--selection", args.scar_selection,
                "--output", args.scar_output,
            ])
            current_mesh = args.scar_output
            continue

        if step == "convert":
            run_step([
                py,
                str(ROOT / "fileConversion" / "ASCIIlegacyToVtkUnstructured.py"),
                "--input", current_mesh,
                "--output", args.converted_output,
            ])
            current_mesh = args.converted_output


if __name__ == "__main__":
    main()
