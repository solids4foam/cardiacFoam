"""
Master pipeline for Niederer benchmark simulation + postprocessing.
Includes a complete CLI interface.

Usage examples:

    python main.py all
    python main.py sim
    python main.py post


    python main.py all --case runs/Sim01
    python main.py post --case postResults
"""
"""
Main pipeline controller for NiedererEtAl2012 benchmark.

This script:
 - Sets up simulation
 - Runs OpenFOAM cases
 - Exports Paraview data
 - Runs line + point post-processing
 - Lets the user control everything from CLI

Configuration is defined here.

================================================================================
Main Setup Pipeline for Niederer Et Al. 2012 Slab Benchmark
================================================================================

This file holds ALL CONFIGURATION for:
- Simulation parameters (dx, dt, ionic models, tissue types)
- Folder paths
- Paraview export script
- run_cases.sh
- Benchmark reference paths

All other scripts import configuration ONLY from here.

================================================================================
"""

from pathlib import Path
import argparse
import subprocess
import importlib.util

# Determine important directories automatically
SETUP_ROOT = Path(__file__).resolve().parent
CASE_ROOT = SETUP_ROOT.parent

CONFIG = {
    # Folders
    "setup_root": SETUP_ROOT,
    "case_root": CASE_ROOT,

    # Simulation tools
    "pvpython_path": "/Applications/ParaView-5.13.3.app/Contents/bin/pvpython",
    "paraview_script_path": SETUP_ROOT / "simulation" / "exportParaview_Niederer.py",
    "run_cases_script": SETUP_ROOT / "simulation" / "run_cases.sh",
    # Input files
    "points_file": SETUP_ROOT / "simulation" / "Niederer_Points.txt",
    # Benchmark reference
    "excel_reference": SETUP_ROOT / "postProcessing" / 
                       "Niederer_graphs_webplotdigitilizer_points_slab" /
                       "WebPlotDigitilizerdata.xlsx",
}

# ===============================================================
# PARAMETER SWEEP CONFIGURATION
# ===============================================================
DX_VALUES = [0.5, 0.2, 0.1]                     # mm ------ This is a map for a structured blockMesh following rectangle dimensions.git   
DX_VALUES = [0.5, 0.2]                                   
DT_VALUES = [0.05, 0.01, 0.005]                 # ms
DT_VALUES = [0.05] 
TISSUE_TYPES = ["epicardialCells"]              
IONIC_MODEL = ["TNNP"] 
SOLVER = ["explicit", "implicit"]
OUTPUT_FOLDER = CASE_ROOT / "testSolvers"
                        


# ===============================================================
# Utility: dynamic module loader
# ===============================================================
def load_module(path):
    spec = importlib.util.spec_from_file_location("mod", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ===============================================================
# SIMULATION STAGES
# ===============================================================

def stage_simulations():
    """Run all Δx, Δt, ionic, tissue simulations by calling the sweep script."""
    print("\n--- Running simulation parameter sweep ---\n")

    sweep_script = Path(CONFIG["case_root"]) / "setupNiedererEtAl2012" / "simulation" / "setup_multiple_simulations_Niederer.py"

    subprocess.run(["python", str(sweep_script)], check=True)


def stage_postprocess():
    """Run line + point postprocessing."""

    print("\n--- Running Line Postprocessing ---\n")
    line_script = load_module(Path(CONFIG["case_root"])/ "setupNiedererEtAl2012" / "postProcessing" / "line_postProcessing.py")
    line_script.plot_line_csvs(
        folder=Path(CONFIG["case_root"]) / OUTPUT_FOLDER,
        excel_path=Path(CONFIG["case_root"]) / "setupNiedererEtAl2012" /CONFIG["excel_reference"]
    )

    print("\n--- Running Points Postprocessing ---\n")
    points_script = load_module(Path(CONFIG["case_root"])/ "setupNiedererEtAl2012" / "postProcessing" / "points_postProcessing.py")
    points_script.plot_3d_points_and_grid(
        folder=Path(CONFIG["case_root"]) / OUTPUT_FOLDER
    )


# ===============================================================
# CLI
# ===============================================================

def build_cli():
    parser = argparse.ArgumentParser(description="Niederer Et Al. 2012 Automation")

    parser.add_argument(
        "action",
        choices=["sim", "post", "all"],
        help="Which stage to execute"
    )

    parser.add_argument(
        "--case",
        help="Override case_root folder",
        default=None
    )

    return parser


def main():
    parser = build_cli()
    args = parser.parse_args()

    if args.case is not None:
        CONFIG["case_root"] = args.case

    if args.action == "sim":
        stage_simulations()

    elif args.action == "post":
        stage_postprocess()

    elif args.action == "all" or "":
        stage_simulations()
        stage_postprocess()

    print("\n🎉 Completed successfully!\n")


if __name__ == "__main__":
    main()
