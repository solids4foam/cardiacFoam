"""
Master pipeline for the singleCell workflow.

Runs:
  - setup_multiple_simulations_singleCell.py
  - singleCellinteractivePlots.py

Usage:
    python main_setup_singleCell.py all
    python main_setup_singleCell.py sim
    python main_setup_singleCell.py post
"""

import subprocess
import sys
from pathlib import Path




# Ionic models and their tissue types
IONIC_MODELS = ["TNNP", "Gaur", "Courtemanche", "BuenoOrovio"]


# PATHS (relative to this script's folder)
HERE = Path(__file__).resolve().parent
SCRIPT_SETUP = HERE / "setup_multiple_simulations_singleCell.py"
SCRIPT_POST = HERE / "singleCellinteractivePlots.py"

PYTHON = "python"




# ---------------FUNCTIONS----------------------
def run_simulation():
    """Run the full OpenFOAM + setup pipeline."""
    print("\n🚀 Running single-cell simulation setup\n")

    cmd = [PYTHON, str(SCRIPT_SETUP)]

    subprocess.run(cmd, check=True)
    print("✔ Simulation setup completed.\n")


def run_postprocessing():
    """Run post-processing and plotting."""
    print("\n📊 Running single-cell post-processing\n")

    cmd = [PYTHON, str(SCRIPT_POST)]

    subprocess.run(cmd, check=True)
    print("✔ Post-processing completed.\n")


def run_all():
    """Clean → simulate → post-process."""
    print("\n=======================================")
    print(" 🚀 SINGLE-CELL FULL PIPELINE")
    print("=======================================\n")

    run_simulation()
    run_postprocessing()

    print("🎉 All tasks completed successfully!\n")




# ----------------MAIN----------------------
def main():
    if len(sys.argv) < 2:
        print("❌ ERROR: Please specify a command: all | sim | post\n")
        sys.exit(1)

    command = sys.argv[1].lower()

    if command == "all":
        run_all()
    elif command == "sim":
        run_simulation()
    elif command == "post":
        run_postprocessing()
    else:
        print(f"❌ Unknown command: {command}")
        print("Valid commands: all | sim | post\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
