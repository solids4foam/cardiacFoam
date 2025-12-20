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
import sys
from pathlib import Path

from setup_multiple_simulations_s1s2 import run_s1s2_simulation
from postProcessing_restCurves import restitution_curves





# ================= USER PARAMETERS =================
IONIC_MODELS   = ["BuenoOrovio", "TNNP", "Gaur", "Courtemanche"]
IONIC_MODELS   = ["BuenoOrovio"]  # For quick testing
TISSUE_MAP = {
    "TNNP":         ["mCells","endocardialCells","epicardialCells"],
    "BuenoOrovio":  ["mCells","endocardialCells","epicardialCells"],
    "Gaur":         ["myocyte"],
    "Courtemanche": ["myocyte"],
             }

STIM_AMPLITUDE_MAP = {
    "TNNP":  60.0,
    "BuenoOrovio": 0.4,
    "Gaur":   60.0,
    "Courtemanche": 70.0,
}

S1_INTERVAL  = 2000
S2_INTERVAL  = [1500, 1200, 1000, 800, 600, 400, 390, 380, 370, 360, 350, 340, 330, 320, 310, 300, 290, 280, 270, 260, 250]
OUTPUT_FOLDER= "BuenoOrovio_S1S2protocol"



# This variables plots all of the Voltage traces during post-processing. It is time consuming for large datasets.
SEE_PLOTS_DURING_POSTPROCESSING = True

SETUP_ROOT = Path(__file__).resolve().parent
CASE_ROOT = SETUP_ROOT.parent
CONFIG = { "run_cases_script": SETUP_ROOT / "run_cases.sh"}





# ================= PIPELINE =================
def run_simulation():
    print("\n🚀 Running S1–S2 simulations\n")

    for model in IONIC_MODELS:
        try:
            stim_amplitude = STIM_AMPLITUDE_MAP[model]
            mapped_tissue = TISSUE_MAP[model]
        except KeyError:
            raise ValueError(
                f"No stim_amplitude defined for ionic model '{model}'"
            )
        for tissue in mapped_tissue:
            run_s1s2_simulation(
                base_dir=CASE_ROOT,
                ionic_model=model,
                tissue_type=tissue,
                s1_interval=S1_INTERVAL,
                s2_intervals=S2_INTERVAL,
                stim_amplitude=stim_amplitude,
                run_cases_script=CONFIG["run_cases_script"],
                output_folder= OUTPUT_FOLDER
            )

    print("\n✔ Simulations completed\n")


def run_postprocessing():

    print("\n📊 Running S1–S2 post-processing\n")
    for model in IONIC_MODELS:
        mapped_tissue = TISSUE_MAP[model]
        print("\n➡ Analyzing restitution_curves with:")
        print(f"   ionic_model      = {model!r} (type: {type(model).__name__})")
        print(f"   tissue_types     = {mapped_tissue!r} (type: {type(mapped_tissue).__name__})")


        restitution_curves(
            base_dir=CASE_ROOT,
            output_folder=OUTPUT_FOLDER,
            ionic_model=model,          
            tissue_types=mapped_tissue,
            show_plots=SEE_PLOTS_DURING_POSTPROCESSING,
        )

print("\n✔ Post-processing completed\n")

def run_all():
    print("\n=======================================")
    print(" 🚀 S1–S2 FULL PIPELINE")
    print("=======================================\n")

    run_simulation()
    run_postprocessing()

    print("\n🎉 All tasks completed successfully!\n")


# ================= CLI =================
def main():
    if len(sys.argv) < 2:
        print("❌ ERROR: Please specify a command: sim | post | all")
        sys.exit(1)

    command = sys.argv[1].lower()

    if command == "sim":
        run_simulation()
    elif command == "post":
        run_postprocessing()
    elif command == "all":
        run_all()
    else:
        print(f"❌ Unknown command: {command}")
        print("Valid commands: sim | post | all")
        sys.exit(1)


if __name__ == "__main__":
    main()
