import subprocess
from pathlib import Path
import pandas as pd
import shutil
from singleCellinteractivePlots import post_processing_single_cell

# -------------------------------------------------------------------
# Paths / Configuration
# -------------------------------------------------------------------

# Folder where this Python script lives (setupSingleCell/)
SCRIPT_DIR = Path(__file__).resolve().parent

# Base OpenFOAM case directory (parent of setupSingleCell/)
BASE_DIR = SCRIPT_DIR.parent

CONFIG = {
    # run_cases.sh is in setupSingleCell/
    "run_cases_script": SCRIPT_DIR / "run_cases.sh",

    # where to collect all .txt outputs
    "output_folder": BASE_DIR / "singleCellOutputs",
}

# -------------------------------------------------------------------
# Ionic models and their tissue types
# (Make sure these names match what your OpenFOAM dictionaries expect!)
# -------------------------------------------------------------------


# Ionic models and their tissue types
IONIC_MODELS = ["TNNP", "Gaur", "Courtemanche", "BuenoOrovio"]
IONIC_MODEL_TISSUE_MAP = {
    "TNNP": ["epicardialCells", "mCells", "endocardialCells"],
    "BuenoOrovio": ["epicardialCells", "mCells", "endocardialCells"],
    "Gaur": ["myocyte"],
    "Courtemanche": ["myocyte"]}


STIMULUS_MAP = {
    "TNNP": 60,
    "BuenoOrovio": 0.4,
    "Gaur": 60,
    "Courtemanche": 60 }


# Name of the properties file inside constant/
# Change this to "electroActivationProperties" if that's the actual filename.
PROPERTIES_FILE_NAME = "cardiacProperties"
STIMULUS_FILE_NAME = "stimulusProtocol"


# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

def update_electro_properties(electro_file: Path, tissue: str, ionic_model: str) -> None:
    """
    Update the tissue and ionicModel entries in the properties file.

    It assumes lines like:
        tissue      something;
        ionicModel  something;
    and replaces them with the selected tissue / ionic model.
    """
    if not electro_file.exists():
        raise FileNotFoundError(f"Properties file not found: {electro_file}")

    lines = electro_file.read_text().splitlines(keepends=True)

    with electro_file.open("w") as f:
        for line in lines:
            stripped = line.strip()

            if stripped.startswith("tissue") and not stripped.startswith("//"):
                f.write(f"tissue      {tissue};\n")
            elif stripped.startswith("ionicModel") and not stripped.startswith("//"):
                f.write(f"ionicModel  {ionic_model};\n")
            else:
                f.write(line)

def update_stimulus_protocol(stimulus_file: Path, ionic_model: str) -> None:
    """
    Update the stim_amplitude entry in the stimulus protocol file
    based on the ionic model using STIMULUS_MAP.

    It assumes a line like:
        stim_amplitude   value;
    """

    if ionic_model not in STIMULUS_MAP:
        raise KeyError(
            f"Ionic model '{ionic_model}' not found in STIMULUS_MAP. "
            f"Available models: {list(STIMULUS_MAP.keys())}"
        )

    if not stimulus_file.exists():
        raise FileNotFoundError(f"Stimulus protocol file not found: {stimulus_file}")

    stim_amp = STIMULUS_MAP[ionic_model]

    lines = stimulus_file.read_text().splitlines(keepends=True)

    with stimulus_file.open("w") as f:
        for line in lines:
            stripped = line.strip()

            if stripped.startswith("stim_amplitude") and not stripped.startswith("//"):
                f.write(f"stim_amplitude   {stim_amp};\n")
            else:
                f.write(line)


def copy_filter_output(output_folder: Path, extension: str = ".txt") -> None:
    """
    Copy all files with given extension from the base case directory
    into 'output_folder', then delete originals.
    """
    parent_dir = BASE_DIR
    output_folder.mkdir(parents=True, exist_ok=True)

    print(f"\n📁 Scanning base directory for '{extension}' files: {parent_dir}")
    files = list(parent_dir.glob(f"*{extension}"))

    if not files:
        print(f"⚠️ No '{extension}' files found in {parent_dir}")
        return

    for f in files:
        dest = output_folder / f.name
        try:
            shutil.copy2(f, dest)
            f.unlink()
            print(f"✅ Moved: {f.name} → {dest}")
        except Exception as e:
            print(f"❌ Failed to move {f.name}: {e}")

    print(f"\n✅ All {extension} files moved to {output_folder}")


# -------------------------------------------------------------------
# Main driver
# -------------------------------------------------------------------

def run_single_cell_cases() -> None:
    print(f"Base (OpenFOAM case) directory: {BASE_DIR}")

    # Properties file inside the base case directory
    tissue_file = BASE_DIR / "constant" / PROPERTIES_FILE_NAME
    stimulus_file = BASE_DIR / "constant" / STIMULUS_FILE_NAME
    print(f"Using properties file: {tissue_file}")

    # Build (tissue, model) combinations dynamically
    pairs = []
    for model in IONIC_MODELS:
        for tissue in IONIC_MODEL_TISSUE_MAP[model]:
            pairs.append((tissue, model))

    print(f"\nRunning {len(pairs)} simulations...\n")

    for i, (tissue, model) in enumerate(pairs, start=1):
        print("==========================================")
        print(f"▶ Simulation {i}/{len(pairs)}")
        print(f"   Tissue type : {tissue}")
        print(f"   Ionic model : {model}")
        print("==========================================\n")

        # 1) Update constant/<properties file> with current tissue + model
        update_electro_properties(tissue_file, tissue, model)
        update_stimulus_protocol(stimulus_file, model)

        # 2) Call run_cases.sh with BASE_DIR (which contains Allclean/Allrun)
        try:
            subprocess.run(
                ["bash", "-l", str(CONFIG["run_cases_script"]), str(BASE_DIR)],
                check=True
            )
        except subprocess.CalledProcessError as e:
            print("❌ Error while running run_cases.sh")
            print("   Command:", e.cmd)
            print("   Return code:", e.returncode)
            # Optional: re-raise if you want the script to stop
            raise

    print("\n✅ All single-cell runs completed.")


# -------------------------------------------------------------------

if __name__ == "__main__":
    run_single_cell_cases()
    copy_filter_output(CONFIG["output_folder"])
    post_processing_single_cell(CONFIG["output_folder"])
