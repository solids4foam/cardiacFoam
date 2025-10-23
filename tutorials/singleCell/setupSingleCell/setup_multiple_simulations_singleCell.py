import subprocess
from pathlib import Path
import pandas as pd
from itertools import product
import shutil
from singleCellinteractivePlots import post_processing_single_cell


# Configuration
CONFIG = {
    "run_cases_script": Path.cwd() / "run_cases.sh",
     "output_folder": Path.cwd().parent / "outputSingleCell"
}

# Parameter sets
TISSUE_TYPES = ["mCells", "endocardialCells", "epicardialCells", "myocyte"]
IONIC_MODELS = ["BuenoOrovio", "TNNP", "Gaur", "Courtemanche"]
IONIC_MODELS = ["BuenoOrovio", "TNNP"]
TISSUE_TYPES = ["mCells"]


def update_electro_properties(electro_file: Path, tissue: str, ionic_model: str) -> None:
    """Update electroActivationProperties with tissue and ionic model."""
    lines = electro_file.read_text().splitlines(keepends=True)
    with electro_file.open("w") as f:
        for line in lines:
            if line.strip().startswith("tissue") and not line.strip().startswith("//"):
                f.write(f"tissue  {tissue};\n")
            elif line.strip().startswith("ionicModel") and not line.strip().startswith("//"):
                f.write(f"ionicModel  {ionic_model};\n")
            else:
                f.write(line)


import shutil
from pathlib import Path

def copy_filter_output(output_folder_name: str, extension: str = ".txt"):
    """
    Copies all files with the given extension from the parent directory
    into a subfolder (created if needed) inside the parent directory,
    then deletes the originals after successful copy.

    Args:
        output_folder_name (str): Name of the destination folder in the parent directory.
        extension (str): File extension filter (default: '.txt').
    """
    parent_dir = Path.cwd().parent
    output_folder = parent_dir / output_folder_name

    # Ensure the destination folder exists
    output_folder.mkdir(parents=True, exist_ok=True)

    print(f"\n📁 Scanning parent directory: {parent_dir}")
    files = list(parent_dir.glob(f"*{extension}"))

    if not files:
        print(f"⚠️ No '{extension}' files found in {parent_dir}")
        return

    for f in files:
        dest = output_folder / f.name
        try:
            shutil.copy2(f, dest)  # preserves metadata
            f.unlink()  # delete the original file
            print(f"✅ Moved: {f.name} → {dest}")
        except Exception as e:
            print(f"❌ Failed to move {f.name}: {e}")

    print(f"\n✅ All {extension} files moved to {output_folder}")


def run_single_cell_cases():
    base_dir = Path.cwd().parent
    print(f"Base directory: {base_dir}")
    tissue_path = base_dir / "constant" / "electroActivationProperties"

    combinations = list(product(TISSUE_TYPES, IONIC_MODELS))
    total = len(combinations)

    for i, (tissue, ionic_model) in enumerate(combinations, start=1):
        print("\n==========================================")
        print(f"▶ Simulation {i}/{total}")
        print(f"   Tissue type : {tissue}")
        print(f"   Ionic model : {ionic_model}")
        print("==========================================\n")

        update_electro_properties(tissue_path, tissue, ionic_model)

        # Run case
        subprocess.run(
            ["bash", "-l", str(CONFIG["run_cases_script"]), str(base_dir)],
            check=True
        )

    print("✅ All single-cell runs completed.")

if __name__ == "__main__":
    #run_single_cell_cases()
    #copy_filter_output(CONFIG["output_folder"])
    post_processing_single_cell(CONFIG["output_folder"])
    



