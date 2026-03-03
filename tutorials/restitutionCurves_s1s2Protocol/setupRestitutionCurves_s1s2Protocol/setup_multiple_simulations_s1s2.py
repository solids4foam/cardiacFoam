import subprocess
from pathlib import Path
from itertools import product
import shutil


# ================= FILE NAMES =================
ELECTRO_PROPERTIES_FILE_NAME = "electroProperties"


# ================= UPDATE FILES =================
def update_electro_properties(electro_file: Path, tissue: str, ionic_model: str) -> None:
    """
    Update the tissue and ionicModel entries in the properties file.
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


def update_s1s2_protocol(
    electro_file: Path,
    s1_interval: int,
    s2_interval: int,
    stim_amplitude: float,
    n_s1: int = 10,
    n_s2: int = 2,
) -> None:
    """
    Update the S1–S2 protocol entries in the electroProperties file.
    """
    lines = electro_file.read_text().splitlines(keepends=True)

    with electro_file.open("w") as f:
        for line in lines:
            stripped = line.strip()

            if stripped.startswith("stim_period_S1") and not stripped.startswith("//"):
                f.write(f"stim_period_S1  {s1_interval};\n")
            elif stripped.startswith("nstim1") and not stripped.startswith("//"):
                f.write(f"nstim1  {n_s1};\n")
            elif stripped.startswith("stim_period_S2") and not stripped.startswith("//"):
                f.write(f"stim_period_S2  {s2_interval};\n")
            elif stripped.startswith("nstim2") and not stripped.startswith("//"):
                f.write(f"nstim2  {n_s2};\n")
            elif stripped.startswith("stim_amplitude") and not stripped.startswith("//"):
                f.write(f"stim_amplitude   {stim_amplitude};\n")  # 👈 HERE
            else:
                f.write(line)

def copy_filter_output(output_folder_name: str, extension: str = ".txt"):
    """
    Copies all files with the given extension from the parent directory
    into a subfolder (created if needed) inside the parent directory,
    then deletes the originals after successful copy.
    """
    parent_dir = Path.cwd().parent
    output_folder = parent_dir / output_folder_name

    output_folder.mkdir(parents=True, exist_ok=True)

    print(f"\n📁 Scanning parent directory: {parent_dir}")
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


# ================= MAIN SETUP ENTRY =================
def run_s1s2_simulation(
    base_dir: Path,
    tissue_type: str,
    ionic_model: str,
    s1_interval: int,
    s2_intervals: list[int],
    stim_amplitude: float,
    run_cases_script: Path,
    output_folder: str,
):
    """
    Runs S1–S2 OpenFOAM single-cell simulations and collects outputs.
    OpenFOAM-safe: does NOT change cwd, case layout, or solver behavior.
    """

    print(f"\n🚀 Base (OpenFOAM case) directory: {base_dir}")

    electro_file = base_dir / "constant" / ELECTRO_PROPERTIES_FILE_NAME

   
    total = len(s2_intervals)
    for i, s2_interval in enumerate(s2_intervals, start=1):
        print("\n==========================================")
        print(f"▶ Simulation {i}/{total}")
        print(f"   Tissue type : {tissue_type}")
        print(f"   Ionic model : {ionic_model}")
        print(f"   S1–S2       : {s1_interval} ms – {s2_interval} ms")
        print("==========================================\n")


        update_electro_properties(electro_file, tissue_type, ionic_model)
        update_s1s2_protocol(electro_file, s1_interval, s2_interval, stim_amplitude)

        subprocess.run(
            ["bash", "-l", str(run_cases_script), str(base_dir)],
            check=True
        )

    print("\n✅ All S1–S2 single-cell runs completed.")

    # Collect outputs using your trusted copy logic
    copy_filter_output(output_folder)


# ================= SAFETY GUARD =================
if __name__ == "__main__":
    raise SystemExit(
        "This module is meant to be imported and called from mainS1-S2protocol.py"
    )
