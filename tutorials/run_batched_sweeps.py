import os
import shutil
import subprocess
import sys
from pathlib import Path

# Path setup
TUTORIALS_ROOT = Path(__file__).parent.resolve()
DRIVERFOAM_DIR = TUTORIALS_ROOT.parent / "applications" / "scripts" / "driverFoam"
sys.path.insert(0, str(DRIVERFOAM_DIR))

from openfoam_driver.ionic_model_catalog import BATCHED_MODELS, IONIC_MODEL_CATALOG

TEMPLATE_DIR = TUTORIALS_ROOT / "NiedererEtAl2011BatchedTemplate"
RUNS_DIR = TUTORIALS_ROOT / "comparisonResults" / "batched_runs"

if not RUNS_DIR.exists():
    RUNS_DIR.mkdir(parents=True)

print(f"Starting batched sweeps for {len(BATCHED_MODELS)} models using full matrix (cpu, euler, rl, soa)...")

results = []

for model in BATCHED_MODELS:
    case_dir = RUNS_DIR / f"NiedererEtAl2011Batched_{model}"
    
    # 1. Clean existing and copy template
    if case_dir.exists():
        shutil.rmtree(case_dir)
    
    shutil.copytree(TEMPLATE_DIR, case_dir)
    
    # Look up the first compatible tissue for this model
    catalog_entry = IONIC_MODEL_CATALOG.get(model)
    if not catalog_entry:
        print(f"Failed to find {model} in catalog.")
        results.append((model, "failed: missing from catalog"))
        continue
    
    tissue = catalog_entry.compatible_tissues[0] if catalog_entry.compatible_tissues else "myocyte"
    parent_model = model.replace("Batched", "")
    
    print(f"\n========================================================")
    print(f"--- Setting up {model} full sweep (tissue: {tissue}) ---")
    print(f"========================================================")
    
    # 2. Modify all electroProperties.* files in the cloned directory
    constant_dir = case_dir / "constant"
    for prop_file in constant_dir.glob("electroProperties.*"):
        content = prop_file.read_text()
        
        # Replace tissue
        content = content.replace("epicardialCells", tissue)
        
        # Replace model
        if prop_file.name == "electroProperties.cpu":
            # The CPU reference uses the non-batched parent model
            content = content.replace("BuenoOrovio", parent_model)
        else:
            # All other batched runs use the batched model
            content = content.replace("BuenoOrovioBatched", model)
            
        prop_file.write_text(content)
    
    # 3. Launch the full comparison suite script inside the case_dir
    print(f"--- Launching runComparisons.sh for {model} ---")
    try:
        run_cmd = ["bash", "runComparisons.sh"]
        
        # We explicitly set the environment variables so that any leaked
        # variables from the user's terminal (like a custom MODES array missing 'cpu')
        # don't break the scripts!
        env = os.environ.copy()
        env["MODES"] = "cpu batched_euler batched_rl batched_soa"
        env["SUBSTEPS"] = "50 20 10 5"
        
        # We stream output to stdout so you can watch it run since it takes a while
        proc = subprocess.run(run_cmd, cwd=case_dir, env=env)
        if proc.returncode != 0:
            print(f"Execution failed for {model} (exit code {proc.returncode})")
            results.append((model, f"failed with code {proc.returncode}"))
        else:
            print(f"Execution successful for {model}.")
            results.append((model, "complete"))
            
    except Exception as e:
        print(f"Failed to launch {model}: {e}")
        results.append((model, f"failed: {e}"))

print("\n--- Summary ---")
for model, status in results:
    print(f"{model}: {status}")


