from pathlib import Path
import subprocess
import sys

# Try to import myokit directly for metadata extraction
try:
    import myokit
    HAS_MYOKIT_LIB = True
except ImportError:
    HAS_MYOKIT_LIB = False

from mapping_engine import run_mapping
from transformer import run_transformation
import discovery

# Ordered pipeline stages
STAGES = ["cellml", "mmt", "ansic", "openfoam"]


def detect_start_stage(input_path: Path) -> str:
    if input_path.suffix == ".cellml":
        return "cellml"
    if input_path.suffix == ".mmt":
        return "mmt"
    if input_path.name == "sim.c":
        return "ansic"

    raise ValueError(f"Cannot infer stage from input: {input_path}")


def validate_pipeline(start: str, end: str):
    if start not in STAGES:
        raise ValueError(f"Invalid start stage: {start}")
    if end not in STAGES:
        raise ValueError(f"Invalid end stage: {end}")
    if STAGES.index(start) > STAGES.index(end):
        raise ValueError(f"Invalid pipeline: {start} → {end}")


def resolve_pipeline(start: str, end: str):
    i0 = STAGES.index(start)
    i1 = STAGES.index(end)
    return STAGES[i0 : i1 + 1]


# cellml → mmt  (Myokit)
# ------------------------------------------------------------
def run_cellml_to_mmt(cellml: Path, verbose=False) -> Path:
    if not cellml.exists():
        raise FileNotFoundError(cellml)

    mmt = cellml.with_suffix(".mmt")

    if verbose:
        print("    Running Myokit (cellml → mmt):")
        print(f"      input  : {cellml}")
        print(f"      output : {mmt}")

    cmd = [
        "myokit",
        "import",
        "cellml",
        str(cellml),
        str(mmt),
    ]

    subprocess.run(cmd, check=True)
    return mmt


# mmt → ansic (Myokit)
# ------------------------------------------------------------
def run_mmt_to_ansic(mmt: Path, verbose=False) -> Path:
    if not mmt.exists():
        raise FileNotFoundError(mmt)

    ansic_dir = Path("ansic")
    ansic_dir.mkdir(exist_ok=True)

    if verbose:
        print("    Running Myokit (mmt → ansic):")
        print(f"      input  : {mmt}")
        print(f"      output : {ansic_dir}/sim.c")

    cmd = [
        "myokit",
        "export",
        "ansic",
        str(mmt),
        str(ansic_dir),
    ]

    subprocess.run(cmd, check=True)

    sim_c = ansic_dir / "sim.c"
    if not sim_c.exists():
        raise FileNotFoundError("Myokit export failed: ansic/sim.c not found")

    return sim_c


# Metadata Extraction (Automated state_map and Discovery)
# ------------------------------------------------------------
def extract_metadata_from_mmt(mmt_path: Path) -> tuple[dict, dict]:
    """
    Automates the state_map and component discovery using Myokit.
    Returns (mapping, discovered_vars)
    """
    if not HAS_MYOKIT_LIB:
        print("    [Warning] Myokit Python library not found. Falling back to manual workflows.")
        return None, {}

    try:
        model = myokit.load_model(str(mmt_path))
        
        # 1. State Mapping
        mapping = {
            i: s.name().replace('.', '_') 
            for i, s in enumerate(model.states())
        }
        
        # 2. Component Discovery (Vm, Iion, Istim)
        discovered = discovery.discover_all(model)
        
        return mapping, discovered
        
    except Exception as e:
        print(f"    [Error] Failed to extract metadata from {mmt_path}: {e}")
        return None, {}


# Stage: ansic → openfoam (mapping_variables.py)
# ------------------------------------------------------------
def run_pipeline(
    input_path: Path,
    start: str,
    end: str,
    model: "Optional[str]" = None,
    verbose: bool = False,
):
    validate_pipeline(start, end)
    stages = resolve_pipeline(start, end)

    current = input_path
    auto_mapping = None
    discovered_vars = {}

    # Track mmt path if we pass through that stage to extract metadata
    mmt_ref = None
    if input_path.suffix == ".mmt":
        mmt_ref = input_path

    for i in range(len(stages) - 1):
        src = stages[i]
        dst = stages[i + 1]

        if verbose:
            print(f"      {src} → {dst}")

        if src == "cellml" and dst == "mmt":
            current = run_cellml_to_mmt(current, verbose)
            mmt_ref = current

        elif src == "mmt" and dst == "ansic":
            # Before exporting to C, capture the state ordering and discover components
            auto_mapping, discovered_vars = extract_metadata_from_mmt(mmt_ref)
            current = run_mmt_to_ansic(current, verbose)

        elif src == "ansic" and dst == "openfoam":
            if model is None:
                raise ValueError("Model name required for OpenFOAM generation. Use --model <ModelName_Year>")

            if "_" not in model:
                raise ValueError("--model must be ModelName_Year")

            name, year = model.rsplit("_", 1)
            
            # 1. Apply Python-based source-to-source rewrites (Replaces Coccinelle)
            if verbose:
                print("    Applying Python-based source-to-source rewrites (replacing Coccinelle)")
            run_transformation(current, verbose=verbose)

            # 2. Generate OpenFOAM-ready code
            if verbose:
                print(f"    Generating OpenFOAM code (discovered: {discovered_vars})")
            
            output_h = f"{name}_{year}.H"
            run_mapping(
                str(current), 
                output_h, 
                mapping=auto_mapping, 
                discovered=discovered_vars, # Pass discovered components
                verbose=verbose
            )
            return {
                "model": name,
                "year": year,
                "header": f"{name}_{year}Names.H",
                "source": output_h
            }

    return stages
