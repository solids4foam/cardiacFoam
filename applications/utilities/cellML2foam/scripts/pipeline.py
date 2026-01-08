from pathlib import Path
import subprocess

from mapping_engine import run_mapping

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


# Stage: ansic → ansic-renamed (Coccinelle)
# ------------------------------------------------------------
def run_coccinelle(sim_c: Path, cocci_file: Path, verbose=False):
    if not sim_c.exists():
        raise FileNotFoundError(sim_c)

    if verbose:
        print("    Running Coccinelle:")
        print(f"      sim.c : {sim_c}")
        print(f"      cocci : {cocci_file}")

    cmd = [
        "spatch",
        "--in-place", 
        str(cocci_file),
        str(sim_c),
    ]

    result = subprocess.run(
        cmd,
        capture_output=not verbose,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "Coccinelle failed:\n" + (result.stderr or "")
        )


# Stage: ansic→ openfoam (mapping_variables.py)
# ------------------------------------------------------------
def run_mapping_variables(sim_c: Path, model: str, year: str, verbose=False):
    state_map = Path("state_map.txt")
    if not state_map.exists():
        raise FileNotFoundError(
            "state_map.txt not found.\n"
            "Renaming variables requires a state_map.txt file "
            "in the working directory."
        )

    output_h = f"{model}_{year}.H"

    if verbose:
        print("    Running mapping_variables.py:")
        print(f"     sim.c : {sim_c}")
        print(f"     state_map     : {state_map}")
        print(f"     output        : {output_h}")

    run_mapping(str(sim_c), output_h, verbose=verbose)


# ------------------------------------------------------------
# Pipeline driver
# ------------------------------------------------------------
def run_pipeline(
    input_path: Path,
    start: str,
    end: str,
    model: str | None = None,
    verbose: bool = False,
):
    validate_pipeline(start, end)
    stages = resolve_pipeline(start, end)

    script_dir = Path(__file__).parent
    cocci_dir = script_dir / "cocci_executable"
    rename_cocci = cocci_dir / "renameVectors.cocci"

    current = input_path

    for i in range(len(stages) - 1):
        src = stages[i]
        dst = stages[i + 1]

        if verbose:
            print(f"      {src} → {dst}")

        # cellml → mmt
        if src == "cellml" and dst == "mmt":
            current = run_cellml_to_mmt(current, verbose)

        elif src == "mmt" and dst == "ansic":
            current = run_mmt_to_ansic(current, verbose)

        elif src == "ansic" and dst == "openfoam":
            state_map = Path("state_map.txt")
            if not state_map.exists():
                print(
                    f"""
                ============================================================
                STATE MAP REQUIRED
                ============================================================

                To continue to OpenFOAM code generation, you must provide a
                'state_map.txt' file.

                The state map defines the ordering of state variables and
                MUST be created by the user based on the Myokit model output.

                Example:

                If your .mmt (or generated code) contains initial values like:

                    membrane.v               = -89.74808
                    CaMK.CaMKt               =  1.095026e-2
                    intracellular_ions.nai   = 12.39736
                    intracellular_ions.nass  = 12.3977
                    intracellular_ions.ki    = 147.7115
                    ...

                Then your state_map.txt should look like:

                    0   V
                    1   CaMKt
                    2   Nai
                    3   Nass
                    4   Ki
                    ...

                Steps to proceed:

                1. Open the generated .mmt or ansic/sim.c file
                2. Locate the list of initial state values
                3. Assign indices in the same order (starting from 0)
                4. Create 'state_map.txt' in this directory
                5. Re-run the tool using:

            
                ============================================================
                
                Please name the output H file using the following format:

                    <ModelName>_<Year>.H

                Example:
                    TorORD_dynCl_2023.H

                From this filename, the tool will automatically infer:
                - Ionic model name
                - Year of publication
                - Names of generated symbols and header files

                The corresponding header will be written as:
                    TorORD_dynCl_2023Names.H

                Notes:
                - Do not use spaces in the model name
                - The year must be numeric

                ============================================================
                
                To finish the process, please write the command:
                
                ./cellML2foam --from ansic --to openfoam --model <ModelName_Year> ansic/sim.c 
                ./cellML2foam --from mmt --to openfoam --model <ModelName_Year>  <ModelName>.mmt

                Example:

                ./cellML2foam --from ansic --to openfoam --model ToRORd_2023 ansic/sim.c

                """
                )
                return None

            if model is None:
                raise ValueError(
                    "Model name required for OpenFOAM generation.\n"
                    "Use --model <ModelName_Year>"
                )

            if "_" not in model:
                raise ValueError(
                    "--model must be of the form <ModelName>_<Year>"
                )

            name, year = model.rsplit("_", 1)
            if not year.isdigit():
                raise ValueError("Year in --model must be numeric")

            # --------------------------------------------------------
            # Always apply Coccinelle rewrite BEFORE mapping
            # --------------------------------------------------------
            if verbose:
                print("    Applying source-to-source rewrites")

            run_coccinelle(
                sim_c=current,
                cocci_file=rename_cocci,
                verbose=verbose,
            )

            # --------------------------------------------------------
            # Generate OpenFOAM-ready code
            # --------------------------------------------------------
            result = run_mapping_variables(
                sim_c=current,
                model=name,
                year=year,
                verbose=verbose,
            )
            return result

    return stages
