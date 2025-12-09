"""
The setup and analysis of convergence tests for simulations using a manufactured solution
for verification of the Electrophysiology solver with the input of dx, dt and dimension.

Main workflow:
--------------
1. User selects parameter combinations from the GUI.
2. For each combination:
   - Update `blockMeshDict."dimension"` with the corresponding mesh resolution.
   - Update `controlDict` with the chosen Δt.
   - Update 'cardiacProperties' with the selected dimension.
3. Runs the case using a helper shell script (`run_cases.sh`).
4. Solves the errors in the solvers and outputs .dat file.

"""

import subprocess
from itertools import product
import sys
import time
import pandas as pd
from pathlib import Path






def type_text(text: str, delay: float = 0.02) -> None:
    """Print text with a typing effect."""
    for char in text:
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)
    print()



def select_simulation_parameters(
    NUMBER_CELLS, DT_VALUES, piecewise: bool = True) -> pd.DataFrame:
    """
    Generate simulation parameter combinations without GUI.

    Returns a DataFrame with columns:
    ["Δt[ms]", "Cells"]
    """
    if piecewise:
        combos = [
            (dt, dx)
            for (dt, dx) in zip(DT_VALUES, NUMBER_CELLS)
        ]
        mode = "Piecewise (1-to-1 pairing)"
    else:
        combos = list(product(DT_VALUES, NUMBER_CELLS))
        mode = "All combinations"

    df = pd.DataFrame(combos, columns=["Δt[ms]","Cells"])
    print(f"⚙️  Generated simulation combinations ({mode}):")
    print(df)
    return df


def update_blockMeshDict_dx(blockMeshDict_path: Path, dx: float, dimension: str) -> None:
    """
    Replace the hex block line in blockMeshDict according to dx value.
    Automatically generates the hex line based on dx.
    """
    lines = blockMeshDict_path.read_text().splitlines(keepends=True)
    replaced = False

    if dimension == "1D":
        hex_line_template = "hex (0 1 2 3 4 5 6 7) ({dx} 1 1) simpleGrading (1 1 1)\n"
    elif dimension == "2D":
        hex_line_template = "hex (0 1 2 3 4 5 6 7) ({dx} {dx} 1) simpleGrading (1 1 1)\n"
    elif dimension == "3D":
        hex_line_template = "hex (0 1 2 3 4 5 6 7) ({dx} {dx} {dx}) simpleGrading (1 1 1)\n"
    else:
        raise ValueError(f"Unsupported dimension type: {dimension}")

    with blockMeshDict_path.open('w') as f:
        for line in lines:
            stripped = line.strip()
            if not replaced and stripped.startswith('hex (0 1 2 3 4 5 6 7)') and not stripped.startswith('//'):
                f.write(hex_line_template.format(dx=int(dx)))
                replaced = True
            else:
                f.write(line)


def update_solver(numericalPropertiesPath: Path, solverType: str) -> None:
    # Decide flag values
    if solverType == "explicit":
        solveExplicit_flag = "yes"
    else:
        solveExplicit_flag = "no"

    lines = numericalPropertiesPath.read_text().splitlines(keepends=True)

    with numericalPropertiesPath.open('w') as f:
        for line in lines:
            stripped = line.strip()
            # Update solveExplicit
            if stripped.startswith("solveExplicit") and not stripped.startswith("//"):
                indent = line[:len(line) - len(line.lstrip())]
                f.write(f"{indent}solveExplicit    {solveExplicit_flag};\n")
                continue
            else:
                f.write(line)


def update_controlDict_dt(controlDict_path: Path, dt_ms: float) -> None:
    """
    Update the deltaT value in controlDict file with the given dt (in milliseconds).
    """
    dt_seconds = dt_ms   # convert ms to seconds
    lines = controlDict_path.read_text().splitlines(keepends=True)
    with controlDict_path.open('w') as f:
        for line in lines:
            if 'deltaT' in line and not line.strip().startswith('//'):
                f.write(f"deltaT          {dt_seconds};\n")
            else:
                f.write(line)


def update_dimension(dimension_path: Path, dimension: str) -> None:
    """
    Update or insert the "dimension" entry in electroActivationProperties.
    Preserves indentation and ignores commented lines.
    """
    lines = dimension_path.read_text().splitlines(keepends=True)
    with dimension_path.open('w') as f:
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("dimension ") and not stripped.startswith("//"):
                indent = line[:len(line) - len(line.lstrip())]
                f.write(f"{indent}dimension \"{dimension}\";\n")
            else:
                f.write(line)


def run_simulation_cases(
    selected_combinations: pd.DataFrame,
    CONFIG,
    DIMENSIONS,
    SOLVER,
) -> None:
    """Run the simulation cases based on selected parameters.

    selected_combinations: DataFrame with at least two columns:
        - Δt (time step)
        - number of cells (or an equivalent spatial parameter)

    DIMENSIONS: list of strings, e.g. ["1D"] or ["1D", "2D", "3D"]
    SOLVER: list of solver names as strings, e.g. ["explicit", "implicit"]
    """

    base_dir = CONFIG["case_file"].parent
    print(f"Base case directory: {base_dir}")

    # Common paths (do not depend on dt or dx)
    controlDict_path = base_dir / "system" / "controlDict"
    dimension_path   = base_dir / "constant" / "cardiacProperties"
    solver_path      = base_dir / "constant" / "numericalsolverProperties"
    run_cases_path   = CONFIG["run_cases_script"]

    # List of (dt, cells) from the DataFrame
    param_list = selected_combinations.values.tolist()

    # Total number of runs = dt/dx combinations × dimensions × solvers
    total = len(param_list) * len(DIMENSIONS) * len(SOLVER)
    case_counter = 0

    # Loop over each dimension and solver, then over dt/dx pairs
    for dimension in DIMENSIONS:
        # BlockMesh dict depends on the dimension
        blockMeshDict_path = base_dir / "system" / f"blockMeshDict.{dimension}"

        for solver_val in SOLVER:
            # solver_val is a *string*: "explicit" or "implicit"
            for (dt_val, dx_val) in param_list:
                case_counter += 1
                print(f"\n🔄 Running simulation ({case_counter}/{total})")
                print(f"   ➤ Dimension : {dimension}")
                print(f"   ➤ Cells     : {dx_val}")
                print(f"   ➤ Δt        : {dt_val}")
                print(f"   ➤ Solver    : {solver_val}", flush=True)

                # --- Update OpenFOAM dictionaries ---
                update_blockMeshDict_dx(blockMeshDict_path, dx_val, dimension)
                update_controlDict_dt(controlDict_path, dt_val)
                update_dimension(dimension_path, dimension)
                update_solver(solver_path, solver_val)

                # --- Run the OpenFOAM case ---
                subprocess.run(
                    ["bash", "-l", str(run_cases_path), str(base_dir), str(dimension)],
                    check=True,
                )

    print("✅ All simulations done.", flush=True)
