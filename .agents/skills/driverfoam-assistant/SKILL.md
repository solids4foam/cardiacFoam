---
name: driverfoam-assistant
description: >
  Use this skill to help a user build, validate, and run OpenFOAM parameter sweeps and optimization cases using the driverFOAM orchestrator. Trigger this when a user asks to "build a case", "run a sweep", "optimize an OpenFOAM simulation", or "set up driverFOAM".
---

# DriverFOAM End-User Assistant Skill

As an AI agent, your goal is to help users bridge the gap between their custom OpenFOAM simulation ideas and the `driverFOAM` automation engine. OpenFOAM cases are highly complex; you will use `driverFOAM`'s strict diagnostic planner to automatically ensure the physical correctness of the user's setup before running it.

**Architecture Status Context:**
- **Phase 1 (CI green):** Complete.
- **Phase 2 (core/plugin decoupling):** Complete (as of 2026-08-18). The `--plugin` flag is stable and the architecture is solver-agnostic.
- **Phase 3 (deterministic experiments & hybrid solvers):** In progress.
- **Catalog Maintenance:** `dict_key_allowlist.json` and the ionic model catalog are auto-gated by CI (failures on stale keys/renames). The `--plugin` selection affects the `capability_manifest` (allowed commands and samplable fields).

Follow this standard workflow when assisting a user with a new or existing case:

## 1. Case Scaffolding Workflow (Hybrid Approach)

**Do NOT build an OpenFOAM case from absolute scratch.** OpenFOAM requires a complex interplay of dictionaries (`fvSchemes`, `fvSolution`, `blockMeshDict`, boundary fields in `0/`).

When a user asks you to build a new case:
1. **Find a Base Tutorial:** Identify the closest existing tutorial in `tutorials/` or `applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/tutorials/` (e.g., `niederer2012` for the closest Bidomain equivalent, `singleCell` as the simplest Monodomain entry point).
2. **Copy the Scaffold:** Copy that tutorial folder to the user's requested location.
3. **Mutate the Scaffold:** Use your code editing tools to modify the `constant/electroProperties`, `system/controlDict`, or boundary conditions to match the user's specific request.

## 2. Sweep Generation

The user will usually want to run a parameter sweep (e.g., testing 3 different ionic models, or 5 different conductivity values). Entry-mode sweeps (using a registered tutorial's `make_spec`) differ from generic (from-scratch) sweeps, but the brain phase handles these differences.

A `sweep.json` has two top-level blocks:
- `base`: the fixed selectors that stay the same across every case (e.g. which
  registered `entry` to target, or `electro_selectors`/`physics_selectors`).
- `sweep`: `mode` (`cross_product` or `zip`) plus `independent`, a flat map of
  axis name → list of values to vary.

1. Create a `sweep.json` file inside the user's case directory.
2. Example of a `sweep.json` overriding the ionic model on the `singleCell` entry
   (`ionic_model` and `tissue` must be given together, and each must be a value
   the target entry's `make_spec()` actually accepts — check
   `IONIC_MODEL_TISSUE_MAP` for valid pairs, or use `driverFoam describe` on
   the entry to see its accepted kwargs):
   ```json
   {
       "base": {
           "entry": "singleCell"
       },
       "sweep": {
           "mode": "zip",
           "independent": {
               "ionic_model": ["Courtemanche", "TNNP", "BuenoOrovio"],
               "tissue": ["myocyte", "epicardialCells", "epicardialCells"]
           }
       }
   }
   ```
3. The sweep is not run as part of `driverFoam plan` — it is run separately
   with `driverFoam sweep-run` (see Section 4).

## 3. The Strict Diagnostics Loop (Auto-Repair)

This is your superpower. Before running the actual simulation, you MUST validate the physics and dictionaries using the `driverFOAM` strict planner.

1. **Run the Planner:**
   Run the following command from the terminal:
   ```bash
   driverFoam plan --strict --entry <registered_entry_name>
   ```
   Note: `--spec` does not belong with `plan` — it is only used with
   `sweep-plan`/`sweep-run` (see Section 4).
2. **Parse the Diagnostics:**
   The strict planner will output a structured JSON report. It will check if the chosen solver supports the chosen ionic model, if required fields like `defaultFieldValues` are present, and if the dictionary groups are complete.
3. **Auto-Repair:**
   If the plan fails, **do not just show the error to the user.**
   Read the JSON error output, open the user's dictionary (e.g., `constant/electroProperties`), and add or fix the missing entries yourself using the active plugin's dictionary catalog as your reference.
4. **Loop until Green:**
   Re-run `driverFoam plan --strict` until the case passes 100% of the diagnostics.

## 4. Execution

Once the strict plan passes, execute the sweep:

1. **Run the Sweep:**
   ```bash
   driverFoam sweep-run --spec sweep.json --output-dir ./sweep_output/
   ```
   Note: `--entry` is not valid with sweep actions (the target entry lives in
   `sweep.json`'s `base` block instead); `--output-dir` is required.
2. **Wait for Terminal State:**
   The post-processing phase (Step 5) will be executed once the sweep reaches a terminal state and all cases are completed. Artifact tracking (where OpenFOAM puts artifacts vs where they are realized) will be verified by the "brain".

## 5. Post-Processing Phase (Brain + Module)

The execution engine hands off to the post-processing phase once a sweep or DAG finishes. This is deliberately split into two independent pieces:

### The Brain
The brain (`build_sweep_context` in `postprocess_phase.py`) grounds the execution data. It:
- Reads `sweep_manifest.json` (started at, finished at, status, axis values).
- Verifies output against the actual disk state.
- Resolves case output directories (entry-mode vs. generic-mode differences are handled automatically).
- Extracts the `setup_root` from each case's `run_document.json`.
- Outputs a grounded `SweepContext` which serves as the single source of truth.

### The Post-Processing Module
The post-processing module (`run_postprocessing_module`) receives the `SweepContext` and a specific reasoning task. **It never re-reads the manifest or re-derives file locations.** 
- If any cases failed during execution, you should explain why and skip post-processing.
- The module lists available scripts using `list_postprocess_scripts()` and dispatches the task to the appropriate script based on its parsed description.

**Available Query Functions:**
You can query the brain for deeper reasoning without bypassing its verification:
- `read_case_workflow_state(context, case_id)`: Read full per-step status, diagnostics, and artifacts for a case.
- `read_case_output_file(context, case_id, relative_path)`: Safely read file content (only paths already verified by the brain).

### PostprocessingProtocol & Script Authoring
Every tutorial post-processing script must expose a `run_postprocessing` function matching the `PostprocessingProtocol` signature:

```python
from openfoam_driver.postprocessing import PostprocessingProtocol

def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    **kwargs: object,
) -> list[dict]:
    # Extract description from docstring for agent discovery
    '''Loads summary CSVs and plots them using the default palette.'''
    pass
```

The catalog extracts the script's `description` directly from this docstring as a `PostprocessScriptInfo` object, allowing the reasoning agent to decide if the script applies to the task.

**Available Utilities (`openfoam_driver.postprocessing`):**
- **Plotly Declarative Traces:** `PlotSpec`, `TraceSpec`, `build_line_traces`
- **Data Loading:** `load_csv_folder` (bulk CSV reading)
- **Styling:** `apply_plotly_layout`, `write_plotly_html` (Plotly) and `configure_matplotlib_defaults`, `finalize_matplotlib_figure` (Matplotlib).
- **Colors:** `DEFAULT_PALETTE`, `GroupShadedColors`
