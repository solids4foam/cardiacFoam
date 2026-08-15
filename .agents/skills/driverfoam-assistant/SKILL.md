---
name: driverfoam-assistant
description: >
  Use this skill to help a user build, validate, and run OpenFOAM parameter sweeps and optimization cases using the driverFOAM orchestrator. Trigger this when a user asks to "build a case", "run a sweep", "optimize an OpenFOAM simulation", or "set up driverFOAM".
---

# DriverFOAM End-User Assistant Skill

As an AI agent, your goal is to help users bridge the gap between their custom OpenFOAM simulation ideas and the `driverFOAM` automation engine. OpenFOAM cases are highly complex; you will use `driverFOAM`'s strict diagnostic planner to automatically ensure the physical correctness of the user's setup before running it.

Follow this standard workflow when assisting a user with a new or existing case:

## 1. Case Scaffolding Workflow (Hybrid Approach)

**Do NOT build an OpenFOAM case from absolute scratch.** OpenFOAM requires a complex interplay of dictionaries (`fvSchemes`, `fvSolution`, `blockMeshDict`, boundary fields in `0/`).

When a user asks you to build a new case:
1. **Find a Base Tutorial:** Identify the closest existing tutorial in `tutorials/` or `applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/tutorials/` (e.g., `niederer2012` for the closest Bidomain equivalent, `singleCell` as the simplest Monodomain entry point).
2. **Copy the Scaffold:** Copy that tutorial folder to the user's requested location.
3. **Mutate the Scaffold:** Use your code editing tools to modify the `constant/electroProperties`, `system/controlDict`, or boundary conditions to match the user's specific request.

## 2. Sweep Generation

The user will usually want to run a parameter sweep (e.g., testing 3 different ionic models, or 5 different conductivity values).

A `sweep.json` has two top-level blocks:
- `base`: the fixed selectors that stay the same across every case (e.g. which
  registered `entry` to target, or `electro_selectors`/`physics_selectors`).
- `sweep`: `mode` (`cross_product` or `zip`) plus `independent`, a flat map of
  axis name → list of values to vary.

1. Create a `sweep.json` file inside the user's case directory.
2. Example of a `sweep.json` overriding the ionic model on the `singleCell` entry
   (`ionic_model` and `tissue` must be given together, and each must be a value
   the target entry's `make_spec()` actually accepts — check
   `IONIC_MODEL_TISSUE_MAP` for valid pairs, or use `foamctl describe` on the
   entry to see its accepted kwargs):
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
3. The sweep is not run as part of `foamctl plan` — it is run separately with
   `foamctl sweep-run` (see Section 4).

## 3. The Strict Diagnostics Loop (Auto-Repair)

This is your superpower. Before running the actual simulation, you MUST validate the physics and dictionaries using the `driverFOAM` strict planner.

1. **Run the Planner:**
   Run the following command from the terminal:
   ```bash
   foamctl plan --strict --entry <registered_entry_name>
   ```
   Note: `--spec` does not belong with `plan` — it is only used with
   `sweep-plan`/`sweep-run` (see Section 4).
2. **Parse the Diagnostics:**
   The strict planner will output a structured JSON report. It will check if the chosen solver supports the chosen ionic model, if required fields like `defaultFieldValues` are present, and if the dictionary groups are complete.
3. **Auto-Repair:**
   If the plan fails, **do not just show the error to the user.**
   Read the JSON error output, open the user's dictionary (e.g., `constant/electroProperties`), and add or fix the missing entries yourself using the active plugin's dictionary catalog as your reference.
4. **Loop until Green:**
   Re-run `foamctl plan --strict` until the case passes 100% of the diagnostics.

## 4. Execution and Summary

Once the strict plan passes, execute the sweep:

1. **Run the Sweep:**
   ```bash
   foamctl sweep-run --spec sweep.json --output-dir ./sweep_output/
   ```
   Note: `--entry` is not valid with sweep actions (the target entry lives in
   `sweep.json`'s `base` block instead); `--output-dir` is required.
2. **Analyze Artifacts:**
   After the run completes, read the `artifacts_manifest.json` file generated in the output directory.
3. **Present to the User:**
   Summarize the results for the user. Tell them exactly where their VTK files, ECG traces, or CSV summaries were generated based on the artifact manifest, and highlight any interesting findings if applicable.
