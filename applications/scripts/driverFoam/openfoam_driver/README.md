# OpenFOAM tutorial driver architecture

This package provides a shared Python automation engine for tutorial sweeps,
OpenFOAM execution, and post-processing.

## Package structure

```text
openfoam_driver/
├── cli.py                          # CLI entrypoint (sim/post/all)
├── core/
│   ├── runtime/
│   │   ├── models.py              # TutorialSpec, CaseConfig contracts
│   │   ├── registry.py            # tutorial name -> make_spec factory
│   │   └── engine.py              # shared simulation/postprocess engine
│   └── defaults/                  # per-tutorial default parameters
├── specs/
│   ├── common.py                  # mutators/path helpers
│   └── tutorials/                 # tutorial-specific make_spec modules
├── postprocessing/                # shared postprocess runner + artifact manifest
├── scripts/run_case.sh            # OpenFOAM shell runner (Allclean/Allrun + options)
└── tests/                         # contract and regression tests for architecture
```

## TutorialSpec contract

Each tutorial module in `specs/tutorials/` builds a `TutorialSpec` with:

- `build_cases()`
- `apply_case(case_root, case)`
- `run_case(case_root, setup_root, case)`
- optional `collect_outputs(case_root, output_dir)`
- optional `postprocess(setup_root, output_dir)`

This keeps all tutorial workflows on one engine while allowing per-tutorial sweep logic.

## Registered tutorials

- `singleCell`
- `niederer2012`
- `manufacturedFDA`
- `restitutionCurves`

Aliases are handled in `core/runtime/registry.py`.

## Install and run

From repository root:

```bash
python3 -m pip install -e applications/scripts/driverFoam
```

Run examples:

```bash
# installed entrypoints
driverFoam sim --tutorial singleCell
foamctl all --tutorial niederer2012
foamctl sim --tutorial manufacturedFDA --dry-run

# module invocation
python3 -m openfoam_driver all --tutorial singleCell
```

## CLI actions

- `sim` : apply and run all planned cases
- `post`: run only post-processing hook
- `all` : `sim` then `post`

Useful flags:

- `--dry-run`
- `--continue-on-error`
- `--config <json>`
- `--tutorials-root <path>`

## Config override model

`--config` accepts either:

- top-level map keyed by tutorial name, or
- direct object with `make_spec(...)` keyword args for the selected tutorial.

Common override keys across tutorials:

- `case_dir_name`
- `setup_dir_name`
- `output_dir_name`
- `run_script_relpath`
- `postprocess_strict_artifacts`

Defaults live in `core/defaults/*.py`.

## Runtime artifacts

- `run_manifest.json`: written by the engine for every run.
- `plots.json`: written by postprocess runner, includes declared artifact metadata.

## Setup-folder dependencies

Expected per tutorial setup assets:

- `singleCell/setupSingleCell/singleCellinteractivePlots.py`
- `manufacturedFDA/setupManufacturedFDA/post_processing_manufactured.py`
- `NiedererEtAl2012/setupNiedererEtAl2012/postProcessing/{cache_postProcessing.py,line_postProcessing.py,points_postProcessing.py}`
- `restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing_restCurves.py`

## Architecture tests

Important contract tests:

- `tests/test_tutorial_architecture_contract.py`
- `tests/test_single_cell_contract.py`

These ensure registry coverage and required `make_spec(...)` keyword contract consistency.
