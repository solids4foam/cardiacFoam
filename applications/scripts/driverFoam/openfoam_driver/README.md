# OpenFOAM tutorial driver architecture

This package provides a shared Python automation engine for tutorial sweeps,
workflow-case execution, and post-processing.

## Package structure

```text
openfoam_driver/
├── cli.py                          # CLI entrypoint (plan/describe/sim/post/all)
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
- `manufacturedFDABidomain`
- `manufacturedFDABathBidomain`
- `restitutionCurves`

In addition to these curated specs, the driver can also run:

- `genericCase` / `randomCase` with `case_dir_name` supplied in config
- any existing case folder directly, for example `foamctl sim --entry ECG`

## Install and run

From repository root:

```bash
python3 -m pip install -e applications/scripts/driverFoam
```

Run examples:

```bash
# installed entrypoints
foamctl all --entry niederer2012
foamctl sim --entry manufacturedFDA --dry-run
foamctl sim --entry manufacturedFDABidomain --dry-run
foamctl plan --strict --entry singleCell
foamctl run --strict --entry singleCell
foamctl step --strict --entry singleCell --step solve
driverFoam sim --entry singleCell

# module invocation
python3 -m openfoam_driver all --entry singleCell
python3 -m openfoam_driver describe --entry singleCell

# repo-local wrapper
applications/scripts/driverFoam/bin/driverFoam sim --entry ECG --dry-run
```

## CLI actions

- `sim` : apply and run all planned cases
- `post`: run only post-processing hook
- `all` : `sim` then `post`
- `describe` : print machine-readable entry/spec metadata as JSON
- `plan` : print a non-mutating strict machine-readable launch contract as JSON
- `step` : execute exactly one normalized strict-plan workflow step
- `run` : execute normalized strict-plan workflow steps until completion or failure

Useful flags:

- `--dry-run`
- `--continue-on-error`
- `--config <json>`
- `--strict` for `plan`, `step`, and `run`
- `--step <id>` for `step`
- `--tutorials-root <path>`

The `describe` action resolves the requested entry and prints:

- the `make_spec(...)` parameter schema and defaults
- resolved case/setup/output paths
- the planned cases for the current configuration
- the grouped dict-entry catalog for `physicsProperties` and `electroProperties`
- the launch plan for `sim`, `post`, and `all`, including the exact driver
  command and expected manifest path

The `plan --strict` action resolves the requested entry, validates the planned
RunDocument v2, checks dict-key catalog coverage, predicts data artifacts, and
exits non-zero if any machine-readable contract is incomplete.

In strict-plan output, `workflow_dag.steps[*]` is normalized for a future step
runner: `command` contains only the executable name, `args` contains argv
arguments, `cwd` is case-relative, `depends_on` is validated against known step
ids, and `produces` lists expected artifact ids when they can be attributed.
The companion `workflow_state` is an initial, non-executed state snapshot:
all steps are `pending`, attempts are `0`, logs and exit codes are `null`, and
`current_step_id` points at the first dependency-free step.

The low-level runner API
`openfoam_driver.core.runtime.workflow_runner.run_workflow_step(...)` executes
one normalized step, writes stdout/stderr logs, and returns an updated
`workflow_state`. It deliberately does not implement multi-step orchestration,
resume, or retry loops yet.

The `step --strict` action is the CLI wrapper around that low-level runner. It
runs only the requested step, writes `workflow_state.json` and
`workflow_logs/<step>.attempt<N>.*.log` under the strict-plan output directory,
prints the final state JSON, and exits non-zero if strict planning or the step
execution fails.

The `run --strict` action reads the same `workflow_state.json` if present and
executes the next runnable step until the workflow completes or a step fails.
It does not retry a failed saved state automatically; use `step --strict` for
explicit manual reruns.

## Config override model

`--config` accepts either:

- top-level map keyed by entry name, or
- direct object with `make_spec(...)` keyword args for the selected entry.

Canonical example config:

- `applications/scripts/driverFoam/openfoam_driver/spec_overrides.example.json`

Common override keys across tutorials:

- `case_dir_name`
- `setup_dir_name`
- `output_dir_name`
- `run_script_relpath`
- `electro_property_overrides`
- `physics_property_overrides`
- `postprocess_strict_artifacts`

Dictionary overrides accept either:

- a mapping of dotted paths to values, for example
  `$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S1: 750`
- or a list of `{key, value, scope}` objects

For `electro_property_overrides`, the special scope token
`$ELECTRO_MODEL_COEFFS` resolves to the currently selected
`<solver>Coeffs` sub-dictionary in `electroProperties`.

The driver override layer is generic: it does not hardcode a whitelist of
allowed keys. Any existing entry in `physicsProperties` or `electroProperties`
can be updated through the dotted-path override API.

Source-backed catalog of repository-known dict entries:

- `applications/scripts/driverFoam/openfoam_driver/dict_entries.py`

That catalog was assembled from the current `src/` readers and groups the
known override paths for:

- `physicsProperties.type`
- `electroProperties.myocardiumSolver`
- common `<solver>Coeffs` keys such as `ionicModel`, `tissue`,
  `solutionAlgorithm`, `writeAfterTime`, and `outputVariables`
- `singleCellStimulus`
- `externalStimulus`
- eikonal-diffusion keys
- ECG keys
- active-tension keys

The catalog also carries machine-readable metadata per entry:

- `value_kind`
- `enum_values`
- `dynamic_path`

Examples:

- `physics_property_overrides: { "type": "electroMechanicalModel" }`
- `electro_property_overrides: { "$ELECTRO_MODEL_COEFFS.ionicModel": "TNNP" }`
- `electro_property_overrides: { "$ELECTRO_MODEL_COEFFS.activeTensionModel.activeTensionModel": "GoktepeKuhl" }`
- `electro_property_overrides: { "$ELECTRO_MODEL_COEFFS.ecgDomains.ECG.electrodePositions.V1": "(1 2 3)" }`

Notes:

- Overrides update existing entries only. They do not insert brand new keys.
- For dimensioned scalars, vectors, and tensors, pass the full OpenFOAM literal
  as a string, for example
  `"[0 -3 0 0 0 1 0] 50000"` or
  `"[-1 -3 3 0 0 2 0] (0.133418 0 0 0.0176062 0 0.0176062)"`.
- The catalog lists repository-backed keys. Additional OpenFOAM ODE-solver
  pass-through keys may also be valid because the dictionaries are forwarded to
  `ODESolver::New(...)`.

Generic-case sweeps can also use a `cases` array. Each case may override
`electro_property_overrides`, `physics_property_overrides`, `dimension`,
`parallel`, `touch_case_foam`, and `openfoam_bashrc` on top of the spec defaults.

Defaults live in `core/defaults/*.py`.

## Runtime artifacts

- `run_manifest.json`: written by the engine for every run.
- `run_report.md`: human-readable summary written alongside the manifest.
- `plots.json`: written by postprocess runner, includes declared artifact metadata.

`run_manifest.json` is the current machine-facing run-state file for local-app
integration. The current schema includes:

- `schema_version`
- `run_id`
- `requested_action`
- `entry`
- `entry_kind`
- `entry_path`
- `source_type`
- `workflow_family`
- `status`
- `postprocess_status`
- `current_case_id`
- `started_at_utc`
- `updated_at_utc`
- `finished_at_utc`
- `total_cases`
- `planned_cases`
- `completed_cases`
- `failed_cases`
- `error`
- `plots_manifest_path`
- `human_report_path`
- `results`

Each `results` item currently reports:

- `case_id`
- `status`
- `duration_s`
- `params`
- `error`
- `index`
- `total_cases`
- `started_at_utc`
- `finished_at_utc`

## Setup-folder dependencies

Expected per tutorial setup assets:

- `singleCell/setupSingleCell/singleCellinteractivePlots.py`
- `manufacturedSolutions/monodomainPseudoECG/setupManufacturedFDA/post_processing_manufactured.py`
- `NiedererEtAl2011/NiedererEtAl2011verification/setupNiedererEtAl2011/postProcessing/{cache_postProcessing.py,line_postProcessing.py,points_postProcessing.py}`
- `restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing_restCurves.py`

## Architecture tests

Important contract tests:

- `tests/test_tutorial_architecture_contract.py`
- `tests/test_single_cell_contract.py`
- `tests/test_mutators.py`

These ensure registry coverage and required `make_spec(...)` keyword contract consistency.
