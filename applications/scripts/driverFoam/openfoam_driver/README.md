# OpenFOAM tutorial driver architecture

This package provides a shared Python automation engine for tutorial sweeps,
workflow-case execution, and post-processing.

## Package structure

```text
openfoam_driver/
├── cli.py                          # CLI entrypoint (plan/step/run/describe/sweep-plan/sweep-run)
├── strict_planning.py              # strict preflight contract report
├── core/
│   ├── runtime/
│   │   ├── run_model.py           # RunDocument v3 + explicit v1/v2 migration
│   │   ├── workflow.py            # workflow DAG normalization/validation
│   │   ├── workflow_state.py      # persisted step state model
│   │   ├── workflow_runner.py     # one-step strict subprocess runner
│   │   ├── workflow_orchestrator.py # runs a workflow_dag to completion
│   │   ├── models.py              # TutorialSpec, CaseConfig contracts
│   │   └── registry.py            # tutorial name -> make_spec factory
│   ├── defaults/                  # per-tutorial default parameters
│   ├── plugin_capabilities.py     # focused internal solver capability adapters
│   ├── compatibility.py           # behavior-preserving legacy boundary
│   └── COMPATIBILITY.md           # fallback reasons, activation, and tests
├── specs/
│   ├── common.py                  # mutators/path helpers
│   ├── mesh_geometry.py           # mesh geometry checks for 1D 3D
│   └── tutorials/                 # tutorial-specific make_spec modules
├── postprocessing/                # shared postprocess runner + artifact manifest
├── scripts/                       # OpenFOAM runner, catalog scanners, allowlists
│   ├── run_case.sh                # Allclean/Allrun wrapper
│   └── dict_key_allowlist.json    # strict scanner parser/pass-through allowlist
└── tests/                         # contract and regression tests for architecture
```

## Standalone Usage

`driverFOAM` can be used completely independently of the `cardiacFoam` C++ repository (e.g. cloned into a `/tmp` directory or installed as a standalone `pip` package).

This mode is designed for CI pipelines or users wanting to drive their own existing OpenFOAM simulations without keeping the entire `cardiacFoam` source tree around.

When run standalone (where `src/` and `tutorials/` siblings do not exist):

- `driverFOAM` automatically falls back to bundled template fixtures (like `electroProperties`) for schema parsing and planning.
- The C++ drift guards (tests that assert Python configurations match C++ `Names.H` headers) are gracefully skipped.
- Regression equivalence and verification tests (which expect physical tutorials on disk) are skipped.
- `strict` and `generic` planning work out of the box using your local simulation paths.

## TutorialSpec contract

Each tutorial module in `specs/tutorials/` builds a `TutorialSpec` with:

- `build_cases()`
- `apply_case(case_root, case)`
- `run_case(case_root, setup, case)`
- optional `collect_outputs(case_root, output_dir)`
- optional `postprocess(setup, output_dir)`

This keeps all tutorial workflows on one engine while allowing per-tutorial sweep logic.

## Solver injection boundary

The public `SolverPlugin` interface remains backward compatible. Internally, a
per-operation `DriverContext` adapts it into focused tutorial, dictionary,
validation, artifact, case-compatibility, C++ mapping, mesh-policy, RunDocument,
and sweep capabilities. Core code consumes those capabilities and does not
reach through the plugin object directly.

Omitting a context still selects cardiacFoam exactly as before. Legacy case,
sweep, mutation, and planning decisions are named in
[`core/COMPATIBILITY.md`](core/COMPATIBILITY.md); Plan 1 moves their ownership
without changing their activation or output.

## Registered tutorials

- `singleCell`
- `cable1DCVConvergence`
- `niederer2012`
- `manufacturedMonodomainPseudoECG`
- `manufacturedBidomain`
- `manufacturedBathBidomain`
- `restitutionCurves`

In addition to these curated specs, the driver can also run:

- `genericCase` / `randomCase` with `case_dir_name` supplied in config
- any existing case folder directly, for example `driverFoam run --strict --entry ECG`

## Install and run

From repository root:

```bash
python3 -m pip install -e applications/scripts/driverFoam
```

Run examples:

```bash
# installed entrypoints
driverFoam plan --strict --entry singleCell
driverFoam run --strict --entry singleCell
driverFoam step --strict --entry singleCell --step solve
driverFoam run --strict --entry niederer2012

# module invocation
python3 -m openfoam_driver run --strict --entry singleCell
python3 -m openfoam_driver describe --entry singleCell

# repo-local wrapper
applications/scripts/driverFoam/bin/driverFoam run --strict --entry ECG
```

## CLI actions

- `plan` : print a non-mutating strict machine-readable launch contract as JSON
- `step` : execute exactly one normalized strict-plan workflow step
- `run` : execute normalized strict-plan workflow steps until completion or failure
- `describe` : print machine-readable entry/spec metadata as JSON
- `sweep-plan` : dry-run every case in a `sweep.json` and report readiness
- `sweep-run` : materialize and run every case in a `sweep.json`

Useful flags:

- `--dry-run`
- `--continue-on-error`
- `--config <json>`
- `--strict` for `plan`, `step`, and `run`
- `--step <id>` for `step`
- `--tutorials-root <path>`

Before planning or running the cardiacFoam plugin, configure the backend and
build identity once for the host. Copy
`applications/scripts/driverFoam/driverfoam-runtime.example.yaml` to a
user-local file and set:

```bash
export DRIVERFOAM_RUNTIME_CONFIG=/absolute/path/driverfoam-runtime.yaml
driverFoam plan --strict --entry manufacturedMonodomainPseudoECG
```

The plugin declares `lightweight` and `full` backends in
`plugins/cardiacfoam/plugin.yaml`. The local selection chooses exactly one,
supplies a solids4foam root for `full`, and points to the build manifest
generated by the top-level `Allwmake`. driverFOAM validates the selected
backend against that manifest, including the OpenFOAM root, linked libraries,
and artifact hashes. It does not auto-discover a different solids4foam
checkout. The OpenFOAM base installation remains separate and is selected with
`--openfoam-bashrc` or `OPENFOAM_BASHRC`.

The `describe` action resolves the requested entry and prints:

- the `make_spec(...)` parameter schema and defaults
- resolved case/setup/output paths
- the resolved cases for the current configuration
- the grouped dict-entry catalog for `physicsProperties` and `electroProperties`
- `strict_launch`: the `run --strict --entry ...` command for this entry,
  including its case/output paths

## Strict autonomous contract

The `plan --strict` action resolves the requested entry, validates the resolved
RunDocument v3, checks dict-key catalog coverage, predicts data artifacts,
normalizes the workflow DAG, and exits non-zero if any machine-readable
contract is incomplete. It does not mutate case files.

Programmatic callers use:

```python
from openfoam_driver.core.strict_planning import strict_plan

report = strict_plan("singleCell")
payload = report.to_json()
```

The JSON report contains:

- `status`: `ok` or `failed`
- `resolved_entry`: selected entry, case root, setup root, output directory,
  entry kind, source type, and workflow family
- `readiness_score`: weighted 0-100 run-readiness score with blocked/warning
  stages
- `simulation_audit`: per-stage explanation of how simulations are created and
  prepared: `build_cases()`, required OpenFOAM files, dictionary resolution,
  workflow DAG normalization, artifact prediction, environment preflight, and
  mesh geometry
- `validation_diagnostics`: RunDocument/config validation
- `workflow_diagnostics`: DAG shape, dependency, command, cwd, and cycle
  diagnostics
- `catalog_coverage_errors`: strict C++ dict-key scanner failures
- `artifact_diagnostics`: unknown solvers, unknown workflow commands, missing
  utility `produces`, and empty predictions caused by missing catalog coverage
- `launch`: exact launch command and expected manifest path
- `workflow_dag`: normalized executable workflow steps
- `workflow_state`: initial or resumed state snapshot
- `expected_artifacts`: predicted data artifacts
- `run_document`: canonical RunDocument v3 payload

In strict-plan output, `workflow_dag.steps[*]` is normalized for the strict step
runner: `command` contains only the executable name, `args` contains argv
arguments, `cwd` is case-relative, `depends_on` is validated against known step
ids, and `produces` lists expected artifact ids when they can be attributed.
The normalizer rejects duplicate step ids, missing dependencies, cycles, unsafe
case-relative `cwd` paths, shell metacharacters in command names, and unknown
step status vocabulary.
The companion `workflow_state` is an initial, non-executed state snapshot:
all steps are `pending`, attempts are `0`, logs and exit codes are `null`, and
`current_step_id` points at the first dependency-free step.

The low-level runner API
`openfoam_driver.core.runtime.workflow_runner.run_workflow_step(...)` executes
one normalized step, writes stdout/stderr logs, and returns an updated
`workflow_state`. It accepts only `pending` or explicitly failed steps, checks
dependencies, confines `cwd` under the case root, records the attempt count,
and writes state atomically when a state path is supplied.

The `step --strict` action is the CLI wrapper around that low-level runner. It
runs only the requested step. If `workflow_state.json` already exists under the
strict-plan output directory, it loads that saved state before execution;
otherwise it starts from the strict-plan initial state. It writes
`workflow_state.json` and `workflow_logs/<step>.attempt<N>.*.log`, prints the
final state JSON, and exits non-zero if strict planning or the step execution
fails.

The `run --strict` action reads the same `workflow_state.json` if present and
executes the next runnable step until the workflow completes or a step fails.
It does not retry a failed saved state automatically; use `step --strict` for
explicit manual reruns.

Strict execution records claimed artifact ids from
`workflow_dag.steps[*].produces` when a step exits successfully. Artifact
realization is reported by the legacy artifact manifest files described below.

## Restitution-curves workflow

The `restitutionCurves` entry runs the
`electrophysiologyProtocols/restitutionCurves_s1s2Protocol` tutorial across ionic
models, tissues, and S2 intervals.

Default supported restitution ionic models and tissues:

- `BuenoOrovio`: `mCells`, `endocardialCells`, `epicardialCells`
- `TNNP`: `mCells`, `endocardialCells`, `epicardialCells`
- `Gaur`: `myocyte`
- `Courtemanche`: `myocyte`
- `Stewart`: `myocyte`

Per-case collection writes trace text files and voltage-trace MP4 files under:

```text
<output_dir>/<ionicModel>/<tissue>/<ionicModel>_<tissue>_S1_*_S2_<interval>.txt
<output_dir>/<ionicModel>/<tissue>/<ionicModel>_<tissue>_S1_*_S2_<interval>.mp4
```

## Config override model

`--config` accepts either:

- top-level map keyed by entry name, or
- direct object with `make_spec(...)` keyword args for the selected entry.

Canonical example config:

- `applications/scripts/driverFoam/openfoam_driver/spec_overrides.example.json`

Common override keys across tutorials:

- `case_dir_name`
- `setup`
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
`$ELECTRO_MODEL_COEFFS` resolves to the selected
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
- `electro_property_overrides: { "$ELECTRO_MODEL_COEFFS.activeTensionModel": "GoktepeKuhl" }`
- `electro_property_overrides: { "$ELECTRO_MODEL_COEFFS.ecgDomains.ECG.electrodePositions.V1": "(1 2 3)" }`

Notes:

- Overrides update existing entries only. They do not insert brand new keys.
- For dimensioned scalars, vectors, and tensors, pass the full OpenFOAM literal
  as a string, for example
  `"[0 -3 0 0 0 1 0] 50000"` or
  `"[-1 -3 3 0 0 2 0] (0.133418 0 0 0.0176062 0 0.0176062)"`.
- The catalog lists repository-backed keys. Additional OpenFOAM ODE-solver
  pass-through keys may also be valid; the dictionaries are forwarded to
  `ODESolver::New(...)`.

Generic-case sweeps can also use a `cases` array. Each case may override
`electro_property_overrides`, `physics_property_overrides`, `dimension`,
`parallel`, `touch_case_foam`, and `openfoam_bashrc` on top of the spec defaults.

Defaults live in `plugins/cardiacfoam/defaults/*.py`.

## Runtime artifacts

- `workflow_state.json`: written by strict `step`/`run`; contains the current
  workflow status, current/failed step id, completed steps, attempt counts,
  logs, exit codes, diagnostics, and claimed produced artifact ids.
- `workflow_logs/<step>.attempt<N>.stdout.log`: stdout for a strict step
  attempt.
- `workflow_logs/<step>.attempt<N>.stderr.log`: stderr for a strict step
  attempt.
- `remediation_history.jsonl`: append-only record of applied `--apply`
  overrides, written by `remediation_audit.append_remediation_record` --
  present only if an override was ever applied for this case.

`workflow_state.json` is the current machine-facing state file for strict
autonomous execution. Status vocabulary is:

- `pending`
- `running`
- `completed`
- `failed`
- `skipped`

Each strict step state reports:

- `step_id`
- `status`
- `attempt`
- `command`
- `args`
- `cwd`
- `started_at`
- `finished_at`
- `exit_code`
- `stdout_log`
- `stderr_log`
- `produced_artifacts`
- `diagnostics`

## Setup-folder dependencies

Expected per tutorial setup assets:

- `singleCell/setup/singleCellinteractivePlots.py`
- `manufacturedSolutions/monodomainPseudoECG/setup/post_processing_manufactured.py`
- `NiedererEtAl2011/NiedererEtAl2011verification/setup/postProcessing/{cache_postProcessing.py,line_postProcessing.py,points_postProcessing.py}`
- `restitutionCurves_s1s2Protocol/setup/postProcessing_restCurves.py`

## Architecture tests

Important contract tests:

- `tests/test_tutorial_architecture_contract.py`
- `tests/test_single_cell_contract.py`
- `tests/test_mutators.py`
- `tests/test_strict_planning.py`
- `tests/test_workflow_contract.py`
- `tests/test_workflow_state.py`
- `tests/test_workflow_runner.py`
- `tests/test_cli_step.py`
- `tests/test_agent_readme_contract.py`
- `tests/test_utility_catalog_contract.py`
- `tests/test_rtst_enum_contract.py`
- `tests/test_ionic_catalog_contract.py`
- `tests/test_mesh_geometry.py`
- `tests/test_mesh_geometry_contract.py`

These ensure registry coverage, required `make_spec(...)` keyword contract
consistency, strict RunDocument validation, strict workflow execution behavior,
README legacy-token gates, utility manifest coverage, and catalog drift checks.

Strict dict-key scanner:

```bash
python3 applications/scripts/driverFoam/scripts/scan-dict-keys.py --strict
```

This scanner fails when new uncatalogued C++ dict keys appear, stale catalog
paths remain, or allowlist entries in
`openfoam_driver/plugins/cardiacfoam/dict_key_allowlist.json` become unused.
