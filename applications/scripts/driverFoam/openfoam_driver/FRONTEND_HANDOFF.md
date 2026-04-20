# Frontend Handoff For `openfoam_driver`

This document summarizes the current `openfoam_driver` implementation, the
backend-oriented GUI preparation work completed so far, the main architectural
strengths and gaps, and the recommended next step for a frontend or local-app
implementation.

The intent is to give a frontend-focused agent a stable handoff point without
requiring it to rediscover the current backend design.

## Scope

The work completed so far is backend and contract preparation, not frontend
implementation.

The driver is now prepared to support a local GUI through:

- entry discovery and execution
- generic dictionary override support
- source-backed editable dict-entry metadata
- machine-readable entry introspection
- machine-readable launch planning
- machine-readable run-state manifests

The frontend itself has not been implemented yet.

## What Was Added Or Changed

### 1. Entry execution

Previously the driver mainly supported a fixed set of curated tutorial specs.
It now also supports:

- registered tutorials
- `genericCase` and `randomCase` aliases
- workflow templates and workflow reference cases
- any existing case folder discoverable under `tutorials/`

Main files:

- `core/runtime/registry.py`
- `specs/tutorials/generic_case.py`

This means the driver no longer depends only on a hardcoded tutorial registry.
It can resolve a direct case folder or workflow entry and run it through a generic spec when appropriate.

### 2. Generic dict override model

The driver now updates both `electroProperties` and `physicsProperties` through
one generic override path model instead of tutorial-specific setter functions.

Main file:

- `specs/common.py`

Key points:

- supports dotted-path overrides
- supports `$ELECTRO_MODEL_COEFFS` as a dynamic scope token
- works for both curated specs and generic cases
- supports nested paths such as active tension and ECG blocks

Important limitation:

- overrides update existing entries only
- they do not create new dictionary keys automatically

### 3. Curated specs migrated to the generic override path

The curated specs were refactored so they use the same dict override mechanism
as generic cases.

Main files:

- `specs/tutorials/single_cell.py`
- `specs/tutorials/niederer_2012.py`
- `specs/tutorials/manufactured_fda.py`
- `specs/tutorials/restitution_curves.py`

This matters because the frontend can now treat curated and generic runs
through the same dict-editing model.

### 4. Source-backed dict-entry catalog

A catalog of repository-known editable entries was added for:

- `physicsProperties`
- `electroProperties`

Main file:

- `dict_entries.py`

This catalog was assembled from the current cardiacFoamv2 `src/` code and is
intended to drive UI field generation and validation hints.

The catalog currently includes:

- `physicsProperties.type`
- `electroProperties.myocardiumSolver`
- common `<solver>Coeffs` entries
- single-cell stimulus entries
- monodomain external-stimulus entries
- eikonal-diffusion entries
- ECG entries
- active-tension entries
- common ODE-solver pass-through entries

Each entry now exposes frontend-relevant metadata:

- `driver_path`
- `description`
- `source_refs`
- `notes`
- `value_kind`
- `ui_control`
- `enum_values`
- `examples`
- `dynamic_path`

This is the current backend source of truth for dict-editing UI hints.

### 5. Introspection API for GUI consumption

A machine-readable introspection layer was added.

Main file:

- `introspection.py`

It is exposed through:

```bash
python3 -m openfoam_driver describe --entry <name-or-path> [--entry-kind <kind>] [--config overrides.json]
```

The `describe` payload currently includes:

- entry resolution details
- entry catalog and workflow catalog data
- registered tutorials
- discovered case directories
- `make_spec(...)` parameter schema and defaults
- resolved case/setup/output paths
- planned cases
- grouped dict-entry metadata
- GUI route and view-model suggestions
- launch metadata for `sim`, `post`, and `all`

This is the main API a frontend should use to build its initial screens, including the first-step choice between a tutorial/example and a real workflow case.

### 6. Launch planning layer

A launch description layer was added.

Main file:

- `launch.py`

It provides a machine-readable description of how the driver would be launched
for each action:

- `sim`
- `post`
- `all`

Each launch item provides:

- exact command
- command display string
- expected `run_manifest.json` path
- expected `plots.json` path when relevant
- resolved case/setup/output paths
- whether post-processing is available

This allows a local app to launch the backend without re-implementing registry
and path-resolution logic.

### 7. Richer run manifests

The execution engine now writes richer run manifests intended for local-app
polling.

Main file:

- `core/runtime/engine.py`

The manifest now includes:

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
- `results`

Each per-case result currently includes:

- `case_id`
- `status`
- `duration_s`
- `params`
- `error`
- `index`
- `total_cases`
- `started_at_utc`
- `finished_at_utc`

This is now the main run-state contract a frontend should poll.

### 8. CLI alignment

The CLI now correctly forwards the requested action to the engine so manifests
report the true action.

Main file:

- `cli.py`

This matters for GUI correctness because `post`, `sim`, and `all` are now
distinguishable in the manifest.

### 9. GUI route/view-model planning layer

A lightweight GUI planning schema was added.

Main file:

- `gui_schema.py`

This is not frontend code. It is a backend-authored recommendation for the
first screen structure and the initial data/view-model boundaries.

It currently suggests routes for:

- tutorial catalog
- tutorial overview
- tutorial configuration
- planned cases
- runs index
- run detail

## Current Code Architecture

The current architecture is in a good state for a frontend handoff because the
driver has a clearer separation between execution, mutation, discovery, and UI
contract generation.

### A. Runtime layer

Files:

- `core/runtime/models.py`
- `core/runtime/registry.py`
- `core/runtime/engine.py`

Responsibilities:

- `models.py`
  Defines `TutorialSpec` and `CaseConfig`, which are the core execution
  contracts.
- `registry.py`
  Resolves tutorial names, aliases, and direct case-folder matches.
- `engine.py`
  Executes simulations and post-processing and writes run manifests.

Assessment:

- This is the correct abstraction boundary for execution.
- `TutorialSpec` remains the main unit of orchestration.
- The frontend should not need to interact with this layer directly.

### B. Spec layer

Files:

- `specs/tutorials/*.py`
- `specs/common.py`

Responsibilities:

- each tutorial builds a `TutorialSpec`
- `common.py` now contains the shared dict override logic and path helpers
- curated specs and generic specs both use the same mutator model

Assessment:

- This layer is much cleaner than before because dict edits are no longer
  hardcoded in tutorial-specific helper APIs.
- The generic spec is especially important for future GUI flexibility because
  it allows running arbitrary case folders.

### C. Introspection and GUI-contract layer

Files:

- `dict_entries.py`
- `introspection.py`
- `launch.py`
- `gui_schema.py`
- `GUI_CONTRACT.md`

Responsibilities:

- expose editable dict metadata
- expose tutorial parameter schema
- expose launch plans
- expose route/view-model planning

Assessment:

- This is the key GUI-preparation layer.
- It converts internal Python objects into JSON-friendly, machine-readable
  payloads.
- It gives the frontend a stable read model without forcing it to understand
  spec internals.

### D. CLI layer

Files:

- `cli.py`
- `__main__.py`

Responsibilities:

- parse CLI actions
- load config overrides
- call `describe`
- call the engine for `sim`, `post`, and `all`

Assessment:

- This is still thin enough, which is good.
- For a local app, this can either remain the launch target or be wrapped by a
  thin desktop-backend service.

## Why This Is Good GUI Preparation

The backend work done so far addresses the main failure mode of many GUI
projects: a frontend starting before the backend has a stable machine-facing
contract.

The preparation is useful for the GUI because it now answers these questions in
a stable way:

- What tutorials or case folders are runnable?
- What inputs are configurable at the tutorial level?
- What dict entries are editable?
- What kind of values do those entries expect?
- What cases will the current configuration expand into?
- What command will actually run?
- Where should the frontend watch for run progress?
- What artifacts should be visible after completion?

That is enough to start a first useful GUI without inventing a second backend
contract.

## Architectural Strengths

### 1. Good separation between execution and introspection

The runtime engine is not coupled to the GUI metadata generation. This is good
because it keeps the execution path simple while still allowing the frontend to
read a structured view of the system.

### 2. One dict override model for both curated and generic cases

This is one of the most important improvements. It reduces special cases and
means the frontend can treat most editing as one generic problem.

### 3. Direct case-folder support

This makes the driver more future-proof. The frontend can expose existing
tutorial case folders even before a dedicated curated spec exists.

### 4. Source-backed editable entry catalog

The dict catalog is not arbitrary UI metadata. It is based on the repository
source, which makes it much more defensible as a UI contract.

### 5. Run manifests are now useful for a local app

Without the richer manifest, the frontend would have needed to scrape stdout or
reverse-engineer state. The manifest is now enough for simple polling-based run
monitoring.

## Current Gaps And Possible Improvements

The backend is ready enough for a first frontend pass, but it is not finished.
These are the main areas that could still be improved.

### 1. Structured log streaming

Current status:

- run state is file-backed through `run_manifest.json`
- logs are still stdout/stderr based

Impact:

- a frontend can show progress and final state
- it cannot yet consume structured incremental logs without wrapping process IO

Suggested improvement:

- add a line-oriented log file or event file per run
- optionally add structured events such as `case_started`, `case_finished`,
  `postprocess_started`, `postprocess_failed`

### 2. Validation could be stricter

Current status:

- the dict override path model is generic
- the dict-entry catalog contains typing hints
- there is still limited hard validation before execution

Impact:

- the frontend can guide the user, but the backend may still only fail once a
  run is attempted

Suggested improvement:

- add preflight validation against known dict paths
- validate required config combinations earlier
- report field-level errors in a machine-readable format

### 3. Run discovery is still implicit

Current status:

- the frontend can poll a manifest path once it already launched a run
- there is no dedicated backend API for enumerating prior runs

Impact:

- a local app can monitor active runs it launched
- showing a robust historical run list would need additional directory scanning

Suggested improvement:

- add a `list_runs(...)` helper that discovers manifests under known output
  directories

### 4. Dict catalog is intentionally incomplete for some pass-through solver keys

Current status:

- common ODE-solver pass-through keys are included
- the backend intentionally documents that additional keys may exist

Impact:

- frontend coverage is strong for common paths
- it is not a mathematically complete schema for every possible OpenFOAM
  pass-through key

Suggested improvement:

- expand source analysis or add an expert-mode raw override editor in the GUI

### 5. The launch path is still CLI-centric

Current status:

- `launch.py` gives the exact CLI command
- there is not yet a dedicated long-lived local service API

Impact:

- this is fine for a first local app
- a more mature desktop app might want a backend process API rather than only
  shelling out to `python -m openfoam_driver`

Suggested improvement:

- if needed later, add a thin backend service layer instead of changing the
  execution engine itself

## Recommended Next Step For The Frontend

The next step is frontend implementation, but it should stay focused on
visibility and control rather than styling.

The first frontend should be a thin client over the current backend contract.

### Recommended visible screens

1. Tutorial catalog

Show:

- registered tutorials
- discovered case folders
- generic-case entrypoint

Data source:

- `describe_tutorial(...).registered_tutorials`
- `available_tutorials`
- `case_directories`

2. Tutorial overview

Show:

- requested tutorial
- resolution mode
- resolved name
- `case_root`
- `setup_root`
- `output_dir`
- metadata summary

Data source:

- `describe_tutorial(...)`

3. Configuration screen

Show:

- `make_spec.parameters`
- common override keys
- grouped dict-entry editors

Data source:

- `make_spec.parameters`
- `common_override_keys`
- `dict_entries`

Notes:

- render fields from `value_kind` and `ui_control`
- keep a raw JSON override view as an expert fallback

4. Planned cases screen

Show:

- case count
- expanded case list
- per-case parameter payload

Data source:

- `spec.cases`

5. Run launch controls

Show:

- `sim`
- `post`
- `all`
- dry-run toggle where valid
- continue-on-error toggle where valid
- exact launch command

Data source:

- `launch`

6. Run monitor

Show:

- current status
- current case id
- counts for planned/completed/failed cases
- per-case result list
- error state

Data source:

- `run_manifest.json`

7. Results/artifacts view

Show:

- `output_dir`
- `plots_manifest_path`
- success/failure summary

Data source:

- `run_manifest.json`
- `plots.json` when present

## Recommended Frontend Data Flow

For a first implementation:

1. On entry selection, call `describe`.
2. Build forms from `make_spec.parameters` and `dict_entries`.
3. Persist current UI state as JSON compatible with `--config`.
4. Use `describe` again after config changes to refresh:
   - resolved paths
   - planned cases
   - launch metadata
5. Launch the driver using the returned launch command or an equivalent local
   backend wrapper.
6. Poll the returned manifest path.
7. Show artifacts once the run reaches a terminal state.

This is enough for a useful local app without inventing a second orchestration
layer.

## Recommended Frontend Boundaries

The frontend should avoid re-implementing backend logic for:

- tutorial resolution
- case discovery
- path derivation
- launch command construction
- dict path semantics
- run progress inference

Instead, it should treat the current backend as the source of truth and focus
on:

- visibility
- form generation
- input persistence
- process launch
- manifest polling
- result presentation

## Suggested First Frontend Milestone

The first milestone should be intentionally modest:

- select a tutorial or case folder
- edit top-level parameters
- edit dict overrides
- preview cases
- launch `sim`
- poll `run_manifest.json`
- show results and errors

That would already provide real value and would validate whether the current
backend contract is sufficient before more advanced UI work.

## Tests And Stability

The backend preparation work is covered by driver tests, including:

- introspection tests
- dict-entry catalog tests
- GUI-schema tests
- launch description tests
- run manifest tests

At the time of writing, the driver-side test suite passes with:

```bash
PYTHONPATH=applications/scripts/driverFoam python3 -m unittest discover applications/scripts/driverFoam/openfoam_driver/tests
```

## Final Assessment

The backend is now in a good state for frontend work.

The important point is that the next step should not be another backend
refactor unless the frontend immediately exposes a missing contract. The best
next move is to build a thin local UI against the existing `describe`,
`launch`, and `run_manifest.json` surfaces, then refine the backend only where
the first real frontend pass proves it necessary.
