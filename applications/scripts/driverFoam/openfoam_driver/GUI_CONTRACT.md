# GUI Contract

This document defines the intended machine-facing contract between a future GUI
and `openfoam_driver`.

## Entry point

Use:

```bash
python3 -m openfoam_driver describe --entry <name-or-path> [--entry-kind <kind>] [--config overrides.json]
```

The command prints a single JSON object. The GUI should treat that JSON as the
source of truth for:

- entry resolution
- `make_spec(...)` parameters and defaults
- resolved case/setup/output paths
- planned cases for the current configuration
- editable dict-entry metadata

## Top-level JSON shape

The `describe` payload currently contains:

- `requested_entry`
- `requested_tutorial`
- `entry`
- `entry_kind`
- `resolution`
  Values:
  `registered`, `case_folder`, `generic_alias`
- `resolved_name`
- `entry_catalog`
- `workflow`
- `workflow_catalog`
- `is_runnable`
- `registered_tutorials`
- `special_tutorial_aliases`
- `available_tutorials`
- `case_directories`
- `common_override_keys`
- `make_spec`
- `factory_overrides`
- `spec`
- `dict_entries`
- `ionic_model_catalog`
- `active_tension_catalog`
- `gui_schema`
- `launch`
- `config_schema`
- `manifest_schema`

## `make_spec` schema

`make_spec.parameters` is a map keyed by parameter name. Each parameter reports:

- `kind`
- `required`
- optional `annotation`
- optional `default`

The GUI can use this to auto-build entry-level forms.

## `spec` payload

`spec` reports:

- `name`
- `case_root`
- `setup_root`
- `output_dir`
- `metadata`
- `cases`

`cases` contains:

- `count`
- `items`

Each case item contains:

- `case_id`
- `params`

This is intended for previewing the currently planned sweep before execution.

## `dict_entries` payload

`dict_entries.physicsProperties` is a flat list.

`dict_entries.electroProperties` is grouped by category:

- `top_level`
- `common_model_coeffs`
- `ode_solver_passthrough`
- `single_cell_stimulus`
- `monodomain`
- `eikonal_diffusion`
- `ecg`
- `active_tension`

Each dict entry reports:

- `driver_path`
- `description`
- `source_refs`
- `notes`
- `value_kind`
- `ui_control`
- `enum_values`
- `examples`
- `dynamic_path`
- `required`
  Boolean. `true` if the parameter must be supplied for a valid simulation; `false` if it has a safe default or is optional.
- `constraints`
  List of strings. Each string describes a constraint on this parameter's value or relationship with another parameter. Empty list if unconstrained.
- `unit`
  SI unit string for the parameter (e.g. `"S/m"`, `"1/m"`, `"F/m²"`, `"s"`, `"A/m³"`). Empty string for enum, boolean, or dimensionless parameters.
- `typical_value`
  A canonical value from a working tutorial as a string. Empty string for enum or boolean parameters. For `dimensioned_scalar_literal` or `dimensioned_tensor_literal` entries, this is the numeric magnitude only (the full OpenFOAM literal must still include the dimension set when writing the file).

## `ionic_model_catalog`

`ionic_model_catalog` is a static machine-readable catalog of cardiac ionic models
and solver compatibility rules. Agents use this key to plan ionic-model
configuration without requiring an OpenFOAM environment at planning time.

Top-level keys:

- `schema_version`
  String. Currently `"1.0"`. Increment when the catalog schema changes in a
  breaking way.
- `ionic_models`
  Map from model name to `IonicModelEntry`.
- `solver_compatibility`
  List of solver compatibility rules.

### `IonicModelEntry` fields

- `states` — list of state variable names (from `stateVariableNames()`)
- `algebraic` — list of algebraic variable names (from `algebraicVariableNames()`)
- `constants` — list of constant parameter names (from `constantVariableNames()`)
- `recommended_exports` — minimum useful set of variable names for `outputVariables.ionic.export`
- `compatible_tissues` — valid `tissue` values for this model (e.g. `"epicardialCells"`, `"myocyte"`)
- `compatible_solvers` — valid `myocardiumSolver` values (e.g. `"monodomainSolver"`, `"singleCellSolver"`)
- `species` — species this model was derived from (e.g. `"human"`, `"pig"`)
- `cardiac_region` — cardiac region (e.g. `"ventricle"`, `"atrium"`, `"purkinje"`)
- `model_type` — `"phenomenological"`, `"ionic"`, or `"manufactured"`
- `description` — human-readable summary including biological context, typical use case, computational cost signal, and key limitations
- `notes` — additional implementation notes (may be empty)
- `aliases` — list of alternative names and abbreviations users may use when referring to this model (e.g. `"ten Tusscher"`, `"TT06"`, `"TT2006"` for `TNNP`); useful for semantic matching and in-context model selection
- `recommended_ode_step` — typical `initialODEStep` in seconds for this model; dimensionless models use a larger step (e.g. `1e-3`), ionic models typically `1e-5`
- `recommended_stimulus_duration` — typical stimulus duration in seconds; `null` for dimensionless or manufactured models where the unit is not seconds
- `recommended_stimulus_intensity` — typical stimulus current density in A/m³; `null` for dimensionless or manufactured models

### `solver_compatibility` rules

Each rule is a dict with:

- `myocardium_solver` — solver name or `"*"` (wildcard)
- `purkinje_solver` — solver name or `"*"` (wildcard)
- `required_coupler` — required coupler name, or `null` if the combination is invalid
- `valid` — boolean
- `reason` — string explanation (present on invalid rules only)

Valid combinations:

| Myocardium solver | Purkinje solver | Coupler |
|---|---|---|
| `monodomainSolver` | `monodomain1DSolver` | `reactionDiffusionPvjCoupler` |
| `eikonalSolver` | `eikonalSolver` | `eikonalPvjCoupler` |

All other cross-physics combinations are invalid. `bidomainSolver` and
`singleCellSolver` do not support Purkinje network coupling.

**Note:** `solutionAlgorithm` (PDE↔ODE coupling loop) and
`electrophysicsAdvanceScheme` (myocardium↔Purkinje multi-domain coupling) are
**independent** parameters — there is no constraint between them.

## `active_tension_catalog`

`active_tension_catalog` is a static machine-readable catalog of active tension
models. Agents use this key to plan electromechanical coupling and active
tension output without conflating those models with ionic-model selection.

Top-level keys:

- `schema_version`
  String. Currently `"1.0"`. Increment when the catalog schema changes in a
  breaking way.
- `active_tension_models`
  Map from model name to `ActiveTensionModelEntry`.

### `ActiveTensionModelEntry` fields

- `states`, `algebraic`, `constants` — same semantics as ionic models
- `rates` — list of rate variable names (from `rateVariableNames()`)
- `recommended_exports` — minimum useful set of variable names for active tension output
- `description`, `notes` — same as ionic models
- `aliases` — alternative names and abbreviations for this active tension model

## `gui_schema`

`gui_schema` is a machine-readable frontend planning layer produced by
`openfoam_driver.gui_schema.describe_gui_schema()`.

It currently contains:

- `routes`
- `view_models`

### `routes`

Each route reports:

- `id`
- `path`
- `title`
- `purpose`
- `primary_data`
- `actions`

This is not a router implementation. It is the recommended application
structure for the first GUI pass.

Current recommended routes:

- `/entries`
- `/entries/:entryId`
- `/entries/:entryId/config`
- `/entries/:entryId/cases`
- `/runs`
- `/runs/:runId`

### `view_models`

Each view-model spec reports:

- `id`
- `description`
- `fields`
- `source`

Current view-models cover:

- entry catalog
- workflow family metadata
- entry overview
- spec parameter form fields
- dict-entry editor fields
- planned cases
- run-launch payload
- run summary
- per-case run result

## `launch`

`launch` is a map keyed by driver action:

- `sim`
- `post`
- `all`

Each launch item reports:

- `action`
- `entry`
- `entry_kind`
- `tutorial`
- `resolved_name`
- `resolution`
- `command`
- `command_display`
- `manifest_path`
- `plots_manifest_path`
- `case_root`
- `setup_root`
- `output_dir`
- `dry_run`
- `continue_on_error`
- `postprocess_available`

This lets a local app present the exact driver command it is about to run
without re-implementing registry or path resolution logic.

## `value_kind`

Current value kinds:

- `enum`
  Fixed word choices.
- `boolean`
  OpenFOAM boolean-like toggles handled as a GUI checkbox.
- `integer`
  Whole-number numeric fields.
- `scalar`
  Floating-point numeric fields.
- `word_list`
  Token lists such as `(Vm Iion)`.
- `scalar_list`
  Numeric lists.
- `vector3`
  A three-component vector like `(x y z)`.
- `vector3_list`
  A list of vectors.
- `dimensioned_scalar_literal`
  Requires an OpenFOAM literal string including dimensions when relevant.
- `dimensioned_tensor_literal`
  Requires a full OpenFOAM tensor literal string.
- `openfoam_literal`
  Fallback type when a more specific GUI type is not yet assigned.

## `ui_control`

These are recommendations, not hard requirements:

- `select`
- `checkbox`
- `number`
- `text`
- `textarea`
- `token_list`
- `vector3`

The GUI may use a richer widget set internally, but it should preserve the
same value semantics.

## Dynamic paths

When `dynamic_path` is `true`, the path contains a placeholder segment rather
than a fixed dictionary key.

Current example:

- `$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.electrodePositions.<name>`

The GUI should not assume that `<name>` can be created arbitrarily. The current
driver mutator updates existing entries only; it does not create new keys.

## `config_schema`

`config_schema` documents the JSON format accepted by `--config`. It contains:

- `description` — prose explanation of the format
- `top_level_shapes` — the two valid shapes (`flat` and `wrapped`)
- `section_fields` — each key accepted inside a tutorial section
- `worked_example` — a ready-to-use minimal config for this tutorial

### Config file shape

**Flat** — one tutorial, parameters at the top level:

```json
{
  "ionic_models": ["TNNP"],
  "electro_property_overrides": {
    "$ELECTRO_MODEL_COEFFS.chi": "140000"
  }
}
```

**Wrapped** — multi-tutorial file, each tutorial in its own named section
(key matching is case-insensitive):

```json
{
  "singleCell": {
    "ionic_models": ["TNNP"],
    "electro_property_overrides": {
      "$ELECTRO_MODEL_COEFFS.chi": "140000"
    }
  }
}
```

### `section_fields`

| Key | Type | Description |
|-----|------|-------------|
| *(spec parameters)* | varies | Parameters accepted by `make_spec()` for this tutorial. Listed in `config_schema.section_fields.spec_parameters.available_keys`. |
| `electro_property_overrides` | mapping or list | Overrides for `constant/electroProperties`. |
| `physics_property_overrides` | mapping or list | Overrides for `constant/physicsProperties`. |
| `case_dir_name` | string | Override the case directory name. |
| `setup_dir_name` | string | Override the setup directory (default: `<case_dir_name>_setup`). |
| `output_dir_name` | string | Override where the manifest and outputs are written. |
| `run_script_relpath` | string | Path to an alternative run script. |

### `electro_property_overrides` format

**Shorthand (recommended)** — a mapping from `driver_path` to value string.
Use exactly the `driver_path` values from `dict_entries.electroProperties`.
The `$ELECTRO_MODEL_COEFFS` token resolves automatically to the correct solver
coeffs dict (e.g. `monodomainSolverCoeffs`):

```json
{
  "$ELECTRO_MODEL_COEFFS.chi": "140000",
  "$ELECTRO_MODEL_COEFFS.cm": "0.01",
  "$ELECTRO_MODEL_COEFFS.initialODEStep": "1e-5",
  "$ELECTRO_MODEL_COEFFS.ionicModel": "TNNP"
}
```

**Explicit** — a list of `{key, scope, value}` objects, for cases where you need
to target a sub-dictionary by its literal name:

```json
[
  {"key": "chi",       "scope": ["monodomainSolverCoeffs"], "value": "140000"},
  {"key": "ionicModel","scope": ["monodomainSolverCoeffs"], "value": "TNNP"}
]
```

---

## `manifest_schema`

`manifest_schema` documents `run_manifest.json` — the run-state source of truth.

### File location

```
<output_dir>/run_manifest.json
```

The exact path is also available as `launch.<action>.manifest_path`.

### Polling contract

Poll every **15–30 seconds**. Stop when `status` is in the terminal state list.
The file is written atomically; reading it at any time is safe.

### Companion file: `action_events.jsonl`

Written to the same directory. Append-only JSONL; one event per line.
Events: `sim_started`, `case_started`, `case_finished`, `sim_finished`,
`postprocess_started`, `postprocess_finished`, `all_started`, `all_finished`.
Each line has `event`, `timestamp_utc`, and `run_id` fields plus event-specific payload.

### Top-level fields

| Field | Type | Description |
|-------|------|-------------|
| `schema_version` | string | Currently `"2.1"` |
| `run_id` | string | Unique run identifier |
| `requested_action` | string | `"sim"`, `"post"`, or `"all"` |
| `entry` | string | Selected entry name |
| `entry_kind` | string\|null | Entry classification |
| `entry_path` | string\|null | Resolved relative path under `tutorials/` |
| `source_type` | string\|null | Resolution source such as `spec_factory` or `workflow_reference_case` |
| `workflow_family` | string\|null | Workflow family name when applicable |
| `case_root` | string | Absolute path |
| `setup_root` | string | Absolute path |
| `output_dir` | string | Absolute path |
| `dry_run` | boolean | |
| `continue_on_error` | boolean | |
| `status` | string | See status values below |
| `postprocess_status` | string | See postprocess status values below |
| `current_case_id` | string\|null | Case currently executing; null between cases |
| `started_at_utc` | string\|null | ISO 8601 |
| `updated_at_utc` | string | ISO 8601 — timestamp of last write |
| `finished_at_utc` | string\|null | ISO 8601; null until terminal state |
| `total_cases` | integer | |
| `planned_cases` | integer | Non-zero only in dry-run mode |
| `completed_cases` | integer | Cases with status `"ok"` |
| `failed_cases` | integer | Cases with status `"failed"` |
| `error` | string\|null | Top-level error message if run failed early |
| `plots_manifest_path` | string\|null | Path to `plots.json` if postprocess produced plots |
| `human_report_path` | string | Path to `run_report.md` |
| `results` | array | Per-case results — see below |

### `status` values

| Value | Terminal? | Meaning |
|-------|-----------|---------|
| `running` | no | Simulation in progress |
| `completed` | **yes** | All cases finished successfully |
| `completed_with_failures` | **yes** | All cases ran; at least one failed |
| `failed` | **yes** | A case failed and `continue_on_error=false` |
| `postprocessing` | no | Simulations done; postprocessing running |
| `postprocess_failed` | **yes** | Postprocessing raised an exception |
| `planned` | **yes** | Dry-run complete; no simulation was run |

### `postprocess_status` values

`not_started` → `running` → `completed` or `failed`.
`skipped` is set in dry-run mode.

### Per-case result fields (`results` array)

| Field | Type | Description |
|-------|------|-------------|
| `case_id` | string | Unique case identifier |
| `status` | string | `"ok"`, `"failed"`, or `"planned"` |
| `duration_s` | float | Wall-clock seconds |
| `params` | object | Parameter values for this case |
| `error` | string\|null | Exception message if `status="failed"` |
| `index` | integer | 1-based position in the case list |
| `total_cases` | integer | |
| `started_at_utc` | string\|null | ISO 8601 |
| `finished_at_utc` | string\|null | ISO 8601 |

---

## Current limitations

- The driver override layer updates existing dict entries only.
- Some OpenFOAM values are intentionally represented as literal strings because
  preserving the exact OpenFOAM syntax matters.
- ODE solver pass-through options are only partially typed. The repository
  forwards the dictionary to `ODESolver::New(...)`, so additional solver-specific
  keys may exist outside the current catalog.
- `run_manifest.json` is the run-state source of truth, but log streaming is
  still stdout/stderr based rather than a structured event channel.

## Recommended GUI strategy

1. Call `describe`.
2. Build entry-level controls from `make_spec.parameters`.
3. Build dict editors from `dict_entries`.
   - Use `required: true` entries to mark mandatory fields.
   - Display `constraints` strings as inline validation hints.
4. Read `ionic_model_catalog` to plan ionic model configuration:
   a. Pick an `ionicModel` → look up `ionic_models[name].recommended_exports`.
   b. Set `outputVariables.ionic.export` to the `recommended_exports` value.
   c. Verify the chosen `tissue` is in `compatible_tissues`.
   d. Verify the chosen solver is in `compatible_solvers`.
   e. Check `solver_compatibility` before combining a myocardium solver with a
      Purkinje solver.
5. Read `active_tension_catalog` to plan active tension configuration:
   a. Pick an `activeTensionModel` only when electromechanical coupling is intended.
   b. Set `outputVariables.activeTension.export` from
      `active_tension_models[name].recommended_exports`.
6. Use `gui_schema.routes` and `gui_schema.view_models` as the initial app structure.
7. Preview `spec.cases`.
8. Build the `--config` JSON using `config_schema` as the template.
   - Start from `config_schema.worked_example.json`.
   - Set spec parameters from `config_schema.section_fields.spec_parameters.available_keys`.
   - Set dict entries via `electro_property_overrides` using driver_path keys from `dict_entries`.
9. Launch through `launch.<action>.command` or an equivalent backend wrapper.
10. Poll `launch.<action>.manifest_path` every 15–30 s.
   - Stop when `manifest_schema.terminal_states` contains the current `status`.
   - Read `results[].status` to identify which cases succeeded or failed.
   - `action_events.jsonl` (same directory) provides a fine-grained event stream.
11. Persist overrides as JSON compatible with `--config`.
