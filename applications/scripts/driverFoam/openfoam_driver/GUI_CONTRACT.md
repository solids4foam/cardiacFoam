# GUI Contract

This document defines the intended machine-facing contract between a future GUI
and `openfoam_driver`.

## Entry point

Use:

```bash
python3 -m openfoam_driver describe --tutorial <name> [--config overrides.json]
```

The command prints a single JSON object. The GUI should treat that JSON as the
source of truth for:

- tutorial resolution
- `make_spec(...)` parameters and defaults
- resolved case/setup/output paths
- planned cases for the current configuration
- editable dict-entry metadata

## Top-level JSON shape

The `describe` payload currently contains:

- `requested_tutorial`
- `resolution`
  Values:
  `registered`, `case_folder`, `generic_alias`
- `resolved_name`
- `registered_tutorials`
- `special_tutorial_aliases`
- `available_tutorials`
- `case_directories`
- `common_override_keys`
- `make_spec`
- `factory_overrides`
- `spec`
- `dict_entries`
- `gui_schema`
- `launch`

## `make_spec` schema

`make_spec.parameters` is a map keyed by parameter name. Each parameter reports:

- `kind`
- `required`
- optional `annotation`
- optional `default`

The GUI can use this to auto-build tutorial-level forms.

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

- `/tutorials`
- `/tutorials/:tutorialId`
- `/tutorials/:tutorialId/config`
- `/tutorials/:tutorialId/cases`
- `/runs`
- `/runs/:runId`

### `view_models`

Each view-model spec reports:

- `id`
- `description`
- `fields`
- `source`

Current view-models cover:

- tutorial catalog
- tutorial overview
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

- `$ELECTRO_MODEL_COEFFS.ECG.<ecgModel>Coeffs.electrodes.<name>`

The GUI should not assume that `<name>` can be created arbitrarily. The current
driver mutator updates existing entries only; it does not create new keys.

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
2. Build tutorial-level controls from `make_spec.parameters`.
3. Build dict editors from `dict_entries`.
4. Use `gui_schema.routes` and `gui_schema.view_models` as the initial app structure.
5. Preview `spec.cases`.
6. Launch through `launch.<action>.command` or an equivalent backend wrapper.
7. Poll `launch.<action>.manifest_path` for run-state updates.
8. Persist overrides as JSON compatible with `--config`.
