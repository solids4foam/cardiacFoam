# OpenFOAM Driver (Python Specs)

Shared automation engine for OpenFOAM tutorials, with per-tutorial Python specs.

## Install

From repository root:

```bash
python3 -m pip install -e applications/scripts/driverFoam
```

## Run

From repository root:

```bash
# OpenFOAM-style wrapper script
applications/scripts/driverFoam/bin/driverFoam sim --tutorial singleCell
applications/scripts/driverFoam/bin/driverFoam post --tutorial singleCell

# after adding scripts to PATH
export PATH="$PWD/applications/scripts/driverFoam/bin:$PATH"
driverFoam sim --tutorial singleCell

# installed entrypoints (after pip install -e applications/scripts/driverFoam)
driverFoam sim --tutorial singleCell
foamctl sim --tutorial singleCell
foamctl post --tutorial singleCell
foamctl all --tutorial singleCell
foamctl sim --tutorial niederer2012 --dry-run
foamctl sim --tutorial manufacturedFDA --dry-run
foamctl sim --tutorial singleCell --tutorials-root tutorials

# or without installing entrypoint
python applications/scripts/driverFoam/run_openfoam_driver.py sim --tutorial singleCell
python applications/scripts/driverFoam/run_openfoam_driver.py post --tutorial singleCell
python applications/scripts/driverFoam/run_openfoam_driver.py all --tutorial singleCell

# module invocation (after `pip install -e applications/scripts/driverFoam`)
python3 -m openfoam_driver sim --tutorial singleCell

# with JSON overrides
python -m openfoam_driver sim --tutorial singleCell \
  --config applications/scripts/driverFoam/spec_overrides.example.json
```

## Dry-run

```bash
python applications/scripts/driverFoam/run_openfoam_driver.py sim --tutorial singleCell --dry-run
```

## Spec Overrides (JSON)

The `--config` file maps directly to each tutorial `make_spec(...)` parameters.
Units are **ms** for time and **mm** for space.

### Common options

| Key | Type | Description |
| --- | --- | --- |
| `case_dir_name` | `string` | Case folder under `tutorials/` (required in most custom configs). |
| `tutorials_root` | `string` | Root folder containing tutorial cases (default `<repo>/tutorials`). |
| `setup_dir_name` | `string` | Setup folder inside the case. Default: `setup<CaseDirName>`. |
| `output_dir_name` | `string` | Output folder inside the case. Default: `postProcessing`. |
| `run_script_relpath` | `string` | Run script path. Default: `applications/scripts/driverFoam/openfoam_driver/scripts/run_case.sh` (relative to repository root). |
| `postprocess_strict_artifacts` | `bool` | If true, fail post-processing when declared artifacts are missing on disk. Default: `false`. |

### `singleCell` options

| Key | Type | Description |
| --- | --- | --- |
| `ionic_models` | `string[]` | Ionic models to sweep. |
| `ionic_model_tissue_map` | `object` | Map: model -> list of tissues. |
| `stimulus_map` | `object` | Map: model -> stimulus amplitude. |
| `electro_properties_scope` | `string` | Parent sub-dictionary scope in `electroProperties` (default `singleCellElectroCoeffs`). Single-cell stimulus keys are written under `singleCellStimulus` inside this scope. |
| `output_glob` | `string` | Pattern for files collected into output dir (default `*.txt`). |
| `electro_properties_relpath` | `string` | Case-relative path to electro dictionary (default `constant/electroProperties`). |
| `postprocess_script_relpath` | `string` | Setup-relative Python postprocess script. |
| `postprocess_function_name` | `string` | Function name called from postprocess script. |

### `manufacturedFDA` options

| Key | Type | Description |
| --- | --- | --- |
| `number_cells` | `int[]` | Cells per sweep point. |
| `dt_values` | `float[]` | Time-step values (**ms**). |
| `dimensions` | `string[]` | Dimensions to run: `1D`, `2D`, `3D`. |
| `solver_types` | `string[]` | Solvers to run: `explicit`, `implicit`. |
| `piecewise_sweep` | `bool` | If true, pairs `dt_values[i]` with `number_cells[i]`. |
| `electro_properties_scope` | `string` | Sub-dictionary scope in `electroProperties` (default `monoDomainElectroCoeffs`). |
| `control_dict_relpath` | `string` | Case-relative path to `controlDict`. |
| `electro_properties_relpath` | `string` | Case-relative path to electro dictionary (default `constant/electroProperties`). |
| `block_mesh_dict_template` | `string` | Case-relative mesh dict template (`{dimension}` placeholder). |
| `postprocess_script_relpath` | `string` | Setup-relative Python postprocess script. |

### `niederer2012` options

| Key | Type | Description |
| --- | --- | --- |
| `ionic_models` | `string[]` | Ionic models to sweep. |
| `ionic_model_tissue_map` | `object` | Map: model -> list of tissues. |
| `dt_values` | `float[]` | Time-step values (**ms**). |
| `dx_values` | `float[]` | Spatial resolution values (**mm**). |
| `solvers` | `string[]` | Solver list (`explicit`, `implicit`). Default in this spec is `implicit` for stability. |
| `slab_size_mm` | `float[3]` | Fixed slab size `(x, y, z)` in **mm** used to derive hex resolution from `dx`. |
| `end_time_by_dx` | `object` | Map: `dx` -> end time (**s**). |
| `electro_properties_scope` | `string` | Sub-dictionary scope in `electroProperties` (default `monoDomainElectroCoeffs`). |
| `control_dict_relpath` | `string` | Case-relative path to `controlDict`. |
| `block_mesh_dict_relpath` | `string` | Case-relative path to `blockMeshDict`. |
| `electro_properties_relpath` | `string` | Case-relative path to electro dictionary (default `constant/electroProperties`). |
| `points_function_object_name` | `string` | Function object name used for point probes (default `Niedererpoints`). |
| `line_function_object_name` | `string` | Function object name used for line probes (default `Niedererlines`). |
| `sampled_field` | `string` | Sampled field name (default `activationTime`). |
| `sampled_points` | `array` | Point definitions `[(label, x, y, z), ...]` written to points CSV. |
| `line_start` | `float[3]` | Line start coordinate in meters. |
| `line_end` | `float[3]` | Line end coordinate in meters. |
| `line_n_points` | `int` | Number of probe samples along the line. |
| `line_postprocess_relpath` | `string` | Setup-relative line postprocess script. |
| `points_postprocess_relpath` | `string` | Setup-relative points postprocess script. |
| `excel_reference_relpath` | `string` | Setup-relative Excel reference path. |
| `line_postprocess_function_name` | `string` | Function name for line postprocess. |
| `points_postprocess_function_name` | `string` | Function name for points postprocess. |

### `restitutionCurves` options

| Key | Type | Description |
| --- | --- | --- |
| `ionic_models` | `string[]` | Ionic models to sweep. |
| `ionic_model_tissue_map` | `object` | Map: model -> list of tissues. |
| `stimulus_map` | `object` | Map: model -> stimulus amplitude. |
| `s1_interval_ms` | `int` | S1 pacing interval (**ms**). |
| `s2_intervals_ms` | `int[]` | S2 intervals (**ms**). |
| `n_s1` | `int` | Number of S1 beats. |
| `n_s2` | `int` | Number of S2 beats. |
| `electro_properties_scope` | `string` | Parent sub-dictionary scope in `electroProperties` (default `singleCellElectroCoeffs`). Single-cell stimulus keys are written under `singleCellStimulus` inside this scope. |
| `electro_properties_relpath` | `string` | Case-relative path to electro dictionary (default `constant/electroProperties`). |
| `control_dict_relpath` | `string` | Case-relative path to `controlDict`. |
| `postprocess_script_relpath` | `string` | Setup-relative Python postprocess script. |

## Setup Folder Dependencies

The driver run script is centralized at `applications/scripts/driverFoam/openfoam_driver/scripts/run_case.sh`.  
Remaining setup-folder dependencies are postprocess/export assets:

### `singleCell/setupSingleCell`

- `singleCellinteractivePlots.py`

### `manufacturedFDA/setupManufacturedFDA`

- `post_processing_manufactured.py`

### `NiedererEtAl2012/setupNiedererEtAl2012`

- `postProcessing/line_postProcessing.py`
- `postProcessing/points_postProcessing.py`
- `postProcessing/Niederer_graphs_webplotdigitilizer_points_slab/WebPlotDigitilizerdata.xlsx`

### `restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol`

- `postProcessing_restCurves.py`
- `animate_trace.py`

## Current Specs

- `singleCell`
- `niederer2012`
- `manufacturedFDA`
- `restitutionCurves`

## Notes

- Runtime lives in `openfoam_driver/core/runtime/`.
- Specs live in `openfoam_driver/specs/tutorials/`.
- Shared spec helpers live in `openfoam_driver/specs/common.py`.
- Central defaults live in `openfoam_driver/core/defaults/`.
- Example overrides file: `applications/scripts/driverFoam/spec_overrides.example.json`.
- The engine writes `run_manifest.json` after simulation planning/execution.

## Postprocess Contract

All driver-managed postprocess scripts expose:

```python
def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **kwargs)
```

Expected return value:
- `None`, or
- a path-like value, or
- a dict with `path` (and optional `kind`, `format`, `label`), or
- a list of the entries above.

Declared artifact entries are collected into `plots.json`.

## Artifact Manifest (`plots.json`)

The postprocess runner writes `<output_dir>/plots.json` with schema version `1.0`.

Top-level keys:
- `schema_version`
- `tutorial`
- `output_dir`
- `generated_at_utc`
- `artifact_count`
- `artifacts`

Per-artifact keys:
- `path` (relative to output dir when possible)
- `absolute_path`
- `exists`
- `kind` (for example `plot`, `table`)
- `format` (for example `html`, `png`, `csv`)
- `label`
- `task_index`
- `module`
- `function`
