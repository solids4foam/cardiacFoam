# Phase 1 Migration (Non-Breaking)

Phase 1 keeps old entrypoints working while defining the new canonical architecture.

## Canonical Entry Points

- `diffusivity/diffusivity.py`
- `scar/scar.py`
- `purkinje/slab/slab.py`
- `tagging/tagEndoEpi_fromUVC.py`
- `conversion/ASCIIlegacyToVtkUnstructured.py`
- `conduction_system_generation.py` (pipeline CLI)

## Canonical Engine API

- `engine/src/cardiac_engine/runtime.py` (`CardiacPreprocEngine`)
- `engine/src/cardiac_engine/conduction_runner.py` (pipeline orchestration)
- `engine/src/cardiac_engine/paths.py` (single source of repository/data paths)

## Core Domain Code (kept in cardiac_preproc)

- `cardiac_preproc/src/cardiac_preproc/steps/*`
- `cardiac_preproc/src/cardiac_preproc/io/*`
- `cardiac_preproc/src/cardiac_preproc/lib/*`
- `cardiac_preproc/src/cardiac_preproc/pipeline/*`

## Compatibility Layer Status

- `cardiac_preproc/scripts/{diffusivity_vtk.py,scar_vtk.py,purkinje_slab.py,purkinje_fractal.py}`: removed in cleanup.
- `cardiac_preproc/scripts/tagEndoEpi_fromUVC.py`: moved to `tagging/tagEndoEpi_fromUVC.py`.
- `cardiac_preproc/configs/*`: moved to top-level product-owned config folders.
- `cardiac_preproc/src/cardiac_preproc/cli/*`: removed in cleanup.
- `cardiac_preproc/src/cardiac_preproc/engine.py`: removed in cleanup.

Canonical entrypoints should be used directly.

## Data Policy

All VTK and `inputs/outputs` assets live under:

- `filesOrganize/`

Legacy folder locations keep symlinks for compatibility during Phase 1.

## Phase 2 (Breaking cleanup)

1. Remove compatibility wrappers and shims.
2. Remove deprecated import paths.
3. Keep only canonical `engine` + product CLIs + core domain packages.
