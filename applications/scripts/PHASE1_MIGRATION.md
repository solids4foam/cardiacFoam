# Phase 1 Migration (Non-Breaking)

Phase 1 keeps old entrypoints working while defining the new canonical architecture.

## Canonical Entry Points

- `diffusivity/diffusivity.py`
- `scar/scar.py`
- `purkinje/slab/slab.py`
- `cardiac_preproc/conduction_system_generation.py` (pipeline CLI)

## Canonical Engine API

- `engine/src/cardiac_engine/runtime.py` (`CardiacPreprocEngine`)
- `engine/src/cardiac_engine/conduction_runner.py` (pipeline orchestration)
- `engine/src/cardiac_engine/paths.py` (single source of repository/data paths)

## Core Domain Code (kept in cardiac_preproc)

- `cardiac_preproc/src/cardiac_preproc/steps/*`
- `cardiac_preproc/src/cardiac_preproc/io/*`
- `cardiac_preproc/src/cardiac_preproc/lib/*`
- `cardiac_preproc/src/cardiac_preproc/pipeline/*`

## Compatibility Layer Kept in Phase 1

- `cardiac_preproc/scripts/*.py`
- `cardiac_preproc/src/cardiac_preproc/cli/*.py`
- `cardiac_preproc/src/cardiac_preproc/engine.py` (shim)

These files forward execution to canonical entrypoints and exist to avoid breaking existing user scripts.

## Data Policy

All VTK and `inputs/outputs` assets live under:

- `filesOrganize/`

Legacy folder locations keep symlinks for compatibility during Phase 1.

## Phase 2 (Breaking cleanup)

1. Remove compatibility wrappers and shims.
2. Remove deprecated import paths.
3. Keep only canonical `engine` + product CLIs + core domain packages.
