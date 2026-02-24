# Upper Architecture

## Top-Level Layout

- `purkinje/`: Purkinje product surface with three child modules:
  - `density_mapper/`
  - `fractal_3d/`
  - `slab/`
- `diffusivity/`: product-facing diffusivity entrypoint layer.
- `scar/`: product-facing scar entrypoint layer.
- `tagging/`: product-facing endo/epi tagging entrypoint + config.
- `conversion/`: product-facing VTK conversion entrypoint + config.
- `cardiac_core/`: shared core domain/step implementations and pipeline.
- `anatomicalSegmentationModels/`: anatomical tagging dependency layer.
- `PHASE1_MIGRATION.md`: non-breaking migration contract and deprecation scope.

## Layering Model

1. Product Layer (`purkinje/*`)
- User-facing workflows, wrappers, examples.
- Stable product-oriented entrypoints.
- Includes sibling product wrappers: `diffusivity/`, `scar/`.

2. Engine Runtime Layer (`cardiac_core/engine`)
- High-level facade and orchestration API surface.

3. Core Domain Layer (`cardiac_core`)
- Domain steps: diffusivity, scar, purkinje_slab.
- VTK I/O and post-processing.
- Pipeline registration and orchestration.

4. Domain-Dependency Layer
- Anatomical segmentation and mesh tagging support.

## API Orientation

The preferred execution API is step-oriented and engine-centered:

- `cardiac_core.engine.CardiacPreprocEngine`
- `cardiac_core.steps.diffusivity.run_diffusivity(...)`
- `cardiac_core.steps.scar.run_scar(...)`
- `cardiac_core.steps.purkinje_slab.run_purkinje_slab(...)`
- `conduction_system_generation.py` as top-level chained pipeline API.

Product CLIs live in:
- `diffusivity/diffusivity.py`
- `scar/scar.py`
- `purkinje/slab/slab.py`
- `tagging/tagging.py`
- `conversion/conversion.py`

They delegate into shared `cardiac_core.steps.*` engine functions.

## Naming Rules

- Use `cardiac_core` for core/reusable implementation modules.
- Use `purkinje_*` naming only inside Purkinje product packages.
- Keep wrapper scripts thin and delegate to engine modules.
