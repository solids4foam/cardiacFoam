# Upper Architecture

## Top-Level Layout

- `engine/`: upper-layer shared runtime facade (`cardiac_engine` package).
- `purkinje/`: Purkinje product surface with three child modules:
  - `density_mapper/`
  - `fractal_3d/`
  - `slab/`
- `diffusivity/`: product-facing diffusivity entrypoint layer.
- `scar/`: product-facing scar entrypoint layer.
- `cardiac_preproc/`: shared core domain/step implementations and pipeline.
- `anatomicalSegmentationModels/`: anatomical tagging dependency layer.
- `PHASE1_MIGRATION.md`: non-breaking migration contract and deprecation scope.

## Layering Model

1. Product Layer (`purkinje/*`)
- User-facing workflows, wrappers, examples.
- Stable product-oriented entrypoints.
- Includes sibling product wrappers: `diffusivity/`, `scar/`.

2. Engine Runtime Layer (`engine/src/cardiac_engine`)
- High-level facade and orchestration API surface.

3. Core Domain Layer (`cardiac_preproc/src/cardiac_preproc`)
- Domain steps: diffusivity, scar, purkinje_slab.
- Compatibility CLIs: `cardiac_preproc.cli.*` (delegate to upper product CLIs).
- VTK I/O and post-processing.
- Pipeline registration and orchestration.

4. Domain-Dependency Layer
- Anatomical segmentation and mesh tagging support.

## API Orientation

The preferred execution API is step-oriented and engine-centered:

- `cardiac_engine.CardiacPreprocEngine`
- `cardiac_preproc.steps.diffusivity.run_diffusivity(...)`
- `cardiac_preproc.steps.scar.run_scar(...)`
- `cardiac_preproc.steps.purkinje_slab.run_purkinje_slab(...)`
- `cardiac_preproc.conduction_system_generation.main()` for chained pipelines.

Product CLIs live in:
- `diffusivity/diffusivity.py`
- `scar/scar.py`
- `purkinje/slab/slab.py`

They delegate into shared `cardiac_preproc.steps.*` engine functions.

## Naming Rules

- Use `cardiac_preproc` for core/reusable implementation modules.
- Use `purkinje_*` naming only inside Purkinje product packages.
- Keep wrapper scripts thin and delegate to engine modules.
