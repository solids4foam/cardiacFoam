# First-Draft Refactor Architecture

Scope: `cardiac_preproc`, `purkinje/density_mapper`, `purkinje/fractal_3d`  
Constraint: keep fractal core where it is (no relocation yet).

## Core 1: VTK Handler Layer

Responsibilities:
- read/write mesh files
- enforce ASCII saves
- post-process VTK output (array conversion, blank-line cleanup)
- field validation and introspection

Current placement:
- `cardiac_preproc/src/cardiac_preproc/io/`
  - `vtk_mesh.py` for shared read/write
  - `postprocess.py` for cleanup
  - `field_checks.py` for required field validation
- `purkinje/density_mapper/src/purkinje_density/mesh_context.py`
  - mesh + segment context assembly
  - shared ASCII write helper

## Core 2: Domain Logic Layer

Responsibilities:
- independent business logic per domain:
  - diffusivity
  - scar
  - purkinje slab
- pure step execution with explicit options and deterministic outputs

Current placement:
- `cardiac_preproc/src/cardiac_preproc/steps/`
  - `diffusivity.py`
  - `scar.py`
  - `purkinje_slab.py`

Pipeline registration:
- `cardiac_preproc/src/cardiac_preproc/pipeline/__init__.py`

## Core 3: Purkinje Generation Strategies

Two tracks:
1. Slab-based tagging (in preprocessing step system).
2. Fractal generation (active development track).

Current boundary:
- Slab stays integrated in `cardiac_preproc` steps.
- Fractal remains in:
  - `purkinje/fractal_3d/Costabal2015_purkinjeNetwork/`

No fractal-core move in this draft.

## CLI Layer

Scripts should stay thin wrappers:
- parse args
- build options
- call step/domain modules

Current state:
- `diffusivity/diffusivity.py`, `scar/scar.py`, and `purkinje/slab/slab.py`
  follow this pattern.

## Next Refactor Iteration (Suggested)

1. Extract a shared VTK domain package used by both `cardiac_preproc` and `purkinje_density`.
2. Move endocardial area computation into a dedicated reusable module (not inside UI/workbench files).
3. Add integration tests for full step chains:
   - `diffusivity -> purkinje_slab -> scar -> convert`
4. Introduce versioned data contracts for generated fields (`Diffusivity`, `Scar`, `purkinjeLayer`, area tags).
