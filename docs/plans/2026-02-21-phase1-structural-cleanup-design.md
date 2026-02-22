# Phase 1: Structural Cleanup — Design

**Scope:** `applications/scripts/conductionSystem_preProcessing` only.
**No changes to:** `anatomicalSegmentationModels/`, `PurkinjeFractal3D/`, any `lib/`, `steps/`, `pipeline/`, or business logic.

## Changes

### 1. Delete dead code
- Remove `branch_original.py` — exact copy of `branch.py` without debug features, unused.
  - Path: `src/conduction_preproc/purkinje_network/Costabal2015_purkinjeNetwork/src/src/fractal_tree/branch_original.py`

### 2. Fix missing `__init__.py`
- Add empty `src/conduction_preproc/tagging/__init__.py`.
- Without it, `from conduction_preproc.tagging import ...` raises `ModuleNotFoundError`.

### 3. Collapse double CLI wrapper layer (Approach A)

Each `scripts/*.py` currently calls an intermediate wrapper inside `src/` which calls the core.
Merge the intermediate wrapper's argparse `main()` into `scripts/*.py`, then delete the intermediate file.

| Keep | Absorbs | Delete |
|------|---------|--------|
| `scripts/diffusivity_vtk.py` | `diffusivity/diffusionTensor_vtk.py` | `diffusivity/diffusionTensor_vtk.py` |
| `scripts/scar_vtk.py` | `scar/scar_vtk.py` | `scar/scar_vtk.py` |
| `scripts/purkinje_slab.py` | `purkinje_network/purkinje_fractal/purkinje_slab.py` | `purkinje_network/purkinje_fractal/purkinje_slab.py` |

Import chain after change:
- `scripts/diffusivity_vtk.py` → `steps.diffusivity.run_diffusivity`
- `scripts/scar_vtk.py` → `steps.scar.run_scar`
- `scripts/purkinje_slab.py` → `purkinje_network.purkinje_fractal.slab.add_purkinje_layer`

## Non-goals
- No field name constants module (deferred to later phase).
- No `src/src/` path flattening (vendored third-party project structure).
- No logic changes anywhere.
