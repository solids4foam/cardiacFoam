# Conduction System Preprocessing

User-facing tools for Purkinje preprocessing (slab + fractal), diffusivity tensors,
and scar tagging.

## Structure

- `conduction_system_generation.py` — main driver (uses `--steps` CLI)
- `scripts/` — individual entrypoints:
  - `purkinje_slab.py`
  - `purkinje_fractal.py`
  - `diffusivity_vtk.py`
  - `scar_vtk.py`
  - `tagEndoEpi_fromUVC.py`
- `configs/` — user defaults and run plan
- `inputs/meshes/` — input meshes
- `src/conduction_preproc/` — implementation code
- `outputs/` — generated results (created at runtime)
- `videoEditorPurkinje/render_purkinje_frames.py` — ParaView pvpython rendering + GIFs

## Quick Start

Run the main driver:

```bash
python3 conduction_system_generation.py --steps diffusivity purkinje_slab
```

Run a single tool:

```bash
python3 scripts/purkinje_slab.py
python3 scripts/diffusivity_vtk.py
python3 scripts/scar_vtk.py
python3 scripts/purkinje_fractal.py
python3 scripts/purkinje_fractal.py --input outputs/ASCIIlegacy_DiffusionTensor_endo_epi_surface.vtk
python3 scripts/tagEndoEpi_fromUVC.py
```
By default the fractal run opens an interactive viewer after VTUs are written.
To skip it:

```bash
python3 scripts/purkinje_fractal.py --no-view-final
```

## Configuration

Use `conduction_system_generation.py --steps ...` to control which steps run.
Each solver reads its default input/output paths from its own config in `configs/`.

If you want to understand the detailed logic for a step, see the corresponding
script in `scripts/` or the implementation in `src/conduction_preproc/`.

## Pipeline Order (Why it matters)

When running multiple steps, follow:

1) `diffusivity`
2) `purkinje_slab`
3) `scar`
4) `convert`

Reason: `purkinje_slab` scales the `Diffusivity` tensor for tagged cells (using
the multiplier in `configs/purkinjeSlab_config.py`), and `scar` can further
scale the tensor inside scar cells using `DIFFUSIVITY_SCAR_MULTIPLIER` from
`configs/scar_config.py`. Running out of order means those tensors won't be
modified as intended.

## Surface Tagging (Endo/Epi)

`scripts/tagEndoEpi_fromUVC.py` tags LV/RV endocardium and epicardium using
`uvc_transmural` and `uvc_intraventricular` point data. The core implementation
lives in `src/conduction_preproc/tagging/tag_endo_epi_surface.py`. The surface
output is triangulated and saved as an unstructured surface; the volume-mapped
output stays as a volume mesh.

Defaults for surface/volume outputs and the optional seed point live in
`configs/tag_endo_epi_surface_config.py`. If you see boundary artifacts, try
setting a `SEED_POINT` (or `--seed-point`) to grow the inside-shared region
from a known location.
