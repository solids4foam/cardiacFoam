# Cardiac Preprocessing Engine

Core engine for preprocessing workflows (diffusivity, scar, tagging, and
Purkinje slab/fractal orchestration).

Upper runtime facade lives in:

- `../engine/src/cardiac_engine/`

## Structure

- `conduction_system_generation.py` — main driver (uses `--steps` CLI)
- `scripts/` — individual entrypoints:
  - `purkinje_slab.py`
  - `purkinje_fractal.py`
  - `diffusivity_vtk.py`
  - `scar_vtk.py`
  - `tagEndoEpi_fromUVC.py`
- `configs/` — compatibility configs for legacy scripts
- `../filesOrganize/` — centralized data assets (inputs/outputs/VTK files)
- `src/cardiac_preproc/` — implementation code
  - `cli/` — core CLI modules (single source of argument parsing)
  - `pipeline/` — step registry/context/runner
  - `steps/` — internal step implementations (`diffusivity`, `scar`)
  - `io/` — shared field checks + VTK post-processing
- `../filesOrganize/cardiac_preproc/outputs/` — generated results
- `videoEditorPurkinje/render_purkinje_frames.py` — ParaView pvpython rendering + GIFs

Product-facing Purkinje wrappers live in:

- `../purkinje/slab/`
- `../purkinje/fractal_3d/`
- `../purkinje/density_mapper/`

Sibling wrappers for non-Purkinje modules:

- `../diffusivity/`
- `../scar/`

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
Each solver reads defaults from colocated product configs:
- `../diffusivity/Diffusivity_config.py`
- `../scar/scar_config.py`
- `../purkinje/slab/purkinjeSlab_config.py`

If you want to understand the detailed logic for a step, see the corresponding
script in `scripts/` or the implementation in `src/cardiac_preproc/`.

## Pipeline Order (Why it matters)

When running multiple steps, follow:

1) `diffusivity`
2) `purkinje_slab`
3) `scar`
4) `convert`

Reason: `purkinje_slab` scales the `Diffusivity` tensor for tagged cells (using
the multiplier in `../purkinje/slab/purkinjeSlab_config.py`), and `scar` can further
scale the tensor inside scar cells using `DIFFUSIVITY_SCAR_MULTIPLIER` from
`../scar/scar_config.py`. Running out of order means those tensors won't be
modified as intended.

## Surface Tagging (Endo/Epi)

`scripts/tagEndoEpi_fromUVC.py` tags LV/RV endocardium and epicardium using
`uvc_transmural` and `uvc_intraventricular` point data. The core implementation
lives in `src/cardiac_preproc/tagging/tag_endo_epi_surface.py`. The surface
output is triangulated and saved as an unstructured surface; the volume-mapped
output stays as a volume mesh.

Defaults for surface/volume outputs and the optional seed point live in
`configs/tag_endo_epi_surface_config.py`. If you see boundary artifacts, try
setting a `SEED_POINT` (or `--seed-point`) to grow the inside-shared region
from a known location.
