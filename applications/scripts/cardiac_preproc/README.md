# Cardiac Preprocessing Engine

Core engine for preprocessing workflows (diffusivity, scar, tagging, and
Purkinje slab/fractal orchestration).

Upper runtime facade lives in:

- `../engine/src/cardiac_engine/`

## Structure

- `../conduction_system_generation.py` — top-level main driver (uses `--steps` CLI)
- `../tagging/` — top-level endo/epi tagging entrypoint + config:
  - `tagging.py`
  - `config_tagging.py`
- `../conversion/` — top-level conversion entrypoint + config:
  - `conversion.py`
  - `config_conversion.py`
- `../filesOrganize/` — centralized data assets (inputs/outputs/VTK files)
- `src/cardiac_preproc/` — implementation code
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
python3 ../conduction_system_generation.py --steps diffusivity purkinje_slab
```

Run a single tool:

```bash
python3 ../purkinje/slab/slab.py
python3 ../diffusivity/diffusivity.py
python3 ../scar/scar.py
python3 ../purkinje/fractal_3d/fractal.py
python3 ../purkinje/fractal_3d/fractal.py --input ../filesOrganize/cardiac_preproc/outputs/ASCIIlegacy_DiffusionTensor_endo_epi_surface.vtk
python3 ../tagging/tagging.py
```
By default the fractal run opens an interactive viewer after VTUs are written.
To skip it:

```bash
python3 ../purkinje/fractal_3d/fractal.py --no-view-final
```

## Configuration

Use `../conduction_system_generation.py --steps ...` to control which steps run.
Each solver reads defaults from colocated product configs:
- `../diffusivity/config_diffusivity.py`
- `../scar/config_scar.py`
- `../purkinje/slab/config_slab.py`

If you want to understand the detailed logic for a step, see the corresponding
canonical module or the implementation in `src/cardiac_preproc/`.

## Pipeline Order (Why it matters)

When running multiple steps, follow:

1) `diffusivity`
2) `purkinje_slab`
3) `scar`
4) `convert`

Reason: `purkinje_slab` scales the `Diffusivity` tensor for tagged cells (using
the multiplier in `../purkinje/slab/config_slab.py`), and `scar` can further
scale the tensor inside scar cells using `DIFFUSIVITY_SCAR_MULTIPLIER` from
`../scar/config_scar.py`. Running out of order means those tensors won't be
modified as intended.

## Surface Tagging (Endo/Epi)

`../tagging/tagging.py` tags LV/RV endocardium and epicardium using
`uvc_transmural` and `uvc_intraventricular` point data. The core implementation
lives in `src/cardiac_preproc/tagging/tag_endo_epi_surface.py`. The surface
output is triangulated and saved as an unstructured surface; the volume-mapped
output stays as a volume mesh.

Defaults for surface/volume outputs and the optional seed point live in
`../tagging/config_tagging.py`. If you see boundary artifacts, try
setting `SHARED_BOUNDARY_SEED` (or `--seed-point`) to grow the inside-shared region
from a known location.
