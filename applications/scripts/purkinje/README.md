# Purkinje Product Architecture

This folder is the top-level Purkinje product surface and is split into three
child domains:

- `density_mapper/`: probability engine + bullseye tooling + VTK mapping.
- `fractal_3d/`: fractal network generation workflow and visualization tooling.
- `slab/`: slab tagging entrypoint and user-facing wrapper.

The compute engine remains centralized in:

- `cardiac_preproc/src/cardiac_preproc/`

Design rule:

- `purkinje/*` folders are product-facing modules and wrappers.
- `cardiac_preproc` contains reusable engine/domain logic and orchestration.
- `engine/` contains the upper shared runtime facade.
