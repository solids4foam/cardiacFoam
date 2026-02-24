# Costabal 2015 Purkinje Network Generator

This folder contains a standalone Python script that generates left- and
right-ventricular Purkinje fractal trees (Costabal 2015 style) on a tagged
surface mesh.

## Main script

- `fractalTree_purkinjeNetwork.py`: loads a triangulated surface mesh (VTK/VTU
  or gmsh `.msh`), picks LV/RV seed points and initial directions, then generates
  two fractal trees using the local `fractal_tree` package.

## Inputs

- Mesh: a triangulated surface with LV/RV tags (VTK/VTU recommended). Default:
  `inputs/meshes/ASCIIlegacy_biventricular_endo_epi_surface_triangulate.vtk`.
- Config: `../purkinjeFractalTree_config.py` (Python file with parameters).

The mesh must contain surface tags named (mapped in config):
- `ENDO_LV` (required)
- `ENDO_RV` (required)
- `EPI` (optional; used only for plotting)

## Configuration parameters (`../purkinjeFractalTree_config.py`)

- `lv_filename`, `rv_filename`: output file name bases.
- `lv_seed`, `rv_seed`: initial seed points (nearest mesh vertex is used).
- `lv_init_dir`, `rv_init_dir`: initial direction vectors.
- `lv_init_length`, `rv_init_length`: initial segment length.
- `lv_N_it`, `rv_N_it`: number of growth iterations.
- `lv_length`, `rv_length`: branch segment length.
- `lv_save_frames`, `rv_save_frames`: write intermediate frames for videos.

## Outputs

Written under `outputs/`:
- `*_lines.txt`: line connectivity for the tree.
- `*_xyz.txt`: node coordinates.
- `*_endnodes.txt`: end-node list.
- `<input_stem>_geometry.vtu`: ParaView-friendly mesh.

If `*_save_frames` is enabled, per-iteration frames are also saved.

## Usage

```bash
python3 fractalTree_purkinjeNetwork.py --input <surface.vtk> --config ../purkinjeFractalTree_config.py
```

Optional arguments:

```bash
python3 fractalTree_purkinjeNetwork.py --input ../../../../filesOrganize/cardiac_preproc/inputs/meshes/biv_ellipsoid.msh --config ../purkinjeFractalTree_config.py
```
By default the script opens an interactive viewer after the final VTUs are written.
Use `--no-view-final` to skip it.
