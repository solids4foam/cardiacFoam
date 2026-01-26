# conductionSystem_preProcessing

Utilities for conduction system pre-processing (Purkinje, diffusivity, scar, conversion). The focus is on:
- tagging a Purkinje layer from UVC fields
- adding diffusion tensors
- converting FIELD arrays to SCALARS/VECTORS/TENSORS and cleaning blank lines

# Requirements
- Python packages: `pyvista`, `numpy`
- Input VTK must be legacy ASCII

# Folder Map
- `lib/` core algorithms (Purkinje slab tagging, diffusion tensor, scar)
- `utils/` VTK parsing, inspection, conversion helpers
- `output/` default outputs
- `purkinje_network/purkinje_slab/purkinje_vtk.py`, `diffusivity/diffusionTensor_vtk.py`, `../scar_creator/scar_vtk.py`, `fileConversion/ASCIIlegacyToVtkUnstructured.py` command-line tools

# CLIs

## conductionSystem_Generation.py
Run the full pipeline (Purkinje slab → diffusivity → scar → conversion).

```
python conductionSystem_Generation.py
```

Common overrides:
```
python conductionSystem_Generation.py --input ASCIIlegacy.vtk --purkinje-output output/purkinjeLayer.vtk --diffusivity-output output/Diffusion_purkinjeLayer.vtk
```


## purkinje_vtk.py
Create a Purkinje layer using UVC fields.

```
python purkinje_network/purkinje_slab/purkinje_vtk.py --mode layer --input ASCIIlegacy.vtk
python purkinje_network/purkinje_slab/purkinje_vtk.py --mode layer --input ASCIIlegacy.vtk --output output/purkinjeLayer.vtk
```
Output: VTK with `purkinjeLayer` in `CELL_DATA`.

Inputs:
- `POINT_DATA`: `uvc_transmural`, `uvc_intraventricular`

## diffusionTensor_vtk.py
Add diffusion tensors (optionally scaled in Purkinje cells), then convert + clean.

```
python diffusivity/diffusionTensor_vtk.py --input ASCIIlegacy.vtk 
python diffusivity/diffusionTensor_vtk.py --input output/purkinjeLayer.vtk --output output/Diffusion_purkinjeLayer.vtk --purkinje-mult 2.0
```
Output: VTK with `Diffusivity` in `CELL_DATA`, converted and cleaned for OpenFOAM.

Inputs:
- `CELL_DATA`: `fiber`, `sheet`
- Optional `CELL_DATA`: `purkinjeLayer` (defaults to zeros if missing)


## scar_vtk.py
Tag scar cells using a selection mesh and change the scar area for 0 zero diffusivity in scar cells.

```
python ../scar_creator/scar_vtk.py --full-mesh output/purkinjeLayer_Diffusivity_IDsGlobal.vtk --selection output/scar_tissue_region.vtu --output output/purkinjeLayer_Diffusivity_scar.vtk --scar-value 1.0 --diffusivity-scale 0.1
```
Output: VTK with `Scar` in `CELL_DATA` (and `Diffusivity` scaled in scar cells when present).

Inputs:
- `CELL_DATA`: `GlobalCellIds` on full mesh and selection
- Optional `CELL_DATA`: `Diffusivity` (scaled in scar cells)


## ASCIIlegacyToVtkUnstructured.py
Inspect and convert FIELD arrays, then remove blank lines.

```
python fileConversion/ASCIIlegacyToVtkUnstructured.py --input ASCIIlegacy.vtk --output output/converted.vtk
```
Output: converted VTK + printed field summary.

# Config Defaults
`config_Diffusivity.py` and `config_purkinjeSlab.py` store default parameters used by the CLIs.
Override per-run with CLI flags where available.

# Recommended Workflow
1) `purkinje_network/purkinje_slab/purkinje_vtk.py` (optional)
2) `diffusivity/diffusionTensor_vtk.py`
3) `../scar_creator/scar_vtk.py`
4) `fileConversion/ASCIIlegacyToVtkUnstructured.py` (if you need a separate conversion/inspection step)
