# heartSim3D-1D Tutorials

This directory contains advanced tutorials modeling cardiac electrophysiology on realistic heart geometries with a coupled 1D Purkinje network.

The tracked tutorials are:
- `eikonalHeart`
- `monodomain3D-eikonal1D`
- `monodomainHeart`

> [!WARNING]
> **Missing Assets Notice**
> 
> The configuration dictionaries, Python post-processing scripts, and case `README` files are all up to date and tracked in version control. However, **the physical 3D heart meshes and the 1D `purkinjeGraph` data files are currently NOT tracked in this repository.**
> 
> If you execute the `./Allrun` scripts on a fresh clone, the solver will attempt to fall back to generating a mesh via `blockMesh`. Please note that the included `blockMeshDict` files only generate a small 20x3x7 mm slab (identical to the Niederer benchmark), rather than an actual heart geometry, and the runs will fail when attempting to read the missing `purkinjeGraph_healthy` files.
> 
> You must manually supply the `constant/polyMesh` and the `constant/purkinjeGraph_healthy` files before running these cases.
