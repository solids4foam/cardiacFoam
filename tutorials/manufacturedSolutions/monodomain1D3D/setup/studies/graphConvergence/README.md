# graphConvergence - monodomain1D3D

## Purpose

Mesh-refinement convergence of the 1D Purkinje graph alone, with no
myocardium coupling: the case's own README "Graph-Only Diagnostic"
(`blockMesh -dict system/blockMeshDict.3D`, then `runPurkinjeGraph`),
reached through the record's `graphOnly` route (`"solver": "graphOnly"`)
instead of by hand. Traces convergence across the six committed
`constant/purkinjeGraph.nodes*` refinements (003 through 161 nodes) at the
native mesh and time-step defaults.

## Execution

From the repository root:

```bash
[omnidriver command to run]
```

## Tracking & Outputs

All generated outputs are saved to the local `results/` folder, which is
explicitly ignored by git. Do not commit generated OpenFOAM data.
