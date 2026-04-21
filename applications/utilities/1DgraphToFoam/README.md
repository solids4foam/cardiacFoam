# 1DgraphToFoam

Converts a legacy ASCII VTK graph made of `VTK_LINE`, `VTK_POLY_LINE`, or
`POLYDATA LINES` into a Foam dictionary at `constant/purkinjeGraph`.

Recommended VTK data layout:

- Store graph coordinates in `POINTS`.
- Store graph connectivity as line elements only. Avoid standalone
  `VTK_VERTEX` cells for root or terminal markers.
- Store node roles in `POINT_DATA`, for example `nodeRole`:
  `0 = internal`, `1 = root`, `2 = terminal`.
- Store edge material data in `CELL_DATA` attached to line cells, for example
  `diffusivity`, `radius`, `fiberGroup`, or `branchId`.

Example:

```bash
1DgraphToFoam purkinje_network_graph.vtk -case <case> -name purkinjeGraph
```

The output contains:

- `rootNode`: from `POINT_DATA nodeRole == 1`, or node 0 if absent.
- `pvjNodes`: from `POINT_DATA nodeRole == 2`, or inferred as degree-one endpoints except the root.
- `pvjLocations`: coordinates of `pvjNodes`, used by the PVJ mapper.
- `conductionEdges`: solver-ready `(nodeA nodeB length conductance)` entries.
- `points`: VTK point coordinates.
- `edges`: retained line connectivity.
- `edgeVtkCellMap`: map from each retained edge to the original VTK cell id.
- `edgeLength`: geometric line length.
- `endpointNodes`: degree-one nodes inferred from connectivity.
- `pointFields`: preserved `POINT_DATA`.
- `edgeFields`: preserved line `CELL_DATA`, mapped through `edgeVtkCellMap`.

For `conductionEdges`, conductance is read from the first available scalar
`CELL_DATA` field in this order: `conductance`, `diffusivity`, `D`, `sigma`,
`conductivity`. If none is present, the utility writes `1.0`.
