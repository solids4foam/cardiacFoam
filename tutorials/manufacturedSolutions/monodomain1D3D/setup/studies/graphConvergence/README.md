# graphConvergence - monodomain1D3D

Spatial convergence of the Purkinje graph alone (`Allrun -domain graph`,
`runPurkinjeGraph`, no junction): the six graphs at a fixed time step, and the
same ladder at a 16 times smaller step as a control of the temporal error.

## Ladder

| graph | nodes | edge length h | dt (fixedDt) | dt (smallDt) |
|---|---|---|---|---|
| nodes003 | 3 | 0.5 | 0.000140174 | 8.760875e-06 |
| nodes011 | 11 | 0.1 | 0.000140174 | 8.760875e-06 |
| nodes021 | 21 | 0.05 | 0.000140174 | 8.760875e-06 |
| nodes041 | 41 | 0.025 | 0.000140174 | 8.760875e-06 |
| nodes081 | 81 | 0.0125 | 0.000140174 | 8.760875e-06 |
| nodes161 | 161 | 0.00625 | 0.000140174 | 8.760875e-06 |

Every run ends at `endTime 0.200028`, 1427 steps of 0.000140174. `Allrun`
also builds the default `20^3` cube, which `runPurkinjeGraph` needs but does
not solve on.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<fixedDt|smallDt>_nodes<k>`, the names that `../summarize.py` reads. From
this folder:

```bash
mkdir -p results
for k in 003 011 021 041 081 161
do
    rsync -a --exclude setup ../../../ results/fixedDt_nodes$k/
    results/fixedDt_nodes$k/Allrun -domain graph -graph nodes$k -dt 0.000140174 -endTime 0.200028
    rsync -a --exclude setup ../../../ results/smallDt_nodes$k/
    results/smallDt_nodes$k/Allrun -domain graph -graph nodes$k -dt 8.760875e-06 -endTime 0.200028
done
python3 ../summarize.py results
```

## Measured behaviour

Results of the last run of this study (`results/`).

`fixedDt`:

| nodes | h | Vm1D L2 | Vm1D Linf | u1 L2 | order Vm1D L2 | order Vm1D Linf | order u1 L2 |
|---|---|---|---|---|---|---|---|
| 3 | 0.5 | 6.1657e-03 | 7.7755e-03 | 9.0583e-04 | | | |
| 11 | 0.1 | 2.3340e-04 | 3.2055e-04 | 2.5902e-05 | 2.03 | 1.98 | 2.21 |
| 21 | 0.05 | 5.1270e-05 | 7.3210e-05 | 2.9048e-07 | 2.19 | 2.13 | 6.48 |
| 41 | 0.025 | 7.1165e-06 | 1.3698e-05 | 6.5135e-06 | 2.85 | 2.42 | -4.49 |
| 81 | 0.0125 | 5.1435e-06 | 1.0198e-05 | 8.0386e-06 | 0.47 | 0.43 | -0.30 |
| 161 | 0.00625 | 7.5989e-06 | 1.4148e-05 | 8.4021e-06 | -0.56 | -0.47 | -0.06 |

`smallDt`:

| nodes | h | Vm1D L2 | Vm1D Linf | u1 L2 | order Vm1D L2 | order Vm1D Linf | order u1 L2 |
|---|---|---|---|---|---|---|---|
| 3 | 0.5 | 6.1741e-03 | 7.7904e-03 | 9.1593e-04 | | | |
| 11 | 0.1 | 2.4118e-04 | 3.3505e-04 | 3.4294e-05 | 2.01 | 1.95 | 2.04 |
| 21 | 0.05 | 5.8862e-05 | 8.3242e-05 | 7.9874e-06 | 2.03 | 2.01 | 2.10 |
| 41 | 0.025 | 1.4183e-05 | 2.0099e-05 | 1.5687e-06 | 2.05 | 2.05 | 2.35 |
| 81 | 0.0125 | 3.1521e-06 | 4.5784e-06 | 1.7381e-08 | 2.17 | 2.13 | 6.50 |
| 161 | 0.00625 | 4.4004e-07 | 8.5593e-07 | 4.0300e-07 | 2.84 | 2.42 | -4.54 |

The graph error falls at second order while the spatial error dominates. At
the fixed step the spatial and temporal errors meet near 21 nodes: the `u1`
error passes through zero there and from 41 nodes on the error stays near
`1e-5`, the temporal error of `dt = 0.000140174`. With a 16 times smaller step
the same pattern moves two levels finer (second order to 41 nodes, the `u1`
sign change at 81). This ladder therefore measures the spatial order of the
graph operator only on its coarse levels; the coupled ladder
([coupledConvergence](../coupledConvergence/README.md)) refines the time step
with `h^2`.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
