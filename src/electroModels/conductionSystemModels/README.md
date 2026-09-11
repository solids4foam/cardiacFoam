# conductionSystemModels

The Purkinje-network solvers, used by `conductionSystemDomain`.

## What's available

- `monodomain1DSolver`: the cable equation on the network graph, with an ionic model at each node.
- `eikonalSolver1D`: activation times only, from edge lengths and a conduction speed.
- `restitutionEikonalSolver1D`: activation times with conduction-velocity restitution. It can re-excite, tracks refractoriness, and flags conduction block and wavebreak.

## Folders

```text
src/electroModels/conductionSystemModels/
├── monodomain1DSolver/
├── eikonalSolver1D/
└── restitutionEikonalSolver1D/
```
