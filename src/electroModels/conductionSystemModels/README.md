# conductionSystemModels

The Purkinje-network solvers, used by `conductionSystemDomain`.

## What's available

- `monodomain1DSolver`: the cable equation on the network graph, with an ionic model at each node.
- `eikonalSolver1D`: activation times only, from edge lengths and a conduction speed.
- `restitutionEikonalSolver1D`: activation times with conduction-velocity restitution. It can re-excite,
  tracks refractoriness, and flags conduction block and wavebreak. Recovery is a fixed-duration
  surrogate: a node accepts a new activation once its activation interval reaches
  `apdNominal + minimumDI90`. It stores no repolarization time, so the `DI` it reports is an
  activation interval minus a constant, not a measured DI90. The capture boundary (`minimumDI90`) and
  the CV table's domain are separately calibrated quantities.

## Folders

```text
src/electroModels/conductionSystemModels/
├── monodomain1DSolver/
├── eikonalSolver1D/
└── restitutionEikonalSolver1D/
```
