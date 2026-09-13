# electroModels/core

The brain of the electrophysiology: it builds the system of domains a case needs and advances it each timestep.

## What's available

- `electroModel`: the top-level electrophysiology model that the solver advances.
- `electrophysiologyModel`: the tissue model, available as `monodomainSolver`, `bidomainSolver` and `eikonalSolver`.
- `electrophysicsSystem`: holds the domains and couplers; `electrophysicsSystemBuilder` creates them.
- `staggeredElectrophysicsAdvanceScheme`: the advance scheme. It advances each domain in turn, with the couplers in between.
- The interfaces the domains implement: `electroDomainInterface`, `electroStateProvider`, `electroStateDomain`, `electroVolumeFieldDomain`.
- The abstract verifier bases, which the concrete verifiers in [src/verificationModels](../../verificationModels/README.md) inherit from.

## Folders

```text
core/
├── system/                  # electrophysicsSystem and its builder
├── advanceSchemes/          # the staggered advance scheme
├── electrophysiologyModel/  # the tissue model
├── verificationModels/      # abstract verifier bases
└── *.H                      # the top-level model and the domain interfaces
```

**Deep dive:** [ARCHITECTURE.md](ARCHITECTURE.md) explains how the core is built inside.
