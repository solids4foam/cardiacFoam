# electroMechanicalModels

The electromechanics part of the core: it couples the electrophysiology to a solid, using solids4foam as the solid backend. Built only in full mode.

## What's available

- `sequentialElectroMechanical`: each timestep it advances the electrophysiology, computes the active tension `Ta` with the chosen active-tension model, and advances the solid, in one pass with no outer correctors. `electroMechanicalModel` is its base class; it follows solids4foam's `fluidSolidInterface` pattern, which is also the intended route to fluid–structure interaction.

The solid needs the fibre fields `f0` and `f0f`. [setFibreField](../../applications/utilities/setFibreField/README.md) writes both; if `f0` comes from elsewhere, [interpolateFibreField](../../applications/utilities/interpolateFibreField/README.md) derives `f0f` from it.

## Folders

```text
src/electroMechanicalModels/
├── electroMechanicalModel/       # base class
└── sequentialElectroMechanical/  # the sequential, weakly coupled scheme
```

## What this does not own

- The active-tension models: [activeTensionModels](../activeTensionModels/README.md).
- The electrophysiology: [electroModels](../electroModels/README.md).
- The solid solver and its constitutive law (`electroMechanicalLaw`): solids4foam.
