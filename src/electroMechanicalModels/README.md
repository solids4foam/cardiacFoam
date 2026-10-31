# electroMechanicalModels

The electromechanics part of the core: it couples the electrophysiology to a solid, using solids4foam as the solid backend. Built only in full mode.

## What's available

- `sequentialElectroMechanical`: each timestep it advances the electrophysiology, computes the active tension `Ta` with the chosen active-tension model, and advances the solid, in one pass with no outer correctors. `electroMechanicalModel` is its base class; it follows solids4foam's `fluidSolidInterface` pattern, which is also the intended route to fluid–structure interaction.

The active tension reads `Vm` or `Cai` from the electrophysiology's ionic model. With `activeTensionSignal templates;` it reads them instead from single-cell templates, one per `ionicHeterogeneity` region, evaluated at `t - activationTime` in each cell. This drives the mechanics from an eikonal activation: the activation is solved once in the first step and the controlDict `deltaT` drives the run. `activeTensionTemplates` takes the same entries as `eikonalECG`'s `personalizedTemplates` ([eikonalECG README](../electroModels/ecgModels/eikonalECG/README.md)):

```cpp
sequentialElectroMechanicalCoeffs
{
    activeTensionModel   LandNiedererTWorld;
    activeTensionSignal  templates;   // electro (default) | templates
    activeTensionTemplates
    {
        ionicModelConfig { ionicModel TNNP; tissue epicardialCells; solver RKF45; singleCellStimulus { ... } }
        nBeats     1;
        duration   0.5;
        dt         1e-5;
    }
}
```

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
