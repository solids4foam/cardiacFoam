# electroMechanicalModels

This library builds `libelectroMechanicalModels`, the electro-mechanics coupling
layer for full solids4foam builds. It is **not compiled** in lightweight
electro-only builds (see `etc/resolveSolids4Foam.sh`).

## Current contents

```text
src/electroMechanicalModels/
├── electroMechanicalModel/       # Abstract base class and runtime selection
├── sequentialElectroMechanical/  # Concrete sequential (weakly coupled) scheme
└── Make/
```

## Core class: `Foam::electroMechanicalModel`

Defined in `electroMechanicalModel/electroMechanicalModel.H`.

- Inherits `physicsModel` (the top-level cardiacFoam/solids4foam entry point).
- Owns an `electroModel` and a `solidModel` in separate mesh regions.
- Declares the runtime selection table for concrete coupling schemes.
- Follows the `fluidSolidInterface` pattern from solids4foam.

Key virtual interface:

- `evolve()` — advance one time step (pure virtual, implemented by derived classes)
- `writeFields(const Time&)` — write fields for both sub-models
- `setDeltaT(Time&)` — propagate time-step updates
- `end()` — cleanup

## Available coupling scheme

- **`sequentialElectroMechanical`** (registered as `sequentialElectroMechanical`)
  - Weakly coupled: no outer correctors per time step.
  - Sequence: electro `evolve()` → compute active tension from ionic Cai signal
    → inject `Ta` field into solid mesh → solid `evolve()`.
  - `Ta` is computed as a linear function of Cai above a threshold:
    `Ta = kTa * max(Cai - CaiThreshold, 0)`.
  - `kTa` and `CaiThreshold` are read from the solver dictionary.

## Build mode

This library is compiled only when `etc/resolveSolids4Foam.sh` finds a valid
solids4foam installation. In electro-only builds, `libelectroMechanicalModels`
is not built and `electroMechanicalModel` is not available as a `physicsModel`
type.

## Selected by

```text
physicsProperties:
    type    electroMechanicalModel;
```

and then:

```text
electroMechanicalProperties:
    electroMechanicalModel  sequentialElectroMechanical;
```
