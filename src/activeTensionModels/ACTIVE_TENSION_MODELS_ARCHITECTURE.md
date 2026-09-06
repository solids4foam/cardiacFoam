# Active Tension Models Architecture

## Overview

`src/activeTensionModels` contains runtime-selectable active-tension models
built around the base class `Foam::activeTensionModel`.

The current implementation is:

- integration-point based
- scalar-state based
- driven by an upstream `ElectromechanicalSignalProvider`

## Base class: `activeTensionModel`

The base class lives in:

- `activeTensionModel/activeTensionModel.H`
- `activeTensionModel/activeTensionModel.C`

### Main responsibilities

- runtime selection via `activeTensionModel::New(...)`
- own the model dictionary and integration-point count
- discover required upstream signals through `Requirements`
- bind an optional `ElectromechanicalSignalProvider`
- run the shared `calculateTension(...)` loop
- provide generic export and write helpers through `activeTensionIO`

### Signal-side contract

The base class depends on:

- `ElectromechanicalSignalProvider`

This is how the model asks for signals such as:

- `Vm` (canonical mV)
- `Cai` (canonical mM)

`driveSignalScaleFactor()` is applied in the common scalar and batched
signal-read paths. It defaults to `1.0`; a model can override it when its
published equations require a different input unit without changing the
provider contract. The original `LandNiederer` variants use `1000.0` to
consume `Cai` in µM.

The `LandNiedererTWorld` variants deliberately leave it at `1.0`: their
maths header takes `Cai` in mM and performs the mM→µM conversion itself. The
two Land families therefore convert in different places, and both are
correct — check the input-convention block at the top of a model's maths
header before adding or removing a scale factor.

### Time-side contract

`timeScaleFactor()` converts OpenFOAM time (seconds) into the model's own
time unit, and `scaledTime()` applies it. Unlike `driveSignalScaleFactor()`
it is **pure virtual**, mirroring `batchedIonicModel`: a defaulted `1.0`
cannot distinguish "this model runs on the OpenFOAM clock" from "nobody
thought about it", so every model states its own scale even when that scale
is `1.0`.

| Model | `timeScaleFactor()` |
| --- | --- |
| `NashPanfilov`, `NashPanfilovBatched` | `1000/12.9` (Aliev-Panfilov dimensionless time) |
| `LandNiederer`, `LandNiedererBatched` | `1000` (ms) |
| `LandNiedererTWorld`, `LandNiedererTWorldBatched` | `1000` (ms) |
| `GoktepeKuhl`, `GoktepeKuhlBatched` | `1.0` |
| `ManufacturedElectromechanics` | `1.0` |

`GoktepeKuhl` shares its rate equation and constants with `NashPanfilov` but
does not rescale time. Whether that is correct is a question about the
CellML/paper provenance of those constants, not about this hook.

## Concrete models

### `GoktepeKuhl`

- runtime name: `GoktepeKuhl`
- integration-point ODE model
- uses the shared base-class export and write machinery

### `NashPanfilov`

- runtime name: `NashPanfilov`
- integration-point ODE model
- follows the same provider and I/O pattern as `GoktepeKuhl`

### `LandNiederer`

- runtime name: `LandNiederer`
- original seven-state intact-human Land et al. (2017) model
- reports active (`AV_Ta`), passive (`AV_Tp`), and total (`AV_T`) tension;
  only active tension is supplied to the active-stress interface
- batched CPU/GPU runtime name: `LandNiedererBatched`

### `LandNiedererTWorld`

- runtime name: `LandNiedererTWorld`
- six-state contraction subsystem extracted from TWorld
- GPU-batched runtime name: `LandNiedererTWorldBatched`

For the folder-level overview, see [`README.md`](./README.md).
