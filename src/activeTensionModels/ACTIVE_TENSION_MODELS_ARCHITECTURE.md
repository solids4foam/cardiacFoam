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
| `ManufacturedElectromechanics` | `1.0` |

This table previously carried a `GoktepeKuhl` row at `1.0` with an open
question about whether that was correct. It was not: the model's `eps`
switch was parameterised in Göktepe & Kuhl's mV/ms units but evaluated on
the normalised `u`, leaving it 100x slow with its switch effectively
disabled. The model has been removed rather than repaired — see
[`future/GOKTEPE_KUHL_REMOVAL.md`](../../future/GOKTEPE_KUHL_REMOVAL.md) for
the analysis and for what to do differently if it is ever re-added.

## Concrete models

### `NashPanfilov`

- runtime name: `NashPanfilov`
- integration-point ODE model
- uses the shared base-class export and write machinery
- faithful to Nash & Panfilov (2004) Eq. (22c)/(23) and Table 1

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
