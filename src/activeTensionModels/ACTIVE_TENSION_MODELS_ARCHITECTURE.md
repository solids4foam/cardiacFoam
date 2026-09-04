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

- `Vm`
- `Cai`

Current concrete models request `Vm`.

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

### `LandNiedererTWorld`

- runtime name: `LandNiedererTWorld`
- six-state contraction subsystem extracted from TWorld
- GPU-batched runtime name: `LandNiedererTWorldBatched`

For the folder-level overview, see [`README.md`](./README.md).
