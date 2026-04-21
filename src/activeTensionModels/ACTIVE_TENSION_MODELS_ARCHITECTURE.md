# Active Tension Models Architecture

This note describes the active-tension layer that actually exists in this tree.
It is not a generic electromechanics survey.

## Overview

`src/activeTensionModels` contains runtime-selectable active-tension models
built around the base class `Foam::activeTensionModel`.

The current implementation is:

- integration-point based
- scalar-state based
- driven by an upstream `ElectromechanicalSignalProvider`

It does not implement a full tensor-mechanics constitutive framework by itself.

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

## Folder responsibilities

### Owned here

- active-tension model runtime selection
- per-integration-point active-tension state
- signal-provider requirements
- model-side export and trace writing hooks

### Not owned here

- myocardium or ionic-model time integration
- solid-mechanics constitutive law infrastructure
- a separate `electroMechanicalModel` implementation in this source tree

## Dependency relationship

The active-tension layer sits alongside electrophysiology rather than inside
`electroModels`:

- upstream signals come from `ionicModel` through
  `ElectromechanicalSignalProvider`
- active-tension output is evolved by `activeTensionModel`
- output helpers come from `activeTensionIO`

For the folder-level overview, see [`README.md`](./README.md).
