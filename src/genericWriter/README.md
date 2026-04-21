# genericWriter

This folder builds `libgenericWriter`, the shared output and stimulus helper
library.

## Current contents

```text
src/genericWriter/
├── ionicModelIO.{H,C}
├── ionicVariableCompatibility.{H,C}
├── stimulusIO.{H,C}
├── activeTensionIO.{H,C}
├── ecgModelIO.{H,C}
├── purkinjeModelIO.{H,C}
├── Make/
└── README.md

```

## Purpose

This library centralizes small reusable pieces that would otherwise be repeated
inside ionic, electro, ECG, Purkinje, or active-tension code.

## Main components

- `ionicModelIO`

  General ionic-model export, trace writing, field-export planning, and
  selected-variable handling

- `ionicVariableCompatibility`

  Shared variable-name compatibility and lookup rules used by ionic exports and
  signal discovery

- `stimulusIO`

  Shared stimulus dictionary parsing and evaluation helpers

- `activeTensionIO`

  Active-tension output helpers

- `ecgModelIO`

  ECG output helpers

- `purkinjeModelIO`

  Purkinje/conduction time-series output helpers

## What this folder does not own

This folder does not own the physics models themselves. It only provides helper
logic used by:

- `ionicModels`

- `electroModels`

- `activeTensionModels`
