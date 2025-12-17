# singleCellElectroActivationFoam tutorial

This tutorial demonstrates how to run **single-cell electrophysiology
simulations** using `singleCellElectroActivationFoam`.

The case is designed to test and analyse a **single integration point**
of a chosen ionic model, allowing inspection of the transmembrane voltage
and internal state variables as functions of time.

---

## Overview

Unlike tissue-level solvers (e.g. `electroActivationFoam`), this case does
**not** solve a spatial voltage PDE. Instead:

- The transmembrane voltage `Vm` is integrated directly as part of the
  ionic model ODE system,
- A single cell (one integration point) is advanced in time,
- Results are written as a time series.

This setup is useful for:

- Testing and validating ionic model implementations,
- Producing single-cell action potentials,
- Debugging ODE solvers and stimulus protocols,
- Exploring parameter sensitivity at the cell level.

---

## Case structure

```
singleCell/
├── constant
│   ├── cardiacProperties
│   ├── electroActivationProperties
│   ├── stimulusProtocol
│   └── timeIntegrationProperties
│   └── sweepCurrents

├── system
│   ├── blockMeshDict
│   └── controlDict
├── plotVoltage
├── setupSingleCell
│   ├── main_singleCell.py
│   ├── run_cases.sh
│   ├── setup_multiple_simulations_singleCell.py
│   └── singleCellinteractivePlots.py
├── Allrun
└── Allclean
```

---

## Model configuration

### Ionic model

The ionic model is selected in `constant/cardiacProperties`, for example:

```cpp
ionicModel  BuenoOrovio;
```

Any model registered with the `ionicModels` library can be used - write `banana` (or any other word) and run the solver to see a list of available models.

---

### Time integration

The file `constant/timeIntegrationProperties` controls:

- The ODE solver (e.g. `RKF45`, `Euler`),
- Initial ODE step size,
- Solver tolerances,
- Time-integration behaviour.

Since this is a single-cell case, the voltage equation is solved entirely
within the ionic model ODE system.

---

### Stimulus protocol

The external stimulus is defined in `constant/stimulusProtocol`, allowing
specification of:

- Stimulus start time,
- Duration,
- Amplitude,
- Pacing protocols (S1/S2).

This enables standard single-cell pacing and restitution studies.

---
### SweepCurrents

The dictionary to acess to run the sweepCurrent utility besides ionicModel and tissue type, the user decides: 
- Voltage interval (typicaly -80 to + 40),
- number of points to interpolate.

This enables ionicModel checks and comparisons of variables that define some currents.

-  steady state gating,
-  time constant tau.

Note that BuenoOrovio has dimensionless u - i.e 0 to +1.5

---

## Running the case

From the case directory:

```bash
./Allrun
```

This will:

1. Run `singleCellElectroActivationFoam`,
2. Automatically invoke `plotVoltage` after the simulation completes.

To clean generated files:

```bash
./Allclean
```

---

## Post-processing

### Voltage trace

The script `plotVoltage` plots **transmembrane voltage vs time** from the
solver output.

This provides a quick visual check of:

- Action potential morphology,
- Depolarisation and repolarisation timing,
- Stimulus response.

---

### Advanced analysis

The `setupSingleCell` directory contains Python utilities for:

- Running multiple single-cell simulations,
- Automating parameter sweeps,
- Generating interactive plots of voltage and state variables.

These scripts are useful for systematic model analysis and development.

---

## Purpose of this tutorial

This case is intended for:

- Single-cell model verification,
- Development and debugging of ionic models,
- Educational demonstrations of cardiac action potentials.

It is **not** intended for tissue-level or patient-specific simulations.

---

## Notes

- Results depend strongly on the chosen ionic model and stimulus protocol.
- Time-step sensitivity should be checked when changing ODE solver settings.
- For spatial propagation and activation timing, use `electroActivationFoam`
  instead.
