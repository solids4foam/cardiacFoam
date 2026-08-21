# Restart state with contained propagation

## Summary

Add a dormant, opt-in restart interface at the ionic and active-tension base
classes, then activate it only for scalar `TNNP` and `LandNiederer` in the
sequential electromechanics path. No other model, solver, or dictionary
changes behavior in this change.

## Propagation and containment

- Add non-pure virtual hooks to `ionicModel` and `activeTensionModel`:
  - `bool supportsRestartState() const { return false; }`
  - `bool readRestartState(const fvMesh&) { return false; }`
  - `void writeRestartState(const fvMesh&) const {}`
- Keep the base defaults inert. This preserves every scalar, batched/GPU,
  manufactured, eikonal, and single-cell model unless it explicitly opts in.
- Override the hooks only in `TNNP` and `LandNiederer`; share only the binary
  file framing/validation utility, not generic state restoration through
  `ioStatesPtr()`.
- Add narrow electro-model forwarding methods with no-op defaults.
  `electrophysiologyModel` forwards only its owned myocardium ionic model;
  `singleCellSolver`, eikonal, conduction, and other electro paths remain
  unchanged.
- Override `sequentialElectroMechanical::writeFields()` to invoke:
  - TNNP bundle write through the full electrophysiology path using the
    electro mesh.
  - LandNiederer bundle write using the solid mesh.
  - Existing field-writing behavior unchanged.
- During sequential-model startup, load the active-tension bundle after
  `D`/`f0` and the ionic provider exist. Only a successful Land load
  suppresses preconditioning.

## Model-specific restoration

- `TNNPState` stores exactly 17 states per local cell. On load, validate
  metadata, restore state arrays, synchronize Vm from the restarted `Vm`
  field, and recompute rates, algebraics, current, and Ca_i without
  integration.
- `LandNiedererState` stores six states plus `prevLambda`. On load, validate
  and restore both; recompute fibre stretch from restarted `D`, then refresh
  algebraics, rates, Ca_i-driven tension, and `Ta` without changing
  `prevLambda`.
- Missing state file means fresh start. Invalid, wrong-model, wrong-size, or
  non-finite state file is fatal rather than silently resetting.
- Bundles are binary, versioned, per processor, and written only at normal
  output times. No new case-dictionary keys or `0/` fields are added.

## Early-step test plan

- Use the coupled TNNP + LandNiederer tutorial with a tiny test overlay:
  output after an early step and stop after only a few steps.
- Compare:
  1. an uninterrupted run over those first steps;
  2. a run through the first output, restarted there, then completed to the
     same few steps.
- At each resulting early output, compare Vm, Ca_i, Ta, D/displacement, and
  stress within solver/output tolerance.
- Assert TNNP/Land bundles exist and that the restarted log skips Land
  preconditioning.
- Run serial first; add a small decomposed variant that verifies rank-local
  bundle validation and the same restart result.

## Assumptions

- Version 1 intentionally excludes all batched/GPU, manufactured, non-TNNP
  ionic, non-Land active-tension, eikonal, and single-cell restart persistence.
- Restart with a changed mesh or decomposition is unsupported and fails
  validation.
