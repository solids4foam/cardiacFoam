# setIonicRestartState

Starts a myocardium from converged ionic states instead of the ionic model's
default initial state. Every cell gets the state of its `ionicHeterogeneity`
region, taken from a single-cell run of that region.

## Why

A single cell paced for many beats settles its slow variables (intracellular
Na⁺, K⁺, SR Ca²⁺). A 3D heart cannot afford those beats. This utility copies
each region's settled state into every cell of that region, so the 3D run
starts close to steady state.

## Inputs

1. The case's own `electroProperties` (`constant/` or `constant/<region>/`):
   `ionicModel` and the `ionicHeterogeneity` block, which must use
   `mode namedRegions`. Cells are assigned to regions exactly as the solver
   assigns the region constants; cells in a `blend` transition get the same
   weighted mix of states.
2. The region field named by `ionicHeterogeneity.field` (for example
   `uvc_transmural`), read from the start time.
3. A state map, `system/setIonicRestartStateDict` by default:

   ```
   regionStates
   {
       endoPig  "/path/to/seedStates/gaur_endo/2000.02/GaurState";
       midPig   "/path/to/seedStates/gaur_mid/2000.02/GaurState";
       epiPig   "/path/to/seedStates/gaur_epi/2000.02/GaurState";
   }
   ```

   Each file is the `<Model>State` restart file a single-cell run writes at
   its output times. Relative paths are relative to the case root.

## Output

`<startTime>/<Model>State`, or `<startTime>/<region>/<Model>State` with
`-region`. cardiacFoam reads it at start-up in place of the default state.

## Usage

```bash
setIonicRestartState                          # serial
setIonicRestartState -region electro          # EM case, electro region
mpirun -np 8 setIonicRestartState -parallel   # after decomposePar
```

In parallel each processor writes the file for its own cells, so run it after
`decomposePar`.

## Errors

It stops when a region has no state file, a map entry is not a region, a
state file belongs to another ionic model, holds more than one cell, or has a
different number of states from the others, when the region field leaves
[0, 1], or when the myocardium uses a `cellZone` subset.
