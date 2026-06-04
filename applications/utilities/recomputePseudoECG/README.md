# recomputePseudoECG

Recomputes the pseudo-ECG signals from stored `Vm` and conductivity fields
without re-running the full 3-D simulation. Useful when electrode positions
need to be corrected or a new electrode configuration needs to be evaluated.

## What it does

1. Reads `constant/electroProperties` and extracts the electrode positions from
   the first `ecgDomains` entry.
2. Reads the conductivity tensor once from `0/Diffusivity`, or from the active
   solver coefficients when the field is not present on disk.
3. Iterates over every stored time directory (or a selected range).
4. For each time step:
   - Reads the `Vm` field.
   - Computes `fvc::grad(Vm)` via the Green-Gauss scheme.
   - Applies the Gima-Rudy dipole integral over all mesh cells:

         phi(P) = sum_c  (sigma_c . grad(Vm)_c) . (x_c - P) * V_c
                          / |x_c - P|^3

5. Writes `postProcessing/pseudoECG.dat` (header + one row per time step).

The output file is a drop-in replacement for the `pseudoECG.dat` written
during a normal `cardiacFoam` run and is compatible with `plot_pseudo_ecg.py`.

## Usage

```bash
# Recompute all stored time steps
recomputePseudoECG

# Selected time range
recomputePseudoECG -time '0.1:0.5'

# Write to a separate file (keep the original)
recomputePseudoECG -output postProcessing/pseudoECG_corrected.dat

# Parallel
mpirun -np 6 recomputePseudoECG -parallel
```

## Options

| Option | Default | Description |
|---|---|---|
| `-output <path>` | `postProcessing/pseudoECG.dat` | Output file path |
| `-vmField <name>` | `Vm` | Name of the membrane potential field to read |
| `-sigmaField <name>` | `Diffusivity` | Conductivity tensor field name in `0/` |
| `-time <range>` | all | OpenFOAM time selector, e.g. `'0.1:0.5'` or `'latest'` |

## Notes

- The conductivity field is read once from `0/Diffusivity` and assumed
  constant in time. If the field is not present, the utility constructs the
  tensor from the same active `<solver>Coeffs` dictionary used by
  `cardiacFoam`.
- Electrode positions are taken from the first `ecgDomains` sub-dictionary in
  `constant/electroProperties`. Update the positions there and re-run this
  utility to correct any errors without re-running the simulation.
- Runs in parallel; results are reduced across processors at each time step.
- Approximately 50× faster than a full simulation re-run since no ODE solving
  or Purkinje coupling is performed.
