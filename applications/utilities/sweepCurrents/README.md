# sweepCurrents

Sweeps one or more ionic currents over a voltage range and writes the
current (and its gate-variable dependencies) as a function of clamped
membrane voltage `Vm`. Used to verify current–voltage (I–V) relationships
and to calibrate ionic model parameters.

## What it does

1. Reads `constant/sweepCurrents` to determine the ionic model, tissue type,
   target currents, and voltage sweep range.
2. Constructs the ionic model with a single integration point.
3. For each requested current, calls `ionicModel::sweepCurrent(...)` which
   clamps `Vm` at each voltage point and records the current and its
   gate-variable dependencies.
4. Writes one CSV file per current to `postProcessing/`:
   - Filename: `<model>_<tissue>_<current>_sweep.<ext>`
   - Columns: `Vm`, then one column per dependency variable.

## Dictionary: `constant/sweepCurrents`

```foam
ionicModel      TNNP;
tissue          endocardialCells;

Vmin            -100;   // mV
Vmax             60;    // mV
points           100;   // number of voltage points

currents        (INa IKr IKs);   // omit to sweep all available currents

listCurrentsOnly false;    // set true to only print available current names
outputExtension  csv;      // output file extension (default: txt)
printCurrentVariables true; // print dependency variable names to terminal
```

## Usage

```bash
# Sweep all available currents
sweepCurrents -case <caseDir>

# List available current names without sweeping
sweepCurrents -case <caseDir>   # with listCurrentsOnly true in dict

# Specific currents
sweepCurrents -case <caseDir>   # with currents (INa IKr); in dict
```

## Notes

- Runs in serial only (`noParallel`).
- `listCurrentsOnly true` (or `listCurrents true`) prints available names
  and exits without sweeping.
- Only currents declared in the model's sweep-current dependency map are
  available. Use `listCurrentsOnly true` to inspect the supported set for
  a given model.
- Output files are appended to `postProcessing/` in the case directory.
