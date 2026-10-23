# redistributeRestartState

Moves the restart state files that cardiacFoam's ionic and active-tension
models write, `<time>/<Model>State`, across the processor boundary.
`decomposePar` and `reconstructPar` handle OpenFOAM fields and leave these
files alone, so without this utility a reconstructed case has no ionic state
to restart from, and a serial state cannot start a parallel run.

## What it does

- Default: for every selected time, gathers
  `processor*/<time>/<Model>State` into `<time>/<Model>State` using each
  processor's `constant/polyMesh/cellProcAddressing`.
- `-decompose`: scatters `<time>/<Model>State` into the processor
  directories the same way.
- `-region <name>`: the files of one mesh region, for example `electro` or
  `solid` in an electromechanics case.

A state file holds one row per cell. A file whose row count is not the cell
count is reported and skipped. Every file written here reads back with the
same `restartStateIO` header the models use.

A Purkinje network's ionic state, `<time>/<Model>State.<network>`, holds one
row per graph node and is written once in the case root, where every
processor reads its own block on restart. It needs no redistribution, and a
network restarts on any number of processors.

## Usage

```bash
reconstructPar
redistributeRestartState                      # all times
redistributeRestartState -latestTime
redistributeRestartState -region electro      # one region of a multi-region case

decomposePar
redistributeRestartState -decompose -time 0.4 # start a parallel run from a serial state
```

The tutorials' `Allrun parallel` branches call it right after
`reconstructPar`, so a reconstructed tutorial is a complete restart point.
