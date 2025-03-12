# Idealised Ventricle

## Overview
This 3-D problem consists of the inflation of an idealisde ventricle, as
proposed in the Land et al. (2015) benchmark article.

## Instructions

### Compile `extractIdealisedVentricleResults` Utilty
The `extractIdealisedVentricleResults` utility is used to extract results from
the idealised ventricle case.
The `extractIdealisedVentricleResults` utility can be compiled with
```bash
wmake extractIdealisedVentricleResults
```

### Run the Cases
The `Allrun` runs the a mesh study for various mesh and solution procedure
configurations. The configurations are defined near the top of the `Allrun` script:
```bash
configs=(
    "BASE=base/snes NAME=poly.hypre PETSC_FILE=petscOptions.hypre"
    "BASE=base/segregated NAME=poly.seg"
e"
)
```
where various flags are used to specify meshing and solution procedure options.
The `Allrun` script is executed as
```bash
./Allrun
```
which creates a directory for the cases called `run_<CPU_NAME>_<DATE_TIME>`, for
example, `run_Apple_M1_Ultra_20250118_151956`. When the `Allrun` script
completes, pdf plots will be available in the directory, if `gnuplot` is
installed.