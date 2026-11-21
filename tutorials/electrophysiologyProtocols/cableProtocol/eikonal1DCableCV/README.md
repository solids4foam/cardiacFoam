# eikonal1DCableCV tutorial

## What it does

Solves the eikonal equation for activation time (`eikonalSolver`,
advection-diffusion form, no ionic model) on a 20 mm cable of
0.1 mm x 0.1 mm cross-section (400 x 1 x 1 cells, 0.05 mm), seeded in the
first 2 mm. The solver sets `deltaT` and `endTime` to 1 and solves once; five
probes sample `activationTime`, from which `setup/extract_cv.py` computes the
conduction velocity (CV).

- **Exercises:** `eikonalSolver` with `eikonalAdvectionDiffusionApproach`
  and a stimulus-box seed on a one-cell-thick cable, and `probes` sampling of
  `activationTime`.
- **Use it for:** tuning `conductivity` or `c0` to a target CV; comparing
  the eikonal CV across mesh resolutions (`cable1DEikonalCV`, `dx` axis) and
  against the monodomain cable (`../monodomain1DCableCV`).

## Regression

Tier: `none`.

None. The eikonal solver with the advection-diffusion form and a
stimulus-box seed is checked by
[`manufacturedSolutions/eikonal/cube`](../../../manufacturedSolutions/eikonal/cube/README.md)
(standard) against an exact solution. What this case adds is the CV
measurement protocol.

## Folder structure

```text
tutorials/electrophysiologyProtocols/cableProtocol/eikonal1DCableCV/
├── 0/
│   └── activationTime
├── constant/
│   ├── electroProperties
│   └── physicsProperties
├── system/
│   ├── blockMeshDict
│   ├── cableProbes
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   └── fvSolution
├── setup/
│   ├── extract_cv.py
│   └── table_summary.py
├── Allrun
├── Allclean
└── README.md
```

## Mesh and pacing setup

- Cable length: 20 mm
- Cross-section: 0.1 mm x 0.1 mm
- Resolution: 0.05 mm along the cable (`400 x 1 x 1`)
- Stimulus region: first 2 mm of the cable, over the whole cross-section

Five probes are placed at:

- 2 mm
- 5 mm
- 10 mm
- 15 mm
- 18 mm

The post-processing script reports:

- activation time at each probe
- segment-wise CV between consecutive probes
- a central calibration CV between 5 mm and 15 mm

The central 5-15 mm value is usually the one you want for conductivity tuning,
because it avoids both the stimulus source region and the distal boundary.

## Parallel execution

`Allrun` accepts a `parallel` argument and uses:

- `decomposePar`
- `runParallel cardiacFoam`
- `reconstructPar`
- `redistributeRestartState`

The decomposition is defined in `system/decomposeParDict`.

## Run

```bash
./Allrun
./Allrun parallel
```

`Allrun` will:

1. generate the mesh with `blockMesh`
2. run `cardiacFoam` in serial or parallel

To extract CV from a manual run, call `python3 setup/extract_cv.py` separately.
The script is under `setup/` for post-processing workflows.

## omnidriver entry

This case is the `cable1DEikonalCV` record; its `dx` axis (metres) rewrites the
block mesh:

```bash
omnidriver --plugin cardiacfoam --repo . plan --strict \
    --entry cable1DEikonalCV --scratch-dir <scratch>
```

## Main calibration knob

The calibration workflow is:

1. run the cable
2. measure CV
3. adjust `conductivity` (or `c0`) in `constant/electroProperties`
4. rerun until the central-cable CV matches the target physiology

The primary CV tuning parameter in this case is:

```cpp
eikonalSolverCoeffs
{
    conductivity  ...
}
```

Increasing `conductivity` increases CV; decreasing it lowers CV.

`c0`, `chi`, `cm` and the spatial resolution also affect the measured value,
so keep those fixed while calibrating conductivity. The time step does not:
the solver sets `deltaT` and `endTime` to 1 and solves once.

## Outputs

- `activationTime` in time directory `1` (the eikonal solver does not write `Vm`)
- probe traces in `postProcessing/cableProbes/`
- CV summary in `postProcessing/cv_summary.txt`
- `setup/table_summary.py` walks a sweep's own `<output_dir>/cases/case_NNNN/
  postProcessing/*_cv_summary.json` (one per case) and writes
  `cable1DCVConvergence_summary.csv`/`.html` plus dx/dt convergence plots at
  the sweep's own output_dir

## Typical use

For a new ionic model or tissue phenotype:

1. calibrate AP shape and APD in a single-cell case first
2. use this 1D cable to tune conductivity to the target CV
3. move the calibrated setup into larger tissue or organ cases

## Notes on dictionary entries

The conductivity matches the monodomain scaling. The eikonal conduction
velocity is `CV = c0 * sqrt(sigma / (chi * cm))`; with `sigma = 0.1334`,
`chi = 140000` and `cm = 0.01`, `sqrt(sigma / (chi * cm)) = 0.00976`, so
`c0 = 72` gives `CV = 0.702 m/s`. `0/activationTime` is initialised to `-1`
("not yet activated") and overwritten by the solve.
