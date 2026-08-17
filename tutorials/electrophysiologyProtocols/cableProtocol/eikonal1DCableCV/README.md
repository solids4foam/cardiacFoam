# monodomain1DCableCV tutorial

This tutorial is the minimal tissue-scale conduction-velocity calibration case
for `cardiacFoam`.

- Electro model: `monodomainSolver`
- Geometry: 20 mm cable with a single-cell cross-section
- Purpose: launch one paced wave, measure activation times at several stations,
  and compute a conduction velocity suitable for conductivity calibration

## Why this case exists

Before a 2D sheet, 3D slab, ventricle, or Purkinje-coupled simulation can be
trusted, the tissue propagation speed has to be checked in isolation. This case
gives the smallest practical setup for that task:

- one paced wavefront
- no geometric curvature effects
- no Purkinje coupling
- direct access to activation-time-based CV measurements

The intended workflow is:

1. choose an ionic model and tissue type
2. run the cable
3. measure CV
4. adjust `conductivity` in `constant/electroProperties`
5. rerun until the central-cable CV matches the target physiology

## Folder structure

```text
tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/
├── constant/
│   ├── electroProperties
│   └── physicsProperties
├── system/
│   ├── blockMeshDict
│   ├── cableProbes
│   ├── controlDict
│   ├── fvSchemes
│   └── fvSolution
├── setup/
│   └── extract_cv.py
├── Allrun
├── Allclean
└── README.md
```

## Mesh and pacing setup

- Cable length: 20 mm
- Cross-section: 0.1 mm x 0.1 mm
- Resolution: 0.1 mm along the cable (`200 x 1 x 1`)
- Stimulus region: first 0.5 mm of the cable

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

The decomposition is defined in `system/decomposeParDict`.

## Run

```bash
./Allrun
./Allrun parallel
```

`Allrun` will:

1. generate the mesh with `blockMesh`
2. run `cardiacFoam` in serial or parallel
3. compute the CV summary from the probe outputs

## Driver convergence entry

This case is also exposed as a registered driver sweep:

```bash
foamctl all --entry cable1DCVConvergence
```

The default sweep config is stored in:

```text
setup/driver_config.json
```

That workflow mutates `blockMeshDict`, `controlDict`, and `electroProperties`
per sweep case, runs in parallel by default, exports one CV summary per case,
and then writes per-ionic-model convergence CSVs plus convergence plots under
model-specific output folders such as `outputsCVConvergence/BuenoOrovio/` or
`outputsCVConvergence/TWorld/`.

## Main calibration knob

The primary CV tuning parameter in this case is:

```cpp
monodomainSolverCoeffs
{
    conductivity  ...
}
```

Increasing `conductivity` increases CV; decreasing it lowers CV.

`chi`, `cm`, spatial resolution, and time step also affect the measured value,
so keep those fixed while calibrating conductivity.

## Outputs

- field writes in time directories (`Vm`, `activationTime`, ...)
- probe traces in `postProcessing/cableProbes/`
- CV summary in `postProcessing/cv_summary.txt`
- sweep artifacts grouped by ionic model inside `outputsCVConvergence/<ionicModel>/`

## Typical use

For a new ionic model or tissue phenotype:

1. calibrate AP shape and APD in a single-cell case first
2. use this 1D cable to tune conductivity to the target CV
3. move the calibrated setup into larger tissue or organ cases
