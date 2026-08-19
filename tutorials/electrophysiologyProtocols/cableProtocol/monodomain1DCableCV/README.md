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
tutorials/coreProtocols/cableProtocol/monodomain1DCableCV/
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

## Driver restitution entry

This case is also exposed as a newly normalized S1-S2 restitution sweep:

```bash
foamctl run --entry cableRestitutionCurves --strict
```

Or you can sweep custom restitution intervals using a JSON config:

```bash
foamctl sweep-run --spec sweep_restitution.json --output-dir validation_run
```

The sweep logic and default S2 pacing intervals are fully centralized in:

```text
applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/defaults/cable_restitution_curves.py
```

*Note:* Because multiple wavefronts are generated, the normalized post-processing step (`setup/postProcessing_cableRestitution.py`) parses the raw voltage traces to isolate the CV of the second (S2) wavefront.

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

## Restitution calibration for restitutionEikonal solver

This case has been calibrated to generate a functional S1-S2 restitution curve for use in dynamic restitution-aware 3D eikonal sweeps (`restitutionEikonalSolver1D`).

### Initial calibration (Resting CV)

The `conductivity` in `constant/electroProperties` was tuned to achieve a baseline conduction velocity of **3.0 m/s**. Using the Stewart ionic model, a conductivity of **2.3 S/m** yields a CV of exactly **3.03 m/s** when paced from a mathematically perfect resting state (the very first beat at t=0).

### S1-S2 protocol automation

A smart branching approach leverages OpenFOAM's native restart capabilities to generate the restitution curve efficiently:

1. **Phase 1 (S1 Drive Train):** 5 S1 beats at a Basic Cycle Length (BCL) of 1000 ms, simulated from t=0 to t=4.25 s across 6 cores. The OpenFOAM field state was saved at t=4.25 s.
2. **Phase 2 (S2 Branches):** For each Diastolic Interval (DI) tested, the t=4.25 s checkpoint was restored and the simulation resumed, injecting the S2 premature beat and simulating only the brief ~25 ms window required for wave propagation.

### Electrophysiological phenomena observed

#### Supernormal conduction (velocity peaking)

When pacing at 1 Hz, the resting membrane potential ($V_m$) does not perfectly return to its absolute minimum before the next beat arrives due to ionic memory (e.g., slight extracellular $K^+$ accumulation). Because $V_m$ sits slightly higher (less negative), the membrane is closer to the excitation threshold, requiring less depolarizing current to trigger adjacent cells. This results in **supernormal conduction**:

- **Beat 1:** 3.03 m/s
- **Beat 2:** 3.22 m/s
- **Beat 5:** 3.17 m/s
- **S2 (DI = 0.700 s):** 3.37 m/s

#### APD prolongation and ERP shift

The nominal resting Action Potential Duration (APD) of the Stewart model is approximately 290 ms. Pacing 5 times at 1.0 Hz caused the APD to physiologically lengthen to approximately 450 ms. Because the APD prolonged, the Effective Refractory Period (ERP) pushed significantly outward. Any premature S2 beats with a DI below 0.330 s fell inside the Absolute Refractory Period and naturally failed to propagate.

### Final restitution curve

By sweeping S2 intervals that successfully propagated outside the ERP, the steep gradient of the restitution curve was isolated. These functional values are hardcoded into `src/electroModels/conductionSystemModels/restitutionEikonalSolver1D/restitutionTemplates.H`:

**Diastolic Intervals (s):**

```
{ 0.330, 0.350, 0.400, 0.450, 0.500, 0.700 }
```

**Conduction Velocities (m/s):**

```
{ 2.03,  2.43,  2.95,  3.18,  3.29,  3.37 }
```

These values are actively used by the eikonal solver during 3D Purkinje sweeps.

## Typical use

For a new ionic model or tissue phenotype:

1. calibrate AP shape and APD in a single-cell case first
2. use this 1D cable to tune conductivity to the target CV
3. if using dynamic restitution in eikonal solves, run an S1-S2 protocol to extract functional restitution values
4. move the calibrated setup into larger tissue or organ cases
