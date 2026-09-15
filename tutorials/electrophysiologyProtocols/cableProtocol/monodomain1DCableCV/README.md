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
- a straight cable geometry
- conduction confined to the cable myocardium
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
- Resolution: 0.2 mm along the cable (`100 x 1 x 1`)
  - The frozen reference protocol below overrides this to 0.1 mm; the
    committed default is the coarser mesh.
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
driverFoam run --strict --entry cable1DCVConvergence
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
driverFoam run --entry cable1DRestitution --strict
```

Or you can sweep custom restitution intervals using a JSON config:

```bash
driverFoam sweep-run \
    --spec tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/sweep.json \
    --output-dir .tmp/driverfoam/cable-restitution
```

The sweep logic and default S2 pacing intervals are fully centralized in
the driverFOAM add-on's `cable_1d_restitution` cardiacFoam plugin defaults.

*Note:* Because multiple wavefronts are generated, the normalized post-processing step (`setup/postProcessing_cableRestitution.py`) parses the raw voltage traces to isolate the CV of the second (S2) wavefront.

### Frozen Stewart true-DI90 reference

`sweep_stewart_true_di90_dt1e-6.json` is runnable as the frozen
Stewart/myocyte reference protocol: five S1 stimuli at a 1 s BCL, one S2, a
0.1 mm cable resolution, `deltaT = 1e-6 s`, and the configured conductivity.
Its `reference_repolarization90_s = 4.304086260869566` was measured at the
proximal reference probe for exactly that conditioned setup.  The driver uses
it to schedule `t(S2) = t(repolarization90_S1) + requestedDI90`.

Do not reuse that timestamp after changing the ionic model, tissue,
conductivity, stimulus protocol, mesh, time step, solver, or reference probe.
It is not yet regenerated automatically: a future model/protocol calibration
must first run and measure its own conditioning reference before launching its
true-DI90 sweep.  In particular, a future Trovato protocol requires a separate
reference measurement rather than this Stewart value.

### Restitution postprocessing contract

`setup/postProcessing_cableRestitution.py` reads the applied schedule from the
`.cardiacfoam_protocol.json` sidecar the driver writes, associates each
detected activation with a labelled stimulus, and classifies the branch into
exactly one `protocol_outcome`:

| Outcome | Meaning |
| --- | --- |
| `captured_and_propagated` | S2 captured and reached every probe |
| `competing_spontaneous_wave` | S2 propagated, but an unforced beat followed it |
| `local_capture_propagation_block` | S2 captured at some probes only |
| `stimulus_no_capture` | no probe responded to S2 |
| `automaticity_before_s2` | an unforced beat pre-empted the S2 stimulus |
| `automaticity_observed` / `no_spontaneous_activation_observed` | no-S2 control branches |
| `analysis_failure` | the trace could not be analysed; the reason is recorded |

Blocked, censored and failed branches are outcome data, not process failures:
every case writes `<case_id>_event_summary.json`, `<case_id>_events.csv` and
the long-form `<case_id>_restitution.csv` regardless of outcome.

Each restitution row keeps four quantities that are **not** interchangeable —
`requested_di90_s`, `stimulus_coupling_interval_s`, `activation_interval_s` and
the measured `measured_di90_s` — alongside S1/S2 APD50/70/90 and the central
CV. `measured_di90_s` is a per-probe quantity: repolarization90 varies by about
10 ms along this cable, so the DI90 axis is defined at the reference probe
(`is_reference_probe`) while CV is measured across the central segment. Both
choices are recorded in `measurement_settings`.

Stimulus-to-activation association uses a per-probe window,
`latency allowance + probe distance / CV floor`, rather than a flat one. A flat
50 ms window bound the measured automaticity beat at 5.203520 s to an S2
stimulus 49.4 ms earlier and reported a fabricated capture; the derived window
is bounded above by that collision.

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

## Restitution calibration for restitutionEikonalSolver1D

This case is also the measurement apparatus for the constants compiled into
`restitutionEikonalSolver1D`. The full record -- protocol, conditioned state,
the CV and APD90 tables, the capture boundary, and the mesh and time-step
convergence behind them -- is in
[Purkinje_S1_S2_Calibration.md](Purkinje_S1_S2_Calibration.md).

It is deliberately not repeated here. An earlier version of this README carried
a second copy of those numbers, and the two drifted: both asserted a conditioned
APD near 450 ms and a capture boundary at DI 0.330 s, neither of which
reproduces.

In brief, measured on this cable with Stewart/myocyte at 2.3 S/m:

- conditioned APD90 is **303 ms**, and 1 Hz pacing *shortens* it from the
  from-rest value rather than prolonging it;
- conduction velocity runs **1.40 m/s at DI90 62 ms to 3.33 m/s at 401 ms**,
  roughly thirty times steeper at the short end than the long end;
- a premature beat captures at DI90 52.85 ms and fails at 25 ms;
- APD90 varies only 7.7% across the whole capturable range, which is why the
  solver carries a constant `apdNominal` and no APD restitution curve.

The sweep specifications in this directory drive the protocol through the
external orchestration add-on. `requestedDI90` is an input and `measuredDI90`
is a result; the calibration table is indexed on the measured value.

## Typical use

For a new ionic model or tissue phenotype:

1. calibrate AP shape and APD in a single-cell case first
2. use this 1D cable to tune conductivity to the target CV
3. if using dynamic restitution in eikonal solves, run an S1-S2 protocol to extract functional restitution values
4. move the calibrated setup into larger tissue or organ cases
