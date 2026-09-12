# Stewart restitution normalization and recalibration plan

## Objective

Establish a reproducible Stewart 2009 Purkinje reference for:

- APD90 restitution;
- conduction-velocity restitution as a function of measured DI90;
- refractory capture/block behavior; and
- spontaneous activation caused by Stewart automaticity.

Use that evidence to make the one-dimensional graph terminology and equations
consistent. Do not alter the graph's scientific behavior until the reference
measurements and their numerical uncertainty have been recorded.

This plan deliberately separates three boundaries:

1. behavior-preserving naming and provenance work;
2. Python protocol and postprocessing corrections; and
3. C++ scientific-model changes, which require approval after calibration.

## Definitions that all artifacts must use

All times are local to the measurement site unless explicitly labelled as
stimulus times.

| Quantity | Definition | Unit |
| --- | --- | --- |
| `stimulusCouplingInterval` | `t(S2 stimulus) - t(last S1 stimulus)` | s |
| `activationInterval` | `t(S2 activation) - t(last S1 activation)` | s |
| `APD90_S1` | `t(S1 repolarization90) - t(S1 activation)` | s |
| `APD90_S2` | `t(S2 repolarization90) - t(S2 activation)` | s |
| `requestedDI90` | requested delay from measured S1 repolarization90 to S2 stimulus | s |
| `measuredDI90` | `t(S2 activation) - t(S1 repolarization90)` | s |
| `stimulusToActivationLatency` | activation time minus its associated stimulus time | s |
| `automaticityCycleLength` | activation interval for an unforced spontaneous beat | s |
| `CV_S2` | probe separation divided by S2 activation-time difference | m/s |

Do not use the unqualified names `DI`, `APD`, `S2 interval`, `ERP`, or
`beatInterval` in new data products. `ERP` is a protocol result (the capture
boundary), not the sum of two constants.

## Current evidence and blockers

The retained calibration document reports five S1 beats at 1 s BCL, an APD of
approximately 0.450 s after conditioning, successful DI90 points from 0.330 to
0.700 s, and automatic activation near the late end of the protocol. The raw
branch traces and `run_smart_restitution.sh` are not retained, so those numbers
cannot currently be regenerated from repository artifacts.

The current graph instead computes:

```text
effectiveDI = activationInterval - apdNominal
minActivationInterval = apdNominal + tableMinimum
```

with `apdNominal = 0.290 s`. The embedded voltage template has fixed morphology
and is independent of this value. The current graph therefore does not yet
implement measured DI90.

The generic single-cell restitution workflow is not a faithful replacement for
the original Stewart calibration:

- its default model is Bueno-Orovio;
- the tutorial-local JSON currently selects TWorld;
- its normal default conditioning BCL is 2 s, not 1 s;
- it applies a coupling interval, not an adaptively scheduled true DI90; and
- automatic beats are not represented as a first-class outcome.

## Phase 0: prerequisites and frozen baseline

### 0.1 Environment

Before launching any case:

1. Source the intended OpenFOAM/cardiacFoam environment.
2. Record the OpenFOAM version, compiler, MPI implementation, git commit, host,
   and build options in a run manifest.
3. Build `cardiacFoam` and the conduction-system library from the same tree.
4. Confirm that `cardiacFoam`, `blockMesh`, `decomposePar`, `reconstructPar`,
   and `mpirun` resolve from that environment.

`mpirun` is visible in the current shell, but `cardiacFoam` is not. No scientific
run should be claimed from the current unsourced shell.

### 0.2 Preserve the existing implementation

Before code changes:

1. Run the smallest existing graph case and save its activation, voltage, DI,
   block-count, and wavebreak fields.
2. Run the current cable case once in serial and once with the intended MPI
   decomposition.
3. Save dictionaries and probe outputs with hashes.
4. Record the current hard-coded table and template as baseline artifacts.

This establishes whether later naming-only changes are behavior preserving.

## Phase 1: make the Stewart protocol reproducible

Create a dedicated Stewart calibration workflow rather than changing the
generic restitution defaults.

### 1.1 Conditioning run

Use the calibration protocol described by the retained audit:

- ionic model: `Stewart`;
- tissue: `myocyte`;
- five applied S1 stimuli;
- S1 BCL: 1.000 s;
- existing cable geometry, conductivity, stimulus region, probe locations,
  spatial discretization, time step, ODE solver, and tolerances initially
  unchanged.

The workflow must write the exact applied stimulus schedule. Do not infer the
number of pulses from `nstim1` because the legacy stimulus implementation uses
inclusive pulse indices (`k <= nstim1`). Prefer explicit stimulus-time lists for
the cable calibration.

Save a restart after the last S1 activation but before its complete
repolarization, together with sufficient probe output to measure the final S1
APD90. The restart time is an implementation detail; it must not be used as a
surrogate for repolarization time.

### 1.2 Measure the conditioned reference

From the unbranched conditioning trace, calculate and retain:

- every S1 stimulus and activation time;
- APD90, APD70, and APD50 for each S1 beat;
- stimulus-to-activation latency at each probe;
- S1 conduction velocity at each interior segment; and
- the final S1 repolarization times at each probe.

Export the final conditioned action-potential waveform as a versioned numeric
artifact. This waveform, rather than an undocumented resting beat, is the
candidate graph voltage template.

### 1.3 Automaticity control

Branch from the conditioned state with no S2 stimulus and run until at least one
unforced activation occurs or until a documented upper time bound is reached.
Record:

- first spontaneous activation time and location;
- local automaticity cycle length;
- propagation of that spontaneous beat through the cable; and
- whether more than one site begins an independent spontaneous wave.

This measurement determines the Stewart automaticity boundary. Do not assume
that the graph's current `escapeInterval = 1.1 s` is exact.

### 1.4 True-DI90 S2 branches

For each requested DI90, schedule the S2 stimulus from the measured final-S1
repolarization90 time at the stimulus/proximal reference site:

```text
tS2_stimulus = tRepolarization90_S1_reference + requestedDI90
```

Use the retained successful points (`0.330, 0.350, 0.400, 0.450, 0.500,
0.700 s`) as reproduction targets, not assumed truth. Add points immediately
below and above the observed capture boundary, and points approaching the
measured automaticity boundary. Choose the final refinement spacing only after
the first coarse sweep.

Use one applied S2 stimulus per branch. A later spontaneous beat is an observed
outcome, not a second S2.

Each branch must end in exactly one status:

- `captured_and_propagated`;
- `local_capture_propagation_block`;
- `stimulus_no_capture`;
- `automaticity_before_s2`;
- `competing_spontaneous_wave`; or
- `analysis_failure`.

The late automaticity cases are censored restitution points. They are not
refractory failures.

## Phase 2: Python changes

### 2.1 Shared event representation

Introduce one tested representation for:

- stimulus events (`S1`, `S2`);
- detected activations;
- repolarization events at 90/70/50%;
- spontaneous activations; and
- case outcome.

Every summary row must include the case ID, model, tissue, probe, requested
DI90, stimulus coupling interval, activation interval, measured DI90, S1/S2
APD90, activation latency, capture status, automaticity status, and CV where
applicable. Failed and censored cases must remain in the CSV/JSON output.

### 2.2 Correct beat association

Refactor `postProcessing_restCurves.py` so that:

1. the duplicate `detect_beats` implementation is removed;
2. threshold crossings remain linearly interpolated;
3. the applied stimulus schedule is read from case metadata, not reconstructed
   only from filenames;
4. each activation is associated with a stimulus only inside a documented
   latency window;
5. an activation without a matching stimulus is labelled spontaneous;
6. a spontaneous beat before S2 produces `automaticity_before_s2`; and
7. S1 and S2 are selected by event identity, not merely as adjacent detected
   beats.

Preserve APD90, APD70, and APD50 outputs, but make the repolarization definition
and baseline/peak construction explicit in the summary metadata.

### 2.3 Correct cable CV extraction

Refactor `postProcessing_cableRestitution.py` and retire the assumptions in
`extract_cv_s2.py` that the last crossing is necessarily S2. The extractor must:

- associate crossings with the labelled S2 event;
- distinguish local capture from propagation to the distal probes;
- interpolate crossing times consistently at every probe;
- calculate segment CV and central-cable CV;
- calculate measured DI90 locally where voltage data permit; and
- emit a structured status instead of treating expected physiological block as
  a successful process with only a printed message.

### 2.4 Calibration orchestration

Add a dedicated Python orchestration entry point that:

1. creates the S1 checkpoint once;
2. obtains the measured reference repolarization time;
3. creates one restart branch per requested DI90;
4. writes explicit absolute stimulus times;
5. runs branches serially or through the existing MPI case runner;
6. invokes postprocessing; and
7. writes a combined provenance manifest, long-form CSV, and JSON summary.

Do not hard-code a presumed `APD90 = 0.450 s` into branch scheduling.

### 2.5 Python tests

Add synthetic trace tests for:

- interpolated activation and repolarization times;
- true DI90 calculation;
- normal S2 capture;
- stimulus artifact without capture;
- local capture followed by distal block;
- spontaneous activation before S2;
- spontaneous activation after a valid S2;
- competing spontaneous and stimulated waves;
- missing/truncated probe data; and
- preservation of failed/censored rows in output.

Add workflow tests proving that requested DI90 is converted into the correct
absolute S2 time and that generated cases contain exactly five applied S1
stimuli and one S2 stimulus.

## Phase 3: numerical checks before fitting the graph

Run a bounded convergence/reproducibility study at:

- one fully recovered point;
- one point on the steep CV region;
- one point on each side of the capture boundary; and
- one case near the automaticity boundary.

For those cases compare:

- baseline and refined time steps;
- baseline and refined axial mesh resolution; and
- serial versus the intended MPI rank count.

Use the observed variation to set acceptance tolerances. Do not select APD, CV,
or activation-time tolerances before this check.

The calibrated table must contain measured values and uncertainty/provenance,
not rounded values copied manually from console output.

## Phase 4: C++ normalization (behavior preserving)

This batch may proceed before the scientific-model change, provided baseline
fields remain identical.

1. Rename local variables `beatInterval` to `activationInterval`.
2. Describe the existing `DI_` calculation as `effectiveDI`, not measured DI90.
3. Describe `minBeatInterval_` as `minimumActivationInterval_`, not ERP.
4. Describe `escapeInterval_` as a timer-based automaticity surrogate.
5. Rename internal table symbols from generic `purkinjeDI`/`purkinjeCV` to
   Stewart-specific, provenance-bearing names without changing values.
6. Keep existing dictionary keys and diagnostic field names as compatibility
   aliases during this batch.
7. Mark first-beat/fully-recovered DI as undefined in diagnostics rather than
   silently presenting the table maximum as a measured DI.
8. Correct README, architecture, template comments, and dictionary-catalog
   claims to match executable behavior.

Acceptance: graph activation times, voltage, block counters, and CV must match
the frozen baseline for the same input, apart from an explicitly added
undefined/fully-recovered diagnostic.

## Phase 5: C++ scientific-model correction

Do this only after reviewing the new Stewart evidence. The intended target is a
true-DI90 graph model.

### 5.1 Stewart restitution profile

Create one Stewart profile containing, with provenance:

- the measured DI90-to-CV table;
- conditioned voltage-template samples;
- the template's measured APD90;
- the measured capture boundary; and
- the measured automaticity cycle length.

The profile must state the conditioning BCL, number of S1 beats, probe used for
DI, cable conductivity, mesh/time-step resolution, and source artifact hash.

### 5.2 Per-node recovery state

Replace the surrogate subtraction with explicit state:

```text
lastRepolarization90Time[node]
measuredDI90 = candidateActivationTime - lastRepolarization90Time[node]
```

For the fixed-template graph, repolarization90 time can be derived from the
accepted activation time plus the calibrated template APD90. This remains a
fixed-APD surrogate; it must not be called APD restitution.

### 5.3 Propagation semantics

Use a documented two-part rule:

1. evaluate outgoing edge CV from the activated source node's measured DI90;
2. test target recovery at the candidate arrival time, not at source departure
   time.

Conceptually:

```text
sourceDI90 = sourceActivation - sourcePreviousRepolarization90
candidateArrival = sourceActivation + edgeLength / CV(sourceDI90)
targetDI90 = candidateArrival - targetPreviousRepolarization90
accept target only if targetDI90 >= calibratedMinimumDI90
```

This removes the current premature block where a target is refractory at
departure but recovered by arrival. It also avoids an implicit circular solve
for target-DI-dependent CV. If a different edge law is desired, that is a
separate model decision.

### 5.4 Automaticity semantics

Retain automaticity as an explicit surrogate, but rename the setting to
`automaticityCycleLength`. Keep `escapeInterval` as a deprecated dictionary
alias for existing cases.

Test whether per-node automaticity produces simultaneous/competing sources.
The no-S2 cable control determines whether the graph needs an explicit
pacemaker-origin policy. Do not introduce that policy without evidence.

### 5.5 Diagnostics and compatibility

Add explicit fields/status for:

- `measuredDI90` for reactivated nodes;
- `activationInterval`;
- `fullyRecovered` or `DI90Defined`;
- stimulated, propagated, and automatic activation origin;
- recovery block; and
- automaticity competition.

Keep the old `DI` field for one compatibility window only if it is clearly
labelled as the legacy effective-DI quantity. Do not silently change the
meaning of an existing output field.

### 5.6 C++ tests

Add focused tests for:

- table interpolation and endpoint clamping;
- first activation/fully recovered behavior;
- true DI90 after a second activation;
- exact capture-boundary behavior;
- target recovery during edge transit;
- automatic activation timing;
- collision of an incoming wave with a scheduled automatic event;
- restart equivalence; and
- deterministic serial/MPI graph results after distributed-graph work lands.

## Phase 6: end-to-end acceptance

Compare the graph against the monodomain Stewart reference for the same protocol.
The final report must include:

1. APD90 and waveform comparison for the conditioned S1 template;
2. measured DI90 versus CV, including numerical uncertainty;
3. capture/block boundary;
4. automaticity-censored boundary;
5. graph versus monodomain activation-time and CV residuals;
6. behavior at DI values between table samples; and
7. serial/MPI consistency.

The existing hard-coded numbers may be retained only if the rerun reproduces
them within tolerances derived from the numerical checks. Otherwise update the
profile and calibration document together from the generated artifacts.

## Proposed implementation batches

| Batch | Scope | Scientific behavior change? |
| --- | --- | --- |
| A | Python event types, synthetic tests, explicit statuses | No |
| B | Stewart checkpoint/branch orchestrator and provenance outputs | No model change |
| C | Reproduce Stewart single-cell/cable data and numerical checks | Experiment only |
| D | C++ terminology and documentation cleanup | No, except new diagnostics |
| E | Review calibration and approve the target graph equations | Decision gate |
| F | True-DI90 recovery, arrival-time gating, Stewart profile | Yes |
| G | End-to-end graph/monodomain validation | Validation |

Do not combine batches D and F. A behavior-preserving terminology diff must be
reviewable independently from the change in restitution and recovery equations.

