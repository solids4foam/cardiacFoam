# Stewart restitution C++ reconciliation plan

## Scope and evidence boundary

This plan reconciles `restitutionEikonalSolver1D` with the reproducible Stewart
cable measurements without silently changing the graph model. The boundary
refinement sweep must finish before any scientific constant, table, gating
equation, voltage template, or public diagnostic field is changed.

The reproducible reference uses Stewart/myocyte, five S1 stimuli at 1.000 s,
`deltaT = 1e-6 s`, `dx = 0.1 mm`, conductivity 2.3 S/m, and six MPI ranks.
The applied schedule is stored in `.cardiacfoam_protocol.json`; activations and
repolarizations are interpolated from the cable probes rather than inferred
from case names.

## Comparison with the retained calibration

| Quantity | Retained calibration | Reproducible measurement | Consequence |
| --- | ---: | ---: | --- |
| Conditioned proximal APD90 | approximately 450 ms | 303.037 ms | The retained APD claim is not reproduced. |
| Conditioned APD90 across probes | not retained | 303.037--309.987 ms | `apdNominal=0.290` is a surrogate, not the measured conditioned APD90. |
| First spontaneous proximal activation | described near the late protocol end | 5.203520 s | The event is explicit and unforced. |
| Proximal automaticity cycle after final S1 activation | implicitly about 1.1 s | 1.202470 s | `escapeInterval=1.1` is not the measured cycle. |
| Spatial automaticity pattern | not retained | 0.222 ms crossing spread over 16 mm | Stewart reactivation is near-synchronous across the cable, not a normal end-launched wave. |
| CV at nominal DI 330 ms | 2.03 m/s | 3.357960 m/s at proximal measured DI90 331.071 ms | Old and new axes/results are not interchangeable. |
| CV at nominal DI 350 ms | 2.43 m/s | 3.363985 m/s at proximal measured DI90 351.065 ms | Same mismatch. |
| CV at nominal DI 400 ms | 2.95 m/s | 3.366801 m/s at proximal measured DI90 401.054 ms | Same mismatch. |
| CV at nominal DI 450 ms | 3.18 m/s | 3.356133 m/s at proximal measured DI90 451.047 ms | Same mismatch. |
| CV at nominal DI 500 ms | 3.29 m/s | 3.335530 m/s at proximal measured DI90 501.043 ms | Values become close only by coincidence. |
| CV at nominal DI 700 ms | 3.37 m/s | 3.184866 m/s at proximal measured DI90 701.052 ms | The reproducible curve declines during late diastolic depolarization. |

The retained table rises steeply from 2.03 to 3.37 m/s. The reproducible
true-DI90 curve is approximately flat near 3.36 m/s from 330--450 ms, then
declines to 3.18 m/s at 700 ms. The retained branch traces and branching script
are absent, and the available 4.25 s field checkpoint does not contain the
Stewart internal state arrays, so its claimed restart equivalence cannot be
verified.

### Concrete label mismatch in the retained 0.500 s case

The fine-step retained case schedules the final S1 stimulus at 4.0 s and S2 at
4.5 s. At the proximal probe it measures:

```text
stimulusCouplingInterval = 4.500000000 - 4.000000000 = 0.500000000 s
activationInterval       = 4.501140255 - 4.001049728 = 0.500090527 s
measuredDI90             = 4.501140255 - 4.304086261 = 0.197053994 s
legacyGraphEffectiveDI   = 0.500090527 - 0.290000000 = 0.210090527 s
```

These four quantities are not interchangeable. The retained table abscissae
0.330--0.700 have the same numeric form as the stimulus coupling intervals used
by the old workflow. The original script/raw branches are missing, so this is
an inference rather than proof, but the available schedule and trace strongly
indicate that the old axis was a stimulus coupling interval mislabeled as DI,
not measured DI90. The boundary sweep's true-DI90=200 ms case is the direct
comparison for the legacy coupling=500 ms case.

## Agreed working interpretation

Until the boundary sweep either confirms or rejects it, use the following
interpretation consistently:

1. The retained CV ordinates `{2.03, 2.43, 2.95, 3.18, 3.29, 3.37}` are
   provisionally credible measurements and must not be discarded merely
   because their axis label was wrong.
2. Their retained abscissae `{0.330, 0.350, 0.400, 0.450, 0.500, 0.700}` are
   most likely S1--S2 stimulus coupling intervals, not measured DI90.
3. The retained approximately 450 ms conditioned APD claim is unsupported and
   conflicts with the reproducible 303--310 ms APD90 measurement.
4. `apdNominal=0.290 s` is close to, but not equal to, the measured proximal
   conditioned APD90. It is currently a fixed duration surrogate in the graph
   equation, not an ERP and not a recorded repolarization event.
5. The new true-DI90=330--700 ms results do not invalidate the old CV values,
   because they sample a later recovery range than the old coupling-labelled
   branches.

Using proximal APD90 `0.303036533 s` and approximately 1 ms S2 latency gives
this provisional axis reinterpretation:

| Retained coupling label | Approximate measured DI90 | Retained CV |
| ---: | ---: | ---: |
| 0.330 s | 0.027 s | 2.03 m/s |
| 0.350 s | 0.047 s | 2.43 m/s |
| 0.400 s | 0.097 s | 2.95 m/s |
| 0.450 s | 0.147 s | 3.18 m/s |
| 0.500 s | 0.197 s | 3.29 m/s |
| 0.700 s | 0.397 s | 3.37 m/s |

These shifted abscissae are validation targets, not final constants. Replace
the approximation with the measured local activation and repolarization times
from each recovered branch before generating a C++ table. In particular, do
not implement the conversion by subtracting one APD from every filename or
case label at runtime.

The compatibility and true-DI90 paths are intentionally separate:

- **Legacy-compatible correction:** preserve the old CV ordinates, correct
  their provenance/axis, and use an explicitly named fixed recovery-duration
  approximation.
- **True-DI90 model:** store/derive repolarization90 per node and query a table
  whose abscissae are measured DI90 values.

Both paths require the existing behavior to remain available until regression
and graph-level comparisons are complete.

## What can change immediately without changing behavior

### Batch A: characterization tests

Freeze the current graph behavior before editing names:

1. Test the six legacy table points and both clamped extrapolation branches.
2. Test the current minimum activation interval:
   `0.290 + 0.330 = 0.620 s`.
3. Test activation acceptance immediately below, at, and above 0.620 s.
4. Test spontaneous scheduling at the current 1.1 s escape interval.
5. Test the current public diagnostic names and values: `DI`, `minDI`,
   `blockCount`, and `wavebreakCount`.
6. Save one small serial graph baseline and verify its activation, voltage,
   block-count, and wavebreak fields byte-for-byte or numerically as
   appropriate.

### Batch B: internal semantic names

Make a behavior-preserving rename only after Batch A passes:

- `DI_` -> `effectiveDI_`;
- `minDI_` -> `minEffectiveDI_`;
- `minBeatInterval_` -> `minimumActivationInterval_`;
- local `DIact` -> `effectiveDIAtActivation`;
- local `DIj` -> `targetEffectiveDI`.

Define the existing quantity explicitly in comments:

```text
effectiveDI = activationInterval - apdNominal
```

It is not measured DI90 because the solver stores activation time but no
repolarization time. Preserve the dictionary key `apdNominal` and public fields
`DI`/`minDI` during this batch. Renaming those external contracts requires a
separate compatibility design.

### Batch C: documentation and catalog correction

Without changing defaults:

1. Mark the retained CV table as a legacy coupling-interval calibration whose
   raw traces are unavailable and whose CV ordinates remain provisionally
   usable after axis reconciliation.
2. Stop describing `apdNominal` as a constant ERP. It is the fixed duration
   subtracted from activation interval by the current surrogate equation.
3. Stop describing `apdNominal + tableMinimum` as a measured ERP. It is the
   current graph's minimum accepted activation interval.
4. State that `escapeInterval` is a graph scheduling parameter, not yet the
   measured Stewart automaticity cycle.
5. Correct the restart documentation: Vm/current fields alone do not preserve
   Stewart gating and concentration state.
6. Apply the same wording to the driverFOAM dictionary catalog and template
   comments while retaining keys and defaults.

## Scientific changes that must wait for the refinement sweep

Do not make the following changes merely from the six completed points:

- replace `purkinjeDI` or `purkinjeCV`;
- change `purkinjeAPDnominal`;
- change `purkinjeEscapeInterval`;
- use the first table abscissa as the capture boundary;
- rename/remove public `DI` or `minDI` output fields;
- replace the voltage template;
- claim an ERP from APD plus a table minimum.

Also do not erase or overwrite the retained CV ordinates. The refinement sweep
must first test whether the shifted-axis interpretation recovers them.

The 0/100/200/300 ms cases determine early capture behavior. The
850/900/950 ms cases distinguish valid late S2 capture from automaticity that
pre-empts S2. Expected physiological block and automaticity are outcome data,
not process failures.

## Planned scientific C++ design after calibration

### Batch D: separate model concepts

Introduce distinct state/configuration for:

- `lastActivationTime`;
- `lastRepolarization90Time`;
- `measuredDI90 = candidateActivationTime - lastRepolarization90Time`;
- `minimumMeasuredDI90` or another explicitly calibrated capture boundary;
- `automaticityCycleLength`;
- the CV-versus-measured-DI90 table.

The capture boundary must not be inferred automatically from the first CV table
abscissa. The table domain and the refractory boundary are separate calibrated
objects.

Implement the legacy-compatible fixed-duration path and the true-DI90 path as
distinct modes. A change to the meaning of the existing default mode must not
occur implicitly through a renamed variable or a replaced array.

### Batch E: event propagation rule

For an accepted source activation:

```text
sourceMeasuredDI90 = sourceActivation - sourcePreviousRepolarization90
edgeCV = CV(sourceMeasuredDI90)
candidateArrival = sourceActivation + edgeLength/edgeCV
targetMeasuredDI90 = candidateArrival - targetPreviousRepolarization90
```

Accept the target only when its independently calibrated recovery rule passes.
Record whether rejection came from recovery, zero conductance, or event
competition. This requires a scientific decision on whether the calibrated CV
belongs to the source site, target site, or a segment reference; the cable
summary must make that choice explicit.

### Batch F: APD/repolarization representation

Choose one of two explicit models:

1. fixed-template APD90 derived from the selected voltage template; or
2. an APD90-restitution model calibrated from the cable/single-cell protocol.

Do not retain a constant called `apdNominal` while claiming true DI90 unless it
actually defines the repolarization event used by the algorithm. Preserve the
legacy mode as a selectable compatibility path for existing cases.

### Batch G: automaticity

Use the no-S2 measurements to decide whether automaticity is:

- scheduled at every graph node, matching the near-synchronous Stewart cable;
- restricted to configured pacemaker nodes; or
- represented by a spatially varying cycle.

Define deterministic priority when automaticity and an external/stimulated
event occur in the same time window. Add `automaticity_before_s2` and competing
event tests before changing the current 1.1 s default.

### Batch H: voltage template provenance

Export a conditioned Stewart waveform with its activation and APD50/70/90
metadata. Version the numeric artifact and generate the C++ array from it.
Keep waveform morphology separate from recovery/CV parameters; replacing the
array is a scientific change even if event timing is unchanged.

### Batch I: complete ionic restart state

The Stewart wrapper initializes `STATES_`, `ALGEBRAIC_`, and `RATES_` in its
constructor, while the current 4.25 s case directory contains only Vm, Vm_0,
activation time, and current fields. Implement restart serialization/import for
all state variables needed to continue the ODE exactly. Prove a continuous run
and a restart branch agree before driverFOAM reuses conditioning checkpoints.

### Batch J: distributed graph execution

Integrate the true-DI90 state only after the graph ownership/halo/event-queue
design in `DISTRIBUTED_1D_GRAPH_PARALLELIZATION_PLAN.md` is validated. MPI
acceptance must compare serial and distributed activation times, recovery
state, event sources, block counts, wavebreak counts, and graph voltage for
normal capture, block, automaticity, and simultaneous competing events.

## Acceptance gate for changing defaults

Changing default scientific behavior requires all of the following:

1. boundary refinement completed with structured outcomes;
2. one combined provenance-bearing CSV/JSON table;
3. chosen reference probe/segment documented;
4. time-step and mesh uncertainty stated;
5. serial/MPI equivalence for the calibration case;
6. legacy graph behavior preserved behind an explicit compatibility mode;
7. true-DI90 unit and graph integration tests passing;
8. tutorial, template, catalog, and C++ comments updated together.
