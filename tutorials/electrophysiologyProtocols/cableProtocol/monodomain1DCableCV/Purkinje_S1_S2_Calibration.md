# Stewart Purkinje S1-S2 calibration

Measurement record for the constants compiled into
`src/electroModels/conductionSystemModels/restitutionEikonalSolver1D/restitutionTemplates.H`.

Measured 2026-09-15. This document **replaces** an earlier calibration note
whose values could not be reproduced; see "What changed, and why" below.

## Protocol

| | |
| --- | --- |
| Ionic model / tissue | Stewart / myocyte |
| Cable | 20 mm, 0.1 x 0.1 mm cross-section |
| Conductivity | 2.3 S/m |
| Drive train | five S1 stimuli at 1000 ms BCL |
| Premature beat | one S2 per branch |
| S2 scheduling | `t(S2) = t(repolarization90 of the last S1) + requestedDI90` |
| Reference repolarization90 | 4.303753683083511 s at the proximal probe |
| Probes | 2, 5, 10, 15, 18 mm |
| Reported CV | central segment, 5 mm to 15 mm |
| Discretisation | dx 0.2 mm, deltaT 1e-5 s, six MPI ranks |

S2 is scheduled from a **measured repolarization time**, so the swept axis is a
diastolic interval rather than a stimulus coupling interval that resembles one.
Stimulus times are emitted as explicit absolute values at twelve significant
figures; at this reference a six-figure format loses 3.74 us, roughly four
steps at the finest time step used here.

Every quantity below is read from the probe traces. `requestedDI90` is an
input; `measuredDI90 = t(S2 activation) - t(S1 repolarization90)` is the
result, and the two differ by the stimulus-to-activation latency of about
1.4 ms. The table is indexed on the measured value.

## Conditioned state

APD90 at the proximal probe across the drive train:

```text
beat 1 (from rest)  315.74 ms
beat 2              301.17 ms
beat 3              301.77 ms
beat 4              302.27 ms
beat 5 (final S1)   302.70 ms
```

Pacing at 1 Hz **shortens** APD90 from its from-rest value and then re-lengthens
it slightly. The conditioned value is 303.037 ms at dx 0.1 mm / deltaT 1e-6 s,
independently reproduced as 302.70 ms at the coarser setting above -- a 0.34 ms
spread, so this quantity is effectively resolution independent.

With no S2 applied, the first unforced activation follows the final S1 by
**1202.470 ms**. That reactivation is near-synchronous across the cable
(0.222 ms of crossing spread over 16 mm) rather than a wave launched from one
end, so it is a pacemaker-like escape rather than a conduction result.

## Restitution

| requested DI90 | measured DI90 | outcome | CV [m/s] | S2 APD90 [ms] |
| ---: | ---: | --- | ---: | ---: |
| 25 ms | -- | `stimulus_no_capture` | -- | -- |
| 50 ms | 52.85 ms | captured, decremental | not defined | 279.50 |
| 60 ms | 62.24 ms | captured and propagated | 1.4042 | 280.71 |
| 70 ms | 71.85 ms | captured and propagated | 1.7999 | 281.63 |
| 80 ms | 81.64 ms | captured and propagated | 2.1521 | 282.39 |
| 90 ms | 91.51 ms | captured and propagated | 2.3423 | 283.15 |
| 100 ms | 101.43 ms | captured and propagated | 2.4992 | 283.96 |
| 200 ms | 201.16 ms | captured and propagated | 3.1542 | 292.04 |
| 400 ms | 401.06 ms | captured and propagated | 3.3280 | 299.10 |

### Conduction velocity

CV falls 44% between DI90 101 ms and 62 ms. Local slopes:

```text
62 -> 101 ms    0.028   (m/s)/ms
101 -> 201 ms   0.0066  (m/s)/ms
201 -> 401 ms   0.00087 (m/s)/ms
```

The short end is roughly thirty times steeper than the long end. Any table
that stops near 100 ms and clamps below it is wrong by a wide margin in exactly
the region where conduction block is decided.

### APD90

Across the whole capturable range APD90 moves 279.5 to 302.7 ms: **23 ms, or
7.7%**, with a maximum slope near 0.09. That is an order of magnitude below the
slope-1 threshold associated with alternans, and it is why the solver uses a
single constant `apdNominal` rather than an APD restitution curve. The constant
is a fixed duration, not a measured repolarization event and not an ERP.

### Capture boundary

Capture fails at DI90 25 ms and succeeds at 52.85 ms, so the boundary lies in
(25, 52.85] ms. `purkinjeMinimumDI90` is set to 50 ms, inside that bracket and
slightly toward the conservative end. Refining it further is what the
seven-case boundary sweep (`sweep_stewart_di90_boundaries_dt1e-6.json`) is for.

This boundary is calibrated **separately from the CV table's domain**. The two
were previously coupled, with the acceptance threshold derived from the table's
lowest abscissa; that is what made the graph refuse premature beats the cable
conducts.

### Where the eikonal assumption fails

At DI90 52.85 ms the S2 beat reaches all five probes, yet segment velocities
run 0.462, 0.943, 2.108 m/s and then invert -- the 18 mm probe activates
0.18 ms *before* the 15 mm probe, against an expected 1.4 ms transit. The wave
is decremental and then recovers, and the far end depolarizes essentially
independently of the arriving front.

A well-defined local wavefront velocity is the eikonal model's core assumption,
and here it does not exist. The postprocessor reports no CV for this branch,
which is correct. This is not a failed measurement; it is a measurement of the
assumption breaking, which is the condition the screening solver exists to
flag.

## Numerical uncertainty

At DI90 ~ 101 ms, holding the other axis fixed:

| dx [mm] | CV [m/s] | | deltaT [s] | CV [m/s] |
| ---: | ---: | --- | ---: | ---: |
| 0.2 | 2.4992 | | 1e-5 | 2.5075 |
| 0.1 | 2.5075 | | 5e-6 | 2.5134 |
| 0.05 | 2.5112 | | 2e-6 | 2.5186 |

Mesh convergence is roughly first order; Richardson extrapolation gives about
2.514 m/s. The tabulated ordinates sit approximately 0.8% below the fully
resolved value. The whole table is recorded at one resolution rather than
mixing corrected and uncorrected points.

APD90 moves 0.26 ms and measured DI90 moves 0.28 ms across a fourfold mesh
refinement, so those quantities are converged well inside their own tolerances.

## What changed, and why

The previous calibration note reported a conditioned APD of approximately
450 ms, an automaticity cycle near 1.1 s, a capture boundary at DI 0.330, and
the CV table `{2.03, 2.43, 2.95, 3.18, 3.29, 3.37}` on abscissae
`{0.330, 0.350, 0.400, 0.450, 0.500, 0.700}`.

None of it reproduces. The conditioned APD is 303 ms, not 450; the automaticity
cycle is 1.202 s; and premature beats capture and propagate far below DI 0.330.
The CV discrepancy survives a fourfold mesh refinement and a fivefold time-step
refinement, so it is not numerical.

The retained table also never sampled the steep region: its lowest ordinate,
2.03 m/s, is higher than the 1.40 m/s measured here at DI90 62 ms, while
conduction fails entirely by 25 ms. Whatever axis those values belong to, they
do not describe restitution of this cable.

Its branch traces and generating script were not retained, so the cause cannot
be established from repository artifacts. The values are replaced rather than
reconciled.

## Reproducing this

The sweep specifications in this directory are runnable through the external
orchestration add-on; each writes a protocol sidecar recording the applied
schedule, and each case emits an event summary, a long-form restitution row and
a per-beat event CSV whether or not it captures. Blocked, censored and
decremental branches are outcome data and are retained in the output.
