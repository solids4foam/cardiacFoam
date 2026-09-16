# Electrophysiology protocols

Standard electrophysiology protocols and sanity checks: the things a cardiac
EP practitioner would run before committing to a large 3D or coupled run —
restitution, conduction velocity, transmural heterogeneity, and re-entry
behaviour — on small, fast, meshless-or-1D/2D cases. Contrast
`manufacturedSolutions/` (verifies the numerics against a known solution) and
`idealizedHeart/` (a full anatomical case using these same protocols).

```text
electrophysiologyProtocols/
├── singleCell/                 single-point AP runs and ionic-model sweeps
├── ionicHeterogeneity/         transmural heterogeneity probe (meshless)
├── restitutionCurves_s1s2Protocol/  S1-S2 APD restitution curves
├── cableProtocol/               1D conduction-velocity calibration
│   ├── monodomain1DCableCV/
│   └── eikonal1DCableCV/
├── rotorInstability/            sustained re-entry / spiral-wave check
└── purkinjeRestitution2D/       multi-beat Purkinje network coupled to a 2D slab
```

## `singleCell/`

Single integration-point electrophysiology: runs an ionic model with no
spatial PDE. The baseline sanity check for any ionic-model change, and the
starting point for the heterogeneity and restitution protocols below.

## `ionicHeterogeneity/`

Runs the meshless `ionicHeterogeneityProbe` utility to check that the
`ionicHeterogeneity` transmural blend itself (endo/M-cell/epi bands, with
smoothed transitions) is stable — the action potential varies smoothly
across the transmural coordinate, with no shape breakdown at a band
boundary — before that same `ionicHeterogeneity` configuration is trusted
in a spatial myocardium case. It checks the blending mechanism itself; for a
specific ionic model's physiology or a disease state, see "Pre-checking
tissue types and pathologies" below.

## `restitutionCurves_s1s2Protocol/`

Single-cell S1-S2 pacing sweeps that generate APD/conduction restitution
curves — the standard check for alternans risk and rate-dependent behaviour
before a tissue-scale run.

## `cableProtocol/`

Minimal 1D tissue-scale conduction-velocity calibration, once per myocardium
solver (`monodomain1DCableCV`, `eikonal1DCableCV`), used to tune conductivity
against a target CV before that solver's conductivity is trusted in a larger
case.

## `rotorInstability/`

A 2D slab workflow that initiates and sustains re-entrant spiral waves
(rotors) with multiple stimuli — checks that the solver and numerical
scheme can sustain self-perpetuating re-entry without spuriously damping or
blowing up.

## `purkinjeRestitution2D/`

A small Purkinje tree coupled to a 2D monodomain slab, paced over several
beats. Checks the network's escape rhythm, capture and block at short
coupling intervals, retrograde activation from the tissue, and the
reaction-diffusion network coupling in parallel.

## Pre-checking tissue types and pathologies

`singleCell` and `ionicHeterogeneity` share a role: confirm the typical
action-potential shape/duration for a tissue type or an ionic
parametrisation (e.g. an `ionicConstantOverrides` block, as used by
`../idealizedHeart/pathos/ionicPathology`) at single-cell or meshless scale
before a tissue-scale run. `singleCell` supports comparisons of this kind
as a study under `setup/studies/` — descriptive figures and a metrics
table, not a regression gate (see its README's `tworldVsGaur` study).

## Regression

`singleCell`, `ionicHeterogeneity`, `rotorInstability` and
`purkinjeRestitution2D` are wired into
`tutorials/Alltest-regression`; see `../README.md` for the full canonical
table and coverage status of the rest.
