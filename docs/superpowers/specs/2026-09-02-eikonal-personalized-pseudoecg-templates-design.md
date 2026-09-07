# Personalized pseudoECG templates for eikonalECG

## Problem

`eikonalECG` (`src/electroModels/ecgModels/eikonalECG/`) reconstructs a
pseudo-ECG from an eikonal activation-time field by the surrogate
`Vm(x,t) = U(t - psi(x))`, where `U` is a Vm(t) action-potential template.
Today `U` comes from three **hardcoded, compile-time** arrays
(`tissueTemplates.H`) — fixed TNNP endo/mid/epi curves generated once from a
reference `singleCell` run and baked into the binary. `calculateTransmuralWeights`
also only accepts `ionicHeterogeneity.mode == transmuralBands` and fatals on
anything else.

Full monodomain, meanwhile, already supports per-case ionic personalization:
`ionicHeterogeneityOrchestrator` resolves a full per-cell ionic-constants
field by composing `transmuralBands` / `namedRegions` / `cellZoneRegions`
(the base mode) with an optional `apexBaseBands` continuous overlay (e.g. a
longitudinal-distance-driven exponential scale on named constants). The
eikonal ECG path ignores all of this — a case with drug-scaled GKr, or a
patient-specific transmural gradient, gets the same generic ECG shape as any
other case.

**Goal:** let `eikonalECG` build its Vm(t) templates from the *same* ionic
parameterization the case's `ionicHeterogeneity` block would give monodomain,
by running cheap single-cell simulations instead of the full 3D
reaction-diffusion solve.

## Non-goals

- No change to the eikonal activation-time solve itself, the lead-field
  construction, or the manufactured-solution verification path.
- No attempt to give every mesh cell a literally unique template — heterogeneity
  is discretized onto a bounded anchor set (see below), not resolved
  per-cell.
- No new driverFOAM orchestration stage. This is a single-binary, in-process
  pipeline; driverFOAM launches the case exactly as it does today.

## Architecture

Trigger: a new `personalizedTemplates{}` block nested inside the `eikonalECG`
config. Its **absence is the fallback** — existing cases keep using the
compiled generic templates and the current `transmuralBands`-only weighting,
unchanged. Its presence runs a three-step in-process pipeline before the
existing eikonal solve/ECG loop:

1. **Anchor generation** ("run the 1D system"). A component factored out of
   `applications/utilities/ionicHeterogeneityProbe/` (so the standalone CLI
   probe and this in-process path share one implementation) builds an anchor
   set from the case's own `ionicHeterogeneity` block:
   - `transmuralBands`: K evenly-spaced anchors along `t in [0,1]`.
   - `namedRegions`: one anchor per named region, at its reference field
     value.
   - `cellZoneRegions`: one anchor per named cell zone (discrete, no field
     value) — new support; the probe currently refuses this mode.
   - If `apexBaseBands` is also configured (composed on top, exactly as
     `myocardiumDomainInterface.C` already does for monodomain), each base
     anchor expands into M sub-anchors along the apex-base field, giving a
     K x M grid.

   Each anchor's ionic constants are resolved via the same
   `ionicHeterogeneityOrchestrator` code monodomain uses, then integrated to
   steady state with the batched ionic-model machinery
   (`ionicModel::New(dict, nAnchors, dt, true)` — the same "cells" abstraction
   a real tissue run uses, just sized to ~K or K*M instead of millions of
   mesh cells; runs on GPU automatically if the case's ionic model has a
   batched/GPU variant, e.g. `BuenoOrovioBatched`).

2. **Populate the cells.** For every real mesh cell, compute blend weights
   against the anchor grid: the base mode's own weighting (transmural
   smoothstep / namedRegions blend / cellZone hard-assignment) combined with
   linear interpolation across the apex-base axis when present. This
   generalizes today's fixed `wEndo_/wMid_/wEpi_` triple to a `K`- or
   `K*M`-length weight vector per cell.

3. **Proceed to the existing 3D eikonal solve + ECG reconstruction**,
   structurally unchanged — `reconstructGradVm` still evaluates
   `Vm(x,t) = U(t-psi(x))` via the chain rule, just blending across however
   many anchors are active instead of a hardcoded 3.

## Dict schema

```c++
eikonalECG
{
    ecgSolver eikonalECG;
    sampling { start 0; end 0.5; deltaT 0.001; }
    electrodePositions { ... }

    personalizedTemplates
    {
        nAnchors          5;     // anchors along the base mode's field
        nApexBaseAnchors  5;     // anchors along apexBaseBands field, if active
        nBeats            10;    // pace to steady state before capture
        duration          1.0;   // s, capture window
        dt                0.0001;// s, ODE integration step for anchor runs
        singleCellStimulus { ... }  // optional; else reuse case's own protocol
    }
}
```

The case's existing `ionicHeterogeneity` block (mode, `ionicConstantOverrides`,
optional `apexBaseBands`) is read as-is; nothing new to configure there.

## Runtime storage

`tissueTemplates.H`'s compiled static arrays are replaced, on this path, by
in-memory `List<scalarField>` (times/values per anchor) populated once at
`solve()` startup, before the sampling loop. `evaluateTemplateDerivative`
(the existing interpolation helper) is reused unchanged against the dynamic
anchors. The compiled header and its 3-array path remain as the
no-`personalizedTemplates` fallback.

## Validation

Reuse `ionicHeterogeneityProbe`'s existing sanity checks (adjacent-anchor
waveform-RMS smoothness, APD monotonic ordering, APD envelope bounding) as a
gate on generated anchors before the 3D solve starts — `FatalError` on an
unphysiological anchor set rather than silently producing a bad ECG. Cap
total anchor count (K, or K*M when apex-base is stacked) with a sane default
and a hard ceiling that fatals with guidance to lower resolution, since each
anchor is a real (cheap but non-free) paced single-cell run.

## Testing strategy

- **Regression parity:** with `personalizedTemplates` absent, all existing
  `eikonalECG` behavior (compiled templates, `transmuralBands`-only,
  manufactured-solution verifier) must be bit-identical to today.
- **K=3 equivalence:** the generalized weight-blending function, run with
  `transmuralBands` and `nAnchors=3` at today's default interface positions,
  should reproduce the current fixed endo/mid/epi blend to numerical
  precision — a direct regression check that generalization didn't change
  the K=3 case's math.
- **Shared-logic consistency:** the in-process anchor generator and the
  standalone `ionicHeterogeneityProbe` CLI, given equivalent dicts, must
  produce matching Vm traces (they share one implementation, so this is
  mostly a wiring check).
- **Physiological sanity:** a tutorial case with `personalizedTemplates` and
  a non-default `ionicConstantOverrides` (e.g. GKr scaled to prolong APD)
  should show the eikonalECG output shift accordingly (later T-wave timing)
  relative to the generic-template baseline — the same kind of steady-state
  cross-check used previously to verify the hardcoded TWorld templates
  against `singleCell` output.
