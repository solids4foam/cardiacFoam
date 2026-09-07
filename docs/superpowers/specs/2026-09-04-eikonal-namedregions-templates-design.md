# N-region (namedRegions) personalized templates for eikonalECG

## Problem

`eikonalECG`'s `personalizedTemplates` feature (landed 2026-09-02) generates
its three Vm(t) action-potential templates by pacing single cells, but it
only understands `ionicHeterogeneity.mode transmuralBands`: exactly three
fixed templates (`TemplateTriplet{endo, mid, epi}`), blended per cell via
`wEndo_/wMid_/wEpi_` computed by directly calling
`ionicHeterogeneity::transmuralBandWeights()`. Any other `mode` fatals at
first `solve()`.

`ionicHeterogeneity` itself already supports a strictly more general mode,
`namedRegions`: an arbitrary number of named regions tiling `[0,1]` on any
named scalar field (`field <word>;`, default `"t"` — the field name is
already configuration, not hardcoded), each resolving to a baseline tissue
type (`epicardialCells`/`mCells`/`endocardialCells`/`myocyte`) plus optional
`ionicConstantOverrides.<name>`, blended via `namedRegionWeightsAt()` using
the identical smoothstep transition math as `transmuralBandWeights()`,
generalized from 2 interfaces to N-1. Monodomain already uses this for
arbitrary tissue partitions; `eikonalECG` cannot.

**Goal:** let `eikonalECG`'s `personalizedTemplates` consume
`mode namedRegions` — N templates on any named coordinate, not just the
fixed endo/mid/epi triple on a hardcoded transmural axis — using the same
per-region baseline/override resolution monodomain already relies on.

## Non-goals

- No change to `mode transmuralBands`'s *behavior* — its output must remain
  numerically identical to today, verified by regression tests (see Testing
  strategy). Its underlying storage and blend loop are not kept as a
  separate, untouched code path: they are unified with `namedRegions` onto
  one generalized representation (see Architecture), with `transmuralBands`
  becoming the N=3 case of it.
- No support for `mode apexBaseBands`. That mode continuously rescales named
  ionic variables by a smooth function of position; it does not select among
  a small set of pre-characterized cell types. Approximating it with a
  template grid would require sampling the rescaling axis at multiple points
  and characterizing the resulting interpolation error — a separate, harder
  problem with no existing precedent in this codebase, explicitly deferred.
- No support for `mode cellZoneRegions`. Structurally close to
  `namedRegions` (same baseline/override resolution) but selection is by
  cell-zone membership, not a scalar field — no blending, no coordinate.
  Deferred as its own near-trivial follow-up once this lands.
- No reconciliation with the original 2026-09-02 design's broader envisioned
  architecture (an anchor-generation component factored out of
  `applications/utilities/ionicHeterogeneityProbe/`, batched/GPU ionic-model
  integration for anchor pacing). Phase 1 built a simpler, standalone
  single-cell pacer (`eikonalTemplateGenerator`) instead of that shared
  component; this phase deliberately continues extending that actual
  implementation rather than reconciling with the earlier, unbuilt vision.
  That divergence is accepted, not resolved, here.
- No attempt to give every mesh cell a literally unique template —
  heterogeneity is still discretized onto the region set, now sized N
  (from the case's `regions{}` dict) instead of fixed at 3.

## Architecture

Trigger: no new configuration surface. A case using `personalizedTemplates`
behaves identically to today if its `ionicHeterogeneity.mode` is
`transmuralBands`; it gains N-region support if `mode` is `namedRegions`.
Both modes now flow through one generalized representation instead of the
old mode being a separate, fixed-3 special case.

**1. One generalized template/weight representation.** Replace the fixed
`TemplateTriplet{endo, mid, epi}` and the 3 fixed `wEndo_/wMid_/wEpi_`
`scalarField` members with `List<DynamicTemplate>` and `List<scalarField>`
respectively, both index-aligned to an ordered list of region names.
`transmuralBands` and `namedRegions` each populate this same pair of lists
via their own, unchanged, mode-specific parsing/weight functions (below) —
only the storage and the blend loop (`reconstructGradVm`) are shared.

**2. Template generation (`eikonalTemplateGenerator`).** The existing
generator does *not* resolve per-tissue-type ionic constants itself. It
builds one `ionicModel` sized to the anchor count (`ionicModel::New(modelDict, nPoints, dt*1000.0, true)`)
and calls `model.configureIonicHeterogeneity(tPoints, heterogeneityDict)` —
the same virtual function the real monodomain solve calls per-mesh-cell
(`configuredIonicModel.H:151-161`), which forwards to
`ionicHeterogeneityOrchestrator::configureTransmuralBandHeterogeneity()`.
That orchestrator function already dispatches on `heterogeneityDict.mode`
(`ionicHeterogeneityOrchestrator.C:157-178`): `namedRegions` routes to
`configureNamedRegionHeterogeneity()` unchanged, today, with zero code
changes required. So `namedRegions` support for constant resolution already
exists — nothing to build there.

What the generator *does* still need, per mode, is which representative
point each anchor is paced at:
- `transmuralBands` (unchanged): `tPoints = [0.0, tMid, 1.0]`, `tMid`
  computed from `endoMInterface`/`mEpiInterface`/`transitionWidth` exactly
  as today.
- `namedRegions` (new): call `ionicHeterogeneity::parseNamedFieldRegions()`
  on `heterogeneityDict.subDict("regions")` to get the sorted `List<NamedFieldRegion>`
  (already validated for tiling/baseline names/minimum count), and use each
  region's midpoint, `0.5*(rangeMin+rangeMax)`, as its `tPoints` entry — the
  direct generalization of how `tMid` is already the mid-band's own
  midpoint-of-effective-range today.

Both branches then set `nPoints` to the anchor count (3, always, for
`transmuralBands`; N for `namedRegions`), construct the model with that many
integration points, and pace/capture one trace per point with the case's
`singleCellStimulus`, using the existing pacing/capture/validation logic
unchanged — generalized from the current hardcoded `traces[3]` array to a
`List<DynamicTemplate>` sized `nPoints`.

**3. Per-cell weights (`eikonalECG::calculateTransmuralWeights`),
unified via one shared region list.** Read the field named by
`ionicHeterogeneity.field` (default `"t"`) exactly as today. Both modes
now produce a `List<NamedFieldRegion>` and call `namedRegionWeightsAt()`
against it: `namedRegions` calls the existing `parseNamedFieldRegions()`
directly; `transmuralBands` builds its region list via a small new shared
helper, `ionicHeterogeneity::synthesizeTransmuralBandRegions(endoMInterface, mEpiInterface)`,
extracted from the inline 3-region construction
`ionicHeterogeneityOrchestrator::configureTransmuralBandHeterogeneity()`
already does (`ionicHeterogeneityOrchestrator.C:212-224`) — a pure
refactor of existing code into a named, reusable function, called from both
that orchestrator (no behavior change) and `eikonalECG` (new use). This
mirrors the orchestrator's own stated intent ("keeps exactly one code path
for region-based heterogeneity") one level up, into `eikonalECG`.

Both modes then write into the same generalized `List<scalarField>`,
index-aligned to the mode's region list. `namedRegionWeightsAt()` returns a
sparse `List<NamedRegionWeight>` (1 or 2 entries, keyed by region *name*)
per cell; build a `name -> index` lookup once over the small, fixed region
list at setup time, and scatter each cell's sparse result into the aligned
weight fields (inactive regions default to 0, matching `scalarField`'s
zero-initialization).

**4. Blend loop (`eikonalECG::reconstructGradVm`), unified.** Replace the
fixed three-term sum with a loop:
`forAll(regions, i) rawDUds += weights[i][cellI] * evaluateTemplateDerivative(templates[i], localTime);`
For `transmuralBands` this evaluates to exactly the same three terms in the
same order as today, verified by the N=3 equivalence test (see Testing
strategy) rather than by leaving the code path untouched.

**5. Mode gating.** The existing first-`solve()` fatal ("mode other than
transmuralBands is rejected") widens to accept `namedRegions` in addition;
every other mode (`apexBaseBands`, `cellZoneRegions`) keeps fataling with
the existing diagnostic.

## Dict schema

No new keys anywhere. A case enables this purely by setting
`ionicHeterogeneity.mode namedRegions;` with its (already-existing-schema)
`regions{ <name> { range (min max); baseline <word>; } ... }` sub-dict, and
`eikonalECG.personalizedTemplates{...}` (unchanged schema from phase 1). The
generalization is entirely on the consuming side inside `eikonalECG`.

## Validation

- Reuse `parseNamedFieldRegions()`'s existing fatals unchanged (tiling,
  reserved name `global`, valid baseline names, minimum two regions) —
  `eikonalECG` calls the identical, already-validated function monodomain
  uses, so no new parsing-time checks are needed.
- Generalize phase 1's existing "generated trace must be finite and
  non-degenerate" per-template check from a fixed 3-iteration loop to a loop
  over N — same check, same fail-loud diagnostic style, no new validation
  category introduced.
- Phase 1's construction-time fatals (missing `ionicModelConfig`,
  `nBeats < 1`, non-positive `duration`/`dt`, etc.) are mode-independent and
  unchanged.

## Testing strategy

- **Regression parity (now a verified property, not a structural
  guarantee):** since storage and the blend loop are unified across both
  modes, existing `transmuralBands` fixtures (personalized or not) must be
  re-run and produce output numerically identical to pre-change baselines,
  captured before this work starts. This is the primary safety net given
  the shared code is no longer literally untouched.
- **N=3 equivalence (real correctness check, not just a sanity check):** a
  `namedRegions` config with 3 regions whose ranges match a `transmuralBands`
  config's `endoMInterface`/`mEpiInterface` (e.g. ranges `[0,0.3]`,
  `[0.3,0.7]`, `[0.7,1]` against `endoMInterface 0.3; mEpiInterface 0.7;`),
  same `transitionWidth`/`smoothing`/`transitionMode`, must produce
  bit-identical per-cell weights to `transmuralBandWeights()` — both
  functions apply the same smoothstep formula with the same boundary
  convention (`ionicHeterogeneity.C:417-419`).
- **New coverage:** 2-region and 5-region `namedRegions` cases verifying
  per-cell weights sum to 1 everywhere, each generated template's resolved
  constants match directly calling `constantsForRegion()`/
  `initialStatesForRegion()` for that region's name/baseline, and the
  resulting pseudo-ECG shifts sensibly under a region-specific override
  (e.g. a region with a GKr-prolonging override shows a later T-wave) — the
  same style of physiological cross-check phase 1 already used for its
  endo/mid/epi templates.
- **driverFOAM regression:** extend the existing
  `test_eikonal_ecg_personalized_templates_regression.py`-style coverage
  with a `namedRegions` tutorial fixture, matching how phase 1 added its own.

## Completion criteria

- `transmuralBands` (personalized or not) produces output identical to
  today, unchanged.
- `namedRegions` works end-to-end for `personalizedTemplates`: N-region
  cases pace correctly, per-cell weights sum to 1, the blend loop applies
  correctly, and the resulting pseudo-ECG is sane.
- The N=3 `namedRegions`/`transmuralBands` equivalence test passes to
  numerical precision.
- Every generated per-region template is validated before use;
  misconfiguration fails loud with a useful diagnostic, consistent with
  phase 1.
- `README.md` (`src/electroModels/ecgModels/eikonalECG/`) is updated to
  document `namedRegions` support, explicitly listing `apexBaseBands` and
  `cellZoneRegions` as still unsupported and deferred.
