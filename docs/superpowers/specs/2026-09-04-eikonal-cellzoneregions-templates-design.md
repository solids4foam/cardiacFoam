# cellZoneRegions personalized templates for eikonalECG

## Problem

`eikonalECG`'s `personalizedTemplates` now supports `ionicHeterogeneity.mode`
`transmuralBands` (fixed 3 templates) and `namedRegions` (N templates on any
named scalar field, blended smoothly) — see
`docs/superpowers/specs/2026-09-04-eikonal-namedregions-templates-design.md`.
It still fatals on `mode cellZoneRegions`, which the monodomain path already
supports: partitioning tissue by literal mesh `cellZone` membership (e.g.
healthy myocardium / border zone / dense scar in an infarct model), rather
than by a continuous coordinate. A case using `cellZoneRegions` for its real
physics falls back to the generic 3-curve compiled ECG, losing that
region-specific electrophysiology from the pseudo-ECG entirely.

**Goal:** let `eikonalECG`'s `personalizedTemplates` consume
`mode cellZoneRegions` — one template per named cell zone.

## Non-goals

- No change to `transmuralBands` or `namedRegions` — both untouched.
- No blending logic of any kind for `cellZoneRegions`: this mode is a hard,
  crisp per-cell assignment in the monodomain path today
  (`ionicHeterogeneityOrchestrator::configureCellZoneRegionHeterogeneity`
  has no `transitionWidth`/`smoothing`/`transitionMode` parameter at all),
  and stays that way here.
- No new shared reusable "mesh cellZone → per-cell region index" utility
  factored out of `myocardiumDomainInterface.C`'s `readTransmuralDistance()`
  — that logic is duplicated (not extracted), matching how phase 1 and the
  `namedRegions` work already duplicate small pieces of per-mode logic
  between `eikonalECG.C` and the monodomain path rather than reaching into
  `myocardiumDomainInterface.C` (an electroDomains-layer file `eikonalECG`
  does not otherwise depend on).

## Architecture

**1. Template generation (`eikonalTemplateGenerator.C`).** New anchor
function `cellZoneRegionAnchors()`, parallel to `namedRegionAnchors()` but
simpler — no transition-width guard, since there is nothing to blend:

```cpp
void cellZoneRegionAnchors
(
    const dictionary& heterogeneityDict,
    scalarField& tPoints,
    wordList& anchorNames
)
{
    if (!heterogeneityDict.found("regions"))
    {
        FatalErrorInFunction
            << "eikonalTemplateGenerator: ionicHeterogeneity mode "
            << "cellZoneRegions requires a 'regions' sub-dictionary."
            << exit(FatalError);
    }

    const List<ionicHeterogeneity::NamedCellZoneRegion> regions =
        ionicHeterogeneity::parseNamedCellZoneRegions
        (
            heterogeneityDict.subDict("regions")
        );

    const label nRegions = regions.size();
    tPoints.setSize(nRegions);
    anchorNames.setSize(nRegions);

    forAll(regions, i)
    {
        tPoints[i] = scalar(i);
        anchorNames[i] = regions[i].name;
    }
}
```

`generatePersonalizedTemplates()`'s mode dispatch gains an
`else if (mode == "cellZoneRegions")` branch calling this. Passing
`tPoints[i] = scalar(i)` works because
`ionicHeterogeneityOrchestrator::configureCellZoneRegionHeterogeneity`
resolves a cell's region via `round(regionIndices[cellI])` used as a direct
0-based index into its per-region constants array (not an OpenFOAM cellZone
ID needing lookup) — so anchor `i`, given coordinate `scalar(i)`, resolves
to exactly `regions[i]`'s baseline/overrides via the same
`constantsForRegion()`/`initialStatesForRegion()` calls `namedRegions`
already exercises. No new constant-resolution logic.

**2. Per-cell weights (`eikonalECG::calculateTransmuralWeights`).** A new
`cellZoneRegions` branch, alongside the existing `transmuralBands`/
`namedRegions` handling, builds `regionWeights_` directly from real mesh
`cellZones()` (crisp 1.0/0.0, no field read, no `namedRegionWeightsAt()`) —
mirroring `myocardiumDomainInterface.C:57-98`'s inline logic:

```cpp
else // mode == "cellZoneRegions"
{
    const List<ionicHeterogeneity::NamedCellZoneRegion> czRegions =
        ionicHeterogeneity::parseNamedCellZoneRegions
        (
            hetDict.subDict("regions")
        );

    const label nRegions = czRegions.size();
    regionWeights_.setSize(nRegions);
    forAll(regionWeights_, i)
    {
        regionWeights_[i].setSize(mesh.nCells(), 0.0);
    }

    boolList claimed(mesh.nCells(), false);

    forAll(czRegions, regionIndex)
    {
        const word& zoneName = czRegions[regionIndex].cellZone;
        const label zoneId = mesh.cellZones().findZoneID(zoneName);

        if (zoneId < 0)
        {
            FatalErrorInFunction
                << "eikonalECG personalizedTemplates: cellZone '"
                << zoneName << "' (region '" << czRegions[regionIndex].name
                << "') not found in this mesh."
                << exit(FatalError);
        }

        const labelList& zoneCells = mesh.cellZones()[zoneId];

        forAll(zoneCells, i)
        {
            const label cellI = zoneCells[i];

            if (claimed[cellI])
            {
                FatalErrorInFunction
                    << "eikonalECG personalizedTemplates: cell " << cellI
                    << " is claimed by more than one cellZoneRegions entry."
                    << exit(FatalError);
            }

            claimed[cellI] = true;
            regionWeights_[regionIndex][cellI] = 1.0;
        }
    }

    forAll(claimed, cellI)
    {
        if (!claimed[cellI])
        {
            FatalErrorInFunction
                << "eikonalECG personalizedTemplates: cell " << cellI
                << " is not claimed by any cellZoneRegions entry. Every "
                << "mesh cell must belong to exactly one named cell zone."
                << exit(FatalError);
        }
    }
}
```

This branch does not read `ionicHeterogeneity.field` at all (there is
nothing to read a field for) — matching `readTransmuralDistance()`'s own
`cellZoneRegions` special-casing, which likewise skips the named-field
read path entirely for this mode.

**3. Blend loop (`reconstructGradVm`).** Unchanged — it already just sums
`regionWeights_[i][cellI] * evaluateTemplateDerivative(...)` over however
many regions the active mode produced; a crisp 1.0/0.0 weight vector is a
degenerate case of the same sum, not a new code path.

**4. Mode gating.** Both existing fatals (`calculateTransmuralWeights`,
`generatePersonalizedTemplates` member) widen their accepted-mode check
from `{transmuralBands, namedRegions}` to
`{transmuralBands, namedRegions, cellZoneRegions}`, with `cellZoneRegions`
gated the same way `namedRegions` already is — only valid when
`personalizedTemplatesEnabled_`, since the compiled 3-curve fallback still
cannot represent an arbitrary region count.

## Dict schema

No new keys. `ionicHeterogeneity.mode cellZoneRegions;` with its
already-existing-schema `regions{ <name> { cellZone <word>; baseline <word>; } ... }`
sub-dict (`ionicHeterogeneity::parseNamedCellZoneRegions`, unchanged) and
`personalizedTemplates{...}` (unchanged).

## Validation

- Reuse `parseNamedCellZoneRegions()`'s existing fatals unchanged (missing
  `cellZone` key, duplicate cellZone claimed by two regions, reserved name
  `global`, minimum two regions).
- New: fatal if a named `cellZone` doesn't exist in the mesh (construction
  of `eikonalECG`'s per-cell weights, not `parseNamedCellZoneRegions`,
  since zone existence is a mesh property, not a dict property).
- New: fatal if any mesh cell is claimed by more than one region, or by
  none — full, non-overlapping mesh coverage is required, mirroring the
  monodomain path's de facto behavior (an uncovered cell trips
  `configureCellZoneRegionHeterogeneity`'s `round(regionIndices[cellI])`
  out-of-range fatal there; this makes the same requirement an explicit,
  clearly-worded fatal here instead of an indirect one).
- Generalize the existing per-template finite/non-degenerate trace check
  (already looped over N anchors for `namedRegions`) — no change needed,
  it already handles any N.

## Testing strategy

- **Regression parity:** `transmuralBands` and `namedRegions` paths are
  untouched code (new branches only), so both stay exactly as they are.
- **New coverage:** a 2-zone and a 3-zone `cellZoneRegions` case (mirroring
  the 2-region/3-region `namedRegions` coverage already added), verifying:
  each generated template's baseline is resolved correctly (cross-checked
  against the same baseline reached via `namedRegions`/`transmuralBands`
  for the identical baseline name — reusing the cross-check pattern already
  proven this session, where `endocardialCells`/`mCells`/`epicardialCells`/
  `myocyte` baselines produce identical peak Vm regardless of which mode
  reaches them); an unclaimed-cell mesh fatals with the new diagnostic; a
  double-claimed cell fatals; a missing cellZone name fatals.
- **driverFOAM regression:** a `cellZoneRegions` tutorial fixture, built by
  cloning the existing `eikonalECGPersonalizedNamedRegions` tutorial and
  replacing its 3 field-range regions with 3 `cellZone`-based regions
  covering equivalent anatomical volumes — same style of "reproduces the
  same physiology through a different mode" fixture already used twice
  (transmuralBands ↔ namedRegions).

## Completion criteria

- `transmuralBands` and `namedRegions` produce output identical to today,
  unchanged.
- `cellZoneRegions` works end-to-end for `personalizedTemplates`: N-zone
  cases pace correctly, per-cell weights are crisp (exactly one region
  weight is 1.0, the rest 0.0, for every cell), and the resulting pseudo-ECG
  is sane.
- Every generated per-region template is validated before use; mesh-zone
  misconfiguration (missing zone, double-claimed cell, unclaimed cell)
  fails loud with a useful diagnostic, consistent with `transmuralBands`
  and `namedRegions`.
- `README.md` (`src/electroModels/ecgModels/eikonalECG/`) is updated to
  document `cellZoneRegions` support, removing it from the "not supported"
  list (leaving only `apexBaseBands` there).
