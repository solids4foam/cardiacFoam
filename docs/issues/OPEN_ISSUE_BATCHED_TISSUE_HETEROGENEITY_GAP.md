# Open issue: 8 batched ionic models cannot do per-tissue heterogeneity

**Status:** open, deferred. Diagnosed 2026-08-06, no code changed.

## Summary

Eight of the twelve `*compactBatched` ionic models are missing
`constantsForTissue()` and `initialStatesForTissue()`. As a result they cannot
re-derive constants for an arbitrary tissue flag, which is what per-zone /
named-region heterogeneity requires when different cell zones select different
tissues within one run.

Uniform single-tissue cases are **unaffected** — those work correctly on all
twelve.

## The split

| | tissues advertised | `constantsForTissue` | `initConsts` uses `tissueFlag` |
|---|---|---|---|
| AlievPanfilov, Courtemanche, Fabbri, Gaur, Grandi, PerisYague, Stewart, Trovato | `myocyte` only | **absent** | **no (0 uses)** |
| BuenoOrovio, TNNP, TWorld, ToRORd_dynCl | epi, M, endo | present | yes (21, 4, 103, 58 uses) |

The correlation is exact. Whoever did the batched ports appears to have used
"does the codegen branch on `tissueFlag`?" as the criterion for whether the
tissue layer was worth porting. That is the wrong test, because tissue is also
meaningful through the **override/scaling layer** even when the codegen ignores
it (see below).

## Why tissue still matters when `initConsts` ignores it

`ionicModelIO::applyConstantOverrides(..., tissueFlag)` maps the flag to a
dictionary **scope name**:

```cpp
word tissueOverrideScopeName(const label tissueFlag)
{
    return tissueFlag == 1 ? "epicardialCells"
         : tissueFlag == 2 ? "mCells"
         : tissueFlag == 3 ? "endocardialCells"
         : tissueFlag == 4 ? "myocyte" : word();
}
```

and applies that scope's entries, which support `scale` and `set`. So a user
differentiates tissues from the case dictionary:

```c++
ionicConstantOverrides
{
    global           { scale { ... } }
    epicardialCells  { scale { ... } }
    endocardialCells { scale { ... } }
}
```

This is real functionality independent of whether the model equations carry
built-in transmural constants.

## What is NOT broken (two corrections made while diagnosing)

Both of these were asserted during investigation and then disproved. They are
recorded so the same wrong turns are not repeated.

1. **"The 8 never call `applyConstantOverrides` at all."** False. All twelve
   call `applyIonicConstantOverrides()` at construction, which routes to
   `ionicModelIO::applyConstantOverrides(..., tissue_)`. The original grep
   looked for the symbol in each model's own `.C` and missed the inherited
   base-class call. Both `global` and tissue-scoped scaling work on all twelve
   for the model-wide tissue.

2. **"The batched restriction to `myocyte` is the honest one, and the scalar
   list is cosmetic."** False, and backwards. It was reached by inspecting
   `initConsts` alone and stopping before the override layer. The tissue labels
   are meaningful as override scopes, so restricting the list removes reachable
   functionality.

The genuinely missing piece is narrower than either claim: the *heterogeneity*
entry point, not the *override* entry point.

## Impact

- Uniform single-tissue runs: **fine on all twelve models.**
- Per-zone / named-region heterogeneity assigning different tissues to
  different zones: **not possible on the 8**, because the tissue cannot be
  selected and the constants cannot be re-derived per flag.
- Scalar vs batched divergence: the scalar twin of each of the 8 accepts
  `{epicardialCells, mCells, endocardialCells, myocyte}`, so the same case
  dictionary is not portable between a model and its batched variant.

## Corresponding driverFOAM state

The catalog on branch `named-region-ionic-heterogeneity` (commit `693df051`,
"distinguish native vs override-only tissue labels in the ionic catalog")
already encodes this distinction, and its split matches the C++ exactly:

```python
native_tissue_labels        # distinct baseline physiology in the equations
approximate_tissue_labels   # accepted only through overrides or scaling
```

with the batched entry inheriting its parent's labels only when
`heterogeneity_capable`, else collapsing to `("myocyte",)`. That flag is
currently a *true statement about incomplete code* rather than about
physiology.

Note that branch is 217 commits behind `no-frontend-minor-errors` and 0 ahead,
so this work is not on the open PR.

## What a fix involves

1. Add `constantsForTissue(tissueFlag)` and `initialStatesForTissue(tissueFlag)`
   to the 8, modelled on `TNNPBatched`. Each needs its `*Names.H` included so
   `<Model>CONSTANTS_NAMES` is in scope — the 8 do not currently include it and
   have no `ioConstantNames()` override, so this is per-file work, not a
   templated edit.
2. Decide what `supportedTissueTypes()` should advertise. **This is a design
   decision, not a mechanical one.** For these 8 the codegen ignores
   `tissueFlag`, so all four labels would return the same baseline and differ
   only through user-supplied overrides. Widening the list advertises four
   labels backed by one physiology; leaving it advertises one label while the
   scalar twin advertises four. Neither is obviously right.
3. Once (1) lands, `heterogeneity_capable` becomes true for all twelve and the
   special-case branch in the driverFOAM catalog becomes dead code.

## Not yet checked

Whether the same 8 also lack `constantsForRegion` (the named-region / cellZone
path). If so it is a third missing piece and should be fixed in the same pass.
