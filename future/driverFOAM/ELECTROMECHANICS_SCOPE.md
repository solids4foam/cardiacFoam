# Electromechanics — deliberately out of scope (as of 2026-08-18)

EM work is **not** part of the current driverFOAM effort. This file records every
EM entanglement found while normalizing dicts and wiring overrides, so none of
it is rediscovered as a surprise later.

## 1. `monodomainTotalLagrangianEM` blockMeshDict convention — REVERTED, deferred

The blockMeshDict normalization (explicit `-dict system/blockMeshDict.<dim>`,
`.3D` = 20^3) was applied repo-wide *except* here. Changes were made and then
reverted on 2026-08-18 because this is EM work.

What a future pass would need to redo:
- `Allrun:69` — `runApplication blockMesh` -> `runApplication blockMesh -dict system/blockMeshDict.3D`
- `system/blockMeshDict.3D` — currently 80x80x80, but the case actually builds
  40x40x40 via the plain unsuffixed `system/blockMeshDict`. The `.3D` file is an
  orphan resolution nothing reads.
- `system/blockMeshDict` (plain, 40^3) — removal candidate once the above lands.

The driver already selects explicitly: `manufactured_monodomain_total_lagrangian_em.py:300`
uses `f"system/blockMeshDict.{dimensions_list[-1]}"` with `DIMENSIONS=["3D"]`, so
the driver builds 80^3 while the tutorial's own `Allrun` builds 40^3. **These
disagree today.** No regression test exists for this case, so nothing certified
is contradicted — which is also why nothing catches it.

## 2. `electroMechanicalNiedererEtAl2011` — expected regression skip

`tutorials/Alltest-regression` `isExpectedSkip()` skips this case in
`lightweight` build mode (exit 77). It only runs under `with-solids4foam`.
Consequence: every lightweight run leaves EM untested, including CI.

## 3. Active-tension artifact prediction — ONE real defect (singleCell)

`plugins/cardiacfoam/artifacts_predictor.py:198` `_predict_active_tension` is
called unconditionally at line 257, outside the `_SOLVER_HANDLERS` dispatch, and
takes only `case_root` -- it never sees the solver.

**The real defect.** `electrophysiologyProtocols/singleCell` declares
`activeTensionModel LandNiederer` but is a pure ODE run that writes **no time
directories at all**. Active tension goes to
`postProcessing/<tissue>_<protocol>_Ta.txt` (e.g.
`TWorld_endocardialCells_S1_1000_Ta.txt`). The predictor emits a required
`{time}/AV_Ta` (`time_indexed=True`, `optional=False`), so the solver succeeds
and the driver then rejects the run with `missing_artifacts`. This blocks
singleCell under `foamctl run --strict` today. Note this is a *non-EM* case, so
despite living in this document it is not gated on EM work resuming.

**Fix shape (agreed 2026-08-18).** One branch, not full solver dispatch: if the
myocardium solver is `singleCellSolver`, emit the `postProcessing/*_Ta.txt`
artifact as non-time-indexed, using an explicit filename-pattern helper --
the solver hardcodes those names today, so keep the regex in one named function
rather than deriving it. Whether that stays hardcoded is a future discussion,
to be had together with the `export ()` correction (`f5f935c2`), since both are
the same question: how much of the artifact contract the driver *derives*
versus *mirrors*.

**NOT a defect — corrected 2026-08-18.** An earlier draft of this file claimed
`produced_by="sequentialElectroMechanical"` was wrong for a non-EM tissue run
with active tension. That was incorrect. `activeTensionModel::New` has exactly
two construction sites in the whole tree:

    src/electroMechanicalModels/sequentialElectroMechanical/sequentialElectroMechanical.C:71
    src/electroModels/myocardiumModels/singleCellSolver/singleCellSolver.C:189-191

`monodomainSolver` and `bidomainSolver` never build one. A 3D tissue run
therefore *cannot* have active tension without going through
`sequentialElectroMechanical` (the solid part), so that `produced_by` value is
correct wherever a tissue case declares active tension. Only `singleCellSolver`
constructs one standalone, and it guards on
`if (electroProperties().found("activeTensionModel"))`.

**Still open, unverified.** `electroMechanicalNiedererEtAl2011` keeps its
electroProperties at `constant/electro/`, not `constant/electroProperties`, so
`_predict_active_tension`'s `properties.exists()` guard returns `()` early and
predicts **no** active-tension artifact at all -- in the one case where `AV_Ta`
genuinely belongs. Not re-checked since the corrections above; verify before
acting.

## 4. `AV_Ta` is an ALGEBRAIC name, not a field name

`AV_Ta` is an index into the model's ALGEBRAIC array
(`src/activeTensionModels/LandNiederer/LandNiederer_2017Names.H:67`,
`src/ionicModels/TWorld/TWorld_2024Names.H:388`). How it materializes depends
entirely on solver family: a tissue solver turns exported algebraics into
volScalarFields in time directories; `singleCellSolver` writes a `.txt` series.
The predictor assumes the tissue shape unconditionally.

## 5. solids4foam case-introspection fields — long-term deferral

`plugins/cardiacfoam/case_introspection.py` deliberately omits `D`, `DD`,
`sigmaHyd` from `_SOLID_SOLVER_FIELDS`, documented in-code with a `NOT YET ADDED`
comment and verified field-name evidence. Simao, 2026-08-17: *"that solids4foam
adding will not be done soon."* Do not surface as near-term work.

## 6. EM verification model is compiled and live

`src/verificationModels/electromechanicsVerification/` is in the build. The
`manufacturedElectromechanicsVerifier` emits norms like every other verifier, so
EM is covered by the same `postProcessing/*.dat` mechanism — no special handling
needed when EM work resumes.
