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

## 3. Active-tension artifact prediction — two opposite defects

`plugins/cardiacfoam/artifacts_predictor.py:198` `_predict_active_tension` is
called unconditionally at line 257, outside the `_SOLVER_HANDLERS` dispatch:

- **Over-prediction:** any case declaring `activeTensionModel` gets a required
  `{time}/AV_Ta` artifact (`time_indexed=True`, not optional, `produced_by`
  hardcoded to `sequentialElectroMechanical`). `electrophysiologyProtocols/singleCell`
  declares `activeTensionModel LandNiederer` but is a pure ODE run that writes
  **no time directories at all** — active tension goes to
  `postProcessing/<tissue>_<protocol>_Ta.txt`. The solver succeeds; the driver
  then rejects the run with `missing_artifacts`. This blocks singleCell from
  running under `foamctl run --strict`.
- **Under-prediction:** `electroMechanicalNiedererEtAl2011` keeps its
  electroProperties at `constant/electro/`, not `constant/electroProperties`, so
  `_predict_active_tension` finds no file, returns `()`, and predicts **no**
  active-tension artifact at all — in the one case where `AV_Ta` genuinely belongs.

Both stem from the same function ignoring solver family. Neither is fixed.

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
