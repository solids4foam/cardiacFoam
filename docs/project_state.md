# cardiacFoam Project State

Last reviewed: 2026-06-14

This note captures the current architecture assumptions that matter for
development. It is a short repository-level state snapshot, not a full design
document. If it disagrees with code, the code wins.

## What changed relative to older cardiacFoam layouts

- The top-level electrophysiology selector is now `myocardiumSolver` in
  `constant/electroProperties`.
- `singleCellSolver` remains a direct `electroModel` subclass.
- `monodomainSolver`, `bidomainSolver`, and `eikonalSolver` are routed through
  `electrophysiologyModel`, which registers under those names with
  `addNamedToRunTimeSelectionTable(...)`.
- Spatial PDE kernels live in the separate `myocardiumSolver` runtime table.
- Older flat names such as `monoDomainElectro`, `eikonalDiffusionElectro`, and
  `singleCellElectro` are obsolete in this tree.

## Current electro stack

- `physicsModel` is selected from `constant/physicsProperties`.
- Electrophysiology runs enter `electroModel::New(...)`.
- The multi-domain orchestration path is owned by `electrophysiologyModel`.
- Domain assembly is centralized in
  `src/electroModels/core/system/electrophysicsSystemBuilder.C`.
- The builder configures myocardium, conduction, ECG, and optional
  extracellular / bath domains.

## Current ionic model stack

Scalar CPU models compiled in `libionicModels`:

- `AlievPanfilov`
- `BuenoOrovio`
- `Courtemanche`
- `Fabbri`
- `Gaur`
- `Grandi`
- `PerisYague`
- `Stewart`
- `TNNP`
- `ToRORd_dynCl`
- `Trovato`
- `TWorld`

Batched models:

- `AlievPanfilovBatched`
- `BuenoOrovioBatched`
- `CourtemancheBatched`
- `FabbriBatched`
- `GaurBatched`
- `GrandiBatched`
- `PerisYagueBatched`
- `StewartBatched`
- `TNNPBatched`
- `ToRORd_dynClBatched`
- `TrovatoBatched`
- `TWorldBatched`

Manufactured / verification models:

- `monodomainFDAManufactured`
- `bidomainFDAManufactured`
- `bathBidomainFDAManufactured`

Notes:

- Optional CUDA sources are conditionally included through
  `src/ionicModels/Make/files-gpu`.
- Older model summaries that mention `ORd` or `tmanufacturedFDA` are stale for
  this repository state.

## Current build reality

- `Allwmake` resolves `solids4foam` through `etc/resolveSolids4Foam.sh`.
- If a built `solids4foam` is unavailable, the build falls back to the bundled
  lightweight `modules/physicsModel` path.
- `FORCE_LIGHTWEIGHT_PHYSICSMODEL=1` forces EP-only mode.
- `src/Allwmake` currently builds:
  - `genericWriter`
  - `ionicModels`
  - `activeTensionModels`
  - `electroModels`
  - `verificationModels`
  - `electroMechanicalModels` only outside lightweight mode

## Files to read before architectural changes

- `README.md`
- `src/electroModels/README.md`
- `src/ionicModels/README.md`
- `src/electroModels/core/electroModel.C`
- `src/electroModels/core/electrophysiologyModel/electrophysiologyModel.C`
- `src/electroModels/electroDomains/myocardiumDomain/myocardiumSolver.H`
- `src/ionicModels/Make/files`
- `src/Allwmake`
- `etc/resolveSolids4Foam.sh`
