# Duplication and Architecture Audit

Role: Agent 6, Duplication and Architecture Hunter  
Date: 2026-07-11  
Mode: independent read-only discovery (only this audit note was created)

## Scope and standard

I inspected project-owned C++, Python, shell/build code, and the component architecture documents, excluding `modules/solids4foam` as external/submodule content and excluding generated equation bodies as consolidation targets. I only propose consolidation where copied correctness-sensitive policy has demonstrably drifted. Repetition that is model-specific, generated, or merely a few stable lines is not a finding.

## Findings and consolidation proposals

### AD-1 — Active electro-model dictionary selection has multiple incompatible sources of truth

**Severity:** S1 — High  
**Confidence:** high

**Current duplicate locations**

- Canonical runtime selection: `src/electroModels/core/electroModel.C:125-180`, especially `:147-152`, requires `myocardiumSolver` and uses `<myocardiumSolver>Coeffs` through construction at `:109`.
- Builder-side re-derivation: `src/electroModels/core/system/electrophysicsSystemBuilder.C:105-117` reads `myocardiumSolver`, but can infer the type from a dictionary name ending in `Coeffs`.
- Current utility fallback: `applications/utilities/recomputePseudoECG/recomputePseudoECG.C:135-156` reads `myocardiumSolver`, defaults it to `monodomainSolver`, and falls back to the flat root dictionary when the coefficient dictionary is missing.
- Mixed legacy/current utility: `applications/utilities/ionicHeterogeneityProbe/ionicHeterogeneityProbe.C:46-80` first accepts obsolete `electroModel`, then current `myocardiumSolver`, then a flat dictionary.
- Legacy-only utility: `applications/utilities/listCellModelsVariables/listCellModelsVariables.C:44-75` recognizes only `electroModel`; `:106-114` also reports only that legacy selector.

**Concrete evidence of drift**

The repository contract says `myocardiumSolver` is the single canonical entry point, and production construction requires it (`electroModel.C:147-152`). In contrast, `listCellModelsVariables` ignores that key. For a normal current dictionary such as `myocardiumSolver monodomainSolver;` with `monodomainSolverCoeffs`, it returns the entire `electroProperties` root at line 75; its caller then expects `ionicModel` in that returned dictionary (`listCellModelsVariables.C:116-117`). The normal nested contract therefore fails before it can instantiate/report the ionic model. The other utility copies disagree over whether a missing selector defaults, whether legacy `electroModel` is supported, and whether missing coefficients are fatal.

**Meaningful differences and intentionality**

- Reading from a preselected `*Coeffs` dictionary in the builder is intentional internal context.
- Some flat-dictionary tolerance may be intentional for utilities or old cases, but there is no shared compatibility policy and no deprecation boundary.
- The legacy-only behavior in `listCellModelsVariables` is not intentional under the documented current architecture; it is stale copy-and-paste parsing.
- `recomputePseudoECG` defaulting an absent selector to monodomain while the runtime requires the key is a materially different policy and should remain only if explicitly documented as utility-specific.

**Proposed source of truth**

Add one small project-owned resolver in the electroModels configuration layer (for example, an `electroConfiguration` helper near `electroModel`) that returns the selected type and active coefficient dictionary under an explicit compatibility policy. Production runtime selection should remain strict. Utilities should call the same resolver with an intentional compatibility mode rather than implementing precedence themselves. If legacy `electroModel` or flat dictionaries are retained, encode and test their precedence there and issue a consistent diagnostic.

**Migration risk**

Medium. Consolidation can expose old cases that currently work only through one utility's private fallback, and making utility parsing strict too early would be a compatibility break. First capture the accepted legacy matrix, then migrate callers without changing policy; deprecation/removal should be separate.

**Tests required before consolidation**

1. Table-driven C++ tests for current `myocardiumSolver + <type>Coeffs`, legacy `electroModel + <type>Coeffs`, already-selected coefficient dictionaries, flat dictionaries, missing selector, and missing coefficient dictionary.
2. Contract tests asserting production construction remains strict while each utility's deliberately supported compatibility modes are explicit.
3. Integration smoke tests for `listCellModelsVariables`, `ionicHeterogeneityProbe`, and `recomputePseudoECG` against the same current-format case.
4. Negative tests asserting identical errors/diagnostics for malformed selectors.

### AD-2 — eikonal ECG independently parses ionic heterogeneity and has already fallen behind the canonical modes

**Severity:** S2 — Medium  
**Confidence:** high

**Current duplicate locations**

- Canonical parsing/validation and weighting: `src/ionicModels/ionicModel/ionicHeterogeneityOrchestrator.C:138-238`, with mode dispatch at `:158-187`, and shared primitives in `src/ionicModels/ionicModel/ionicHeterogeneity.C:104-158`.
- ECG-side copy: `src/electroModels/ecgModels/eikonalECG/eikonalECG.C:231-331`; it independently locates `ionicHeterogeneity` (`:243-256`), reads defaults (`:268-310`), loads the anatomical field (`:271-297`), and loops over cells (`:317-331`).

**Concrete evidence of drift**

The canonical orchestrator recognizes `transmuralBands`, `namedRegions`, and `cellZoneRegions` (`ionicHeterogeneityOrchestrator.C:158-187`). The ECG copy reads `mode` at `eikonalECG.C:268` but never dispatches on it; it unconditionally reads the three-band interfaces and calls `transmuralBandWeights` at `:307-325`. Thus a valid `namedRegions` configuration is silently interpreted using default three-band boundaries, and `cellZoneRegions` is treated as a field-based three-band configuration. This can make the ECG reconstruction template inconsistent with the ionic heterogeneity used by the myocardium.

**Meaningful differences and intentionality**

The consumer outputs differ intentionally: the ionic orchestrator creates per-cell constants/states, whereas eikonal ECG needs region weights for template blending. The duplicated *configuration interpretation, field acquisition, mode dispatch, and weight construction* are not inherently consumer-specific. Restricting eikonal ECG to three anatomical templates might be an intentional capability limit, but silently accepting other canonical modes as if they were `transmuralBands` is not a safe expression of that limit.

**Proposed source of truth**

Keep `ionicHeterogeneity`/`ionicHeterogeneityOrchestrator` as the owner of heterogeneity schema and weight semantics. Expose a narrow, consumer-neutral API that parses a heterogeneity dictionary and returns validated per-cell named weights (or an explicit unsupported-mode result). Have both ionic constant/state configuration and eikonal ECG consume it. If eikonal ECG supports only exactly three anatomical baselines, validate that constraint explicitly after shared parsing and fail clearly for incompatible named/cell-zone configurations.

**Migration risk**

Medium to high because weights affect scientific output. The refactor must preserve current `transmuralBands` values bit-for-bit or within a stated floating-point tolerance. Adding explicit rejection for configurations currently silently misinterpreted changes behavior, but changes it from incorrect/ambiguous to diagnosed.

**Tests required before consolidation**

1. Golden equivalence tests for hard and blended `transmuralBands`, including exact boundaries, zero transition width, and out-of-range field values.
2. Shared-parser tests for `namedRegions` and `cellZoneRegions`, including gaps, overlaps, duplicate zones, bad baselines, and unsupported smoothing.
3. An integration test that compares ionic-region weights with eikonal ECG template weights cell-by-cell for the same dictionary and field.
4. Negative tests proving eikonal ECG rejects modes/baselines it cannot represent instead of silently applying three-band defaults.
5. A representative ECG regression with an explicit tolerance before and after migration.

### AD-3 — CUDA device selection and fallback policy is copied across all twelve batched model wrappers

**Severity:** S2 — Medium  
**Confidence:** high

**Current duplicate locations**

The same `cudaGetDeviceCount` / `rank % nDevices` / `cudaSetDevice` / `useDevice_` constructor policy appears in:

- `src/ionicModels/AlievPanfilovBatched/AlievPanfilovBatched.C:123-135`
- `src/ionicModels/BuenoOrovioBatched/BuenoOrovioBatched.C:182-194`
- `src/ionicModels/CourtemancheBatched/CourtemancheBatched.C:151-163`
- `src/ionicModels/FabbriBatched/FabbriBatched.C:157-169`
- `src/ionicModels/GaurBatched/GaurBatched.C:163-173`
- `src/ionicModels/GrandiBatched/GrandiBatched.C:153-165`
- `src/ionicModels/PerisYagueBatched/PerisYagueBatched.C:153-165`
- `src/ionicModels/StewartBatched/StewartBatched.C:151-163`
- `src/ionicModels/TNNPBatched/TNNPBatched.C:168-185`
- `src/ionicModels/ToRORd_dynClBatched/ToRORd_dynClBatched.C:202-214`
- `src/ionicModels/TrovatoBatched/TrovatoBatched.C:179-191`
- `src/ionicModels/TWorldBatched/TWorldBatched.C:223-235`

The corresponding `useDevice_` state is also repeated in every model header (for example `AlievPanfilovBatched.H:59`, `TNNPBatched.H:61`, and `ToRORd_dynClBatched.H:67`) even though all models inherit `configuredBatchedIonicModel`.

**Concrete evidence of drift**

Only TNNP emits a warning when CUDA is compiled in but unavailable (`TNNPBatched.C:180-185`). The other eleven copies silently fall back to the host path. Formatting and local `devId` handling also differ, but those are cosmetic; the user-visible fallback diagnostic is a real policy divergence. Device/rank mapping and CUDA-error behavior are correctness-sensitive runtime policy, not ionic-model mathematics.

**Meaningful differences and intentionality**

Kernel launchers, support-buffer sizes, state synchronization, and model-specific memory requirements are intentional differences and should remain in model wrappers. The initial device discovery and rank mapping have no demonstrated model-specific need. TNNP's unique warning appears to be incremental drift rather than a scientifically meaningful exception.

**Proposed source of truth**

Move only device discovery/selection and its diagnostic policy into the batched CUDA infrastructure (`configuredBatchedIonicModel`, `batchedIonicModel`, or a small helper beside `batchedIonicCoreCUDA.H`). Return a selected-device result that wrappers can use; do not attempt to genericize model kernels or allocation layouts. Keep fallback policy configurable only if the repository needs a strict “GPU required” mode.

**Migration risk**

Medium. CUDA initialization order, per-rank device binding, and compilation without CUDA must remain unchanged. Moving ownership of `useDevice_` can affect destruction and synchronization if done too broadly; the first migration should centralize selection while leaving model-specific CUDA buffers and cleanup intact.

**Tests required before consolidation**

1. Compile both with and without `CARDIAC_ENABLE_CUDA`/`HAS_CUDA`.
2. Runtime tests for zero visible devices, one device, multiple devices, and more MPI ranks than devices; assert consistent mapping and diagnostics for every batched model.
3. Inject/observe `cudaGetDeviceCount` and `cudaSetDevice` failures and confirm the chosen fallback/fatal policy.
4. CPU-versus-GPU numerical equivalence smoke tests for at least one simple and one high-state-count model, with stated tolerances.
5. Lifecycle tests that construct/destroy multiple model types sequentially and verify no stale device or buffer state.

## Architecture observations not promoted to findings

- Scalar and batched ionic equation implementations are necessarily different data layouts and many equation headers are generated. A broad scalar/batched unification would carry high scientific and accelerator risk without a demonstrated benefit; existing shared metadata, heterogeneity orchestrator, and batched base classes are the appropriate consolidation level.
- The three myocardium solvers repeat some OpenFOAM solver scaffolding, but their equations, fields, and solution algorithms are materially different. No consolidation is proposed without a specific drift case.
- `Allwmake` and `src/Allwmake` both source `etc/resolveSolids4Foam.sh`, but the latter is a defensive standalone-entry behavior guarded by `_SOLIDS4FOAM_RESOLVED`; the resolver is already the single policy source.
- Repeated runtime-registration macros and `*Coeffs` naming are framework conventions, not inappropriate duplication.

## Priority

1. **AD-1** first: it is an already-broken current-format utility contract and the shared resolver would reduce future selector drift.
2. **AD-2** next, with scientific equivalence tests: it can silently make ECG output inconsistent with configured ionic regions.
3. **AD-3** as a focused CUDA infrastructure change after device-policy tests exist; avoid bundling it with model equation changes.
