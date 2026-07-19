# Agent 10 — Hostile Peer Review

Audit target: `f6b798807e8b8c7b685766d6928c16dc5997e150` on
`no-frontend-minor-errors`. This review covers every proposed finding in the
nine discovery reports. Production code was not edited. Evidence was checked
against tracked source; pre-existing working-tree changes and untracked files
were not used as sole evidence.

## Disposition summary

There are 30 report entries: 12 accepted at their reported severity, 8 accepted
with reduced severity, 3 require investigation, 3 are rejected, and 4 are
duplicates. The apparent count overstates the repair set: build failure masking,
electro-dictionary selector drift, and cellML test coverage were each reported
more than once.

| Finding | Disposition | Final severity | Hostile-review rationale |
|---|---|---:|---|
| OFCPP-01 | Accepted with reduced severity | S2 | Directly confirmed: the utility recognizes only `electroModel` and then passes a canonical current dictionary root to `ionicModel::New`. This breaks a public utility, but not the solver or scientific output, so S1 is exaggerated. |
| OFCPP-02 | Accepted | S1 | `set -e` does not propagate the left side of the two `tee` pipelines; grep is not an exit-status oracle. CI can certify a partial/stale build. |
| OFCPP-03 | Accepted | S2 | Explicit dependency selection checks a source header while auto-discovery checks `lnInclude`; an initialized, unbuilt tree is therefore misclassified as full mode. |
| OFCPP-04 | Investigation required | S3 candidate | The `const_cast` is real, but all current state providers may intentionally expose mutable storage through a read-only metadata interface. No current invalid storage or missed synchronization was shown. Do not add a broad virtual API until model capability and GPU synchronization requirements are established. |
| NC-1 | Investigation required | S1 candidate | Nested calls to `pimple.loop()` on the same object are directly present and highly suspicious. The report assumes, rather than dynamically demonstrates, exact OpenFOAM loop-state behavior and intended corrector ownership. Any fix changes solver ordering and numerical results; require a call-count reproduction on supported OpenFOAM versions before approval. |
| NC-2 | Accepted | S1 | Scalar preconditioning explicitly prevents a spurious global tension transient, the base hook is a no-op, and the batched class does not override it. This is a demonstrated paired-model initialization contract gap with numerical consequences. |
| NC-3 | Investigation required | S2 candidate | The scalar clamp and batched absence are certain, but the intended scientific policy is not. Adding or removing the clamp changes results. Measure realistic rates and obtain maintainer approval before harmonization. |
| PY-01 | Accepted | S1 | `--outdir` is parsed but never passed to the pipeline; generated headers land in CWD and `sort_folder` then looks under `outdir`. The advertised non-default path deterministically fails after partial output. |
| PY-02 | Accepted | S2 | Independently confirmed that `run_mapping(output_c)` writes the names header by basename in CWD. This is the lower-level path defect and should be fixed with PY-01, but is a distinct API contract. |
| PY-03 | Accepted | S2 | README and optional extra advertise a dashboard action/module that do not exist. Maintainer must choose removal or restoration; restoration is not justified without source/history evidence. |
| PY-04 | Duplicate of TC-04 | — | Both identify the same missing cellML regression/CI boundary. Keep TC-04 as the control finding; tests belong in the PY-01/PY-02 patch. |
| F1 | Accepted | S2 | README publishes `*Batched` selector strings while registrations and the Python catalogue use `*compactBatched`. Copying documented values fails runtime selection. Do not add aliases casually. |
| F2 | Accepted with reduced severity | S3 | Registration, build, dictionary enum, and tutorial evidence are direct. Missing introspection is misleading, but the report does not demonstrate a failed MMS execution or sampling request, so S2 is too high. |
| F3 | Accepted | S2 | The manifest inserts biological species labels into field names and infers a solid region from any non-single-cell EP solver. This expands strict planning's accept-surface with nonexistent fields. |
| F4 | Rejected as unsupported | — | Registry/display metadata enumerate known tutorials; no established contract says every registered tutorial is runnable in the current build mode. Root/tutorial documentation already states the full-mode limitation. Availability metadata could be a feature, not a confirmed consistency defect. |
| AD-1 | Duplicate of OFCPP-01 | — | The demonstrated defect is the stale `listCellModelsVariables` selector. A shared multi-mode resolver is a broader architectural proposal with migration risk, not required to repair it. Keep narrow utility fixes and separately document supported legacy precedence. |
| AD-2 | Accepted | S2 | `eikonalECG` reads `mode` but never dispatches on it and always computes transmural-band weights, while canonical parsing supports `namedRegions` and `cellZoneRegions`. Silent reinterpretation can alter ECG output. Prefer explicit rejection first; shared abstraction is optional. |
| AD-3 | Accepted with reduced severity | S3 | Twelve CUDA policy copies and the unique TNNP warning are confirmed, but no different mapping, fallback result, or numerical behavior was shown. Centralization is risky and not yet warranted; first align/define diagnostics and test policy. |
| RS-01 | Duplicate of OFCPP-02 | — | Same root `Allwmake` pipeline defect. |
| RS-02 | Accepted | S1 | Fourteen executable cleaners use caller-relative destructive patterns. Absolute invocation from another case can delete that case's times, processors, mesh, and logs. Self-directory anchoring is a small, non-numerical safety fix. |
| RS-03 | Accepted | S2 | Unquoted `cd ${0%/*}` word-splits checkout paths with spaces in both root entrypoints. |
| RS-04 | Accepted | S2 | Root clean ignores subordinate failures and normally exits with the final `find` status, allowing stale build products after reported success. |
| TC-01 | Duplicate of OFCPP-02 | — | Same build-status failure and regression requirement. |
| TC-02 | Accepted with reduced severity | S3 | CI omits tracked non-orthogonal and electromechanical MMS workflows, but aggregate MMS workflows may be intentionally too expensive for every CI row. The tracked `monodomain1D3D` case is also incomplete as a committed runnable case (only setup files are tracked). Add reduced quantitative gates only after runtime budgeting. |
| TC-03 | Accepted | S2 | Any exit 77 is accepted in every mode and final output says all tests passed. This can silently erase the only EM regression even in full mode. |
| TC-04 | Accepted | S2 | No cellML tests or path-triggered CI job exist despite compiled-code generation; PY-01/PY-02 demonstrate the consequence. |
| TC-05 | Accepted | S3 | Real `foamDictionary` mutation tests skip in Python CI and are not run in the sourced OpenFOAM matrix. This is a genuine integration gap, appropriately low severity. |
| DOC-001 | Accepted with reduced severity | S3 | README reverses the opt-in write policy, advertises nonexistent `-noScale`, omits valid options, and gives stale thresholds. It is serious documentation drift, but executable behavior remains safe. |
| DOC-002 | Accepted with reduced severity | S3 | `-nSteps` and `-deltaT` are not registered; `controlDict` owns both. The defect blocks README-derived commands but does not affect existing execution. |
| DOC-003 | Accepted | S3 | Exact committed tutorial paths disagree with both indexes and the Purkinje README. |

## Deduplication map

- **Build status:** OFCPP-02 is canonical; RS-01 and TC-01 are duplicates.
- **Electro dictionary selection:** OFCPP-01 is the concrete defect; AD-1 is a
  duplicate plus an unapproved broad consolidation proposal.
- **cellML test boundary:** TC-04 is canonical; PY-04 is a duplicate. PY-01 and
  PY-02 remain distinct behavioral defects addressed by the same patch set.
- **Related but not duplicate:** TC-02 is general MMS coverage; NC-1 needs a
  specific loop-ownership reproduction. F1 and DOC-003 are separate
  documentation contracts. F2 and F3 affect different catalogue/manifest data.

## Minimal approved repair set

1. **Build/clean safety:** add pipeline-status propagation to root `Allwmake`;
   propagate component build failures; quote root self-directory changes;
   propagate `Allwclean` failures; anchor the 14 affected tutorial cleaners to
   their own directory. Do not reformat or rewrite deletion patterns.
2. **Current dictionary utility:** update `listCellModelsVariables` to read
   `myocardiumSolver` plus `<type>Coeffs`, with an explicitly tested legacy path
   only if compatibility is desired. Do not introduce a repository-wide
   resolver in this patch.
3. **Dependency resolver:** require the same built-tree `lnInclude` predicate
   for explicit and auto-discovered solids4foam paths.
4. **cellML paths:** pass `outdir` through `run_pipeline`; make both generated
   headers siblings of the requested output; add temporary-directory tests and
   a focused CI trigger. Avoid broader Python modernization.
5. **Driver contracts:** remove species labels from `samplable_fields`; expose
   solid fields only from positive electromechanical context; correct batched
   selector documentation; add manufactured active-tension catalogue metadata.
   Resolve the dashboard by removing stale docs/extra unless restoration is
   separately approved.
6. **Eikonal heterogeneity:** immediately reject unsupported modes clearly.
   Preserve `transmuralBands` bit-for-bit. A shared weight abstraction is a
   later change only if equivalence tests justify it.
7. **Batched Land–Niederer preconditioning:** implement only after a regression
   captures scalar resting-state behavior and an approved tolerance. Treat
   reference changes as numerical baseline changes.
8. **CI skip semantics and docs:** encode expected skips by build mode; fail on
   unexpected exit 77; run the two real dictionary-mutation tests in a sourced
   OpenFOAM job; correct the three documented CLI/path surfaces.

## Required tests and numerical/public-contract cautions

| Area | Minimum validation |
|---|---|
| Build and clean | Stub child exit with token-free output; assert nonzero and no later stage. Test root paths with spaces. Invoke every affected cleaner by absolute path from a sentinel CWD. Test explicit built/unbuilt dependency fixtures. Then clean full and forced-lightweight builds. |
| Utility selection | Current monodomain and single-cell dictionaries; explicitly supported legacy form; missing selector/coefficients negatives; build on oldest/newest supported OpenFOAM. This changes diagnostic/public configuration handling, not solver numerics. |
| cellML | Neutral CWD plus nested `--outdir`; both headers co-located; include resolves; no leaked CWD files; subprocess failure remains nonzero; CI path trigger. This is a public CLI/API correction. |
| Manifest/catalogue/docs | Negative assertions that species/region labels and plain-EP `Ta`/`lambda` are absent; positive full-EM case; registration-to-catalogue and README runtime-name set equality; dashboard packaging smoke test matching the chosen support decision. |
| Eikonal heterogeneity | Golden `transmuralBands` weights at boundaries and transition widths with zero behavioral tolerance where feasible; explicit failures for `namedRegions`/`cellZoneRegions` until supported. Any later support requires ECG reference comparison. |
| Land–Niederer | Zero-stimulus resting-Cai scalar/batched state and tension comparison; multi-cell uniformity; CPU/OpenMP/CUDA parity; approved new EM references. This explicitly changes numerical initialization. |
| Regression runner | Synthetic expected and unexpected exit-77 cases in full/lightweight modes; summaries must distinguish expected skips. Run real dictionary mutation tests with a zero-skip assertion. |
| Documentation | Compare built `-help` output with READMEs; boundary meshes for `checkMeshGeometry`; `controlDict`-driven Purkinje run; verify every indexed tutorial path exists. |

## Deferred or rejected changes

- Do not change PIMPLE loop ownership until call counts reproduce the defect on
  supported OpenFOAM releases. If confirmed, approve it as a numerical S1 with
  coupled 1D–3D, implicit monodomain/bidomain MMS, serial/decomposed, and
  reference-tolerance review.
- Do not harmonize the Land–Niederer stretch-rate clamp until the intended
  scientific policy and realistic rate distribution are established.
- Do not add a mutable ionic-state API merely to remove a `const_cast` until a
  current incompatible model or synchronization contract is identified.
- Do not centralize all electro dictionary or CUDA policy logic as part of the
  narrow repairs. Tests and an explicit compatibility/fallback policy must
  precede those refactors.
- Do not add build-mode metadata to tutorial registry/schema solely from F4;
  treat it as a separately scoped product feature if maintainers want
  environment-aware discovery.

## Validation performed for this review

Static source inspection confirmed each cited control flow, registration,
dictionary lookup, output path, workflow condition, and documentation mismatch.
The repository commit and dirty state were recorded. No production file was
edited, no OpenFOAM build was claimed, and no numerical reference was changed.
The environment did not establish a sourced OpenFOAM runtime, so the numerical
and cross-version checks above remain mandatory for repair agents/CI.
