# cardiacFoam consistency and correctness audit — final report

## Executive summary

The hostile review and post-repair integration review accepted **24 unique
findings**: **5 S1**, **12 S2**, and **7 S3**. There are no accepted S0 or S4 findings. The hostile review's prose
summary says 20 accepted findings, but its per-finding disposition table contains
22; this report follows the table and its explicit deduplication map, plus two
confirmed findings from the continued integration analysis.

The working tree contains narrow repairs or preventive controls for 21 findings.
Three accepted findings are wholly deferred: CF-012 (CUDA policy consolidation),
CF-016 (expanded MMS gates), and CF-024
(PIMPLE loop ownership). CF-011's
unsafe silent interpretation is repaired by explicit rejection, while its broader
shared-parser refactor is intentionally deferred. No generated ionic equations,
constants, tolerances, solver ordering, or numerical references were changed.

The audit began from commit `f6b798807e8b8c7b685766d6928c16dc5997e150`
on `no-frontend-minor-errors`. The checkout was substantially dirty before the
audit, including a modified external `modules/solids4foam` submodule and many
untracked research/tutorial artifacts. Those are pre-existing user changes and
are not audit repairs. The submodule was not edited by the audit.

## Convention map

- Maintained code and build files are authoritative over READMEs and project
  memory. Nearby equivalent components define style; no global reformatting was
  performed.
- `myocardiumSolver` plus `<type>Coeffs` is the current electro-properties
  selection contract. Runtime names, registrations, catalogues, dictionaries,
  and documentation must use identical selector strings.
- Runtime-selectable C++ follows local OpenFOAM ownership, dictionary, fatal
  diagnostic, and registration idioms. Compatibility spans OpenFOAM v2312-v2512.
- Full mode requires a built solids4foam tree; lightweight mode uses the bundled
  physicsModel fallback. Source presence alone does not establish a full install.
- Python CLI arguments must propagate through pipeline APIs, and paired generated
  headers must remain co-located. Public package extras/actions must correspond to
  executable modules.
- Shell entrypoints are callable from any CWD, quote derived paths, fail closed,
  and propagate child status. Destructive cleaners must anchor to their case.
- Large ionic equation headers are generated code; generators/contracts should be
  repaired instead of normalizing outputs. `modules/solids4foam` is external.
- Scientific changes require an explicit hypothesis, tolerances, regression or
  equivalence tests, and maintainer approval before reference changes.

## Validated findings

Ordered by severity; canonical details and exact evidence are in
`validated-findings.json`.

### S1 — High

| ID | Finding | Disposition |
|---|---|---|
| CF-002 | Top-level build masks child failures | Repaired: pipeline and component failure propagation |
| CF-004 | LandNiedererBatched omits resting-Cai preconditioning | Repaired with batched resting-state conditioning; full EM/backend parity remains pending |
| CF-005 | `cellML2foam --outdir` is not propagated | Repaired with temporary-directory tests |
| CF-013 | Tutorial cleaners operate on the caller directory | Repaired by anchoring identified cleaners |
| CF-024 | Nested ownership consumes strong-coupling PIMPLE correctors | **Accepted, deferred** pending numerical regressions |

### S2 — Medium

| ID | Finding | Disposition |
|---|---|---|
| CF-001 | Utility reads obsolete electro selector contract | Repaired; OpenFOAM runtime validation pending |
| CF-003 | Explicit unbuilt solids4foam tree selects full mode | Repaired with consistent built-tree predicate |
| CF-006 | `run_mapping` separates paired outputs | Repaired with output-path tests |
| CF-007 | driverFoam advertises absent dashboard API | Repaired by removing unsupported docs/extra |
| CF-008 | README publishes unregistered batched selectors | Repaired in documentation |
| CF-010 | Capability manifest mixes labels/fields and invents mechanics fields | Repaired with negative/positive tests |
| CF-011 | eikonalECG silently reinterprets heterogeneity modes | Repaired by explicit rejection; shared refactor deferred |
| CF-014 | Root entrypoints fail in whitespace paths | Repaired by quoted fail-closed directory changes |
| CF-015 | Top-level clean masks subordinate failures | Repaired by fail-fast propagation |
| CF-017 | Regression runner accepts arbitrary skips | Repaired with mode-specific expected skips |
| CF-018 | cellML2foam lacks tests and a CI trigger | Repaired with focused tests/workflow |
| CF-023 | Executable scripts have content before their shebang | Repaired by placing `#!` at byte zero |

### S3 — Low

| ID | Finding | Disposition |
|---|---|---|
| CF-009 | ManufacturedElectromechanics absent from introspection | Repaired with catalogue entry/tests |
| CF-012 | CUDA fallback diagnostics diverge across wrappers | **Accepted, deferred** pending policy and GPU tests |
| CF-016 | Scientific MMS workflows absent from CI regression gate | **Accepted, deferred** pending complete fixtures/runtime budget |
| CF-019 | Real foamDictionary mutation tests lack a CI home | CI control added; live sourced run pending |
| CF-020 | checkMeshGeometry README contradicts CLI/safety default | Documentation repaired |
| CF-021 | runPurkinjeGraph documents nonexistent time options | Documentation repaired |
| CF-022 | Tutorial indexes name nonexistent paths | Documentation repaired |

## Rejected, duplicate, and uncertain findings

- OFCPP-04: investigation only. The `const_cast` is real, but no incompatible
  state provider or synchronization failure was demonstrated. A mutable state API
  would be a broad unapproved change.
- NC-1 is no longer uncertain: official OpenFOAM v2512 source confirms that every
  `pimpleControl::loop()` call increments the shared corrector counter. It is now
  accepted as CF-024, but remains unrepaired because changing ownership alters
  solver ordering and requires coupled numerical regressions.
- NC-3: investigation only. Scalar and batched stretch-rate policies differ, but
  the intended scientific policy is unresolved; harmonization could change results.
- F4: rejected. No established contract requires every registered tutorial to be
  runnable in every build mode; environment-aware discovery is a product feature.
- Duplicates: AD-1 folds into CF-001; RS-01 and TC-01 fold into CF-002; PY-04
  folds into CF-018. The proposed repository-wide electro resolver and broad CUDA
  centralization were not approved as parts of narrow fixes.

## Patch summary mapped to finding IDs

- **Build and clean safety — CF-002, CF-003, CF-013, CF-014, CF-015:** root and
  component scripts propagate errors and quote self-directory changes; dependency
  resolution consistently requires a built tree; identified cleaners anchor locally.
- **Executable portability — CF-023:** four tracked scripts now place their
  shebang at byte zero so direct `execve` launchers can run them.
- **C++ configuration contracts — CF-001, CF-011:** the utility handles canonical
  selector/coefficient dictionaries with a labelled legacy path; eikonalECG rejects
  unsupported heterogeneity modes instead of silently changing their meaning.
- **cellML generation — CF-005, CF-006, CF-018:** outdir reaches the pipeline,
  paired headers are siblings, focused tests were added, and a path-triggered CI
  workflow was introduced.
- **Driver contracts — CF-007, CF-009, CF-010, CF-019:** stale dashboard surfaces
  were removed; manufactured tension metadata was added; field capability inference
  was corrected; sourced foamDictionary tests were added to the build matrix.
- **Regression semantics — CF-017:** CI passes build mode and the runner fails
  unexpected exit-77 skips with an accurate summary.
- **Active-tension initialization — CF-004:** `LandNiedererBatched` now applies
  configurable resting-Cai preconditioning through its existing batched hot
  path and resets transient/stretch history before physical time begins.
- **Documentation — CF-008, CF-020, CF-021, CF-022:** runtime names, CLI behavior,
  safety defaults, configuration ownership, and tutorial paths now match code/tree.
- **No production repair — CF-012, CF-016, CF-024:** these remain accepted and
  deferred for the reasons recorded above.

Files changed for unrelated research, tutorial development, generated results,
the external submodule, and pre-existing working-tree work are not included in
this mapping and must not be attributed to the audit.

## Validation performed

Passed locally:

- `python3 -m json.tool .audit/validated-findings.json`
- `python3 -m pytest applications/scripts/cellML2foam/tests/test_output_paths.py -q`
  — **3 passed**
- `python3 -m pytest applications/utilities/listCellModelsVariables/tests/ -q`
  — **3 passed**
- Complete driverFoam validation — **719 passed, 3 skipped, 126 subtests**
- `PYTHONPATH=applications/scripts/driverFoam python3 -m pytest
  applications/scripts/driverFoam/openfoam_driver/tests/test_capability_manifest.py
  applications/scripts/driverFoam/openfoam_driver/tests/test_audit_contract_repairs.py -q`
  — **15 passed**
- Strict driver dictionary-key scanner — clean
- `bash -n` on root/component build and clean scripts, the regression runner,
  and solids4foam resolver
- byte-zero shebang assertions for every tracked script touched by CF-023; an
  equivalent pre-fix fixture reproduced direct-execution errno 8
- `git diff --check`

Unavailable locally because no OpenFOAM environment is sourced and `wmake` and
`foamVersion` are absent:

- full and forced-lightweight builds and cleans;
- v2312/v2412/v2512 compilation and runtime-selection validation;
- utility `-help` and mesh/graph runtime checks;
- sourced foamDictionary integration tests;
- OpenFOAM tutorial, ECG, and electromechanical numerical regressions;
- CUDA compilation/runtime and CPU/OpenMP/CUDA parity.

## Remaining risks

- CF-004 changes the batched initial state as intended, but scalar/batched,
  OpenMP/CUDA, and full electromechanical reference parity still require a
  sourced build and numerical review.
- CF-024 is a confirmed strong-coupling loop defect, while NC-3 remains an
  investigation; neither should be changed without reference-tolerance review.
- CF-011 now fails clearly for unsupported modes, but named/cell-zone ECG template
  support remains absent.
- Full/lightweight, cross-version, destructive-cleaner sentinel, and CUDA checks
  have not run in this host environment; CI or a sourced developer environment
  must supply them.
- CF-016 leaves coupled, non-orthogonal, and electromechanical MMS coverage gaps.
- The dirty checkout makes file-level attribution important. Review audit-mapped
  hunks separately from unrelated user work before committing.

## Demonstrated preventive controls

- A path-filtered cellML2foam workflow exercises paired-output placement.
- Capability tests reject species labels, plain-EP mechanics fields, and invalid
  single-cell solid inference while retaining a positive EM case.
- Build-mode-aware regression skip policy turns unexpected skips into failures.
- The sourced OpenFOAM matrix requires foamDictionary for real mutation tests,
  preventing silent test skips.
- Build/clean orchestration now carries child failure status; destructive cleaners
  anchor to their own cases.
- Runtime/configuration drift is constrained by focused selector-contract and
  driver contract tests. Registration/catalogue/documentation set-equality remains
  a recommended next control where not yet implemented.
