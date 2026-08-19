# Publication Roadmap

The roadmap is ordered by evidentiary dependency. Do not begin the expensive
agent study before the deterministic software and numerical baselines are
stable; otherwise failures cannot be attributed to the agent or the driver.

Effort estimates are person-weeks of focused work and are deliberately rough.

## Roadmap summary

| Phase | Objective | Key outputs | Exit gate | Effort |
|---|---|---|---|---:|
| 0 | Freeze scope and claims | Paper charter, claim/evidence matrix | Authors agree on cardiac software-paper scope and prohibited overclaims | 0.5 |
| 1 | Make the artifact clean | Green CI, reviewed fixture drift, fixed paths, release checklist | Zero unexplained failures/skips in required matrix | 1–2 |
| 2 | Stabilize generic contracts | Explicit readiness semantics, generic RunDocument design, compatibility telemetry | Non-cardiac v2 context uses no cardiac fallback | 2–5 |
| 3 | Produce deterministic evidence | Numerical equivalence, fault corpus, recovery, overhead, provenance results | E1–E4 and E7 complete with frozen analysis | 4–8 |
| 4 | Demonstrate OpenFOAM generality | Independent solids4foam plugin and case ladder | E6 passes with zero core changes | 3–6 |
| 5 | Demonstrate reproducibility | Tagged artifacts, container/env recipe, blinded external replay | E5 passes and deviations are documented | 1–3 |
| 6 | Optional human/agent study | User or LLM benchmark with baselines and ablations | E8/E9 completed under preregistered protocol | 4–10 |
| 7 | Write, release, and submit | Manuscript, archive DOI, data/code availability, response-ready supplement | Internal skeptical review passes | 3–5 |

## Phase 0 — Scope and claim freeze

| ID | Task | Deliverable | Completion test |
|---|---|---|---|
| P0.1 | Select near-term paper type | One-page paper charter | States venue class, audience, primary RQ, primary endpoint, and excluded claims |
| P0.2 | Adopt terminology | Glossary in manuscript repository | Separates plan-valid, execution-ready, numerically verified, scientifically validated |
| P0.3 | Freeze architecture baseline | Tagged commit and environment manifest | All reported code/test counts resolve to the tag |
| P0.4 | Assign evidence owners | Owner column added to this roadmap | Every required experiment and artifact has one accountable owner |

Recommended decision: write the cardiac/OpenFOAM research-software paper first.
Treat agent evaluation as optional or a later paper unless resources permit a
proper controlled benchmark.

## Phase 1 — Clean and releasable software

| ID | Priority | Task | Acceptance criterion |
|---|---:|---|---|
| P1.1 | P0 | Review the failing tutorial-characterization digests | Explain every changed section; update fixture only when changes are intended; full test passes |
| P1.2 | P0 | Fix stale `buildAndTest.yml` test paths | Binary CI reaches and passes the intended real-`foamDictionary` tests |
| P1.3 | P0 | Classify all 11 local skips | Required solver/OpenFOAM tests run in a suitable CI job; optional skips have documented reasons |
| P1.4 | P0 | Run driver regression-equivalence phase 2 in built-solver CI | Every supported mapped case is executed through driverFOAM and compared to committed reference values |
| P1.5 | P0 | Make CI status unambiguous | Required PR and release workflows are green on the tagged baseline |
| P1.6 | P1 | Add coverage reporting | Publish line/branch coverage by core, cardiac plugin, and integration layer; do not use one aggregate percentage alone |
| P1.7 | P1 | Complete release metadata | License, citation file, changelog, contribution guide, code of conduct/governance as appropriate, version policy |
| P1.8 | P1 | Build clean-install tests | Wheel/sdist install in a new environment; console scripts, schemas, profiles, and entry points work without editable checkout |

## Phase 2 — Contract and architecture hardening

| ID | Priority | Task | Acceptance criterion |
|---|---:|---|---|
| P2.1 | P0 | Resolve plan/readiness status semantics | Truth-table tests cover structural, environment, warning, and execution states; all CLI paths use one launch predicate |
| P2.2 | P0 | Design a generic RunDocument version | Core `run-document.json` defines config as an open object (`additionalProperties: true`). Core engine dynamically evaluates `jsonschema.validate()` against a plugin-provided schema, yielding structured `StrictDiagnostic`s on failure to enable agent self-healing. |
| P2.3 | P0 | Provide v2→new-version migration | Existing archived cardiac documents migrate deterministically with schema tests |
| P2.4 | P0 | Instrument compatibility fallbacks | Tests can assert which fallback was called; explicit non-cardiac v2 contexts call none |
| P2.5 | P1 | Move cardiac detection/override semantics into plugin | Generic packages contain no cardiac dictionary vocabulary outside versioned compatibility/migration code |
| P2.6 | P1 | Generalize legacy generic case mutation API | Generic dictionary overrides are path/schema based; old cardiac names remain deprecated aliases only |
| P2.7 | P1 | Make audit text plugin-neutral | Non-cardiac reports contain no electro/ionic/cardiac terms |
| P2.8 | P1 | Remove duplicate schema drift risk | One schema source of truth; any second copy is generated and checked |
| P2.9 | P1 | Test trust boundary end to end | RunDocument ingestion, command/cwd/path checks, case-script caveat, and allowed-root behavior are regression gated |

Phase 2.2 is mandatory before claiming a fully project-neutral configuration
model. It is not mandatory for a narrowly scoped cardiacFoam paper if the
limitation is stated explicitly.

## Phase 3 — Deterministic experiments

| ID | Experiment | Work package | Exit artifact |
|---|---|---|---|
| P3.1 | E1 | Freeze six representative cases, input trees, observables, and tolerances | `equivalence_protocol.yaml` plus pilot justification |
| P3.2 | E1 | Execute direct, minimal-Python, and driver conditions | Raw runs, comparison CSVs, numerical-equivalence figure |
| P3.3 | E2 | Build fault corpus with independently authored hold-out subset | Versioned mutation corpus and adjudicated labels |
| P3.4 | E2 | Run all validators/baselines | Confusion matrices, localization results, invalid-launch counts |
| P3.5 | E3 | Define interruption/staleness schedule | Repeatable fault-injection harness |
| P3.6 | E3 | Measure recovery and recomputation | Recovery table and state/provenance audit |
| P3.7 | E4 | Implement benchmark harness | Machine-readable timing/memory/output-size records |
| P3.8 | E4 | Run sweep sizes and workload classes | Scaling/overhead plots with uncertainty |
| P3.9 | E7 | Perturb inputs, plugin, workflow, solver binary/library, environment | Provenance sensitivity and resume-decision table |
| P3.10 | All | Archive analysis code and raw data | One command regenerates every manuscript table/figure |

### Existing 16-experiment migration

| Task | Target |
|---|---|
| Convert legacy paths | Move the nine `legacy_bash` experiment executions to driver-backed workflows where scientifically equivalent |
| Preserve native baselines | Keep the original Bash runners as comparison conditions, not silently replace them |
| Materialize results | Produce all 16 declared result CSVs from the frozen code/environment |
| Validate references | Record provenance and rationale for all 13 declared references; explain three experiments without a reference |
| Add acceptance checks | Contract tests must verify result values/criteria, not only paths and JSON structure |
| CI stratification | Small numerical gates on PR; complete matrix on scheduled/release runs |

## Phase 4 — solids4foam external plugin

| ID | Task | Acceptance criterion |
|---|---|---|
| P4.1 | Create separate plugin repository and package | Installs through `driverfoam.plugins` against a tagged driverFOAM release |
| P4.2 | Implement profile, command, case, and dictionary capabilities | Patch-test case describes/plans without cardiac fallbacks or vocabulary |
| P4.3 | Execute case ladder S1–S5 | Native and driver conditions satisfy the same regression/benchmark criteria |
| P4.4 | Add solids4foam fault mutations | Project-specific detection precision/recall is reported |
| P4.5 | Demonstrate sweep and resume | One material/mesh sweep and one interrupted nonlinear/coupled run recover correctly |
| P4.6 | Add OpenFOAM version matrix | Declared supported variants pass; unsupported variants fail clearly |
| P4.7 | Freeze integration-cost evidence | Core changes = 0; plugin LOC, hours, limitations, and unsupported cases reported |

See [`SOLIDS4FOAM_CASE_STUDY.md`](SOLIDS4FOAM_CASE_STUDY.md) for the detailed
case and capability design.

## Phase 5 — Independent reproduction

| ID | Task | Acceptance criterion |
|---|---|---|
| P5.1 | Create tagged source and data releases | Persistent DOI; exact code/data/environment identifiers |
| P5.2 | Provide container or exact environment recipe | Clean machine can install/build without undocumented local paths |
| P5.3 | Select blinded reproducer | Person was not involved in implementation or experiment authoring |
| P5.4 | Run reproduction protocol | Reproducer obtains expected outputs within tolerance; all deviations logged |
| P5.5 | Incorporate feedback without hiding failures | Final supplement reports initial failures and corrective documentation changes |

Use ACM terminology precisely: same-team/same-setup repeatability is weaker than
independent reproduction and replication [R7].

## Phase 6 — Optional user and agent studies

### User study

- Obtain institutional ethics determination if necessary.
- Use a randomized crossover design.
- Stratify novice and experienced OpenFOAM users.
- Predeclare the primary endpoint, preferably successful completion within a
  fixed time or time to first valid result.
- Retain interaction/event logs rather than relying only on questionnaires.

### Agent study

- Freeze hidden tasks and contamination controls.
- Use the same model, prompt budget, and stopping criteria in raw-shell and
  driverFOAM conditions.
- Include deterministic driverFOAM as a non-agent control.
- Use multiple model families and repeated trials.
- Score execution, numerical correctness, safety, cost, and intervention
  separately.
- Run ablations for strict planning (including JSON schema validation), introspection, provenance/resume, and structured diagnostic repair.

## Phase 7 — Manuscript and release

| ID | Task | Acceptance criterion |
|---|---|---|
| P7.1 | Write methods before interpreting results | Protocol, exclusions, thresholds, and statistics match frozen plans |
| P7.2 | Generate figures/tables from archived data | No hand-copied result values |
| P7.3 | Add limitations and threat model | Explicitly rejects physical-correctness and sandbox claims |
| P7.4 | Perform internal skeptical review | Reviewer can map every abstract claim to a result/table/test |
| P7.5 | Produce release and DOI | Manuscript cites exact software and data releases |
| P7.6 | Select venue and reformat | Scope and evidence meet venue requirements, not only template formatting |

## Critical path

```text
scope freeze
  → green/reproducible CI
  → driver numerical-equivalence gate
  → frozen deterministic benchmark
  → generic contract decision
  → solids4foam external plugin
  → independent reproduction
  → manuscript/release
  → optional agent study
```

## Stop/go gates

- **After Phase 1:** if the binary driver equivalence cannot be made stable,
  stop publication work and fix execution semantics.
- **After E1:** if driver and native results differ beyond explained numerical
  tolerances, investigate before measuring productivity or agents.
- **After E2 pilot:** if false positives are high, narrow the “strict” claim and
  improve diagnostic calibration before the final benchmark.
- **After solids4foam S1:** if core changes are required, classify each as a
  generic missing seam or project-specific exception; do not hide them inside
  the plugin claim.
- **Before submission:** if independent reproduction fails, publish the failure
  analysis internally and repeat only after the artifact is corrected.

References are listed in [`REFERENCES.md`](REFERENCES.md).

