# driverFOAM architecture pass: sweep, execution, and state unification

**Date:** 2026-08-19
**Status:** Design, not approved, not implemented
**Relationship:** companion to
`2026-08-19-driverfoam-plugin-seam-documentation-design.md`, which stays
deliberately narrow. This document holds the structural findings that spec
excludes.
**Roadmap position:** Phase 2 is complete; these are Phase 3-adjacent. None is
a blocker for the decoupling-evidence claim.

---

## 0. Provenance of these findings

Every claim below was verified against the code on 2026-08-19. Where an
external review's claim proved stale, that is recorded, because a stale finding
repeated confidently is worse than no finding.

**Verified stale — already fixed, no action needed:**

| claim | actual state |
|---|---|
| Artifact predictor conflates unknown with known-empty | Fixed. `detection.py` returns the tuple directly; `None` only when no `export` block matches. `artifacts_predictor.py:23` guards with `is not None`. Regression test added 2026-08-19 (`TestEmptyExportListIsKnownEmpty`), mutation-checked. |
| `_replace_blockmesh_resolution` duplicated across 5 tutorial specs | Fixed, and better than proposed. One definition, in `specs/mesh_provisioning.py`. Three of the five named files no longer exist; the two survivors have zero raw `hex (...)` edits. All 18 tutorials moved to `plugins/cardiacfoam/tutorials/`. |
| RunDocument remains cardiac-shaped, blocks neutral config | Fixed. `schemas/run-document.json` `config` is `additionalProperties: true`, "Core imposes no required keys or shape." This is roadmap P2.2 (P0). |

**Verified true — the subject of this document:** §1, §2, §3, §4.

---

## 1. Two sweep paths with opposite isolation semantics

### 1.1 Finding

`core/runtime/sweep_runner.py` documents this against itself:

> "Entry-based sweeps target an existing registered tutorial whose
> `apply_case()`/`build_cases()` mutate that tutorial's own shared `case_root`
> in place … rather than writing an isolated per-case directory the way
> `build_and_launch` does for generic case_folder sweeps."

So one CLI verb dispatches to two materialization models:

| | entry-based | case_folder |
|---|---|---|
| seam | `TutorialSpec.build_cases` / `apply_case` | `SweepMaterializerCapability.route` / `materialize` |
| workspace | one shared `case_root`, mutated in place | isolated directory per case |
| isolation | none between cases | full |
| reset between cases | `_clean_stale_time_directories` only | fresh directory |

### 1.2 What is *not* wrong

A common misreading is that `route`/`materialize` are redundant hooks that
should be replaced by a pure/impure `plan`/`apply` pair. They already are that
pair. Verified by AST scan of `plugins/cardiacfoam/sweep.py`:

```
route_case_values  write-ops: NONE      -> pure, returns a routed dict
materialize_case   write-ops: write_text, chmod, build_and_launch
```

Core also already owns axis expansion: `sweep_expansion.py` is "Pure axis
expansion for parameter sweeps (cross product / zip)", with `compute_case_count`,
zip-length validation, and `expand_sweep`. The `ParameterSpace → RunSpec →
Materializer` layering exists. Do not "introduce" it.

The redundancy is the **third** surface — `TutorialSpec.build_cases`/
`apply_case` — not the second.

### 1.3 The real hazard: cross-invocation state leak

Because entry-based cases share one `case_root`, correctness requires every
`apply_case` to be **total over the keys it mutates**: a key set on one path
must be set or removed on every other path. Nothing enforces, documents, or
tests this invariant.

It has already been violated and hand-patched.
`plugins/cardiacfoam/tutorials/manufactured_bath_bidomain.py:337`:

> "Symmetric cleanup: a prior electrodePair case sharing this `case_root` may
> have left `surfaceCurrentPatches.xMin` behind, which would collide with
> `groundPatches.xMin` below the same way."

Scope of exposure, verified:

- Sweeps are **strictly sequential** — no concurrency primitives in
  `sweep_runner.py`, no parallel/jobs flag in the CLI. This is *not* a race.
- The audited conditional mutations (`manufactured_bath_bidomain.py:315,353`,
  `cable_1d_cv_convergence.py:119`, `manufactured_monodomain_pseudo_ecg.py:255,259`)
  all branch on **spec-construction parameters** (`ecg_enabled`,
  `fda_bath_variant`, `conductivity`, `electro_properties_scope`), never on
  `case.params`. Within one sweep the branch is therefore constant, so no
  within-sweep leak exists today.
- The leak is **cross-invocation**: a `case_root` left by a previous run with
  different flags carries keys the current `apply_case` never overwrites.

So: real, bounded, currently latent, and mitigated by per-tutorial vigilance in
exactly one place. The risk is that the invariant is invisible — the next
conditional mutation someone adds will not come with a cleanup branch.

### 1.4 Options

**A. Enforce the invariant, keep both paths.** Add a helper expressing
"set this key or remove it," and a test that materializes a tutorial twice under
differing spec parameters and asserts the resulting dictionaries are identical
to a from-clean materialization. Cheapest; leaves the duplication.

**B. Give entry-based sweeps isolated workspaces.** Materialize each entry case
into its own directory, as case_folder sweeps already do. Removes the invariant
entirely rather than policing it. Cost: tutorials assume a stable `case_root`;
output collection and resume both key off it.

**C. Collapse to one path.** Entry-based becomes a `Materializer` over the same
`route`/`materialize` seam; `build_cases`/`apply_case` become compatibility
wrappers and then go. Correct end state, largest change, touches all 18
tutorials.

**Recommendation: A now, C as the eventual target, B only as part of C.** A is
days and removes the live risk. B without C creates a third semantics. C should
not start before the state-file convergence in §3, or it will be done twice.

---

## 2. `run_case` as a second execution pathway — DONE 2026-08-20 (63a1c6e6)

Removed entirely. The framing below was wrong in one way worth recording: it
called `run_case` "live, not vestigial: set by 12 tutorials", which measured
*assignments*, not calls. Zero production invocations existed; only two tests
called it. So the proposed audit ("what does each do that a DAG step cannot
express?") assumed load-bearing behaviour that was not there. Net -2210 lines.

### Original framing

`TutorialSpec` carries `build_cases`, `apply_case`, `run_case`,
`collect_outputs`. `run_case` is live, not vestigial: set by 12 tutorials plus
`core/runtime/generic_case.py:350`.

Once a case is materialized, execution should be a DAG. A cardiac tutorial, a
solids4foam case, and any future non-OpenFOAM program should reach one executor
through their workflow DAG. A parallel `run_case` path means every execution
concern — resume, retry, telemetry, provenance, artifact capture — either gets
implemented twice or silently works on only one path.

**Proposal:** treat `run_case` as deprecated. Audit what each of the 12 does
that a DAG step cannot express; that list is the actual work item. Do not remove
it before that audit — the reasons it exists are undocumented and some are
likely legitimate (parallel decomposition, multi-stage solves).

---

## 3. Two runtime state worlds — DONE 2026-08-20 (5cfb8f46)

Also under-described below. This was not two worlds; it was one world plus a
ghost that `describe` advertised to agents as the source of truth, with
polling guidance, while nothing had ever written it — and six stale copies
were committed under `tutorials/`, so a compliant agent could read a fossil
and report a run that never happened. Retired in favour of
`workflow_state.json`.

### Original framing

`workflow_state.json` is read by 8 core modules (`cli.py`, `introspection.py`,
`postprocess_phase.py`, `workflow_orchestrator.py`, `fresh.py`,
`output_collection.py`, `sweep_runner.py`, `artifacts.py`).
`run_manifest.json` is read by 4 (`introspection.py`, `run_discovery.py`,
`execution_context.py`, `artifacts.py`).

`introspection.py` and `artifacts.py` read **both**.

This is the highest-leverage item here, and the reason is compounding: every new
capability risks acquiring a strict version and a legacy version. The plugin
capability work, `run_case` retirement, and sweep unification all touch state,
so each one done before convergence is done twice.

**Proposal:** before further capability work, decide whether `run_manifest.json`
is derived from `workflow_state.json` or retired. Write that decision down as a
migration note with a schema test, in the manner of P2.3.

---

## 4. Residual cardiac vocabulary in core

`openfoam_driver/specs/fixtures/template/constant/electroProperties` is a
cardiac dictionary inside core `specs/`, read only by
`tests/plugins/cardiacfoam/test_template_contract.py`. It belongs in
`plugins/cardiacfoam/`.

Trivial, but it is exactly the P2.5 acceptance criterion: "Generic packages
contain no cardiac dictionary vocabulary outside versioned
compatibility/migration code."

---

## 4.5 FIXED 2026-08-19: manufacturedBathBidomain could not apply its own case

**Status: FIXED. Three stacked faults, each masked by the previous one.**
Regression test: `tests/plugins/cardiacfoam/test_bath_bidomain_variant_apply.py`
(4 tests, covering both variants, the default path, and round-trip idempotence).

| # | fault | fix |
|---|---|---|
| 1 | block remover (`remove_electro_property_dict`) used to delete the scalar `xMin` | added `remove_foam_entry` to `mutators.py` + `remove_electro_property_entry` wrapper |
| 2 | `manufacturedBidomain.fdaBathVariant` over-deleted from the hex template by `91a3debb` | restored (the tet template kept it, which is why tet tests passed) |
| 3 | `groundPatches {}` empty, but `update_foam_entry` updates rather than inserts | bath patch entries now upsert via `ensure_electro_property_entry` |

Faults 1 and 3 are one gap: the mutation API was block-shaped
(`remove_foam_dict`, `ensure_foam_dict`) while this tutorial needed
entry-shaped remove *and* upsert. The strict "key must exist" default is
preserved everywhere else — only the bath patch entries upsert, because those
are the keys whose presence legitimately varies by variant.

Fault 2 was silent rather than loud in C++: `bathECGManufacturedVerifier.C:233`
uses `lookupOrDefault(..., "groundElectrode")`, so a missing key would have made
the ECG verifier read `groundElectrode` while `verificationModel` said
`electrodePair` — the two halves of the verification disagreeing without any
error.

**Still open (not a bug, a decision):** `make_spec` defaults to
`groundElectrode` while the committed template is in the `electrodePair` state.
Applying now works either way, but whoever runs the tutorial with defaults
rewrites the tracked file — the shared-`case_root` drift of §1.3.

**Original finding, retained for context:**

```
registered spec, defaults; 12 cases; first = 1D_10_cells_implicit_DT0p00892857
  apply_case -> KeyError: "Dictionary 'xMin' has no opening brace"
```

`make_spec()` with no arguments, on a pristine copy of the checked-in
`tutorials/manufacturedSolutions/bathBidomain`. Verified at the spec level
(`make_spec` -> `apply_case`); NOT verified end to end through `foamctl`, so if
some CLI path supplies `electrodePair` the crash stays latent there.

**Nothing is invalid.** `groundElectrode` is a valid C++ variant and the
committed `electroProperties` is valid OpenFOAM. This is a driverFOAM Python
fault, raised before any solver runs. The message means "you asked me to delete
a dictionary named `xMin`, but `xMin` is not a dictionary" -- not "your
dictionary is malformed".

Two faults compound:

1. **Wrong tool for the shape.** Both bath-variant cleanups in
   `plugins/cardiacfoam/tutorials/manufactured_bath_bidomain.py` (lines ~321
   and ~340) call `remove_electro_property_dict` on `xMin`, which is a scalar
   entry (`xMin -0.01;`), not a block. `remove_foam_dict` scans forward for an
   opening brace, finds none, and raises. `missing_ok=True` does not help: it
   wraps only scope resolution, so a key that is *present but the wrong shape*
   passes straight through it.

   This is an API gap, not carelessness. The entire removal surface is
   `remove_foam_dict_via_foamDictionary`, `remove_foam_dict`, and
   `remove_electro_property_dict` -- all block removers. No scalar-entry
   remover exists, so the author reached for the only tool available.

2. **Committed template contradicts the default.** `fda_bath_variant` defaults
   to `groundElectrode` in both `_apply_case` and `make_spec`, but the tracked
   `electroProperties` is in the electrodePair state (`groundPatches {}` empty,
   `surfaceCurrentPatches` holding `xMin`/`xMax`). That asymmetry is the
   shared-`case_root` mechanism's fingerprint: the committed file records
   whichever variant ran last.

**Why no test catches it:** the tet suite passes `mesh_family="tet"` and builds
against `tmp_path`, never exercising the default hex path against the real
tracked case.

**Fix when taken up:** add `remove_foam_entry` beside `remove_foam_dict`
(delete the matched line when it has no opening brace instead of raising),
point both cleanups at it, and decide whether the committed template should
become `groundElectrode` or the default should become `electrodePair` -- they
currently disagree, which is what makes the cleanup fire at all.

**This is the concrete argument for §1.4 option A.** The proposed invariant
test -- materialize a tutorial twice under differing spec parameters, assert
the dictionaries match a from-clean materialization -- fails on today's code
for exactly this reason. It is not hypothetical hygiene.

---

## 5. Standing design rules extracted from this pass

Two rules earned by real defects, worth applying to new seams:

1. **UNKNOWN ≠ KNOWN_EMPTY.** Never let one optional value carry both "I could
   not determine this" and "I determined this, and it is empty." The artifact
   predictor cost a real debugging session to this. Audit any new
   `X | None` return that a consumer branches on.
2. **Shared mutable workspace requires total mutation.** If two operations write
   the same workspace, every key one sets, the other must set or remove.
   Prefer isolation over policing; if isolation is impractical, make totality
   explicit and tested.

---

## 6. Sequencing

```
converge state files (§3)
   -> enforce mutation invariant (§1.4 option A)
   -> audit and retire run_case (§2)
   -> collapse sweep paths (§1.4 option C)

independent, any time: move the cardiac fixture (§4)
```

§3 first because it is the one everything else touches. §4 can be done in a
minute by anyone.

---

## 7. Explicit non-goals

- Introducing `StudySpec`/`ParameterSpace`/`RunSpec` as new named types. The
  concepts exist under different names; renaming is churn without a consumer.
- A separate post-processing orchestration system. Post-processing is already a
  DAG step and should stay one.
- Any change justified only by an out-of-tree plugin that does not yet exist.
  Per the roadmap's own gate: build the solids4foam plugin, then classify each
  core change it forces as a genuine missing seam or a project-specific
  exception.
