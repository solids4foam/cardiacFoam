# driverFOAM: agent–solver parity, provenance, and agnosticity

Date: 2026-08-10
Status: approved design, not yet implemented

## Problem

driverFOAM's strict contract (`plan --strict` → RunDocument v2 → normalized
workflow DAG → command allowlist → artifact reconciliation) is a sound agent
interface and is not changed by this design. Four gaps sit on top of it.

1. **The agent never sees the solver's numerical state.**
   `workflow_runner.run_workflow_step` executes `subprocess.run` and writes
   stdout to a file that nothing parses. `failure_classification.classify_failure`
   decides retryable-vs-fatal purely from diagnostic codes
   (`RETRYABLE_CODES = {"workflow_step_timeout"}`). Consequences: a solve that
   diverges to NaN but exits 0 is reported as `completed`; a solve that stops
   early at `writeInterval` boundaries is indistinguishable from one that
   reached `endTime`; and the agent's only evidence for a real failure is
   `--tail-lines` of raw log text.

2. **Artifacts assert existence, never validity.**
   Strict planning predicts that a file will appear and reconciles that
   prediction after each step. Nothing asserts the contents are finite, in
   range, or converging at the expected order — which for an MMS-heavy
   repository is the question that matters. `verification_contracts.py` has
   `observables`/`aggregation`, but it is a separate experiments layer and is
   not wired into the strict run path.

3. **Resume has no input provenance.**
   `workflow_state.json` records step status but no digest of the inputs that
   produced it. A `completed` state is replayed as success even when the case
   dictionaries or the solver binary have changed since. This is documented at
   length in `AGENT_GUIDE.md` and has already caused a near-miss: a sweep
   re-run after a solver code change reported the previous day's numbers as
   fresh, caught only because two code versions cannot agree to six
   significant figures. The mitigation today is a human remembering `--fresh`.

4. **Solver-specific knowledge still leaks out of the plugin seam.**
   The `DriverContext` / `PluginCapabilities` work correctly removed the
   process-global active plugin, but several top-level modules still hardcode
   cardiacFoam: `capability_manifest.py` (`_ELECTRO_SOLVER_FIELDS`,
   `_SOLID_SOLVER_FIELDS`, and a literal `constant/electroProperties` path),
   `tutorial_contracts.REQUIRED_FILES`, `introspection.py`'s embedded
   `$ELECTRO_MODEL_COEFFS` schema text, and `utility_catalog.py`'s fixed
   `applications/utilities/` scan root. `plugin_api_version` exists but nothing
   enforces it, and plugin loading is limited to a trusted
   `--plugin module:Class` import.

## Non-goals

- **No C++ changes.** The C++ dictionary read-surface stays a hand-maintained
  Python catalog (`plugins/cardiacfoam/dict_entries_catalog.py`) cross-checked
  by the approximate scanner and its allowlist. The decision is deliberate: the
  Python catalog is updated as the source of truth in this project's workflow.
  Deriving the catalog from the solver (a `-dumpContract` mode, or per-model
  manifest sidecars extending the `utility.manifest.toml` pattern) remains
  available as a future spec but is explicitly out of scope here.
- **No streaming telemetry.** Logs are parsed after a step completes, not
  tailed live. Live progress reporting and early-abort on divergence are a
  later, separate change.
- **No new predicate language.** Observable assertions use a closed
  vocabulary, for the same reason `report_catalog.py`'s `applicable_when` is
  restricted to flat key-equality: a richer v2 must not be able to silently
  re-interpret v1 documents.
- **No change to the public strict loop.** `plan --strict`, `run --strict`,
  `step --strict`, `sweep-plan`, `sweep-run`, and RunDocument execution keep
  their existing shape and exit semantics. Every addition below is additive.

## Sequencing

Four phases, ordered by ascending structural cost, each independently
committable and each leaving the driver in a working state.

Phase 0 is a prerequisite, not part of this design: the branch currently
carries the in-flight plugin-seam refactor uncommitted (25 modified files, 12
untracked, ~750 insertions). Verify the suite against it and commit it on its
own so each later phase has a clean diff and a rollback point.

Phases 1–3 are plugin-free by construction. That is a design property, not an
accident: it keeps them independent of Phase 4's protocol changes, and it
means Phase 4's acceptance test can use them as instrumentation.

---

## Phase 1 — Input-provenance digest

**Goal.** A resumed run whose inputs no longer match the recorded inputs is
refused, rather than replaying a stale `completed` state as fresh.

**New module:** `core/runtime/provenance.py`

```
compute_input_digest(case_root, *, commands) -> InputDigest
```

`InputDigest` is a frozen dataclass with `algorithm`, `digest` (a
`sha256:`-prefixed hex string), and `components` (an ordered mapping of what
went into it, for diagnostics).

Digest inputs, in order:

- Every file under `<case_root>/system/` and `<case_root>/constant/`, walked in
  sorted path order. Files below a 1 MiB threshold contribute their content
  hash. Files at or above it — meshes, imported geometry — contribute
  `(relative_path, size, mtime_ns)` instead. Rationale: content-hashing a
  multi-gigabyte `polyMesh` on every resume is not affordable, and any real
  mesh change moves both size and mtime.
- The resolved executable for each command in the workflow DAG, contributing
  `(resolved_path, size, mtime_ns)`. Binaries are never content-hashed:
  `wmake` rewrites them, so mtime is already a reliable signal, and this keeps
  the digest cheap.
- Time directories (`0/`, `0.001/`, …), `postProcessing/`, `workflow_logs/`,
  and `processor*/` are excluded — they are outputs, not inputs, and including
  them would make every digest self-invalidating.

**Storage.** `WorkflowRunState` gains one optional field, `input_digest: str |
None = None`, serialized into `workflow_state.json`. `SweepManifest` gains
`case_input_digest` per case entry alongside the existing `override_hash`.
Both default to `None` so an existing on-disk state loads unchanged.

**Enforcement.** On resume, `run --strict` / `step --strict` recompute the
digest and compare:

| Recorded | Computed | Behaviour |
|---|---|---|
| `None` (pre-Phase-1 state) | any | Proceed, emit warning-level `provenance_unrecorded` |
| matches | matches | Proceed as today |
| differs | differs | **Refuse**: error-level `stale_inputs` diagnostic, non-zero exit |

The `stale_inputs` diagnostic names which component classes changed
(`system/`, `constant/`, or `executables`) so the agent can tell a dictionary
edit from a recompile without diffing by hand.

**Refuse, never auto-delete.** A mismatch does not wipe the output directory.
Deletion remains an explicit `--fresh`, consistent with the existing safety
posture. `--fresh` already deletes the directory, so it trivially satisfies
the check.

`sweep-run` applies the same rule per case. A case recorded as `completed`
whose digest no longer matches is not skipped; its manifest status becomes a
new terminal value `stale`, and the sweep refuses to proceed — the same posture
as the existing spec-hash mismatch, which this parallels at case granularity.
`stale` is not cleared by `--retry-failed` (that flag means "rerun failures",
not "accept changed inputs"); it is cleared by `--fresh` or a new
`--output-dir`.

**Testing.** Digest stability across repeated computation on an unchanged
case; digest change on a dictionary edit, on a large-file mtime change, and on
a touched executable; resume refused on mismatch; resume permitted on a
`None` digest with the warning present; `--fresh` clearing the refusal.

---

## Phase 2 — Cost model and sweep parallelism

Two independent additions that share a phase because both are mechanical.

### 2a. Deterministic cost estimate

`StrictPlanReport` gains a `cost_estimate` block, computed without executing
anything:

- `n_cells` — from the `polyMesh/owner` header `note` field when a mesh is
  present, otherwise the product of the `blockMeshDict` block counts.
- `n_timesteps` — `ceil(endTime / deltaT)` from `system/controlDict`.
- `n_writes` — derived from `writeControl` / `writeInterval`.
- `cell_steps` — `n_cells × n_timesteps`, the single comparable magnitude.

Each field is `None` when its input cannot be resolved; a partial estimate is
reported rather than suppressed, and an unresolvable estimate is never an
error on its own.

`cell_steps` above a warn threshold emits a warning-level `high_cost_estimate`
diagnostic. A new `--max-cell-steps N` promotes that to an error-level refusal.
Default is warn-only, so no existing plan starts failing. The flag applies to
`plan`, `run`, `sweep-plan`, and `sweep-run`; in sweep actions it is evaluated
per case at plan time, and a single case over the cap refuses the whole sweep
before anything executes, consistent with how `--max-cases` already gates.
An unresolvable estimate never triggers the refusal — only a computed value
over the cap does.

This is not a runtime predictor — it deliberately does not model per-cell cost,
ionic model complexity, or hardware. It exists so an agent can compare two
candidate configurations, and so a 10⁹-step accident is visible before launch
rather than after.

### 2b. `sweep-run --jobs N`

The loop in `sweep_runner.sweep_run` becomes a bounded thread pool of size `N`
(default 1, preserving today's behaviour exactly). Threads rather than
processes: each case is already executed as a subprocess, so the pool only
needs to wait on them.

Constraints:

- **Refused when `base.entry` is set.** Entry-mode materialization mutates the
  tutorial's shared `case_root`; `AGENT_GUIDE.md` already documents that these
  sweeps must never be parallelized. `--jobs > 1` with `base.entry` is an
  error before any case runs, not a warning.
- **Manifest writes are serialized** behind a single lock; the manifest stays
  atomically rewritten as it is today, so an interrupted parallel sweep resumes
  exactly like an interrupted sequential one.
- Generic-mode cases each materialize into their own `<output_dir>/<case_id>/`
  and share no state, which is what makes this safe.
- `--jobs` does not interact with MPI decomposition; a parallel-decomposed case
  still consumes `numberOfSubdomains` ranks, and oversubscription is the
  caller's responsibility. Documented, not enforced.

**Testing.** Cost estimate correctness against a known mesh and controlDict;
partial estimate when inputs are missing; warn vs refuse at the threshold;
`--jobs > 1` refused for entry-mode; a parallel sweep producing the same
manifest as the sequential one; resume-after-interrupt of a parallel sweep.

---

## Phase 3 — Run telemetry and result contract

### 3a. Telemetry

**New module:** `core/runtime/telemetry.py`, a pure parser over OpenFOAM's
`solverPerformance` log output. The format is stock OpenFOAM, not
cardiacFoam-specific — verified against a real log:

```
Time = 5e-05
diagonal:  Solving for Vm, Initial residual = 0, Final residual = 0, No Iterations 0
ExecutionTime = 0.04 s  ClockTime = 0 s
```

The parser is therefore solver-agnostic and takes no plugin dependency.

After a step completes, its stdout log is parsed into
`workflow_logs/<step>.attempt<N>.telemetry.jsonl`, one JSON record per
`Time =` block:

```json
{"time": 5e-05,
 "execution_time_s": 0.04,
 "clock_time_s": 0.0,
 "solves": [{"field": "Vm", "solver": "diagonal",
             "initial_residual": 0.0, "final_residual": 0.0, "iterations": 0}]}
```

Unparseable lines are skipped, never fatal. A log with no recognizable time
blocks yields an empty telemetry file and a `telemetry_unavailable`
warning — a utility step such as `blockMesh` legitimately produces none.

A summary is attached to `WorkflowStepState` as `telemetry_summary`:
`n_timesteps`, `last_time`, `max_final_residual`, `diverged`, `stalled`.

### 3b. Evidence-driven failure classification

Three diagnostic codes, derived from the summary:

- `solver_diverged` — any residual or reported field value is NaN or infinite.
- `solve_incomplete` — `last_time` is short of the case's `endTime` beyond a
  relative tolerance, while the process exited 0.
- `solver_stalled` — final residual fails to decrease across a configurable
  trailing window of timesteps while iteration count is at its ceiling.

All three **fail the step**, consistent with the status-driven-not-exit-code-driven
contract that `missing_artifacts` already establishes: a NaN is not a success.
All three are classified *fatal*, never retryable — rerunning an identical
diverged configuration cannot succeed.

Thresholds (stall window length, `endTime` relative tolerance) are per-step
`retry_policy`-adjacent configuration with documented defaults, so a case with
a legitimately flat residual can opt out without disabling the whole check.

`candidate_remediations` gains evidence-keyed hints for these codes — still
suggestions, still agent-applied, but now grounded in the telemetry rather
than in a static code→hint table.

### 3c. Expected observables

RunDocument v2 gains an optional `expectedObservables` array. Each entry:

```json
{"id": "vm_finite",
 "source": {"kind": "telemetry"},
 "assertion": {"type": "finite"}}
```

`source.kind` is `telemetry` or `artifact` (with `artifact_id`). The assertion
vocabulary is closed:

| `type` | Parameters | Meaning |
|---|---|---|
| `finite` | — | No NaN/inf in the referenced values |
| `range` | `min`, `max` (either optional) | All values within bounds |
| `monotonic` | `direction` | Values non-increasing / non-decreasing |
| `reached_end_time` | `rtol` | Final time within `rtol` of `endTime` |
| `convergence_order` | `expected`, `tolerance` | Observed order from a sweep's error series within tolerance |

Anything outside this set is rejected at schema-validation time, before
execution. Adding a type is a deliberate contract change.

Observables are evaluated at terminal status and written to
`observables_realized.json` beside `artifacts_realized.json`, each entry
carrying `id`, `status` (`passed` / `failed` / `not_evaluated`), the observed
value, and the assertion. A failed *required* observable fails the run; an
observable marked `optional` reports and does not.

`convergence_order` is the one assertion evaluated at sweep scope rather than
run scope, since it needs an error series across cases. It is evaluated by
`sweep-run` against the completed manifest.

**Testing.** Parser fixtures for a clean run, a diverged run, a truncated run,
and a non-solver utility log; each diagnostic code raised on its fixture and
absent on the clean one; all three classified fatal; every assertion type
passing and failing; unknown assertion type rejected at validation; a run with
no `expectedObservables` behaving exactly as today.

---

## Phase 4 — Complete the agnosticity seam

**Protocol additions.** `SolverPlugin` gains the members that the leaking
modules currently hardcode:

- `get_samplable_fields(resolved) -> Mapping[str, tuple[str, ...]]` — region
  name to field names. Replaces `_ELECTRO_SOLVER_FIELDS` / `_SOLID_SOLVER_FIELDS`
  and the `has_solid_region` inference in `capability_manifest.py`.
- `resolve_case_models(case_root)` moves behind the plugin: the literal
  `constant/electroProperties` path and the three `detect_*` calls are
  cardiacFoam's knowledge, not the core's.
- `get_required_case_files() -> tuple[str, ...]` — replaces
  `tutorial_contracts.REQUIRED_FILES`.
- `get_utility_roots() -> tuple[Path, ...]` — replaces `utility_catalog.py`'s
  fixed repo-relative scan root, so a second plugin can ship its own utilities.
- `get_override_schema()` — the `$ELECTRO_MODEL_COEFFS` schema text currently
  embedded in `introspection.py`.

The generic plugin implements every one of these with empty or minimal
returns, which is what makes `--plugin none` a real mode rather than a
degraded one.

**API version enforcement.** Core declares `SUPPORTED_PLUGIN_API_VERSIONS`.
`driver_context()` rejects a plugin outside that range with a clear error
naming both versions. Today `plugin_api_version` is validated as a non-empty
string and otherwise ignored.

**Discovery.** Plugins are discovered through an `importlib.metadata` entry
point group (`driverfoam.plugins`). `--plugin module:Class` is retained and
keeps its existing label as an unsafe local-development path. `--plugin` may
name a discovered plugin by id. Ambiguity between a discovered id and an
import target is resolved in favour of the discovered plugin, and the
`module:Class` form is recognized only by its colon.

**Compatibility.** The existing top-level shims (`dict_entries.py`,
`sweep_routing.py`, `sweep_materialize.py`) keep re-exporting from
`plugins.cardiacfoam.*`, and `COMPATIBILITY.md` is extended with each newly
delegated surface. No public import path breaks.

**Acceptance test — run it, do not read it.** The phase is complete when a
real case executes end-to-end under `--plugin none`: `plan --strict` produces
a plan with an empty-but-valid capability manifest, `run --strict` executes
its workflow DAG, telemetry is parsed, and the input digest is recorded. A
test that only asserts the imports resolve does not demonstrate agnosticity;
one that runs a case does.

**Testing.** Every new protocol member exercised on both the cardiac and
generic plugins; a plugin with an unsupported `api_version` rejected; entry-point
discovery finding a fixture plugin; `--plugin module:Class` still working; the
end-to-end generic run above.

---

## Cross-cutting decisions

**Everything is additive and defaults to today's behaviour.** New state fields
default to `None`, new diagnostics for absent data are warnings, `--jobs`
defaults to 1, cost estimation defaults to warn-only, and a RunDocument
without `expectedObservables` behaves exactly as it does now. An agent written
against the current contract keeps working.

**The driver reports; the agent decides.** Telemetry, cost estimates, and
observable results are evidence. The driver's only new judgements are the
three narrow, mechanically-defined failure codes in 3b — it still does not
choose a remediation, tune a parameter, or interpret physics.

**Digests and observables are separate concerns.** Provenance answers "did
this actually run with these inputs"; observables answer "is the answer any
good". Conflating them would make a legitimate re-run look like a bad result.

## Risks

- **Stall detection false positives** (Phase 3b) could block a legitimate run
  whose residual is genuinely flat — the case in this repository being an
  already-converged MMS solve. Mitigated by configurable thresholds and by
  keeping `solver_stalled` distinct from `solver_diverged`, so the unambiguous
  signal is not weakened by the heuristic one.
- **Digest cost on large meshes** (Phase 1). Mitigated by the size threshold;
  worth measuring against the largest torso mesh in the tutorials before
  fixing the 1 MiB cut-off.
- **Parallel sweep resource contention** (Phase 2b). `--jobs N` multiplied by
  MPI ranks can oversubscribe a workstation. Documented; not enforced, because
  the driver cannot know what else the machine is doing.
- **Phase 4 touches every consumer of the leaking modules.** Mitigated by
  doing it last, behind commits for Phases 1–3, and by keeping the shims.
