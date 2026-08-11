# driverFOAM: revised next phases

Date: 2026-08-10  
Status: Phase 1 IMPLEMENTED (565b25de..HEAD); Phases 2-6 not started  
Plan 1 baseline: `ebeaf868`

## Decision

The direction is sound, but the original four-phase sequence is not safe to
implement as written. Provenance is necessary for trustworthy resume, and the
remaining agnosticity work is necessary for driverFOAM's core goal. Telemetry
and observable validation are valuable managed-mode features. Cost estimation
and concurrent sweeps are useful, but are not prerequisites for running a
normal OpenFOAM case.

Ordinary terminal execution remains outside this contract. A user can continue
to run `Allrun`, an OpenFOAM solver, or any case-authored script directly.
driverFOAM is an optional managed path that plans, records, resumes, and checks
the same case. Generic Allrun execution with `--plugin none` already has a real
end-to-end test; generic cases must not be restricted to cardiac layouts.

## Corrections to the original proposal

1. `0/` is normally an input and must not be blanket-excluded from provenance.
   The selected start-time directory, case scripts, the normalized workflow,
   plugin identity, and case-local includes also affect execution.
2. A single digest recorded at the beginning self-invalidates when an earlier
   workflow step legitimately creates `constant/polyMesh` or changes another
   downstream input. Provenance must be checkpoint-based.
3. Storing only an aggregate digest cannot produce a diagnostic that says
   whether `system/`, `constant/`, or executables changed. Structured component
   digests must be persisted.
4. A legacy completed state with no digest cannot both be accepted unchanged
   and be protected from stale replay. Safe behavior requires refusal unless
   the user explicitly accepts an unprovenanced state.
5. `step --apply` is an authorized input change. It needs a provenance lineage
   and downstream-state invalidation path; treating every change as stale would
   break the existing remediation loop.
6. Entry-mode sweeps mutate one shared tutorial root. Their old completed cases
   cannot be validated later by hashing the root, which contains only the most
   recently materialized case.
7. OpenFOAM `runApplication` commonly redirects solver output to `log.<app>`.
   Parsing only driverFOAM's captured stdout misses generic Allrun solves.
8. `solve_incomplete` is meaningful only for a solve step with an explicit
   target and compatible `stopAt` policy. `solver_stalled` is heuristic and
   must not be fatal by default in its first release.
9. Artifact assertions need a variable/column selector and a format reader.
   An `artifact_id` alone does not define which numeric values to inspect.
10. Sweep convergence order is a sweep-level assertion, not a per-run
    RunDocument assertion.
11. A thread pool plus a manifest lock does not handle Ctrl-C or terminate
    child process groups. Parallel execution needs explicit process lifecycle
    management and a single manifest-owning coordinator.
12. The remaining plugin leak includes command authorization and utility
    artifact lookup, not just samplable fields and utility scan roots. The core
    allowlist still names `cardiacFoam`, and artifact planning still reads the
    global `UTILITY_CATALOG`.
13. Adding `expectedObservables` to RunDocument v2 is not wire-compatible with
    older readers because the v2 schema has `additionalProperties: false`.
    The result contract belongs in RunDocument v3, with v2 loading retained.
14. The Plan 1 capability bundle should be extended. Core consumers should not
    return to calling a larger monolithic public plugin object directly.

## Invariants for all phases

- Direct OpenFOAM and Allrun terminal workflows remain unchanged.
- Public commands keep their names; behavior changes are explicit in JSON,
  diagnostics, state schema versions, and release notes.
- Core owns execution safety, normalized DAGs, state transitions, generic
  OpenFOAM parsing, and atomic persistence.
- Plugins own solver vocabulary, case semantics, extra commands, utility
  manifests, telemetry hints, artifact readers, and additional provenance
  inputs.
- The driver records evidence. It never silently tunes physical parameters.
- Each phase lands as an independently revertible commit and passes the full
  characterization suite from Plan 1.

## Phase 1 — Close the agnosticity and plugin-version seams

Do this before adding instrumentation so provenance, telemetry, and
observables do not acquire new cardiac assumptions.

### Internal capabilities

Extend `PluginCapabilities` with focused capabilities:

- `CommandAuthorizationCapability`: plugin solver commands, plugin utility
  manifests, and utility roots. The core retains only solver-neutral OpenFOAM,
  MPI, and case-script commands.
- `CaseIntrospectionCapability`: resolved case models and samplable fields.
- `CaseFileContractCapability`: structured required/conditional case-file
  rules, rather than one cardiac tuple.
- `OverrideSchemaCapability`: plugin-authored configuration and mutation help.
- `RuntimeEvidenceCapability`: declared solver-step roles, telemetry source
  globs, extra runtime dependencies, and artifact value readers.

Move these remaining consumers behind the bundle:

- `capability_manifest.py` fixed electro/solid fields and
  `constant/electroProperties` parsing.
- `tutorial_contracts.py` cardiac required files.
- `introspection.py` electro/physics override prose and examples.
- `utility_catalog.py`, strict planning, and artifact prediction's global
  utility catalog.
- `workflow.py`'s cardiac solver and utility command names.

### Public plugin compatibility

- Keep `SolverPlugin` API v1 loadable through the Plan 1 adapter.
- Define a distinct `SolverPluginV2` contract for the new capabilities; do not
  silently add required methods to v1.
- Set `SUPPORTED_PLUGIN_API_VERSIONS = {"1", "2"}` and reject any other
  version before catalog or case code runs.
- Migrate the built-in cardiac and generic plugins to v2.
- Keep v1 cardiac-shaped fallbacks only in `core.compatibility`, covered by
  tests and documented for later removal.
- Discover installed plugins through `importlib.metadata` group
  `driverfoam.plugins`. A colon continues to mean trusted `module:Class` local
  import. A discovered ID wins only when the argument contains no colon.
- Record entry-point distribution name/version in plugin identity provenance.

### Exit gate

- The Plan 1 characterization fixture is unchanged except for intentionally
  versioned capability fields.
- A real temporary Allrun case plans and executes under `--plugin none`.
- Its plan contains no cardiac command, field, required-file, utility, override,
  or artifact semantics.
- Cardiac and generic plugins exercise every v2 capability; a v1 fixture still
  loads through compatibility; an unsupported version is rejected.

## Carried into Phase 2 from Phase 1

Known agnosticity residual, deliberately not closed in Phase 1: `TutorialSpec`
carries `electro_properties_relpath` / `physics_properties_relpath`, and core's
generic-case factory populates them even under `--plugin none`. It predates
Phase 1 and closing it renames fields that flow into `spec.metadata` ->
`resolved_entry` -> the RunDocument, moving the `entry_resolution` and
`run_document_v2` digests for all twelve tutorials for reasons unrelated to
versioned capability fields — which would have destroyed Phase 1's ability to
demonstrate its whole delta was `api_version` + `capability_digest`. It belongs
in its own commit with its own before/after evidence.

It is pinned by
`tests/core/test_generic_plan_has_no_cardiac_semantics.py::test_known_residual_tutorialspec_carries_cardiac_field_names`,
which asserts the leak as a *known* fact. Closing it will fail that test
loudly; delete the test as part of the fix.

Also carried forward, from the Phase 1 whole-branch review: `resume_guard.py`
style/unused-parameter items (moot once Phase 2 deletes the module), the
missing integration test for `cli.py`'s failed-status resume payload,
`specs/function_object_fields.py` hardcoding cardiac regions, and
`utility_catalog.UTILITY_CATALOG` now having zero production consumers while
still being eagerly built from a hardcoded path.

## Phase 2 — Checkpointed input provenance

### Model

Add `core/runtime/provenance.py` with frozen, JSON-safe types:

- `ProvenanceComponent(kind, path, method, digest, size, mtime_ns, strength)`
- `ProvenanceSnapshot(schema_version, aggregate_digest, components,
  workflow_digest, plugin_identity)`

`WorkflowRunState` gains optional `resume_provenance`. Each
`WorkflowStepState` gains optional `input_provenance_digest`. Old JSON still
loads because both default to `None`.

The resume snapshot represents the filesystem immediately after the last
persisted transition, not only the original case. Before executing a step,
record its input digest. After any terminal step transition, recompute and
atomically persist the next resume snapshot. This permits `blockMesh` and
other producer steps to create legitimate downstream inputs without causing
self-invalidating state.

### Inputs

Canonical inputs include:

- normalized workflow DAG and selected plugin identity;
- `system/**` and `constant/**`;
- the actual initial-condition directory selected by `startFrom`/`startTime`,
  including `0/` when appropriate;
- case-local workflow contracts and the exact Allrun-family scripts referenced
  by the DAG;
- resolved executables using the same cwd, PATH, environment, and command
  resolution used by the executor, including the MPI payload executable;
- plugin-declared extra dependency paths and runtime libraries.

Do not follow an escaping symlink silently. Record the link and its resolved
external target fingerprint, or emit an explicit weak/unavailable provenance
component.

Small files use full SHA-256. The large-file policy must be benchmarked against
the largest supported case before fixing the threshold. Metadata-only entries
are allowed for large meshes, but must carry `strength="metadata"`; the driver
must not describe such a snapshot as content-complete.

### Resume policy

| Existing state | Comparison | Action |
|---|---|---|
| no state file | n/a | compute and record before first execution |
| pending state, no completed steps, no provenance | unrecorded | adopt current snapshot and warn |
| legacy state with completed/failed/running work, no provenance | unrecorded | refuse unless `--accept-unprovenanced-state` |
| recorded snapshot | match | resume |
| recorded snapshot | mismatch | refuse with `stale_inputs` and component diff |

`--fresh` remains the ordinary safe way to start a new lineage. It never
deletes the case root unless that root is already the explicitly safety-gated
output directory.

For `step --apply`, first verify that the pre-mutation state matches its saved
checkpoint, apply the mutation transactionally, invalidate the selected step
and all transitive dependents, record a provenance event, then establish the
new snapshot. Unknown pre-existing edits still refuse.

### Sweeps

- Generic isolated cases persist the structured snapshot digest in each
  manifest entry and compare it before a completed case is skipped.
- A mismatch marks the case `stale` and blocks the sweep. `--retry-failed`
  does not accept it.
- Entry-mode completed-case skipping is `unverifiable` while cases share one
  mutable root. Refuse that resume unless explicitly accepted, or rerun the
  case. Do not claim per-case provenance until isolated entry materialization
  exists.
- Preserve the original manifest `created_at` on resume.

### Exit gate

Tests cover unchanged stability; config, initial field, Allrun, workflow,
binary, MPI payload, and plugin dependency changes; legitimate changes made by
an earlier workflow step; legacy-state policy; `step --apply`; generic sweep
skip/refusal; and entry-mode unverifiability.

## Phase 3 — Deterministic cost evidence

Add a pure `core/runtime/cost_estimate.py`. Report evidence and confidence,
not a runtime prediction.

`CostEstimate` contains:

- `n_cells`, its source, and confidence;
- `start_time`, `end_time`, fixed/adaptive timestep classification;
- exact, lower-bound, or unknown `n_timesteps`;
- exact or unknown `n_writes`, honoring `writeControl`, `writeInterval`, and
  `writeAtEnd` where resolvable;
- exact/bounded `cell_steps` and an overall quality label.

Cell count precedence:

1. a valid global `constant/polyMesh/owner` header note;
2. a decomposed-mesh estimate explicitly labelled as such;
3. the sum of cell products across all parseable `hex` blocks in
   `blockMeshDict`;
4. unknown.

Use `(endTime - startTime) / deltaT`, not `endTime / deltaT`. Adaptive time
stepping, `latestTime`, non-`endTime` stop policies, expressions, and unsupported
write controls produce bounds or `None`, never invented precision.

Add `cost_estimate` to `StrictPlanReport`, not RunDocument v2. A default high
cost threshold may add a dedicated warning block; it must not silently change
the existing readiness score. `--max-cell-steps` is an opt-in refusal and is
recomputed at execution rather than trusting an agent-authored document.

Sweep planning evaluates every materialized case and refuses the whole sweep
before launch when a computed exact/upper estimate exceeds the explicit cap.
Unknown estimates do not fail.

## Phase 4 — Telemetry as evidence first

### Capture and parsing

Add a pure OpenFOAM log parser plus a source collector:

- parse driver-captured stdout/stderr;
- for Allrun-family steps, parse only declared or newly-created/modified
  `log.*` files so stale logs are not reused;
- attach source path, application segment, time, execution/clock time, solve
  field, solver, residuals, and iteration count to every record;
- tolerate MPI prefixes, multi-region output, restarts, repeated time blocks,
  and partial final records;
- write bounded/streamed JSONL without loading a solver log into memory.

`WorkflowStepState.telemetry_summary` defaults to `None`. It reports parsed
source files, record count, last time, maximum finite residual, non-finite
evidence, and parser quality. `telemetry_unavailable` is emitted only when a
step is explicitly a solve step or requests telemetry; utilities do not receive
noise by default.

### Gating rollout

Land telemetry collection as evidence-only first. In a separate commit:

- `solver_diverged` is fatal only for mechanically unambiguous non-finite
  residual/time evidence or an explicit fatal marker.
- `solve_incomplete` is fatal only when the workflow declares a solve target
  and the case uses a compatible `stopAt` policy.
- `solver_stalled` starts as opt-in or warning-only. Promotion to a fatal
  default requires a repository-wide evaluation fixture demonstrating no
  false positives for converged MMS and zero-iteration solves.

Remediation hints remain advisory and include the telemetry evidence that
triggered them. No automatic parameter mutation is added.

## Phase 5 — Versioned result and observable contracts

Introduce RunDocument v3 rather than changing the closed v2 schema. Continue
to load/execute v2 and provide an explicit v2-to-v3 migration. Do not reinterpret
old verification contracts automatically.

Run-level observable sources are typed:

- telemetry source plus a named metric/field selector;
- artifact source with `artifact_id`, variable/column selector, and optional
  time/case reduction.

Core provides bounded readers for explicitly supported generic formats such as
JSON and delimited tables. Plugins provide readers for solver/domain formats.
Unknown format/selector combinations become `not_evaluated` with a reason,
not an implicit pass.

Run-level assertions are the closed set `finite`, `range`, `monotonic`, and
`reached_end_time`. Each realization records assertion, observed value/summary,
status, reason, and evidence source. Only an explicitly required failed
observable fails the run.

`convergence_order` belongs in a versioned sweep result contract because it
requires a cross-case error series. It defines the mesh/time abscissa,
artifact/column source, fit method, minimum number of points, expected order,
and tolerance. It is evaluated only after the required cases complete.

Write `observables_realized.json` atomically and expose it as a core artifact.
Test every assertion and parser on passing, failing, missing, malformed, NaN,
and very large inputs.

## Phase 6 — Safe generic sweep concurrency

Add `sweep-run --jobs N`, defaulting to 1. Initially allow `N > 1` only when
every case has an isolated case root and the active materialization capability
declares concurrent use safe. Continue to reject shared-root entry mode before
any mutation.

Use a coordinator/worker design:

- workers materialize and launch one isolated case;
- only the coordinator mutates and atomically writes the manifest;
- the manifest retains expansion order regardless of completion order;
- counters and timestamps are coordinator-owned;
- existing `created_at` is retained on resume;
- subprocesses run in controllable process groups;
- Ctrl-C stops scheduling, terminates active groups with a bounded escalation,
  records interrupted/cancelled cases, and leaves a resumable manifest.

Do not share a plugin materializer between workers unless its capability says
it is thread-safe; otherwise create isolated plugin/context instances or keep
materialization serialized. MPI oversubscription remains documented and may
emit an advisory estimate, but is not silently enforced.

Acceptance compares sequential and parallel manifests after removing volatile
timestamps, exercises mixed success/failure, interrupts a live parallel sweep,
checks for orphaned children, and resumes to the same final semantic state.

## Priority and necessity

- **Required for trustworthy managed execution:** Phase 1 and Phase 2.
- **Strongly recommended for agent–solver parity:** Phase 4 and Phase 5.
- **Useful safety/ergonomics, independently deferrable:** Phase 3 and Phase 6.
- **Not required for ordinary OpenFOAM users:** all phases. Direct terminal
  execution remains the simplest path when planning, provenance, resume, and
  machine-readable evidence are not needed.

## Deferred

- C++ `-dumpContract` or generated model manifests.
- Live log streaming and early process termination on divergence.
- A general predicate language.
- Automatic remediation or physical parameter tuning.
- Parallel execution for shared-root entry sweeps.
- Mandatory full-content hashing of multi-gigabyte meshes.
