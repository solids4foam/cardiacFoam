# driverFOAM security model

driverFOAM lets a semi-trusted agent author a `RunDocument` (config + workflow
DAG + launch paths) that the strict executor runs as subprocesses. This is the
threat model the code is written against. It assumes a local, single-tenant
host.

## Trust boundaries

- **Trusted:** the `openfoam_driver` package; binaries on `PATH` and under
  `$FOAM_APPBIN` / `$FOAM_USER_APPBIN`; the ambient process environment
  (including `PATH` and `$FOAM_*BIN`).
- **Agent-provided (validated at ingestion):** the `RunDocument`. Validation
  happens once, at ingestion — `RunDocument` load + `build_execution_inputs` +
  the CLI `run`/`step` path:
  - `config` via `validate_run`.
  - `workflowDag` normalized, then the command allowlist
    (`validate_workflow_commands`) — a known OpenFOAM/driver command, an
    `Allrun`-family case script, a registered utility, or an installed OpenFOAM
    app. Absolute-path and arbitrary `./script` commands are rejected.
  - `launch.caseRoot` must be an existing, runnable OpenFOAM case
    (`registry._case_is_runnable`); `caseRoot`/`outputDir` are resolved to
    canonical absolute paths; when `DRIVERFOAM_ALLOWED_RUNS_ROOT` is set, both
    must resolve under it.
  This is the only path untrusted document content reaches execution.
- **Case-authored (untrusted, unsandboxed by design):** the contents of
  `Allrun`-family scripts. Running a case runs its scripts.

## Output-location contract

OpenFOAM outputs are `caseRoot`-relative: serial/reconstructed time
directories `caseRoot/<time>/<field>`, decomposed `caseRoot/processor<N>/<time>/<field>`,
function-object outputs `caseRoot/postProcessing/…`, config in `constant/` and
`system/`. The artifact-existence gate resolves patterns under the canonical
`caseRoot` and accepts either the reconstructed or the decomposed location for
time-indexed artifacts. `outputDir` is the driver-bookkeeping base (manifests,
workflow state, logs, plots); it defaults under `caseRoot` but may be a separate
results directory — it is not forced under `caseRoot`.

## Mitigations

- Single command allowlist owner (`validate_workflow_commands`), enforced once
  at ingestion.
- No command shadowing: bare names resolve via `PATH` only; only `Allrun`-family
  resolve case-locally.
- Workflow `cwd` cannot escape `caseRoot`.
- `caseRoot` must be a runnable OpenFOAM case; `caseRoot`/`outputDir` resolved to
  canonical paths; opt-in `DRIVERFOAM_ALLOWED_RUNS_ROOT` containment.
- Steps run argv-style (no shell).
- Override / spec **values** are rejected at the `update_foam_entry` write
  path if they are directive- or entry-terminating-shaped. The command
  allowlist gates *what binary runs*, not the *content* of the dicts it
  reads, so this is enforced at the write path instead: `mutators._format_value`
  (tier 1 — the path almost every override takes) raises `ValueError` on any
  value containing `;`, a newline, or `#`, before the value is written. The
  foamlib tier (tier 2, the line-scanner fallback) has its own explicit
  `_reject_directive_shaped` guard mirroring the same rule — it does not rely
  on foamlib's incidental type-strictness, which is narrower (foamlib only
  objects to a string that would read back as a different type, so it lets
  `'#includeEtcFuncs'`, a bare `'#'`, and `'PCG#calc'` through unconverted).
  This closes the gap previously recorded here as documented-but-unenforced:
  a value carrying `;` can no longer append a second dictionary entry, and a
  `#codeStream` / `#calc` / coded-function-object value can no longer reach a
  dict file through `update_foam_entry` to be compiled and executed by the
  solver at run time. This guard covers `update_foam_entry` specifically —
  see the dict-regeneration gap noted below, which it does not cover.

## Explicitly NOT mitigated

- Arbitrary code inside an invoked `Allrun` (running a case is running its code).
- No rlimit / output-size bounds (local DoS).
- Trusts the ambient `PATH` and `$FOAM_*BIN`.
- Assumes a single-tenant host.
- `run_workflow_step` is a trusted low-level primitive: a Python caller that
  invokes it directly with an unvalidated `case_root` / `log_dir` / `state_path`
  / command bypasses path and command validation. Untrusted document content
  never reaches the runner except through validated ingestion.
- The dict-regeneration/synthesis path is not screened for injection. When an
  override targets a selector key (e.g. `myocardiumSolver`),
  `specs/apply_overrides.py` routes sibling `$TOKEN.`-scoped override values
  in the same call to `RegenerationScope.regenerate`, which serializes them
  via `specs/dict_builder._openfoam_value_token` /
  `_serialize_block` in `plugins/cardiacfoam/dict_builder.py`. That path
  performs no injection screening at all — a `#codeStream` payload routed
  through it reaches the written dict file unrejected. Only
  `update_foam_entry`'s direct write path (both tiers, above) is guarded.
  This is a pre-existing gap; `dict_builder.py` synthesizes dicts rather than
  mutating them and was never migrated to the guarded write path.
