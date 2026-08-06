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

## Explicitly NOT mitigated

- Arbitrary code inside an invoked `Allrun` (running a case is running its code).
- Override / spec **values** (not just commands) are written verbatim into case
  dictionaries. A value containing an OpenFOAM coded entry — `#codeStream`,
  `#calc`, or a coded function object — is compiled and executed by the solver
  at run time. The command allowlist gates *what binary runs*, not the *content*
  of the dicts it reads; dict values are arbitrary code at solve time by design.
  This is acceptable under the local/single-tenant trust model (the same actor
  authoring values could run the solver directly) but is NOT a sandbox against a
  malicious value channel.
- No rlimit / output-size bounds (local DoS).
- Trusts the ambient `PATH` and `$FOAM_*BIN`.
- Assumes a single-tenant host.
- `run_workflow_step` is a trusted low-level primitive: a Python caller that
  invokes it directly with an unvalidated `case_root` / `log_dir` / `state_path`
  / command bypasses path and command validation. Untrusted document content
  never reaches the runner except through validated ingestion.
