# driverFOAM Agent Guide

This is the agent contract for planning, launching, polling, and inspecting
cardiacFoam runs through the Python `openfoam_driver` package. Read this once
before driving the orchestrator.

## What the agent can do

| Action | Function | Module |
|---|---|---|
| Discover tutorials, dict keys, ionic models, utilities | `describe_tutorial(...)` | `openfoam_driver.core.introspection` |
| Build a non-mutating strict launch contract | `strict_plan(...)` | `openfoam_driver.core.strict_planning` |
| Execute an agent-authored RunDocument | `driverFoam run/step --run-document <file>`; `build_execution_inputs(...)` | `openfoam_driver.core.runtime.run_document_exec` |
| Execute one strict workflow step | `run_workflow_step(...)` | `openfoam_driver.core.runtime.workflow_runner` |
| Read/write strict workflow state | `workflow_state_from_json(...)`, `WorkflowRunState.to_json()` | `openfoam_driver.core.runtime.workflow_state` |
| Validate RunDocument v3 or migrate v1/v2 explicitly | `RunDocument.from_json(...)`, `RunDocument.migrate_v1(...)`, `RunDocument.migrate_v2(...)` | `openfoam_driver.core.runtime.run_model` |
| Validate a configuration before launching | `validate_run(run, *, entries=None)` | `openfoam_driver.core.specs.validation` |
| Synthesize a fresh `electroProperties` / `physicsProperties` | `build_electro_properties(...)`, `build_physics_properties(...)` | `openfoam_driver.plugins.cardiacfoam.dict_builder` |
| Parse an existing `electroProperties` back to selectors + overrides | `parse_electro_properties(path)` | `openfoam_driver.plugins.cardiacfoam.dict_builder` |
| Build + launch a one-shot run (runs through the strict executor) | `build_and_launch(...)` | `openfoam_driver.plugins.cardiacfoam.dict_builder` |
| Locate predicted outputs | `strict_plan(...)`'s `expected_artifacts` field (also in `driverFoam plan --strict` JSON) | `openfoam_driver.core.strict_planning` |
| Verify outputs vs predictions | `artifact_reconciliation` in `run --strict`/`step --strict` JSON output | `openfoam_driver.core.runtime.reconciler` |
| List past runs | `list_runs(root)` | `openfoam_driver.core.runtime.run_discovery` |
| Plan/run a parameter sweep | `driverFoam sweep-plan/sweep-run --spec sweep.json --output-dir <dir>` | `openfoam_driver.core.runtime.sweep_runner` |

## Preferred strict agent loop

The solver-injection refactor does not change this public loop. A single
per-operation driver context now supplies focused solver capabilities
internally, while omitted contexts, RunDocument v3, legacy fallbacks, commands,
diagnostics, and artifacts retain their established behavior.

Use strict planning before launching. It is the only path that tells an agent
whether the run is machine-readable, validated, catalog-covered, artifact
predictable, and workflow-addressable before execution starts.

For the cardiacFoam plugin, configure
`applications/scripts/driverFoam/driverfoam-runtime.example.yaml` once per
host and expose it through `DRIVERFOAM_RUNTIME_CONFIG`. The plugin declares
the `lightweight` and `full` physics backends in its `plugin.yaml`; the local
file selects one backend, its OpenFOAM bashrc, the full-mode solids4foam root,
and the generated `cardiacFoam.build.json` manifest. The manifest is not a
build step you run yourself: `runtime_profile.py` generates or refreshes it
automatically, on the fly, whenever it is missing or older than the compiled
`cardiacFoam` solver — by inspecting the solver's actual linked libraries
(`otool -L`/`ldd`) to infer which backend was compiled, never by trusting an
asserted flag. driverFOAM rejects an unset, invalid, unbuilt, or
compiled-metadata-mismatched selection instead of letting a shell resolver
silently select another checkout. This runtime file is separate from
case/sweep overrides and applies to all cardiacFoam entries.

```bash
driverFoam plan --strict --entry singleCell
driverFoam run --strict --entry singleCell
```

The `plan --strict` command is non-mutating. It prints JSON with:

- `status`: `ok` or `failed`
- `entry`: the raw entry identifier as requested (pre-resolution)
- `resolved_entry`: case/spec identity and paths
- `readiness_score`: weighted 0-100 score summarising whether the driver has
  enough concrete case-generation and run-preparation evidence to execute
- `simulation_audit`: scored stages showing exactly how simulations are created
  and prepared: `build_cases()`, required OpenFOAM files, dictionary
  resolution, workflow DAG normalization, artifact prediction, environment
  preflight, and mesh geometry
- `validation_diagnostics`: RunDocument and configuration validation results
- `workflow_diagnostics`: normalized workflow-DAG validation results (command
  allowlist, DAG structure)
- `catalog_coverage_errors`: strict dict-key coverage failures
- `artifact_diagnostics`: solver/utility/artifact prediction coverage failures
- `environment_diagnostics`: missing executables, unsourced OpenFOAM env, missing MPI launcher
- `mesh_geometry_diagnostics`: mesh-scale / geometry sanity checks
- `workflow_dag`: normalized executable steps
- `workflow_state`: initial pending step state
- `expected_artifacts`: predicted machine-readable artifacts
- `launch`: exact launch command and output paths

The `run --strict` command executes normalized steps until completion or
failure. It writes:

- `workflow_state.json` under the strict-plan output directory
- `workflow_logs/<step>.attempt<N>.stdout.log`
- `workflow_logs/<step>.attempt<N>.stderr.log`

If `workflow_state.json` already exists, `run --strict` resumes from that
state. If the saved state is `failed`, it exits non-zero and does not retry the
failed step automatically. Use `step --strict` for an explicit manual rerun:

```bash
driverFoam step --strict --entry singleCell --step solve
```

**Resuming can silently replay stale results.** If `workflow_state.json`
already says `completed` — e.g. a leftover case directory from a previous
session, code change, or experiment — `run --strict`/`step --strict` report
success and exit 0 without invoking the solver at all; there is no warning.
This was hit in practice: a sweep re-run after a solver code change reported
the previous day's numbers as fresh, caught only because the "new" errors
matched the old ones to six significant figures — two different code
versions cannot agree that precisely, so identical numbers meant identical
(non-)execution, not agreement. Any re-run intended as a genuine before/after
comparison after a code or config change MUST pass `--fresh`, which deletes
the resolved output directory before running so the workflow executes
exactly as it would on a first run:

```bash
driverFoam run --strict --entry singleCell --fresh
```

`--fresh` refuses to delete anything that doesn't look like driverFOAM's own
output (no `workflow_state.json`/`sweep_manifest.json`/`run_document.json`
found), the filesystem root, your home directory, or a path outside
`DRIVERFOAM_ALLOWED_RUNS_ROOT` when that's set — but it does not prompt for
confirmation, so treat any `--output-dir`/case directory you point it at as
fully disposable and copy out anything you want to keep first.

`--max-total-attempts <N>` caps the total number of step executions across the
whole run (a retry-storm guard on top of each step's per-step `max_attempts`).
It defaults to unbounded, preserving prior behavior.

For `sweep-run`, `--case-timeout-s <seconds>` sets a wall-clock timeout per case
subprocess; a case that exceeds it is recorded as failed (with a `timeout_error`
in its summary) and the sweep continues to the next case rather than hanging.
Defaults to no timeout.

Programmatic planning uses the same contract:

```python
from openfoam_driver.core.strict_planning import strict_plan

report = strict_plan("singleCell")
payload = report.to_json()
if payload["status"] != "ok":
    raise RuntimeError(payload)
print(payload["workflow_state"]["current_step_id"])
```

### Executing an agent-authored RunDocument

`plan --strict` emits a complete `run_document` (RunDocument v3) in its JSON
output. An agent can persist that document, edit it (e.g. tune `config`, add or
reorder `workflowDag` steps, set per-step `retry_policy`), and execute the
edited document directly — the driver runs *your* document instead of
regenerating one from `--entry`:

```bash
# 1. Plan and capture the run document the planner produced.
driverFoam plan --strict --entry singleCell > plan.json
python3 -c "import json; json.dump(json.load(open('plan.json'))['run_document'], open('run.json','w'))"

# 2. (optional) edit run.json — config, workflowDag, retry_policy, expectedArtifacts.

# 3. Execute the document. No --entry; --strict is implied by the document.
driverFoam run  --run-document run.json
driverFoam step --run-document run.json --step solve   # single step
```

`--run-document` is mutually exclusive with `--entry` (and with
`--config`/`--entry-kind`/`--tutorials-root`). Before executing, the driver:

1. Loads and schema-validates the document (a `version: "1"` document is
   migrated to v2 automatically).
2. Runs `validate_run` on its `config`.
3. Re-normalizes the supplied `workflowDag` and enforces the **command
   allowlist**: each step's command must be a known OpenFOAM/driver core
   command, a recognized case script (`Allrun`-family), a `UTILITY_CATALOG`
   entry, or an executable installed under `$FOAM_APPBIN`/`$FOAM_USER_APPBIN`
   (any core OpenFOAM app or your own compiled utility). Arbitrary non-OpenFOAM
   commands are rejected before anything runs. Note: when OpenFOAM is not
   sourced, only the core set + case scripts + `UTILITY_CATALOG` are accepted.
4. Requires `launch.caseRoot` and `launch.outputDir`.

If any of these produce an error-level diagnostic, the command prints
`{"status": "failed", "diagnostics": [...]}` and exits non-zero **without
executing anything**. Otherwise execution, `workflow_state.json` resume,
retry/backoff, and `failure_context` behave exactly as for the `--entry` path.

**Command-boundary guarantees.** Steps run argv-style (no shell). A step's
working directory cannot escape `caseRoot`. Bare command names resolve via
`PATH` only — a case directory **cannot shadow** a trusted binary such as
`cardiacFoam`. Only recognized case scripts (`Allrun`-family), named bare
(`Allrun`) or as `./Allrun`, resolve to case-local files; arbitrary
`./script` and absolute-path commands are rejected by the allowlist. Note
this does **not** sandbox the code *inside* an invoked `Allrun` — running a
case means running its scripts, which is arbitrary case-authored code by
design. The trust model is local/single-tenant: it assumes `PATH` and the
`$FOAM_*BIN` variables are not attacker-controlled.

See [`SECURITY.md`](SECURITY.md) for the full trust model, output-location
contract, and the explicit list of what is and is not mitigated. For the
plugin-boundary compatibility fallbacks (optional-hook defaults, legacy
shims), see `openfoam_driver/core/compatibility.py`.

## Compatibility one-shot loop

The legacy `build_and_launch(...)` path remains supported for existing scripts,
but it is not the preferred autonomous path because it mutates and launches in
one call instead of first emitting a strict contract.

```python
from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch

result = build_and_launch(
    electro_selectors={
        "myocardiumSolver": "singleCellSolver",
        "ionicModel": "AlievPanfilov",
        "tissue": "myocyte",
    },
    physics_selectors={"type": "electroModel"},
    case_dir="/tmp/my_run/case",
)
print(result)  # {"case_dir": ..., "status": "complete", "workflow_state": {...}}
```

That single call:

1. Calls `build_electro_properties(...)` and `build_physics_properties(...)`.
2. Runs the validator on both — raises `ValueError` if your selectors break a structured constraint.
3. Writes `case/constant/electroProperties` and `case/constant/physicsProperties`.
4. Constructs a `generic_case` spec pointing at the case directory, whose
   `workflow_dag` mirrors `pre_solve_commands`/`solver_command` exactly.
5. Runs that workflow_dag to completion through the same strict executor
   `run --strict` uses (`run_workflow`), validating every command against
   the allowlist first.
6. Returns the final `workflow_state` (per-step status, logs, exit codes).

## Sweeping a parameter grid

For running many cases off one parameter grid, use `sweep-plan`/`sweep-run`
instead of hand-looping `build_and_launch`. A `sweep.json` has two top-level
objects:

- `"base"`: fixed values applied to every case — `electro_selectors`,
  `physics_selectors`, `electro_overrides`, `physics_overrides`, `delta_t`,
  `end_time` (same shapes as `build_and_launch`'s kwargs).
- `"sweep"`: `"mode"` (`"cross_product"` or `"zip"`), `"independent"` (axis
  name to list of values, routed into `build_and_launch`'s parameters per
  below), and `"dependent"` (a list of `{"name", "derive", "of"}` entries for
  derived *labels only* — not routed through the selector rules below —
  currently the only registered `derive` function is `case_id_template`,
  which joins the named `of` values into a `caseId` label).

```json
{
  "base": {
    "electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells"},
    "physics_selectors": {"type": "electroModel"}
  },
  "sweep": {
    "mode": "cross_product",
    "independent": {"ionicModel": ["TNNP", "BuenoOrovio"], "deltaT": [1e-6, 2e-6]},
    "dependent": [{"name": "caseId", "derive": "case_id_template", "of": ["ionicModel", "deltaT"]}]
  }
}
```

Each resolved case's axis values route automatically into `build_and_launch`'s
parameters: `myocardiumSolver`/`ionicModel`/`tissue` go to `electro_selectors`,
`type` goes to `physics_selectors`, `deltaT`/`endTime` go to the dedicated
`delta_t`/`end_time` kwargs, `dx` goes to the dedicated `dx` kwarg (mesh
resolution in mm, see below), any other `system/controlDict` key is rejected
outright, and any key that isn't a recognized electroProperties/
physicsProperties driver_path is rejected outright too (it would otherwise
have no effect on the generated case). Everything recognized falls through
to `electro_overrides`.

Every case is *materialized* fresh: `build_and_launch(..., dry_run=True)`
writes its dict files, and the sweep runner additionally writes a generated
`Allrun` script and a `workflow_contract.json`, into `<output_dir>/<case_id>/`.
This is not a registered-tutorial lookup; each case is its own on-disk
`case_folder` entry.

### Mesh provisioning for from-scratch cases

A freshly materialized `case_folder` has no author-supplied mesh, so
`build_and_launch` provisions one based on `myocardiumSolver`:

- `singleCellSolver` (no real geometry): a bundled static 1-cell polyMesh is
  copied into `constant/polyMesh/` directly — no `blockMesh` step needed.
- `monodomainSolver`/`bidomainSolver`/`eikonalSolver` (need real geometry): a
  generic default `system/blockMeshDict` is written (a small cubic slab,
  "walls" patch — **not** tuned to any specific tutorial's science), and the
  generated `Allrun` runs `blockMesh` before `cardiacFoam`. Sweep this mesh's
  resolution with the `dx` axis (**metres**, isotropic cell size — note this
  differs from `niederer_2012.py`'s own `DX_VALUES`, which are in
  millimetres; the two are unrelated mechanisms, see below). `dx` derives
  the cell count for the fixed default slab size via
  `specs/mesh_provisioning.py::cell_counts_from_dx`, which raises
  `ValueError` if `dx` does not evenly divide the slab size — deliberately
  no silent rounding, matching the same rigor
  `niederer_2012.py::_replace_blockmesh_resolution` already established for
  its own (different, millimetre, non-cubic) slab; both now share the
  `cell_counts_from_dx` calculation, differing only in how the result gets
  written (`mesh_provisioning.py` generates a fresh file from its own
  template; `niederer_2012.py` patches an existing author-provided file).
  `dx` is meaningless for `singleCellSolver` (no geometry to resolve) and
  raises `ValueError` rather than silently having no effect. `dx` also has
  nothing to do with real anatomical meshes imported via
  `vtkUnstructuredToFoam` (most real tutorials) — those are unstructured
  meshes with no cell-size concept, and this mechanism never touches them.
- A mesh already present under `constant/polyMesh/` or `system/blockMeshDict`
  is never clobbered by a repeat `build_and_launch` call, regardless of that
  call's own `overwrite` flag — this protects a hand-authored custom mesh
  from being silently replaced by the generic default.
If the sweep declares a `caseId` dependent entry, it becomes the case's
directory name (validated for uniqueness and path-safety); otherwise cases
are named `case_0001`, `case_0002`, ... in expansion order.

Both actions enforce a safety cap of 200 expanded cases by default (override
with `--max-cases`), checked before any case is expanded or materialized:

```bash
driverFoam sweep-plan --spec sweep.json --output-dir .tmp/driverfoam/sweeps/my_sweep/
driverFoam sweep-run --spec sweep.json --output-dir .tmp/driverfoam/sweeps/my_sweep/
```

`sweep-plan` materializes and strict-plans every case without launching
anything. `sweep-run` additionally launches each case and is resumable:
re-invoking it against the same `--output-dir` skips cases already recorded
as `completed` in `sweep_manifest.json`, leaves `failed` cases alone unless
`--retry-failed` is passed, and refuses to proceed at all if `sweep.json` has
changed since that output directory's manifest was created (a spec-hash
mismatch) — use a fresh `--output-dir` or resolve the mismatch first.

`--fresh` applies here too, and matters more: a solver/code change
invalidates every case in the sweep equally, so `sweep-run --fresh` deletes
the *entire* `--output-dir` (not just individual cases) before re-running
everything from scratch — this also sidesteps the spec-hash-mismatch refusal
above, since there's no old manifest left to compare against. Mutually
exclusive with `--retry-failed` (resume-only-failures vs. wipe-everything are
contradictory intents).

See `openfoam_driver/core/runtime/sweep_runner.py` for the full implementation.

### Sweeping an existing registered tutorial (`base.entry`)

The generic mode above always materializes a fresh, from-scratch `case_folder`
via `build_and_launch`. Some tutorials (`niederer2012`, `manufacturedMonodomainPseudoECG`, and
others under `openfoam_driver/specs/tutorials/`) instead expose their own
`make_spec(**kwargs)` with tutorial-specific parameters (e.g. `niederer2012`'s
`dx_values`/`dt_values`/`end_time_by_dx`, in millimetres/milliseconds;
`manufacturedMonodomainPseudoECG`'s `dimensions`/`number_cells`/`dt_values`). To sweep one of
these instead of a from-scratch case, set `base.entry` to the tutorial's
registered name:

```json
{
  "base": {
    "entry": "niederer2012",
    "solvers": ["implicit"],
    "end_time_by_dx": {"0.5": 0.2, "0.2": 0.08, "0.1": 0.055}
  },
  "sweep": {
    "mode": "zip",
    "independent": {
      "dx_values": [[0.5], [0.2], [0.2], [0.2], [0.1], [0.1], [0.1]],
      "dt_values": [[0.01], [0.01], [0.005], [0.001], [0.01], [0.005], [0.001]]
    },
    "dependent": [
      {"name": "output_dir_name", "derive": "output_dir_name_template", "of": ["dx_values", "dt_values"]}
    ]
  }
}
```

Every axis value is forwarded verbatim as a keyword argument to that
tutorial's own `make_spec(**overrides)` — there is no fixed vocabulary the way
generic mode has (`electro_selectors`/`dx`/etc.); `make_spec` validates its
own keyword arguments and an unrecognized one is a normal `TypeError`,
reported as that case's `materialization_error`, same as any other per-case
failure. Values fixed across every case in the sweep (like `solvers`/
`end_time_by_dx` above) go in `base`; per-case values come from
`independent`/`dependent` and win on conflict.

**One case per resolved combination, and why.** Several of these tutorials'
own `apply_case()` methods patch `system/controlDict`/`system/blockMeshDict*`
directly instead of writing an isolated per-case directory. driverFOAM stages
a fresh copy under the disposable workspace before applying those mutations,
but each resolved axis combination must still collapse to exactly one case —
if it doesn't (e.g. a config that still fans out internally because a
constraining kwarg like `solvers` is missing),
`sweep-plan`/`sweep-run` reports that case as `failed` with a clear
`materialization_error` rather than silently applying only the first of
several. In practice this means giving `dt`/`dx`-style axes their own
dedicated sweep row (`"zip"` mode with per-case single-element lists, as
above) instead of relying on the tutorial's own internal multi-value fan-out.
Entry-mode sweeps remain serial because each case owns its staged case tree
and post-processing boundary.

**Entry-mode execution is staged and disposable.** `sweep-plan` and
`sweep-run` copy the registered tutorial into
`<output_dir>/cases/<case_id>/` before calling `apply_case()`; the source under
`tutorials/` is never the mutable execution root. If `--output-dir` is omitted,
the default is `<repo>/.tmp/driverfoam/sweeps/<spec-name>`. All generated
meshes, processor/time directories, logs, workflow state, manifests,
post-processing output, and archives must stay below this repository-local
`.tmp/driverfoam/` workspace. Keep a failed workspace when diagnosing a run;
cleanup is an explicit, disposable-output action.

Everything else — the manifest, `--retry-failed`, `--case-timeout-s`,
`--max-cases`, resumability — is identical to generic mode.

## Polling a long-running run

For a run that takes minutes, prefer the async-friendly polling pattern,
against the real run-state file:

```python
import json
import time
from pathlib import Path

state_path = Path("<output_dir>/workflow_state.json")
while True:
    state = json.loads(state_path.read_text())
    if state["status"] in {"completed", "failed", "skipped"}:
        break
    time.sleep(15)
```

`workflow_state.json` is written by the strict workflow orchestrator and
updated after every step, so the read above is safe at any instant.

## Post-processing phase (brain + module)

The execution engine hands off to the postprocessing phase once a workflow or sweep reaches a terminal state. This is split into two independent pieces:

1. **The brain (`build_sweep_context`)**: Reads the sweep's own record (`sweep_manifest.json`), verifies it against what is actually on disk (resolving entry-mode vs generic-mode output directory differences), and returns a single grounded `SweepContext`.
2. **The postprocessing module (`run_postprocessing_module`)**: A separate function that receives the `SweepContext` and a task. **It never re-reads the manifest or re-derives file locations.** It lists each case's postprocessing script catalog via `list_postprocess_scripts()`.

If an agent needs deeper reasoning than the flat summary, it must use the brain's query functions:

- `read_case_workflow_state(context, case_id)`
- `read_case_output_file(context, case_id, relative_path)`

These query functions raise clearly on an unknown case ID and safely restrict reads to files the brain has already verified.

### Authoring postprocessing scripts

Every postprocessing script in a tutorial's `setup/` directory must expose a `run_postprocessing` function matching the `PostprocessingProtocol` signature:

```python
def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **kwargs: object) -> list[dict]: ...
```

The script's docstring is statically extracted as its `description`, allowing reasoning agents to decide if the script applies to a task. Reusable plotting and styling utilities are exposed under `openfoam_driver.postprocessing`.

## Verifying outputs

Strict planning predicts artifacts before launch and assigns artifact ids to
workflow steps when catalog coverage is available. The strict step/run path
now reconciles claimed artifact ids against on-disk files after each step. If an
expected artifact is missing, the step automatically fails with a `missing_artifacts` code.

### Reading a failed strict step

When `run --strict` fails or a `step --strict` ends `failed`, the printed JSON
carries a top-level `failure_context` object for the failed step:

- `step_id`, `attempt`, `exit_code` — identity of the failed attempt. Note
  `exit_code` may be `0` even on failure (e.g. `missing_artifacts`): the
  contract is **status-driven**, never exit-code-driven.
- `diagnostics` — the diagnostic codes the runner emitted.
- `stdout_log` / `stderr_log` — paths, for a full read.
- `stdout_tail` / `stderr_tail` — the last `--tail-lines` lines (default 200) of
  each log, bounded to 64 KiB.
- `stdout_truncated` / `stderr_truncated` — whether content was dropped.

The driver surfaces raw tails and status only. It does **not** judge convergence
or pick a fix — interpretation and remediation are the agent's job. The loop is:
read `failure_context` → edit the case dict (e.g. via `build_electro_properties`
or `mutators.py`) → `step --strict --step <id>` reruns the failed step (the
`attempt` counter increments).

To shorten that loop, `failure_context` also carries a
`candidate_remediations` array — **suggestions only**, the agent applies them.
Each entry has `diagnostic_code`, `driver_path`, `change` (a human-readable
transform, descriptive), `rationale`, `source` (`"static"`), and `confidence`. A
hint with an empty `driver_path` is advisory. The ladder emits static,
diagnostic-code-keyed hints; when one matches, it points the agent straight at
the failure. When none matches, the array is empty and the agent reasons from
`failure_context` and the catalog. For numerical control such as `deltaT`, the
per-ODE stability limit is the anchor: around `1e-6` s for biophysical
(Hodgkin-Huxley-style) ionic models and around `2e-5` s for phenomenological
models.

To apply a chosen fix mechanically, write an overrides file
(`[{"driver_path": "...", "value": "..."}]`) and run:

```
driverFoam step --strict --step <id> --apply overrides.json
```

This validates each override for *applyability*, applies it via the dict mutators
(resolving `$ELECTRO_MODEL_COEFFS.*` to the case's solver-specific coeffs block),
reruns the step (`attempt++`), and appends one record to `remediation_history.jsonl`
under the output directory. The driver accepts three forms of overrides:

1. `$ELECTRO_MODEL_COEFFS.*`: Catalog-addressable entries. For a `dynamic_path`,
   replace each template placeholder with the concrete instance name in the
   `driver_path` (for example,
   `$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.global.scale.myChannel`). The
   concrete key must already exist in the generated dictionary.
2. `system/path/to/dict:entry_path`: Explicit overrides for any OpenFOAM dictionary (e.g., `system/fvSolution:solvers/V/tolerance`). The file path must be strictly inside `system/`. If the case uses multiple regions (e.g., electromechanics), check `constant/physicsProperties` to determine if you need to target `system/electro/fvSolution` or the top-level `system/fvSolution`.
   - **Note on entry paths**: `/` traverses nested blocks. OpenFOAM lets a sub-dictionary be keyed by a quoted regex instead of a literal name (e.g. a solver block declared as `"Vm|VmFinal|u|uFinal"`); mutators.py resolves an ordinary member name (`solvers/Vm/tolerance`) against such a pattern automatically, so you do **not** need to know the pattern or spell it out in quotes — just use the field name you actually mean (e.g. `system/electro/fvSolution:solvers/Vm/tolerance`). An exact literal key always wins over a pattern match if both exist.
3. Flat string paths (e.g., `deltaT`): Routed to `system/controlDict` for backward compatibility.

Invalid overrides are rejected **before** any mutation or rerun.

**Derived constants are not overridable.** Some models expose constants that are
*computed* from other (user-facing) constants at `initConsts` — e.g. the
Land-Niederer active-tension transition rates (`AC_k_uw`, `AC_k_ws`, `AC_k_wu`,
`AC_k_su`, `AC_cds`, `AC_cdw`, `AC_ktm_block`, `AC_A`, `AC_XSSS`, `AC_XWSS`,
`AC_fPKA_TnI`, `AC_PKAForceMultiplier`). The active-tension catalog deliberately
omits these from its `constants` list, and overriding one has no effect (the
solver recomputes it from its inputs). To *change* such a quantity, override the
user-facing constants it derives from. An agent may still **reason about** derived
values (e.g. predict how halving `AC_dr` shifts `AC_k_su`) — just don't try to set
them directly. (Note: the ionic catalog, which is auto-generated from the full C++
constant enum, *does* list derived constants; the same rule applies there — listed
≠ overridable.)

`run --strict` and `step --strict` already compare predicted artifacts to
on-disk reality after every step/run -- no separate file to read. The
printed JSON payload carries an `artifact_reconciliation` object, built by
`reconcile_artifacts()` (`core/runtime/reconciler.py`):

```python
payload = json.loads(run_strict_stdout)  # the JSON run --strict prints
reconciliation = payload["artifact_reconciliation"]
print(reconciliation["matched_count"], "/", reconciliation["predicted_count"])
for artifact in reconciliation["artifacts"]:
    if artifact["status"] == "missing" and not artifact["optional"]:
        print("  warning missing required:", artifact["artifact_id"])
```

Missing-but-optional artifacts are not errors. They only appear under specific
configurations, for example probes that were not enabled.

## Discovering what's valid

Three layers of discovery:

1. **What tutorials exist?** `from openfoam_driver.core.introspection import describe_launch_matrix; describe_launch_matrix()` returns every registered entry.
2. **What dict keys can I set?** Iterate `openfoam_driver.dict_entries.ELECTRO_PROPERTY_ENTRY_GROUPS` and `PHYSICS_PROPERTY_ENTRIES` for case-physics entries. For time-control use `openfoam_driver.dict_entries.CONTROL_DICT_ENTRIES` (`deltaT`, `endTime`). Each entry carries `driver_path`, `value_kind`, `enum_values`, `unit`, `typical_value`, and structured constraints (`applicable_when`, `forbidden_when`, `required_when`, `mutually_exclusive_with`).
3. **What ionic models can I pick?** `from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG`. Each entry carries `states`, `algebraic`, `compatible_solvers`, `compatible_tissues`, `species`, `cardiac_region`, `recommended_exports`.
4. **What utilities are known?** `from openfoam_driver.core.utility_catalog import UTILITY_CATALOG`. Strict planning fails when a workflow command has missing required `produces` metadata.
5. **What dict keys have parser limitations?** Read `openfoam_driver/plugins/cardiacfoam/dict_key_allowlist.json`. Strict dict-key scanning fails when new uncatalogued keys appear, stale catalog paths remain, or allowlist entries become unused.
6. **What commands may a workflow step run, and what fields may a function object sample?** Read the `capability_manifest` block emitted by both `describe --entry <name>` and `plan --strict --entry <name>` (and `describe_entry(...)` / `strict_plan(...).to_json()` programmatically). It is the authoritative, machine-readable accept-surface: `allowed_commands` (`core`, `case_scripts`, `utilities`, plus the `$FOAM_APPBIN` note) mirrors the command allowlist exactly, and `samplable_fields` lists the field names the *resolved* model exposes,
keyed by region. **Both blocks are plugin-dependent.** For cardiacFoam the
regions are `electro` / `solid`; under `--plugin none` neither key is
present (only `note`), so read the keys that are there rather than
assuming a fixed set. Author `workflowDag` commands and `functions{}` field lists against this instead of guessing — a command outside `allowed_commands` is rejected before execution, and a field outside `samplable_fields` is dropped silently by the solver (see below).

**Hand-built case directories need both an `Allrun` and a `workflow_contract.json`.**
A directory resolved as `entry_kind="case_folder"` (any case directory under
`tutorials_root` that isn't a registered tutorial) needs an executable
`Allrun` script *and* a `workflow_contract.json` whose `"steps"` array is
non-empty. Without a populated `"steps"` array, the registry silently sets
the resolved entry's workflow DAG to `None` — there is no diagnostic that
names `workflow_contract.json` or `Allrun` specifically, so `strict_plan`
just blocks at the `workflow_preparation` stage with a generic "workflow DAG
is missing or invalid" error and no pointer to the actual cause. This was
undocumented until it was hit directly while building the sweep feature
(`sweep_materialize.py` writes both files for exactly this reason).

## What the validator catches

`validate_run(run)` runs seven families of checks:

- **Required fields** — every `required` entry has a value.
- **Enum membership** — values for enum-typed entries are in `enum_values`.
- **Structured constraints** — `applicable_when` / `forbidden_when` / `required_when` / `mutually_exclusive_with`.
- **Solver coupling** — pairings like (`singleCellSolver`, any Purkinje) reject with the table's stated reason.
- **Block references** — `domainCouplings.<name>.conductionNetworkDomain` must point at a declared block.
- **Tissue heterogeneity** — `ionicHeterogeneity` requires a supported `ionicModel` and `endoMInterface < mEpiInterface`.
- **Tissue compatibility** — `tissue` must be in the `ionicModel`'s `compatible_tissues`.

If the dict builder rejects your input with `ValueError`, the message lists every violation. Fix the selectors or overrides and call again.

## Function objects (probes, sampling, sets, …)

Function objects are **OpenFOAM's, not driverFOAM's.** Anything you put in a
case's `controlDict` `functions { … }` block is defined by the OpenFOAM
documentation, not by this driver — so there is no driver catalog, builder, or
helper for them, and there shouldn't be. Author them the normal OpenFOAM way:

- **Reuse OpenFOAM's shipped library.** `functions { #includeFunc probes(...) }`
  pulls a ready-made, documented object from `$FOAM_ETC/caseDicts/postProcessing/`.
  `ls "$FOAM_ETC/caseDicts/postProcessing"` lists what is available — that
  directory *is* the reference; do not re-derive these from tutorials.
- **Or write a full typed block** (`type probes; libs (...); fields (...);
  probeLocations (...);`) exactly as the OpenFOAM docs specify. To attach it to
  an existing case, write the fragment into `system/<Name>` and `#include` it
  from a `functions{}` entry, or set it through the `system/<dict>:<entry>`
  override form (the `system/path/to/dict:entry_path` form documented above).

**The only parts you can't get from OpenFOAM docs — because they are
cardiacFoam-specific:**

- **Sample-able field names.** The *object* is OpenFOAM's; the *fields* it can
  sample are this solver's: membrane voltage `Vm`, `activationTime`, total
  ionic current `Iion`; active tension `Ta` and fibre stretch `lambda`;
  bidomain potentials `phiE` / `phiI`; per-ionic-model species (e.g. `Ca_i`).
  The authoritative, model-specific list is the catalogs already noted under
  "Discovering what's valid" (`IONIC_MODEL_CATALOG` states / algebraic /
  `recommended_exports`, `ACTIVE_TENSION_MODEL_CATALOG`). Sample only names that
  exist for your chosen model, or the solver drops them.
- **Regions (multi-region cases only).** Electromechanical cases split fields
  across two regions: `electro` (`Vm`, `Ca_i`, ionic state) and `solid` (`Ta`,
  `lambda`, mechanics). A function object on such a case must carry
  `region electro;` or `region solid;` accordingly. Single-region electro cases
  take no `region` entry. But see the electromechanics note below before
  driving such a case at all.

> ### Electromechanics is not currently working
>
> **Do not select an electromechanical entry, and do not try to fix one.**
> `manufacturedMonodomainTotalLagrangianEM` is registered but fails
> `plan --strict` with `missing_electro_properties`, `myocardiumSolver is
> required` and `empty_artifact_prediction`. The cause is known: these cases
> lay their dicts out per region (`constant/electro/electroProperties`,
> `constant/solid/solidProperties`) while the planner looks for
> `constant/electroProperties`.
>
> This is a deliberately deferred gap, not a defect to discover. If you are
> here because a plan failed on that entry, the correct response is to report
> that electromechanics is unsupported and stop -- not to repair the planner,
> move the dicts, or work around it. Every other registered entry plans clean;
> if one of *those* fails, that is a real regression worth investigating.

Outputs land where OpenFOAM puts them:
`postProcessing/<functionObjectName>/<time>/<field>`.

**Strict planning now checks sampled field names.** `plan --strict` parses each
`controlDict` `functions{}` sub-dict's `fields (...)` list and emits a
**warning-level** `unknown_sampled_field` diagnostic (in the report's
`function_object_diagnostics`) for any field the resolved model does not expose —
`region solid;` blocks are checked against the mechanics fields, everything else
against the electro fields (`capability_manifest.samplable_fields`). This is
**non-blocking**: it never fails a plan, because the catalog can lag the C++
solver and a false positive must not block a run — but it turns the solver's
otherwise-silent field drop into a visible signal. `#includeFunc` shorthands are
not parsed (their field lists live in `$FOAM_ETC/caseDicts`). Set
`SKIP_FUNCTION_OBJECT_DIAGNOSTICS=1` to bypass the check entirely.

## Common patterns

### Override a single dict key

```python
build_electro_properties(
    selectors={"myocardiumSolver": "monodomainSolver",
               "ionicModel": "TNNP",
               "tissue": "epicardialCells"},
    overrides={
        "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude": "60",
        "$ELECTRO_MODEL_COEFFS.solutionAlgorithm": "implicit",
    },
)
```

Override paths use the full `$ELECTRO_MODEL_COEFFS.<key>` form. Top-level keys (like `myocardiumSolver`) live in `selectors`, not `overrides`.

**Block-gated families.** Some groups only appear once you configure them.
`singleCellStimulus.*` is one: override any key under it — as above — and the
rest of the family fills from its typical values, so the four keys the solver
requires together (`stim_start`, `stim_period_S1`, `stim_duration`,
`stim_amplitude`) are never written half-complete. Override none of them and
**no stimulus block is generated at all**, which is deliberate: a run without a
stimulus is legal, and inventing one from defaults would silently pace a case
that asked for nothing. `bathPotentialDomain.*`, `ecgDomains.*` and
`conductionNetworkDomains.*` behave the same way.

### Configure a bath bidomain run

```python
build_electro_properties(
    selectors={"myocardiumSolver": "bidomainSolver",
               "ionicModel": "bathBidomainFDAManufactured"},
    overrides={
        "$ELECTRO_MODEL_COEFFS.bathPotentialDomain.bathCellZones": "(bath organ)",
    },
)
```

Declaring any `bathPotentialDomain.*` override auto-enables the bath block — the bath leaves typical-value default unless overridden.

### Read back an existing dict

```python
from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties

parsed = parse_electro_properties("/path/to/case/constant/electroProperties")
# {"selectors": {"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP", ...},
#  "overrides": {"$ELECTRO_MODEL_COEFFS.solutionAlgorithm": "explicit", ...}}
```

Pass the result directly to `build_electro_properties` to round-trip:

```python
from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties, parse_electro_properties

parsed = parse_electro_properties(existing_path)
text = build_electro_properties(parsed["selectors"], overrides=parsed["overrides"] or None)
```

Only non-default values appear in `overrides`. Entries matching the catalog's
`typical_value` are omitted. `dynamic_path` entries and keys outside the
catalog are silently ignored by the parser, but strict planning and the strict
dict-key scanner are the contract gates for new generated plans.

### Run a smoke test before a full sweep

```python
build_and_launch(
    electro_selectors={
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "TNNP",
        "tissue": "epicardialCells",
    },
    physics_selectors={"type": "electroModel"},
    case_dir="/path/to/case",
    end_time=0.001,   # 1 ms — just enough to verify the case launches
    delta_t=0.0001,
)
```

If the call returns without raising, the case structure, boundary conditions, and property files are consistent enough to run. Then widen `end_time` for production. Requires an existing `system/controlDict` in the case directory — `build_and_launch` patches it in-place.

### Run with pre-solve commands

```python
build_and_launch(
    electro_selectors={...},
    physics_selectors={"type": "electroModel"},
    case_dir="/tmp/my_run/case",
    pre_solve_commands=["blockMesh", "setTorsoOrganConductivityField"],
    openfoam_bashrc="/opt/openfoam/etc/bashrc",
)
```

Each entry in `pre_solve_commands` runs in `case_dir` before `cardiacFoam`. Strings are shell-split; lists are passed directly. When `openfoam_bashrc` is set every command is sourced into the OpenFOAM environment.

### Parsing Complex OpenFOAM Dictionaries

`mutators.py` mutates dictionaries in two tiers.

**Tier 1** is a line-based reader/writer. It is the primary path because it
returns and writes values *verbatim* -- `5e-6` stays `5e-6`. It handles block
comments, `#include`, and `#calc` correctly.

**Tier 2** is `foam_backend.py`, backed by foamlib. It is consulted only when
tier 1 cannot locate the target -- most commonly a brace inside a quoted value,
which defeats brace counting. foamlib parses in process and never evaluates
`#calc` or `#codeStream`.

Reads never reach tier 2: `read_foam_entry` and `read_foam_dict_block` return
verbatim source text, and foamlib returns typed values.

driverFOAM no longer shells out to the `foamDictionary` binary. Behaviour no
longer depends on whether OpenFOAM is sourced. If you are writing tools that
query these dictionaries, use the `mutators.py` API -- not `grep` or `sed`.

### Find past runs

```python
from openfoam_driver.core.runtime.run_discovery import list_runs
for manifest in list_runs("/path/to/runs/dir"):
    print(manifest["run_id"], manifest["status"], manifest["_manifest_path"])
```

## Known gaps

These are real limitations; the agent must not assume them:

- **Automatic retry** is mechanical and bounded. `run --strict` retries a step
  whose failure is classified *retryable* (currently `workflow_step_timeout`) up
  to its `retry_policy.max_attempts` (or the run's `default_max_attempts`), with
  exponential backoff (`retry_policy.backoff_seconds`). Fatal failures
  (`missing_artifacts`, exec errors, generic nonzero exit / FOAM FATAL ERROR) are
  never retried. Between retryable attempts the persisted `workflow_state.json` is
  kept resumable, so a crash during backoff resumes into another retry. It does
  **not** read logs to reclassify failures or mutate configuration between
  attempts (that is deferred). A *terminal*-failed saved state is still refused by
  `run --strict`; use `step --strict` to rerun it manually.

- **Environment preflight** is command-aware but not exhaustive. Strict planning
  derives the executables your plan will run from its `workflow_dag` steps and
  errors (`environment_diagnostics`, a first-class report field) if any are missing
  from `PATH`, if `WM_PROJECT_DIR` is unset, or if the plan is parallel but no
  `mpirun`/`mpiexec` is found. It warns on a partially-sourced environment
  (`WM_PROJECT_VERSION` / `FOAM_USER_LIBBIN` unset). It does **not** yet check free
  disk space or output-directory writability. Set `SKIP_ENV_DIAGNOSTICS=1` to bypass
  the gate (used by the test suite).

- **Active-tension models beyond NashPanfilov and GoktepeKuhl** are not in `active_tension_catalog.py`. Future C++ models must be registered there before artifact prediction will cover their state variables.

If your agent depends on any of these, expect failure and consider a workaround (e.g. starting from an existing tutorial template and overriding deltas rather than constructing from scratch).

## Where to read further

- `applications/scripts/driverFoam/openfoam_driver/dict_entries.py` — every dict key with its constraints
- `applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/ionic_model_catalog.py` — every ionic model
- `applications/scripts/driverFoam/openfoam_driver/core/utility_catalog.py` — every utility's CLI surface and outputs
- `applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/solver_coupling.py` — cross-domain coupler rules
- `applications/scripts/driverFoam/openfoam_driver/core/strict_planning.py` — strict preflight report and RunDocument v3 assembly
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/run_model.py` — RunDocument v3 model and explicit v1/v2 migration
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/workflow.py` — workflow DAG normalization and validation
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/workflow_state.py` — persisted step state model
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/workflow_runner.py` — low-level strict step executor
- `applications/scripts/driverFoam/schemas/run-document.json` — canonical RunDocument v3 JSON Schema

## Plugin selection (Phase 1)

`--plugin` accepts an installed plugin id from the `driverfoam.plugins`
entry-point group, a trusted `module.path:PluginClass` local-development
import (a colon always selects this form), or `none` for generic OpenFOAM.
The `capability_manifest` accept-surface is plugin-dependent:
`allowed_commands.core` lists solver-neutral OpenFOAM commands plus the
active plugin's own, so it changes with `--plugin`.

---

## Adding a New cardiacFoam Tutorial

This is for adding one more registered tutorial (a manufactured-solution
case, a benchmark, a new sweep-able configuration) that uses the cardiacFoam
solver family already wired up — not for adding support for a different
solver binary. For that, see "Plugin Guide — Adding a New Solver" below.

**The registry is the single source of truth:**
`openfoam_driver/plugins/cardiacfoam/tutorials/registry.py` holds
`SPEC_FACTORIES` (id, and its lowercase alias, → factory function) and
`REGISTERED_TUTORIALS` (the canonical id tuple). Both are exported through
`CardiacFoamPlugin.get_tutorial_catalog()`.

1. **Add an id** to the `CardiacTutorialID` enum in
   `openfoam_driver/plugins/cardiacfoam/tutorials/ids.py`, e.g.
   `MY_NEW_CASE = "myNewCase"`.
2. **Write `tutorials/my_new_case.py`** with a `make_spec(...) -> TutorialSpec`
   factory. `TutorialSpec` (`core/runtime/models.py`) needs `name`,
   `case_root`/`setup_root`/`output_dir` (build via the shared
   `resolve_spec_paths(...)` helper), `build_cases` (returns
   `list[CaseConfig]`), `apply_case` (mutates the case's dict files per
   `CaseConfig`, typically via `apply_electro_property_overrides`/
   `apply_physics_property_overrides` from `plugins/cardiacfoam/overrides.py`),
   and a `metadata` dict with at least a `workflow_dag` (a `solve` step at
   minimum). There is no per-tutorial output-collection callback to wire up
   — output discovery globs the case's actual on-disk files instead.
   `single_cell.py` is the smallest complete worked example of this shape.
3. **Do not accept any of the four dead postprocess-selector parameter
   names** — `cv_extract_script_relpath`, `postprocess_script_relpath`,
   `postprocess_function_name`, `table_summary_relpath`. The runtime cannot
   discover through them (see `test_tutorial_postprocessing_contract.py`,
   which fails the whole suite if any factory's signature advertises one).
4. **Register it** in `registry.py`: import your `make_spec`, add
   `CardiacTutorialID.MY_NEW_CASE.value` (and its `.lower()` form) to
   `SPEC_FACTORIES`, add the id to `REGISTERED_TUTORIALS`.
5. **Add a display entry** in `tutorials/display.py` (a `TutorialDisplay(...)`)
   — an exporter cross-checks this one-to-one against `REGISTERED_TUTORIALS`,
   so a tutorial without a display card (or a display card without a
   factory) fails.
6. **Extend the characterization fixture** —
   `openfoam_driver/tests/core/test_cardiac_tutorial_characterization.py`
   iterates every id in `REGISTERED_TUTORIALS` against
   `tests/fixtures/cardiac_tutorial_characterization.json`. Regenerate it
   once your factory exists (there's no dedicated CLI for this — write
   `{"tutorials": _current_characterization()}` back to the fixture path,
   `json.dumps(..., indent=2, sort_keys=True)`, matching the existing
   formatting).

Once registered, drive it exclusively through `driverFoam`
(plan/run/sweep) per `CLAUDE.md` — never a bespoke shell script.

---

## Plugin Guide — Adding a New Solver to driverFOAM

This section is for **plugin authors** — developers or AI agents who need to
add support for a new OpenFOAM solver to driverFOAM. End-users running existing
solvers do not need to read this section.

> **Quickest path:** Follow the dedicated skill at
> `../../../../.agents/skills/driverfoam-plugin-builder/SKILL.md`
> (relative to this file), which contains a complete step-by-step workflow,
> a worked `ShallowWaterPlugin` example, and a troubleshooting table.

### What a plugin is

A driverFOAM plugin is a Python class that implements the `SolverPlugin`
contract defined in `openfoam_driver/core/plugin_interface.py`. It creates a
clean boundary between the generic execution engine and all solver-specific
knowledge.

Two Protocol classes define the contract:

| Class | Members | Required when |
|---|---|---|
| `SolverPlugin` | 27 | Always |
| `SolverPluginOptionalHooks` | 14 (probe-based) | Never required; enable capabilities |

### Mandatory files

| File | Purpose |
|---|---|
| `my_solver_plugin.py` | Python class implementing the contract |
| `plugin.yaml` | Manifest: identity, case file rules, optional C++ roots |
| `pyproject.toml` entry-point | `[project.entry-points."driverfoam.plugins"]` |

### Required Members (all plugins)

```python
plugin_name             # str — human display name
plugin_id               # str — reverse-DNS id, must match plugin.yaml
plugin_version          # str — plugin semantics version
plugin_api_version      # str — "2", the only supported contract version
get_profile()           # PluginProfile from load_plugin_profile("plugin.yaml")
get_dict_entries()      # tuple[DictEntry, ...] — globally unique driver_paths
get_dictionary_catalog() # DictionaryCatalog — entries by document name
get_dict_groups()       # dict[str, tuple[DictEntry, ...]] — by logical group
get_capabilities()      # CapabilityManifest via build_capability_manifest()
get_tutorial_catalog()  # dict with spec_factories, registered_tutorials
get_tutorial_displays() # tuple[TutorialDisplay, ...]
validate_configuration(spec)   # tuple[StrictDiagnostic, ...]
validate_run_semantics(context) # tuple[...]
predict_data_artifacts(case_root, spec) # tuple[DataArtifact, ...]
get_solver_commands()           # frozenset[str] — artifact-producing binaries
get_auxiliary_commands()        # frozenset[str] — meshers, decomposers
get_utility_manifests()         # dict[str, Any]
get_utility_roots()             # tuple[Path, ...]
resolve_case_models(case_root)  # dict — best-effort, never raise
get_samplable_fields(resolved)  # dict[str, tuple[str, ...]] — by region
get_override_schema(tutorial, info) -> dict
get_run_document_config_schema() -> dict  # JSON Schema
get_dict_entry_catalog()        # dict — entries by document name (unserialized)
get_solve_step_commands()       # frozenset[str]
get_telemetry_source_globs(command) # tuple[str, ...]
get_extra_provenance_paths(case_root) # tuple[RuntimeDependency, ...]
get_artifact_value_reader(format)    # Any | None
```

`GenericOpenFOAMPlugin` (`core/generic_plugin.py`) is the canonical scaffold —
copy it and fill in identity properties; every required member already has a
neutral implementation to start from.

### Key Optional Hooks (`SolverPluginOptionalHooks`, probed with `getattr`)

Two have no neutral fallback — sweeps fail if they are absent:

| Hook | If absent |
|---|---|
| `route_sweep_case_values(...)` | **Sweeps refused by name** |
| `materialize_sweep_case(...)` | **Sweeps refused by name** |
| `has_case_marker(case_root)` | `False` |
| `is_nondimensional_case(spec)` | `False` (SI mesh checks on) |
| `build_run_document_config(spec)` | `({}, ())` |
| `get_override_scopes()` | `()` |
| `get_regeneration_scopes()` | `()` |
| `get_report_catalog()` | `()` |
| `get_named_catalogs()` | `{}` |

### Entry-point registration

```toml
[project.entry-points."driverfoam.plugins"]
mysolver = "my_package.my_solver_plugin:MySolverPlugin"

[tool.setuptools.package-data]
"my_package" = ["plugin.yaml"]
```

### Validation commands

```bash
# Verify entry-point is discoverable
python -c "from importlib.metadata import entry_points; \
           print(list(entry_points(group='driverfoam.plugins')))"

# Load and validate
python -c "
from openfoam_driver.core.plugin_interface import load_plugin_context
ctx = load_plugin_context('mysolver')
print('OK:', ctx.identity)
"

# Strict plan
driverFoam --plugin mysolver plan --strict --entry <tutorial_or_case_path>
```

### `validate_plugin()` cross-validation rules

- `profile.plugin_id` **must equal** `plugin.plugin_id`
- `profile.api_version` **must equal** `plugin.plugin_api_version`
- All `DictEntry.driver_path` values must be **globally unique**

### Common errors

| Error | Cause |
|---|---|
| `KeyError: 'mysolver'` | Wrong entry-point group or not installed |
| `TypeError: missing required members: X` | Missing v1 methods |
| `TypeError: missing v2 contract; missing: X` | Missing v2 callables |
| `TypeError: profile id does not match plugin_id` | YAML id ≠ class property |
| `TypeError: duplicate paths: X` | Two `DictEntry` share same `driver_path` |
| Sweep refused: `does not implement route_sweep_case_values` | Implement sweep hooks |

### See also

- `openfoam_driver/core/plugin_interface.py` — full Protocol definitions
- `openfoam_driver/core/generic_plugin.py` — minimal v2 scaffold to copy
- `openfoam_driver/core/generic-plugin.yaml` — annotated `plugin.yaml` template
- `openfoam_driver/plugins/cardiacfoam_plugin.py` — full v2 reference
- `KEY_FILES.md` — navigational map for all reader types
