# driverFOAM Agent Guide

This is the agent contract for planning, launching, polling, and inspecting
cardiacFoam runs through the Python `openfoam_driver` package. Read this once
before driving the orchestrator.

## What the agent can do

| Action | Function | Module |
|---|---|---|
| Discover tutorials, dict keys, ionic models, utilities | `describe_tutorial(...)` | `openfoam_driver.introspection` |
| Build a non-mutating strict launch contract | `strict_plan(...)` | `openfoam_driver.strict_planning` |
| Execute an agent-authored RunDocument | `foamctl run/step --run-document <file>`; `build_execution_inputs(...)` | `openfoam_driver.core.runtime.run_document_exec` |
| Execute one strict workflow step | `run_workflow_step(...)` | `openfoam_driver.core.runtime.workflow_runner` |
| Read/write strict workflow state | `workflow_state_from_json(...)`, `WorkflowRunState.to_json()` | `openfoam_driver.core.runtime.workflow_state` |
| Validate RunDocument v2 or migrate v1 explicitly | `RunDocument.from_json(...)`, `RunDocument.migrate_v1(...)` | `openfoam_driver.core.runtime.run_model` |
| Validate a configuration before launching | `validate_run(run, *, entries=None)` | `openfoam_driver.specs.validation` |
| Synthesize a fresh `electroProperties` / `physicsProperties` | `build_electro_properties(...)`, `build_physics_properties(...)` | `openfoam_driver.specs.dict_builder` |
| Parse an existing `electroProperties` back to selectors + overrides | `parse_electro_properties(path)` | `openfoam_driver.specs.dict_builder` |
| Build + launch a one-shot run (runs through the strict executor) | `build_and_launch(...)` | `openfoam_driver.specs.dict_builder` |
| Locate predicted outputs | Read `artifacts_manifest.json` (sidecar, atomic) | `<output_dir>/` |
| Verify outputs vs predictions | Read `artifacts_realized.json` (written at terminal status) | `<output_dir>/` |
| List past runs | `list_runs(root)` | `openfoam_driver.core.runtime.run_discovery` |
| Plan/run a parameter sweep | `foamctl sweep-plan/sweep-run --spec sweep.json --output-dir <dir>` | `openfoam_driver.core.runtime.sweep_runner` |

## Preferred strict agent loop

The solver-injection refactor does not change this public loop. A single
per-operation driver context now supplies focused solver capabilities
internally, while omitted contexts, RunDocument v2, legacy fallbacks, commands,
diagnostics, and artifacts retain their established behavior.

Use strict planning before launching. It is the only path that tells an agent
whether the run is machine-readable, validated, catalog-covered, artifact
predictable, and workflow-addressable before execution starts.

```bash
foamctl plan --strict --entry singleCell
foamctl run --strict --entry singleCell
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
foamctl step --strict --entry singleCell --step solve
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
foamctl run --strict --entry singleCell --fresh
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
from openfoam_driver.strict_planning import strict_plan

report = strict_plan("singleCell")
payload = report.to_json()
if payload["status"] != "ok":
    raise RuntimeError(payload)
print(payload["workflow_state"]["current_step_id"])
```

### Executing an agent-authored RunDocument

`plan --strict` emits a complete `run_document` (RunDocument v2) in its JSON
output. An agent can persist that document, edit it (e.g. tune `config`, add or
reorder `workflowDag` steps, set per-step `retry_policy`), and execute the
edited document directly — the driver runs *your* document instead of
regenerating one from `--entry`:

```bash
# 1. Plan and capture the run document the planner produced.
foamctl plan --strict --entry singleCell > plan.json
python3 -c "import json; json.dump(json.load(open('plan.json'))['run_document'], open('run.json','w'))"

# 2. (optional) edit run.json — config, workflowDag, retry_policy, expectedArtifacts.

# 3. Execute the document. No --entry; --strict is implied by the document.
foamctl run  --run-document run.json
foamctl step --run-document run.json --step solve   # single step
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
plugin-boundary compatibility fallbacks (v1 plugin support, legacy shims),
see `openfoam_driver/core/compatibility.py`.

## Compatibility one-shot loop

The legacy `build_and_launch(...)` path remains supported for existing scripts,
but it is not the preferred autonomous path because it mutates and launches in
one call instead of first emitting a strict contract.

```python
from openfoam_driver.specs.dict_builder import build_and_launch

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
foamctl sweep-plan --spec sweep.json --output-dir sweeps/my_sweep/
foamctl sweep-run --spec sweep.json --output-dir sweeps/my_sweep/
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
via `build_and_launch`. Some tutorials (`niederer2012`, `manufacturedFDA`, and
others under `openfoam_driver/specs/tutorials/`) instead expose their own
`make_spec(**kwargs)` with tutorial-specific parameters (e.g. `niederer2012`'s
`dx_values`/`dt_values`/`end_time_by_dx`, in millimetres/milliseconds;
`manufacturedFDA`'s `dimensions`/`number_cells`/`dt_values`). To sweep one of
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
own `apply_case()` mutate that tutorial's *shared* `case_root` in place
(confirmed for `niederer2012` and `manufacturedFDA`: they patch
`system/controlDict`/`system/blockMeshDict*` directly rather than writing an
isolated per-case directory). So each resolved axis combination must collapse
to exactly one case — if it doesn't (e.g. a config that still fans out
internally because a constraining kwarg like `solvers` is missing),
`sweep-plan`/`sweep-run` reports that case as `failed` with a clear
`materialization_error` rather than silently applying only the first of
several. In practice this means giving `dt`/`dx`-style axes their own
dedicated sweep row (`"zip"` mode with per-case single-element lists, as
above) instead of relying on the tutorial's own internal multi-value fan-out.
Because materialization mutates shared state, entry-mode sweeps must never be
parallelized across cases — `sweep_run`'s plain sequential loop already
guarantees this.

**`sweep-plan` is not fully non-mutating here.** Unlike generic mode (which
only ever writes into a fresh directory under `--output-dir`), entry-mode's
`apply_case()` mutates the tutorial's real, shared `case_root` — including
during `sweep-plan`. Don't run either action against a tutorial whose
`case_root` holds results you care about without first checking what's there
(or testing against a scratch copy with a `tutorials_root` override in `base`
pointed elsewhere).

Everything else — the manifest, `--retry-failed`, `--case-timeout-s`,
`--max-cases`, resumability — is identical to generic mode.

## Polling a long-running legacy run

For legacy engine runs that take minutes, prefer the async-friendly polling
pattern:

```python
import json
import time
from pathlib import Path

manifest_path = Path("<output_dir>/run_manifest.json")
while True:
    manifest = json.loads(manifest_path.read_text())
    if manifest["status"] in {"completed", "completed_with_failures", "failed",
                              "postprocess_failed"}:
        break
    time.sleep(15)
```

`run_manifest.json` is rewritten atomically (`os.replace` of a `.tmp` sibling),
so the read above is safe at any instant. Do not implement polling that opens
`.tmp` files directly.

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
foamctl step --strict --step <id> --apply overrides.json
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
   - **Note on entry paths**: `foamDictionary` uses `/` to traverse nested blocks. If a block name contains special characters (like the `Vm|VmFinal|u|uFinal` solver block), you **must** wrap that specific block name in quotes within the path: e.g., `system/electro/fvSolution:solvers/"Vm|VmFinal|u|uFinal"/tolerance`.
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

Legacy engine runs already write `artifacts_realized.json` at terminal status.
After such a run reaches a terminal status, agents should compare predicted
artifacts to on-disk reality:

```python
manifest = json.loads(Path("<output_dir>/run_manifest.json").read_text())
realized_path = manifest["artifacts_realized_path"]
realized = json.loads(Path(realized_path).read_text())

for case in realized["cases"]:
    print(case["case_id"], case["matched_count"], "/", case["predicted_count"])
    for artifact in case["artifacts"]:
        if artifact["status"] == "missing" and not artifact["optional"]:
            print("  warning missing required:", artifact["artifact_id"])
```

Missing-but-optional artifacts are not errors. They only appear under specific
configurations, for example probes that were not enabled.

## Discovering what's valid

Three layers of discovery:

1. **What tutorials exist?** `from openfoam_driver.introspection import describe_launch_matrix; describe_launch_matrix()` returns every registered entry.
2. **What dict keys can I set?** Iterate `openfoam_driver.dict_entries.ELECTRO_PROPERTY_ENTRY_GROUPS` and `PHYSICS_PROPERTY_ENTRIES` for case-physics entries. For time-control use `openfoam_driver.dict_entries.CONTROL_DICT_ENTRIES` (`deltaT`, `endTime`). Each entry carries `driver_path`, `value_kind`, `enum_values`, `unit`, `typical_value`, and structured constraints (`applicable_when`, `forbidden_when`, `required_when`, `mutually_exclusive_with`).
3. **What ionic models can I pick?** `from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG`. Each entry carries `states`, `algebraic`, `compatible_solvers`, `compatible_tissues`, `species`, `cardiac_region`, `recommended_exports`.
4. **What utilities are known?** `from openfoam_driver.utility_catalog import UTILITY_CATALOG`. Strict planning fails when a workflow command has missing required `produces` metadata.
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
from openfoam_driver.specs.dict_builder import parse_electro_properties

parsed = parse_electro_properties("/path/to/case/constant/electroProperties")
# {"selectors": {"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP", ...},
#  "overrides": {"$ELECTRO_MODEL_COEFFS.solutionAlgorithm": "explicit", ...}}
```

Pass the result directly to `build_electro_properties` to round-trip:

```python
from openfoam_driver.specs.dict_builder import build_electro_properties, parse_electro_properties

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

### Parsing Complex OpenFOAM Dictionaries (foamDictionary)

`mutators.py` implements a hybrid parsing architecture for all OpenFOAM dictionary mutations (`read_foam_entry`, `update_foam_entry`, `ensure_foam_dict`, `remove_foam_dict`).

If a target dictionary uses complex OpenFOAM C++ syntax (e.g., `#include` macros, `/* block comments */`, nested scopes, `#calc`), the naive Python regex parser may fail with `KeyError: unbalanced braces` or `KeyError: not found`.

To handle this, all mutator functions automatically fallback to using OpenFOAM's native `foamDictionary` C++ executable if it exists in the `PATH`.
If you are writing custom bash scripts or tools that need to query values from these complex dictionaries, do not rely on `grep` or `sed`. Instead, use the native CLI or the `mutators.py` API:

```bash
# Safely extract a value, ignoring comments and expanding macros
foamDictionary system/controlDict -entry functions/myFunction/type -value

# Safely modify a value inline
foamDictionary system/controlDict -entry startTime -set 0.0
```

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
- `applications/scripts/driverFoam/openfoam_driver/utility_catalog.py` — every utility's CLI surface and outputs
- `applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/solver_coupling.py` — cross-domain coupler rules
- `applications/scripts/driverFoam/openfoam_driver/strict_planning.py` — strict preflight report and RunDocument v2 assembly
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/run_model.py` — RunDocument v2 model and v1 migration
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/workflow.py` — workflow DAG normalization and validation
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/workflow_state.py` — persisted step state model
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/workflow_runner.py` — low-level strict step executor
- `applications/scripts/driverFoam/schemas/run-document.json` — canonical RunDocument v2 JSON Schema

## Plugin selection (Phase 1)

`--plugin` accepts an installed plugin id from the `driverfoam.plugins`
entry-point group, a trusted `module.path:PluginClass` local-development
import (a colon always selects this form), or `none` for generic OpenFOAM.
The `capability_manifest` accept-surface is plugin-dependent:
`allowed_commands.core` lists solver-neutral OpenFOAM commands plus the
active plugin's own, so it changes with `--plugin`.
