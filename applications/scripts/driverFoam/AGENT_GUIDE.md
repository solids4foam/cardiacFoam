# driverFOAM Agent Guide

This is the agent contract for planning, launching, polling, and inspecting
cardiacFoam runs through the Python `openfoam_driver` package. Read this once
before driving the orchestrator.

## What the agent can do

| Action | Function | Module |
|---|---|---|
| Discover tutorials, dict keys, ionic models, utilities | `describe_tutorial(...)`, `describe_launch_matrix()` | `openfoam_driver.introspection` |
| Build a non-mutating strict launch contract | `strict_plan(...)` | `openfoam_driver.strict_planning` |
| Execute one strict workflow step | `run_workflow_step(...)` | `openfoam_driver.core.runtime.workflow_runner` |
| Read/write strict workflow state | `workflow_state_from_json(...)`, `WorkflowRunState.to_json()` | `openfoam_driver.core.runtime.workflow_state` |
| Validate RunDocument v2 or migrate v1 explicitly | `RunDocument.from_json(...)`, `RunDocument.migrate_v1(...)` | `openfoam_driver.core.runtime.run_model` |
| Validate a configuration before launching | `validate_run(run, *, entries=None)` | `openfoam_driver.specs.validation` |
| Synthesize a fresh `electroProperties` / `physicsProperties` | `build_electro_properties(...)`, `build_physics_properties(...)` | `openfoam_driver.specs.dict_builder` |
| Parse an existing `electroProperties` back to selectors + overrides | `parse_electro_properties(path)` | `openfoam_driver.specs.dict_builder` |
| Build + launch a legacy one-shot run | `build_and_launch(...)` | `openfoam_driver.specs.dict_builder` |
| Run a registered tutorial through the legacy engine | `DriverEngine(spec=..., requested_action=...).run_simulations()` | `openfoam_driver.core.runtime.engine` |
| Poll legacy engine progress | Read `run_manifest.json` (atomic), tail `action_events.jsonl` (one JSON per line) | `<output_dir>/` |
| Locate legacy predicted outputs | Read `artifacts_manifest.json` (sidecar, atomic) | `<output_dir>/` |
| Verify legacy outputs vs predictions | Read `artifacts_realized.json` (written at terminal status) | `<output_dir>/` |
| List past runs | `list_runs(root)` | `openfoam_driver.core.runtime.run_discovery` |

## Preferred strict agent loop

Use strict planning before launching. It is the only path that tells an agent
whether the run is machine-readable, validated, catalog-covered, artifact
predictable, and workflow-addressable before execution starts.

```bash
foamctl plan --strict --entry singleCell
foamctl run --strict --entry singleCell
```

The `plan --strict` command is non-mutating. It prints JSON with:

- `status`: `ok` or `failed`
- `resolved_entry`: case/spec identity and paths
- `validation_diagnostics`: RunDocument and configuration validation results
- `catalog_coverage_errors`: strict dict-key coverage failures
- `artifact_diagnostics`: solver/utility/artifact prediction coverage failures
- `environment_diagnostics`: missing executables, unsourced OpenFOAM env, missing MPI launcher
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

Programmatic planning uses the same contract:

```python
from openfoam_driver.strict_planning import strict_plan

report = strict_plan("singleCell")
payload = report.to_json()
if payload["status"] != "ok":
    raise RuntimeError(payload)
print(payload["workflow_state"]["current_step_id"])
```

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
print(result)  # {"case_dir": ..., "status": "complete", "results": [...]}
```

That single call:

1. Calls `build_electro_properties(...)` and `build_physics_properties(...)`.
2. Runs the validator on both — raises `ValueError` if your selectors break a structured constraint.
3. Writes `case/constant/electroProperties` and `case/constant/physicsProperties`.
4. Constructs a `generic_case` spec pointing at the case directory.
5. Launches `DriverEngine.run_simulations()`.
6. Returns the per-case results.

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
3. **What ionic models can I pick?** `from openfoam_driver.ionic_model_catalog import IONIC_MODEL_CATALOG`. Each entry carries `states`, `algebraic`, `compatible_solvers`, `compatible_tissues`, `species`, `cardiac_region`, `recommended_exports`.
4. **What utilities are known?** `from openfoam_driver.utility_catalog import UTILITY_CATALOG`. Strict planning fails when a workflow command has missing required `produces` metadata.
5. **What dict keys have parser limitations?** Read `openfoam_driver/scripts/dict_key_allowlist.json`. Strict dict-key scanning fails when new uncatalogued keys appear, stale catalog paths remain, or allowlist entries become unused.

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

### Configure a bath bidomain run

```python
build_electro_properties(
    selectors={"myocardiumSolver": "bidomainSolver",
               "ionicModel": "bathBidomainFDAManufactured"},
    overrides={
        "$ELECTRO_MODEL_COEFFS.bathPotentialDomain.bathCellZones": "(bath organ)",
        "$ELECTRO_MODEL_COEFFS.bathPotentialDomain.heartCellZone": "myocardium",
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
- `applications/scripts/driverFoam/openfoam_driver/ionic_model_catalog.py` — every ionic model
- `applications/scripts/driverFoam/openfoam_driver/utility_catalog.py` — every utility's CLI surface and outputs
- `applications/scripts/driverFoam/openfoam_driver/solver_coupling.py` — cross-domain coupler rules
- `applications/scripts/driverFoam/openfoam_driver/strict_planning.py` — strict preflight report and RunDocument v2 assembly
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/run_model.py` — RunDocument v2 model and v1 migration
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/workflow.py` — workflow DAG normalization and validation
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/workflow_state.py` — persisted step state model
- `applications/scripts/driverFoam/openfoam_driver/core/runtime/workflow_runner.py` — low-level strict step executor
- `applications/scripts/driverFoam/schemas/run-document.json` — canonical RunDocument v2 JSON Schema
- `docs/superpowers/plans/2026-05-19-driverfoam-agentic-integration.md` — historical architecture rationale
