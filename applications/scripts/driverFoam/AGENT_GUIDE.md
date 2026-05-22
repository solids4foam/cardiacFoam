# driverFOAM Agent Guide

This is the agent contract for launching, polling, and inspecting cardiacFoam runs through the Python `openfoam_driver` package. Read this once before driving the orchestrator.

## What the agent can do

| Action | Function | Module |
|---|---|---|
| Discover tutorials, dict keys, ionic models, utilities | `describe_tutorial(...)`, `describe_launch_matrix()` | `openfoam_driver.introspection` |
| Validate a configuration before launching | `validate_run(run, *, entries=None)` | `openfoam_driver.specs.validation` |
| Synthesize a fresh `electroProperties` / `physicsProperties` | `build_electro_properties(...)`, `build_physics_properties(...)` | `openfoam_driver.specs.dict_builder` |
| Build + launch a one-shot run | `build_and_launch(...)` | `openfoam_driver.specs.dict_builder` |
| Run a registered tutorial | `DriverEngine(spec=..., requested_action=...).run_simulations()` | `openfoam_driver.core.runtime.engine` |
| Poll progress | Read `run_manifest.json` (atomic), tail `action_events.jsonl` (one JSON per line) | `<output_dir>/` |
| Locate predicted outputs | Read `artifacts_manifest.json` (sidecar, atomic) | `<output_dir>/` |
| Verify outputs vs predictions | Read `artifacts_realized.json` (written at terminal status) | `<output_dir>/` |
| List past runs | `list_runs(root)` | `openfoam_driver.core.runtime.run_discovery` |

## Minimum-viable agent loop

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

## Polling a long-running run

For runs that take minutes, prefer the async-friendly polling pattern:

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

`run_manifest.json` is rewritten atomically (`os.replace` of a `.tmp` sibling), so the read above is safe at any instant. Do not implement your own polling that opens `.tmp` files directly.

## Verifying outputs

After the run reaches a terminal status, agents should compare predicted artifacts to on-disk reality:

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

Missing-but-optional artifacts are not errors — they only appear under specific configurations (e.g. probes that weren't enabled).

## Discovering what's valid

Three layers of discovery:

1. **What tutorials exist?** `from openfoam_driver.introspection import describe_launch_matrix; describe_launch_matrix()` returns every registered entry.
2. **What dict keys can I set?** Iterate `openfoam_driver.dict_entries.ELECTRO_PROPERTY_ENTRY_GROUPS` and `PHYSICS_PROPERTY_ENTRIES`. Each entry carries `driver_path`, `value_kind`, `enum_values`, `unit`, `typical_value`, and structured constraints (`applicable_when`, `forbidden_when`, `required_when`, `mutually_exclusive_with`).
3. **What ionic models can I pick?** `from openfoam_driver.ionic_model_catalog import IONIC_MODEL_CATALOG`. Each entry carries `states`, `algebraic`, `compatible_solvers`, `compatible_tissues`, `species`, `cardiac_region`, `recommended_exports`.

## What the validator catches

`validate_run(run)` runs five families of checks:

- **Required fields** — every `required` entry has a value.
- **Enum membership** — values for enum-typed entries are in `enum_values`.
- **Structured constraints** — `applicable_when` / `forbidden_when` / `required_when` / `mutually_exclusive_with`.
- **Solver coupling** — pairings like (`bidomainSolver`, any Purkinje) reject with the table's stated reason.
- **Block references** — `domainCouplings.<name>.conductionNetworkDomain` must point at a declared block.

If the dict builder rejects your input with `ValueError`, the message lists every violation. Fix the selectors or overrides and call again.

## Common patterns

### Override a single dict key

```python
build_electro_properties(
    electro_selectors={"myocardiumSolver": "monodomainSolver",
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
    electro_selectors={"myocardiumSolver": "bidomainSolver",
                       "ionicModel": "bathBidomainFDAManufactured"},
    overrides={
        "$ELECTRO_MODEL_COEFFS.bathPotentialDomain.bathCellZones": "(bath organ)",
        "$ELECTRO_MODEL_COEFFS.bathPotentialDomain.heartCellZone": "myocardium",
    },
)
```

Declaring any `bathPotentialDomain.*` override auto-enables the bath block — the bath leaves typical-value default unless overridden.

### Find past runs

```python
from openfoam_driver.core.runtime.run_discovery import list_runs
for manifest in list_runs("/path/to/runs/dir"):
    print(manifest["run_id"], manifest["status"], manifest["_manifest_path"])
```

## What's not yet supported

These are known gaps; the agent must not assume them:

- **Active-tension model variables** are predicted for `NashPanfilov` and `GoktepeKuhl` when an `activeTensionModel` block is declared in `electroProperties`. Future C++ models must be added to `active_tension_catalog.py` first.
- **Reverse parsing** of an existing `electroProperties` is available via `parse_electro_properties(path)` in `openfoam_driver.specs.dict_builder`. Returns `{"selectors": {...}, "overrides": {...}}` that round-trips through `build_electro_properties`.
- **Bidomain + Purkinje coupling** is currently rejected by the validator — the `bidomainPvjCoupler` C++ class does not exist yet.
- **Per-step retry / checkpointing through the workflow DAG** is not implemented — the DAG is descriptive metadata.

If your agent depends on any of these, expect failure modes and consider the workaround (e.g. start from an existing tutorial template and override the deltas, rather than constructing from scratch).

## Where to read further

- `applications/scripts/driverFoam/openfoam_driver/dict_entries.py` — every dict key with its constraints
- `applications/scripts/driverFoam/openfoam_driver/ionic_model_catalog.py` — every ionic model
- `applications/scripts/driverFoam/openfoam_driver/utility_catalog.py` — every utility's CLI surface and outputs
- `applications/scripts/driverFoam/openfoam_driver/solver_coupling.py` — cross-domain coupler rules
- `docs/superpowers/plans/2026-05-19-driverfoam-agentic-integration.md` — the architecture rationale behind all of the above
