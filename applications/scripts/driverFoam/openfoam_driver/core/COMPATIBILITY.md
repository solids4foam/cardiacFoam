# Plan 1 compatibility boundaries

Plan 1 changes dependency ownership without changing behavior. The items below
remain active, emit no new warnings, and retain their current precedence. Each
has a named boundary so a later behavior-changing plan can address it without
rediscovering the historical reason.

| Boundary | Why and exact activation | Preserved by | Plan 2 seam |
|---|---|---|---|
| `legacy_default_driver_context` | Public CLI/Python operation receives no context; select built-in cardiacFoam. | plugin-context, CLI matrix, strict-plan, validation tests | Default plugin selection and deprecation |
| `legacy_case_marker` / `legacy_case_runnable_without_workflow` | Plugin lacks case hooks; recognize `electroProperties*`, then require physics and standard system dictionaries when no contract/Allrun establishes runnability. | case compatibility matrix, filesystem ingest, RunDocument execution | Generic pipeline discovery |
| `legacy_generic_case_mutation` | Direct core `make_spec` uses electro/physics override arguments. | cardiac generic-case and template tests | Explicit generic mutation contract |
| `legacy_route_sweep_case` / `legacy_materialize_sweep_case` | Plugin lacks sweep hooks; use the historical build-and-launch cardiac sweep shape. | sweep routing, materialization, runner, manifest tests | Solver-neutral sweep schema |
| Allrun workflow fallback | A generic spec has no explicit solver command; normalize one `Allrun` step. | generic-case and filesystem-ingest tests | Manual/managed pipeline precedence |
| RunDocument v1 migration | Loader receives `version: "1"`; map conservatively to v2 fields. | RunDocument model/execution tests | RunDocument v3 and mandatory replanning |
| `build_and_launch` | Legacy caller requests one-shot cardiac dictionary synthesis and execution. | dict-builder, workflow, and agent-guide contracts | Canonical materialization plan |
| Flat controlDict override | Override path contains neither a document prefix nor `:`; validate and apply it as a controlDict entry. | apply-overrides and remediation tests | Explicit document/path mutations |
| Native-to-text dictionary fallback | `foamDictionary` is absent or its guarded operation fails; use the existing limited parser/update path. | mutator backend and integrity tests | Backend pinning and failure policy |
| Typical-value population | Applicable catalog entry is absent from explicit selectors/overrides and supplies `typical_value`. | dict-builder default-precedence tests | Default provenance/strict-explicit mode |
| Repository-layout search | Monorepo, tutorials sibling, or packaged fixture is needed and the preferred root is unavailable. | path, standalone, scanner, and verification tests | Explicit roots in agent plans |
| Historical cardiac artifact layout | Canonical output is absent; tutorial artifact locator checks processor/legacy filenames. | cardiac artifact and manufactured-tutorial tests | Minimum supported solver version |
| Verification `legacy_bash` | Experiment declares `legacy_bash`; execute the repository reproduction script. | verification-contract tests | Canonical workflow-DAG migration |

Only `openfoam_driver.core.compatibility` may import cardiac implementation
modules from core. Other core consumers use `DriverContext.capabilities`.
