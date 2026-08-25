# driverFOAM changelog

This file starts with the Phase 2 core-decoupling branch. Earlier history is
not reconstructed here — use `git log` for anything before it.

Only changes that break a **Python import surface** or an **on-disk document
format** are listed. Behaviour-preserving refactors, new tests, and
documentation are out of scope; `git log` is the record for those.

## Unreleased — ROADMAP Phase 2: core/plugin decoupling (branch `ep-work-onto-main`)

Plan: `docs/superpowers/plans/2026-08-16-driverfoam-roadmap-phase2-core-decoupling.md`.
Ledger: `.superpowers/sdd/progress.md`, section "driverFOAM ROADMAP Phase 2".

The theme: cardiacFoam-specific vocabulary moved out of the generic
`openfoam_driver` core and into `openfoam_driver.plugins.cardiacfoam`. Every
individual change was verified to have zero remaining in-tree importers at the
time it landed; this entry exists for **out-of-tree** consumers, for whom none
of that verification applies.

### Dictionary mutation and dependencies

Plan: `docs/superpowers/plans/2026-08-20-driverfoam-foamlib-tier2-backend.md` and
`docs/superpowers/plans/2026-08-21-driverfoam-foamlib-parser-migration-phase2.md`
(both gitignored, local-only — not the ROADMAP Phase 2 plan this section is
nested under; noted here since the two pieces of work landed on the same
branch at the same time).

- Dictionary mutation no longer shells out to `foamDictionary`. The complex-syntax
  fallback is now foamlib, in process. Written bytes no longer depend on whether
  OpenFOAM is sourced.
- **`update_foam_entry_via_foamDictionary`, `remove_foam_dict_via_foamDictionary`,
  and `ensure_foam_dict_via_foamDictionary` are deleted.** Any out-of-tree importer
  of these three names breaks; there is no replacement to import — the fallback
  they provided is now internal to `mutators.py`/`foam_backend.py`, not a
  separately callable public function.
- **`update_foam_entry` now raises `ValueError` on override values it previously
  accepted.** A value containing `;`, a newline, or `#` is rejected before it
  reaches a case dictionary (closes a documented value-channel injection risk —
  see `SECURITY.md`). A caller that was relying on such a value being written
  verbatim will now get an exception instead.
- **Provenance note:** in a *sourced* environment, mutated dictionary bytes change
  from `foamDictionary`'s re-serialisation to the line tier's form, so sha256
  provenance digests move for anyone who previously ran sourced. No golden digest
  baselines are committed, so nothing in-tree breaks; archived Paper I provenance
  JSON will not match a fresh re-run.
- Python floor raised to 3.11 (foamlib requirement). `numpy` is now a core
  dependency rather than a `[post]` extra.
- Four more hand-rolled OpenFOAM parsers (`specs/function_object_fields.py`,
  `plugins/cardiacfoam/detection.py`, `specs/mesh_geometry.py`,
  `plugins/cardiacfoam/run_document_config.py`) now read through foamlib instead
  of manual brace/paren-depth counting. No import-surface or on-disk format
  change; noted here because it's part of the same migration.

### Document format

- **`RunDocument` is now version `"3"`.** `RunDocument.from_json()` rejects
  both `"1"` and `"2"` documents with a `ValueError` instead of silently
  accepting them. Migrate explicitly: `RunDocument.migrate_v1(data)` or the
  new `RunDocument.migrate_v2(data)`. `to_json()` emits only `"3"`.
- **`schemas/run-document.json`'s `config` is now open**
  (`{"type": "object", "additionalProperties": true}`). Core no longer
  declares the `anatomy`/`physics`/`stimulus`/`solver` phases or the
  `physicsSlice` shape; a plugin declares its own config schema via
  `SolverPlugin.get_run_document_config_schema()`. The cardiac declaration
  lives in `plugins/cardiacfoam/config_schema.py`. Core validates the config
  against it on both the plan-emission and the document-ingestion path,
  reporting a `plugin_config_schema_violation` diagnostic.
- `openfoam_driver/schemas/run-document.json` is now **generated** from the
  hand-authored `schemas/run-document.json` by
  `schemas/generate_run_document_schema.py`. Edit the hand-authored copy only.

### Modules moved (import path changed, contents otherwise unchanged)

| Was | Now |
|---|---|
| `openfoam_driver.specs.detection` | `openfoam_driver.plugins.cardiacfoam.detection` |
| `openfoam_driver.specs.overrides` | `openfoam_driver.plugins.cardiacfoam.overrides` |
| `openfoam_driver.specs.system_templates` | `openfoam_driver.plugins.cardiacfoam.system_templates` |
| `openfoam_driver.strict_planning` | `openfoam_driver.core.strict_planning` |
| `openfoam_driver.planning_types` | `openfoam_driver.core.planning_types` |
| `openfoam_driver.capability_manifest` | `openfoam_driver.core.capability_manifest` |
| `openfoam_driver.introspection` | `openfoam_driver.core.introspection` |
| `openfoam_driver.report_catalog` | `openfoam_driver.core.report_catalog` |
| `openfoam_driver.utility_catalog` | `openfoam_driver.core.utility_catalog` |
| `openfoam_driver.tutorial_contracts` | `openfoam_driver.core.tutorial_contracts` |
| `openfoam_driver.tutorials_display` | `openfoam_driver.core.tutorials_display` |
| `openfoam_driver.sweep_derivation_catalog` | `openfoam_driver.core.sweep.sweep_derivation_catalog` |
| `openfoam_driver.sweep_expansion` | `openfoam_driver.core.sweep.sweep_expansion` |

These nine modules were solver-agnostic already; this only finishes moving
them into `core/` alongside the rest of the generic engine. `dict_entries.py`,
`sweep_materialize.py`, and `sweep_routing.py` deliberately stay at the
package root instead of joining them: each carries a lazily-imported (PEP 562
`__getattr__`, or a `_*_legacy` function-body import) fallback into
`plugins.cardiacfoam`, and `core/` is enforced by
`test_plugin_dependency_boundary.py` to import cardiac code only through
`compatibility.py`. `cli.py` also stays at the package root — the entry
surface, not engine internals.

The whole `specs/` package also moved under `core/` as `core/specs/`, keeping
its own identity (like `core/runtime/`, `core/contracts/`, `core/sweep/`)
rather than scattering its twelve modules flat:

| Was | Now |
|---|---|
| `openfoam_driver.specs.apply_overrides` | `openfoam_driver.core.specs.apply_overrides` |
| `openfoam_driver.specs.common` | `openfoam_driver.core.specs.common` |
| `openfoam_driver.specs.dict_builder` | `openfoam_driver.core.specs.dict_builder` |
| `openfoam_driver.specs.function_object_fields` | `openfoam_driver.core.specs.function_object_fields` |
| `openfoam_driver.specs.mesh_geometry` | `openfoam_driver.core.specs.mesh_geometry` |
| `openfoam_driver.specs.mesh_provisioning` | `openfoam_driver.core.specs.mesh_provisioning` |
| `openfoam_driver.specs.paths` | `openfoam_driver.core.specs.paths` |
| `openfoam_driver.specs.spatial_pacing` | `openfoam_driver.core.specs.spatial_pacing` |
| `openfoam_driver.specs.tet_mesh_provisioning` | `openfoam_driver.core.specs.tet_mesh_provisioning` |
| `openfoam_driver.specs.utils` | `openfoam_driver.core.specs.utils` |
| `openfoam_driver.specs.validation` | `openfoam_driver.core.specs.validation` |
| `openfoam_driver.specs.validation_types` | `openfoam_driver.core.specs.validation_types` |

The bundled fixture also moved: `openfoam_driver/specs/fixtures/template/...`
is now `openfoam_driver/core/specs/fixtures/template/...` (pyproject.toml
package-data entry updated to match). While moving `paths.py`, its dormant
Tier-3 `parents[3]` fallback (only reachable for a standalone install outside
this monorepo, so never exercised by this repo's own test suite) turned out
to have been off by one at the pre-move depth already — the move happened to
land it on the correct index by coincidence; see the code comment for detail.
ARCHITECTURE.md's stale claim that `apply_overrides.py` still imports
`plugins.cardiacfoam.detection` directly was also corrected while moving it:
it only imports `core/compatibility.py`.

The cardiac halves of three more modules were split out, leaving the generic
half at the original path:

| Symbol(s) | Was | Now |
|---|---|---|
| `resolve_context`, `build_electro_properties`, `parse_electro_properties`, `build_physics_properties`, `build_and_launch`, `_serialize`, `_entry_scope_and_key`, `_COEFFS_PREFIX` | `specs.dict_builder` | `plugins.cardiacfoam.dict_builder` |
| `MESHLESS_SOLVERS`, `BLOCK_MESH_SOLVERS`, `provision_mesh`, the `single_cell_polymesh` fixtures | `specs.mesh_provisioning` | `plugins.cardiacfoam.mesh_provisioning` |
| `read_purkinje_graph_bbox`, `discover_purkinje_graphs` | `specs.mesh_geometry` | `plugins.cardiacfoam.mesh_geometry` |

### Names removed

- **`openfoam_driver.specs.common` dropped 16 re-exports** from its facade and
  its `__all__`: `detect_myocardium_solver_name`,
  `detect_electro_coeffs_scope`, `detect_ionic_model_name`,
  `detect_ionic_export_list`, `electro_properties_has_block`,
  `detect_verification_model_type`, `detect_active_tension_model_name`,
  `detect_active_tension_export_list`, `_resolve_scope_tokens`,
  `normalize_entry_overrides`, `apply_entry_overrides`,
  `apply_electro_property_overrides`, `apply_physics_property_overrides`,
  `remove_electro_property_dict`, `ensure_electro_property_dict`, plus the
  `_IONIC_EXPORT_RE` / `_BLOCK_DECL_RE` / `_AT_EXPORT_RE` regexes. Import them
  from `plugins.cardiacfoam.detection` / `plugins.cardiacfoam.overrides`
  directly. `specs.common` is now genuinely generic (paths + utils only).
- **`openfoam_driver.report_catalog.REPORTS` is gone.** The concrete catalog
  moved to `plugins.cardiacfoam.reports.CARDIAC_REPORTS`; reach the active
  plugin's catalog through
  `driver_context.capabilities.report_catalog.reports()`.
  `report_catalog` keeps only the shared machinery (`ReportDefinition`, the
  `applicable_when` evaluator, the JSON record shape).
  `scripts/export-report-catalog.py` gained an optional `--plugin` flag; its
  default output is unchanged.
- **`openfoam_driver.tutorial_contracts` lost 3 deprecated constants**:
  `CORE_REQUIRED_FILES`, `SOLVER_REQUIRED_FILES`, `CONDITIONAL_FILES`. The
  authoritative source is `driver_context.capabilities.case_files`.
- `openfoam_driver.dict_entries`'s module-level cardiac import became a
  PEP 562 lazy `__getattr__` shim. Attribute access still works; a
  `from ... import *`-style eager expectation may not.

### Signatures changed

- `specs.dict_builder.select_applicable_entries(context, *, entries)` — the
  `entries` keyword is now **required** (was `entries: list[DictEntry] | None
  = None`, defaulting to the cardiac catalog). The cardiac-defaulting wrapper
  of the same name lives at
  `plugins.cardiacfoam.dict_builder.select_applicable_entries`.
- `core.runtime.generic_case.make_spec(...)` — the four cardiac-named
  parameters `electro_property_overrides`, `physics_property_overrides`,
  `electro_properties_relpath`, `physics_properties_relpath` are replaced by
  the neutral `dict_file_overrides` / `dict_file_relpaths` mappings, keyed by
  whatever names a plugin gives its own dictionaries. The old names still work
  as **deprecated aliases**, translated by
  `core.compatibility.legacy_generic_case_dict_file_aliases()`. The same
  rename applies to `TutorialSpec.metadata` keys.

### Output shape changed

- `introspection.describe_tutorial(...)` — the top-level
  `ionic_model_catalog` and `active_tension_catalog` keys are now nested under
  a single `plugin_catalogs` key, sourced from the plugin's
  `get_named_catalogs()` capability. Non-cardiac plugins get an empty mapping
  rather than cardiac-shaped placeholders.
- Strict-plan audit success text and `samplable_fields` defaults no longer
  name cardiac dictionaries or the `electro`/`solid` regions; both are
  plugin-sourced, with neutral empty defaults.

### Known open (not delivered by this branch)

- **P2.5-followup: `$ELECTRO_MODEL_COEFFS` scope resolver.** The sentinel is
  still hardcoded in `specs/validation.py` and `scripts/_dict_keys_scanner.py`
  and re-parsed by `specs/apply_overrides.py`. The plan's own exit-gate grep
  returns 12 hits, not zero. See the correction note at the end of
  `.superpowers/sdd/progress.md`'s Phase 2 section.
