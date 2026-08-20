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

- Dictionary mutation no longer shells out to `foamDictionary`. The complex-syntax
  fallback is now foamlib, in process. Written bytes no longer depend on whether
  OpenFOAM is sourced.
- **Provenance note:** in a *sourced* environment, mutated dictionary bytes change
  from `foamDictionary`'s re-serialisation to the line tier's form, so sha256
  provenance digests move for anyone who previously ran sourced. No golden digest
  baselines are committed, so nothing in-tree breaks; archived Paper I provenance
  JSON will not match a fresh re-run.
- Python floor raised to 3.11 (foamlib requirement). `numpy` is now a core
  dependency rather than a `[post]` extra.

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
