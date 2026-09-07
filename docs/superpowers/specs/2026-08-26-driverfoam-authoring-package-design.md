# driverFOAM: decouple core, then ship a solver-integration package

**Date:** 2026-08-26
**Status:** Phases 0-3 IMPLEMENTED and verified 2026-08-26 (14 commits,
`f6ee44d2..d0ebb9df` on ep-work-onto-main; suite 1705 passed / 9 skipped / 0
failed). Phases 4-5 (the `authoring/` package and its three verbs) still need
their own plan -- see §7-§11.

Evidence for the phases that landed:
- core has zero implicit context resolution (static AST guard, mutation-verified)
- a generic plan resolves the cardiac default zero times, and loads no cardiac
  module at all
- a plugin declaring its own phases gets full required-field and enum
  validation; the silent-skip test is mutation-verified
- cardiacFoam reaches 15/15 optional hooks
- `core/compatibility.py` cardiac references: 21 -> 2, both documented; the 19
  removed branches were instrumented and shown unreachable across the full
  suite before deletion
**Scope:** `applications/scripts/driverFoam/`
**Predecessor:** `2026-08-19-driverfoam-plugin-seam-documentation-design.md` (COMPLETE).
That spec's §5 deliberately deferred *"any plugin-author guide, out-of-tree
reference plugin, or CI portability test."* This spec delivers the first, plus
the core decoupling that makes it honest. The out-of-tree reference plugin
(solids4foam — `future/driverFOAM/SOLIDS4FOAM_CASE_STUDY.md`) stays deferred.

**Every measurement below was taken on 2026-08-26 against `ep-work-onto-main`.
None is inherited from a prior session.**

---

## 1. Problem

Two problems, one cause.

A developer who wants to drive a new OpenFOAM solver can load a plugin and then
do nothing else with the shipped tooling. And core, beneath that tooling, still
defaults to cardiacFoam tens of thousands of times per test run.

### 1.1 The contract works; everything after it does not

Established by experiment, not reading. `core/generic_plugin.py` was copied to a
fresh package, the four identity properties edited exactly as
`.agents/skills/driverfoam-plugin-builder/SKILL.md` instructs:

```
LOADED: org.myproject.shallowwater 2 sha256:3dbf9
```

Then, attempting to produce that plugin's machine-readable catalogs:

```
$ python scripts/export-dict-catalog.py      --out X --plugin my_package:ShallowWaterPlugin
export-dict-catalog.py: error: unrecognized arguments: --plugin
$ python scripts/export-tutorials-catalog.py --out X --plugin ...
export-tutorials-catalog.py: error: unrecognized arguments: --plugin
$ python scripts/export-utility-catalog.py   --out X --plugin ...
export-utility-catalog.py: error: unrecognized arguments: --plugin
```

The three that refuse `--plugin` are the three that bypass the plugin contract:

| script | how it bypasses the contract |
|---|---|
| `export-dict-catalog.py` | imports `plugins.cardiacfoam.ionic_model_catalog` and `.active_tension_catalog` directly; hardcodes `PHASES = ("anatomy","physics","stimulus","solver")` |
| `export-tutorials-catalog.py` | calls `default_driver_context()` |
| `export-utility-catalog.py` | reads `core.utility_catalog.UTILITY_CATALOG`, a module-level constant built from `Path(__file__).resolve().parents[4] / "utilities"` |

Each has a contract member built for the job and unused:
`get_dict_entry_catalog()`, `get_tutorial_catalog()`, `get_utility_manifests()`
/ `get_utility_roots()`.

Two of the three commit anti-patterns the project's own plugin-builder skill
lists by name: the `anatomy/physics/stimulus/solver` tuple is its anti-pattern #1
("cardiac key leak"), and `default_driver_context()` is its anti-pattern #5.

### 1.2 Nothing consumes the exports

No catalog JSON is committed anywhere (`git ls-files`, checked 2026-08-26). The
only consumers of all five exporters are their own conformance tests, which run
them by subprocess and assert on shape. They are write-only tools kept alive by
tests that exist to test them.

### 1.3 The guide teaches a contract that was deleted

`4da7db4c` (2026-08-25 12:34) — *"driverFOAM: collapse SolverPlugin v1/v2 into a
single required contract"*. `SKILL.md` was last modified 17:31 the same day and
still teaches the removed model:

| SKILL.md claims | code |
|---|---|
| `SolverPluginV2`, `_REQUIRED_V2_MEMBERS` | neither symbol exists |
| "14 v1 members" + "13 v2 members" as separate contracts | one flat `_REQUIRED_PLUGIN_MEMBERS`, 27 entries |
| troubleshooting: `Use "1" or "2" only` | `SUPPORTED_PLUGIN_API_VERSIONS = frozenset({"2"})` — following the doc gets you refused |
| error text `SolverPlugin declares plugin_api_version '2' but does not implement the v2 contract; missing: ...` | actual text is `SolverPlugin does not implement the plugin contract; missing: ...` |

The same stale claim appears in `core/generic-plugin.yaml`
(`api_version: "1" or "2"`) and `core/plugin_interface.py:287`
(docstring cites `_REQUIRED_V2_MEMBERS`).

The guide is also not in the repository. `.gitignore:121` ignores `.agents/`;
`driverfoam-assistant/SKILL.md` was force-added and is tracked,
`driverfoam-plugin-builder/SKILL.md` was not. The project's only
solver-integration guide is an untracked file describing a deleted contract.

Scaffold and guide also disagree: `generic_plugin.get_tutorial_catalog()` returns
`{"registered_tutorials": (), "spec_factories": {}}`; the guide's template adds
`"make_generic_case_spec": None`.

### 1.4 `--help` advertises the wrong solver's tutorials

```
$ driverFoam --plugin shallowwater --help
  --entry ENTRY   Entry name or relative workflow/case path to run
                  (singleCell, cable1DCVConvergence, niederer2012, ...)
```

`cli.py build_parser()` interpolates `list_tutorials()` into `--entry`'s help at
parser-construction time, before `--plugin` is parsed.

### 1.5 Phase vocabulary is core-owned, and unrecognised phases fail silently

Phases live in four places. Only the assignment is plugin-owned:

| concern | location | plugin-owned |
|---|---|---|
| vocabulary and order | `core/runtime/run_model.py:44` — `Phase = Literal["anatomy","physics","stimulus","solver"]` | no |
| per-entry shape | `core/contracts/dictionary.py` — `phases: frozenset[str]` | free-form at runtime; the `Literal` is a hint only |
| assignment | `plugins/cardiacfoam/dict_entries_catalog.py` (54 physics, 18 solver, 17 stimulus, 1 anatomy, 4 multi) and `common_dict_entries.py` (9 solver, 1 physics) | yes |
| meaning | `core/specs/validation.py` — `_PHASE_ORDER = get_args(Phase)`, `primary_phase()`, `_slice_value()` | no |

`primary_phase()` returns `None` for any phase outside the core `Literal`, and
its five call sites split two ways:

```
validation.py:353          ph = primary_phase(e); if ph is None: continue   # required-field check SKIPPED
validation.py:371          ph = primary_phase(e); if ph is None: continue   # enum check SKIPPED
validation.py:412          ph = primary_phase(e) or "physics"
dict_builder.py:239        ph = primary_phase(e) or "physics"
run_document_config.py:96  phase = primary_phase(entry_obj) or "physics"
```

A plugin declaring its own phase words gets required-field and enum validation
**silently skipped**, and its entries written into a `"physics"` slice it never
declared. Both failures are silent — the only silent-wrong defect in this set.

`RunDocument.config` is already plugin-neutral: `schemas/run-document.json`
declares it `{"type":"object","additionalProperties":true}` with *"Core imposes
no required keys or shape"*. The envelope concern in
`SOLIDS4FOAM_CASE_STUDY.md` is therefore already resolved; `Phase` is what remains.

### 1.6 Core defaults to cardiacFoam 51,540 times per test run

Every `legacy_*` function in `core/compatibility.py` was wrapped with a counter
and the full suite run (1690 tests). Invocation census:

```
51540  legacy_default_driver_context
   77  legacy_describe_config_resolution
   77  legacy_generic_case_dict_file_relpaths
   38  legacy_case_marker
    5  legacy_nondimensional_case
    5  legacy_generic_case_alias_names
    3  legacy_generic_case_dict_file_aliases
    2  legacy_case_runnable_without_workflow, legacy_route_sweep_case,
       legacy_resolve_case_models, legacy_run_document_config
    1  legacy_dict_regeneration_scopes, legacy_materialize_sweep_case,
       legacy_samplable_fields, legacy_generic_case_mutation,
       legacy_report_catalog
```

The dominant coupling is not the gated fallbacks — it is
`resolve_public_driver_context(None)`, which returns the cardiac context. Core
has roughly 25 call sites of the form
`driver_context = resolve_public_driver_context(driver_context)`, in
`introspection.py`, `strict_planning.py`, `specs/apply_overrides.py`,
`specs/dict_builder.py`, `specs/validation.py`, `runtime/registry.py`,
`runtime/sweep_runner.py`, `runtime/artifacts.py`, `runtime/run_document_exec.py`.
Each silently becomes cardiacFoam when a context is not threaded through.

Structural census of `core/compatibility.py`: 28 functions, **21 import cardiac
code** — 19 gated on `plugin_id == "org.cardiacfoam"`, 2 ungated
(`legacy_default_driver_context`, `legacy_generic_case_mutation`, both known and
deliberate per the predecessor spec).

### 1.7 Two hard dependencies and two optional ones are never imported

AST scan of every `import` across the package:

| declared | where declared | actually imported |
|---|---|---|
| `foamlib>=1.7.5,<2` | hard | core (4 files), plugins (5) |
| `jsonschema>=4.0` | hard | core (3), tests (1) |
| `PyYAML>=6.0` | hard | core (1), plugins (1), tests (1) |
| `gmsh>=4.15.2` | hard | **never** |
| `numpy>=1.24` | hard | **never** |
| `matplotlib`, `pandas`, `plotly` | `post` extra | postprocessing |
| `prompt-toolkit>=3.0` | `post` extra | **never** |
| `openpyxl>=3.1` | `post` extra | **never** |

`gmsh` is a false positive of the worst kind: it appears throughout the codebase
as a *workflow command string* — an external binary invoked by the DAG runner,
exactly like `blockMesh` or `checkMesh`. `core/specs/tet_mesh_provisioning.py:68`
says so explicitly: *"Never invokes gmsh -- this is a pure file render."* The
Python `gmsh` distribution is a ~100MB wheel installed to satisfy a binary that
must be on `PATH` regardless.

### 1.8 Root cause

Every defect above has one shape: **something holds a private copy of, or a
private default for, what the plugin contract already says.** The guide, the
scaffold, and the exporters each copied the contract. Core copied the *choice of
plugin*. Both rot the same way, and adding more hand-maintained copies would rot
identically.

### 1.9 What is not wrong

The contract itself is sound: `SolverPlugin` (27 required members) plus
`SolverPluginOptionalHooks` (14 hooks, 14/14 with prose docstrings), both
introspectable at runtime. The predecessor spec gated the capability fallbacks.
Cardiac code is not misplaced in core — it all lives in `plugins/cardiacfoam/`.
The problem is core's routing defaults and the absence of author-facing tooling,
not the contract and not the code's location.

---

## 2. Approach

Two principles, applied in order.

**Core never guesses which plugin it is driving.** A `DriverContext` is passed
explicitly; the cardiac default survives only at the public edge, where "no
plugin supplied" legitimately means "the built-in one."

**The contract is the single source for author-facing artifacts.** Everything a
plugin author reads is generated from it or diff-checked against it by a test.

### 2.1 Phases

Ordered so each lands on a green suite and the risky work precedes the new
surface.

| phase | content | risk |
|---|---|---|
| 0 | dead code and unused dependencies (§3) | none — nothing imports any of it |
| 1 | explicit `DriverContext` in core (§4) | **highest** — ~25 call sites |
| 2 | plugin-declared phases (§5) | low |
| 3 | cardiac reaches 14/14; retire gated fallbacks (§6) | low |
| 4 | the `authoring/` package (§7–§9) | low — new surface |
| 5 | removals, moves, anti-staleness (§10–§11) | low |

Phase 1 is separable and valuable alone. If it proves larger than estimated, it
may ship on its own and phases 2–5 follow; the reverse order is not available,
because §5 and §6 both assume an explicit context.

---

## 3. Phase 0 — dead code and unused dependencies

No behaviour change. Nothing here is imported by anything.

**`core/capability_manifest.resolve_case_models`** — delete. A shim documented
as deprecated (*"Kept for callers that imported this function directly before it
became plugin-owned"*) that calls `default_driver_context()` directly, a hard
cardiac default inside core. Verified: **zero production callers**; the only two
importers are `tests/core/test_case_introspection_capability.py` and
`tests/plugins/cardiacfoam/test_capability_manifest.py`, both of which test the
shim itself and are deleted with it.

**`pyproject.toml` dependencies** — hard deps 5 → 3, `post` extra 5 → 3:

- drop `gmsh>=4.15.2` (§1.7 — a workflow binary, not a Python import)
- drop `numpy>=1.24`
- drop `prompt-toolkit>=3.0`, `openpyxl>=3.1` from `post`

`gmsh`, `gmshToFoam`, and `checkMesh` remain in `CORE_NEUTRAL_COMMANDS`
(`core/runtime/workflow.py:56`) — they are authorized commands, which is
correct and unaffected.

Acceptance: full suite green, and a clean install in a fresh virtualenv without
the four packages still passes.

---

## 4. Phase 1 — core never guesses its plugin

The largest change and the one that actually decouples core.

Remove `resolve_public_driver_context(driver_context)` from core call sites and
make `driver_context` a required keyword parameter. The affected modules, from
§1.6: `introspection.py`, `strict_planning.py`, `specs/apply_overrides.py`,
`specs/dict_builder.py`, `specs/validation.py`, `runtime/registry.py`,
`runtime/sweep_runner.py`, `runtime/artifacts.py`, `runtime/run_document_exec.py`.

The cardiac default is retained at exactly four public-edge locations, each of
which is a documented compatibility boundary rather than core internals:

| edge | why the default is correct there |
|---|---|
| `openfoam_driver/dict_entries.py` | public module with external callers; already uses PEP 562 lazy re-exports for the same reason |
| `openfoam_driver/sweep_materialize.py` | public helper |
| `openfoam_driver/sweep_routing.py` | public helper |
| `cli.py:847` | "no `--plugin` given" legitimately means the built-in plugin |

`legacy_default_driver_context` and `resolve_public_driver_context` both survive,
called from four named places instead of everywhere.

### 4.1 Method

This is a mechanical refactor with a non-mechanical failure mode: a missed site
does not raise, it silently keeps working *and stays cardiac*. So the guard
comes first.

1. Add a test asserting `legacy_default_driver_context` is never invoked from
   within `core/` — using the existing `track_fallback_calls()` instrumentation,
   which already supports exactly this assertion (its docstring calls it "the
   P2.4 assertion an explicit non-cardiac v2 context should satisfy").
2. Convert call sites, innermost first, letting the type checker and the failing
   guard drive the order.
3. Re-run the §1.6 census. Acceptance is `legacy_default_driver_context`
   attributable only to the four edges.

### 4.2 Explicit non-goal

`--plugin` does not become mandatory on the CLI. `driverFoam plan --entry
singleCell` must keep working unchanged. The default moves; it does not vanish.

---

## 5. Phase 2 — plugin-declared phases

Fixes §1.5, the only silent-wrong defect.

Add `get_phases() -> tuple[str, ...]` to `SolverPluginOptionalHooks` — an
ordered tuple, since the order *is* the semantics (`primary_phase()` returns the
first declared phase an entry claims; the rest are read-only mirrors).

- `_PHASE_ORDER` derives from the active plugin's `get_phases()` instead of
  `get_args(Phase)`.
- Absent hook → `legacy_phases`, gated on `plugin_id == "org.cardiacfoam"`
  returning `("anatomy","physics","stimulus","solver")`, and for every other
  plugin the phases actually present in its `DictEntry` catalog. This follows
  the gating pattern the predecessor spec established.
- `primary_phase()` returning `None` becomes a reported condition, not a silent
  `continue`: the three `or "physics"` sites (§1.5) instead use the plugin's
  first declared phase, and the two `continue` sites emit a diagnostic.
- `Phase` stays as a type alias for documentation, no longer the runtime source
  of order.

Cardiac behaviour must be bit-identical: `CardiacFoamPlugin.get_phases()`
returns the same four in the same order, so `_PHASE_ORDER` is unchanged for it.
Pinned by the existing cardiac validation and dict-builder suites.

---

## 6. Phase 3 — cardiac reaches 14/14, gated fallbacks retire

Measured hook coverage (2026-08-26):

| plugin | optional hooks | absent |
|---|---|---|
| cardiacFoam | 12 / 14 | `get_report_catalog`, `get_config_resolution_description` |
| generic | 3 / 14 | 11 |
| minimal (test fixture) | 0 / 14 | 14 |

`CardiacFoamPlugin` gains both, plus `get_phases()` from §5. Then the
`plugin_id == "org.cardiacfoam"` branches are removed from `legacy_report_catalog`
and `legacy_config_resolution_description`.

**Net effect: `core/compatibility.py` stops importing `plugins/cardiacfoam/reports.py`
entirely** — retiring the exact indirection the dev-tools inventory needed a
page-long case study to explain, and the sole reason `reports.py` looks orphaned
to any static trace.

The remaining gated fallbacks are then dead for cardiac (which implements
everything) and neutral for everyone else. They exist for "v1 plugins", a
population the predecessor spec §1.3 established as zero. Each is deleted only
when the §1.6 census shows it at zero invocations outside tests that exercise it
deliberately; any that still fire are retained with a comment naming the caller.
Deleting on census evidence, not on assumption, is the rule here.

**Behaviour pinning.** §10 merges `test_report_catalog_export.py` away, so it
cannot simply "pass unmodified". Instead: before this phase lands, record the
current cardiac `reports.json` as a committed fixture; afterwards the merged
suite asserts byte-identical output. That pins the result across both the merge
and the fallback removal — the pair of changes that could silently drop a report.

---

## 7. Phase 4 — the `authoring/` package

*§8 and §9 detail this phase's three verbs.*

```
openfoam_driver/authoring/
  __init__.py         public API
  contract.py         introspect SolverPlugin + SolverPluginOptionalHooks -> model
  scaffold.py         model -> a plugin package on disk
  conformance.py      DriverContext -> tuple[StrictDiagnostic, ...]
  catalogs.py         DriverContext -> the four machine-readable catalogs
  pytest_plugin.py    third-party CI entry point
  templates/          .py.in / .yaml.in / .toml.in fragments
```

A sibling to `core/`, `plugins/`, `postprocessing/`, following the precedent
`postprocessing/` sets: maintained code whose consumer is not the `plan`/`run`
path. `core/` runs cases; `authoring/` exists before a plugin does.
`authoring/` imports `core/`; `core/` never imports `authoring/`.

`contract.py` is the keystone: it reads `_REQUIRED_PLUGIN_MEMBERS`,
`SUPPORTED_PLUGIN_API_VERSIONS`, `_PLUGIN_ID_RE`, and
`SolverPluginOptionalHooks`, recovering each member's signature and docstring by
`inspect`. Both `scaffold.py` and the generated `AUTHORING.md` render from this
one model, so neither can describe a contract the code does not have.

### 7.1 Prerequisite: five quote marks

`plugin_interface.py` quotes TYPE_CHECKING annotations inconsistently —
`get_extra_provenance_paths(self, case_root: "Path")` quotes,
`predict_data_artifacts(self, case_root: Path, ...)` does not. Under Python
3.14's deferred annotations these import fine but fail `inspect.signature()`.
Measured: 18/27 resolve, 4 fail only as properties (handled via `.fget`), and 5
genuinely fail — `get_dict_entries`, `get_dict_groups`, `get_tutorial_displays`
(`DictEntry`, `TutorialDisplay`), `validate_configuration` (`TutorialSpec`),
`predict_data_artifacts` (`Path`). Quoting those five makes the contract fully
introspectable and the file internally consistent. Annotation-only; no runtime
effect.

### 7.2 Verbs

`cli.py` uses a flat positional `action` with `choices=[...]`, and
`sweep-plan`/`sweep-run` already establish the hyphenated-namespace idiom. Three
choices are added; no parser restructuring:

```
driverFoam plugin-new    --id ID --name NAME --solver BINARY --out DIR
driverFoam plugin-check  --plugin TARGET [--strict] [--json]
driverFoam plugin-export --plugin TARGET --out DIR
```

Deliberately absent: `plugin-run` / `plugin-plan`. Those verbs exist and already
take `--plugin`. The toolkit stops where the plugin becomes drivable.

---

## 8. Phase 4, cont. — `plugin-new`

```
<out>/
  pyproject.toml                 entry-point stanza + package-data
  src/<pkg>/plugin.py            27 required members, stubbed
  src/<pkg>/plugin.yaml          from generic-plugin.yaml's annotated template
  src/<pkg>/dictionaries.py      DictEntry examples
  tests/test_conformance.py      imports authoring.pytest_plugin
  README.md                      the three commands, in order
```

**Every stub carries its contract docstring, copied by `inspect` at generation
time.** Nothing is transcribed, so no generated file can outlive a contract
change the way `SKILL.md` did.

**All 14 optional hooks are emitted commented-out**, each with its docstring and
documented absent-behaviour verbatim — `"Absent -> sweeps are refused by name"`,
`"Absent -> () ... every unknown file is a required input"`. This is the direct
fix for predecessor §1.0: *"a plugin author reading the public protocol never
learns the hook exists."*

Naming is derived, not asked twice: `--id org.myproject.shallowwater` yields
distribution `driverfoam-shallowwater` and package `driverfoam_shallowwater`,
from the id's last dot-segment. `--name` is the display string; `--solver` seeds
`get_solver_commands()` and `get_solve_step_commands()`. An id failing
`_PLUGIN_ID_RE` is refused before anything is written.

`--out` defaults to the current directory and is **refused if it resolves inside
the driverFOAM package tree**, so the tool cannot produce the in-tree plugin that
would invalidate a future portability claim.

---

## 9. Phase 4, cont. — `plugin-check` and `plugin-export`

### 9.1 `plugin-check`

Emits `StrictDiagnostic` records — the same `level`/`code`/`message` shape
`plan --strict` produces, so agents and existing tooling parse it without a new
vocabulary. Non-zero exit on any `error`; `--strict` promotes `warning`.

| # | family | level | checks |
|---|---|---|---|
| 1 | contract completeness | error | 27 required members present and callable; identity strings non-empty; `plugin_id` matches `_PLUGIN_ID_RE`; api version supported; `plugin.yaml` id/api match the class; `driver_path` uniqueness |
| 2 | fallback reliance | warning | per optional hook: implemented, or resolving to a named `legacy_*` fallback, and what that returns |
| 3 | cardiac reachability | error | drive every capability with this plugin; assert `plugins.cardiacfoam` never imported |
| 4 | vocabulary leakage | error | `DictEntry.phases` outside the plugin's own `get_phases()`; `build_run_document_config()` returning cardiac phase keys the plugin never declared |
| 5 | catalog well-formedness | warning | empty `description`; `enum_values` on a non-enum `value_kind`; `required_when`/`forbidden_when`/`mutually_exclusive_with` naming absent `driver_path`s |
| 6 | implicit-default reliance | warning | run the plugin's operations under `track_fallback_calls()`; report any `legacy_default_driver_context` hit as a place core fell back to the built-in plugin |

Family 1 is `validate_plugin()` + `driver_context()` reported as a complete
diagnostic list rather than a `TypeError` on first failure. Family 3 generalises
`tests/core/test_capability_fallback_neutrality.py` — today asserted only for
`GenericOpenFOAMPlugin` — and ships it to third parties; it is the check that
would have caught the wrong-`Allrun` defect. Family 6 is the residue of §4: after
Phase 1 it should be empty for any well-formed plugin, and it stays as the
regression guard.

### 9.2 `plugin-export`

```
<out>/{dict-catalog,tutorials,utilities,reports}.json + manifest.json
```

Everything resolves through `driver_context.capabilities.*`; the exporter imports
no plugin package. The `PHASES` hardcode is replaced by the plugin's
`get_phases()` (§5). Fan-out semantics are preserved from the current exporter:
one record per `(entry, phase)` pair, each stamped with a single `phase` and
retaining the full sorted `phases` list.

`tutorials.json` keeps the current cross-check — the plugin's
`get_tutorial_displays()` ids must equal `registry.list_tutorials()` — but
reports a mismatch as a diagnostic rather than `SystemExit`.

---

## 10. Phase 5 — removals, moves, and the remaining core fixes

### 10.1 `core/utility_catalog.py`

```python
# line 542, today
UTILITY_CATALOG: Final[dict[str, UtilityManifest]] = load_utility_manifests(
    Path(__file__).resolve().parents[4] / "utilities"
)
```

A core module, eagerly evaluated at import, hardcoded to this repository's
layout. The constant is **deleted outright**, leaving `load_utility_manifests(roots)`
as the only API, called with the active plugin's `get_utility_roots()`.

No compatibility shim, and none should be added: verified, `UTILITY_CATALOG` has
no production consumer — its only non-test importer is
`scripts/export-utility-catalog.py`, which §10.3 deletes. The contract seam
already produces the identical result:

```
CardiacFoamPlugin().get_utility_roots() -> ('<repo>/applications/utilities',)
load_utility_manifests over those roots -> 12 utilities
UTILITY_CATALOG (hardcoded)             -> 12 utilities
identical: True
```

`tests/core/test_mesh_geometry.py` and `tests/core/test_utility_catalog_export.py`
are updated to call `load_utility_manifests()` with the cardiac roots explicitly.

Retaining a lazy `__getattr__` — the precedent `dict_entries.py` sets for
`CONTROL_DICT_ENTRIES` — was considered and rejected. That precedent protects
*external* importers of a public name; this name has none, and keeping it would
hide the repo-layout hardcode rather than remove it.

### 10.2 `cli.py build_parser()`

`--entry`'s help stops interpolating `list_tutorials()` at construction time
(§1.4). The plugin-specific entry list moves to where `--plugin` is known:
`registry.resolve_entry`'s `KeyError` already names valid entries, and `describe`
reports them.

### 10.3 Removals and moves

| path | action | rationale |
|---|---|---|
| `scripts/export-dict-catalog.py` | delete | → `plugin-export`; cardiac imports and `PHASES` hardcode die with it |
| `scripts/export-tutorials-catalog.py` | delete | → `plugin-export`; drops the `default_driver_context()` anti-pattern |
| `scripts/export-utility-catalog.py` | delete | → `plugin-export`, from `get_utility_manifests()` |
| `scripts/export-report-catalog.py` | delete | → `plugin-export`; already contract-correct, merely one of four |
| `scripts/regenerate-ionic-catalog.py` | move → `plugins/cardiacfoam/` | cardiac tooling belongs in the cardiac plugin — the boundary rule applied to the tools |
| `openfoam_driver/scripts/_names_parser.py` | move with it | its only consumer |
| `openfoam_driver/scripts/_rtst_scanner.py` | fold into `tests/drift_guards/test_rtst_enum_contract.py` | 188 lines, one consumer, and that consumer is a test |
| the four export conformance tests | merge → one parametrized suite | four near-identical subprocess-and-assert-shape tests |
| `.coverage` (untracked) | delete | regenerable |
| `applications/scripts/driverFoam/ROADMAP.md` (untracked) | delete | a three-line coverage note, not a roadmap; name collides with the tracked `future/driverFOAM/ROADMAP.md` |
| `devTools_inventort.md` (untracked) | delete | misspelled, untracked, a hand-maintained snapshot of what `plugin-check` reports on demand |

Retained deliberately: `scripts/export-capability-seams.py` and
`schemas/generate_run_document_schema.py` (genuinely core-internal generators),
and `scripts/scan-dict-keys.py` (contract-correct, already accepts `--plugin`,
and its C++ dict-key audit is a different job from plugin conformance).

After this, `scripts/` holds exactly two tracked files.

---

## 11. Anti-staleness

The disease is staleness, so these are load-bearing.

1. **`AUTHORING.md` is generated from `contract.py`**, never hand-written, with a
   `--check` mode and a conformance test — the pattern
   `export-capability-seams.py` already uses for `ARCHITECTURE.md`. A
   v1/v2-style divergence becomes structurally impossible.
2. **Scaffold round-trip in CI**: generate into a temp dir, load, `plugin-check
   --strict`, assert clean. A 28th required member either appears in the
   generated scaffold or breaks this test.
3. **The guide enters the repository.** `.gitignore:121` gains an exception for
   `.agents/skills/driverfoam-plugin-builder/`, matching how
   `driverfoam-assistant/SKILL.md` was force-added. The skill is rewritten as a
   thin pointer to the generated `AUTHORING.md` plus the three commands — it must
   not restate the contract, because a second copy is what rotted.
4. **The boundary guard widens.** `tests/core/test_plugin_dependency_boundary.py`
   scans only `core/`; `tests/core/test_no_top_level_cardiac_imports.py` checks
   exactly one file despite a docstring claiming a general rule. Both widen to
   cover `scripts/` and `authoring/`. `scripts/` sitting outside
   `openfoam_driver/` is why `export-dict-catalog.py` could import cardiac freely.
5. **A dependency guard.** A test asserting every distribution in
   `[project].dependencies` is imported somewhere in the package — so a future
   unused dependency fails CI instead of accumulating (§1.7).
6. **Stale references fixed at source:** `generic-plugin.yaml`'s
   `api_version: "1" or "2"`, `plugin_interface.py:287`'s `_REQUIRED_V2_MEMBERS`
   mention, and `generic_plugin.get_tutorial_catalog()`'s missing
   `make_generic_case_spec` key.

---

## 12. Out of scope

**The out-of-tree reference plugin.** `SOLIDS4FOAM_CASE_STUDY.md` stage S1
remains deferred. This work is its prerequisite, not its substitute.

**Any change to the publication claim.** `ARCHITECTURE.md:18` concedes a
portability claim needs an out-of-tree plugin and end-to-end CI. Neither is
delivered here, so the claim stays "decoupling evidence" — though §4 and §6
materially strengthen the evidence behind it.

**Making `--plugin` mandatory** (§4.2). **`cli.py` decomposition** — 963 lines and
on the tech-debt list, but three new `action` choices do not justify restructuring
it here. **Converting `cli.py` to argparse subparsers.**

---

## 13. Testing

New `openfoam_driver/tests/authoring/`:

- `test_contract_model.py` — the introspected model matches
  `_REQUIRED_PLUGIN_MEMBERS` and `SolverPluginOptionalHooks` exactly; every
  member yields a signature and a docstring
- `test_scaffold_round_trip.py` — generate, load, `plugin-check --strict` clean;
  and generation refuses an `--out` inside the package tree
- `test_conformance_checks.py` — one deliberately-broken fixture plugin per check
  family, each proving the check fires. Mutation-verified, following the
  predecessor spec's method
- `test_catalog_export.py` — parametrized over `{cardiac, generic, minimal}`,
  replacing the four merged export tests
- `test_authoring_doc_current.py` — `AUTHORING.md --check`

New in `tests/core/`: the Phase 1 guard (§4.1), the dependency guard (§11.5), and
the widened boundary scans (§11.4).

**Baseline, measured 2026-08-26**, unsourced OpenFOAM environment,
`.venv/bin/python -m pytest openfoam_driver/tests -q`:

```
1690 passed, 9 skipped, 130 subtests passed in 207.54s   (exit 0)
```

Zero failures. Two notes for the implementer. The predecessor spec closed at 1574
passed with 3 known failures; the suite has since grown by 116 tests and those 3
are gone — do not carry forward a "3 pre-existing failures" expectation, since
any failure during this work belongs to this work. And a characterization suite
absorbs regressions without the failure *count* moving (a lesson already paid for
in this repo), so compare failure content, not totals.

Each phase lands on a green suite before the next begins.

---

## 14. Success criteria

1. No module under `core/` resolves an implicit cardiac context. Proven two
   ways, both mutation-verified: a static AST guard that `core/` never calls
   `resolve_public_driver_context`, and a runtime guard that a plan driven by a
   non-cardiac plugin reaches `legacy_default_driver_context` zero times.

   **Corrected 2026-08-26 (measured).** This criterion originally read *"proven
   by re-running the §1.6 census"*. That is not a valid proof and the census
   must not be used as one: the public `default_driver_context()` routes
   straight through `legacy_default_driver_context()`, so the whole-suite
   counter cannot distinguish the defect (an implicit fallback) from correct
   usage (an explicit call). Phase 1 eliminated the implicit calls in core while
   threading ~60 explicit ones into tests, so the census moved only
   51,540 -> 51,358 despite the defect being gone. The §1.6 census remains a
   useful *discovery* tool for finding where fallbacks fire; it is not an
   acceptance metric.
2. A plugin declaring its own phase vocabulary receives full required-field and
   enum validation, and never has entries written into a `"physics"` slice it did
   not declare.
3. `core/compatibility.py` contains no import of `plugins/cardiacfoam/reports.py`.
4. `plugin-check` reports cardiacFoam at 14/14 optional hooks with zero
   error-level diagnostics.
5. A developer runs `plugin-new`, `plugin-check`, `plugin-export` and holds a
   loading plugin plus its four machine-readable catalogs, without editing
   driverFOAM and without reading `core/`.
6. `plugin-export` produces correct catalogs for `{cardiac, generic, minimal}`
   and for a freshly scaffolded plugin.
7. `driverFoam --plugin <non-cardiac> --help` names no cardiacFoam tutorial.
8. `AUTHORING.md --check` passes, and fails when a contract member is added
   without regeneration.
9. Every declared runtime dependency is imported somewhere in the package,
   asserted by test.
10. `scripts/` contains exactly two tracked files, both core-internal generators.
11. No new hand-maintained description of the plugin contract exists anywhere in
    the tree.
