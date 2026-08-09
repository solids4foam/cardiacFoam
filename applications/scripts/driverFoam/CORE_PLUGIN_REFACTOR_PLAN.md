# driverFOAM Core / Plugin Refactor Plan

## Decision

Make the core orchestrator independent of cardiacFoam **and make plugin
selection explicit, immutable for a run, and recorded in every strict
contract**.  A plugin may describe solver meaning; it may not replace the
workflow security boundary or silently alter an already planned run.

The compatibility target is preserved:

```text
foamctl plan --strict --entry singleCell
foamctl run  --strict --entry singleCell
```

continues to select the built-in cardiacFoam plugin when the CLI is launched
without a plugin selection.  The implicit selection is a CLI bootstrap
compatibility shim only; it must not remain a mutable module-global dependency
inside the core.

## Audit of the starting point

The existing direction is sound, but it is not yet sufficient to prove solver
agnosticism or safe injection.

| Finding | Evidence | Required correction |
|---|---|---|
| Plugin state is process-global. | `core/plugin_interface.py` stores `_ACTIVE_PLUGIN`; `set_active_plugin()` changes it and `get_active_plugin()` lazily imports cardiacFoam. | Resolve a plugin once into a per-operation `DriverContext`; pass that context explicitly. No mutable global plugin state. |
| The loader is arbitrary Python import. | `cli.py --plugin module:Class` imports and instantiates an unchecked class. | Keep this only as an explicitly documented trusted-development escape hatch, or replace it with installed entry points. Validate the object and record its identity. It is not a sandbox. |
| “Generic” case execution is cardiacFoam-owned and assumes cardiac files. | `plugins/cardiacfoam/tutorials/generic_case.py`; `registry.py` detects `electroProperties` and requires `physicsProperties`. | Move a minimal case-folder spec to core; recognize structural OpenFOAM cases without cardiac dictionaries. |
| Validation remains cardiac-specific before plugin validation runs. | `specs/validation.py`, `validation_rules.py`, and `dict_entries.py` import cardiacFoam catalogs/rules. | Split generic RunDocument validation from plugin semantic validation. The core must import no cardiac modules. |
| Strict planning cannot pass for a neutral plugin. | `strict_planning._artifact_diagnostics()` treats an empty prediction as an error. | Core predicts run-state/log artifacts and lets domain predictions be empty. |
| Python/C++ drift checking assumes one monorepo C++ tree and one global catalog. | `strict_planning._catalog_diagnostics()` scans `<repo>/src`; `_dict_keys_scanner.py` imports the active catalog but uses one global allowlist. | Each plugin declares catalog provenance, C++ source roots, scanner settings, and its reviewed allowlist. |
| Builder and sweep routing encode cardiac selectors and materialization. | `dict_builder.py` invokes a plugin tutorial factory but hard-codes `cardiacFoam`; `sweep_routing.py` imports cardiac-shaped catalogs. | Make materialization and routing a plugin capability, with a small core generic-mutation alternative. |
| Current protocol proof is nominal only. | The mock in `test_plugin_architecture.py` does not implement all protocol members; `Protocol` is not runtime-checked. | Add runtime contract validation and a genuine minimal-plugin integration suite. |

These issues mean that merely moving files would preserve hidden cardiacFoam
dependencies, test-order dependence, and non-reproducible execution.

## Boundary and trust model

### Core owns

- CLI parsing and plugin selection policy
- entry/path resolution, generic case-folder detection, and RunDocument schema
- workflow-DAG normalization, command authorization, argv execution, retries,
  state, resume, freshness, and output containment
- generic OpenFOAM structure checks and core-owned bookkeeping artifacts
- generic text/dictionary mutation primitives only
- the diagnostic, artifact, and manifest schemas

### Plugin owns

- named tutorial registry and display metadata
- solver-specific configuration intent, selector vocabulary, catalog entries,
  defaults, and materialization
- solver/domain semantic validation
- domain artifact prediction and field capability descriptions
- C++/Python catalog provenance and drift scanning configuration

### Non-negotiable security rules

1. A plugin is trusted in-process code.  Loading one is equivalent to running
   arbitrary Python; it is not an agent-safety boundary.
2. Plugin hooks return data and diagnostics.  They do not receive a shell,
   subprocess callback, or workflow executor.
3. Only core normalizes and authorizes workflow commands.  A plugin can
   propose a DAG, but cannot expand the command allowlist or bypass path,
   argv, output-root, retry, or `--fresh` checks.
4. `Allrun` remains deliberately case-authored code under the existing
   local/single-tenant trust model.  A generic case folder is therefore not
   safe to execute merely because it is discoverable.
5. An agent must expose the selected plugin id/version and its source before
   planning or running.  It must never infer a plugin from untrusted case
   content.

## Target contracts

Use small, typed, data-oriented contracts.  Avoid one large `SolverPlugin`
interface with nullable methods and untyped dictionaries.

```python
@dataclass(frozen=True)
class PluginIdentity:
    id: str                 # reverse-DNS or distribution-qualified, not display name
    version: str
    api_version: str
    source: str             # installed entry point or explicit trusted import
    capability_digest: str  # canonical JSON SHA-256 of planning-relevant data

@dataclass(frozen=True)
class DriverContext:
    plugin: SolverPlugin | None
    identity: PluginIdentity | None
    mode: Literal["generic", "plugin"]

class SolverPlugin(Protocol):
    def describe(self) -> PluginDescriptor: ...
    def tutorials(self) -> TutorialRegistry: ...
    def catalog(self) -> DictionaryCatalog | None: ...
    def semantic_validators(self) -> tuple[SemanticValidator, ...]: ...
    def materializer(self) -> CaseMaterializer | None: ...
    def predict_artifacts(self, request: ArtifactRequest) -> tuple[DataArtifact, ...]: ...
    def cxx_mapping(self) -> CxxMapping | None: ...
```

The exact names can vary, but these properties are required:

- Every returned value is immutable or treated as immutable by core.
- `PluginDescriptor`, `TutorialRegistry`, `DictionaryCatalog`, and
  `CxxMapping` have versioned JSON forms and deterministic ordering.
- Diagnostics always include `origin: "core" | "plugin:<id>"` and a stable
  code.  Plugin exceptions become a clear error diagnostic; they must never
  produce a traceback-shaped successful plan.
- `None` means the capability is unavailable, not “use cardiacFoam”.
- Plugins receive a read-only request object with resolved paths, config,
  normalized DAG, and identity; they do not mutate `TutorialSpec` or process
  globals.

### Declarative OpenFOAM case profiles: use YAML as data, not as the plugin API

YAML is a good seam for **declarative case evidence and dictionary targets**;
it is the wrong seam for arbitrary solver semantics or executable behaviour.
Use a schema-validated, versioned YAML document (JSON is accepted because it
is a YAML subset) that is loaded with `safe_load`, canonicalized to JSON for
the capability digest, and converted immediately to typed Python objects.
Do not permit YAML tags, anchors/merges, templates, expressions, imports, or
embedded Python/shell.

Split declarative data into two scopes:

1. `driverfoam-case.yaml` is optional, case-local evidence: its workflow
   intent, explicit required/generated files, and explicitly permitted generic
   dictionary mutations.  It replaces the narrow parts of the current
   `workflow_contract.json`, but JSON remains supported during migration.
2. A package-owned `plugin.yaml` describes the solver profile: its catalog,
   dictionary targets, standard OpenFOAM conventions, C++ mapping roots, and
   static artifacts.  It is part of the plugin distribution, never inferred
   from untrusted case files.

For example, a profile may declare the *shape* of a cardiac case without
placing any cardiac name in core:

```yaml
schema_version: 1
plugin: {id: org.cardiacfoam, api_version: "1"}
case_profile:
  dictionaries:
    - id: run-control
      path: system/controlDict
      kind: openfoam_dictionary
      role: openfoam.control_dict
      required: always
    - id: electro-properties
      path: constant/electroProperties
      kind: openfoam_dictionary
      role: plugin.configuration
      required: always
  generated_files:
    - path: constant/polyMesh/points
      produced_by: mesh
catalog:
  entries_file: catalog/electro-properties.yaml
cxx_mapping:
  source_roots: [../../../src]
  reviewed_allowlist: catalog/dict-key-allowlist.json
```

The core interprets only `kind`, path containment, `required`,
`generated_files`, and the generic `openfoam_dictionary` mutation grammar.
`role: openfoam.control_dict` is useful metadata for generic reports and
mutation helpers, but it is not a claim that every OpenFOAM solver requires a
`controlDict`.  A plugin—or a case contract—declares whether it is needed.
`plugin.configuration` is opaque to core.

Keep conditional requirements deliberately small (`always`, `never`, or a
structured data predicate such as `all: [{path: solver, equals: piso}]`).
If a condition needs C++ model knowledge, nested dictionary traversal, ranges,
or an explanation tailored to users, it belongs in the plugin semantic
validator.  In particular, do not encode ionic-model compatibility,
heterogeneity ranges, solver/coupler relationships, or mesh provisioning
algorithms in YAML.

This gives a useful three-level seam:

```text
core OpenFOAM mechanics  →  declarative profile/case evidence  →  plugin Python semantics
paths, dictionary I/O,       files, paths, simple predicates,      model compatibility,
DAG/state/security           static catalog metadata               synthesis, prediction
```

This is deliberately OpenFOAM-aware, not solver-agnostic in the abstract.
The core may own OpenFOAM path/dictionary primitives (`foamDictionary` when
available, read/update, `system/controlDict` as a named convention, and
standard `Allrun` handling).  It must not turn those conventions into
universal requirements or assume that one dictionary identifies every solver.

### Plugin resolution

Implement `resolve_driver_context(selection)`, called once by each public
operation (`plan`, `run`, `step`, sweep, introspection, and Python API).  The
CLI default resolves the built-in cardiac plugin; the Python API defaults to
`generic` unless the caller supplies a context.  This removes surprising
cardiac defaults from reusable library calls.

Preferred selection is an installed package entry point such as
`driverfoam.plugins`, addressed by stable id (`--plugin cardiacfoam`).  During
migration, retain `--plugin module:Class` only with an opt-in flag such as
`--unsafe-plugin-import`, label it as trusted local development, and include
the literal import target in `PluginIdentity.source`.

`RunDocument` v3 gains:

```json
"plugin": {
  "id": "org.cardiacfoam",
  "version": "…",
  "api_version": "1",
  "capability_digest": "sha256:…"
}
```

`run --run-document` resolves that exact plugin id and rejects a missing or
incompatible id/API/digest before executing.  A controlled `--accept-plugin-
drift` may re-plan and write a new document; it must not silently execute an
old plan under new semantics.  V1/V2 migration injects the selected default
identity and emits a migration warning.  Do not pin Python wheel hashes in
this phase; distribution provenance plus capability digest is the useful
semantic reproducibility boundary.

## Migration plan

### Phase 0 — Freeze behavior and establish seams

Before moving production code, add characterization tests for:

- every current named cardiacFoam entry’s normalized DAG, diagnostics shape,
  and artifact ids
- `--plugin` selection, malformed plugins, and no cross-test/plugin leakage
- direct case-folder planning with and without `electroProperties`,
  `physicsProperties`, and `Allrun`
- RunDocument v1/v2 compatibility and resume/fresh behavior
- existing command-boundary guarantees.

Record a fixture corpus of sanitized strict-plan JSON (assert stable fields,
not volatile absolute paths or timestamps).  This is the regression oracle
for all later phases.

**Exit gate:** the current test suite passes and tests can run a plugin A then
plugin B in the same interpreter without affecting each other.

### Phase 1 — Contextual plugin injection and runtime contract validation

Files:

- change `core/plugin_interface.py`; add `core/plugin_resolution.py` and
  `core/plugin_contracts.py`
- change `cli.py`, `strict_planning.py`, `introspection.py`, runtime registry,
  artifact prediction, builders, and sweep entry points to accept a context
- add a core test fixture that constructs a new context per test.

Actions:

1. Add the contracts above and `validate_plugin(plugin)`.  It checks required
   members, JSON serializability, stable id grammar, non-empty version/API,
   duplicate tutorial names, duplicate catalog paths, and artifact-id
   uniqueness.
2. Replace all calls to `get_active_plugin()` with an explicit context passed
   from the public boundary.  Delete `_ACTIVE_PLUGIN`, `set_active_plugin`,
   and the lazy cardiac import.
3. Make context explicit in `strict_plan`, `load_entry_spec`, introspection,
   dictionary building, sweep planning, and RunDocument execution.  Avoid a
   `contextvars` compromise: it still hides a dependency and makes agent/tool
   concurrency harder to reason about.
4. Keep a CLI-only `default_context()` that selects cardiacFoam when no plugin
   is specified.  Add `--plugin none` for truly generic invocation.

**Exit gate:** core package import does not import `plugins.cardiacfoam`; two
different contexts are safe in one process; every output report names the
resolved plugin identity.

### Phase 2 — Core generic case folders

Files:

- add `core/runtime/generic_case.py`
- change `core/runtime/registry.py`, `execution_context.py`, and strict audit
- reduce `plugins/cardiacfoam/tutorials/generic_case.py` to a compatibility
  adapter, then delete it after deprecation.

Define a generic case as a directory with either a valid `workflow_contract`
or an executable permitted case-script convention (`Allrun` / `Allrun.*`).
`system/controlDict`, `fvSchemes`, `fvSolution`, `constant/electroProperties`,
and `physicsProperties` become *diagnostic evidence*, never generic discovery
requirements.  Core chooses on-disk `workflow_contract` steps over the
`Allrun` fallback and makes missing/incomplete structure an explainable strict
diagnostic rather than routing through a cardiac factory.

Core generic support provides only:

- path containment and output location resolution
- a normalized `Allrun` fallback DAG
- the existing generic text/dictionary mutation hook for explicitly supplied
  paths
- core bookkeeping artifacts: `workflow_state.json`, workflow logs,
  `artifacts_manifest.json`, and terminal `artifacts_realized.json`.

It does not synthesize cardiac dictionaries, assume solver names, source a
cardiac run script, or invoke postprocessing modules.

**Exit gate:** `foamctl plan --strict --plugin none --entry <plain-case>` can
produce a valid plan for a minimal `Allrun` case and an explicit failure for a
non-runnable directory, without importing cardiacFoam.

### Phase 3 — Catalog, materialization, and C++↔Python provenance

Files:

- move `DictEntry`, `Phase`, and generic catalog types to
  `core/contracts/dictionary.py`
- move cardiac `PHYSICS_PROPERTY_ENTRIES`, `CONTROL_DICT_ENTRIES`, batched
  models, selector paths, and templates into `plugins/cardiacfoam/`
- change `specs/dict_builder.py`, `sweep_routing.py`, and scanner APIs to use
  `DictionaryCatalog` / `CaseMaterializer`
- change `_dict_keys_scanner.py` into a source-root-parameterized core utility.

`DictionaryCatalog` must carry, per entry:

- logical document and concrete OpenFOAM path template
- selector/override role, value type, phases, default, structured constraints
- C++ provenance: source root id, file, symbol or read expression, and line
  when known
- the mapping confidence (`exact`, `generated`, `reviewed_heuristic`) and a
  human rationale for non-exact mappings.

`CxxMapping` supplies one or more C++ roots and a plugin-owned allowlist.
`strict_plan` invokes scanning only when the active plugin supplies this
mapping.  It reports results under `plugin_cxx_catalog_diagnostics`, tagged
with the plugin id; it no longer scans `<repo>/src` unconditionally.  The
scanner remains a drift detector, not proof that a C++ key is semantically
accepted—dynamic names, macro expansion, aliases, and inheritance need an
explicit reviewed mapping.

The cardiac materializer owns cardiac selector routing, dictionary synthesis,
default mesh policy, templates, and the solver command.  The core can expose
generic safe `set_foam_entry` / copy / file-creation primitives, but must not
know `myocardiumSolver`, `ionicModel`, `dx`, or `cardiacFoam`.

**Exit gate:** `rg 'plugins\.cardiacfoam' openfoam_driver/core
openfoam_driver/specs openfoam_driver/sweep_*.py` has no production hits; a
cardiac catalog scan produces the same reviewed drift report as the baseline.

### Phase 4 — Validation pipeline

Files:

- add `core/runtime/validation_contracts.py`
- split `specs/validation.py` into core RunDocument/schema validation and
  `plugins/cardiacfoam/validation.py`
- move `validation_rules.py` and cardiac model coupling/heterogeneity rules
  under the plugin.

Pipeline order:

```text
resolve context → resolve entry → build/parse RunDocument
  → core structural validation → normalize + authorize workflow
  → plugin semantic validation → core artifact bookkeeping
  → plugin artifact enrichment → final plan
```

Core validates paths, schema, DAG shape, output paths, declared artifact
schema, and generic case structure.  The plugin validates semantic values
such as ionic/tissue compatibility, coupling topology, heterogeneity ranges,
and required cardiac dictionary blocks.  Pass the plugin a normalized
`ValidationRequest` (intent/config plus case evidence), so validation does not
quietly parse one file while materialization validates another representation.

Preserve diagnostic JSON field names where possible.  Add `origin` additively;
do not rename current cardiac diagnostic codes during the extraction.

**Exit gate:** a no-domain plugin receives no cardiac validation errors; all
current cardiac invalid fixtures retain their codes and fail status.

### Phase 5 — Artifact and capability composition

Files:

- change `core/runtime/artifacts.py`, `strict_planning.py`,
  `capability_manifest.py`, and function-object field diagnostics
- retain `plugins/cardiacfoam/artifacts_predictor.py` as the cardiac provider.

Core always declares bookkeeping artifacts.  Plugins may append domain
artifacts, but an empty plugin prediction is valid.  Merge artifacts through
a core function that enforces stable ids, namespace rules (`core.*` and
`plugin.<id>.*`), path containment, and deterministic sort order.  A plugin
cannot overwrite core artifacts.

Replace the cardiac-shaped global capability manifest with:

```json
{
  "core": {"allowed_commands": "…", "run_document_api": "…"},
  "plugin": {"identity": "…", "capabilities": "…"}
}
```

Function-object sampled-field checks run only when a plugin advertises field
capabilities.  Generic runs report them as not-applicable, not as a cardiac
warning.

**Exit gate:** the minimal plugin strict plan succeeds with only core
artifacts; cardiac artifact ids and manifest serialization remain stable.

### Phase 6 — Registries, sweeps, and public API migration

Files:

- change `registry.py`, `introspection.py`, `tutorials_display.py`,
  `dict_builder.py`, `sweep_routing.py`, `sweep_runner.py`, and CLI help
- add plugin entry-point metadata for built-in cardiacFoam.

Keep filesystem discovery in core and named tutorials in the plugin.  Make
introspection explicitly distinguish `generic_case_folder`,
`plugin_registered_tutorial`, and `workflow_contract`.

Split sweep modes:

- generic sweeps mutate only explicitly declared files/paths from a generic
  sweep schema; no solver selector vocabulary
- plugin materialized sweeps route axes through that plugin’s catalog and
  materializer, and pin plugin identity/digest in `sweep_manifest.json`.

Deprecate public helpers that rely on ambient plugin state.  Provide explicit
replacement signatures and a one-release warning shim only if external users
need it.

**Exit gate:** a cardiac sweep still materializes byte-equivalent dicts for
the baseline fixtures; a generic sweep needs no cardiac plugin.

### Phase 7 — Proof, hardening, and removal

Add `tests/plugins/minimal/` with a real plugin (not a mock) that has no
catalog, no semantic validators, no domain artifacts, and no named tutorials.
It must plan and run a tiny permitted `Allrun` fixture, emitting core state and
artifacts.  Add a second fixture plugin with a one-entry catalog to prove
catalog isolation.

Required tests:

- plugin A / plugin B / generic execution in one process, including parallel
  planning threads
- missing, malformed, duplicate-id, and throwing plugins
- plugin identity mismatch and capability-digest drift for RunDocuments and
  sweep manifests
- command allowlist rejection for a plugin-proposed arbitrary command
- no-plugin generic discovery with no cardiac dictionaries
- C++ scanner runs only the active plugin’s declared roots and allowlist
- full cardiac regression equivalence and its existing real-solver integration
  matrix.

Move tests according to ownership: true engine behavior under `tests/core`,
cardiac semantics under `tests/plugins/cardiacfoam`, and cross-plugin
contract/security tests under `tests/plugins` or `tests/core/plugin_*`.
Remove global-plugin compatibility code only after all gates pass.

## Acceptance criteria

The refactor is complete when:

- core and generic/spec/sweep production modules have no direct imports of
  `plugins.cardiacfoam`;
- no mutable active-plugin global remains;
- each plan, run, RunDocument, and sweep manifest identifies the exact plugin
  semantics used;
- a plain OpenFOAM case folder can be strictly planned and run without
  cardiac-specific files or a cardiac plugin;
- a minimal plugin can run without domain validation or domain artifacts;
- cardiacFoam retains its current CLI UX, strict diagnostics, materialized
  dictionary behavior, and artifact manifest format;
- C++↔Python mapping is plugin-scoped, reviewable, and no longer assumes the
  monorepo `src/` tree;
- plugins cannot bypass core workflow authorization; and
- the full regression suite and real-solver verification suite pass.

## Implementation order

Execute phases 0–2 first.  Do not move catalogs or validation until context
injection and generic execution are independently green; otherwise failures
will be impossible to attribute.  Then execute 3–5 as one vertical slice for
the cardiac plugin, followed by registry/sweep migration and the minimal
plugin proof.

The governing rule remains simple: core decides **whether and how a workflow
is safely executed**; a plugin decides **what a solver configuration means**.
