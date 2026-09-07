# driverFOAM plugin seam: capability gating and documentation

**Date:** 2026-08-19
**Status:** COMPLETE and verified 2026-08-19 (Part A, Part B, §4.5 amendment).

**Final suite:** 3 failed / 1567 passed, against a 7 failed / 1455 passed
baseline. The same 3 pre-existing, unrelated failures remain; +112 new tests
(4 fallback-neutrality, 108 seam-documentation conformance). Zero regressions.

**Docstring coverage:** `plugin_interface.py` 54/54. In
`plugin_capabilities.py` all 31 public classes are documented (0 undocumented).
A naive scanner still reports ~100 items there: 40 one-line `Protocol` method
stubs, covered by their class docstring, and 62 private `_*Adapter` members,
which are implementation and deliberately out of scope per §5. That is the
correct end state, not residual debt -- the original "96 undocumented items"
figure counted exactly these.

**Conformance tests are mutation-verified.** Un-gating a fallback, inventing an
`:adapts:` member, and replacing a docstring with `"""TODO."""` were each
introduced deliberately and each caught.

**Verification of Part A (full suite, unsourced OpenFOAM env):**

| | failed | passed |
|---|---:|---:|
| baseline before change | 7 | 1455 |
| after gating + test updates | **3** | **1459** |

The 3 remaining failures are pre-existing and unrelated, present on a clean
checkout: a stale `cardiac_tutorial_characterization.json` digest fixture, and
two `TestBuildAndLaunchDirectRun` tests whose `patch("subprocess.run")` is
over-broad and also intercepts the OpenFOAM bashrc probe. Of the original 7, 4
were this change's own acceptance tests, now green. Zero regressions.

Cardiac behaviour verified identical: `has_case_marker` True,
`is_runnable_without_workflow` True, `materialize` still writes
`'#!/bin/sh\nblockMesh\ncardiacFoam\n'`. The generic plugin now returns
False/False and refuses materialization by hook name.

Three existing tests were updated rather than relaxed, all of which pinned the
pre-gate behaviour:
`test_plugin_capabilities.py::test_legacy_plugin_case_evidence_preserves_pre_capability_behavior`
(inverted and renamed), the RunDocument config assertion in the same file
(cardiac phase vocabulary -> `{}`), and
`test_sweep_routing.py::test_routing_uses_the_selected_plugin_catalog`
(now asserts the stronger refusal). The last was NOT predicted by the
pre-change audit and was found only by running the suite.

**Scope:** `applications/scripts/driverFoam/openfoam_driver/core/`

---

## 1. Problem

A review flagged that driverFOAM's "plugin capabilities are not explained,"
citing 96 undocumented items in `core/plugin_capabilities.py`. Investigation
showed the count is largely an artefact of a docstring scanner treating every
one-line `Protocol` method stub as a separate item, and that the file's real
problem is not missing prose.

Measured docstring coverage:

| file | documented / defs |
|---|---|
| `core/plugin_capabilities.py` | 16 / 134 |
| `core/plugin_interface.py` | 25 / 39 |
| `core/plugin_discovery.py` | 4 / 4 |
| `core/plugin_profile.py` | 2 / 7 |
| `core/generic_plugin.py` | 12 / 31 |

Of the 21 capability `Protocol` classes, 10 carry substantial docstrings and 11
carry none. The split is not arbitrary. The 10 documented ones are the later,
optional hooks whose `compatibility.py` fallback checks
`plugin_id == "org.cardiacfoam"` and returns a neutral empty value for every
other plugin.

The 11 undocumented ones divide further, and the division matters more than the
count. Verified 2026-08-19 by grepping every probed hook name against
`plugin_interface.py`:

- **7 mandatory adapters** (`tutorials`, `dictionaries`, `manifest`,
  `configuration_validator`, `run_semantic_validator`, `artifacts`,
  `cxx_mapping`) call declared v1 `SolverPlugin` members unconditionally. Their
  semantics genuinely are documented upstream; a pointer suffices.
- **4 optional adapters** (`case_compatibility`, `run_document_configuration`,
  `mesh_diagnostic_policy`, `sweep_materializer`) probe hooks that appear in
  **neither `SolverPlugin` nor `SolverPluginV2`**. For these,
  `plugin_capabilities.py` is the only place in the codebase where the
  extension point is described at all.

### 1.0 The three findings are one finding

Those 4 optional adapters are the same 4 that carry ungated cardiac fallbacks
(§1 below). So for exactly these four capabilities, all three problems coincide:

1. the extension point is **undiscoverable** — a plugin author reading the
   public protocol never learns the hook exists;
2. it is **undocumented** — no docstring on the capability either;
3. **not implementing it is silently penalised** — the fallback applies cardiac
   semantics rather than a neutral default.

A third-party author cannot find the hook, is not told it exists, and is
punished for not implementing it. That compound is the honest answer to
"the plugin capabilities are not explained," and it is a much better finding
than a docstring count.

Fourteen probed hooks are duck-typed this way in total (derived by regex
over every `getattr(self.plugin, ...)` probe in the adapters, 2026-08-19). Required-ness is driven
by the explicit `_REQUIRED_PLUGIN_MEMBERS` / `_REQUIRED_V2_MEMBERS` tuples, not
by `Protocol` class membership, so declaring them is purely documentary and
cannot break v1 or v2 plugin loading — see §4.5.

Two conclusions follow.

**The audit targeted the wrong file for its stated purpose.** Line 1 of the
module reads "Internal, focused capability seams for solver plugins." The
public contract an external author implements is `SolverPlugin` in
`plugin_interface.py`. `PluginCapabilities` is core's internal view *over* a
plugin, pointing the opposite direction. Documenting all 21 protocols "for
external plugin developers" documents the wrong side of the boundary.

**The documentation gap and a real defect are the same gap.** Classifying all
25 `legacy_*` fallbacks in `core/compatibility.py` by whether they gate on
plugin identity:

- **13 gated** — import from `plugins/cardiacfoam/` only when the active plugin
  is `org.cardiacfoam`; neutral empty value otherwise.
- **8 ungated** — import from `plugins/cardiacfoam/` regardless of which plugin
  is loaded.
- 4 further `legacy_*` functions import no cardiac code at all.

Six of the 8 ungated are reachable through a `PluginCapabilities` adapter:
`legacy_case_marker`, `legacy_case_runnable_without_workflow`,
`legacy_run_document_config`, `legacy_nondimensional_case`,
`legacy_route_sweep_case`, `legacy_materialize_sweep_case`. They map to exactly
four capability protocols — `CaseCompatibilityCapability`,
`RunDocumentConfigurationCapability`, `MeshDiagnosticPolicyCapability`,
`SweepMaterializerCapability` — and every one of those four is in the
undocumented eleven. The docstring was written precisely when a fallback
decision had to be justified, and skipped where the fallback silently stayed
cardiac.

The remaining 2 ungated functions are out of scope and stay as they are:
`legacy_default_driver_context` encodes "cardiacFoam is the default plugin when
none is supplied," a product decision rather than a leak;
`legacy_generic_case_mutation` serves direct callers of core `make_spec` and is
already documented as a Plan-2 seam.

### 1.1 The cardiac code is not misplaced

It is worth stating explicitly, because it is the natural first reading and it
is wrong: no cardiac logic lives in core. `plugins/cardiacfoam/
case_compatibility.py`, `plugins/cardiacfoam/sweep.py`, and
`plugins/cardiacfoam/planning_policy.py` are all inside the plugin package,
each with a "Why this exists / Activation / Compatibility" header. The defect
is not *where the code lives* but *when core decides to call into it*.

### 1.2 Observable consequence

`GenericOpenFOAMPlugin`, the shipped non-cardiac built-in, implements only one
of the six hooks (`build_run_document_config`). The other five fall through to
cardiac code when the generic plugin is active:

| question core asks | what answers it under the generic plugin |
|---|---|
| does this case belong to my plugin? | scans for `constant/electroProperties*` |
| runnable without a workflow? | additionally demands `constant/physicsProperties` |
| is it non-dimensional? | parses `electroProperties` for `singleCellSolver` |
| route these sweep axes | validates against `electroProperties`/`physicsProperties` vocabulary |
| materialize the sweep case | writes `Allrun` containing `cardiacFoam` |

The first four fail closed: a non-cardiac case has no `electroProperties`, so
they return `False` or reject the axis. The wrong question is asked, but the
answer is safe.

The fifth does not. `materialize_case` in `plugins/cardiacfoam/sweep.py`
hardcodes `"blockMesh\ncardiacFoam\n"` into the generated run script.

Reproduced against the shipped `GenericOpenFOAMPlugin` on 2026-08-19:

```
>>> materialize SUCCEEDED under the generic (non-cardiac) plugin
Allrun contents -> '#!/bin/sh\nblockMesh\ncardiacFoam\n'
files created   -> ['Allrun', 'constant', 'system', 'workflow_contract.json']
```

One qualification, established by the same experiment. Reaching that outcome
requires the routed values to satisfy cardiac validation
(`myocardiumSolver`, a `tissue` from the cardiac enum, a `physics.type` of
`electroModel`/`electroMechanicalModel`). With arbitrary routed values the
call instead raises `ValueError` out of `build_electro_properties` or
`build_physics_properties`. So the capability has two failure modes under a
non-cardiac plugin, and both are defects:

- **common:** a validation error phrased entirely in cardiac vocabulary,
  raised at a user who never selected cardiacFoam — safe but misleading
- **serious:** where the values do validate, a complete cardiac case with an
  `Allrun` invoking the `cardiacFoam` binary

The gate replaces both with one honest error naming the unimplemented hook.

### 1.3 Why the shims exist, and why that reason is empty

Their docstrings state it plainly: Plan-1 compatibility for legacy v1 plugins
written before the hooks existed. But `CardiacFoamPlugin` implements all six
hooks, so the fallbacks are dead code for cardiac; and the only other plugins
in the tree are `GenericOpenFOAMPlugin` and a test fixture, neither of which is
a legacy plugin needing compatibility. The shims protect a population of zero
while mis-serving the one real non-cardiac plugin the project ships.

---

## 2. Publication claim

driverFOAM claims **decoupling evidence**, not extensible platform: the
orchestration core is separated from domain logic, demonstrated by a plugin
protocol, a neutral built-in plugin, and boundary tests.

The stronger claim — that an external project can integrate without core
changes — is stated as roadmap with its outstanding requirements named. This
matches `ARCHITECTURE.md:18`, which already concedes that a publication-strength
portability claim needs an out-of-tree reference plugin and an end-to-end CI
test.

Fixing the ungated fallbacks (§3) is a precondition for even the weaker claim
being told honestly, since the shipped non-cardiac plugin currently inherits
cardiac semantics.

---

## 3. Part A — Gate the six ungated capability fallbacks

Bring the six capability-path ungated fallbacks in line with the 13 gated ones:
return cardiac behaviour only when `plugin_id == "org.cardiacfoam"`, and a
neutral value otherwise.

| fallback | capability | neutral return for non-cardiac |
|---|---|---|
| `legacy_case_marker` | `CaseCompatibilityCapability` | `False` |
| `legacy_case_runnable_without_workflow` | `CaseCompatibilityCapability` | `False` |
| `legacy_nondimensional_case` | `MeshDiagnosticPolicyCapability` | `False` |
| `legacy_run_document_config` | `RunDocumentConfigurationCapability` | empty config, no diagnostics (see below) |
| `legacy_route_sweep_case` | `SweepMaterializerCapability` | raise `SweepValidationError` naming the missing hook |
| `legacy_materialize_sweep_case` | `SweepMaterializerCapability` | raise `SweepValidationError` naming the missing hook |

For the first four, gating changes nothing observable: the cardiac path already
returned the neutral value for a case with no `electroProperties`. Gating makes
the *reason* correct rather than incidental.

`legacy_run_document_config`'s neutral return follows the precedent already set
by its gated sibling `legacy_run_document_config_schema`, which returns
`{"type": "object", "additionalProperties": True}` for non-cardiac plugins —
"this plugin constrains nothing." The builder's matching neutral is an empty
config with no diagnostics: it validates nothing rather than validating against
a cardiac vocabulary the plugin never declared.

For the two sweep hooks, a neutral empty return would silently produce an
unrunnable case, so they must fail loudly instead. The error names the
unimplemented hook and the active plugin, so the message tells a plugin author
what to write. This is the only behavioural change visible to a user, and it
replaces generating an `Allrun` that invokes the wrong solver binary.

Cardiac behaviour is unchanged throughout, because `CardiacFoamPlugin`
implements all six hooks and never reaches a fallback.

### 3.0 Implementation note: four fallbacks cannot currently see the plugin

The 13 already-gated fallbacks take the plugin as their first argument
(`legacy_report_catalog(plugin)`), which is how they test identity. Four of the
six to be gated do not receive it at all:

| fallback | current signature |
|---|---|
| `legacy_case_marker` | `(case_root)` |
| `legacy_case_runnable_without_workflow` | `(case_root)` |
| `legacy_nondimensional_case` | `(spec)` |
| `legacy_run_document_config` | `(spec)` |
| `legacy_materialize_sweep_case` | `(*, case_dir, routed)` |
| `legacy_route_sweep_case` | `(*, base, resolved_axis_values, driver_context)` — has the context, so reaches the plugin already |

So gating requires adding a `plugin` parameter to five of them and updating the
corresponding adapter call sites, which already hold `self.plugin`. This is
safe with respect to the P2.4 instrumentation: `_instrumented` wraps with
`*args, **kwargs` and records only `func.__name__`.

`SweepValidationError` is defined in core (`sweep_expansion.py`, a `ValueError`
subclass), so the sweep fallbacks can raise it without importing from the
plugin.

### 3.1 Tests

1. Per fallback, a direct unit test: cardiac plugin identity yields the cardiac
   result; a foreign identity yields the neutral value or the named error.
2. One boundary test asserting `GenericOpenFOAMPlugin` never reaches
   `plugins.cardiacfoam` through any capability, by driving each capability
   method with the generic plugin active and asserting the module is not
   imported. This is the regression guard: it fails if a future capability adds
   another ungated fallback.
3. Existing cardiac suites must pass unchanged. Any that fail indicate a hook
   whose implementation was assumed to be the fallback.

---

## 4. Part B — Documentation architecture

### 4.1 Audiences and artefacts

No new top-level documents. The project already has the split; the seam
contract lands across what exists.

| audience | artefact | what they get |
|---|---|---|
| agent driving `foamctl` | `AGENT_GUIDE.md`, `describe` JSON | unchanged — agents consume case/config vocabulary, not the seam |
| reviewer, future maintainer | `ARCHITECTURE.md` | the generated seam table |
| reader building on the code | docstrings | public contract in `plugin_interface.py`, seam rationale in `plugin_capabilities.py` |

The direction of the boundary is stated in the table header: `SolverPlugin` is
the public contract a plugin author implements; `PluginCapabilities` is core's
internal view over that plugin.

### 4.2 Docstring rubric

Each of the 21 capability protocols carries prose rationale in the existing
house voice — why the seam exists, what breaks if it is got wrong, as
`CaseFileContractCapability` already does — followed by four fields:

```python
class CaseCompatibilityCapability(Protocol):
    """Whether a case on disk belongs to this plugin, and whether it can
    run without a driverFOAM workflow.

    <prose rationale>

    :adapts: SolverPlugin.has_case_marker, SolverPlugin.is_case_runnable_without_workflow
    :consumed-by: openfoam_driver/core/runtime/registry.py
    :fallback: legacy_case_marker, legacy_case_runnable_without_workflow
    :status: optional
    """
```

- `:adapts:` — the `SolverPlugin` member(s) this renames, or `none`
- `:consumed-by:` — core module path(s) containing real call sites
- `:fallback:` — the `compatibility.py` function(s) used when the hook is
  absent, or `none` for mandatory members
- `:status:` — `mandatory` or `optional`

The 10 protocols with existing prose keep it and gain the field block; the 11
bare ones get both. Part A removes the need for a neutral-vs-cardiac annotation
on `:fallback:`, since after gating every fallback is neutral.

### 4.5 Amendment (proposed, needs approval): declare the duck-typed hooks

**Not in the approved scope. Raised because §1.0 changes what Part B is for.**

Fourteen hooks the adapters probe exist nowhere in `plugin_interface.py`. The
`:adapts:` field would therefore read `none (duck-typed)` for a third of the
table — which records the problem accurately but does not fix it. An author
reading the public contract still cannot discover the extension points.

Proposal: add a `SolverPluginOptionalHooks` Protocol to `plugin_interface.py`
declaring all 14 with signatures and docstrings, explicitly described as
optional and structurally unenforced. Then `:adapts:` names a real, findable
member for every capability.

Safety, verified: `validate_plugin` gates on the explicit
`_REQUIRED_PLUGIN_MEMBERS` and `_REQUIRED_V2_MEMBERS` tuples, never on
`Protocol` membership, and `plugin_interface.py` already documents why those
lists are kept separate ("adding them there would break v1 loading"). Adding a
Protocol class that no list references is inert at runtime — no plugin, v1 or
v2, changes loading behaviour.

Cost: roughly 14 signatures and docstrings. Risk: a future maintainer adds one
of these to a required tuple, which the boundary tests would catch.

Recommended. Without it, Part B documents an undiscoverable API rather than
making it discoverable — and discoverability was the original complaint.

### 4.3 Generated seam table

A "Plugin capability seams" section in `ARCHITECTURE.md`, generated from the
field blocks by a script under `applications/scripts/driverFoam/scripts/`. One
row per capability, columns matching the four fields, ordered as
`PluginCapabilities` declares its fields. The rendered table is committed;
the script regenerates it; the test diff-checks it.

### 4.4 Conformance test

For every field of `PluginCapabilities`:

1. its `Protocol` has a docstring containing all four fields, each non-empty
2. every `:adapts:` target is a real `SolverPlugin` member, or `none`
3. every `:fallback:` names a real `compatibility.py` function, or `none`
4. every `:consumed-by:` module contains at least one call site for that
   capability — subset semantics, so adding a consumer does not break the build
5. the committed `ARCHITECTURE.md` table matches a fresh regeneration

Checks 2–4 are what make the presence test non-vacuous: `"""TODO."""` fails
check 1, and a stale or invented reference fails checks 2–4. All 45 capability
call sites go through `driver_context.capabilities.<field>.<method>`, so check 4
is a mechanical scan.

---

## 5. Scope

**In:**

- `core/compatibility.py` — gate the six capability fallbacks (Part A)
- `core/plugin_capabilities.py` — 21 protocols get rubric plus field block; the
  8 `*Request` dataclasses get a one-line docstring each (parameter objects,
  no field block, their meaning lives in the consuming capability);
  `PluginCapabilities` and `adapt_plugin_capabilities` get full docstrings
- `core/plugin_interface.py` — the 14 undocumented items get prose docstrings,
  no field block; this is the public contract and the first file a reviewer opens
- `ARCHITECTURE.md` — generated seam table section
- new: table generator script, conformance test, Part A tests

**Out, deliberately:**

- the 22 private `_*Adapter` classes — implementation, not contract; the
  conformance test must not require docstrings on them
- `core/generic_plugin.py`, `core/plugin_profile.py`
- `legacy_default_driver_context`, `legacy_generic_case_mutation` — ungated but
  not capability fallbacks, and correct as they stand
- any plugin-author guide, out-of-tree reference plugin, or CI portability test
  — these are the extensible-platform claim's evidence, deferred
- any change to `describe` output or the plugin loading path
- making the six hooks mandatory on `SolverPlugin` — considered and rejected;
  it would require `GenericOpenFOAMPlugin` to implement sweep routing and
  materialization, which is real design work with no current consumer

---

## 6. Non-goals

No capability registry module, no documentation generator beyond the single
table script, no Sphinx or pdoc site, no rewrite of the 10 protocols whose
prose is already good.

---

## 7. Success criteria

1. `GenericOpenFOAMPlugin` reaches no `plugins.cardiacfoam` code through any
   capability, asserted by test.
2. A sweep under a plugin lacking `materialize_sweep_case` fails with a message
   naming the hook, instead of writing an `Allrun` that invokes `cardiacFoam`.
3. All 21 capability protocols carry rubric prose and a valid field block,
   asserted by test.
4. `ARCHITECTURE.md` contains a seam table that a test proves current.
5. Existing cardiac test suites pass unchanged.
6. `ARCHITECTURE.md` needs no new limitation entry about cardiac-shaped
   fallbacks, because after Part A there is none to confess.
