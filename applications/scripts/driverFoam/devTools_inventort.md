# driverFOAM dev-tools inventory

Generated 2026-08-26. Answers one question: for any given file under
`applications/scripts/driverFoam/`, is it (a) reached by driverFOAM at
runtime when you run `driverFoam plan/run/sweep-run`, (b) a maintainer tool
that exists but is never invoked by a simulation, (c) test-only support code,
or (d) not real repo content at all (docs, scratch, generated artifacts)?

## How this was built, and how to re-verify it

This was produced by five parallel research passes (Haiku-model agents, one
per subtree: `core/`, `plugins/cardiacfoam/`, `scripts/`+`postprocessing/`,
`tests/`, repo-root docs/config), then spot-checked by hand against the
actual source before writing this file. **Treat every row as a claim to
re-check, not a fact** — this is a snapshot, not something that stays in
sync with the code. Two things are known to have gone wrong in the first
pass (see "Known correction" below), which is itself evidence this needs
independent verification, not blind trust.

The general re-verification method for any row below:

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors/applications/scripts/driverFoam
# Does anything import this module by dotted path?
grep -rn "module.dotted.path" . --include="*.py"
# Does anything reference this file's distinctive function/class names?
grep -rln "distinctive_function_name" .
```

Grep-based reachability has a specific blind spot: it misses anything
reached through a **fallback/compatibility shim** rather than a direct
import at the call site. See the `reports.py` case study below — that one
file is why "no direct importer" is not proof of "unused."

---

## 1. Runtime-core

Reached (directly or transitively) from `openfoam_driver/cli.py`'s actual
verbs: `plan`, `run`, `step`, `sweep-plan`, `sweep-run`, `describe`. This is
"the real driverFOAM" — code exercised on every case you plan or run.

- Nearly all of `openfoam_driver/core/` (63 of 64 files; the one exception is
  `capability_seams.py`, see §2).
- Nearly all of `openfoam_driver/plugins/cardiacfoam/` (29 of 30 files —
  `reports.py` included, see the case study below).
- `openfoam_driver/cli.py`, `openfoam_driver/__main__.py`, `bin/driverFoam`
  (the shell wrapper that sets `PYTHONPATH` and delegates to
  `python3 -m openfoam_driver`).
- `openfoam_driver/dict_entries.py`, `sweep_materialize.py`,
  `sweep_routing.py` — loose files directly under `openfoam_driver/`, all
  imported by `core/` runtime modules.
- `openfoam_driver/scripts/_dict_keys_scanner.py` — **this is the C++
  dict-key scanner from our earlier conversation.** It is not a dev-only
  tool despite living in `scripts/`. Verified directly:

  ```
  core/strict_planning.py:67   from ..scripts._dict_keys_scanner import catalogued_paths, strict_dict_key_report
  core/strict_planning.py:280  for key in ("unmatched_cxx_reads", "stale_paths", "unmatched_subdicts", "unused_allowlist"):
                                    ...diagnostic("error", f"plugin_dict_key_{key}", ...)
  ```

  Its four output buckets become **hard errors** in `plan --strict`. The
  position-aware checker you preferred (`core/specs/case_dict_keys.py`) is
  wired into the *same function*, a few lines away, but only ever emits
  **warnings** (`uncatalogued_case_dict_key`, warning level, never fails the
  plan). Same call site, two differently-trusted mechanisms — likely the
  root of why the scanner felt confusing to reason about.

- `openfoam_driver/scripts/run_case.sh`, `openfoam_driver/schemas/run-document.json`.

## 2. Wired maintainer tooling (real, invoked by hand, never touched by a simulation)

These exist, are documented, and are exercised by dedicated pytest
conformance tests (in `--check`/drift-detection mode) — but nothing in the
`plan`/`run`/`sweep-run` execution path calls them. No CI config or Makefile
in this repo invokes them automatically either (checked: no `.yml`/`.yaml`
file references any of these paths) — they are **manual commands**, run by
a maintainer after a source change, whose correctness is then pinned by a
test.

| File | What it does | How to call it |
|---|---|---|
| `openfoam_driver/core/capability_seams.py` | Owns parsing/rendering of the plugin-capability seam table | Not called directly — via the CLI script below |
| `scripts/export-capability-seams.py` | Thin CLI over `capability_seams.py`; renders the seam table into `ARCHITECTURE.md` | `python scripts/export-capability-seams.py` to write; `python scripts/export-capability-seams.py --check` to verify `ARCHITECTURE.md` is current without writing (this is what the conformance test runs) |
| `scripts/export-dict-catalog.py` | Exports `dict_entries` + ionic/active-tension catalogs to JSON, one record per `(entry, phase)` pair | `python scripts/export-dict-catalog.py --out /path/to/catalog.json` |
| `scripts/export-report-catalog.py` | Exports the active plugin's post-run report catalog to JSON | `python scripts/export-report-catalog.py --out /path/to/reports.json [--plugin <id\|module:Class\|none>]` (defaults to built-in cardiacFoam) |
| `scripts/export-tutorials-catalog.py` | Exports the tutorial catalog; cross-checks plugin's `get_tutorial_displays()` against `core.runtime.registry.list_tutorials()` so nothing is listed without a backend factory | `python scripts/export-tutorials-catalog.py --out /path/to/tutorials.json` |
| `scripts/export-utility-catalog.py` | Exports the utility-command catalog to JSON | `python scripts/export-utility-catalog.py --out /path/to/utility-catalog.json` |
| `scripts/regenerate-ionic-catalog.py` | Regenerates the ionic-model catalog from OpenFOAM `*_Names.H` headers via `openfoam_driver/scripts/_names_parser.py` | `python scripts/regenerate-ionic-catalog.py` to write; `--check` to verify without writing |
| `scripts/scan-dict-keys.py` | Standalone CLI over the **same** `strict_dict_key_report` used inside `plan --strict` — a human-readable report instead of a plan-embedded diagnostic | `python scripts/scan-dict-keys.py` (default, top 50 rows/section); `--limit N`; `--strict` to exit non-zero on unreviewed drift. Its own `--help` text says outright: *"'driverFoam plan --strict' already runs this same check and reports it as plugin_dict_key_* diagnostics."* |
| `openfoam_driver/scripts/_names_parser.py` | Header parser used only by `regenerate-ionic-catalog.py` | Not called directly |
| `schemas/generate_run_document_schema.py` | Copies the hand-authored `schemas/run-document.json` to the packaged `openfoam_driver/schemas/run-document.json`; a test fails if they drift | `python schemas/generate_run_document_schema.py`, run after any hand-edit to the source copy |

Each exporter/regenerator above has a pytest conformance test (e.g.
`tests/core/test_report_catalog_export.py`,
`tests/core/test_capability_seam_documentation.py`,
`tests/plugins/cardiacfoam/test_ionic_catalog_contract.py`) that runs it via
subprocess in check mode. That means `pytest` catches drift between the
committed output and a fresh render — it does **not** mean the tool runs
automatically outside of a manual invocation or the test suite.

## 3. Runtime code, but for a different consumer than the CLI

`openfoam_driver/postprocessing/` (`plot_builder.py`, `plotting_common.py`,
`style.py`, `table_writer.py`) — real, exported, maintained code. It is used
by external tutorial post-processing scripts *after* a run, not called from
`plan`/`run`/`sweep-run` itself. Not a dev tool, just a different call site
than the simulation path.

## 4. Test-only support (exists only for the test suite)

- `openfoam_driver/plugins/cardiacfoam/ionic_catalog_verification.py` —
  imported only by `test_ionic_catalog_verification.py` and
  `test_ionic_catalog_live_verification.py`.
- `openfoam_driver/scripts/_rtst_scanner.py` — verified: the only reference
  to it anywhere in the repo (besides itself and a comment in
  `_dict_keys_scanner.py` saying "identical pattern to `_rtst_scanner.py`")
  is `tests/drift_guards/test_rtst_enum_contract.py`. A "scanner" module
  reachable from exactly one test file and nothing else.
- `openfoam_driver/scripts/__init__.py` — empty.

## 5. Docs/config — live and actually read by code

`AGENT_GUIDE.md`, `ARCHITECTURE.md`, `CHANGELOG.md`, `KEY_FILES.md`,
`SECURITY.md` — all git-tracked, updated 2026-08-21/25.
`equivalence_protocol.yaml` — read by `dual_run.py`/`test_protocol.py`.
`driverfoam-runtime.yaml` / `.example.yaml` — the `.yaml` is gitignored by
design (per-host copy of the tracked `.example.yaml` template).

## 6. Scratch / generated — not real repo content

- **Root `ROADMAP.md`** (in `applications/scripts/driverFoam/`, *not*
  `future/driverFOAM/`) — confirmed via `git ls-files` that this one is
  **untracked**. Its content is a leftover pytest-coverage note
  ("Generated 92% coverage report... 2026-08-25T21:40Z"), not a document
  anyone wrote. Do not confuse it with `future/driverFOAM/ROADMAP.md`,
  which **is** tracked (`git ls-files future/driverFOAM/` includes it) and
  is the real one.
- `dynamicCode/` — OpenFOAM's `#codeStream` compiled-code cache, regenerable
  per invocation.
- `.coverage`, `.pytest_cache/`, `cardiacfoam_tutorials_driver.egg-info/`,
  `.venv/`, `__pycache__/` — all regenerable, untracked or gitignored, none
  hand-authored.

---

## Case study: the `reports.py` fallback (why grep alone isn't enough)

The first-pass agent for `plugins/cardiacfoam/` called
`plugins/cardiacfoam/reports.py` an orphan: *"defines `CARDIAC_REPORTS` but
never wired; `get_report_catalog()` stub not implemented."* That's wrong,
and wrong in a way worth understanding, because it's exactly the "nested
inside other content" pattern you were asking about — a plugin file reached
through a **core-side fallback**, invisible if you only check "does the
plugin class import this?"

**The mechanism, in order:**

1. `openfoam_driver/core/plugin_interface.py:393` declares
   `get_report_catalog()` as an *optional* method on the `SolverPlugin`
   protocol — plugins may implement it, but don't have to.

2. `CardiacFoamPlugin` (`plugins/cardiacfoam_plugin.py`) does **not**
   implement `get_report_catalog()`. If you only grep that one file for
   `reports`, you find nothing — this is what tricked the first pass.

3. `core/plugin_capabilities.py:961-971` defines `_ReportCatalogAdapter`,
   the thing that actually gets called when driverFOAM needs a plugin's
   report catalog:

   ```python
   @dataclass(frozen=True)
   class _ReportCatalogAdapter:
       plugin: "SolverPlugin"

       def reports(self) -> tuple["ReportDefinition", ...]:
           hook = getattr(self.plugin, "get_report_catalog", None)
           if callable(hook):
               return tuple(hook())
           from .compatibility import legacy_report_catalog

           return legacy_report_catalog(self.plugin)
   ```

   Since `CardiacFoamPlugin` has no `get_report_catalog` attribute, `hook`
   is `None`, `callable(hook)` is `False`, and it falls through to
   `legacy_report_catalog`.

4. `core/compatibility.py:417-428`:

   ```python
   def legacy_report_catalog(plugin) -> tuple:
       """v1 plugins predate get_report_catalog(). Same rule as
       legacy_override_schema: only the built-in cardiac plugin has an
       authored post-run report catalog; other v1 plugins get no reports and
       must declare their own by migrating to v2."""

       if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
           from ..plugins.cardiacfoam.reports import CARDIAC_REPORTS

           return CARDIAC_REPORTS
       return ()
   ```

   This is the only place in the whole repo that imports
   `plugins/cardiacfoam/reports.py`. It's gated on `plugin_id ==
   "org.cardiacfoam"` — a string check, not a class relationship — so a
   static "who imports this class" trace would miss it entirely; you have
   to know to look at the *string-keyed* fallback in a completely different
   package (`core/compatibility.py`) than the file itself
   (`plugins/cardiacfoam/`).

5. The docstring on `ReportCatalogCapability` in `plugin_capabilities.py`
   (line 552-555) actually documents this contract explicitly with
   `:adapts:`, `:consumed-by:`, and `:fallback:` tags — so the design *is*
   documented, just not in the file itself, and not discoverable by
   grepping the plugin package for its own name.

**Why this happened, structurally:** `core/compatibility.py` exists
specifically to let "v1" plugins (predating some newer optional protocol
methods) keep working without every plugin author having to implement every
optional hook — the fallback logic for several such hooks
(`legacy_override_schema`, `legacy_describe_config_resolution`,
`legacy_report_catalog`, `legacy_named_catalogs`) all live there, and all of
them special-case `org.cardiacfoam` by string, reaching into
`plugins/cardiacfoam/*` from outside the plugin package. If you're auditing
reachability into `plugins/cardiacfoam/`, `core/compatibility.py` is a
second root to check, not just `cardiacfoam_plugin.py`.

**How to call it, concretely:** nothing calls `reports.py` directly by name
from outside — you reach it by asking the capability system for the active
plugin's report catalog, e.g. via `scripts/export-report-catalog.py --out
reports.json` (see §2), or internally via
`plugin_capabilities.py`'s `ReportCatalogCapability.reports()` on whatever
capability object driverFOAM resolved for the active plugin.
