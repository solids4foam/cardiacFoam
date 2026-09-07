# driverFOAM: foamlib replaces foamDictionary as the complex-syntax backend

**Date:** 2026-08-20
**Status:** Design, not implemented
**Branch context:** `ep-work-onto-main` (see "Concurrency" below)
**Supersedes:** the first draft of this file, which proposed replacing
`mutators.py` wholesale. That premise was wrong; see "Rejected: full parser
replacement".

## Thesis

`mutators.py` has a three-tier architecture that is not documented as such:

1. a fast, **verbatim** line-based reader/writer;
2. a shell-out to the `foamDictionary` binary, taken **first** on every write
   whenever it is on `PATH`;
3. the line-based path again as the exception fallback.

Tier 2 is the problem. It makes dict bytes depend on whether OpenFOAM happens
to be sourced, it re-serialises values (`0.0` -> `0`, `5.5e-3` -> `0.0055`), and
it *executes* `#calc` / `#codeStream` to produce a value.

**This spec replaces tier 2 with foamlib.** Tier 1 stays exactly as it is.

That is the whole change. It is smaller than the first draft proposed and it is
the only version the evidence supports.

## Rejected: full parser replacement

The first draft proposed rewriting `mutators.py` on `foamlib.FoamFile`,
deleting the line scanner, and collapsing 787 lines to ~200. **This is not
possible.**

`read_foam_entry` returns `str | None` — the *verbatim on-disk token* — and
`mutators.py:216-224` records why: returning parsed values made generated dicts
and provenance digests depend on the environment. foamlib parses to typed
Python values:

```
f["deltaT"]    -> 5e-06          (float; on disk: 5e-6)
f["computed"]  -> ('#calc', ...) (tuple)
```

foamlib exposes **no public API for verbatim source text**. Its public surface
is dict-like only (`as_dict`, `get`, `getall`, `items`, `keys`, `values`,
`dumps`, ...). It tracks byte spans internally (`entry.start` / `entry.end` in
`_files/_parsing/__init__.py`) but does not expose them.

Routing reads through foamlib would therefore re-introduce, on the read side,
the exact defect that `9128ba40..4647b45e` closed. `provenance_inputs.py:106`
and `:116` call `.strip()` on the result; a float raises `AttributeError`.

**Constraint, permanent unless upstream exposes spans: reads stay on tier 1.**

## Rejected: simply deleting tier 2

Deleting the `shutil.which("foamDictionary")` branches without a replacement is
tempting — it is pure deletion and removes the environment dependency. It also
loses real capability.

Measured behaviour of tier 1 against hard syntax:

| dict syntax | tier 1 (line scanner) | foamlib |
|---|---|---|
| `/* block comment { with braces } */` | read+write OK | OK |
| `#include "..."` | read+write OK | OK |
| `#calc "2*3"` | read+write OK | OK |
| `note "value with { brace";` | read -> `None`; write -> `KeyError` | read+write OK, file preserved |

A brace inside a quoted string defeats the brace counter. On read this returns
`None` **silently** — the failure mode that makes an override a no-op. On write
it raises `KeyError: "Scope 'solvers' not found"`.

`AGENT_GUIDE.md:733-747` is wrong on both counts: it blames block comments,
`#include` and `#calc` (all handled), and it describes `foamDictionary` as a
fallback when the code takes it first.

## Evidence

Produced by running foamlib 1.7.5 against this repo, not by reading its docs.

### foamlib facts

- v1.7.5, GPL-3.0-only (matches cardiacFoam — no new licence obligation),
  requires-python `>=3.11`.
- Deps: `numpy`, `rich`, `aioshutil`, `multicollections`, `typing-extensions`.
  Pure Python; no C extension, no pyparsing.
- Hand-written parser (~1531 LOC) with ~2500 LOC of parsing tests including an
  explicit OpenFOAM grammar-tolerance suite. Actively maintained.

### Parse coverage

Every `FoamFile`-headed file under `tutorials/`: **10,830 parsed, 0 failures**
(~326 MB), covering every `system/`, `constant/` and small `0/` dictionary —
the complete set `mutators.py` touches. A 633-file / 5.5 GB tail of `polyMesh/`
and large field data was still running at time of writing; it informs only the
deferred field-I/O question, not this spec.

### Directives are inert and typed

```
pwned      -> ('#codeStream', {'code': '#{ os << system("touch /tmp/PWNED"); #}'})
computed   -> ('#calc', '"2.0*3.0"')
alias      -> '$computed'   (string; not expanded)
codeStream executed?  False
```

Never evaluated on read. Directives come back as a tuple whose head is the
directive name — a structural discriminator for injection screening that
`foamDictionary` cannot provide, because it executes them.

### Writes preserve the file

Editing `endTime` in `tutorials/template/system/controlDict` left the banner,
all 13 comment blocks and every sibling entry byte-identical. foamlib splices
into the original contents using tracked spans.

Two deltas from tier 1's output, both of which the wrapper must normalise:

- **Separator width.** foamlib writes `endTime 250;` (one space); tier 1 writes
  four — `f"{key}    {_format_value(value)};"` at `mutators.py:422`, `:427`,
  `:447`. Five test assertions hard-code the four-space form
  (`test_mutators.py:149, 370, 414, 554, 567`).
- **One inserted blank line** before the edited entry. Verified **idempotent**
  across four successive writes (blank-line count 17, length 2061, unchanged);
  it does not accumulate. OpenFOAM parses the result identically.

### foamlib auto-creates missing keys, silently

```
f["endTimee"] = 999      # typo
-> no error; key created
```

`FoamFile.__setitem__` has no fail-closed mode. `update_foam_entry`'s
`add_if_missing=False` default exists precisely because `foamDictionary -set`
auto-creating bogus keys has bitten this project before
(`specs/apply_overrides.py:241-242`). The wrapper must preserve fail-closed
behaviour or a mistyped override key produces a valid-looking dict, a passing
plan, and quietly wrong results.

### Booleans serialise correctly

`True`/`False` -> `yes`/`no`, matching `_format_value` (`mutators.py:37-40`).

### foamlib is type-strict on write (measured)

Every override value arriving from `sweep.json` or the CLI is a `str`, which
`_format_value` passes through as raw text. foamlib refuses strings that would
read back as another type:

```
'1e-6'            -> ValueError: cannot be stored as a string
'0.0'             -> ValueError
'(1e6 1e6 1e6)'   -> ValueError
'uniform 0'       -> ValueError
'3D'              -> ValueError: invalid string   (digit-prefixed)
'#calc "2*3"'     -> ValueError: invalid string
'1e-6;  rogue  1' -> ValueError: invalid string   (injection payload)
'yes'             -> 'k yes;'  (UserWarning: stored as True)
'PCG'             -> 'k PCG;'
'"Vm|VmFinal"'    -> 'k "Vm|VmFinal";'
'$FOAM_CASE/x'    -> 'k $FOAM_CASE/x;'
```

Two consequences, one favourable and one constraining.

**Favourable:** a value carrying a stray `;` plus a second entry, or a `#calc`
directive, is rejected *at the write call*. This closes the value-channel
injection gap by prevention, not detection — stronger than Step 4 anticipated.
It also mirrors `dict_builder._openfoam_value_token`'s digit-prefix quoting
rule, which is why `'3D'` is refused unquoted.

**Constraining:** the wrapper must coerce str -> proper Python type before
calling foamlib, and foamlib then respells (`"1e-6"` -> `1e-06`). Byte-identity
with tier 1 is therefore **impossible** for numeric-looking string values on the
foamlib path.

This does not sink the design, because foamlib occupies tier 2 only, and tier 2
already respelled — `foamDictionary` turned `5.5e-3` into `0.0055`. Respelling
on the complex-syntax fallback is no worse than the status quo, and tier 1
continues to write verbatim on the normal path. What it does invalidate is a
byte-identity exit criterion for tier-2 cases: no verbatim baseline exists,
because there never was one.

## The five tier-2 branches

```
core/runtime/mutators.py:277   update_control_dict
core/runtime/mutators.py:393   update_foam_entry
core/runtime/mutators.py:495   remove_foam_dict
core/runtime/mutators.py:754   ensure_foam_dict
specs/apply_overrides.py:332   controlDict override route (its own branch)
```

The four in `mutators.py` swallow all exceptions (`except Exception: pass`) and
fall through to tier 1. The one in `apply_overrides.py` has no fallback and
bypasses `update_foam_entry` entirely.

The verification is gated the same way: `test_apply_overrides.py:200` and `:216`
skip unless `REQUIRE_FOAMDICTIONARY=1` (which no CI job sets), and
`test_mutators.py:577` is `@skipUnless`-gated. The divergence between tiers is
asserted nowhere.

## Public API — actual signatures

`mutators.py` is **787 lines** and exports **seven** public functions:

```
read_foam_entry(file_path, key, *, scope=None) -> str | None
update_foam_entry(file_path, key, value, *, scope=None, add_if_missing=False) -> None
remove_foam_entry(file_path, entry_name, *, scope=None, missing_ok=False) -> None
read_foam_dict_block(file_path, dict_name, *, scope=None) -> str | None
ensure_foam_dict(file_path, dict_name, block_text, *, scope=None) -> bool
remove_foam_dict(file_path, dict_name, *, scope=None, missing_ok=False) -> None
update_control_dict(control_dict_path, *, delta_t=None, end_time=None,
                    start_time=None, write_interval=None, write_control=None,
                    write_format=None, purge_write=None) -> None
```

All seven are preserved unchanged. `remove_foam_entry` (`:566`) is used by
`plugins/cardiacfoam/overrides.py:5-10`.

### Verbatim block carry-forward is a hard constraint

`read_foam_dict_block` (`:634-660`) returns raw on-disk text so it can be
replayed verbatim through `ensure_foam_dict`'s `block_text` parameter.
Live users: `plugins/cardiacfoam/dict_builder.py:460` (`_capture_dynamic_containers`)
-> `:490` (`_carry_forward_dynamic_containers`), documented as *"intentionally a
raw, unmodified carry-forward, not a resynthesis"*, guarded by
`tests/plugins/cardiacfoam/test_dict_builder.py:1485`.

**These two stay on raw text splicing regardless of backend.** A
parse -> re-serialise round trip drops the comments and formatting this path
exists to preserve.

## Architecture

### Division of responsibility

| Concern | Owner |
|---|---|
| Fast path: verbatim reads, simple writes | tier 1 (unchanged) |
| Complex-syntax fallback (parse, locate, splice) | **foamlib `FoamFile`** (replaces `foamDictionary`) |
| Fail-closed key policy (`add_if_missing`, `missing_ok`) | wrapper |
| Separator/blank-line normalisation to tier-1 byte form | wrapper |
| Exception mapping | wrapper |
| Verbatim block carry-forward | tier 1, permanently |
| Planning, provenance, manifests, sweeps, process launch | driverFOAM (unchanged) |

### Exception contract

`specs/apply_overrides.py:343` catches a fixed tuple —
`except (OSError, KeyError, ValueError, RuntimeError)` — and re-raises as
`OverrideError`. `plugins/cardiacfoam_plugin.py:283,288` and
`core/runtime/sweep_runner.py:222` depend on `KeyError` specifically.

**No foamlib exception type may escape `mutators.py`.** The wrapper catches
foamlib's hierarchy and re-raises as `FileNotFoundError` / `KeyError` /
`ValueError` per the existing contract.

## Implementation sequence

### Step 0 — Measure string serialisation — DONE

Result recorded under "foamlib is type-strict on write" above. Outcome:
proceed, with a str -> type coercion layer in the wrapper and a revised
verification criterion in Step 2.

### Step 1 — Packaging

- Bump `requires-python` from `>=3.10` to `>=3.11` in `pyproject.toml:10` and
  `uv.lock:3`. Both CI jobs already pin 3.11; the only in-tree guard is
  `utility_catalog.py:104`. This is a declared-contract break, not a functional
  one, but it must be stated: **3.10 support is dropped.**
- Add `foamlib>=1.7.5,<2` to `[project].dependencies`.
- **`numpy` becomes a core dependency** (it is currently only in the `[post]`
  extra). This changes what a minimal install pulls.
- Regenerate `uv.lock` — the standalone CI job runs `uv run` and resolves from
  the lock.

### Step 2 — Replace tier 2 with foamlib

Swap the five branches to route through foamlib instead of `foamDictionary`.
Delete `update_foam_entry_via_foamDictionary`, `remove_foam_dict_via_foamDictionary`,
`ensure_foam_dict_via_foamDictionary`, the four `shutil.which` checks, the
`subprocess` import, and the `apply_overrides.py:332` branch plus its import at
`:55`. Tier 1 and the `re` import are untouched.

Wrapper normalises output to tier-1 byte form where it can: four-space
separator, no inserted blank line. It also coerces `str` override values to
proper Python types before calling foamlib, and re-raises foamlib's type-strict
`ValueError` as the driverFOAM contract's `ValueError` (see Exception contract).

**Test work — this breaks collection, not just assertions:**

- `test_mutators.py:37-43` imports `update_foam_entry_via_foamDictionary` at
  **module level**. Deleting the name is an `ImportError` that fails collection
  of all 21 tests in the module.
- Eleven `mock.patch` targets disappear: `test_mutators.py:142, 181, 324, 494,
  531, 534, 546, 549, 562` and `test_apply_overrides.py:246, 249`. `mock.patch`
  raises `AttributeError` on a missing attribute.
- The dual-path tests these parameterise collapse to single-path tests.
- The `foamDictionary`-gated tests (`test_apply_overrides.py:200, 216`,
  `test_mutators.py:577`) are deleted or made unconditional.

**Verification — differential harness (new code).** The existing
`test_mutators.py` suite contributes **no corpus coverage**: all 21 tests build
synthetic inline fixtures (`path.write_text("a { b 1; }\n")` at `:154`,
`path.write_text("deltaT 1e-06;\n")` at `:528, 543, 559`). It cannot be reused
for this.

A new harness enumerates the `FoamFile`-headed dicts under `tutorials/`, applies
a generated mutation set per file (edit first scalar leaf; edit first nested
leaf; add-if-missing into first sub-block; remove first sub-block) through both
the pre-change and post-change implementations. The pre-change implementation is
retained as a **test-only** reference and deleted in the same PR.

**Two assertion classes, because byte-identity does not hold everywhere:**

- *Tier-1 cases* (the overwhelming majority — dicts the line scanner handles):
  assert **byte equality**. foamlib is not invoked; output must not move.
- *Tier-2 cases* (dicts that defeat the line scanner, e.g. brace-in-string):
  byte-identity is unavailable — the old path was `foamDictionary`, which
  respelled, and is only present when OpenFOAM is sourced. Assert instead:
  the target entry holds the intended value; every sibling entry and comment is
  unchanged; the file re-parses; and no `#calc`/`#codeStream` was evaluated.

**Behaviour change to record:** in a *sourced* environment, output bytes change
— from `foamDictionary`'s re-serialisation to tier 1's form. That is the point
(determinism), but it means provenance digests move for anyone who previously
ran sourced. `provenance.py:88` and `capture_provenance.sh:9-14` sha256 file
bytes. No golden baselines are committed, so nothing in-tree breaks; the
exposure is comparing archived Paper I provenance JSON against a fresh re-run.
Record in `CHANGELOG.md`; never absorb silently.

**Docs:** rewrite `AGENT_GUIDE.md:733-747` (wrong on both direction and cause)
and `specs/mesh_provisioning.py:46` (refers to "`foamDictionary`-based mutation
machinery").

**Exit criteria:** all five branches gone; no `foamDictionary` reference in the
write path; the brace-in-string case reads and writes correctly; differential
harness green under both assertion classes; a `str` override of shape `"1e-6"`
still succeeds end-to-end through `apply_overrides`; reference implementation
deleted; docs corrected.

### Step 3 — Regex-key resolver (independent of backend)

`tutorials/template/system/fvSolution` keys its solvers block as
`"Vm|VmFinal|u|uFinal"`. OpenFOAM resolves `Vm` against that pattern. **Neither
tier 1 nor foamlib does** — `_find_dict_block_bounds` uses `re.escape`, and
foamlib raises `KeyError`. This is a pre-existing gap, not a migration cost.

Resolver: exact literal key wins; otherwise walk `"`-quoted sibling keys in
reverse declaration order and take the first whose regex fully matches.

**This changes behaviour in the user-facing direction.** `solvers/Vm/tolerance`
currently fails closed and writes nothing; after this it resolves and mutates.
That is the intended fix — it removes a real limit on what an agent can express
— but it must land separately so the gain is attributable.

**Exit criteria:** `solvers/Vm/tolerance` resolves on the template `fvSolution`;
literal-over-pattern precedence tested; reverse-declaration tie-break tested; a
name matching no pattern still fails closed.

### Step 4 — Directive rejection

Step 0 showed foamlib already refuses directive-shaped and semicolon-injected
string values at the write call, so the tier-2 path is covered for free. This
step extends the same rule to **tier 1**, which still writes any string
verbatim: reject override values that are or contain OpenFOAM directives before
either tier sees them, using the parsed-directive tuple as the discriminator.

Without this step the protection is accidental and partial — present only on the
complex-syntax fallback, absent on the path almost every override actually
takes. Update `SECURITY.md` to state the control is enforced and by which
mechanism.

## Non-goals

**Parsers not migrated.** After this work the repo still contains four
hand-rolled OpenFOAM parsers. The headline is "tier 2 retired", not "one parser".

| File | Disposition |
|---|---|
| `specs/function_object_fields.py:49-95` — `_strip_comments`, `_balanced_block`, `_iter_subdicts` over `controlDict` | Migrate later. Second parser on a file `mutators.py` writes; feeds the strict planner's `unknown_sampled_field` check |
| `plugins/cardiacfoam/detection.py:20-43` — brace-depth counting over `electroProperties` | Migrate later |
| `plugins/cardiacfoam/runtime_evidence.py:122-148` — paren-depth over `controlDict` `libs (...)` | Migrate later |
| `specs/utils.py:38-73` — `replace_block_mesh_resolutions`, `partition(") (")` on `blockMeshDict` | Never; `dx`-as-sweep-axis string surgery |
| `specs/dict_builder.py:307-347` — `_openfoam_value_token`, `_serialize_block` | Never (synthesis, not mutation). Note: its digit-prefix quoting rule has no counterpart in `_format_value`, so the two writers already disagree |
| `specs/mesh_geometry.py:81-91` — polyMesh binary headers | Never; meshes out of scope |

**Meshes.** foamlib can parse `polyMesh/faces`, but mesh handling stays as-is.

**Field I/O.** `FoamFieldFile` (binary/compressed -> numpy) is a genuine gap in
driverFOAM and a clean fit, but it belongs to result-interpretation and gets its
own spec.

**Case orchestration.** `FoamCase`, `AsyncFoamCase`, `AsyncSlurmFoamCase` are
declined — they overlap driverFOAM's planner, manifests, provenance and sweep
orchestration and are thinner than what exists. Deferred, not precluded:
driverFOAM's execution layer remains the sole owner of process launching, so a
Slurm backend can slot in behind it later.

**Not a performance project.** foamlib parses this corpus at roughly 45 MB/s.
Adequate, not transformative, and not a justification.

## Concurrency

Another session was committing to `ep-work-onto-main` while this spec was
written (commits at 23:40 and 00:14 on 2026-08-19/20, plus 36 staged files).
`mutators.py`, `specs/apply_overrides.py`, `tests/core/test_mutators.py` and
`pyproject.toml` were all clean at the time of the audit. **Re-verify every
line number in this spec before implementing** — they were accurate as of
commit `12f1189a`.

## Open items

- Upstream issue for the inserted blank line; not blocking.
- Large-field tail (633 files / 5.5 GB) — informs the deferred field-I/O spec.
