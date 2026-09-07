# Normalized Pattern for Building a driverFOAM Dict Catalog (Any Solver)

**Status: methodology derived and verified against driverFOAM's actual
internals this session, for cardiacFoam's `electroProperties` catalog.**
Written so whoever builds a catalog for a *different* OpenFOAM solver under
driverFOAM doesn't have to reverse-engineer this from scratch, and doesn't
repeat the container-vs-leaf confusion documented in
`future/DICT_KEY_ALLOWLIST_BACKLOG.md`.

## 1. The one fact that resolves most confusion: the schema is leaf-only

Read `applications/scripts/driverFoam/openfoam_driver/core/contracts/dictionary.py`
(52 lines, all of it) before building any catalog. The entire schema is:

```python
@dataclass(frozen=True)
class DictEntry:
    driver_path: str
    ...
    value_kind: str = "openfoam_literal"
    dynamic_path: bool = False
    ...

def build_group(defaults, entries) -> tuple[DictEntry, ...]:
    """Apply group defaults without mutating plugin-owned entry objects."""
```

There is **no separate "container" or "block" type**. `build_group` is a
pure Python-source convenience (share `phases`/other defaults across a
tuple of entries written together in the file) — it has *nothing* to do
with OpenFOAM dict nesting. `DictEntry` always represents exactly one leaf
write: a single scalar/word/bool/enum/list/etc. value at the end of a
dotted `driver_path` string like:

```
$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.rPvj
```

Every segment before the last one (`domainCouplings`, `<name>`) is
**implicit structure** — it exists only as a substring of that one leaf
entry's path. It is never itself catalogued, never has its own `DictEntry`,
and structurally *cannot*, because the schema has no field for "this is a
container, here's what belongs inside it."

## 2. How the builder proves this — trace it yourself, don't take it on faith

`applications/scripts/driverFoam/openfoam_driver/specs/dict_builder.py`:

- `_set_nested(node, path, value)` (~line 493): walks `path[:-1]`, calling
  `cursor.setdefault(segment, {})` for each — i.e. it **creates every
  intermediate sub-dict automatically and only ever writes a value at the
  final segment**, `path[-1]`.
- `_serialize_block(tree, indent)` (~line 541): recursively walks the
  resulting nested dict; a Python `dict` value becomes an OpenFOAM
  `key\n{\n...\n}` block, anything else becomes `key value;`. Again:
  structure vs. leaf is decided purely by "is this a dict or not," which
  falls straight out of how many `driver_path` segments got written into
  it by `_set_nested` — never from any catalog metadata about the segment
  itself.
- `<name>`-style placeholders (`_PLACEHOLDER_RE`, ~line 511) mark a segment
  where the OpenFOAM dict author chooses an arbitrary key (e.g. the name of
  one particular PVJ coupling, one particular ECG electrode). `dynamic_path
  = True` is set on the **leaf entries underneath** that placeholder
  (`rPvj`, `pvjKernel`, etc.), never on a `DictEntry` for the placeholder
  segment itself — because there's no such thing as an entry for a
  placeholder segment.

**Conclusion**: if you're tempted to write a `DictEntry` whose `driver_path`
ends at a pure container name (`...domainCouplings`, `...ecgDomains`,
`...ionicConstantOverrides`, etc.) — don't. It will never be read back out
by `dict_builder.py`'s value-population pass (`populate_values`,
`_serialize`), because those only ever match `driver_path`s that terminate
in a real leaf.

## 3. Why the strict scanner still flags container names as "drift"

`_dict_keys_scanner.py` works from the *opposite* direction: it regexes
`.C`/`.H` files for `dict.lookup("X")` / `dict.get<T>("X")` / `dict.found("X")`
/ `dict.subDict("X")` call sites, with no notion of "is this the tail of
some catalog entry's path, or an intermediate segment." Two consequences,
both observed directly in cardiacFoam's `dict_key_allowlist.json` this
session (see `future/DICT_KEY_ALLOWLIST_BACKLOG.md` for the full 90-key
audit):

- A pure container name (only ever `.found("X")` + `.subDict("X")`, never a
  typed `.get`/`.lookup`) shows up in `absent_keys` even when every real
  leaf underneath it is fully catalogued — the scanner has no way to know
  the container itself doesn't need its own entry.
- The *same* name can show up in **both** `absent_keys` and
  `unmatched_subdicts` at once (e.g. `ODESolver`, `constants`, `global`,
  `initialStates` in this repo) — two different regex patterns each fire
  on the same code, one classifying it as a possibly-missing value, the
  other as a possibly-missing subdict group. This is scanner noise, not a
  real double-gap. Don't read a name appearing twice as "extra broken" —
  check the actual C++ read pattern (section 4 below) to see which single
  thing is really going on.

Waiving these in `dict_key_allowlist.json` is therefore *correct*, not
lazy, **for genuine pure-container names**. The mistake to avoid (this
session's actual finding) is waiving something that turns out to be a real
un-catalogued leaf value alongside all the genuine containers, without
checking which is which first.

## 4. The decision procedure — apply this per key, don't pattern-match on the name

For every scanner-flagged key name, find its C++ read site(s)
(`grep -rn '"keyname"' src/ --include='*.C' --include='*.H'`) and classify
by the **read pattern**, not by what the name looks like:

| Pattern seen | What it means | Catalog action |
|---|---|---|
| `dict.get<T>("X")`, `dict.lookup("X")`, `dict.lookupOrDefault<T>("X", default)` where `T` is `scalar`/`word`/`Switch`/`label`/`vector`/`fileName`/`labelList`/`scalarList`/etc. | **Leaf value** | Write a `DictEntry`. `value_kind` from `T`; `enum_values` if the word is drawn from a fixed set (grep call sites / `if(word=="a") ... else if(word=="b")` chains); `typical_value` from the literal default. |
| Only `dict.found("X")` followed by `dict.subDict("X")`, never a typed read of `X` itself | **Pure container** | No `DictEntry` for `X`. Catalog the leaves inside it instead, with `driver_path`s that include `X` as an interior segment. |
| A quoted string compared with `==` against a variable (`var == "X"`), or checked via `someWordList.found("X")` where `someWordList` is a hardcoded C++ list, not a `dictionary` | **False positive** — not a dict key at all | Confirm it's genuinely a string/enum-tag comparison, not a mis-grepped dict read, then leave it waived permanently (document why, as the allowlist file's own description already does for `Vm`/`Iion` in this repo). |
| Real leaf value, but the containing `dict` object is a *different* top-level properties file (e.g. `electroMechanicalProperties` vs `electroProperties`) than the one this catalog module documents | **Out of scope for this catalog file, in scope for another one** | Don't catalog here. If that other properties file has (or gets) its own catalog module, it belongs there instead. **Check whether the same C++ class is also reachable from the dict this catalog *does* cover before deciding this** — see the `preconditioningTime` case below, it's a trap. |
| Read via a generic/shared parent object not owned by this repo at all (e.g. OpenFOAM's own `ODESolver::New(*this, dict)` reading `solver`/`absTol`/`relTol`/`maxSteps` from whatever dict it's handed) | **Out of scope, upstream-owned** | Never catalogued here regardless of which properties file the dict came from — it's not this repo's key to document. Goes in `stale_paths` (scanner blind spot), not `absent_keys`. |

### The trap this session actually hit: dual-scope classes

`preconditioningTime` looked like `D` (out of scope, "belongs to
electromechanics") at first glance, same bucket as `TaScale`/`gamma`.
Wrong: `activeTensionModel::New(...)` (which reads it, via `LandNiederer`)
is called from **two** places —
`sequentialElectroMechanical.C` (passes `electroMechanicalProperties()` —
genuinely out of scope) **and** `singleCellSolver.C` (passes
`electroProperties()` directly, when an `activeTensionModel` subdict is
present — genuinely in scope for *this* catalog). Don't classify a key as
out-of-scope from a single call site. Grep **every** call site of the
function/constructor that reads it (`grep -rn 'ClassName::New\|ClassName(' src/`)
before deciding scope.

## 5. What "using an LLM + driverFOAM's own usage" means in practice

This whole methodology — trace every C++ read site, classify by pattern
not by name, check every call site for dual-scope traps, cross-reference
against what tutorials actually set — is exactly what an LLM agent should
be doing systematically for a new solver's catalog, not something to
shortcut by pattern-matching key names or trusting the scanner's raw
`absent_keys`/`stale_paths`/`unmatched_subdicts` output at face value. The
scanner tells you *where to look*; it does not tell you *what's true*.
Concretely, for a new solver:

1. Read `openfoam_driver/core/contracts/dictionary.py` in full (52 lines —
   there is no excuse to skip it) before writing a single `DictEntry`.
2. Read `openfoam_driver/specs/dict_builder.py`'s `_set_nested` and
   `_serialize_block` to see exactly how `driver_path` strings become
   nested OpenFOAM blocks for *this specific solver's* properties file —
   the `$ELECTRO_MODEL_COEFFS.` prefix / `<solver>Coeffs` resolution logic
   is cardiacFoam-specific; a different solver's builder will have its own
   equivalent root-scope convention that needs tracing the same way, not
   assumed to match.
3. For every dict-read call site the scanner (or a plain grep sweep) turns
   up, apply the table in section 4 — read pattern, not name shape.
4. Grep `tutorials/` (or that solver's equivalent case directory) for every
   candidate leaf key to confirm it's genuinely exercised somewhere, not
   just theoretically readable — this is how `nNonOrthogonalCorrectors`
   would have been caught as "not even this solver's dict" (it's read from
   `system/fvSolution`'s `PIMPLE` block, standard OpenFOAM, not
   `electroProperties` at all) rather than misclassified as a rarely-set
   "legacy fallback" of this catalog.
5. Verify every new entry against the strict scanner before considering it
   done:
   ```bash
   cd applications/scripts/driverFoam && source .venv/bin/activate
   python3 -m pytest openfoam_driver/tests/core/test_strict_planning.py \
       -k allowlist_is_current -q
   ```
   If removing a waiver after cataloging still fails, the `driver_path`
   doesn't actually match what `dict_builder.py`/the scanner resolves
   (dynamic segment mismatch, wrong prefix, wrong nesting) — that's real
   signal, not a reason to force the waiver back in.

## 6. Precedent that does NOT generalize: the ionic-model catalog

`applications/scripts/driverFoam/openfoam_driver/scripts/generate_catalog.py`
and `scripts/regenerate-ionic-catalog.py` look like an existing
"auto-generate the catalog from C++" precedent, but they solve a narrower,
different problem: verifying the **flat** `ionic_model_catalog.py` (states/
algebraic-variable names per ionic model) against the runtime output of
`listCellModelsVariables` — no nesting, no scope resolution, no dynamic
placeholders. It's a drift-detector for a flat enum list, not a generator
for a nested, scope-dependent properties catalog like `electroProperties`.
There is currently no automated generator for the nested case — the
methodology in this document (LLM-driven, read-pattern-based, per-solver)
is the only mechanism that exists for it.
