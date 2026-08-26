# driverFOAM is migrating to the omniD workspace

**Status as of 2026-08-27: migration in progress, next round runs soon.**

The Python orchestrator in this directory is being split into a standalone,
three-package workspace so that projects other than cardiacFoam can drive their
own solvers with it. That work is happening in a separate repository, `omniD`
(local checkout: `~/omnidriver`).

**If you are about to write plugin, packaging, or solver-integration code:
write it in omniD, not here.** The package boundary this directory only
approximates already exists there.

## The shape it is moving to

| package | contains | rule |
|---|---|---|
| `omnidriver` | DAG execution, schemas, provenance, the plugin contract | zero OpenFOAM vocabulary, zero physics |
| `omnidriver-openfoam` | `foamlib` mutators, mesh provisioning, OpenFOAM parsing | knows OpenFOAM, knows no cardiology |
| `omnidriver-cardiacfoam` | electrophysiology, ionic models, the cardiac plugin | depends on both |

Here, all three live tangled in `openfoam_driver/`, separated only by
convention and a handful of boundary tests.

## What each repository is for

- **omniD** — the orchestrator. All future plugin and packaging work.
- **this repository** — the cardiacFoam solver, its C++ sources, its utilities,
  and its tutorials and cases. Real physics content changes belong here.

## The porting rule, and why it needs care

**The two repositories share no git history.** omniD begins at its own initial
commit, so there is no common ancestor: `cherry-pick`, `merge`, and `rebase` are
all unavailable. Every port is a deliberate file-level copy plus a re-run of the
tests on the receiving side.

Both trees have moved independently since the split, and the same fix has
already been made twice by accident — the dead `initialODEStep` key was removed
here in `a9ee8462` and again in omniD's `0f4d877`. Before fixing anything in
either tree, check whether the other has already done it.

## Not yet ported to omniD (as of 2026-08-27)

Measured, not assumed. Each item is present here and absent there:

| work | here | omniD |
|---|---|---|
| explicit `DriverContext` in core (no implicit cardiac default) | 21/21 call sites converted | **7 core files still implicit** |
| `get_phases()` — plugin-declared phase vocabulary | present | **absent**, so the silent validation-skip defect is live |
| cardiac-gated `legacy_*` branches in `compatibility.py` | 2, both documented | **20** |
| cardiacFoam optional-hook coverage | 15/15 | not yet measured |
| `from __future__ import annotations` on `plugin_interface.py` | fixed (`6921e4f2`) | **missing — this is what makes omniD's CI red** |
| guard test for TYPE_CHECKING annotations | present | absent |
| unused declared dependencies dropped | done | `gmsh`, `numpy` still declared, neither imported |

Conversely, omniD has already solved things this tree has not: marker-based path
resolution instead of `Path(__file__).parents[N]`, utility manifests shipped as
package data, a core free of `foamlib`/OpenFOAM imports, cross-package
entry-point discovery, and GitHub Actions CI on Python 3.11 and 3.12.

**The work is complementary in both directions.** Neither tree is simply ahead.

## The one that blocks everything

omniD's `core/plugin_interface.py` annotates types imported only under
`if TYPE_CHECKING:` without `from __future__ import annotations`. On any Python
before 3.14 the module cannot be imported at all:

```
NameError: name 'DictEntry' is not defined
```

omniD's CI matrix is Python 3.11 and 3.12 — exactly the affected versions — so
all six jobs fail at collection (reproduced: 38 collection errors). Three files
need the one-line fix. Nothing else in omniD can be validated until it lands.

This was invisible locally because both repos' virtualenvs are Python 3.14,
where PEP 649 defers annotation evaluation and hides it. **Verify against the
`requires-python` floor, not the interpreter you happen to be running.**

## Related

- `~/omnidriver/GITHUB_MIGRATION.md` — the receiving side's view
- `~/omnidriver/MIGRATION_AUDIT_v2.md`, `~/omnidriver/future/` — omniD's own
  decoupling records
