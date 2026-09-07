# omniDriver migration: recorded findings

Date: 2026-09-07
Status: recorded, not scheduled. Nothing here is in progress.

Written when driverFOAM was removed from cardiacFoam. Everything below was
measured at that date against `~/omnidriver` and cardiacFoam's
`ep-work-onto-main`. Verify before acting: several items are about drift, and
drift continues.

## 1. Where agent guidance should live

Measured content of the old `.agents/skills/driverfoam-assistant/SKILL.md`
(139 lines): roughly 76% driver mechanics, 4% cardiacFoam physics, 14% seam.
It is a driver document that used cardiacFoam as its worked example, and its
driver half is already superseded by omniDriver's own `AGENT_GUIDE.md`
(1023 lines).

The recommendation is a three-way split matching the boundary omniDriver
already enforces in CI (`scripts/check-import-boundaries.py`, waiver list
empty):

| Skill | Package | Content |
|---|---|---|
| driver mechanics | `omnidriver` | sweep spec schema, strict planning, run documents, the postprocessing brain/module split, `PostprocessingProtocol` |
| OpenFOAM | `omnidriver-openfoam` | dictionaries, meshing, decomposition |
| cardiac | `omnidriver-cardiacfoam` | ionic models, what to read and change |

The organising principle: a skill goes stale when the thing it describes
changes, so it belongs next to whatever can invalidate it.

**The seam should not be prose.** Facts that depend on both sides -- valid
ionic-model/tissue pairs, solver compatibility, the 161 dictionary entries --
must be discovered at runtime (`describe`, the strict planner's diagnostics,
`get_dictionary_catalog()`), not written down. Prose describing a seam rots
from either side and neither repo's CI can catch it. The old skill got this
right in one place (line 76, "the *active plugin's* dictionary catalog") and
wrong in another (naming `IONIC_MODEL_TISSUE_MAP` directly).

**Load-bearing unknown:** whether the harness discovers skills from an
installed package or only from a project's `.agents/`. If only the latter, a
skill inside the omniDriver wheel is invisible to an agent working in a
cardiacFoam checkout, and the split becomes mandatory rather than preferable.
Check this first; it determines what cardiacFoam's own stub must carry.

## 2. Known drift in omniDriver, as of this date

- `packages/omnidriver-cardiacfoam/src/omnidriver/cardiacfoam/plugin.yaml:80`
  has `cxx_mapping.source_roots: ../../../../../../src`. Six `..` copied
  verbatim across a layout change: it resolved correctly inside cardiacFoam
  and resolves to a nonexistent path from omniDriver's root. Whatever the
  external contract should be -- an env var, a flag, or dropping the mapping
  when no checkout is present -- this is not it.
- `AGENT_GUIDE.md` contains references to the retired flat layout
  (`omnidriver/plugins/cardiacfoam/...`) and none to the real package path,
  despite `CLAUDE.md` claiming every module path in it is import-checked.
  `KEY_FILES.md` was fixed for exactly this on 2026-09-03 and carries a banner
  saying a navigational map with stale paths is worse than no map.
- Roughly 258 tests in the cardiac package are `skip_without_monorepo`-gated,
  so a quarter of its coverage cannot run without a cardiacFoam checkout.
- `AGENT_GUIDE.md` carries a 53-line "Adding a New cardiacFoam Tutorial"
  section -- pure cardiacFoam procedure inside the driver's guide, and the
  prose twin of the plugin-inside-the-driver problem.
- Four tracked files (`AGENT_GUIDE.md`, `KEY_FILES.md`,
  `core/plugin_interface.py`, `core/generic_plugin.py`) referenced
  `.agents/skills/driverfoam-plugin-builder/SKILL.md`, which existed only as
  an untracked file on one machine. That skill is now committed on
  cardiacFoam's `dev` branch; the references still need repointing.

## 3. The plugin is itself on the wrong side of the seam

The cardiac plugin is ~13,400 lines -- comparable to the core it plugs into --
and reaches out of itself in three places: the `cxx_mapping` source root above,
the utility manifests (bundled by copy into omniDriver, with drift already
recorded in its own `future/UTILITY_CATALOG_STANDALONE_GAP.md`), and eleven
`CASE_DIR_NAME` literals naming directories in cardiacFoam's `tutorials/`.
All eleven resolved when checked; they are one tutorial rename away from not
resolving, and an MMS tutorial reorg is already planned.

This is why agent guidance for the seam should not live in the plugin either.

## 4. The tutorial post-processing coupling is not a problem

Measured, and worth not re-deriving: **zero** of the twelve tutorial scripts
that import `openfoam_driver` read driverFOAM output artifacts -- no
`run_report.md`, no `artifacts_manifest.json`, no `workflow_state.json`, no
case records. They read OpenFOAM's own `postProcessing/` CSVs and
tutorial-local JSON.

The three modules they use -- `postprocessing/style.py`,
`plotting_common.py`, `table_writer.py` -- have zero `openfoam_driver`
references between them, total 413 lines, and are generic rather than cardiac
(grepping them for cardiac terms returns only licence headers). The heavy
42-module transitive import comes entirely from `openfoam_driver/__init__.py`
eager imports: import weight, not orchestrator coupling.

So `pip install` makes those scripts work unchanged. Vendoring the 413 lines
into a `tutorials/_common/` package remains available as a fallback and was
deliberately not done.

## 5. fvOptions is a blocker for cardiacFoam's own electromechanics

Not an omniDriver item, recorded here because it is easy to mis-scope as one.

On stock solids4foam `development`, `nonLinGeomTotalLagTotalDispSolid` has no
working fvOptions hook at all: `evolveImplicitSegregated` has none, and the
`formResidual` line is commented out. `linGeomTotalDispSolid` has it in both
paths, with the sign discrepancy that solids4foam issue #283 reports.

cardiacFoam's `monodomainTotalLagrangianEM` drives the solid through
`manufacturedSolidForce`, an `fv::option` on exactly that model. It therefore
works only with the local solids4foam patch, and a user building against stock
solids4foam gets no body force -- silently, since a missing fvOption does not
error.

So solids4foam #323 is not only an upstream courtesy: it is what makes
cardiacFoam's electromechanics MMS reproducible by anyone else. Its acceptance
criteria ask for an MMS demonstrating convergence on both the segregated and
PETSc paths; that case already exists in cardiacFoam and needs porting in a
mechanics-only form, without the electrophysiology stack. #283 also names
`linGeomTotalDispSolid.C:936`, which the local patch does not touch.
