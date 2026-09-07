# cardiacFoam first release: branch split design

Date: 2026-09-07
Status: EXECUTED 2026-09-07. Landed on `ep-work-onto-main` as `aa47e7d2` and
pushed to PR #24. This document is now a record, not a plan; the outcome
section at the end says what actually happened and where it differed.
Scope: the cardiacFoam repository only. Work belonging to omniDriver or
solids4foam is named where it touches this plan, but is not part of it.

## Goal

Ship the first public cardiacFoam release. Its subject is electrophysiology.
Electromechanics ships as compiled source but is not claimed to run. driverFOAM
does not ship at all: it becomes an external add-on.

## Context

1. `main` and `ep-work-onto-main` have **unrelated histories** (no merge base).
   The release is a wholesale replacement of `main`, delivered as PR #24. This
   is a known past mistake, accepted and not re-litigated.
2. Electromechanics depends on four modified files in the `modules/solids4foam`
   submodule that are not upstream. CI never sees them: it checks out with
   `submodules: false` and points `SOLIDS4FOAM_INST_DIR` at the Docker image's
   `/lib/solids4foam`.
3. driverFOAM's real home is the separate `omnidriver` repository, a
   three-package workspace with an enforced import boundary.

## Decisions

### D1 - `main` is the electrophysiology release, built lightweight only

`main` receives the content of `ep-work-onto-main`. The three
`with-solids4foam` CI matrix jobs are dropped; `main` builds and tests
lightweight only, on all three OpenFOAM versions.

The `with-solids4foam` leg exercises a runtime path whose correctness depends
on unlanded solids4foam patches. Gating the release on it would block the
merge; making it advisory would train everyone to ignore a permanently red
check.

### D2 - electromechanics stays, and is honest about what it does

EM source (`activeTensionModels`, `electroMechanicalModels`, `couplingModels`)
ships and compiles: `src/electroMechanicalModels/Make/options` builds it
against the `modules/physicsModel` shim under `USE_LIGHTWEIGHT_PHYSICSMODEL`.

The EM tutorial `NiedererEtAl2011/electroMechanicalNiedererEtAl2011` also
ships. `tutorials/Alltest-regression` already classifies it as an expected skip
in lightweight mode, so it reports as skipped rather than failing.

Nothing EM-related moves off `main`. The only thing `main` cannot do is a real
solid solve, and that is visible as a skip.

### D3 - driverFOAM leaves the repository

`applications/scripts/driverFoam/` (340 tracked files, 3.8 MB) is deleted,
along with its CI and its repo-level documentation. After removal, no reference
to driverFOAM may remain anywhere **except** inside tutorial `README.md` files,
where documenting an optional external add-on is legitimate.

References under `src/` and `applications/utilities/` are explicitly
unacceptable. Both are comments only:

- `src/ionicModels/ionicModel/ionicModelFamilyInfo.H:20`
- `applications/utilities/checkMeshGeometry/checkMeshGeometry.C:96`

### D4 - the tutorial post-processing layer is left exactly as it is

Twelve tutorial Python files carry 21 unguarded module-level
`from openfoam_driver ...` imports. They are **not modified**.

This was originally assumed to be a coupling problem. It is not. Measured:

- **0 of 12** read driverFOAM output artifacts. No `run_report.md`, no
  `artifacts_manifest.json`, no `workflow_state.json`, no case records. They
  read OpenFOAM's own `postProcessing/` CSVs and tutorial-local JSON produced
  by tutorial-local scripts. This is the only category of coupling that
  packaging cannot fix, and it is absent.
- The three modules they import - `postprocessing/style.py`,
  `postprocessing/plotting_common.py`, `postprocessing/table_writer.py` - carry
  **zero** `openfoam_driver` references between them, total 413 lines, and are
  generic rather than cardiac-specific.
- The only plugin dependency is one dict of five electrode coordinate strings.
- The heavy transitive import surface (42 modules) comes entirely from
  `openfoam_driver/__init__.py` eager imports. That is import weight, not
  orchestrator coupling: nothing is constructed and no configuration is read.

So `pip install` of the external package makes these scripts work unchanged.
The intended long-term fix is to link them from omniDriver once that dependency
is connected. Vendoring the 413 lines remains available as a fallback and is
explicitly not being done now.

Accepted consequence: until omniDriver is installed, running the `setup/`
post-processing layer raises `ModuleNotFoundError`. Case execution is
unaffected - every `Allrun`, `Allclean`, `regressionTest.sh` and
`Alltest-regression` is driverFOAM-free.

### D5 - `dev` branches from `ep-work-onto-main` and keeps everything

`dev` is cut from `929abfc8`, before any deletion, and carries driverFOAM, its
two dedicated workflows, and the three driverFOAM steps inside
`buildAndTest.yml`. It must branch from `ep-work-onto-main`, not from today's
`main`: the two have no common ancestor.

Ordering matters. `dev` is cut first, then driverFOAM is stripped on
`ep-work-onto-main`, then the merge happens. Merging first and deleting
afterwards would put the package in the public history of the release branch,
contradicting the premise that driverFOAM is not part of cardiacFoam.

`dev` is also the backup: nothing deleted from `main` is lost while it exists.

### D6 - the solids4foam patches are preserved, then curated

The four modified submodule files are committed to a branch so they survive.
Choosing which become upstream pull requests is a separate, later decision.

Status at time of writing: preserved as `d28c6527` on
`electromechanical-coupling-wip`. Two of the four changes are now upstream as
solids4foam PRs #414 (Jacobian scaling, fixes #337) and #415
(`fixedValueCorrected` write, refs #409). A third is issue #416. The fvOptions
work (#323/#283) is not started.

### D7 - agent guidance: minimal replacement now, proper split later

`.agents/skills/driverfoam-assistant/` is deleted with the package. Measured,
its content is roughly 76% driver mechanics, 4% cardiacFoam physics, and 14%
seam - it is a driver document that uses cardiacFoam as its worked example, and
its driver half is already superseded by omniDriver's own 1023-line
`AGENT_GUIDE.md`.

`CLAUDE.md` currently devotes 26 of its 30 lines to mandating driverFOAM.
Deleting that without replacement leaves a stub with no guidance at all, which
is worse than the status quo. It is replaced with a short section covering:

- how a case is actually run here: `Allrun`, `Allclean`, `regressionTest.sh`
- the surviving warning from the old mandate, which is still true and still
  valuable: do not hand-roll a shell script that mutates tracked dictionaries.
  This caused real, silent damage once - an ad-hoc script flipped a tracked
  `fvSchemes` default and overwrote `box.geo.template`, undetected until a
  manual audit. The warning survives; only its "use driverFOAM instead"
  phrasing goes.
- one line noting that an optional external orchestrator exists and ships its
  own guidance.

The full three-way skill split - one skill per omniDriver package, matching the
import boundary its CI already enforces - is **deferred** and belongs to the
omniDriver migration, not to this release.

## What ships where

| | `main` | `dev` |
|---|---|---|
| Electro core + tutorials | yes | yes |
| `manufacturedSolutions/` | yes | yes |
| EM source | yes, compiled lightweight | yes |
| EM tutorial | yes, skips | yes, runs |
| Tutorial `setup/` post-processing | yes, unmodified | yes |
| driverFOAM package | no | yes |
| driverFOAM CI | no | yes |
| `with-solids4foam` CI | no | yes |
| solids4foam submodule | clean upstream pin | patched pin |

## Execution order

1. Cut and push `dev` from `ep-work-onto-main`. **Done.**
2. Strip driverFOAM on `ep-work-onto-main`:
   - `git rm -r applications/scripts/driverFoam/`
   - `git rm .github/workflows/driverFoamTests.yml
     .github/workflows/driverfoam_standalone.yml`
   - `git rm -r .agents/skills/driverfoam-assistant/`
   - remove the three driverFOAM steps from `buildAndTest.yml`
     (`Real foamlib mutation tests`, `Live catalog and runtime-dependency
     verification`, `Driver regression-equivalence (phase 2)`)
   - drop the three `with-solids4foam` matrix jobs
   - clean the three driverFOAM blocks in `.gitignore`
   - rewrite `CLAUDE.md` per D7; fix `README.md` and `tutorials/README.md`
   - reword the two source comments
3. Verify: `git grep -inI -e driverfoam -e openfoam_driver` returns hits only in
   `tutorials/**/README.md`.
4. Verify `./tutorials/Alltest-regression` still passes. It is now the only
   end-to-end CI gate, since the driver regression-equivalence step is gone.
5. Merge PR #24.

## Non-goals

- Modifying the twelve importing Python files. Deferred to the omniDriver link.
- Rewriting the roughly 24 `setup/studies/*/README.md` whose only documented
  execution path is a `driverFoam sweep-run` line. They stay.
- The three-way agent skill split. Deferred to omniDriver.
- Landing the solids4foam changes upstream beyond what D6 records.
- Any change to `main`'s submodule pointer.

## Defects found, handled separately

These are real and pre-existing. None blocks the strip; all predate it.

- **The bidomain post-processing carries a verbatim copy of the monodomain
  case's entire pseudo-ECG apparatus, and has been unimportable for months.**
  `tutorials/manufacturedSolutions/bidomain/setup/post_processing_manufactured.py:45`
  imports `manufactured_fda`, removed by `44bed361` when tutorial keys were
  renamed to canonical identifiers; `fe3a0f11` synced docs and tests but not the
  tutorial scripts. `monodomainPseudoECG` imports the correct canonical name and
  is unaffected - bidomain is the only broken script, which fits, since it is the
  copy nobody re-ran.

  The apparatus is dead. The bidomain case configures no ECG at all, yet its
  post-processing carries 105 electrode lines across 24 ECG functions - counts
  numerically identical to `monodomainPseudoECG`'s, which indicates a verbatim
  copy rather than parallel development. `bathBidomain`, which legitimately does
  this work, configures ECG at case level and has a proportionate 1 line and 5
  functions.

  The governing rule: bidomain alone has nothing to do with ECG. Extracellular
  potential that could be compared against an ECG belongs to bath-bidomain. So
  all ECG and electrode content is deleted from bidomain, including the
  electrode-geometry figure and its VTP export; only `phiE` survives, as core
  bidomain physics. Queued as its own task.
- **`openfoam_driver/scripts/run_case.sh:27`** falls back to a hardcoded
  `/Volumes/OpenFOAM-v2412/etc/bashrc`. It leaves with the package, so this
  resolves itself here, but the same path must not reappear in omniDriver.
- **`.agents/skills/driverfoam-plugin-builder/SKILL.md`** (593 lines) is
  untracked and gitignored via `.gitignore:121`, yet referenced from four
  tracked files in omniDriver. It exists only on one machine. It needs a home
  in omniDriver before `.agents/` is touched here.
- **`tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/README.md:38`**
  ships a hardcoded `/Users/simaocastro/omnidriver/...` path. Should not go
  public.

## Risks and accepted consequences

- Post-processing raises `ModuleNotFoundError` without the add-on. Accepted per
  D4, temporary.
- EM has no runtime CI anywhere on `main`. Accepted per D1; compile coverage
  remains.
- `Alltest-regression` becomes the only end-to-end CI gate.
- The two `run_cases.sh` in `singleCell` and `restitutionCurves_s1s2Protocol`
  become dead scripts, since they exec a runner inside the deleted package.
  They are thin wrappers around `Allclean` / `blockMesh` / `Allrun` and can be
  inlined later if wanted.

## Preconditions satisfied

- PR #24's three failing checks are fixed and pushed (`ef885952`, `929abfc8`):
  driverFOAM suite 1767 passed / 0 failed both natively and with OpenFOAM hidden
  from discovery; header normalizer a no-op; markdownlint clean across all 96
  tracked files in the workflow globs.
- `dev` cut from `929abfc8` and pushed.
- solids4foam patches preserved as `d28c6527`.

## Open item

The six build-matrix jobs on PR #24 have sat at `pending` with zero duration
since 11:42. This is infrastructure, not code, but the merge cannot be called
green until it is resolved.

## Outcome

Executed 2026-09-07. Differences from the plan as written, and results.

**The strip landed as `aa47e7d2`** -- 356 files, +52 / -72,656. Verified after
landing: zero driverFOAM references outside `tutorials/`; zero tracked files
left under `applications/scripts/driverFoam/`; `Alltest-regression` and the
submodule untouched; `./Allwmake` completed with "There were no build errors".

**The CI matrix kept `mode` as a single-value dimension** (`[lightweight]`)
rather than removing it, so the three `with-solids4foam` shell branches remain
but are dead. Re-enabling the EM leg later is a one-word change instead of
restoring deleted logic.

**Three further files were reworded** beyond the planned list, because the grep
sweep found references the plan had not enumerated:
`future/STEWART_RESTITUTION_CPP_RECONCILIATION_PLAN.md`,
`monodomain1DCableCV/Purkinje_S1_S2_Calibration.md`, and two case-dict comments.

**Three defects were caught in review** and fixed before landing: `CLAUDE.md`
contradicted itself (it forbade shell for running cases, but `Allrun` is shell
-- the distinction is now between a tutorial's committed `Allrun` and a new
ad-hoc script); a dict comment claimed verification by a suite that no longer
verifies it; and a `tutorials/README.md` section documented the driver's
registry semantics, which passed a keyword grep but described a mechanism no
longer present.

**The bidomain defect was fixed separately and landed first** (`5cb6346d`,
`3c6b5e16`), removing 2,617 lines. It went further than the defects section
above anticipated: rather than repointing a broken import, the whole inherited
pseudo-ECG apparatus came out, on the rule that bidomain alone has nothing to
do with ECG and that extracellular-potential-versus-ECG work belongs to
bath-bidomain. `phiE` was preserved; ECG references are at zero. The strip was
rebased onto that work, with no overlap between them.

**`.tmp/` was deleted** -- 23 GB of driverFOAM sweep output. Two sweep specs
that existed nowhere else were recovered onto `dev` first
(`sweep_stewart_true_di90_late_dt1e-6.json`, and a pre-debug-suppression
variant of the automaticity spec). The written analysis it also held --
`results.md`, `RESULTS.md`, input-sha256 provenance -- was not preserved, a
deliberate trade given the runs are to be redone.

**`dev` became more than "main + driverFOAM".** It also carries both agent
skills (one of which existed nowhere in git), 33 notes and design specs, the
two recovered sweep specs, and the omniDriver migration notes. That was not in
the plan; it followed from the decision to lose no local information.
