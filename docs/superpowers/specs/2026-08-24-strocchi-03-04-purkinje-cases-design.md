# Strocchi 03/04 Purkinje case preparation — design

Status: approved by user 2026-08-24. Executing per
`docs/superpowers/plans/2026-08-24-strocchi-03-04-purkinje-cases.md`.

## Discovered during execution (not known when this spec was written)

- **EnSight filename mangling:** VTK's `vtkGenericEnSightReader` silently
  strips the leading digit run from any referenced filename shaped
  `<digits>.<ext>` (`03.geo` → `.geo`), confirmed by direct testing to be
  independent of the `.case` file's own name and of the `TIME` block.
  Filenames with an underscore instead (`03_uvc_transmural.ens`) are
  unaffected. Fixed with a symlinked "sanitized" bundle (originals
  untouched) — see the plan's Task 2/8 Step 1b.
- **`sheet` field silently corrupted on import:** `newVtkUnstructuredToFoam`
  only correctly reconstructs one 3-component CELL_DATA array as a
  `volVectorField` per import (whichever pyvista wrote as VTK's single
  "active vectors" — `fiber`, in this pipeline); any other 3-component
  array lands in a generic legacy-VTK `FIELD` block, and the importer's
  FIELD-block handling flattens multi-component sub-arrays into a bogus
  `volScalarField` instead of reconstructing them as vectors. `sheet` was
  affected; confirmed pyvista's own reader parses it correctly as
  `(n_cells, 3)`, isolating the bug to this importer. Fixed by re-importing
  `sheet` alone through a second VTK file with it promoted to the
  active-vectors slot, then copying just the resulting `0/sheet` over — see
  the plan's Task 3/9 Step 2b. **Always verify `grep class .../0/sheet`
  reads `volVectorField` after any import.**
- **Mesh units:** see the millimetre-vs-metre note in Stage A step 4 below
  — the mesh must stay in millimetres, not SI metres, to match the
  validated `generatePurkinjeTreeDict` growth parameters.
- **`1DgraphToFoam` also assumes SI-metre input:** its `-graphStep` default
  is `3e-4` *metres*, which on this millimetre-scale mesh subdivided a
  ~44k-node tree into ~45 million edges before being caught. Always pass
  `-graphStep 0` (disable subdivision — `generatePurkinjeTree` already
  discretizes at `l_segment` resolution) given the mesh stays in mm.
- **Non-deterministic uninitialized-memory bug in `newVtkUnstructuredToFoam`
  on case 04:** the sheet-field reimport fix (promote `sheet` to VTK's
  active-vectors slot, reimport, copy `0/sheet` over) produced a
  `defaultFaces` boundary value of denormalized garbage
  (`7.9e-323, 2.3e-318, 1.2e-322`) on the first attempt, causing a
  downstream "Imbalanced brackets" parse failure in `setCardiacConductivity`
  — confirmed the internalField itself was 100% well-formed (every one of
  2,404,586 tuples matched a strict numeric-tuple regex; the file's closing
  structure was syntactically perfect) and the source `sheet` array had no
  NaN/Inf/zero-norm entries, isolating the defect to an uninitialized
  variable in the importer that isn't always zeroed. Re-running the exact
  same reimport produced a clean `uniform (0 0 0)` boundary value instead.
  Not reproducible on case 03 (smaller mesh) in this session — re-run the
  reimport and re-verify the boundary block if this recurs, rather than
  assuming a first success/failure is deterministic.
- **LV/RV tag assignment has no written reference in this dataset** — no
  StrocchiData README, no cardiacCoreStandalone doc states which of
  `tags ∈ {1,2}` is LV vs RV. Verified instead by rendering the mesh:
  tag 1 forms the apex (LV's anatomical signature; RV normally doesn't
  reach the apex), tag 2 sits above it with an open vessel-like stump (RV
  outflow-tract-like geometry), and the near-apex `uvc_intraventricular`
  slice shows tag/sign 1 (−1) as an almost-complete ring with 2 (+1)
  reduced to a sliver — consistent with LV dominating near the apex. All
  three signals (3D shape, apex slice, 79:21 volume split matching typical
  LV:RV mass ratio) agree: `uvc_intraventricular < 0` → LV, `> 0` → RV,
  confirming `uvcConventionDict`'s existing convention holds for this
  anatomy too. Repeat this visual check for case 04 rather than assuming it
  carries over.

## Goal

Prepare 4 ready-to-transfer OpenFOAM case directories from raw Strocchi
biventricular anatomy data (`/Users/simaocastro/Documents/Data/StrocchiData`,
cases `03` and `04`), each carrying a converted mesh, initial fields, and a
`constant/purkinjeGraph` — but no solver/`electroProperties` wiring. The
solver ingests the mesh/graph the same way regardless of anatomy, so that
part is explicitly out of scope; the deliverable is data, meant to be copied
to a workstation for the actual run.

Healthy only — no scar code, no ischemia overrides.

## Scope: 4 cases

| Case | Anatomy | Purkinje scenario |
| --- | --- | --- |
| `bivCase_Strocchi03_human` | StrocchiData `03` | human normal: all-leaves, endocardial |
| `bivCase_Strocchi03_pigExtended` | StrocchiData `03` | morphometric, transmural-extended |
| `bivCase_Strocchi04_human` | StrocchiData `04` | human normal: all-leaves, endocardial |
| `bivCase_Strocchi04_pigExtended` | StrocchiData `04` | morphometric, transmural-extended |

Human and pig scenarios for the same case number share the identical
converted mesh and base fields (fiber/sheet/uvc_*) — the Purkinje generator
config is the only thing that differs within a pair, matching the existing
`bivCase` / `bivCase_human_tree_v2` / `bivCase_Pig_Morphometric_Tree`
precedent in `cardiacCoreStandalone`.

**Explicitly not used:** `StrocchiData/03-350um.vtk` and `04-350um.vtk`.
These were checked and are separate ~3 GB whole-torso meshes (`elemTag`,
fiber, sheet; ~59.7M cells; no UVC fields), matching the KCL/Niederer
torso-tag convention from an unrelated project — not the biventricular
source. The real source per case is its `.case`/`.geo`/`.ens` EnSight bundle,
which already carries fiber/sheet (per the user's explicit instruction to use
the data that already has fibers).

## Pipeline

### Stage A — EnSight → OpenFOAM mesh conversion (once per anatomy, i.e. once for 03, once for 04)

1. `caseToASCIIlegacy.py --input 0X.case` (in `StrocchiData/`) merges the
   EnSight blocks (point data: fiber, sheet, tags, uvc_transmural,
   uvc_intraventricular, uvc_longitudinal, uvc_rotational,
   electrode_endo_rv) into `ASCIIlegacy0X.vtk`.
2. `biventricular_filter.py --input ASCIIlegacy0X.vtk` keeps only
   `tags ∈ {1, 2}` → `ASCIIlegacy0X_biventricular.vtk`.
3. `newVtkUnstructuredToFoam ASCIIlegacy0X_biventricular.vtk -case <caseDir>`
   (repo-local tool in `cardiacCoreStandalone/src/newVtkUnstructuredToFoam`,
   already compiled) writes `constant/polyMesh` and imports all point/cell
   arrays as OpenFOAM fields into `0/`.
4. **Mesh stays in millimetres — do not scale to SI metres.** The
   `newVtkUnstructuredToFoam` README recommends `transformPoints -scale
   '(0.001 0.001 0.001)'`, and Stage A originally did this, but it was
   reverted after discovering (mid-execution, on case 03) that the
   validated `generatePurkinjeTreeDict` growth parameters reused from
   `bivCase_human_tree_v2`/`bivCase_Pig_Morphometric_Tree`
   (`l_segment 0.3`, `initLength 6-30`, `length 5-6`, ...) are themselves
   millimetre-scale absolute lengths, not ratios — confirmed by the
   historical `2026-06-11-purkinje-seed-deduction-RESULTS.md`, which states
   the mesh's own endo edge length is "~0.925 mm". Scaling the mesh to
   metres while reusing those parameters unchanged would make a single
   `l_segment` march step (0.3 "units") 300 mm — larger than the whole
   heart — and degenerate the tree immediately, the same failure mode as
   the historical bad-seed stub case. Unit consistency for the eventual
   solver run (conductivity in S/m, capacitance, etc.) is explicitly out of
   scope here per the user ("we can scale only in the solver itself...
   solver agnostic") — only the millimetre-scale Purkinje-generation
   geometry is this pipeline's concern.
5. `checkMeshGeometry` (bounding box sanity) and `checkMesh` (topology) before
   any cardiacCore utility touches the case. `checkMeshGeometry` should
   report "Max dimension ... suggests mesh is in mm, not meters" — that is
   the expected, intended result here, not a problem to fix.
6. All imported fields except `Conductivity` are genuinely dimensionless
   (fiber, sheet, tags, uvc_*) and need no dimension fix. `Conductivity` is
   produced later by `setCardiacConductivity`, not by the importer, so the
   dimension fix (Stage B) applies there, not here.

This produces one converted case per anatomy. Clone it into the two
scenario-specific directories (human, pig) before Stage B, so scenario runs
don't clobber each other's `0/` outputs.

### Stage B — cardiacCore pipeline (per scenario case dir, via `agent/pipeline_runner.py`)

Plan-driven, not raw shell — reuses the existing validated runner in
`cardiacCoreStandalone/agent/`, following the exact step graph already used
by `bivCase_human_tree_v2` / `bivCase_Pig_Morphometric_Tree`
`pipeline_plan.json` files.

Steps, in dependency order:

1. `setCardiacConductivity`: `df 0.1334; ds 0.01761; dn 0.01761; fiberField fiber; sheetField sheet;`
   (monodomain, transversely isotropic). These are not the orthotropic
   `cardiacCoreStandalone` defaults (0.1143/0.052/0.016) — per the user's
   explicit instruction, replaced with the "typical" literature value
   already used elsewhere in this repo's own `electroProperties` fixture
   template (`applications/scripts/driverFoam/openfoam_driver/specs/fixtures/template/constant/electroProperties`),
   which is the harmonic-mean monodomain-equivalent of the widely-cited
   Niederer et al. (2011) bidomain conductivities (σ_il=0.17, σ_el=0.62 S/m
   longitudinal → 0.1334; σ_it=0.019, σ_et=0.24 S/m transverse → 0.01761 —
   verified by direct computation, not assumed).
   Fix `Conductivity` dimensions after this step: `[0 0 0 0 0 0 0]` →
   `[-1 -3 3 0 0 2 0]` (S/m), per `newVtkUnstructuredToFoam`'s documented
   post-import step.
2. `setCardiacAnatomy`: `zApicalMid 0.3333333; zMidBasal 0.6666667; zApexCap 0.08; grooveMode auto;`
3. `setPurkinjeMorphometry` (pig scenario only — human's `allLeaves`
   selection doesn't consume the weight fields): `grooveMode auto;`
4. Seed/direction deduction (see below) → write `generatePurkinjeTreeDict`.
5. `generatePurkinjeTree`.
6. `1DgraphToFoam <purkinje_vtk> -case <caseDir> -name purkinjeGraph` →
   `constant/purkinjeGraph`.

`system/uvcConventionDict` is identical across all 4 cases (fixed repo
convention, not an anatomy fact):

```
uvc { transmuralField "uvc_transmural"; intraventricularField "uvc_intraventricular"; longitudinalField "uvc_longitudinal"; rotationalField "uvc_rotational"; }
intraventricularChambers { LV -1.0; RV 1.0; }
transmural { min 0.0; max 1.0; }
```

### Seed/direction deduction (per anatomy — cannot reuse `bivCase`'s hardcoded XYZ)

`cardiacCoreStandalone/agent/deduce_purkinje_seeds.py` is referenced in
`docs/superpowers/specs/2026-06-11-purkinje-seed-deduction-RESULTS.md` but is
missing from the checkout. Reimplement the documented method as a one-off
analysis script (reads the OpenFOAM `AHA_Segment`/`aha_angle`/`uvc_*` fields
after Stage B step 2, not a cardiacCore C++ change):

- LV seed = nearest LV-endocardial point to the centroid of basal-septal AHA
  segments `{2, 3}`, `uvc_transmural ≤ 0.15`.
- RV seed = nearest RV-endocardial point to the centroid of basal-septal RV
  cells (`uvc_intraventricular > 0`, `|aha_angle| < 0.6`, `longitudinal > 0.6`).
- `hisBundleSeed` = midpoint of LV/RV seeds.
- `lineEnd` per ventricle = `seed + unit apical surface tangent`
  (`-grad(uvc_longitudinal)` projected onto the local surface tangent plane).

Growth numerics are anatomy-independent (validated tuning, reused verbatim
from `bivCase_human_tree_v2` / `bivCase_Pig_Morphometric_Tree` — these two
tutorials use *different* `initLength`/`length`, not just different
N_it/repulsivity, so they are kept fully separate below rather than factored
into a shared block):

- Shared across both scenarios: `growthModel surfaceFollow; l_segment 0.3; branchAngle 0.5;`
- Human (from `bivCase_human_tree_v2`): `initLength 30; length 6; N_it 35;`
  (both ventricles), `repulsivity 0.2` LV / `0.1` RV,
  `terminalSelectionModel allLeaves; terminalModel endocardial;` for both
  ventricles. No `extension` block.
- Pig (from `bivCase_Pig_Morphometric_Tree`): `initLength 6; length 5;`
  (both ventricles), `N_it 24` LV / `27` RV, `repulsivity 0.05` LV / `0.1` RV,
  `transmuralMax 0.4` LV / `0.6` RV, `longitudinalMax 0.9` both,
  `terminalModel transmural` both; LV `terminalSelectionModel weightedField`
  with `terminalCount` set per "matching" rule below; RV stays `allLeaves`
  (weighted RV selection is rejected by the current generator). Extension:
  `type gradientFollow; dMin 0.1; dMax 0.4; stepLen 0.5; maxSteps 100; epiClamp 0.8;`

### "Similar to human, but extended" — matching pig terminal count

Per anatomy: run the human scenario first, read the LV terminal/PMJ count
from the `generatePurkinjeTree` log or `postProcessing/generatePurkinjeTree`
output, then set that anatomy's pig scenario `lv.terminalCount` to the same
number. Same PMJ density as the human tree; transmural placement instead of
plain endocardial is the only structural difference, per the user's intent.

## Known risk

Real-anatomy tree growth is parameter-sensitive: the historical record shows
a bad-seed run degenerating to a 13-terminal stub. After each
`generatePurkinjeTree` run, check terminal/coverage counts and the
`exitedMesh`/`maxStepsExhausted` log counters before calling a case done. If
a tree looks degenerate, iterate seeds/`l_segment`/`N_it` before proceeding
to `1DgraphToFoam`.

## Out of scope

- `constant/electroProperties`, ionic model choice, solver variant selection
  (monodomain/bidomain/eikonal), pseudo-ECG — deferred by the user; "the code
  enters inside the solver the same way" regardless of anatomy.
- Scar/ischemia — explicitly excluded (healthy only).
- Actually running `cardiacFoamEP` — the 4 cases are prepared for transfer
  to a separate workstation, not solved locally.
