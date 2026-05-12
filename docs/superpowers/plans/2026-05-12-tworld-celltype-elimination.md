# TWorld AC_celltype Elimination Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the `CONSTANTS[AC_celltype]` indirection from TWorld, rewriting tissue dispatch to use `int tissueFlag` directly (matching TNNP/ToRORd_dynCl/Gaur/BuenoOrovio), and promoting the one hot-path tissue-dependent scaling to its own per-cell `gnalTissueScale` constant. Closes a latent heterogeneity bug at `TWorld_2024.H:1542`.

**Architecture:** Single C++ refactor across three TWorld files; `NUM_CONSTANTS` is unchanged (one constant removed, one added). Validation is bit-identical-diff against pre-refactor traces of the singleCell tutorial running TWorld at all three supported tissues.

**Tech Stack:** OpenFOAM-v2412 / cardiacFoam, C++ (cellML-generated kernel in headers), bash for build (`Allwmake`) and per-tissue validation runs.

**Spec:** [docs/superpowers/specs/2026-05-12-tworld-celltype-elimination-design.md](docs/superpowers/specs/2026-05-12-tworld-celltype-elimination-design.md)

---

## Files Touched

| File | What changes |
|---|---|
| `src/ionicModels/TWorld/TWorld_2024Names.H` | Remove `AC_celltype,` from `CONSTANTS_INDEX` enum; append `gnalTissueScale,` before `NUM_CONSTANTS`; drop `int tissueFlag,` from `TWorldcomputeVariables` declaration |
| `src/ionicModels/TWorld/TWorld_2024.H` | Remove the `TWorldCellType` helper function; remove `"AC_celltype",` from `TWorldCONSTANTS_NAMES`; append `"gnalTissueScale",`; remove the `CONSTANTS[AC_celltype] = TWorldCellType(tissueFlag);` line; rewrite 9 init-time `(CONSTANTS[AC_celltype] == K.0)` ternaries to use `(tissueFlag == K)` directly; add `CONSTANTS[gnalTissueScale] = (tissueFlag == 1) ? 0.7 : 1.0;`; drop `int tissueFlag,` from `TWorldcomputeVariables` definition signature; replace 1 hot-path branch with branch-free form |
| `src/ionicModels/TWorld/TWorld.C` | Drop `tissue()` from 3 call sites in `TWorldcomputeVariables` calls (lines 143, 172, 245). Keep `tissue()` at line 73 (the `TWorldinitConsts` call) |

**Out of scope (per spec §7):** No clamp-inventory work, no heterogeneity wiring, no SoA shim, no CUDA, no other models, no integrator changes.

---

## Phase 0: Baselines (already complete — context only)

Three baseline traces were captured at HEAD `26e1a50` (= the parent of the design-doc commit `1d9dc8e`) before this plan was written:
- `~/cardiacFoam-baselines-tworld/HEAD.txt` — pinned SHA
- `~/cardiacFoam-baselines-tworld/TWorld_endocardialCells/postProcessing/`
- `~/cardiacFoam-baselines-tworld/TWorld_mCells/postProcessing/`
- `~/cardiacFoam-baselines-tworld/TWorld_epicardialCells/postProcessing/`

These are the comparison reference for Phase 2. Do not delete or overwrite them.

---

## Phase 1: Refactor

Owned by: **implementer agent**.

### Task 1.1: Eliminate AC_celltype, promote gnalTissueScale, drop hot-path tissueFlag

**Files:**
- Modify: `src/ionicModels/TWorld/TWorld_2024Names.H` (lines 627, 725, 733)
- Modify: `src/ionicModels/TWorld/TWorld_2024.H` (lines 626, 813–820, 857, 900, 936, 988, 997, 998, 1012, 1024, 1036, 1042, 1252, 1542; plus one new init line and one new names-array line)
- Modify: `src/ionicModels/TWorld/TWorld.C` (lines 143, 172, 245)

**Critical environment requirement:** Every shell invocation that touches `wmake`, `Allwmake`, `Allrun`, or `cardiacFoam` MUST be prefixed with `source /Volumes/OpenFOAM-v2412/etc/bashrc &&`.

- [ ] **Step 1: Read enum & callers to confirm anchors**

```bash
grep -n "AC_celltype\|NUM_CONSTANTS\|TWorldcompute\|TWorldinit" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TWorld/TWorld_2024Names.H
grep -nE "AC_celltype|TWorldCellType" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TWorld/TWorld_2024.H | head -20
grep -n "::TWorldcomputeVariables\|tissue()" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TWorld/TWorld.C
```

Expected anchor locations:
- `TWorld_2024Names.H:627` — `    AC_celltype,` (the enum entry)
- `TWorld_2024Names.H:725` — `    NUM_CONSTANTS` (the enum end)
- `TWorld_2024Names.H:733` — `TWorldcomputeVariables(...)` declaration with `int tissueFlag`
- `TWorld_2024.H:626` — `    "AC_celltype",` (the names-array entry)
- `TWorld_2024.H:813–820` — `inline Foam::scalar TWorldCellType(...)` helper
- `TWorld_2024.H:857` — `CONSTANTS[AC_celltype] = TWorldCellType(tissueFlag);`
- `TWorld_2024.H:900,936,988,997,998,1012,1024,1036,1042` — 9 init-time branches on `CONSTANTS[AC_celltype]`
- `TWorld_2024.H:1252` — `TWorldcomputeVariables(...)` definition with `int tissueFlag`
- `TWorld_2024.H:1542` — hot-path `AV_GNaL` branch
- `TWorld.C:143,172,245` — three `tissue(),` lines inside `TWorldcomputeVariables` calls
- `TWorld.C:73` — one `tissue(),` line inside the `TWorldinitConsts` call (KEEP THIS)

If the line numbers shift slightly (e.g. because of a small unrelated edit), match by content not by literal line number. The anchors above are the ones to find.

- [ ] **Step 2: Edit `TWorld_2024Names.H` — enum**

Remove the line containing `    AC_celltype,` (around line 627). Then immediately before the `NUM_CONSTANTS` line, insert:

```cpp
    gnalTissueScale,
```

The trailing block should look like:

```cpp
    AC_minCaI,
    AC_steepnessCaI,
    AC_steepnessCaSR,
    AC_tauInact,
    AC_tauInact2,
    gnalTissueScale,
    NUM_CONSTANTS
};
```

- [ ] **Step 3: Edit `TWorld_2024Names.H` — drop tissueFlag from `TWorldcomputeVariables` declaration**

Locate the line:
```cpp
TWorldcomputeVariables(double VOI,double* CONSTANTS,double* RATES,double* STATES,double* ALGEBRAIC,int tissueFlag,bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus);
```

Replace with:
```cpp
TWorldcomputeVariables(double VOI,double* CONSTANTS,double* RATES,double* STATES,double* ALGEBRAIC,bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus);
```

Do NOT touch `TWorldinitConsts` — it keeps `int tissueFlag`.

- [ ] **Step 4: Edit `TWorld_2024.H` — remove names-array entry, add new entry**

Find the line `    "AC_celltype",` in `TWorldCONSTANTS_NAMES` (around line 626). Remove that single line.

Then find the end of the names array, which currently looks like:

```cpp
    "AC_tauInact",
    "AC_tauInact2",
};
```

Change it to:

```cpp
    "AC_tauInact",
    "AC_tauInact2",
    "gnalTissueScale",
};
```

- [ ] **Step 5: Edit `TWorld_2024.H` — remove the `TWorldCellType` helper**

Find this 8-line block (around lines 813–820):

```cpp
inline Foam::scalar TWorldCellType(const int tissueFlag)
{
    return (tissueFlag == 1) ? 1.0
         : (tissueFlag == 2) ? 2.0
         : (tissueFlag == 3) ? 0.0
         : 0.0;
}
```

Delete the entire block (including the trailing blank line after the closing brace if there is one). Preserve any `// * * *` separator lines or `void` declaration that immediately follow.

- [ ] **Step 6: Edit `TWorld_2024.H` — remove the AC_celltype init assignment**

Locate the line (around line 857 inside `TWorldinitConsts`):

```cpp
    CONSTANTS[AC_celltype] = TWorldCellType(tissueFlag);
```

Delete that single line.

- [ ] **Step 7: Edit `TWorld_2024.H` — rewrite 9 init-time branches**

Each branch maps `(CONSTANTS[AC_celltype] == 1.0)` to `(tissueFlag == 1)` and `(CONSTANTS[AC_celltype] == 2.0)` to `(tissueFlag == 2)`. The trailing `else` clause is unchanged in every case (it covers endo, where the original `AC_celltype` was 0.0 and the new condition is `tissueFlag != 1 && tissueFlag != 2`, equivalent to `tissueFlag == 3`).

Apply the following 9 substitutions verbatim:

**Line 900 — `AC_INaCa_celltype_factor`:**

Before:
```cpp
    CONSTANTS[AC_INaCa_celltype_factor] = ((CONSTANTS[AC_celltype] == 1.0) ? 1.1 : ((CONSTANTS[AC_celltype] == 2.0) ? 1.4 : 1.0));
```

After:
```cpp
    CONSTANTS[AC_INaCa_celltype_factor] = ((tissueFlag == 1) ? 1.1 : ((tissueFlag == 2) ? 1.4 : 1.0));
```

**Line 936 — `AC_Vmax_SRCaP`:**

Before:
```cpp
    CONSTANTS[AC_Vmax_SRCaP] = ((CONSTANTS[AC_celltype] == 1.0) ? 1.2 * CONSTANTS[AC_Vmax_SRCaP_b] * CONSTANTS[AC_Jup_multiplier] : CONSTANTS[AC_Vmax_SRCaP_b] * CONSTANTS[AC_Jup_multiplier]);
```

After:
```cpp
    CONSTANTS[AC_Vmax_SRCaP] = ((tissueFlag == 1) ? 1.2 * CONSTANTS[AC_Vmax_SRCaP_b] * CONSTANTS[AC_Jup_multiplier] : CONSTANTS[AC_Vmax_SRCaP_b] * CONSTANTS[AC_Jup_multiplier]);
```

**Line 988 — `AC_IK1_celltype_factor`:**

Before:
```cpp
    CONSTANTS[AC_IK1_celltype_factor] = ((CONSTANTS[AC_celltype] == 1.0) ? 1.1 : ((CONSTANTS[AC_celltype] == 2.0) ? 1.3 : 1.0));
```

After:
```cpp
    CONSTANTS[AC_IK1_celltype_factor] = ((tissueFlag == 1) ? 1.1 : ((tissueFlag == 2) ? 1.3 : 1.0));
```

**Line 997 — `AC_Gto_fast`:**

Before:
```cpp
    CONSTANTS[AC_Gto_fast] = ((CONSTANTS[AC_celltype] == 1.0) ? 0.29856 * CONSTANTS[AC_Itof_multiplier] : ((CONSTANTS[AC_celltype] == 2.0) ? 0.14928 * CONSTANTS[AC_Itof_multiplier] : 0.01276 * CONSTANTS[AC_Itof_multiplier]));
```

After:
```cpp
    CONSTANTS[AC_Gto_fast] = ((tissueFlag == 1) ? 0.29856 * CONSTANTS[AC_Itof_multiplier] : ((tissueFlag == 2) ? 0.14928 * CONSTANTS[AC_Itof_multiplier] : 0.01276 * CONSTANTS[AC_Itof_multiplier]));
```

**Line 998 — `AC_Gto_slow`:**

Before:
```cpp
    CONSTANTS[AC_Gto_slow] = ((CONSTANTS[AC_celltype] == 1.0) ? 0.02036 * CONSTANTS[AC_Itos_multiplier] : ((CONSTANTS[AC_celltype] == 2.0) ? 0.04632 * CONSTANTS[AC_Itos_multiplier] : 0.0721 * CONSTANTS[AC_Itos_multiplier]));
```

After:
```cpp
    CONSTANTS[AC_Gto_slow] = ((tissueFlag == 1) ? 0.02036 * CONSTANTS[AC_Itos_multiplier] : ((tissueFlag == 2) ? 0.04632 * CONSTANTS[AC_Itos_multiplier] : 0.0721 * CONSTANTS[AC_Itos_multiplier]));
```

**Line 1012 — `AC_PCa`:**

Before:
```cpp
    CONSTANTS[AC_PCa] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_PCa_b] * CONSTANTS[AC_ICaLPCa_multiplier] * 1.025 : ((CONSTANTS[AC_celltype] == 2.0) ? CONSTANTS[AC_PCa_b] * CONSTANTS[AC_ICaLPCa_multiplier] * 1.1 : CONSTANTS[AC_PCa_b] * CONSTANTS[AC_ICaLPCa_multiplier]));
```

After:
```cpp
    CONSTANTS[AC_PCa] = ((tissueFlag == 1) ? CONSTANTS[AC_PCa_b] * CONSTANTS[AC_ICaLPCa_multiplier] * 1.025 : ((tissueFlag == 2) ? CONSTANTS[AC_PCa_b] * CONSTANTS[AC_ICaLPCa_multiplier] * 1.1 : CONSTANTS[AC_PCa_b] * CONSTANTS[AC_ICaLPCa_multiplier]));
```

**Line 1024 — `AC_PCa_P`:**

Before:
```cpp
    CONSTANTS[AC_PCa_P] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_PCa_Pb] * 1.025 : ((CONSTANTS[AC_celltype] == 2.0) ? CONSTANTS[AC_PCa_Pb] * 1.1 : CONSTANTS[AC_PCa_Pb]));
```

After:
```cpp
    CONSTANTS[AC_PCa_P] = ((tissueFlag == 1) ? CONSTANTS[AC_PCa_Pb] * 1.025 : ((tissueFlag == 2) ? CONSTANTS[AC_PCa_Pb] * 1.1 : CONSTANTS[AC_PCa_Pb]));
```

**Line 1036 — `AC_IKr_celltype_factor`:**

Before:
```cpp
    CONSTANTS[AC_IKr_celltype_factor] = ((CONSTANTS[AC_celltype] == 1.0) ? 1.25 : ((CONSTANTS[AC_celltype] == 2.0) ? 0.7 : 1.0));
```

After:
```cpp
    CONSTANTS[AC_IKr_celltype_factor] = ((tissueFlag == 1) ? 1.25 : ((tissueFlag == 2) ? 0.7 : 1.0));
```

**Line 1042 — `AC_IKs_celltype_factor`:**

Before:
```cpp
    CONSTANTS[AC_IKs_celltype_factor] = ((CONSTANTS[AC_celltype] == 1.0) ? 1.4 : ((CONSTANTS[AC_celltype] == 2.0) ? 0.5 : 1.0));
```

After:
```cpp
    CONSTANTS[AC_IKs_celltype_factor] = ((tissueFlag == 1) ? 1.4 : ((tissueFlag == 2) ? 0.5 : 1.0));
```

- [ ] **Step 8: Edit `TWorld_2024.H` — add the new `gnalTissueScale` init assignment**

Immediately after the `AC_IKs_celltype_factor` line you just rewrote (around line 1042), insert one new line:

```cpp
    CONSTANTS[gnalTissueScale] = (tissueFlag == 1) ? 0.7 : 1.0;
```

Position rationale: keeps all tissue-dependent CONSTANTS assignments grouped together at the end of the existing tissue-dependence block in `TWorldinitConsts`.

- [ ] **Step 9: Edit `TWorld_2024.H` — drop `int tissueFlag` from `TWorldcomputeVariables` definition**

Locate at line 1252:
```cpp
TWorldcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
```

Replace with:
```cpp
TWorldcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
```

- [ ] **Step 10: Edit `TWorld_2024.H` — replace hot-path branch at line 1542**

Locate:
```cpp
    ALGEBRAIC[AV_GNaL] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_GNaL_b] * CONSTANTS[AC_INaL_multiplier] * (1.0 + ALGEBRAIC[AV_fINaLp]) * 0.7 : CONSTANTS[AC_GNaL_b] * CONSTANTS[AC_INaL_multiplier] * (1.0 + ALGEBRAIC[AV_fINaLp]));
```

Replace with:
```cpp
    // AV_GNaL: epi-only 0.7 factor promoted to CONSTANTS[gnalTissueScale]
    // (1.0 in non-epi). Branch-free; no AC_celltype dependency.
    ALGEBRAIC[AV_GNaL] = CONSTANTS[AC_GNaL_b] * CONSTANTS[AC_INaL_multiplier]
                       * (1.0 + ALGEBRAIC[AV_fINaLp]) * CONSTANTS[gnalTissueScale];
```

- [ ] **Step 11: Edit `TWorld.C` — drop `tissue()` from three `TWorldcomputeVariables` call sites**

Three call sites at lines 143, 172, 245 each have a `            tissue(),` line (with that exact leading whitespace, or close to it — check the actual file). The structure of each call is:

```cpp
        ::TWorldcomputeVariables
        (
            ...,
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        ,
            stimulusProtocol()
        );
```

For each of the three call sites, delete the line containing `tissue(),`. Preserve the rest of the call's formatting verbatim (including the unusual standalone comma between `solveVmWithinODESolver()` and `stimulusProtocol()` that's on its own line — that's pre-existing style across this codebase, intentional or not).

Do NOT touch the `tissue(),` at line 73 — that's inside the `TWorldinitConsts` call, which keeps `tissueFlag`.

- [ ] **Step 12: Verify no remaining `AC_celltype` or `TWorldCellType` references**

```bash
grep -nE "AC_celltype|TWorldCellType" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TWorld/
```

Expected: zero matches. If any remain, the corresponding edit was missed — go back and fix.

- [ ] **Step 13: Verify hot-path no longer references `tissueFlag`**

```bash
grep -nE "tissueFlag" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TWorld/TWorld_2024.H
```

Expected: matches only inside `TWorldinitConsts` body and signature (lines ~823, plus the 9 rewritten ternaries and 1 new `gnalTissueScale` assignment around lines 900–1042). NO matches inside `TWorldcomputeVariables` (lines 1252+).

```bash
grep -nE "tissueFlag" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TWorld/TWorld_2024Names.H
```

Expected: exactly one match — inside `TWorldinitConsts` declaration (around line 730). NOT in the `TWorldcomputeVariables` declaration.

- [ ] **Step 14: Verify TWorld.C call sites no longer pass `tissue()` to `*compute*`**

```bash
grep -nF "tissue()" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TWorld/TWorld.C
```

Expected: exactly one match — line 73 (inside the `TWorldinitConsts` call). No `tissue()` inside any `TWorldcomputeVariables` call.

- [ ] **Step 15: Verify `gnalTissueScale` is properly wired (1 declaration + 1 names entry + 1 init + 1 hot-path use)**

```bash
grep -n "gnalTissueScale" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TWorld/TWorld_2024Names.H /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TWorld/TWorld_2024.H
```

Expected output:
- One match in `TWorld_2024Names.H` (the enum entry between `AC_tauInact2` and `NUM_CONSTANTS`)
- Three matches in `TWorld_2024.H`: the names-array string `"gnalTissueScale"`, the init assignment `CONSTANTS[gnalTissueScale] = ...`, the hot-path use `CONSTANTS[gnalTissueScale]` in the `AV_GNaL` formula

- [ ] **Step 16: Build**

```bash
source /Volumes/OpenFOAM-v2412/etc/bashrc && cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors && ./Allwmake 2>&1 | tee /tmp/build-tworld-celltype.log | tail -10
```

Expected: build completes; final lines show "There were no build errors: enjoy cardiacFoam!" or equivalent. If the build fails, the most likely causes are:
- "use of undeclared identifier `AC_celltype`" → a missed reference. Re-run Step 12's grep and fix.
- "use of undeclared identifier `TWorldCellType`" → the helper removal was incomplete or a reference elsewhere remained.
- "too many arguments to function `TWorldcomputeVariables`" → a call site still has `tissue()`. Re-run Step 14's grep.
- "`tissueFlag` was not declared in this scope" → the hot-path body still references it. Re-run Step 13's grep on the file.
- "`gnalTissueScale` was not declared in this scope" → either the enum entry is missing in `_2024Names.H` or the file isn't including the names header.

- [ ] **Step 17: Smoke test all three tissues**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/singleCellprotocols/singleCell
cp constant/electroProperties constant/electroProperties.orig
source /Volumes/OpenFOAM-v2412/etc/bashrc
for tissue in endocardialCells mCells epicardialCells; do
    python3 -c "
import re
p = 'constant/electroProperties'
s = open(p).read()
s = re.sub(r'(ionicModel\s+)\w+;', r'\1TWorld;',     s)
s = re.sub(r'(tissue\s+)\w+;',     r'\1${tissue};', s)
open(p, 'w').write(s)
"
    ./Allclean 2>/dev/null || true
    CF_SKIP_PLOTS=1 ./Allrun 2>&1 | tail -3
    echo "--- TWorld/${tissue} done ---"
done
mv constant/electroProperties.orig constant/electroProperties
```

Expected: each of the three runs completes (final output line should mention `End` or "Results written to:"). Bit-identical comparison vs baselines is Phase 2.

- [ ] **Step 18: Commit (only the three TWorld files)**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git add src/ionicModels/TWorld/TWorld_2024.H \
        src/ionicModels/TWorld/TWorld_2024Names.H \
        src/ionicModels/TWorld/TWorld.C
git diff --cached --stat
git commit -m "refactor(TWorld): eliminate AC_celltype indirection, drop hot-path tissueFlag

Bring TWorld in line with the rest of the ionic-model family
(TNNP/ToRORd_dynCl/Gaur/BuenoOrovio): tissue dispatch happens via
direct (tissueFlag == K) ternaries inside TWorldinitConsts; the
hot path has zero tissue-related branches.

Removes:
- AC_celltype enum entry from CONSTANTS_INDEX
- 'AC_celltype' string from TWorldCONSTANTS_NAMES
- TWorldCellType helper function
- CONSTANTS[AC_celltype] = TWorldCellType(tissueFlag) assignment
- 9 init-time (CONSTANTS[AC_celltype] == K.0) ternaries (rewritten
  to use tissueFlag directly)
- 1 hot-path (CONSTANTS[AC_celltype] == 1.0) branch on AV_GNaL
- int tissueFlag from TWorldcomputeVariables signature
- tissue() argument from 3 TWorldcomputeVariables call sites in
  TWorld.C

Adds:
- gnalTissueScale enum entry + names string + init assignment
  (epi=0.7, non-epi=1.0); replaces the hot-path branch with
  CONSTANTS[gnalTissueScale] multiplication

Net NUM_CONSTANTS change: 0 (one removed, one added).

Closes a latent heterogeneity bug at TWorld_2024.H:1542 where a
50/50 endo+M blend would mis-dispatch as epi (the blended
AC_celltype value lands at 1.0). With AC_celltype eliminated and
the hot-path branch promoted to a per-cell scaling constant, the
blend is well-defined.

Reproduces existing TWorld single-cell traces bit-identically;
validation in Phase 2 of the plan.

Spec: docs/superpowers/specs/2026-05-12-tworld-celltype-elimination-design.md"
```

Important constraints:
- Do NOT use `git add -A` or `git add .`. Only the three explicitly named files.
- Do NOT touch any file outside `src/ionicModels/TWorld/`.
- Do NOT touch `TWorldinitConsts` signature — it keeps `int tissueFlag`.
- Do NOT touch the `tissue()` call at `TWorld.C:73`.

**Report back (under 350 words):**
- Brief edit summary (paste the `git diff --cached --stat` output)
- The 4 verification command outputs (Steps 12, 13, 14, 15) — pass/fail with offending lines if fail
- Build result (Step 16) — last 5 lines
- Smoke test result per tissue (Step 17) — pass/fail
- Commit SHA
- Status: DONE / DONE_WITH_CONCERNS / NEEDS_CONTEXT / BLOCKED

Note about LSP diagnostics: if you see clang/LSP errors about `'stimulusIO.H' file not found`, `'dictionary.H' file not found`, `'ionicModel.H' file not found`, or `'Foam' undeclared identifier`, those are tooling artifacts from missing OpenFOAM include paths in the LSP, not real compile errors. The actual `Allwmake` build is the source of truth.

---

## Phase 2: Validation against pre-refactor baselines

Owned by: **validation agent** (or done inline by the orchestrator).

### Task 2.1: Bit-identical diff per tissue

**Files:**
- Read: `~/cardiacFoam-baselines-tworld/{TWorld_endocardialCells,TWorld_mCells,TWorld_epicardialCells}/postProcessing/`
- Create: `~/cardiacFoam-postrefactor-tworld/{TWorld_endocardialCells,TWorld_mCells,TWorld_epicardialCells}/`
- Create: `~/cardiacFoam-diffs-tworld/{TWorld_<tissue>}.diff`

- [ ] **Step 1: Confirm we're on the post-refactor commit**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git log --oneline -3
```

Expected: top commit is the Phase 1 commit; second commit is `1d9dc8e` (the design doc); third is `26e1a50` (= the baseline HEAD per `~/cardiacFoam-baselines-tworld/HEAD.txt`).

- [ ] **Step 2: Write the post-refactor runner**

Create `/tmp/run_tworld_postrefactor.sh`:

```bash
#!/bin/bash
set -e
TISSUE=$1
TUT=/Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/singleCellprotocols/singleCell
OUT=~/cardiacFoam-postrefactor-tworld/TWorld_${TISSUE}
mkdir -p "$OUT"
cd "$TUT"
cp constant/electroProperties constant/electroProperties.orig
python3 -c "
import re
p = 'constant/electroProperties'
s = open(p).read()
s = re.sub(r'(ionicModel\s+)\w+;', r'\g<1>TWorld;', s)
s = re.sub(r'(tissue\s+)\w+;',     r'\g<1>${TISSUE};', s)
open(p, 'w').write(s)
"
./Allclean 2>/dev/null || true
CF_SKIP_PLOTS=1 ./Allrun > /dev/null 2>&1
cp -r postProcessing "$OUT/" 2>/dev/null || true
cp log.cardiacFoam   "$OUT/" 2>/dev/null || true
mv constant/electroProperties.orig constant/electroProperties
echo "Post-refactor captured: $OUT"
```

```bash
chmod +x /tmp/run_tworld_postrefactor.sh
mkdir -p ~/cardiacFoam-postrefactor-tworld
```

- [ ] **Step 3: Run all three tissues**

```bash
source /Volumes/OpenFOAM-v2412/etc/bashrc
/tmp/run_tworld_postrefactor.sh endocardialCells
/tmp/run_tworld_postrefactor.sh mCells
/tmp/run_tworld_postrefactor.sh epicardialCells
ls ~/cardiacFoam-postrefactor-tworld/
```

Expected: each invocation prints `Post-refactor captured: ~/cardiacFoam-postrefactor-tworld/TWorld_<tissue>` and exits 0. Final `ls` shows three directories.

- [ ] **Step 4: Bit-identical diff per tissue**

```bash
mkdir -p ~/cardiacFoam-diffs-tworld
for tissue in endocardialCells mCells epicardialCells; do
    pair=TWorld_$tissue
    diff -ur \
        --exclude='log.cardiacFoam' --exclude='*.log' \
        ~/cardiacFoam-baselines-tworld/$pair \
        ~/cardiacFoam-postrefactor-tworld/$pair \
        > ~/cardiacFoam-diffs-tworld/$pair.diff 2>&1 || true
    lines=$(wc -l < ~/cardiacFoam-diffs-tworld/$pair.diff)
    echo "$pair: $lines diff lines"
done
```

Expected: each pair shows `0 diff lines`. If any pair shows >0:
- Small diff (< 50 lines) and entirely numerical at the last digit → likely operator-reordering inside one of the rewritten ternaries (e.g. AC_PCa or AC_Gto_fast had a multi-factor expression, and the rewrite accidentally re-grouped factors). Inspect the offending tissue's first few diff lines to identify which ALGEBRAIC or STATE drifted, then trace back to the init-time CONSTANT that was rewritten differently from the original.
- Large diff or every line different → a math change. Most likely a missed branch (`AC_celltype == 2.0` rewritten as `tissueFlag == 1` instead of `tissueFlag == 2`), or the `gnalTissueScale` was set to the wrong value for the wrong tissue.
- Recovery: the design rationale (spec §4) is that every old branch maps to its new branch with identical operand pairs. A non-zero diff means that mapping was broken somewhere. Re-read the offending init line in both the original (via `git show 26e1a50:src/ionicModels/TWorld/TWorld_2024.H | sed -n '<line>p'`) and the post-refactor file, and identify the difference.

- [ ] **Step 5: Final leftover-grep**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
echo "=== AC_celltype anywhere in TWorld? ==="
grep -nE "AC_celltype|TWorldCellType" src/ionicModels/TWorld/ || echo "(none — pass)"

echo ""
echo "=== tissueFlag past TWorldcomputeVariables start? ==="
cv=$(grep -n "^TWorldcomputeVariables" src/ionicModels/TWorld/TWorld_2024.H | head -1 | cut -d: -f1)
echo "TWorldcomputeVariables starts at line $cv"
grep -nE "tissueFlag" src/ionicModels/TWorld/TWorld_2024.H | awk -F: -v cv=$cv '$1+0 >= cv'
echo "(no output above ^ = pass)"
```

Expected:
- `AC_celltype | TWorldCellType` grep: zero matches (prints `(none — pass)`).
- `tissueFlag` past the `TWorldcomputeVariables` start: zero lines (prints only the `(no output above ^ = pass)` confirmation).

- [ ] **Step 6: Report (do NOT commit anything in this phase)**

Report to the user:
- Phase 1 commit SHA + 1-line message
- Per-tissue diff result (PASS/FAIL with line count)
- Build status (PASS — already verified in Phase 1 Step 16)
- Leftover-grep status (PASS/FAIL with offending lines if FAIL)
- Overall: PASS / FAIL

Phase 2 is purely verification — no code changes, no commits.

---

## Success criteria (from spec §8)

This refactor is complete when **all** of the following hold:

1. ✅ `Allwmake` builds clean ("There were no build errors: enjoy cardiacFoam!")
2. ✅ Bit-identical traces (`0 diff lines`) for all three tissues vs `~/cardiacFoam-baselines-tworld/`
3. ✅ `grep -nE "AC_celltype|TWorldCellType" src/ionicModels/TWorld/` returns zero matches
4. ✅ `grep -nE "tissueFlag" src/ionicModels/TWorld/TWorld_2024.H` matches only inside `TWorldinitConsts`
5. ✅ `git log --oneline` shows exactly one new commit on `no-frontend-minor-errors` for this work, touching only the three TWorld files

---

End of plan.
