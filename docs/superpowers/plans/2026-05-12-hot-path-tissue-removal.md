# Hot-Path Tissue-Flag Removal Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the discrete `int tissueFlag` parameter from the hot-path evaluators of three ionic models (TNNP, ToRORd_dynCl, Gaur) by promoting tissue-dependent values into the per-cell `CONSTANTS` array — preparing the codebase for per-cell GPU dispatch matching the existing BuenoOrovio heterogeneity pattern.

**Architecture:** For each branch on `tissueFlag` inside `*computeRates` / `*computeVariables`, find a unified algebraic form covering both branches (typically by zeroing one term's amplitude in non-active branches), promote the literal numbers to named entries in the model's `CONSTANTS_INDEX` enum, and assign per-tissue values inside the existing `*initConsts` tissue branch. The `*initConsts` function keeps `tissueFlag` because it runs once per tissue class at startup. The hot-path functions lose the parameter entirely.

**Tech Stack:** OpenFOAM-v2412 / cardiacFoam, C++ (cellML-generated kernels in headers), bash for build (`Allwmake`) and per-tutorial regression runs.

**Spec:** [docs/superpowers/specs/2026-05-12-hot-path-tissue-removal-design.md](docs/superpowers/specs/2026-05-12-hot-path-tissue-removal-design.md)

---

## Files Touched

| File | What changes |
|---|---|
| `src/ionicModels/TNNP/TNNP_2004Names.H` | +10 enum entries; drop `int tissueFlag` from 2 declarations |
| `src/ionicModels/TNNP/TNNP_2004.H` | +10 per-tissue assignments in `TNNPinitConsts`; replace 4 hot-path branches; drop `int tissueFlag` from 2 signatures |
| `src/ionicModels/TNNP/TNNP.C` | drop `tissue()` from 4 call sites |
| `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023Names.H` | +4 enum entries; drop `int tissueFlag` from 1 declaration |
| `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H` | +4 per-tissue assignments in `ToRORd_dynClinitConsts`; replace 3 hot-path branches; drop `int tissueFlag` from 1 signature |
| `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl.C` | drop `tissue()` from 3 call sites |
| `src/ionicModels/Gaur/Gaur_2021Names.H` | drop `int tissueFlag` from 1 declaration |
| `src/ionicModels/Gaur/Gaur_2021.H` | drop `int tissueFlag` from 2 signatures (`GaurinitConsts` keeps it) |
| `src/ionicModels/Gaur/Gaur.C` | drop `tissue()` from 3 call sites |

**Out of scope (explicit anti-goals from spec §9):** No heterogeneity wiring, no SoA shim files, no CUDA scaffolding, no BuenoOrovio changes, no integrator changes, no clamp implementation, no profiling, no GPU build.

---

## Task ordering rationale

1. **Phase 0** captures pre-refactor traces — must run before any code change so we have a comparison baseline.
2. **Phase 1 (Gaur)** is signature-only: warm-up for the implementer/reviewer agents.
3. **Phase 2 (TNNP)** is the main payload — 3 sub-tasks per the spec's section-boundary discipline.
4. **Phase 3 (ToRORd_dynCl)** repeats Phase 2's pattern with smaller scope.
5. **Phase 4** is final integrated validation across all 7 `(model, tissue)` pairs.
6. **Phase 5** is the parallel clamp-inventory report (fully independent — can start any time after Phase 0).

---

## Phase 0: Capture pre-refactor baseline traces

Owned by: **validation agent**, in its own worktree.

### Task 0.1: Capture seven baseline traces from current `main`

**Files:**
- Read: `tutorials/singleCellprotocols/singleCell/constant/electroProperties`
- Read: `tutorials/singleCellprotocols/singleCell/Allrun`
- Create (outside repo): `~/cardiacFoam-baselines/<model>_<tissue>/` for each of seven pairs

**Pairs to run:**
| Model | Tissue |
|---|---|
| `TNNP` | `endocardialCells` |
| `TNNP` | `mCells` |
| `TNNP` | `epicardialCells` |
| `ToRORd_dynCl` | `endocardialCells` |
| `ToRORd_dynCl` | `mCells` |
| `ToRORd_dynCl` | `epicardialCells` |
| `Gaur` | `myocyte` |

- [ ] **Step 1: Confirm clean working tree on `main` (or designated baseline ref)**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git status --short
git rev-parse HEAD
```
Expected: working tree may have unrelated untracked files but no staged or modified tracked files in `src/ionicModels/{TNNP,ToRORd_dynCl,Gaur}` or `tutorials/singleCellprotocols/singleCell/`. Record HEAD SHA — write it to a note like `~/cardiacFoam-baselines/HEAD.txt`.

- [ ] **Step 2: Build the toolchain on the baseline**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
./Allwmake 2>&1 | tee ~/cardiacFoam-baselines/build-baseline.log
```
Expected: build completes without errors. If it fails, stop and report — we cannot baseline against an unbuildable tree.

- [ ] **Step 3: Create baseline storage**

```bash
mkdir -p ~/cardiacFoam-baselines
```

- [ ] **Step 4: Write the per-pair baseline runner script**

Create `/tmp/run_baseline.sh`:

```bash
#!/bin/bash
# Args: <model> <tissue>
set -e
MODEL=$1
TISSUE=$2
TUT=/Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/singleCellprotocols/singleCell
OUT=~/cardiacFoam-baselines/${MODEL}_${TISSUE}
mkdir -p "$OUT"

cd "$TUT"

# Snapshot original electroProperties so we can restore it
cp constant/electroProperties constant/electroProperties.orig

# Patch ionicModel and tissue keys in the singleCellSolverCoeffs block
python3 -c "
import re
p = 'constant/electroProperties'
s = open(p).read()
s = re.sub(r'(ionicModel\s+)\w+;', r'\g<1>${MODEL};', s)
s = re.sub(r'(tissue\s+)\w+;',     r'\g<1>${TISSUE};', s)
open(p, 'w').write(s)
"

# Clean previous outputs and run
./Allclean 2>/dev/null || true
CF_SKIP_PLOTS=1 ./Allrun

# Capture all top-level numerical outputs (the tutorial dumps Vm and ionic
# state traces under postProcessing or top-level time directories)
cp -r postProcessing "$OUT/" 2>/dev/null || true
cp log.cardiacFoam   "$OUT/" 2>/dev/null || true
for d in [0-9]*; do
    [ -d "$d" ] && cp -r "$d" "$OUT/" 2>/dev/null || true
done

# Restore original electroProperties
mv constant/electroProperties.orig constant/electroProperties
echo "Baseline captured: $OUT"
```

```bash
chmod +x /tmp/run_baseline.sh
```

- [ ] **Step 5: Run the baseline script for all seven pairs**

```bash
/tmp/run_baseline.sh TNNP         endocardialCells
/tmp/run_baseline.sh TNNP         mCells
/tmp/run_baseline.sh TNNP         epicardialCells
/tmp/run_baseline.sh ToRORd_dynCl endocardialCells
/tmp/run_baseline.sh ToRORd_dynCl mCells
/tmp/run_baseline.sh ToRORd_dynCl epicardialCells
/tmp/run_baseline.sh Gaur         myocyte
```
Expected: each run prints `Baseline captured: ~/cardiacFoam-baselines/<model>_<tissue>` and exits 0. If any run fails (e.g. `Gaur` does not support `mCells` — that's why it's only run with `myocyte`), record the failure but do not block — that pair is excluded from validation.

- [ ] **Step 6: Verify baseline corpus exists**

```bash
ls ~/cardiacFoam-baselines/
```
Expected output: 7 directories, one per pair. Each contains at minimum a `log.cardiacFoam` file plus time-step directories or `postProcessing/`.

- [ ] **Step 7: Report**

Report to the orchestrator: HEAD SHA recorded, build OK, all 7 baselines captured. Include the directory listing in the report. **Do not commit anything to the repo at this stage.** Implementation is allowed to begin.

---

## Phase 1: Gaur (signature-only refactor)

Owned by: **implementer agent** (one task), then **reviewer agent**.

### Task 1.1: Drop `tissueFlag` from Gaur hot-path signatures and callers

**Files:**
- Modify: `src/ionicModels/Gaur/Gaur_2021Names.H` (forward declaration)
- Modify: `src/ionicModels/Gaur/Gaur_2021.H` (function definition)
- Modify: `src/ionicModels/Gaur/Gaur.C:134, 167, 241` (three call sites)

- [ ] **Step 1: Read current state to confirm location of edits**

```bash
grep -n "GaurcomputeVariables\|tissueFlag" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/Gaur/Gaur_2021.H /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/Gaur/Gaur_2021Names.H
```
Expected: shows `GaurcomputeVariables(...)` declaration in `Names.H` and definition in `Gaur_2021.H`, both with `int tissueFlag` parameter; `GaurinitConsts(...)` separately retains its `tissueFlag` parameter.

- [ ] **Step 2: Edit `Gaur_2021Names.H`**

Find the line declaring `GaurcomputeVariables` and remove the `int tissueFlag,` parameter only. The exact substring to remove from that one line is `int tissueFlag, ` (with trailing space and comma). Do not touch the `GaurinitConsts` declaration.

- [ ] **Step 3: Edit `Gaur_2021.H`**

Find the function definition `GaurcomputeVariables(...)` (one occurrence). Remove the `int tissueFlag,` parameter from its signature. Do not touch the function body — `GaurcomputeVariables` does not reference `tissueFlag` anywhere in its body, so no body changes are needed. Do not touch `GaurinitConsts`.

- [ ] **Step 4: Edit `Gaur.C` — all three call sites**

In `src/ionicModels/Gaur/Gaur.C`, three call sites currently look like:

```cpp
::GaurcomputeVariables
(
    /* time/VOI */,
    CONSTANTS_.data(),
    /* ... */,
    ALGEBRAICI.data(),
    tissue(),
    solveVmWithinODESolver()
,
    stimulusProtocol()
);
```

For each of the three call sites (around lines 134, 167, 241), remove the line `        tissue(),` (or however it is formatted in the local context — drop the entire `tissue()` argument and its trailing comma). Verify the call still has the correct number of arguments matching the new signature.

- [ ] **Step 5: Verify with grep**

```bash
grep -n "GaurcomputeVariables" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/Gaur/Gaur.C
grep -n "tissue\(\)" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/Gaur/Gaur.C
```
Expected: three `::GaurcomputeVariables` matches; only `tissue()` matches that remain are inside `GaurinitConsts` calls or other unrelated context (NOT inside any `GaurcomputeVariables` call).

- [ ] **Step 6: Build**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
./Allwmake 2>&1 | tee /tmp/build-gaur.log | tail -30
```
Expected: build completes without errors. If the build fails, fix the call site count or signature mismatch and rebuild before continuing.

- [ ] **Step 7: Quick smoke test**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/singleCellprotocols/singleCell
cp constant/electroProperties constant/electroProperties.orig
python3 -c "
import re
p = 'constant/electroProperties'
s = open(p).read()
s = re.sub(r'(ionicModel\s+)\w+;', r'\1Gaur;',     s)
s = re.sub(r'(tissue\s+)\w+;',     r'\1myocyte;', s)
open(p, 'w').write(s)
"
./Allclean 2>/dev/null || true
CF_SKIP_PLOTS=1 ./Allrun 2>&1 | tail -5
mv constant/electroProperties.orig constant/electroProperties
```
Expected: solver completes (no FPE, no crash); last lines of output mention `End` or similar OpenFOAM completion message. Bit-identical comparison happens in Phase 4.

- [ ] **Step 8: Commit**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git add src/ionicModels/Gaur/Gaur_2021.H \
        src/ionicModels/Gaur/Gaur_2021Names.H \
        src/ionicModels/Gaur/Gaur.C
git commit -m "refactor(Gaur): drop tissueFlag from GaurcomputeVariables signature

The Gaur model never branches on tissueFlag in its hot path —
the parameter was carried in the signature but unused. Drop it
from GaurcomputeVariables and the three TNNP-style call sites in
Gaur.C. GaurinitConsts retains tissueFlag (still used at startup
for tissue-class constant baking).

Part of the hot-path tissue-flag removal refactor (spec
docs/superpowers/specs/2026-05-12-hot-path-tissue-removal-design.md)."
```

### Task 1.2: Reviewer verification of Phase 1

Owned by: **reviewer agent**.

- [ ] **Step 1: Verify only the expected files changed**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git diff HEAD~1 --stat
```
Expected: exactly three files (`Gaur_2021.H`, `Gaur_2021Names.H`, `Gaur.C`), no others.

- [ ] **Step 2: Verify the diff has no unintended content**

```bash
git diff HEAD~1 src/ionicModels/Gaur/
```
Verify: only deletions of `int tissueFlag,` (or `tissue(),`); no other lines changed; `GaurinitConsts` still has its `tissueFlag` parameter; no body code changed.

- [ ] **Step 3: Confirm build still passes**

```bash
./Allwmake 2>&1 | tail -5
```
Expected: success messages.

- [ ] **Step 4: Confirm leftover-grep is correct for Gaur**

```bash
grep -nE "tissueFlag" src/ionicModels/Gaur/
```
Expected: matches exist only inside `GaurinitConsts` declaration/definition. Specifically:
- `Gaur_2021.H:460` (`GaurinitConsts` definition signature)
- `Gaur_2021Names.H` (one line — `GaurinitConsts` declaration)
No other matches.

- [ ] **Step 5: Sign off or send back to implementer**

If all four checks pass, report PASS to the orchestrator. If any fails, report the specific check that failed with the offending lines, and the implementer agent re-runs Task 1.1 with corrections.

---

## Phase 2: TNNP (constants promotion + hot-path refactor)

Owned by: **implementer agent** (three tasks, with reviewer between each), then **reviewer agent**.

### Task 2.1 [Section a]: Add 10 enum entries; drop `tissueFlag` from hot-path declarations

**Files:**
- Modify: `src/ionicModels/TNNP/TNNP_2004Names.H`

- [ ] **Step 1: Read the existing enum to identify the insertion point**

```bash
grep -n "tau_fCa\|NUM_CONSTANTS" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TNNP/TNNP_2004Names.H
```
Expected: `tau_fCa` at line 66 (last entry before `NUM_CONSTANTS`); `NUM_CONSTANTS` at line 68.

- [ ] **Step 2: Insert 10 new enum entries between `tau_fCa,` and `NUM_CONSTANTS`**

Modify `src/ionicModels/TNNP/TNNP_2004Names.H` line 66-68 from:

```cpp
    tau_fCa,

    NUM_CONSTANTS};
```

to:

```cpp
    tau_fCa,

    // Tissue-dependent algebraic parameters (promoted from hot-path branches
    // in TNNPcomputeRates / TNNPcomputeVariables — see spec
    // docs/superpowers/specs/2026-05-12-hot-path-tissue-removal-design.md §3.1)
    sInfPrefactor,
    sInfShift,
    sInfScale,
    tauSGaussAmp,
    tauSGaussShift,
    tauSGaussWidth,
    tauSSigmoidAmp,
    tauSSigmoidShift,
    tauSSigmoidScale,
    tauSOffset,

    NUM_CONSTANTS};
```

- [ ] **Step 3: Drop `int tissueFlag` from `TNNPcomputeRates` and `TNNPcomputeVariables` declarations**

In `src/ionicModels/TNNP/TNNP_2004Names.H`, the file currently has these two declarations around lines 173 and 175:

```cpp
TNNPcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus);

TNNPcomputeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus);
```

Change each to remove the `int tissueFlag,` parameter:

```cpp
TNNPcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus);

TNNPcomputeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus);
```

Do **not** touch the `TNNPinitConsts` declaration on line 171 — it keeps `tissueFlag`.

- [ ] **Step 4: Verify with grep**

```bash
grep -nE "tissueFlag|NUM_CONSTANTS" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TNNP/TNNP_2004Names.H
```
Expected: `tissueFlag` only inside the `TNNPinitConsts` declaration (one match). `NUM_CONSTANTS` reachable via the enum closing brace.

- [ ] **Step 5: Do NOT build or commit yet**

The function definitions in `TNNP_2004.H` still have `tissueFlag` — building now will fail. Continue to Task 2.2.

### Task 2.2 [Section b]: Add 10 per-tissue assignments in `TNNPinitConsts`

**Files:**
- Modify: `src/ionicModels/TNNP/TNNP_2004.H` (inside `TNNPinitConsts`, after line 354)

- [ ] **Step 1: Read the existing init code to identify insertion point**

```bash
sed -n '345,360p' /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TNNP/TNNP_2004.H
```
Expected: shows the existing tissue branches for `g_to` ending at line 354, followed by `STATES[s] = 1;` at line 356.

- [ ] **Step 2: Insert 10 per-tissue assignments after the existing `g_to` branch**

In `src/ionicModels/TNNP/TNNP_2004.H`, locate the lines that look like:

```cpp
//Conditional tissue conductance Gto
CONSTANTS[g_to] = (tissueFlag == 3)
          ? 0.073
          : 0.294;

STATES[s] = 1;
```

Insert this block immediately AFTER the `CONSTANTS[g_to] = ...;` statement (3 lines, ending at the semicolon) and BEFORE `STATES[s] = 1;`:

```cpp

// Tissue-dependent algebraic parameters (promoted from hot-path branches
// in TNNPcomputeRates / TNNPcomputeVariables — see spec
// docs/superpowers/specs/2026-05-12-hot-path-tissue-removal-design.md §3.1).
// Original branch: only tissueFlag == 1 (endo) selects the first set;
// M-cells and epi share the second set.
{
    const bool isEndo = (tissueFlag == 1);

    CONSTANTS[sInfPrefactor]    = isEndo ?    1.00000 :    1.10000;
    CONSTANTS[sInfShift]        = isEndo ?   20.00000 :   28.00000;
    CONSTANTS[sInfScale]        = isEndo ?    5.00000 :    6.00000;

    CONSTANTS[tauSGaussAmp]     = isEndo ?   85.00000 : 1000.00000;
    CONSTANTS[tauSGaussShift]   = isEndo ?   45.00000 :   67.00000;
    CONSTANTS[tauSGaussWidth]   = isEndo ?  320.00000 : 1000.00000;

    // Sigmoid term: amp == 0 in non-endo collapses the whole term;
    // scale must remain non-zero to avoid 1/0 in IEEE math even when
    // amp is zero (the result is multiplied by 0, but only if the
    // intermediate doesn't blow up to NaN first).
    CONSTANTS[tauSSigmoidAmp]   = isEndo ?    5.00000 :    0.00000;
    CONSTANTS[tauSSigmoidShift] = isEndo ?  -20.00000 :    0.00000;
    CONSTANTS[tauSSigmoidScale] = isEndo ?    5.00000 :    1.00000;

    CONSTANTS[tauSOffset]       = isEndo ?    3.00000 :    8.00000;
}

```

- [ ] **Step 3: Verify with grep**

```bash
grep -n "sInfPrefactor\|tauSGaussAmp\|tauSOffset" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TNNP/TNNP_2004.H
```
Expected: each constant appears at least twice — once in the new init block (this task) and once in the future hot-path edits (Task 2.3 will add those references). Right now only the init-block references exist.

- [ ] **Step 4: Do NOT build or commit yet**

Hot-path branches still reference `tissueFlag`. Continue to Task 2.3.

### Task 2.3 [Section c]: Replace 4 hot-path branches with unified formulas; drop `tissueFlag` from signatures

**Files:**
- Modify: `src/ionicModels/TNNP/TNNP_2004.H` lines 386, 426-431, 516, 523-529

- [ ] **Step 1: Drop `int tissueFlag` from `TNNPcomputeRates` definition signature (line 386)**

Locate this line in `src/ionicModels/TNNP/TNNP_2004.H`:

```cpp
TNNPcomputeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
```

Change to:

```cpp
TNNPcomputeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
```

- [ ] **Step 2: Replace the `s_inf` and `tau_s` branches inside `TNNPcomputeRates` (lines 426-431)**

Locate this block:

```cpp
//Conditional tissue flag for the Ito current
ALGEBRAIC[s_inf] = (tissueFlag == 1)
    ? 1.00000 / (1.00000 + exp((STATES[V] + 20.0000) / 5.00000))
    : 1.10000 / (1.00000 + exp((STATES[V] + 28.0000) / 6.00000));
ALGEBRAIC[tau_s] = (tissueFlag == 1)
     ? 85.0000*exp(- pow(STATES[V]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[V] - 20.0000)/5.00000))+3.00000
     : 1000.0000*exp(- pow(STATES[V]+67.0000, 2.00000)/1000.000)+8.00000;
```

Replace with:

```cpp
// Ito (tissue dependence now lives in CONSTANTS[sInf*] and CONSTANTS[tauS*],
// set per tissue in TNNPinitConsts — see spec §3.1).
ALGEBRAIC[s_inf] = CONSTANTS[sInfPrefactor]
    / (1.00000 + exp((STATES[V] + CONSTANTS[sInfShift]) / CONSTANTS[sInfScale]));
ALGEBRAIC[tau_s] =
      CONSTANTS[tauSGaussAmp]
      * exp(- pow(STATES[V] + CONSTANTS[tauSGaussShift], 2.00000)
            / CONSTANTS[tauSGaussWidth])
    + CONSTANTS[tauSSigmoidAmp]
      / (1.00000 + exp((STATES[V] + CONSTANTS[tauSSigmoidShift])
                       / CONSTANTS[tauSSigmoidScale]))
    + CONSTANTS[tauSOffset];
```

**Critical:** preserve the operator structure of the original endo branch verbatim inside the new formula — same `pow(.., 2.00000)/width` form, same `1.00000+exp(...)` denominators, same number formatting (`5.00000`, etc.). The non-endo branch's `tauSSigmoidAmp = 0` collapses the second term so the math reduces exactly to the original non-endo formula.

- [ ] **Step 3: Drop `int tissueFlag` from `TNNPcomputeVariables` definition signature (line 516)**

Locate:

```cpp
TNNPcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
```

Change to:

```cpp
TNNPcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
```

- [ ] **Step 4: Replace the `s_inf` and `tau_s` branches inside `TNNPcomputeVariables` (lines 523-529)**

Locate this block (note: it has a slightly different leading-whitespace pattern than the Rates block — preserve column alignment of the original):

```cpp
//Conditional tissue flag for the Ito current
ALGEBRAIC[s_inf] = (tissueFlag == 1)
    ? 1.00000 / (1.00000 + exp((STATES[V] + 20.0000) / 5.00000))
    : 1.10000 / (1.00000 + exp((STATES[V] + 28.0000) / 6.00000));

ALGEBRAIC[tau_s] = (tissueFlag == 1)
     ? 85.0000*exp(- pow(STATES[V]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[V] - 20.0000)/5.00000))+3.00000
     : 1000.0000*exp(- pow(STATES[V]+67.0000, 2.00000)/1000.000)+8.00000;
```

Replace with the same unified formulas as Step 2 (identical text — copy-paste from Step 2's replacement):

```cpp
// Ito (tissue dependence now lives in CONSTANTS[sInf*] and CONSTANTS[tauS*],
// set per tissue in TNNPinitConsts — see spec §3.1).
ALGEBRAIC[s_inf] = CONSTANTS[sInfPrefactor]
    / (1.00000 + exp((STATES[V] + CONSTANTS[sInfShift]) / CONSTANTS[sInfScale]));
ALGEBRAIC[tau_s] =
      CONSTANTS[tauSGaussAmp]
      * exp(- pow(STATES[V] + CONSTANTS[tauSGaussShift], 2.00000)
            / CONSTANTS[tauSGaussWidth])
    + CONSTANTS[tauSSigmoidAmp]
      / (1.00000 + exp((STATES[V] + CONSTANTS[tauSSigmoidShift])
                       / CONSTANTS[tauSSigmoidScale]))
    + CONSTANTS[tauSOffset];
```

- [ ] **Step 5: Drop `tissue()` from all four call sites in `TNNP.C`**

The four call sites are around lines 146 (`TNNPcomputeVariables`), 158 (`TNNPcomputeRates`), 191 (`TNNPcomputeRates`), 250 (`TNNPcomputeVariables`). Each currently looks like:

```cpp
::TNNPcomputeVariables
(
    tEnd,
    CONSTANTS_.data(),
    RATESI.data(),
    STATESI.data(),
    ALGEBRAICI.data(),
    tissue(),
    solveVmWithinODESolver()
,
    stimulusProtocol()
);
```

For each of the four call sites, delete the line `    tissue(),` (preserve the rest of the call's formatting verbatim). After the edit, the call should look like:

```cpp
::TNNPcomputeVariables
(
    tEnd,
    CONSTANTS_.data(),
    RATESI.data(),
    STATESI.data(),
    ALGEBRAICI.data(),
    solveVmWithinODESolver()
,
    stimulusProtocol()
);
```

Do this for all four call sites (the `tEnd` argument and method name vary — `TNNPcomputeVariables` vs `TNNPcomputeRates`, `tEnd` vs `0.0` vs `t` — but the `tissue(),` line removal is the same in all four).

- [ ] **Step 6: Verify with grep — no remaining hot-path `tissueFlag`**

```bash
grep -nE "tissueFlag" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TNNP/TNNP_2004.H /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TNNP/TNNP_2004Names.H
```
Expected: matches only inside `TNNPinitConsts` (declaration in Names.H, definition signature line 310 plus the existing `g_Ks` / `g_to` / new constants ternaries in `TNNP_2004.H`). No matches in `TNNPcomputeRates` or `TNNPcomputeVariables` bodies or signatures.

```bash
grep -n "tissue\(\)" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/TNNP/TNNP.C
```
Expected: matches only inside `TNNPinitConsts` calls (around line 70 in `TNNP.C`); no `tissue()` inside any `TNNPcompute*` call.

- [ ] **Step 7: Build**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
./Allwmake 2>&1 | tee /tmp/build-tnnp.log | tail -30
```
Expected: build completes without errors. Common failure modes if it doesn't:
- "too many arguments to function `TNNPcomputeVariables`" → a call site still has `tissue()`.
- "`tissueFlag` was not declared in this scope" → a hot-path body still references it (you missed one of the 4 occurrences).
- "expected `;` before ..." → bracket/parenthesis mismatch from the multi-line edits.

- [ ] **Step 8: Quick smoke test (one tissue per pair, just to confirm the solver runs)**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/singleCellprotocols/singleCell
cp constant/electroProperties constant/electroProperties.orig
for tissue in endocardialCells mCells epicardialCells; do
    python3 -c "
import re
p = 'constant/electroProperties'
s = open(p).read()
s = re.sub(r'(ionicModel\s+)\w+;', r'\1TNNP;',     s)
s = re.sub(r'(tissue\s+)\w+;',     r'\1${tissue};', s)
open(p, 'w').write(s)
"
    ./Allclean 2>/dev/null || true
    CF_SKIP_PLOTS=1 ./Allrun 2>&1 | tail -3
    echo "--- TNNP/${tissue} done ---"
done
mv constant/electroProperties.orig constant/electroProperties
```
Expected: each of the three runs completes (`End` appearing in the OpenFOAM log). Bit-identical comparison is in Phase 4.

- [ ] **Step 9: Commit**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git add src/ionicModels/TNNP/TNNP_2004.H \
        src/ionicModels/TNNP/TNNP_2004Names.H \
        src/ionicModels/TNNP/TNNP.C
git commit -m "refactor(TNNP): promote tissue-dependent algebraics to CONSTANTS

Add 10 new entries to TNNP CONSTANTS_INDEX (sInfPrefactor/Shift/Scale,
tauSGaussAmp/Shift/Width, tauSSigmoidAmp/Shift/Scale, tauSOffset),
populated per tissue in TNNPinitConsts. Replace the four hot-path
ternary branches on tissueFlag (s_inf and tau_s, in both
TNNPcomputeRates and TNNPcomputeVariables) with unified formulas
that read those CONSTANTS — non-endo's zeroed sigmoid amplitude
collapses the extra term back to the original non-endo math.

Drop int tissueFlag from TNNPcomputeRates and TNNPcomputeVariables
signatures and from the four call sites in TNNP.C. TNNPinitConsts
keeps tissueFlag (still used at startup for tissue-class baking).

NUM_CONSTANTS goes from 42 to 52. Reproduces existing single-cell
Vm traces bit-identically (validation in Phase 4 of the plan).

Part of the hot-path tissue-flag removal refactor (spec
docs/superpowers/specs/2026-05-12-hot-path-tissue-removal-design.md)."
```

### Task 2.4: Reviewer verification of Phase 2

Owned by: **reviewer agent**.

- [ ] **Step 1: Diff against pre-Phase-2 commit**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git diff HEAD~1 --stat
```
Expected: exactly three files (`TNNP_2004.H`, `TNNP_2004Names.H`, `TNNP.C`). No other files.

- [ ] **Step 2: Verify enum has all 10 new entries in the right position**

```bash
sed -n '24,80p' src/ionicModels/TNNP/TNNP_2004Names.H
```
Verify: 10 new entries (`sInfPrefactor` through `tauSOffset`) appear after `tau_fCa,` and before `NUM_CONSTANTS};`.

- [ ] **Step 3: Verify each new constant is set inside `TNNPinitConsts`**

```bash
for c in sInfPrefactor sInfShift sInfScale tauSGaussAmp tauSGaussShift tauSGaussWidth tauSSigmoidAmp tauSSigmoidShift tauSSigmoidScale tauSOffset; do
  count=$(grep -c "CONSTANTS\[$c\]\s*=" src/ionicModels/TNNP/TNNP_2004.H)
  echo "$c: assigned $count times in TNNP_2004.H (expected: 1)"
done
```
Expected: every constant assigned exactly once (inside the `TNNPinitConsts` body). If any shows 0, the init-block assignment is missing for that constant — fail.

- [ ] **Step 4: Verify each new constant is read in the hot-path bodies**

```bash
for c in sInfPrefactor sInfShift sInfScale tauSGaussAmp tauSGaussShift tauSGaussWidth tauSSigmoidAmp tauSSigmoidShift tauSSigmoidScale tauSOffset; do
  count=$(grep -c "CONSTANTS\[$c\]" src/ionicModels/TNNP/TNNP_2004.H)
  echo "$c: total references $count in TNNP_2004.H (expected: 3 — 1 init + 2 hot-path)"
done
```
Expected: every constant appears 3 times total (1 in init, 1 in `TNNPcomputeRates` body, 1 in `TNNPcomputeVariables` body). Anything <3 means a hot-path use is missing.

- [ ] **Step 5: Verify hot-path tissue branches are gone**

```bash
grep -nE "tissueFlag\s*[=!,)]" src/ionicModels/TNNP/TNNP_2004.H
```
Expected: matches only at:
- The `TNNPinitConsts` signature (line 310 region)
- Inside `TNNPinitConsts` body (the `g_Ks`, `g_to`, and new 10 ternaries)
No matches inside `TNNPcomputeRates` or `TNNPcomputeVariables` bodies (lines 386–522 for Rates, 516+ for Variables).

- [ ] **Step 6: Verify build still passes**

```bash
./Allwmake 2>&1 | tail -5
```
Expected: success.

- [ ] **Step 7: Sign off or send back**

If steps 1–6 all pass: report PASS. Otherwise, report the specific failed check with offending lines, and the implementer agent re-runs the relevant sub-task.

---

## Phase 3: ToRORd_dynCl (constants promotion + hot-path refactor)

Owned by: **implementer agent** (three tasks, with reviewer between each), then **reviewer agent**.

### Task 3.1 [Section a]: Add 4 enum entries; drop `tissueFlag` from hot-path declaration

**Files:**
- Modify: `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023Names.H`

- [ ] **Step 1: Find `NUM_CONSTANTS` and the last existing enum entry**

```bash
grep -n "NUM_CONSTANTS\|enum CONSTANTS" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023Names.H
```
Expected: shows the `enum CONSTANTS_INDEX { ... NUM_CONSTANTS };` block. Note the last entry before `NUM_CONSTANTS`.

- [ ] **Step 2: Insert 4 new entries before `NUM_CONSTANTS`**

In `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023Names.H`, locate `NUM_CONSTANTS` inside the `CONSTANTS_INDEX` enum and insert immediately before it:

```cpp

    // Tissue-dependent algebraic parameters (promoted from hot-path branches
    // in ToRORd_dynClcomputeVariables — see spec
    // docs/superpowers/specs/2026-05-12-hot-path-tissue-removal-design.md §3.2)
    deltaEpiAmp,
    deltaEpiShift,
    deltaEpiScale,
    jrelTissueScale,

```

- [ ] **Step 3: Drop `int tissueFlag` from `ToRORd_dynClcomputeVariables` declaration**

Find the forward declaration of `ToRORd_dynClcomputeVariables` in `ToRORd_dynCl_2023Names.H` and remove the `int tissueFlag,` parameter (with comma and trailing space). Do not touch `ToRORd_dynClinitConsts`.

- [ ] **Step 4: Verify**

```bash
grep -nE "tissueFlag|deltaEpi|jrelTissueScale" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023Names.H
```
Expected: 4 new constants present in the enum; `tissueFlag` only inside `ToRORd_dynClinitConsts` declaration.

- [ ] **Step 5: Do NOT build or commit yet — definitions still mismatch**

### Task 3.2 [Section b]: Add 4 per-tissue assignments in `ToRORd_dynClinitConsts`

**Files:**
- Modify: `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H` (inside `ToRORd_dynClinitConsts`, after the existing per-tissue scaling block)

- [ ] **Step 1: Find a good insertion point inside `ToRORd_dynClinitConsts`**

```bash
sed -n '760,795p' /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H
```
Expected: shows `CONSTANTS[AC_cmdnmax] = ...` (line 769) followed by the `tissueSelect`-driven `STATES[...]` initializations (line 782+). Insert the new block between `AC_cmdnmax` assignment and the first `STATES[V] = tissueSelect(...)` line.

- [ ] **Step 2: Insert 4 per-tissue assignments**

Find the line:

```cpp
    CONSTANTS[AC_cmdnmax] = ((tissueFlag == 1) ? CONSTANTS[AC_cmdnmax_b] * 1.3 : CONSTANTS[AC_cmdnmax_b]);
```

Immediately after it (and before the `STATES[V] = tissueSelect(...)` block), insert:

```cpp

    // Tissue-dependent algebraic parameters (promoted from hot-path branches
    // in ToRORd_dynClcomputeVariables — see spec
    // docs/superpowers/specs/2026-05-12-hot-path-tissue-removal-design.md §3.2).
    {
        const bool isEndo = (tissueFlag == 1);
        const bool isMCell = (tissueFlag == 2);

        // AV_delta_epi: only the endo branch produces a sigmoid term.
        // M-cells and epi reduce to constant 1.0 via amp = 0.
        CONSTANTS[deltaEpiAmp]   = isEndo ? 0.95 : 0.0;
        CONSTANTS[deltaEpiShift] = isEndo ? 70.0 : 0.0;
        CONSTANTS[deltaEpiScale] = isEndo ?  5.0 : 1.0;  // must be non-zero

        // AV_Jrel_inf and AV_Jrel_infp scale by 1.7 only in M-cells.
        CONSTANTS[jrelTissueScale] = isMCell ? 1.7 : 1.0;
    }

```

- [ ] **Step 3: Verify**

```bash
grep -n "deltaEpi\|jrelTissueScale" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H
```
Expected: each of the 4 constants appears at least once in the init block (this task) and will appear again in hot-path bodies (Task 3.3).

- [ ] **Step 4: Do NOT build or commit yet — hot-path branches still reference `tissueFlag`**

### Task 3.3 [Section c]: Replace 3 hot-path branches with unified formulas; drop `tissueFlag` from signature

**Files:**
- Modify: `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H` line 848 (signature), 948 (delta_epi), 1128–1129 (Jrel)
- Modify: `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl.C` lines 134, 167, 241 (call sites)

- [ ] **Step 1: Drop `int tissueFlag` from `ToRORd_dynClcomputeVariables` definition signature**

Locate at line 848:

```cpp
ToRORd_dynClcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
```

Change to:

```cpp
ToRORd_dynClcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
```

- [ ] **Step 2: Replace `AV_delta_epi` branch (line 948)**

Locate:

```cpp
    ALGEBRAIC[AV_delta_epi] = ((tissueFlag == 1) ? 1.0 - 0.95 / (1.0 + exp((STATES[V] + CONSTANTS[AC_EKshift] + 70.0) / 5.0)) : 1.0);
```

Replace with:

```cpp
    // AV_delta_epi: tissue dependence promoted to CONSTANTS[deltaEpi*] —
    // amp == 0 in non-endo collapses the sigmoid term back to constant 1.0.
    ALGEBRAIC[AV_delta_epi] = 1.0 - CONSTANTS[deltaEpiAmp]
        / (1.0 + exp((STATES[V] + CONSTANTS[AC_EKshift] + CONSTANTS[deltaEpiShift])
                     / CONSTANTS[deltaEpiScale]));
```

- [ ] **Step 3: Replace `AV_Jrel_inf` and `AV_Jrel_infp` branches (lines 1128-1129)**

Locate:

```cpp
    ALGEBRAIC[AV_Jrel_inf] = ((tissueFlag == 2.0) ? ALGEBRAIC[AV_Jrel_inf_b] * 1.7 : ALGEBRAIC[AV_Jrel_inf_b]);
    ALGEBRAIC[AV_Jrel_infp] = ((tissueFlag == 2.0) ? ALGEBRAIC[AV_Jrel_infp_b] * 1.7 : ALGEBRAIC[AV_Jrel_infp_b]);
```

Replace with:

```cpp
    // Jrel scalings: M-cell-only 1.7x factor promoted to CONSTANTS[jrelTissueScale].
    ALGEBRAIC[AV_Jrel_inf]  = ALGEBRAIC[AV_Jrel_inf_b]  * CONSTANTS[jrelTissueScale];
    ALGEBRAIC[AV_Jrel_infp] = ALGEBRAIC[AV_Jrel_infp_b] * CONSTANTS[jrelTissueScale];
```

- [ ] **Step 4: Drop `tissue()` from all three call sites in `ToRORd_dynCl.C`**

The three call sites are at lines 134, 167, 241 of `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl.C`. Each has the same structure as the TNNP call sites:

```cpp
::ToRORd_dynClcomputeVariables
(
    /* time */,
    CONSTANTS_.data(),
    /* ... */,
    ALGEBRAICI.data(),
    tissue(),
    solveVmWithinODESolver()
,
    stimulusProtocol()
);
```

For each of the three call sites, delete the line `    tissue(),` (preserve the rest of the call's formatting verbatim).

- [ ] **Step 5: Verify with grep — no remaining hot-path `tissueFlag`**

```bash
grep -nE "tissueFlag" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H | head -10
```
Expected: matches only inside `ToRORd_dynClinitConsts` (signature line 571 region and the existing per-tissue ternaries lines 619–790, plus the new 4 ternaries you just added). NO matches at lines 848 (signature), 948, 1128, or 1129.

```bash
grep -n "tissue\(\)" /Users/simaocastro/noFrontendCardiacFoam_minor_errors/src/ionicModels/ToRORd_dynCl/ToRORd_dynCl.C
```
Expected: matches only inside `ToRORd_dynClinitConsts` calls; none inside `ToRORd_dynClcomputeVariables` calls.

- [ ] **Step 6: Build**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
./Allwmake 2>&1 | tee /tmp/build-torord.log | tail -30
```
Expected: build completes without errors.

- [ ] **Step 7: Quick smoke test**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/singleCellprotocols/singleCell
cp constant/electroProperties constant/electroProperties.orig
for tissue in endocardialCells mCells epicardialCells; do
    python3 -c "
import re
p = 'constant/electroProperties'
s = open(p).read()
s = re.sub(r'(ionicModel\s+)\w+;', r'\1ToRORd_dynCl;', s)
s = re.sub(r'(tissue\s+)\w+;',     r'\1${tissue};',    s)
open(p, 'w').write(s)
"
    ./Allclean 2>/dev/null || true
    CF_SKIP_PLOTS=1 ./Allrun 2>&1 | tail -3
    echo "--- ToRORd_dynCl/${tissue} done ---"
done
mv constant/electroProperties.orig constant/electroProperties
```
Expected: each of the three runs completes without crashing.

- [ ] **Step 8: Commit**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git add src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H \
        src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023Names.H \
        src/ionicModels/ToRORd_dynCl/ToRORd_dynCl.C
git commit -m "refactor(ToRORd_dynCl): promote tissue-dependent algebraics to CONSTANTS

Add 4 new entries to ToRORd_dynCl CONSTANTS_INDEX (deltaEpiAmp/Shift/
Scale, jrelTissueScale), populated per tissue in ToRORd_dynClinitConsts.
Replace the three hot-path ternary branches on tissueFlag (AV_delta_epi
plus AV_Jrel_inf and AV_Jrel_infp, all in ToRORd_dynClcomputeVariables)
with unified formulas that read those CONSTANTS — non-endo's zeroed
deltaEpiAmp collapses the sigmoid back to constant 1.0; jrelTissueScale
of 1.0 in non-M-cells leaves AV_Jrel_inf/infp unscaled.

Drop int tissueFlag from ToRORd_dynClcomputeVariables signature and
from the three call sites in ToRORd_dynCl.C. ToRORd_dynClinitConsts
keeps tissueFlag.

Reproduces existing single-cell Vm traces bit-identically (validation
in Phase 4 of the plan).

Part of the hot-path tissue-flag removal refactor (spec
docs/superpowers/specs/2026-05-12-hot-path-tissue-removal-design.md)."
```

### Task 3.4: Reviewer verification of Phase 3

Owned by: **reviewer agent**.

- [ ] **Step 1: Diff scope check**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git diff HEAD~1 --stat
```
Expected: exactly three files (`ToRORd_dynCl_2023.H`, `ToRORd_dynCl_2023Names.H`, `ToRORd_dynCl.C`).

- [ ] **Step 2: Verify enum has all 4 new entries**

```bash
grep -n "deltaEpi\|jrelTissueScale" src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023Names.H
```
Expected: 4 lines showing the new entries.

- [ ] **Step 3: Verify each new constant is set inside `ToRORd_dynClinitConsts`**

```bash
for c in deltaEpiAmp deltaEpiShift deltaEpiScale jrelTissueScale; do
  count=$(grep -c "CONSTANTS\[$c\]\s*=" src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H)
  echo "$c: assigned $count times (expected: 1)"
done
```
Expected: each constant assigned exactly once.

- [ ] **Step 4: Verify each new constant is read in the hot path**

```bash
for c in deltaEpiAmp deltaEpiShift deltaEpiScale jrelTissueScale; do
  count=$(grep -c "CONSTANTS\[$c\]" src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H)
  echo "$c: total references $count (expected: 2 — 1 init + 1 hot-path)"
done
```
Expected: each constant appears exactly 2 times.

- [ ] **Step 5: Verify hot-path tissue branches are gone**

```bash
grep -n "tissueFlag" src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H | awk -F: '$1+0 > 800 {print}'
```
Expected: no output (no `tissueFlag` references past line 800, i.e., outside `ToRORd_dynClinitConsts` which ends well before line 800).

- [ ] **Step 6: Verify build passes**

```bash
./Allwmake 2>&1 | tail -5
```
Expected: success.

- [ ] **Step 7: Sign off**

If all checks pass, report PASS. Otherwise, send the failed check + offending lines back to the implementer agent.

---

## Phase 4: Final integrated validation

Owned by: **validation agent**, in the same worktree as Phase 0 baselines.

### Task 4.1: Run all 7 pairs and bit-identical-diff against baseline

**Files:**
- Read: `~/cardiacFoam-baselines/<pair>/` (Phase 0 outputs)
- Read: `tutorials/singleCellprotocols/singleCell/`
- Create: `~/cardiacFoam-postrefactor/<pair>/` (this task's outputs)
- Create: `~/cardiacFoam-validation-report.md` (final report)

- [ ] **Step 1: Confirm refactor is at HEAD**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git log --oneline -5
```
Expected: top three commits are the Gaur, TNNP, ToRORd_dynCl refactors. The fourth commit should match the SHA recorded in `~/cardiacFoam-baselines/HEAD.txt`.

- [ ] **Step 2: Build with refactor applied**

```bash
./Allwmake 2>&1 | tee /tmp/build-postrefactor.log | tail -10
```
Expected: build completes without errors. Compare against `~/cardiacFoam-baselines/build-baseline.log` — same set of libraries should build, no new warnings.

- [ ] **Step 3: Set up post-refactor output directory**

```bash
mkdir -p ~/cardiacFoam-postrefactor
```

- [ ] **Step 4: Reuse the baseline runner with a different OUT path**

Create `/tmp/run_postrefactor.sh` (copy of `/tmp/run_baseline.sh` with the `OUT` path changed):

```bash
#!/bin/bash
set -e
MODEL=$1
TISSUE=$2
TUT=/Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/singleCellprotocols/singleCell
OUT=~/cardiacFoam-postrefactor/${MODEL}_${TISSUE}
mkdir -p "$OUT"

cd "$TUT"
cp constant/electroProperties constant/electroProperties.orig

python3 -c "
import re
p = 'constant/electroProperties'
s = open(p).read()
s = re.sub(r'(ionicModel\s+)\w+;', r'\g<1>${MODEL};', s)
s = re.sub(r'(tissue\s+)\w+;',     r'\g<1>${TISSUE};', s)
open(p, 'w').write(s)
"

./Allclean 2>/dev/null || true
CF_SKIP_PLOTS=1 ./Allrun

cp -r postProcessing "$OUT/" 2>/dev/null || true
cp log.cardiacFoam   "$OUT/" 2>/dev/null || true
for d in [0-9]*; do
    [ -d "$d" ] && cp -r "$d" "$OUT/" 2>/dev/null || true
done

mv constant/electroProperties.orig constant/electroProperties
echo "Post-refactor captured: $OUT"
```

```bash
chmod +x /tmp/run_postrefactor.sh
```

- [ ] **Step 5: Run all 7 pairs**

```bash
/tmp/run_postrefactor.sh TNNP         endocardialCells
/tmp/run_postrefactor.sh TNNP         mCells
/tmp/run_postrefactor.sh TNNP         epicardialCells
/tmp/run_postrefactor.sh ToRORd_dynCl endocardialCells
/tmp/run_postrefactor.sh ToRORd_dynCl mCells
/tmp/run_postrefactor.sh ToRORd_dynCl epicardialCells
/tmp/run_postrefactor.sh Gaur         myocyte
```
Expected: each completes successfully. If any fails to run that succeeded in Phase 0, report immediately — that's a bigger issue than a numerical drift.

- [ ] **Step 6: Bit-identical diff per pair**

For each pair, diff the time-step directories and any `postProcessing/` content. The simplest robust diff is `diff -ur` over the captured directories (will catch any byte-level difference in the output files). The expected outcome is that the only differences are timestamps inside `log.cardiacFoam` (which contains run-time info and should be excluded).

```bash
mkdir -p ~/cardiacFoam-diffs
for pair in TNNP_endocardialCells TNNP_mCells TNNP_epicardialCells \
            ToRORd_dynCl_endocardialCells ToRORd_dynCl_mCells ToRORd_dynCl_epicardialCells \
            Gaur_myocyte; do
    diff -ur \
        --exclude='log.cardiacFoam' \
        --exclude='*.log' \
        ~/cardiacFoam-baselines/$pair \
        ~/cardiacFoam-postrefactor/$pair \
        > ~/cardiacFoam-diffs/$pair.diff 2>&1 || true
    lines=$(wc -l < ~/cardiacFoam-diffs/$pair.diff)
    echo "$pair: $lines diff lines"
done
```
Expected: every pair shows `0 diff lines`. If any pair shows >0 lines:
- If diff is small (< 20 lines) and entirely within numerical fields → likely a documented operator-reordering issue (see spec §6 final paragraph). Investigate before reporting failure.
- If diff is large or affects the structure → real bug. Report immediately with the diff content for the failing pair.

- [ ] **Step 7: Run leftover-grep success criterion**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git grep -nE "tissueFlag" src/ionicModels/TNNP src/ionicModels/ToRORd_dynCl src/ionicModels/Gaur
```
Expected matches **only** in:
- `TNNP_2004Names.H` line referring to `TNNPinitConsts`
- `TNNP_2004.H` `TNNPinitConsts` definition signature plus the per-tissue ternaries inside its body
- `ToRORd_dynCl_2023Names.H` line referring to `ToRORd_dynClinitConsts`
- `ToRORd_dynCl_2023.H` `ToRORd_dynClinitConsts` definition signature plus the per-tissue ternaries inside its body
- `Gaur_2021Names.H` line referring to `GaurinitConsts`
- `Gaur_2021.H` `GaurinitConsts` definition signature

NO matches inside any `*computeRates`, `*computeVariables`, or `*Compact` body or signature.

- [ ] **Step 8: Write validation report**

Create `~/cardiacFoam-validation-report.md`:

```markdown
# Hot-path tissue-flag removal — validation report

**Baseline HEAD:** $(cat ~/cardiacFoam-baselines/HEAD.txt)
**Post-refactor HEAD:** $(git -C /Users/simaocastro/noFrontendCardiacFoam_minor_errors rev-parse HEAD)

## Build
- Baseline: PASS / FAIL (with link to ~/cardiacFoam-baselines/build-baseline.log)
- Post-refactor: PASS / FAIL (with link to /tmp/build-postrefactor.log)

## Bit-identical diff per (model, tissue) pair

| Pair | Diff lines | Status |
|---|---|---|
| TNNP / endocardialCells | N | PASS / FAIL |
| TNNP / mCells | N | PASS / FAIL |
| TNNP / epicardialCells | N | PASS / FAIL |
| ToRORd_dynCl / endocardialCells | N | PASS / FAIL |
| ToRORd_dynCl / mCells | N | PASS / FAIL |
| ToRORd_dynCl / epicardialCells | N | PASS / FAIL |
| Gaur / myocyte | N | PASS / FAIL |

## Leftover-reference grep
- Output of `git grep -nE "tissueFlag" src/ionicModels/{TNNP,ToRORd_dynCl,Gaur}`
- Verified all matches are within `*initConsts` only: PASS / FAIL

## Overall
PASS / FAIL — explain any failure in detail.
```

Fill in actual numbers and PASS/FAIL values, then report the contents of this file back to the orchestrator.

- [ ] **Step 9: Do NOT commit anything in this phase**

Phase 4 is purely verification — no code changes.

---

## Phase 5: Clamp inventory report (parallel, independent)

Owned by: **clamp-checker agent**, can run any time after Phase 0.

### Task 5.1: Enumerate clamp requirements per model for the eventual SoA-Euler path

**Files:**
- Read: `src/ionicModels/TNNP/TNNP_2004.H`, `TNNP_2004Names.H`
- Read: `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H`, `ToRORd_dynCl_2023Names.H`
- Read: `src/ionicModels/Gaur/Gaur_2021.H`, `Gaur_2021Names.H`
- Create: `docs/superpowers/specs/clamp-inventory.md`

- [ ] **Step 1: For each of the three models, list every state in the `STATES_INDEX` enum and classify it**

For each model open `<Model>_<year>Names.H` and find the `STATES_INDEX` enum. For each state, classify it as one of:

- **gate** — variable representing a probability or fraction; mathematical range `[0, 1]`. Common names contain `m`, `h`, `j`, `d`, `f`, `s`, `r`, `Xr1`, `Xr2`, `Xs`, `g`, `xs`, `xr`, etc.
- **concentration** — ion concentration; mathematical range `≥ 0` (typically clamped to `≥ SMALL` to avoid `log(0)` in subsequent calls). Common names: `Na_i`, `K_i`, `Ca_i`, `Ca_SR`, `Cli`, `Clss`, `Nai`, `Nass`, `Ki`, `Kss`, etc.
- **voltage** — `V` or `Vm`; no clamp.
- **other / structural** — for example CaMK trapping fractions (gate-like, [0,1]), Markov-chain Hodgkin–Huxley state probabilities (must sum to 1, but each in [0,1]); list with rationale.

- [ ] **Step 2: For each state classified as `gate` or `concentration`, decide the clamp**

- gate → clamp to `[0.0, 1.0]`
- concentration → clamp to `>= SMALL` (where `SMALL` is the OpenFOAM small-positive constant; the agent does not need to know its numerical value, just the symbol)
- voltage / other → no clamp

If any state is genuinely ambiguous (e.g. `f_inf` in TNNP could be argued for either category), default to `gate` and add a note flagging it for the implementer to decide during shim factoring.

- [ ] **Step 3: Write `docs/superpowers/specs/clamp-inventory.md`**

```markdown
# SoA-Euler clamp inventory — TNNP, ToRORd_dynCl, Gaur

**Generated:** 2026-05-12
**Purpose:** Per-model list of `(stateName, clampType, bound)` tuples that the next-session SoA-shim factoring will consume. Per `validation_strategy.md §2.2` of the upstream cardiacFoam GPU port: clamps are essential for detailed CellML-generated models because explicit Euler can drift gates outside `[0,1]` or concentrations negative within one substep, and the next call's `log(neg)` or `1/(1+exp(...))` will then NaN.

## TNNP (17 states)

| State | Classification | Clamp |
|---|---|---|
| V | voltage | none |
| K_i | concentration | `>= SMALL` |
| Na_i | concentration | `>= SMALL` |
| Ca_i | concentration | `>= SMALL` |
| ... (fill in all 17) ... | | |

## ToRORd_dynCl (~40 states)

| State | Classification | Clamp |
|---|---|---|
... (fill in all states) ...

## Gaur (28 states)

| State | Classification | Clamp |
|---|---|---|
... (fill in all 28) ...

## Notes / ambiguous cases

- (any per-state notes the agent flagged in Step 2)
```

- [ ] **Step 4: Commit the inventory**

```bash
cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors
git add docs/superpowers/specs/clamp-inventory.md
git commit -m "docs(superpowers): clamp inventory for TNNP/ToRORd_dynCl/Gaur

Per-model list of (stateName, clampType, bound) tuples for the
eventual SoA-Euler path, generated alongside the hot-path tissue-flag
removal refactor. Consumed by the next-session SoA shim factoring.

Sibling deliverable to docs/superpowers/specs/2026-05-12-hot-path-
tissue-removal-design.md."
```

- [ ] **Step 5: Report**

Report to the orchestrator: clamp inventory committed at `<sha>`. List total counts per model: gates / concentrations / other, and any ambiguous cases the implementer should know about for the shim work.

---

## Success criteria (from spec §10)

The plan is complete when **all** of the following are true:

1. `git grep -nE "tissueFlag\s*[=!,)]" src/ionicModels/{TNNP,ToRORd_dynCl,Gaur}` returns matches only inside `*initConsts` function bodies and signatures.
2. `Allwmake` builds the touched libraries and dependent solvers without errors or new warnings.
3. The validation report from Task 4.1 Step 8 reports PASS for all 7 `(model, tissue)` pairs.
4. The clamp inventory at `docs/superpowers/specs/clamp-inventory.md` is committed.
5. The four refactor commits (Gaur, TNNP, ToRORd_dynCl, plus the clamp-inventory) are on the working branch and the design doc commit is the parent of the first refactor commit.

---

End of plan.
