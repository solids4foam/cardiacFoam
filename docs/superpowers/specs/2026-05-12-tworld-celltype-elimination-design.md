# TWorld AC_celltype elimination + hot-path tissue removal — design

**Status:** draft, awaiting user review
**Date:** 2026-05-12
**Scope:** TWorld only
**Out of scope:** other ionic models; heterogeneity wiring; SoA shim work; CUDA build
**Baseline ref:** commit `26e1a50` (= current HEAD, the IKr renormalisation commit)
**Baseline traces:** `~/cardiacFoam-baselines-tworld/{TWorld_endocardialCells, TWorld_mCells, TWorld_epicardialCells}/postProcessing/` captured before any code change

---

## 1. Goal

Bring TWorld in line with the rest of the ionic-model family (TNNP, ToRORd_dynCl, Gaur, BuenoOrovio): tissue dispatch happens via direct `(tissueFlag == K)` ternaries inside `TWorldinitConsts`; the hot path has zero tissue-related branches; tissue-dependent values live entirely in the produced per-cell `CONSTANTS` array.

This eliminates a **latent heterogeneity bug**: TWorld currently routes tissue dispatch through `CONSTANTS[AC_celltype]` (a double-valued `1.0/2.0/0.0`). The 9 init-time branches on `AC_celltype` are safe in the discrete-tissue case because each tissue template is built with a single `tissueFlag`. But the one hot-path branch at `TWorld_2024.H:1542` re-reads `AC_celltype` *after* per-cell heterogeneity blending would happen — and a 50/50 endo (0.0) + M-cell (2.0) blend produces `AC_celltype = 1.0`, which the hot-path branch would mis-dispatch as "epi" and apply a spurious 0.7 GNaL scaling.

The fix has two parts that must land together:
1. Eliminate the `AC_celltype` indirection (Family 2 → Family 1 pattern alignment).
2. Promote the one hot-path scaling factor to its own per-cell constant.

---

## 2. The pattern after the change

```
TWorldinitConsts(tissueFlag):
  // 11 init-time tissue-dependent constants, set via direct tissueFlag ternaries:
  CONSTANTS[AC_INaCa_celltype_factor] = (tissueFlag == 1) ? 1.1 : (tissueFlag == 2) ? 1.4 : 1.0;
  // ... 10 more like this ...
  // 1 new constant for the hot-path scaling:
  CONSTANTS[gnalTissueScale] = (tissueFlag == 1) ? 0.7 : 1.0;

TWorldcomputeVariables(...):  // no tissueFlag parameter anymore
  // line 1542 becomes branch-free:
  ALGEBRAIC[AV_GNaL] = CONSTANTS[AC_GNaL_b] * CONSTANTS[AC_INaL_multiplier]
                     * (1.0 + ALGEBRAIC[AV_fINaLp]) * CONSTANTS[gnalTissueScale];
```

`TWorldCellType` helper, `AC_celltype` enum entry, and `"AC_celltype"` name string are all removed.

---

## 3. Concrete changes

### 3.1 `src/ionicModels/TWorld/TWorld_2024Names.H`

- **Remove** `AC_celltype,` (line 627) from `CONSTANTS_INDEX` enum.
- **Append** `gnalTissueScale,` immediately before `NUM_CONSTANTS` (preserving the convention that newly-added constants go at the end so existing index values shift only by the removal).
- **Drop** `int tissueFlag,` from the `TWorldcomputeVariables` forward declaration. Keep it on `TWorldinitConsts`.

Net `NUM_CONSTANTS` change: -1 (remove `AC_celltype`) + 1 (add `gnalTissueScale`) = **0**.

### 3.2 `src/ionicModels/TWorld/TWorld_2024.H`

- **Remove** the `TWorldCellType` inline helper function (lines 813–820).
- **Remove** the `"AC_celltype",` entry from the `TWorldCONSTANTS_NAMES` array (line 626).
- **Append** `"gnalTissueScale",` to the `TWorldCONSTANTS_NAMES` array, immediately before the closing `};`, in the same enum position as the new entry (i.e. at the end).
- **Remove** the line `CONSTANTS[AC_celltype] = TWorldCellType(tissueFlag);` (line 857).
- **Rewrite** the 9 init-time branches (lines 900, 936, 988, 997, 998, 1012, 1024, 1036, 1042): replace `(CONSTANTS[AC_celltype] == 1.0)` → `(tissueFlag == 1)`, `(CONSTANTS[AC_celltype] == 2.0)` → `(tissueFlag == 2)`, and the implicit endo case stays as the trailing `else`. (Exhaustive `grep` confirms there are no other init-time `AC_celltype` reads in the file.)
- **Add** the new `gnalTissueScale` assignment alongside the existing tissue-dependent block in `TWorldinitConsts`:
  ```cpp
  CONSTANTS[gnalTissueScale] = (tissueFlag == 1) ? 0.7 : 1.0;
  ```
- **Drop** `int tissueFlag,` from the `TWorldcomputeVariables` definition signature (line 1252).
- **Replace** the hot-path branch at line 1542 with the branch-free unified form:
  ```cpp
  ALGEBRAIC[AV_GNaL] = CONSTANTS[AC_GNaL_b] * CONSTANTS[AC_INaL_multiplier]
                     * (1.0 + ALGEBRAIC[AV_fINaLp]) * CONSTANTS[gnalTissueScale];
  ```

### 3.3 `src/ionicModels/TWorld/TWorld.C`

- **Drop** `tissue()` from the three `::TWorldcomputeVariables(...)` call sites (around lines 143, 172, 245). The fourth occurrence — line 73, inside the `TWorldinitConsts` call — stays.

### 3.4 No other files touched

- `Make/files`: TWorld already registered (`TWorld/TWorld.C` at line 17).
- BuenoOrovio, TNNP, ToRORd_dynCl, Gaur: untouched.
- Tutorial dictionaries: untouched (the singleCell tutorial works for TWorld today).

---

## 4. Validation

**Baseline** (already captured in this session, before any code change):
- `~/cardiacFoam-baselines-tworld/HEAD.txt` — pinned to commit `26e1a50`.
- `~/cardiacFoam-baselines-tworld/TWorld_{endocardialCells,mCells,epicardialCells}/postProcessing/` — three trace files (~428K each), produced by `tutorials/singleCellprotocols/singleCell/` with `ionicModel TWorld` and `tissue {endo,M,epi}` swapped via the dict.

**Post-refactor procedure:**
1. Run the build: `source /Volumes/OpenFOAM-v2412/etc/bashrc && ./Allwmake`. Must report "There were no build errors: enjoy cardiacFoam!".
2. Re-run the singleCell tutorial for the same three pairs into `~/cardiacFoam-postrefactor-tworld/TWorld_{tissue}/`.
3. `diff -ur --exclude='log.cardiacFoam' --exclude='*.log'` baseline vs post-refactor per pair.
4. Expect **0 diff lines** per pair (bit-identical).

**Bit-identity argument:** The init-time math reduces algebraically. Old `(CONSTANTS[AC_celltype] == 1.0) ? a : b` is true iff `tissueFlag == 1` (because `TWorldCellType(tissueFlag)` returns `1.0` exactly when `tissueFlag == 1`). Same for `== 2.0` ↔ `tissueFlag == 2`. So every old branch maps to the new branch with identical operand pairs and identical evaluation — IEEE result must be the same bit-for-bit.

Hot-path branch: old `(CONSTANTS[AC_celltype] == 1.0) ? K * 0.7 : K` (where K is the rest of the expression). New `K * CONSTANTS[gnalTissueScale]` with `gnalTissueScale = 0.7` for epi, `1.0` for non-epi. These are mathematically identical; the FP result in the epi branch is `K * 0.7` either way. In the non-epi branch the old code returns exactly `K`; the new code returns `K * 1.0`, which is bit-identical to `K` in IEEE.

**Failure mode:** if a diff appears, the most likely causes are (a) operator-reordering inside one of the rewritten init ternaries (e.g. associativity change in the multi-factor expressions for `AC_PCa`), or (b) a missed branch reference still reading `AC_celltype` after removal (would be a build error). Recovery is to keep the original endo/M/epi branch operand structure verbatim, just replace the ternary condition.

---

## 5. Leftover-grep success criteria

After the refactor:

```bash
grep -nE "AC_celltype|TWorldCellType" src/ionicModels/TWorld/
```
Expected: zero matches.

```bash
grep -nE "tissueFlag" src/ionicModels/TWorld/TWorld_2024.H
```
Expected: matches only inside `TWorldinitConsts` (signature line ~823, plus the 12 ternaries and the new `gnalTissueScale` assignment in its body). Zero matches in `TWorldcomputeVariables` body or signature.

```bash
grep -nF "tissue()" src/ionicModels/TWorld/TWorld.C
```
Expected: one match — the `TWorldinitConsts` call at line 73. No `tissue()` argument in any `TWorldcomputeVariables` call.

```bash
source /Volumes/OpenFOAM-v2412/etc/bashrc && ./Allwmake 2>&1 | tail -5
```
Expected: "There were no build errors: enjoy cardiacFoam!".

---

## 6. What does NOT change

- No new ionic-model registration in `Make/files` — TWorld is already there.
- No tutorial dictionary edits.
- `TWorldinitConsts` keeps its `int tissueFlag` parameter — same as every other model.
- `NUM_CONSTANTS` stays at its current value (one removed, one added).
- The 9 init-time tissue-dependent CONSTANTS (`AC_INaCa_celltype_factor`, `AC_Vmax_SRCaP`, `AC_IK1_celltype_factor`, `AC_Gto_fast`, `AC_Gto_slow`, `AC_PCa`, `AC_PCa_P`, `AC_IKr_celltype_factor`, `AC_IKs_celltype_factor`) keep their values — only the *condition* in the ternary changes from `AC_celltype` to `tissueFlag`.

---

## 7. Out of scope (explicit anti-goals)

- **No** clamp-inventory work in this spec. TWorld extends the inventory; that's a separate follow-up after this refactor lands.
- **No** `configureIonicHeterogeneity` for TWorld. Same rationale as the previous spec — heterogeneity wiring requires consumer-side context not in scope.
- **No** SoA shim, `*Compact.H`, `*Batch.H`, or `*_cuda.cu` for TWorld. Same scope boundary as the previous refactor.
- **No** changes to other ionic models.
- **No** integrator changes, no clamp implementation.

---

## 8. Success criteria

This refactor is complete when:

1. Build clean under `Allwmake`.
2. Bit-identical traces vs `~/cardiacFoam-baselines-tworld/` for all three tissues.
3. Leftover-grep checks in §5 all pass.
4. `git log` shows one or two new commits on `no-frontend-minor-errors` containing the three file changes (TWorld_2024.H, TWorld_2024Names.H, TWorld.C) — and no other files touched by this work.

---

End of design.
