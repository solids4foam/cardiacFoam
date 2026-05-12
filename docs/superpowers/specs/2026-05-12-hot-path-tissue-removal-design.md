# Hot-path tissue-flag removal — design

**Status:** draft, awaiting user review
**Date:** 2026-05-12
**Scope:** TNNP, ToRORd_dynCl, Gaur (signature-only)
**Out of scope:** TNNP/ToRORd_dynCl heterogeneity wiring; BuenoOrovio (already done); other ionic models

---

## 1. Goal

Make every per-cell ionic-model property — including tissue-class differences — live entirely in the per-cell `CONSTANTS` vector, so the SoA evaluator's hot path has zero tissue-related branching. This unblocks per-cell tissue dispatch on a single mesh (the existing `BuenoOrovio` heterogeneity pattern) and makes the SoA kernel signature uniform across models.

This is **GPU-prep work**, not feature work. It does not add heterogeneity wiring for any new model. It changes the shape of the legacy CPU evaluators so that any future heterogeneity work — for TNNP, ToRORd_dynCl, or any extension of the BuenoOrovio scheme — becomes a localised addition that does not touch the hot path again.

---

## 2. The pattern (one paragraph)

Today, three model headers branch on a scalar `int tissueFlag` parameter inside their hot-path functions (`computeRates` / `computeVariables`). For each branch, find a unified algebraic form that covers both sides — typically `prefactor / (1 + exp((V + shift) / scale))` plus an offset, where one branch's `prefactor = 0` collapses the term entirely. Promote the literal numbers to named entries in the `CONSTANTS` enum, set per-tissue values inside the existing `initConsts` tissue branch (which keeps `tissueFlag` because it runs once per tissue class at startup), and drop `tissueFlag` from every hot-path function signature. After this, the SoA evaluator's signature matches `BuenoOrovioComputeVariablesBatch`: just `CONSTANTS_SoA`, `STATES_SoA`, `RATES_SoA`, `SUPPORT_SoA` plus the time/stimulus arguments.

---

## 3. Per-model changes

### 3.1 TNNP

Hot-path branches (4 occurrences in 2 functions):

| Location | What branches |
|---|---|
| [TNNP_2004.H:426-431](src/ionicModels/TNNP/TNNP_2004.H:426) inside `TNNPcomputeRates` | `s_inf`, `tau_s` |
| [TNNP_2004.H:523-529](src/ionicModels/TNNP/TNNP_2004.H:523) inside `TNNPcomputeVariables` | `s_inf`, `tau_s` |

**`s_inf` unification** — both branches share the same shape, just promote the three literals:

```cpp
ALGEBRAIC[s_inf] = CONSTANTS[sInfPrefactor]
    / (1.0 + exp((STATES[V] + CONSTANTS[sInfShift]) / CONSTANTS[sInfScale]));
```

**`tau_s` unification** — endo has 3 terms (gaussian + sigmoid + offset), non-endo has 2 (gaussian + offset). Use the 3-term superset and zero the sigmoid amplitude in non-endo:

```cpp
ALGEBRAIC[tau_s] =
      CONSTANTS[tauSGaussAmp]
      * exp(-pow(STATES[V] + CONSTANTS[tauSGaussShift], 2.0)
            / CONSTANTS[tauSGaussWidth])
    + CONSTANTS[tauSSigmoidAmp]
      / (1.0 + exp((STATES[V] + CONSTANTS[tauSSigmoidShift])
                   / CONSTANTS[tauSSigmoidScale]))
    + CONSTANTS[tauSOffset];
```

**New `CONSTANTS_INDEX` entries (10), inserted before `NUM_CONSTANTS`:**

| Constant | Endo (tissueFlag == 1) | M / Epi (tissueFlag != 1) |
|---|---|---|
| `sInfPrefactor` | 1.00 | 1.10 |
| `sInfShift` | 20.0 | 28.0 |
| `sInfScale` | 5.0 | 6.0 |
| `tauSGaussAmp` | 85.0 | 1000.0 |
| `tauSGaussShift` | 45.0 | 67.0 |
| `tauSGaussWidth` | 320.0 | 1000.0 |
| `tauSSigmoidAmp` | 5.0 | **0.0** (collapses term) |
| `tauSSigmoidShift` | -20.0 | 0.0 (irrelevant when amp=0) |
| `tauSSigmoidScale` | 5.0 | **1.0** (must be ≠0) |
| `tauSOffset` | 3.0 | 8.0 |

(The original TNNP code branches only on `tissueFlag == 1` for these algebraics — M-cells and epi share the same formula. The init code's separate `g_Ks` and `g_to` branches discriminate further between M and epi at init time, and stay as-is.)

`NUM_CONSTANTS` grows from 42 → 52.

### 3.2 ToRORd_dynCl

Hot-path branches (3 occurrences, all in `ToRORd_dynClcomputeVariables`):

| Location | What branches |
|---|---|
| [ToRORd_dynCl_2023.H:948](src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H:948) | `AV_delta_epi` (endo-only sigmoid term) |
| [ToRORd_dynCl_2023.H:1128-1129](src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H:1128) | `AV_Jrel_inf`, `AV_Jrel_infp` (M-cell-only 1.7× scaling) |

**`AV_delta_epi` unification** — same trick as TNNP `s_inf`, with `AC_EKshift` already in `CONSTANTS`:

```cpp
ALGEBRAIC[AV_delta_epi] = 1.0 - CONSTANTS[deltaEpiAmp]
    / (1.0 + exp((STATES[V] + CONSTANTS[AC_EKshift] + CONSTANTS[deltaEpiShift])
                 / CONSTANTS[deltaEpiScale]));
```

**Jrel scalings unification** — both lines share the same M-cell 1.7× factor, so a single new constant suffices:

```cpp
ALGEBRAIC[AV_Jrel_inf]  = ALGEBRAIC[AV_Jrel_inf_b]  * CONSTANTS[jrelTissueScale];
ALGEBRAIC[AV_Jrel_infp] = ALGEBRAIC[AV_Jrel_infp_b] * CONSTANTS[jrelTissueScale];
```

**New `CONSTANTS_INDEX` entries (4):**

| Constant | Endo | M-cell | Epi |
|---|---|---|---|
| `deltaEpiAmp` | 0.95 | 0.0 | 0.0 |
| `deltaEpiShift` | 70.0 | 0.0 | 0.0 |
| `deltaEpiScale` | 5.0 | 1.0 | 1.0 |
| `jrelTissueScale` | 1.0 | 1.7 | 1.0 |

### 3.3 Gaur

Zero hot-path branches. `GaurcomputeVariables` (the only hot-path function in Gaur — there is no separate `GaurcomputeRates`) takes `int tissueFlag` in its signature but never references it. Action: drop the parameter from the signature and callers. No new constants.

---

## 4. Files touched

This repo currently has only the legacy CPU wrapper (`<Model>.C`/`<Model>.H`), the cellML-generated math (`<Model>_<year>.H`), and the enum file (`<Model>_<year>Names.H`) per ionic model. There are no `*Compact.H`, `*Batch.H`, or `*Batched/` directories yet — those are part of the next-session SoA-shim work, not this refactor. The file lists below cover every file that currently exists and is affected.

### 4.1 TNNP
- `src/ionicModels/TNNP/TNNP_2004Names.H` — add 10 enum entries before `NUM_CONSTANTS`; drop `int tissueFlag` from forward declarations of `TNNPcomputeRates`, `TNNPcomputeVariables`.
- `src/ionicModels/TNNP/TNNP_2004.H` — add 10 per-tissue assignments inside the `TNNPinitConsts` tissue branch (after the existing `g_Ks` / `g_to` lines); replace 4 hot-path branches with the unified formulas; drop `int tissueFlag` from `TNNPcomputeRates` and `TNNPcomputeVariables` signatures.
- `src/ionicModels/TNNP/TNNP.C` — drop the `tissue()` argument from all four call sites (lines 146, 158, 191, 250 in current `main`).

### 4.2 ToRORd_dynCl
- `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023Names.H` — add 4 enum entries before `NUM_CONSTANTS`; drop `int tissueFlag` from forward declaration of `ToRORd_dynClcomputeVariables`.
- `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023.H` — add 4 per-tissue assignments inside the `ToRORd_dynClinitConsts` tissue branch; replace 3 hot-path branches with the unified formulas; drop `int tissueFlag` from `ToRORd_dynClcomputeVariables` signature.
- `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl.C` — drop the `tissue()` argument from the three call sites (lines 134, 167, 241).

### 4.3 Gaur
- `src/ionicModels/Gaur/Gaur_2021Names.H`, `src/ionicModels/Gaur/Gaur_2021.H` — drop `int tissueFlag` from `GaurcomputeVariables` signature/declaration.
- `src/ionicModels/Gaur/Gaur.C` — drop the `tissue()` argument from the three call sites (lines 134, 167, 241).

---

## 5. What does NOT change

- `*initConsts` functions keep `int tissueFlag` — they run once per tissue class at startup and use it to bake a per-tissue constants vector. This is the same shape as `BuenoOrovio::constantsForTissue(tissueFlag, dict())` and is the input to any future heterogeneity blender.
- `BuenoOrovio` source files. Reference implementation, untouched.
- No `TNNPBatched::configureIonicHeterogeneity` or `ToRORd_dynClBatched::configureIonicHeterogeneity` — heterogeneity wiring is deliberately out of scope (we don't have all the consumer-side context for it yet).
- Per-cell `CONSTANTS_SoA` layout, `batchedIonicCore`, `BatchedKernelExecutor`, dictionary keys, the SoA shim itself. The shim factoring is the next session's work; this session only changes what the shim *will* call.

---

## 6. Validation

Run by a dedicated agent (see §8 below). The single-cell tutorial that already exists in this repo is `tutorials/singleCellprotocols/singleCell/`, parameterised by `ionicModel` and `tissue` keys in `constant/electroProperties` (today set to `BuenoOrovio` + `endocardialCells`). The agent uses this tutorial as the validation harness, and `tutorialsTest-regression/singleCellprotocols/singleCell/` as the baseline-comparison source if it carries pre-refactor outputs.

For each `(model, tissue)` pair below the agent must demonstrate:

| Model | Tissues to sweep |
|---|---|
| TNNP | `endocardialCells`, `mCells`, `epicardialCells` |
| ToRORd_dynCl | `endocardialCells`, `mCells`, `epicardialCells` |
| Gaur | `myocyte` (single supported tissue per `Gaur::supportedTissueTypes()`) |

The validation procedure per pair:

1. **Capture pre-refactor trace.** Before any code changes land, run `tutorials/singleCellprotocols/singleCell/` once per `(model, tissue)` pair on the current `main`, save `Vm` (and any other exported field) trace to a baseline file outside the repo or in a worktree. Seven runs total.
2. **Run post-refactor.** After the refactor is implemented, run the same seven `(model, tissue)` pairs.
3. **Bit-identical diff.** Compare baseline vs post-refactor `Vm` trace point-for-point. Tolerance: zero ULPs in the typical case; ≤ 1 ULP acceptable only if the agent can demonstrate the difference comes from a documented operator-reordering case in §6 of this doc.
4. **Build clean.** `Allwmake` of the three touched ionic-model libraries and any dependent solver/utility (`cardiacFoam`, `ionicHeterogeneityProbe`) completes without errors and without new warnings.
5. **Leftover-reference grep.** `git grep -nE "tissueFlag" src/ionicModels/{TNNP,ToRORd_dynCl,Gaur}` returns matches only inside `*initConsts` function bodies and signatures; no occurrences in `*computeRates` or `*computeVariables`.

The bit-identity argument: with no heterogeneity dict configured, `HETEROGENEOUS_CONSTANTS_` stays empty (it is empty by default; only BuenoOrovio populates it today) and the wrapper class passes the homogeneous `CONSTANTS_` to the legacy function. That `CONSTANTS_` is populated by `*initConsts(tissueFlag)` and now contains the same numerical values that the hot-path branches previously evaluated to. So at the IEEE level, the new evaluator should produce the same RATES and ALGEBRAIC values for the same inputs.

If any trace fails to match bit-identically, the most likely cause is operator-reordering inside the unified formula (e.g. `pow((V+shift), 2.0) / width` vs `pow((V+shift)/sqrt(width), 2.0)` — algebraically equal, FP-distinct). The fix is to keep the operator structure of the original endo branch and let non-endo's zeroed amplitude collapse the term — never restructure the surviving branch.

---

## 7. Trade-off acknowledged (explicit)

This refactor commits the codebase to **parameter-level** blending semantics for transmural transition bands, not value-level. In a transition band:

- BuenoOrovio (already): `CONSTANTS_blended = w_endo * CONSTANTS_endo + w_nonE * CONSTANTS_nonE`, then evaluate once. For its model this is equivalent to value-level blending because every tissue-dependent quantity enters as a linear coefficient.
- TNNP / ToRORd_dynCl (after this change): same `CONSTANTS_blended`, then evaluate once. For their models this is **not** equivalent to value-level blending — `tau_s` and `AV_delta_epi` contain nonlinear functions of the tissue-blended parameters. A 50/50 cell will get a unimodal `tau_s` curve whose gaussian peak shifts and broadens nonlinearly with the blend, rather than the pointwise mean of the two pure-tissue curves.

This is the standard cardiac-modelling interpretation: a cell in a transition band is one cell with continuum-interpolated phenotype, not a probabilistic mixture of two distinct cell types. For typical `transitionWidth = 0.1` over a 1 mm band this is physiologically defensible. The alternative — value-level blending — would require keeping a per-cell weight buffer, computing both branches per cell, and combining values — measurably more expensive per step and inconsistent with the BuenoOrovio precedent. We are explicitly choosing the BuenoOrovio convention.

---

## 8. Implementation roles

The implementation plan (next session) will dispatch:

- **Implementer agent.** Per model in order (TNNP → ToRORd_dynCl → Gaur), edit the files in §4, advancing through three logical sections per model: (a) `*Names.H` enum + signature drops, (b) `*initConsts` per-tissue assignments, (c) `*computeRates` / `*computeVariables` unified formulas. Stop at each section boundary.
- **Reviewer agent.** Fires at each section boundary (not per keystroke). Verifies: every removed branch is replaced by the unified formula reading from a CONSTANTS entry that exists in the enum and is set in `initConsts` for every tissue; the operator structure of the surviving branch is preserved; no caller still passes `tissue()`; the file compiles in isolation.
- **Validation agent.** Runs §6 procedure end-to-end: captures pre-refactor traces from the `singleCellprotocols/singleCell` tutorial for all seven `(model, tissue)` pairs **before any code lands**, then re-runs after each model's implementation completes, diffs bit-identically, runs the build, and runs the leftover-reference grep. Reports per-pair pass/fail with traces. Owns its own working directory (worktree) so baseline runs are not contaminated by in-progress edits.
- **Clamp-checker agent (parallel, host-only).** Independent of the tissue refactor: enumerate which states in TNNP, ToRORd_dynCl, and Gaur need clamps for the eventual SoA-Euler path (gates → [0,1], concentrations → ≥SMALL). Output is a per-model list of `(stateName, clampType, bound)` tuples that the next-session shim factoring will consume. This agent does not modify code; it produces a markdown report at `docs/superpowers/specs/clamp-inventory.md`.

The implementation plan will spell out the agent-handoff cadence, the exact sub-task each agent owns, and the merge criteria. That document is written in the next session via `superpowers:writing-plans`.

---

## 9. Out of scope (explicit anti-goals)

To prevent scope creep during implementation:

- **No** `configureIonicHeterogeneity` for TNNP, ToRORd_dynCl, or any model. Wiring per-cell heterogeneity to the consumer side (myocardiumDomain, transmural distance fields, dictionary contracts) is a separate design that requires more codebase context than was available for this brainstorm.
- **No** new `*Compact.H`, `*Batch.H`, or `*Batched/` derived class for any model. Building the SoA shim layer is the next session's work, enabled by but distinct from this refactor.
- **No** new `*Batched_cuda.cu` placeholder kernels. CUDA scaffolding lands together with the SoA shims.
- **No** changes to BuenoOrovio. Reference implementation.
- **No** changes to ionic models outside the four GPU-port targets (`AlievPanfilov`, `Courtemanche`, `Fabbri`, `Grandi`, `ORd`, `Stewart`, `Trovato`, `tmanufacturedFDA`, `bidomainFDAManufactured`, `monodomainFDAManufactured`, `PerisYague`, `TWorld`).
- **No** integrator changes, no clamp implementation (the clamp-checker agent only reports; clamps land with the shim factoring next session).
- **No** profiling, no benchmarking, no GPU build (`HAS_CUDA` stays off for this change).

---

## 10. Success criteria

This refactor is complete when:

1. `git grep -nE "tissueFlag\s*[=!,)]" src/ionicModels/{TNNP,ToRORd_dynCl,Gaur}` returns matches only inside `*initConsts` function bodies and signatures.
2. `Allwmake` builds the touched libraries and dependent solvers without errors or new warnings.
3. The validation agent reports bit-identical traces vs current `main` for all seven supported `(model, tissue)` pairs (TNNP × {endo, M, epi} + ToRORd_dynCl × {endo, M, epi} + Gaur × myocyte) using `tutorials/singleCellprotocols/singleCell/` parameterised via the `electroProperties` dict.
4. The design doc is committed; the implementation plan that consumes it is committed.

---

End of design.
