# coupledPVJ Magnitude Jump After uniform -> linear Kernel Switch

**Status: open, not yet investigated.** This document exists so a future
session doesn't have to re-derive the context below from git history alone.

## 1. Background — the bug this was a side effect of fixing

`purkinjeNiedererEtAl2011`'s regression test (Phase 1, `coupledPVJ` checks)
was non-deterministic across identical runs: `pvj6`'s "final" value bounced
between `10449.716` / `11110.503` / `13648.917` (and possibly other values)
depending on run, while every other PVJ (`pvj0`-`pvj5`, `pvj7`) was bit-
reproducible to ~1e-4.

Root cause (confirmed, not just suspected): the tutorial's mesh is a uniform
0.2mm hex grid, and `pvj6`'s terminal location `(0.017, 0.0015, 0.0015)` sits
exactly on a cell face in x and exactly at cell centres in y/z. With
`pvjRadius = 3e-4` m, several of `pvj6`'s "corner" cells (offset half a cell
in x, one full cell in y, one full cell in z) sit at
`sqrt(0.1e-3^2 + 0.2e-3^2 + 0.2e-3^2) = 3.0000...e-4` m — **exactly** equal to
`pvjRadius`. `pvjMapper.C`'s coupling-sphere membership test
(`r <= radius_`, see `src/electroModels/electroCouplers/pvjCoupler/pvjMapper.C`
around line 121) is a hard cutoff evaluated per-cell on the **locally
decomposed** mesh; `decomposeParDict` uses `method scotch`, whose partition
of the mesh is not guaranteed bit-reproducible run to run, and a boundary
cell reconstructed via a different point/face summation order can flip
`r <= radius_` by a single ULP. Combined with `pvjKernel uniform` (weight is a
hard 0/1, not tapered), a single boundary cell flipping in/out changes
`sphereVolumes_[6]` — and since the whole coupling is apparently running near
a stiffness/blow-up threshold for every PVJ (all "final" values are large,
by design or not — see open question below), that conductance change was
enough to send the explicit-scheme blow-up to a different magnitude each run.

`pvj7`, despite being only 0.5mm from `pvj6` (their coupling spheres
genuinely overlap: `pvjRadius=3e-4`, separation `5e-4 < 2*3e-4`), has a
different tie pattern (`x` and `z` both on cell faces, not `x` alone) whose
worst-case corner distance is `2.449e-4` m — comfortably inside the radius,
no exact tie, hence fully reproducible.

## 2. Fix applied

Switched `pvjKernel` from `uniform` to `linear` in
`tutorials/NiedererEtAl2011/purkinjeNiedererEtAl2011/constant/electroProperties.monodomain`.
`linear` tapers weight to ~0 as `r -> radius_`, so a boundary cell flipping
membership now changes the weighted sum by ~0 instead of a full `1.0` —
this doesn't remove the underlying decomposition-dependent last-bit
non-determinism, it makes its *consequence* negligible.

Verified deterministic: 5 fresh runs (fresh `cp -a` copy + full regression
each time, in an isolated `/private/tmp` clone built from scratch — see
session history) all produced bit-identical `coupledPVJ` values. Reference
file `tutorials/NiedererEtAl2011/purkinjeNiedererEtAl2011/regression/purkinjeSlab.reference`
was updated to match. `regressionTest.sh` now passes cleanly and
repeatably.

## 3. Open question — why did every PVJ's magnitude jump ~2-4x?

Switching `uniform` -> `linear` didn't just stabilize `pvj6`; **all 8** PVJ
final values shifted substantially, uniformly upward:

| PVJ | uniform (old ref) | linear (new ref) | ratio |
|---|---|---|---|
| pvj0 | 35699.608 | 155914.66 | 4.37x |
| pvj1 | 38596.205 | 105094.24 | 2.72x |
| pvj2 | 33803.03  | 155914.64 | 4.61x |
| pvj3 | 44286.159 | 166335.47 | 3.76x |
| pvj4 | 44286.162 | 166335.47 | 3.76x |
| pvj5 | 37918.556 | 161511.86 | 4.26x |
| pvj6 | 11110.503 (unstable) | 46464.865 | 4.18x |
| pvj7 | 14183.262 | 46300.933 | 3.26x |

That's a bigger, more uniform shift than "a few boundary cells changed
weight" would predict on its own, and is worth treating as suspicious rather
than just accepted because the tests pass. Two suspects, neither checked yet:

1. **Normalization mismatch in `pvjMapper.C`.** `sphereVolumes_[i]` is a
   weighted sum (`w * cellVolumes[cellI]`, weight either `1.0` for uniform
   or `1 - r/radius` for linear). Swapping kernels shrinks the *average*
   weight per included cell (linear's mean weight over a sphere is less than
   uniform's), so `sphereVolumes_[i]` shrinks. `volumetricSource` divides
   `couplingCurrent[i] / sphereVolumes_[i]` — if `couplingCurrent[i]` itself
   isn't rescaled to match, this ratio inflates. Whether that's intended
   (each kernel is a genuinely different physical deposition profile) or a
   bug (linear should conserve total deposited current the way uniform does,
   just reshape it spatially) has not been checked.
2. Whether the entire "coupledPVJ final" quantity is even meant to reach
   values in the 1e4-1e5 range at all — nothing in this investigation
   confirmed what physical quantity `postProcessing/purkinjeNetwork.dat`
   column 8/9 (etc.) actually represents in absolute terms, or whether a
   blow-up to O(1e4-1e5) by `t=0.02` (endTime) is expected model behavior
   or itself a symptom of the `explicit` `pvjCouplingScheme` being stiff at
   this `rPvj`/`pvjRadius`/mesh-resolution combination. Both `pvj6`/`pvj7`
   time series (see session history / `postProcessing/purkinjeNetwork.dat`)
   show a jump of ~1e8x in the *final* timestep alone (from ~1e-4 at
   `t=0.015` to ~1e4 at `t=0.02`), which independently smells like a
   near-instability regardless of kernel choice.

## 4. `explicit` vs `implicit` pvjCouplingScheme — tested, hypothesis falsified

Set `pvjCouplingScheme implicit` (in the same `electroProperties.monodomain`
`domainCouplings.couplingA` block) and reran 3x fresh. Result: **also fully
deterministic** (identical across all 3 runs — the kernel fix's determinism
holds regardless of coupling scheme), but the magnitude did **not** collapse
toward something physiologically modest the way an explicit-instability
artifact should have:

| PVJ | linear + explicit | linear + implicit | change |
|---|---|---|---|
| pvj0 | 155914.66 | 124743.74 | -20% |
| pvj1 | 105094.24 | 82539.617 | -21% |
| pvj2 | 155914.64 | 124752.8  | -20% |
| pvj3 | 166335.47 | 153383.09 | -8% |
| pvj4 | 166335.47 | 153383.12 | -8% |
| pvj5 | 161511.86 | 136942.39 | -15% |
| pvj6 | 46464.865 | 46630.119 | +0.4% |
| pvj7 | 46300.933 | 46406.468 | +0.2% |

Going implicit only trims 8-21% off most values and barely moves `pvj6`/
`pvj7` at all (well within noise of the kernel-shape difference alone). If
the O(1e4-1e5) magnitude were an explicit-scheme stability artifact,
switching to an unconditionally-stable implicit coupling should have
produced values an order of magnitude smaller, not a ~10-20% trim. **This
rules out "explicit-scheme instability" as the explanation for the large
magnitude.** The large values are very likely a genuine (if still
unexplained) property of the model/units as currently configured, not a
numerical blow-up.

Decision: kept `pvjCouplingScheme implicit` (it is at least as sound
numerically as explicit and no worse for determinism), reference file
updated to match the implicit values above (right-hand column), committed.

## 5. Suggested next steps (magnitude question still open)

- Trace `couplingCurrentAtPvjs` / `volumetricSource` /
  `depositCoupling` in `reactionDiffusionPvjCoupler.C` and `pvjMapper.C`
  end-to-end for one PVJ, by hand, with `debugCoupling true` (already on in
  this tutorial's dict) to see exactly which term produces the O(1e4-1e5)
  magnitude, and whether that's expected given `rPvj=10000`,
  `pvjRadius=3e-4`, and this mesh's cell volumes.
- Check what physical quantity `postProcessing/purkinjeNetwork.dat` column
  8/9 (`coupledPVJ pvjN final`) actually represents in absolute terms —
  this was never confirmed. Given `implicit` ruled out instability, the
  magnitude itself (not its stability) is now the open question.
- Decide whether `sphereVolumes_`/`couplingCurrent` normalization in
  `pvjMapper.C` needs a fix so `uniform` and `linear` kernels are
  physically comparable (same total deposited current, just different
  spatial spread), given the kernel switch alone changed magnitudes by
  3-4x (section 3) independent of the explicit/implicit question.

## 6. Where things stand in the repo

- `tutorials/NiedererEtAl2011/purkinjeNiedererEtAl2011/constant/electroProperties.monodomain`
  — `pvjKernel` changed `uniform` -> `linear`; `pvjCouplingScheme` set to
  `implicit`.
- `tutorials/NiedererEtAl2011/purkinjeNiedererEtAl2011/regression/purkinjeSlab.reference`
  — all 8 `coupledPVJ` values updated to the deterministic `linear` +
  `implicit` results (right-hand column in section 4's table).
- `regressionTest.sh` for this tutorial passes cleanly and reproducibly
  with the above. Committed to `ep-work-onto-main`.
