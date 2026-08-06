# SBDF2 Order Reduction — Investigation Record

**Status: root cause found and fixed.** The scheme was missing one of the two
PDE↔ODE couplings that SBDF2 requires. See
[SBDF2_TIME_STEPPING_RATIONALE.md](SBDF2_TIME_STEPPING_RATIONALE.md) for the
scheme as it now stands; this document records how the defect was found, and
in particular records that the *first* diagnosis was incomplete and why its
supporting evidence was misleading.

## 1. Original symptom

Temporal convergence on the manufactured monodomain problem
(`monodomainPseudoECG`, diagonal conductivity, hex mesh) showed **both**
`godunov` and `sbdf2` converging at p ≈ 1.0, near-identically. `sbdf2` was
expected to reach p ≈ 2.0.

| dt | Godunov Vm L2 (1D) | SBDF2 Vm L2 (1D) |
|---|---|---|
| 0.003125 | 0.000141962 | 0.000142288 |
| 0.0015625 | 7.07047e-05 | 7.07879e-05 |
| 0.00078125 | 3.52375e-05 | 3.52585e-05 |
| 0.000390625 | 1.75445e-05 | 1.75498e-05 |

The two schemes differed by ~0.2% — i.e. the SBDF2 source extrapolation was
having almost no effect at all.

## 2. Ruled out early

- **OpenFOAM's `backward` ddtScheme.** Traced through `backwardDdtScheme.C`;
  for constant dt beyond the first step it is genuine BDF2,
  $(3V^{n+1}-4V^n+V^{n-1})/(2\Delta t)$. Correct.
- **The extrapolation formula and history bookkeeping.** Verified with
  temporary debug prints: `IionOld_`/`IionOldOld_` correctly fall back to
  Godunov for the first two steps, then hold real values; the correction
  $r(I^n - I^{n-1})$ tracks $\Delta t\,dI/dt$ to the right order of magnitude.
- **Dict-key scoping.** `timeCouplingScheme` lands in `monodomainSolverCoeffs`
  and is read at construction.

Two genuine bugs were found and fixed along the way (a non-null `IionOldPtr()`
causing a first-step doubling artifact, and an off-by-one history-readiness
gate that should be `timeIndex() >= 3`, not `>= 2`). Neither explained the
order.

## 3. First diagnosis — real, but not the root cause

The first diagnosis identified that `solveIonicCurrent` for step $k$ integrates
the gating variables to $t_k$ but with $V_m$ held at $V_m(t_{k-1})$, so the
value labelled "$I_{ion}$ at $t_k$" is really

$$I_{ion}\big(\mathbf{w}(t_k),\,V_m(t_{k-1})\big)$$

The fix was to refresh $I_{ion}$ after the diffusion solve, at the correct
$V_m(t_k)$. This is a real defect and the fix is retained (it is
**connection 2** in the rationale document).

The measured result after that fix looked conclusive — clean p ≈ 2.0 over the
first three comparisons, then estimates rising to 2.2–2.8, then collapsing.
The collapse was attributed to a **spatial-error floor**, and an N=1280 rerun
was presented as independent confirmation because the plateau moved down ~4.4×
in 2D, matching the theoretical $1/N^2$ prediction.

**That conclusion was wrong**, for a reason that is worth recording.

## 4. Why the supporting evidence was misleading

The refresh fixes the $V_m$ argument of $I_{ion}$, but **not the gating
variables themselves**. `solveODE` still integrated $\mathbf{w}$ across the
whole step with $V_m$ frozen — which is precisely Godunov splitting on the
reaction side, and is $O(\Delta t)$.

Three independent checks established this.

**(a) Analytically.** For the manufactured model, with
$W = u_1 + u_3 - V$ and the exact solution $W = (1+t)G$,
$u_2^2 = G^{-1}(1+t)^{-2}$, $V = \sqrt{1+t}\,F$:

$$
\frac{\partial}{\partial V}\!\left(\frac{du_1}{dt}\right)
  = -\frac{3}{2(1+t)} - \frac{F}{2G(1+t)^{3/2}},
\qquad
\frac{\partial}{\partial V}\!\left(\frac{du_2}{dt}\right) = u_2^3
$$

Both non-zero on the exact solution. Substituting $V^n$ for $V(t)$ therefore
gives local error $O(\Delta t^2)$ with a smooth, non-oscillating coefficient,
which accumulates to **global $O(\Delta t)$** — the textbook splitting error.

**(b) 0-D numerically.** Integrating the manufactured gating ODEs alone with a
high-accuracy RK4 sub-integrator, comparing frozen-$V$ against
linearly-varying-$V$:

| | $e_{u_1}$ order | $e_{u_2}$ order | $e_{I_{ion}}$ order |
|---|---|---|---|
| frozen $V$ (as implemented) | 1.000 | 1.000 | 1.000 |
| linear $V$ (proposed) | 2.000 | 2.000 | 2.000 |

Order 1.000 to three digits, over five refinements.

**(c) A full replica of the discrete scheme.** A 1-D Python replica of the
exact cardiacFoam step (BDF2 diffusion, AB2 source, per-cell gating
integration, same history rotation) reproduces the post-first-fix numbers to
**three significant figures** (4.7467e-4 vs the measured 4.7305e-4;
1.2177e-4 vs 1.2139e-4; 3.0240e-5 vs 3.0082e-5), confirming the replica is
faithful.

In that replica the reference solution is computed **on the same spatial
grid**, so it contains *exactly zero* spatial error — yet it shows the
identical order drift (2.06 → 2.14 → 2.27) that was attributed to a spatial
floor. A spatial floor cannot be the explanation for a drift that survives the
complete removal of spatial error.

**The decisive test needs no reference solution at all.** Comparing the two
scheme variants directly:

| dt | ‖V(first-fix) − V(both-fixes)‖₂ | p |
|---|---|---|
| 0.2/2³ | 1.474041e-04 | — |
| 0.2/2⁴ | 3.720427e-05 | 1.986 |
| 0.2/2⁵ | 8.737859e-06 | 2.090 |
| 0.2/2⁶ | 1.842867e-06 | 2.245 |
| 0.2/2⁷ | 3.580247e-07 | 2.364 |
| 0.2/2⁸ | 1.526090e-07 | 1.230 |
| 0.2/2⁹ | 9.967751e-08 | 0.614 |
| 0.2/2¹⁰ | 5.768409e-08 | 0.789 |
| 0.2/2¹¹ | 3.095187e-08 | 0.898 |

The difference decays at **p → 1**. Two second-order schemes cannot differ at
first order. The first-fix-only scheme therefore carries a genuine $O(\Delta t)$
term that the fixed one does not.

The apparent p ≈ 2 originally measured was a **pre-asymptotic window**: the
$O(\Delta t^2)$ term dominates at coarse dt, and the residual $O(\Delta t)$
term only emerges below dt ≈ 8e-4 — exactly where the original tables showed
the order "collapsing".

## 5. Root cause, confirmed against the literature

The reference scheme is Ethier & Bourgault (2008) eq. (2.14), which is a
**pair** of equations — the gating variable $v$ is advanced by the same BDF2
formula with the same AB2 extrapolation as $u$:

$$
M\,\frac{\tfrac{3}{2}v^{n+1} - 2v^{n} + \tfrac{1}{2}v^{n-1}}{\Delta t}
  = \epsilon\Big(2G(u^{n},v^{n}) - G(u^{n-1},v^{n-1})\Big)
$$

The paper never freezes $u$ during the gating update. cardiacFoam implemented
only the first of the two equations. That missing second equation is the root
cause.

## 6. Fix

Both couplings, as described in the rationale document:

1. **PDE → ODE** (*this was the missing one, and it is order-critical*):
   the ionic ODE integration is handed $dV_m/dt$ and carries $V_m$ along its
   linear extrapolant across the step instead of holding it frozen.
2. **ODE → PDE** (the first diagnosis, retained): $I_{ion}$ is re-evaluated
   algebraically at the just-solved $V_m^{n+1}$ before the history rotation.

The `dt*1.0e-6` hack that originally stood in for the algebraic evaluation has
been **deleted**. It was a magic constant that perturbed the gating state on
every step and relied on that perturbation staying negligible. It is replaced
by a real non-advancing evaluate (`evaluateIonicCurrent`), and models that
cannot provide one now cause a hard error under `sbdf2` rather than silently
degrading.

`evaluateIonicCurrent` duplicates `solveODE`'s algebraic tail, so it was
checked directly for dropped terms: on the first step, where no rate is
supplied and both paths see the same $V_m$, `maxAbsDiff = 0` against
`maxAbsIion = 0.677` — bit-identical.

## 7. Results

Real solver, 1-D, N=1280, endTime 0.2:

| dt | L1 | L2 | L∞ | p (L2) |
|---|---|---|---|---|
| 0.025 | 2.7906e-04 | 3.2760e-04 | 6.0481e-04 | — |
| 0.0125 | 7.1875e-05 | 8.4549e-05 | 1.5653e-04 | 1.954 |
| 0.00625 | 1.8185e-05 | 2.1421e-05 | 3.9730e-05 | 1.981 |
| 0.003125 | 4.5202e-06 | 5.3351e-06 | 9.9220e-06 | 2.005 |
| 0.0015625 | 1.0732e-06 | 1.2755e-06 | 2.3929e-06 | 2.064 |
| 0.00078125 | 2.0763e-07 | 2.5711e-07 | 5.0091e-07 | 2.311 |

Two cross-checks on these numbers:

- They match the Python replica's fixed variant to **four significant figures**
  (3.2768e-4, 8.4634e-5, 2.1506e-5, 5.4203e-6), independently confirming the
  C++ implements the analysed scheme.
- They are uniformly ~30% **below** the first-fix-only N=1280 results
  (4.7316e-4, 1.2150e-4, 3.0194e-5, 7.2248e-6) at every dt.

### 2-D (N=640, 409,600 cells)

| dt | L1 | L2 | L∞ | p (L2) | L2 (first-fix only) | gain |
|---|---|---|---|---|---|---|
| 0.025 | 3.4469e-04 | 4.3568e-04 | 1.1091e-03 | — | 6.7027e-04 | 1.54× |
| 0.0125 | 8.7618e-05 | 1.1096e-04 | 2.8429e-04 | 1.973 | 1.7101e-04 | 1.54× |
| 0.00625 | 2.1504e-05 | 2.7295e-05 | 7.0471e-05 | 2.023 | 4.1873e-05 | 1.53× |
| 0.003125 | 4.7421e-06 | 6.0660e-06 | 1.6044e-05 | 2.170 | 9.3360e-06 | 1.54× |
| 0.0015625 | 5.2542e-07 | 7.5053e-07 | 2.3136e-06 | 3.015 | 1.3359e-06 | 1.78× |
| 0.00078125 | 5.4049e-07 | 6.7043e-07 | 1.5190e-06 | 0.163 | 5.9151e-07 | 0.88× |

Second order confirmed independently in 2-D (p = 1.973, 2.023), with a uniform
**1.54× error reduction** against the first-fix-only scheme at every dt in the
clean region.

#### The 2-D floor is spatial — confirmed

The 2-D floor sits at ~6.7e-7 at N=640. Repeating the fine end of the ladder
at N=1280 (1,638,400 cells) tests the $1/N^2$ prediction directly:

| dt | L2 (N=640) | L2 (N=1280) | ratio |
|---|---|---|---|
| 0.003125 | 6.0660e-06 | 6.6083e-06 | 0.92 |
| 0.0015625 | 7.5053e-07 | 1.2644e-06 | 0.59 |
| 0.00078125 | 6.7043e-07 | 1.4187e-07 | **4.73** |

At the floor point the error drops **4.73×** for a 2× mesh refinement —
matching the theoretical 4× for a second-order spatial scheme. **The 2-D floor
is spatial discretisation error.**

This also vindicates the original document's spatial-floor explanation *for
2-D* — it measured 4.40× on the same test, independently reproduced here as
4.73×. That explanation was only wrong where it was generalised to explain the
order deficit itself, and where the analogous 1-D floor is concerned (below).

The ratios below 1 at coarser dt are the expected sign-cancellation: at N=640
the spatial error partially cancels the temporal error, flattering that
resolution. Removing most of it at N=1280 exposes more of the true temporal
error, which is why N=1280 looks *worse* at dt = 0.003125 and 0.0015625. The
clean-region orders (p = 1.973, 2.023 at coarse dt, where temporal error
dominates) are unaffected by this.

### The 1-D tail: a characterised but unexplained floor

The clean second-order region is **dt = 0.025 → 0.0015625** — four
comparisons, p = 1.954, 1.981, 2.005, 2.064, spanning more than two decades of
error (3.3e-4 down to 1.3e-6). Below that the apparent order rises above 2 and
the error eventually turns around:

| dt | L2 (N=1280) |
|---|---|
| 0.00078125 | 2.5711e-07 |
| 0.000390625 | 3.0158e-08 |
| 0.0001953125 | 7.0925e-08 |

The error passes through a near-zero crossing around dt ≈ 4e-4 and then grows.
That is the signature of a dt-independent error term of opposite sign, roughly
6e-8 in magnitude. **Four candidate explanations were tested; all four were
falsified.**

1. **Spatial discretisation error** (in 1-D). At the re-grown point dt = 0.0001953125:
   N=1280 → 7.0925e-08, N=2560 → 5.7716e-08, N=5120 → 6.0355e-08. Refining
   the mesh 4× does not reduce the floor; it plateaus near 6e-8. Not spatial.
2. **ODE solver tolerance.** Re-running with `absTol 1e-12`, `relTol 1e-10`
   (default `relTol` is 1e-4) reproduces the loose-tolerance result
   *bit-identically* — 3.015790e-08 vs 3.01579e-08. Not the ODE tolerance.
3. **Linear solver tolerance.** Already `tolerance 1e-15`, `relTol 0` with
   PCG/DIC. Not the linear solve.
4. **A spatial floor moving with N as 1/N²** — the explanation given in the
   original version of this document. Falsified by (1) *in 1-D*. Note this is
   the one place where 1-D and 2-D genuinely differ: the same test in 2-D
   gives a clean 4.73× drop and confirms a spatial floor there.

So a dt-, N-, and tolerance-independent error of ~6e-8 remains, and **its
origin is not known**. It is recorded here as an open question rather than
given a plausible-sounding explanation, because assigning one prematurely is
exactly the mistake that let the incomplete first fix stand.

An earlier reading of this data — that the N=1280-vs-N=2560 agreement at
dt = 0.003125/0.0015625/0.00078125 proved the drift was non-spatial — was
itself over-read: those points sit at or before the crossing, where the
comparison reflects the crossing position rather than the floor magnitude. The
discriminating comparison is the one in (1), at the re-grown point.

None of this affects the second-order claim. The floor sits at 6e-8, four
orders of magnitude below where the clean region begins, and the Python
replica — which has no spatial error, no ODE-solver tolerance and no linear
solve — converges to p = 2.000 at every refinement.

## 8. Scope

- Godunov is unchanged and remains the first-order control: with no rate
  supplied `activeVmRate()` is zero, reproducing the frozen $V_m$ exactly, and
  the refresh is skipped.
- No prior `godunov`-based MMS result is affected.
- Bidomain shares `myocardiumDomain::advance()` and so inherits both
  connections, but has not been exercised by a convergence sweep here.
- The 11 non-batched CellML models do not yet implement the two capabilities
  and will hard-error under `sbdf2`; their batched twins work.

## 9. Run note — sweeps silently resume stale results

`driverFoam sweep-run` resumes from an existing `workflow_state.json`. If a
case directory left over from an earlier session reports `completed`, the
sweep re-reports it as `completed` / `fresh`, **runs no solver, and exits 0**.
There is no warning.

This nearly invalidated the 2-D check here: the first 2-D sweep after the fix
returned the previous session's numbers, and the only tell was that all six
values were bit-identical to the pre-fix table to six significant figures. Two
different schemes cannot agree to six digits — identical numbers mean identical
execution, not agreement.

**Any sweep re-run after a code change must `rm -rf` the case directories
first**, otherwise it reports the old code's results as if they were new. This
applies to every before/after verification in this repo, not just SBDF2.

## 10. Build note

Adding the capability virtuals changes vtable layouts across library
boundaries. A plain `Allwmake` is **not sufficient** — it reported "no build
errors" while producing a `libelectroModels.dylib` that exported
`electroModel::read()`'s thunk at base offset 384 while two of its own objects
referenced offset 368, so the solver died at load with
`symbol not found in flat namespace '__ZThn368_N4Foam12electroModel4readEv'`.
`electroModel.o` was stale and wmake's dependency tracking did not catch the
header change. An explicit `wclean` of the affected libraries before
`Allwmake` is required.
