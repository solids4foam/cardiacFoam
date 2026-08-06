# SBDF2 (Semi-Implicit BDF2) Time-Stepping in cardiacFoam

## 1. The problem being discretised

The monodomain model couples a parabolic PDE for the transmembrane potential
$V_m$ to a system of ODEs for the cell state $\mathbf{w}$ (gating variables,
concentrations):

$$
\chi C_m \frac{\partial V_m}{\partial t}
  = \nabla\!\cdot\!(\boldsymbol{\sigma}\nabla V_m)
  - \chi C_m\, I_{ion}(V_m,\mathbf{w})
  + S
\qquad\text{(PDE)}
$$

$$
\frac{d\mathbf{w}}{dt} = \mathbf{g}(V_m,\mathbf{w})
\qquad\text{(ODE)}
$$

The coupling runs **both ways**: $I_{ion}$ carries the cell state into the PDE,
and $V_m$ enters the ODE right-hand side. Any time-stepping scheme has to be
second-order accurate in *both* directions to be second order overall. This is
the single most important fact in this document, and it is what the original
implementation got wrong.

The bidomain case adds the elliptic $\phi_e$ equation but the reaction coupling
is structurally identical, so everything below applies unchanged.

## 2. The reference scheme

The scheme is **SBDF2** — second-order semi-implicit backward differentiation,
an IMEX method: diffusion implicit, reaction explicit with second-order
extrapolation.

The canonical statement is Ethier & Bourgault (2008), equation (2.14). Written
in their notation ($u$ = transmembrane potential, $v$ = ionic state variable,
$F$ = ionic current term, $G$ = gating right-hand side, $A_i$ = the diffusion
operator, $M$ = mass matrix):

$$
M\,\frac{\tfrac{3}{2}u^{n+1} - 2u^{n} + \tfrac{1}{2}u^{n-1}}{\Delta t}
  = \frac{1}{\epsilon}\Big(2F(u^{n},v^{n}) - F(u^{n-1},v^{n-1})\Big)
  - A_i\big(u^{n+1} + u_e^{n+1}\big)
$$

$$
M\,\frac{\tfrac{3}{2}v^{n+1} - 2v^{n} + \tfrac{1}{2}v^{n-1}}{\Delta t}
  = \epsilon\Big(2G(u^{n},v^{n}) - G(u^{n-1},v^{n-1})\Big)
$$

**SBDF2 is a pair of equations, not one.** Both unknowns are advanced by the
same BDF2 formula, and both explicit terms use the same AB2 extrapolation
$2(\cdot)^n - (\cdot)^{n-1}$. The gating variable $v$ is *not* obtained by
sub-integrating with $u$ held at $u^n$ — the paper never freezes $u$.

Compare the implicit Gear scheme, their eq. (2.15), which differs only in
evaluating $F$ and $G$ at $(u^{n+1}, v^{n+1})$ and therefore requires a
nonlinear solve. SBDF2 exists precisely to avoid that Newton iteration while
keeping second order.

## 3. What cardiacFoam actually does, and why

cardiacFoam cannot implement the second equation literally. Its architecture
runs a **per-cell adaptive stiff ODE integrator** (RKF45 and friends, with
substepping) for the cell model, which is far more robust for stiff ionic
models than a single BDF2 step and is shared with the single-cell solver.
Replacing it with a one-step BDF2/AB2 update would be a rewrite and would lose
that robustness.

So the second equation is realised *equivalently* rather than literally. The
requirement it encodes is:

> over $[t^n, t^{n+1}]$, the gating integration must see a **second-order
> accurate** $V_m$.

The natural device — and the same one already used for the source term — is
**linear extrapolation from the two known voltage levels**:

$$
V_m(t) \;\approx\; V_m^{n} + (t - t^{n})\,
  \frac{V_m^{n} - V_m^{n-1}}{\Delta t_{n-1}},
\qquad t \in [t^{n},\,t^{n+1}]
$$

### Why linear, and why this is the right order

- **It is the matching order.** The extrapolant's error is
  $O(\Delta t^2)$ uniformly on the step. Feeding an $O(\Delta t^2)$-accurate
  $V_m$ into the ODE gives local error $O(\Delta t^3)$, hence global
  $O(\Delta t^2)$ — exactly what BDF2 needs, and no more. A higher-order
  extrapolant would be wasted work; a constant (frozen $V_m$) is $O(\Delta t)$
  and destroys the order.
- **It is literally the same device as the source extrapolation.** The AB2
  term $2I^n - I^{n-1}$ *is* linear extrapolation of $I_{ion}$ through the two
  known levels, evaluated at $t^{n+1}$. Extrapolating $V_m$ the same way is
  that identical idea applied to the other direction of the coupling. The two
  connections are symmetric.
- **It costs nothing.** Every CellML-generated model already pins
  `RATES[V] = 0.0` when `solveVmWithinODESolver` is false. Writing the
  extrapolant slope into that row instead makes the existing adaptive
  integrator carry $V_m$ along the extrapolant. The slope is constant over the
  step, so that row is integrated *exactly* — no accuracy is lost in the
  device itself.

### Variable time steps

Both extrapolations generalise with the step ratio
$r = \Delta t_n / \Delta t_{n-1}$:

$$
I^{\text{extrap}} = I^{n} + r\,(I^{n} - I^{n-1}),
\qquad
\left.\frac{dV_m}{dt}\right|_{\text{step }n} = \frac{V_m^n - V_m^{n-1}}{\Delta t_{n-1}}
$$

For constant $\Delta t$, $r = 1$ and the source reduces to $2I^n - I^{n-1}$,
matching eq. (2.14).

## 4. The two connections in code

| # | Direction | Requirement | Implementation |
|---|---|---|---|
| 1 | PDE → ODE | Gating integrated against a 2nd-order-accurate $V_m$ | `myocardiumDomain` computes `VmRate_ = (Vm_ - VmPrev_)/deltaT0` and hands it to the ionic model via `setVmRate()`; models return it as the $V_m$ row of their derivative vector |
| 2 | ODE → PDE | $I^n$ in the extrapolation must be $I_{ion}$ at time level $n$ | after the diffusion solve produces $V_m^{n+1}$, `evaluateIonicCurrent()` recomputes $I_{ion}$ algebraically at that $V_m$ without moving the gating state |

**Connection 1 is order-critical.** Without it the scheme is asymptotically
first order no matter what is done to the source term — freezing $V_m$ during
the reaction step *is* Godunov splitting, and it is what caps Godunov at
$O(\Delta t)$ in the first place.

**Connection 2 is not order-critical once connection 1 is in place** (the
residual $V_m$ error in $I^n$ is then already $O(\Delta t^2)$), but it is kept
because it makes `IionOld_` genuinely be $I_{ion}$ at time level $n$ rather
than at some extrapolated voltage, and its correctness does not depend on how
good the extrapolant happened to be. It is also what the exported `Iion` field
should mean.

### Ordering within a step

```
prepareTimeStep():   IionOldOld_ = IionOld_;  IionOld_ = Iion_
                     VmRate_ = (Vm_ - VmPrev_)/deltaT0;  VmPrev_ = Vm_
advance():           setVmRate(VmRate_)
                     solveIonicCurrent(t0, dt)          -> gating to t^{n+1}
                     clearVmRate()
                     solveDiffusionImplicit(...)        -> Vm(t^{n+1})
                     refreshIonicCurrent(t0 + dt)       -> Iion at (w^{n+1}, Vm^{n+1})
```

The refresh must land in `Iion_` *before* the next step's `prepareTimeStep`
rotates it into `IionOld_`; that rotation is what feeds the extrapolation.

### Startup

`IionOldOld_` only holds a genuinely computed value from the third step
(`timeIndex() >= 3`), so the source extrapolation falls back to the
frozen-endpoint value before then. The $V_m$ rate needs only two levels and so
switches on at `timeIndex() >= 2`. Both fallbacks affect $O(1)$ steps and
contribute $O(\Delta t^2)$ globally, which does not degrade the order.

## 5. Model support and the capability gate

Both connections need something from the ionic model that the base `ionicModel`
interface did not previously expose:

- `supportsVmRateCoupling()` — does `solveODE` honour a supplied $dV_m/dt$?
- `supportsIonicCurrentEvaluation()` / `evaluateIonicCurrent()` — can $I_{ion}$
  be recomputed at a given $V_m$ *without advancing the state*?

Both default to `false`/`NotImplemented`, so unmodified models keep exactly
their historical frozen-$V_m$ behaviour. `myocardiumDomain`'s constructor
**hard-fails** if `timeCouplingScheme sbdf2` is selected with a model that does
not implement both, naming the missing capability. There is deliberately no
silent fallback: a scheme that quietly degrades to first order while being
labelled SBDF2 is what caused this whole investigation.

Currently implemented in:

- `monodomainFDAManufactured` — the MMS verification model.
- `batchedIonicModel` — written once in the shared base, so all 12 batched
  models (`TNNPBatched`, `ToRORd_dynClBatched`, `BuenoOrovioBatched`, …) get
  both connections with no per-model changes.

The 11 non-batched CellML models still return `false` and will error under
`sbdf2`; use `godunov`, or their batched twins. Extending them is one file at a
time and needs no new machinery — each already calls
`<Model>computeVariables(...)` internally, which is exactly the non-advancing
evaluate, and already pins `RATES[V] = 0.0`, which is exactly the row the rate
goes into.

### Verification of the evaluate hook

`evaluateIonicCurrent` duplicates the algebraic tail of `solveODE` and so could
silently drop a term (the `/CONSTANTS_[Cm]` scaling, the manufactured-source
correction). This was checked directly rather than assumed: on the first time
step, where no $V_m$ rate is supplied and `solveODE` therefore ends with
`S[V] == Vm[i]`, the two paths must agree exactly. Measured
`maxAbsDiff = 0` against `maxAbsIion = 0.677` — bit-identical, not merely
machine precision.

## 6. Required case settings

**`constant/electroProperties`**

```c++
electrophysicsAdvanceScheme staggeredElectrophysicsAdvanceScheme;

monodomainSolverCoeffs
{
    timeCouplingScheme  sbdf2;   // default: godunov
}
```

**`system/fvSchemes`**

```c++
ddtSchemes
{
    default  backward;   // BDF2; Euler here drops the scheme to first order
}
```

`backward` is genuine textbook BDF2 for constant $\Delta t$ beyond the first
step, $(3V^{n+1} - 4V^n + V^{n-1})/(2\Delta t)$ — this was traced through
OpenFOAM's `backwardDdtScheme.C` and confirmed.

## 7. Verified order

1-D manufactured monodomain (`monodomainPseudoECG`), diagonal conductivity,
hex mesh, $N = 1280$, endTime $0.2$:

| $\Delta t$ | $L_1$ | $L_2$ | $L_\infty$ | $p$ ($L_2$) |
|---|---|---|---|---|
| 0.025 | 2.7906e-04 | 3.2760e-04 | 6.0481e-04 | — |
| 0.0125 | 7.1875e-05 | 8.4549e-05 | 1.5653e-04 | 1.954 |
| 0.00625 | 1.8185e-05 | 2.1421e-05 | 3.9730e-05 | 1.981 |
| 0.003125 | 4.5202e-06 | 5.3351e-06 | 9.9220e-06 | 2.005 |
| 0.0015625 | 1.0732e-06 | 1.2755e-06 | 2.3929e-06 | 2.064 |
| 0.00078125 | 2.0763e-07 | 2.5711e-07 | 5.0091e-07 | 2.311 |

2-D, N=640 (409,600 cells), same case and dt ladder:

| $\Delta t$ | $L_1$ | $L_2$ | $L_\infty$ | $p$ ($L_2$) |
|---|---|---|---|---|
| 0.025 | 3.4469e-04 | 4.3568e-04 | 1.1091e-03 | — |
| 0.0125 | 8.7618e-05 | 1.1096e-04 | 2.8429e-04 | 1.973 |
| 0.00625 | 2.1504e-05 | 2.7295e-05 | 7.0471e-05 | 2.023 |
| 0.003125 | 4.7421e-06 | 6.0660e-06 | 1.6044e-05 | 2.170 |
| 0.0015625 | 5.2542e-07 | 7.5053e-07 | 2.3136e-06 | 3.015 |
| 0.00078125 | 5.4049e-07 | 6.7043e-07 | 1.5190e-06 | 0.163 |

Second order is confirmed independently in both dimensions, with errors a
uniform 1.5× (2-D) and 1.3× (1-D) below the pre-fix scheme at every dt in the
clean region.

Both ladders meet a floor at fine $\Delta t$, but of different origin:

- **2-D (~6.7e-7 at N=640): spatial discretisation error, confirmed.**
  Doubling to N=1280 drops the floor 4.73×, matching the theoretical 4× for a
  second-order scheme. Refine the mesh to push it down.
- **1-D (~6e-8 at N=1280): origin unknown.** Spatial resolution (tested to
  N=5120), ODE tolerance and linear-solver tolerance were each tested and
  ruled out. Recorded as an open question in the finding document.

For order demonstrations use dt ≥ 1.5e-3 (1-D at N=1280) or dt ≥ 6e-3
(2-D at N=640); refining N extends the clean range in 2-D.

Godunov remains the first-order control and is unchanged by this work: with no
rate supplied `activeVmRate()` is zero, which reproduces the historical frozen
$V_m$ exactly, and the refresh is skipped.

See [SBDF2_ORDER_REDUCTION_FINDING.md](SBDF2_ORDER_REDUCTION_FINDING.md) for
the investigation that established connection 1 was missing, including the
independent evidence chain.

## 8. References

1. **Ethier, M., & Bourgault, Y. (2008).** *Semi-implicit time-discretization
   schemes for the bidomain model.* SIAM Journal on Numerical Analysis,
   46(5), 2443–2468. DOI 10.1137/070680503.
   [PDF](https://mysite.science.uottawa.ca/ybourg/publications/SINUM_EthierBourgault08.pdf)
   — eq. (2.14) is the SBDF2 scheme quoted in §2; eq. (2.15) is the implicit
   Gear scheme it avoids; §3.2.2 gives the stability analysis.

2. **Pathmanathan, P., et al. (2010).** *Verification of solver algorithms for
   the cardiac monodomain equation.* International Journal for Numerical
   Methods in Biomedical Engineering, 26(9), 1076–1089. — on how the treatment
   of the ionic current integration governs the observed temporal order.

3. **Roy, T., Bourgault, Y., & Pierre, C. (2020).** *Analysis of time-stepping
   methods for the monodomain model.* Computational and Applied Mathematics,
   39, 230.
   [link](https://link.springer.com/article/10.1007/s40314-020-01254-z)

4. **Cervi, J., & Spiteri, R. J. (2018).** *High-order operator splitting for
   the bidomain and monodomain models.* SIAM Journal on Scientific Computing,
   40(2), A769–A786.
   [link](https://epubs.siam.org/doi/10.1137/17M1137061)
