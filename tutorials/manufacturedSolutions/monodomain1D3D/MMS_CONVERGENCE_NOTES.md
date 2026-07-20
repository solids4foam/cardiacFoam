# MMS Convergence Verification — Coupled 1D-3D Monodomain

## Overview

Manufactured Method of Solutions (MMS) for the coupled 1D Purkinje / 3D myocardium
monodomain system.

**Manufactured solutions**

| Domain | V_exact | F(x,y,z) |
|--------|---------|-----------|
| 1D graph | `√(1+t) · cos(πx)` | dim=1: x along branch |
| 3D myocardium | `√(1+t) · cos(πx)cos(2πy)cos(3πz)` | dim=3 |

Ionic gates: `u1_exact = u2_exact = V_exact/2`, with corresponding forcing terms.

**Current graph placement.** The checked-in graph family now uses
`y=1/6, z=1/3`, with PVJ terminals at `(0, 1/6, 1/3)` and `(1, 1/6, 1/3)`.
These terminals are still on the `x=0` and `x=1` tissue boundary faces, where
the manufactured 3D normal flux is zero.  The old placement
`y=0.5, z=1/3` made `cos(2πy)cos(3πz)=1`, so `F_3D(terminal)=F_1D`.  The new
placement gives `cos(π/3)cos(π)=-0.5`, so the terminal coupling terms do not
cancel pointwise.

**Mesh pairing** — `h_1D ≈ h_3D` throughout the sweep:

| N³ mesh | Graph nodes | h_3D | h_1D |
|---------|-------------|------|------|
| 10³ | 11 | 0.100 | 0.050 |
| 20³ | 21 | 0.050 | 0.025 |
| 40³ | 41 | 0.025 | 0.013 |
| 80³ | 81 | 0.0125 | 0.006 |

**Time-step scaling** — `dt ~ h²` (anchored at N=80, dt=1.4×10⁻⁴):

| N | dt |
|---|----|
| 10 | 8.97×10⁻³ |
| 20 | 2.24×10⁻³ |
| 40 | 5.61×10⁻⁴ |
| 80 | 1.40×10⁻⁴ |

This matches the standalone MMS driver and ensures temporal error ≪ h².

---

## 1D standalone — O(h²) confirmed

The 1D graph solver converges at exactly second order across all refinements.
The 1D result is **independent of coupling** and serves as a reference throughout.

---

## 3D standalone — O(h²) confirmed

With coupling disabled the 3D implicit monodomain also converges at second order,
matching the standalone `monodomainPseudoECG` driver results.

---

## Historical coupled sweep — rPvj=1, pvjRadius=0.11

### Setup

- PVJ terminals at `(0, 0.5, 1/3)` and `(1, 0.5, 1/3)` — **on the x=0/x=1 boundary faces**
- MMS consistency: `F_3D(terminal) = cos(πx)·cos(π)·cos(π) = cos(πx) = F_1D(x)` ✓
- `rPvj = 1.0`, `pvjRadius = 0.11`, linear kernel

### Bug fixed during this work — V_1D uninitialized at first coupling step

In `staggeredElectrophysicsAdvanceScheme::advance`, the order was:

```
prepareConductionCouplings(t0, dt)   ← reads V_1D = 0 at t=0  ← SPIKE
advanceConductionDomains(t0, dt)     ← preProcess() sets V_1D = V_exact at t=0
```

Fix: added `preInitialize()` virtual to `electroDomainInterface`, overridden in
`conductionSystemDomain` to call `preProcess()`, and called
`system.preInitializeConductionDomains()` before `prepareConductionCouplings` when
`t0 == 0`.  This reduced the first-step coupling source from ±1333 → ±111 (the
remaining level was the uncancelled PVJ residual at the manufactured reference,
not uninitialized state).

### Results at t=0.1

| N | h | Linf_3D_Vm | rate_3D | Linf_1D_Vm | rate_1D | coupling_srcL1 |
|---|---|---|---|---|---|---|
| 10 | 0.100 | 8.89e-02 | — | 4.23e-03 | — | 4.55e-03 |
| 20 | 0.050 | 1.68e-01 | **-0.92** | 1.06e-03 | **2.00** | 1.24e-02 |
| 40 | 0.025 | 1.81e-01 | **-0.10** | 2.65e-04 | **2.00** | 1.31e-02 |
| 80 | 0.013 | 1.84e-01 | **-0.03** | 6.62e-05 | **2.00** | 1.34e-02 |

### Results at t=0.02 vs t=0.1 — coupling accumulation over time

| N | Linf_3D (t=0.02) | Linf_3D (t=0.1) | Δ |
|---|---|---|---|
| 10 | 0.093 | 0.089 | −4% |
| 20 | 0.161 | 0.168 | +4% |
| 40 | 0.176 | 0.181 | +3% |
| 80 | 0.180 | 0.184 | +2% |

The 3D error barely grows from t=0.02 to t=0.1.  The pre-fix coupling-source
contamination was a fixed PVJ residual at the manufactured reference, not a
temporal accumulation effect.

### Historical root cause framing — sphere-averaging mismatch at boundary terminal

The volumetric coupling gathers a volume-weighted average of V_3D over all cells
within the sphere of radius R:

```
V_3D_avg = Σ( kernel(r) · V_3D[cell] · cellVol ) / Σ( kernel(r) · cellVol )
```

**At the exact solution this average ≠ V_1D_terminal** because:

1. The terminals sit on the x=0 / x=1 **boundary faces**.  The sphere can only
   sample cells with x > 0, so the one-sided average is always less than the
   pointwise value at x=0.
2. The 3D manufactured function `cos(πx)cos(2πy)cos(3πz)` varies spatially over
   the sphere — the average over a finite volume is not the value at the centre.

That observation was useful but incomplete. The fixed-radius `pvjMapper` average
does not create an unbounded point-source error under refinement: the gather and
deposit operators converge to well-defined, bounded continuum quantities for the
chosen radius and kernel. The evidence is the coupling source plateau
(`coupling_srcL1 ≈ 1.3e-2` by N=80), not growth without bound.

The standalone manufactured ionic model cancels the monodomain PDE residual. The
coupled MMS additionally adds the exact PVJ source to that manufactured ionic
residual, so the analytical 1D/3D fields satisfy the coupled equations. Active
coupling itself is left raw: the verifier does not alter terminal currents, 3D
source fields, or implicit source coefficients. The coupled verifier also
computes the exact PVJ source after the solve and reports the diagnostic residual
`S_pvj(V_num) - S_pvj(V_exact)`.

For the explicit source path, the diagnostic compares the deposited terminal
current source against the exact terminal current source. For
`pvjCouplingScheme implicit`, the production operator is split into a
network-voltage source and tissue-voltage diagonal sink, so the diagnostic forms
the effective source `source - coeff*Vm` before comparing against the exact
split source. The production coupler does not contain FDA-specific formulas; it
only calls the generic verification hook that supplies the manufactured ionic
source term.

---

## Historical corrected coupled sweep — unidirectional, rPvj=1, pvjRadius=0.11, t=0.1

Configuration: `couplingMode unidirectional`, `rPvj=1.0`, `pvjRadius=0.11`,
linear kernel, terminals on the x-boundary faces using the older cancelling
`y=0.5, z=1/3` placement.

| N | h | Linf_3D_Vm | rate_3D | Linf_1D_Vm | rate_1D | coupling_srcL1 | rate_coupling |
|---|---|---|---|---|---|---|---|
| 10 | 0.1000 | 4.9598e-03 | — | 3.4317e-04 | — | 1.5884e-03 | — |
| 20 | 0.0500 | 1.3812e-03 | **1.8444** | 8.8475e-05 | **1.9556** | 4.1188e-04 | **1.9472** |
| 40 | 0.0250 | 3.4859e-04 | **1.9863** | 2.1956e-05 | **2.0106** | 1.0324e-04 | **1.9962** |
| 80 | 0.0125 | 8.7734e-05 | **1.9903** | 5.5004e-06 | **1.9970** | 2.5993e-05 | **1.9898** |

This table came from an earlier verifier-side manufactured PVJ forcing
experiment. It is kept as historical context only; it should not be used as raw
coupling validation.

---

## Historical corrected coupled sweep — bidirectional, rPvj=1, pvjRadius=0.11, t=0.1

Same case, with `couplingMode bidirectional` to exercise the 1D deposit path that
is cleared in unidirectional mode.

| N | h | Linf_3D_Vm | rate_3D | Linf_1D_Vm | rate_1D | coupling_srcL1 | rate_coupling |
|---|---|---|---|---|---|---|---|
| 10 | 0.1000 | 4.9598e-03 | — | 3.3428e-04 | — | 1.5853e-03 | — |
| 20 | 0.0500 | 1.3812e-03 | **1.8444** | 8.7047e-05 | **1.9412** | 4.1135e-04 | **1.9463** |
| 40 | 0.0250 | 3.4859e-04 | **1.9863** | 2.1777e-05 | **1.9990** | 1.0317e-04 | **1.9954** |
| 80 | 0.0125 | 8.7734e-05 | **1.9903** | 5.4749e-06 | **1.9919** | 2.5984e-05 | **1.9893** |

This table uses the same historical verifier-side manufactured PVJ forcing and
exercises the retrograde 1D deposit path as well as the tissue source path.

### What this test does and doesn't verify

The raw coupled MMS validates the explicit/implicit PVJ schemes as they are
implemented. The verifier writes the exact PVJ source magnitude and the residual
source norms so the source behaviour can be checked directly without changing
the solve.

---

## Coupled sweep — rPvj=1e6 (coupling negligible), t=0.1

Setting rPvj=10⁶ reduces the coupling source to ~2×10⁻⁷ (noise).
Both 1D and 3D recover O(h²) independently:

| N | h | Linf_3D_Vm | rate_3D | Linf_1D_Vm | rate_1D | coupling_srcL1 |
|---|---|---|---|---|---|---|
| 10 | 0.100 | 4.96e-03 | — | 4.23e-03 | — | 1.7e-07 |
| 20 | 0.050 | 1.38e-03 | **1.84** | 1.06e-03 | **2.00** | 2.3e-07 |
| 40 | 0.025 | 3.50e-04 | **1.98** | 2.65e-04 | **2.00** | 2.2e-07 |
| 80 | 0.013 | 8.91e-05 | **1.97** | 6.62e-05 | **2.00** | 2.2e-07 |

The N=10→20 rate of 1.84 is a known coarse-mesh effect; asymptotic rate is ~2.

---

## Coupled sweep — pvjRadius=0.055 (half sphere), rPvj=1.0, t=0.1

N=10 is skipped (cell size h=0.1 > pvjRadius=0.055).
Halving the sphere reduced the pre-fix absolute error level ~4× compared to
R=0.11, but 3D Linf still plateaued because the PVJ residual at the manufactured
reference was still uncancelled:

| N | h | Linf_3D_Vm | rate_3D | Linf_1D_Vm | rate_1D | coupling_srcL1 |
|---|---|---|---|---|---|---|
| 20 | 0.050 | 2.28e-02 | — | 1.06e-03 | — | 7.7e-04 |
| 40 | 0.025 | 4.07e-02 | **-0.84** | 2.65e-04 | **2.00** | 1.4e-03 |
| 80 | 0.013 | 4.27e-02 | **-0.07** | 6.62e-05 | **2.00** | 1.5e-03 |

---

## Coupled MMS scheme matrix — non-cancelling graph, t=0.1

The current coupled MMS uses the non-cancelling graph placement
`y=1/6, z=1/3` and leaves the production PVJ source assembly raw. The exact PVJ
source is added only through the manufactured ionic residual.

For `rPvj=1.0`, `pvjRadius=0.11`, and N=10→20→40→80:

| Case | rate_3D_Vm | rate_1D_Vm | rate_sourceError |
|------|------------|------------|------------------|
| unidirectional, PVJ explicit | 1.8450, 1.9863, 1.9903 | 1.9556, 2.0106, 1.9970 | 2.0806, 2.0233, 2.0015 |
| unidirectional, PVJ implicit | 1.8444, 1.9863, 1.9903 | 1.9556, 2.0106, 1.9970 | 1.8104, 1.9328, 1.9897 |
| bidirectional, PVJ explicit | 1.8450, 1.9863, 1.9903 | 1.9613, 2.0156, 2.0001 | 2.0808, 2.0235, 2.0016 |
| bidirectional, PVJ implicit | 1.8444, 1.9863, 1.9903 | 1.9557, 2.0162, 2.0005 | 1.8111, 1.9318, 1.9893 |

The source magnitude is not expected to converge to zero; it is the applied PVJ
source. The relevant coupling metric is `coupling_sourceErrorL1`, which converges
at approximately second order.

The final 80³ values are:

| Case | Linf_3D_Vm | Linf_1D_Vm | coupling_sourceErrorL1 |
|------|------------|------------|------------------------|
| unidirectional, PVJ explicit | 8.77337e-05 | 5.50037e-06 | 4.81068e-05 |
| unidirectional, PVJ implicit | 8.77337e-05 | 5.50037e-06 | 1.49776e-05 |
| bidirectional, PVJ explicit | 8.77337e-05 | 5.51280e-06 | 4.81107e-05 |
| bidirectional, PVJ implicit | 8.77337e-05 | 5.51498e-06 | 1.49728e-05 |

A stronger but still clean stress point is `rPvj=0.3`, which raises the N=40
exact source magnitude from about 3.03 to about 10.11 while preserving O(h²)
field convergence:

| Case | rate_3D_Vm | rate_1D_Vm | rate_sourceError |
|------|------------|------------|------------------|
| unidirectional, PVJ explicit | 1.8453, 1.9864 | 1.9556, 2.0106 | 2.0593, 2.0167 |
| unidirectional, PVJ implicit | 1.8443, 1.9863 | 1.9556, 2.0106 | 1.7300, 1.8740 |
| bidirectional, PVJ explicit | 1.8453, 1.9864 | 1.9657, 2.0169 | 2.0594, 2.0168 |
| bidirectional, PVJ implicit | 1.8443, 1.9863 | 1.9543, 2.0178 | 1.7308, 1.8730 |

Lowering to `rPvj=0.1` is useful as a stability stress test, but the explicit
PVJ scheme is too strong on the 10³ mesh and the coarse point is not part of the
recommended convergence evidence.

---

## Summary

| Test | 3D rate | 1D rate | Conclusion |
|------|---------|---------|------------|
| 3D standalone | **~2** | — | solver correct |
| 1D standalone | — | **~2** | solver correct |
| Coupled, rPvj=1, R=0.11, pre-fix | plateau/diverges | **~2** | standalone MMS forcing missed PVJ residual |
| Coupled, rPvj=1, R=0.11, MMS-forced unidirectional | **~2** | **~2** | historical forced-source experiment |
| Coupled, rPvj=1, R=0.11, MMS-forced bidirectional | **~2** | **~2** | historical forced-source experiment |
| Current coupled scheme matrix, rPvj=1, R=0.11 | **~2** | **~2** | explicit/implicit and uni/bi converge |
| Current coupled scheme matrix, rPvj=0.3, R=0.11 | **~2** | **~2** | stronger coupling still converges |
| Coupled, rPvj=1e6, R=0.11 | **~2** | **~2** | coupling negligible → both solvers verified |
| Coupled, rPvj=1, R=0.055, pre-fix | plateau/diverges (lower amplitude) | **~2** | smaller sphere reduced but did not cancel the missing PVJ forcing |

**Key takeaways**

- The 1D and 3D solvers individually converge at O(h²).
- Active coupling adds a non-vanishing PVJ source for the non-cancelling graph
  placement.
- The coupled MMS verifier now supplies the exact PVJ source as a manufactured
  ionic term and reports the exact PVJ source and residual source; it does not
  subtract anything from the production coupling source.
- The exact source is computed in the verification layer using the same PVJ
  radius and kernel as the production mapper; the production coupler does not
  depend on the manufactured FDA reference.
