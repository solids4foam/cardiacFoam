# Final Bath-Bidomain Interface Formulation

## Purpose

This document is the authoritative technical summary for the conformal
tetrahedral bath-bidomain verification used in Paper I. Historical hypotheses,
failed controls, and gate-by-gate evidence remain in
`setup/interfaceStudy/STAGE_REPORT.md`.

## Physical interface problem

Let `H` denote myocardium and `B` the bath. The interface must satisfy

```text
phiE_H = phiB
n . Ge grad(phiE_H) = n . sigmaB grad(phiB)
n . Gi grad(Vm + phiE_H) = 0
```

The last identity is intracellular no-leak. The extracellular heart-side
coefficient is `Ge`, never `Gi + Ge`.

## Selected discrete formulation

The selected runtime configuration is

```text
intracellularAssembly               matchedSubmesh;
interfaceConductivityInterpolation  distanceWeightedHarmonic;
```

The global conformal mesh owns one extracellular potential. Its extracellular
operator uses `Ge` in heart cells and `sigmaB` in bath cells. Across a material
face, the normal two-point resistance is represented by

```text
sigma_f = (d_H + d_B)/(d_H/sigma_H + d_B/sigma_B).
```

The intracellular contribution is assembled on the myocardium submesh. The
same submesh operator is used for both intracellular terms:

```text
L_e^(H+B) phiE + P_H^T L_i^H P_H phiE = -P_H^T L_i^H Vm.
```

`P_H` restricts the global field to heart cells and `P_H^T` maps the heart
matrix into the global system. This preserves the discrete identity

```text
L_i^H(Vm + phiE) = L_i^H(Vm) + L_i^H(phiE)
```

and imposes intracellular zero flux once on the artificial heart-submesh
boundary. Internal and processor-coupled matrix coefficients are mapped, so
the formulation is both implicit and parallel-safe.

The library retains `currentSplit` as its backward-compatible default. This
paper tutorial selects `matchedSubmesh` explicitly.

## Why harmonic interpolation remains necessary

Harmonic interpolation was not the dominant cause of the original N=40
failure, but it is not redundant. The coefficient test proves that
distance-weighted harmonic interpolation reproduces exact series resistance
for unequal owner/face distances. The naive `Gi + Ge` heart coefficient is
physically wrong. Weighted and unweighted harmonic formulas happened to give
similar errors on this mesh family; consistent intracellular assembly produced
the major improvement.

## Verification configuration

- Meshes: conformal tetrahedral N=10, 20, 40, 80.
- N=80 mesh: 6,831,714 cells and 29,528 internal material-interface faces.
- Spatial schemes: least-squares gradients and `Gauss linear corrected`
  Laplacians.
- Time step scales approximately with `h^2`.
- Linear solves use zero relative tolerance. At N=80, `phiE` reaches the
  1000-iteration cap with final algebraic residual around `3e-14`, far below
  the reported discretization errors.
- N=80 ran on six local MPI ranks and completed in approximately 91 minutes.
- Same-mesh serial/parallel comparisons agree exactly across all 109 metrics at
  N=10, 20, and 40.

## Final tetrahedral results

The fitted orders use all four nominal resolutions.

| Metric | N=10 | N=20 | N=40 | N=80 | fitted order |
| --- | ---: | ---: | ---: | ---: | ---: |
| heart `phiE` L2 | 6.066e-3 | 1.812e-3 | 5.165e-4 | 1.630e-4 | 1.75 |
| bath `phiE` L2 | 8.214e-3 | 2.445e-3 | 7.086e-4 | 2.241e-4 | 1.74 |
| x=0 reconstructed heart flux L2 | 3.515e-3 | 4.607e-4 | 8.478e-5 | 1.388e-4 | 1.64 |
| x=1 reconstructed heart flux L2 | 3.751e-3 | 4.564e-4 | 8.132e-5 | 1.348e-4 | 1.69 |
| x=0 reconstructed flux jump L2 | 3.558e-3 | 5.493e-4 | 1.191e-4 | 9.040e-5 | 1.81 |
| x=1 reconstructed flux jump L2 | 3.796e-3 | 5.431e-4 | 1.209e-4 | 9.137e-5 | 1.83 |
| x=0 intracellular leakage L2 | 3.509e-3 | 5.345e-4 | 1.327e-4 | 5.368e-5 | 2.01 |
| x=1 intracellular leakage L2 | 3.752e-3 | 5.272e-4 | 1.355e-4 | 5.613e-5 | 2.01 |
| x=0 assembled flux L2 | 4.680e-4 | 2.942e-4 | 3.444e-4 | 2.875e-4 | 0.19 |
| x=1 assembled flux L2 | 4.612e-4 | 3.027e-4 | 3.510e-4 | 2.828e-4 | 0.19 |

Net exterior-current residual is approximately `2.5e-9` at N=80. The
assembled-flux RMS error is approximately 2.8--2.9% of the manufactured current
`alpha = 0.01`.

## Interpretation and limitations

The evidence supports the following conclusions:

1. Heart and bath potentials converge near second order.
2. Intracellular no-leak converges near second order.
3. Independently reconstructed flux continuity converges at approximately
   order 1.8.
4. The algebraic face flux is locally single-valued and global current balance
   is excellent.
5. Constitutive accuracy of the assembled face flux is bounded below 5% but
   converges slowly, with fitted order about 0.19.
6. Reconstructed heart flux is non-monotone from N=40 to N=80, although its
   four-level fitted order is 1.6--1.7. Exact-field reconstruction converges at
   approximately third order over that interval, so this is numerical local-
   gradient sensitivity rather than postprocessor failure.
7. Bath-side reconstructed flux remains near a `1e-4` floor.

The paper must not claim one universal bath-bidomain convergence order. It must
report the order of each quantity and distinguish exact algebraic conservation
from constitutive-flux accuracy.

## Paper-ready claim

> A domain-consistent finite-volume assembly uses a common myocardium operator
> for both intracellular potential contributions while retaining a single
> conservative extracellular heart--bath flux on a conformal mesh. Across four
> tetrahedral refinements, extracellular potentials converge at approximately
> order 1.75, reconstructed interface-current continuity at approximately order
> 1.8, and intracellular no-leak at approximately second order. The assembled
> face current remains globally conservative and below 3% RMS error at the
> finest resolution, but displays low-order, mesh-sensitive constitutive
> accuracy.

## Reproduction

Serial matched sweep:

```bash
RESOLUTIONS="10 20 40" bash setup/run_matched_serial_sweep.sh
```

Parallel N=80:

```bash
ASSEMBLY=matchedSubmesh \
METHODS=distanceWeightedHarmonic \
RESOLUTIONS=80 \
NPROCS=6 \
bash setup/run_parallel_interface_sweep.sh
```

Same-mesh parallel equivalence:

```bash
N=40 NPROCS=6 \
METHOD=distanceWeightedHarmonic \
ASSEMBLY=matchedSubmesh \
bash setup/run_parallel_equivalence.sh
```

Final numerical and exact CSVs are under
`setup/interfaceStudy/matchedSubmesh/distanceWeightedHarmonic/N80/`. Frozen
mesh hashes are under `setup/interfaceMeshBank/`; the reproducible compressed
mesh caches are generated locally and excluded from version control.

## Retained controls

The following are intentionally retained:

- `currentSplit`: backward compatibility and negative control;
- `unweightedHarmonic`: equal-distance/control comparison;
- `naiveLinearSigmaTotal`: deliberately incorrect physics control;
- limited/orthogonal/Gauss scheme scripts: rejected-method evidence;
- exact-field residual diagnostics: proof that measurement reconstruction is
  not the source of the observed local-gradient behavior.

They are not selected production settings for this paper tutorial.
