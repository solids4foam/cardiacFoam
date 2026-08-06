# Manufactured Eikonal ECG

This case verifies the eikonal activation-time solve and the eikonal template
ECG calculation on manufactured unit domains.

The eikonal solver writes `psi` as the activation time. The ECG solver samples
`Vm(x,t) = U(t - psi(x))` internally and writes `postProcessing/eikonalECG.dat`.

Run one 3D case directly with:

```sh
blockMesh -dict system/blockMeshDict.3D
./Allrun
```

Run the available manufactured dimensions through the setup directory with:

```sh
./setup/run_all_dimensions.sh
```

Expected verification outputs:

- `postProcessing/manufacturedEikonalActivationTime.dat`
- `postProcessing/eikonalECG.dat`
- `postProcessing/manufacturedEikonalECG.dat`
- `postProcessing/manufacturedEikonalECGSummary.dat`

The post-processing helper in `setup` collects activation
and ECG summary files into CSV tables for manufactured mesh studies.

After the driver runs, generate the convergence plot with:

```sh
python3 setup/plot_convergence.py postProcessing
```

This writes `postProcessing/convergence_plot.pdf` and `.png` with three panels:
spacing h. Cases that hit the nonlinear solver iteration cap are marked with ×.

## Error Calculation and Quadrature

It is important to clarify how the errors are evaluated for the different fields in this verification suite:

1. **Activation Times ($\tau$)**: No quadrature is used here. Because the manufactured activation time is a simple analytical function $\tau(x) = \exp(k \cdot x)$, it can be evaluated exactly at any point. To calculate the error, the solver performs a point-by-point comparison between the numerical activation time solved at each OpenFOAM mesh cell's center and the exact analytical mathematical formula evaluated at that identical cell center point.

2. **The ECG Computation (Where Quadrature is Used)**: Unlike the activation time, the ECG signal is defined mathematically as a volume integral over the entire domain. To get the "exact" baseline reference to compare OpenFOAM against, we calculate the integral of the manufactured analytical gradient field. Because this specific multidimensional integral does not have a simple closed-form algebraic solution, the exact "analytical" reference integral is computed using a highly accurate Gauss-Legendre Quadrature (e.g. $q=96$), which integrates the continuous analytical function down to machine precision. The ECG error is the comparison between OpenFOAM's numerical mesh integration (which simply sums up the cell values $\times$ cell volumes) against this near-perfect reference integral computed using the Gauss-Legendre quadrature.

## Tetrahedral (unstructured) mesh variant

`setup/mesh/tet/` is an activatable overlay of this same case on a genuinely
unstructured mesh: identical `constant/` and `system/` dicts, except the mesh
generator changes and `setup/mesh/tet/fvSolution` is swapped in for the
duration of a tet run and restored on exit. It was formerly the standalone
`eikonalTetMMS` tutorial, merged in here (the 5 shared dicts were
byte-identical, so this case stayed the canonical hex case unchanged; only the
tet-specific `box.geo.template`, `fvSolution`, and run scripts became the
overlay). The same pattern was later applied to
`monodomainPseudoECG/setup/mesh/tet/` (formerly `monodomainTetMMS`).

Run the gradient-scheme convergence sweep (paper table):

```bash
cd tutorials/manufacturedSolutions/eikonalECG
bash setup/mesh/tet/run_scheme_study.sh
```

Runs `Gauss linear` and `leastSquares` gradient reconstruction across
`N = 10, 20, 40` (Gauss--linear) and `N = 10, 20, 40, 80` (least-squares,
overridable via `RESOLUTIONS`), writes `setup/results/scheme_study.csv`, and
persists the canonical Paper I table:

```bash
python3 applications/scripts/paperI_results/aggregate.py eikonal_tet
```

The mesh-quality-only sweep (`setup/mesh/tet/run_eikonal_tet.sh`) and its
resolution/end-time overrides follow the same convention as
`monodomainPseudoECG`'s tet variant; see that case's README for the full
description of the mesh-generation and effective-spacing methodology, which
this case reuses unchanged (same `box.geo.template`,
`setup/mesh/tet/summarize_tet.py` convergence-order convention).
