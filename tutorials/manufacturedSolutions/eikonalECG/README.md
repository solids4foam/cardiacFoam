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
