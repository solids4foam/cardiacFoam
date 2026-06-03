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
./setupManufacturedEikonalECG/run_all_dimensions.sh
```

Expected verification outputs:

- `postProcessing/manufacturedEikonalActivationTime.dat`
- `postProcessing/eikonalECG.dat`
- `postProcessing/manufacturedEikonalECG.dat`
- `postProcessing/manufacturedEikonalECGSummary.dat`

The post-processing helper in `setupManufacturedEikonalECG` collects activation
and ECG summary files into CSV tables for manufactured mesh studies.

After the driver runs, generate the convergence plot with:

```sh
python3 setupManufacturedEikonalECG/plot_convergence.py postProcessing
```

This writes `postProcessing/convergence_plot.pdf` and `.png` with three panels:
activation-time error, ECG Linf error, and PIMPLE iteration count — all vs mesh
spacing h. Cases that hit the nonlinear solver iteration cap are marked with ×.
