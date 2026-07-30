# Rotated-anisotropy monodomain manufactured solution

This case verifies the three-dimensional monodomain diffusion operator on the
unit cube with a constant, symmetric positive-definite conductivity tensor
whose principal directions are rotated relative to the Cartesian mesh.

The exact voltage is

\[
V(x,y,z,t)=\sqrt{1+t}\,
\sin^2(\pi x)\sin^2(\pi y)\sin^2(\pi z).
\]

Its gradient vanishes on every cube face, so the solver's homogeneous-flux
boundary condition is exact even though the conductivity has non-zero
off-diagonal components. The verification model adds the continuous analytic
source

\[
S=\beta V-\nabla\cdot(K\nabla V)
\]

at each PDE step. The source contains the complete mixed-Hessian contribution
and is not evaluated with a finite-volume operator.

Run the Cartesian smoke case after sourcing OpenFOAM:

```bash
./Allrun
```

Set `N_CELLS` to change the number of cells per direction:

```bash
N_CELLS=20 ./Allrun
```

The final volume-weighted error norms are written under `postProcessing/`.

For the paper convergence control, run the dedicated sweeps:

```bash
bash setup/run_rotated_anisotropy_hex.sh
bash setup/run_rotated_anisotropy_tet.sh
```

They write `setup/results/rotated_anisotropy_hex_convergence.csv` and
`setup/results/rotated_anisotropy_tet_convergence.csv`. The observed rates
match the existing monodomain convergence interpretation: Cartesian hexes are
second order, least-squares Delaunay tetrahedra are near second order, and
Gauss-linear tetrahedra remain low-order.

## Symbolic source audit

The production solver evaluates the analytic expression directly in C++. An
independent SymPy script derives

\[
\nabla\cdot(K\nabla V)
\]

for a constant symmetric tensor and proves that the resulting manufactured
source is algebraically identical to the C++ expression:

```bash
python3 derive_manufactured_source.py
```

The script is an audit utility only; it is not called by the solver and does
not use a discrete finite-volume operator.
